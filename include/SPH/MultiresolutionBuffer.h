#pragma once

#include <climits>
#include <string>
#include "ParticleSet.h"
#include <TNL/Particles/GhostZone.h>
#include "OpenBoundaryBuffers.h"
#include "SPHTraits.h"

//FIXME:
#include "shared/Interpolation.h"
#include "shared/WendlandC2ABFs.h"

namespace TNL {
namespace SPH {

template< typename SPHCaseConfig >
class MassNodes
{
public:

   using DeviceType = typename SPHCaseConfig::DeviceType;
   using SolverTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using IndexType = typename SolverTraitsType::GlobalIndexType;
   using RealType = typename SolverTraitsType::RealType;
   using VectorType = typename SolverTraitsType::VectorType;
   using IndexArrayType = typename SolverTraitsType::IndexArrayType;
   using ScalarArrayType = typename SolverTraitsType::ScalarArrayType;
   using VectorArrayType = typename SolverTraitsType::VectorArrayType;

   void
   setSize( const IndexType& size )
   {
      points.setSize();
      normal.setSize();
      massFlux.setSize();
      mass.setSize();
      particlesToCreate.setSize();
   }

   void
   sort()
   {
      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::sort_by_key( thrustDevice,
                           particlesToCreate.getData(),
                           particlesToCreate.getData() + this->numberOfMassNodes,
                           thrust::make_zip_iterator( thrust::make_tuple( points.getArrayData(),
                                                                          normal.getArrayData(),
                                                                          massFlux.getArrayData(),
                                                                          mass.getArrayData() ) ) );
   }

   IndexType numberOfMassNodes;

   VectorArrayType points;
   VectorArrayType normal;

   ScalarArrayType massFlux;
   ScalarArrayType mass;
   IndexArrayType  particlesToCreate;
};

template< typename ParticlesType,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables,
          typename OpenBoundaryConfig,
          typename SPHDefs >
// TODO: Derive from ParticleSet or OpenBoundaryBuffer?
class MultiresolutionBoundary : public ParticleSet< ParticlesType, SPHCaseConfig, Variables, IntegratorVariables >
{
public:

   using DeviceType = typename SPHCaseConfig::DeviceType;
   using SolverTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using IndexType = typename SolverTraitsType::GlobalIndexType;
   using RealType = typename SolverTraitsType::RealType;
   using VectorType = typename SolverTraitsType::VectorType;
   using IndexArrayType = typename SolverTraitsType::IndexArrayType;

   using ParticleZone = TNL::ParticleSystem::ParticleZone< typename ParticlesType::Config, typename ParticlesType::DeviceType >;
   using MassNodes = MassNodes< SPHCaseConfig >;

   // matrix and vector type for interpolation
   // FIXME: Hard-def ABFs

   using ABFs = Interpolation::WendlandC2ABFs< 2, 2, SPHCaseConfig >;
   using KernelFunction = typename SPHDefs::KernelFunction;
   using MFD = Interpolation::MFD< 2, 2, RealType, ABFs, KernelFunction >;

   using MfdVectorType = typename MFD::BaseVectorType;
   using MfdMatrixType = typename MFD::BaseMatrixType;

   MultiresolutionBoundary() = default;

   void initialize()
   {

   };

   //compute mass fluxes and create new buffer particles
   template< typename FluidPointer, typename ModelParams >
   void
   accumulateMasses( FluidPointer& fluid_neihgbor, ModelParams& modelParams, const RealType dt )
   {
      auto searchInFluid = fluid_neihgbor->getParticles()->getSearchToken( fluid_neihgbor->getParticles() );

      //const auto view_points_massNodes = this->getParticles()->getPoints().getConstView();
      const auto view_points_massNodes = this->massNodes.points.getConstView();
      const auto view_normals_massNodes = this->massNodes.normal.getConstView();
      //auto view_rho_massPoints = this->getVariables()->rho.getView();
      //auto view_v_massPoints = this->getVariables()->v.getView();
      const auto view_m_massPoints = massNodes->getMass().getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const RealType h = modelParams.h;
      const RealType m = modelParams.m;
      const RealType dx = modelParams.dp; //FIXME
      const VectorType v_subdomain = 0.f; // velocity of moving subdomain
      const RealType div_r_trashold = modelParams.div_r_trashold;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius(); //TODO: Is it this one?
      const RealType extrapolationDetTreshold = modelParams.extrapolationDetTreshold;

      auto interpolateFluid = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_x,
            MfdMatrixType* M_x,
            MfdVectorType* brho_x,
            MfdVectorType* bvx_x,
            MfdVectorType* bvy_x,
            MfdVectorType* bvz_x,
            VectorType* div_r_x ) mutable
      {
         const VectorType r_j = view_points_fluid[ j ];
         const VectorType r_xj = r_x - r_j;
         const RealType drs = l2Norm( r_xj );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho_fluid[ j ];
            const RealType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const auto b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            *bvx_x += v_j[ 0 ] * b_x;
            *bvy_x += v_j[ 1 ] * b_x;
            *bvz_x += v_j[ 2 ] * b_x;

            //div r
            const VectorType gradW = r_xj * KernelFunction::F( drs, h );
            *div_r_x += gradW * V_j;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_x = view_points_massNodes[ i ];
         const VectorType normal_x = view_normals_massNodes[ i ];
         MfdMatrixType M_x = 0.f;
         MfdVectorType brho_x = 0.f;
         MfdVectorType bvx_x = 0.f;
         MfdVectorType bvy_x = 0.f;
         MfdVectorType bvz_x = 0.f;
         VectorType div_r_x = 0.f;

         ParticlesType::NeighborsLoop::exec(
               i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bvx_x, &bvy_x, &bvz_x, &div_r_x );

         RealType rho_x;
         RealType vx_x;
         RealType vy_x;
         RealType vz_x;

         if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            vx_x  = Matrices::solve( M_x, vx_x   )[ 0 ];
            vy_x  = Matrices::solve( M_x, vy_x   )[ 0 ];
            vz_x  = Matrices::solve( M_x, vx_x   )[ 0 ];
         }
         else if( M_x( 0, 0 ) > 0.f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            vx_x =  vx_x[ 0 ]   / M_x( 0, 0 );
            vy_x =  vy_x[ 0 ]   / M_x( 0, 0 );
            vz_x =  vx_x[ 0 ]   / M_x( 0, 0 );
         }
         else{
            // not sure what to do here
         }

         RealType m_flux_x;
         const VectorType v_x = { vx_x, vy_x, vz_x };
         const RealType m_x_lessZero = TNL::max( 0, -rho_x * ( ( v_x - v_subdomain ), normal_x ) * dx * dt );
         const RealType m_x_geqZero = -rho_x * ( ( v_x - v_subdomain ), normal_x ) * dx * dt;

         // free surface correction
         if( div_r_x < div_r_trashold )
            m_flux_x = 0.f;

         view_m_massPoints[ i ] = m_flux_x;
      };
      massNodes->getParticles()->forAll( particleLoop );
   }

   //interpolate values to buffer particles
   template< typename FluidPointer, typename ModelParams >
   void
   interpolateVariables( FluidPointer& fluid_neihgbor, ModelParams& modelParams )
   {
      //auto searchInFluid = this->getParticles()->getSearchToken( fluid_neihgbor->getParticles() )
      typename ParticlesType::NeighborsLoopParams searchInFluid( fluid_neihgbor->getParticles() );

      const auto view_points_overlap = this->getParticles()->getPoints().getConstView();
      auto view_rho_overlap = this->getVariables()->rho.getView();
      auto view_v_overlap = this->getVariables()->v.getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const RealType h = modelParams.h;
      const RealType m = modelParams.m;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius(); //TODO: Is it this one?
      const RealType extrapolationDetTreshold = modelParams.extrapolationDetTreshold;

      auto interpolateFluid = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_x,
            MfdMatrixType* M_x,
            MfdVectorType* brho_x,
            MfdVectorType* bvx_x,
            MfdVectorType* bvy_x,
            MfdVectorType* bvz_x ) mutable
      {
         const VectorType r_j = view_points_fluid[ j ];
         const VectorType r_xj = r_x - r_j;
         const RealType drs = l2Norm( r_xj );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho_fluid[ j ];
            const RealType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const auto b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            *bvx_x += v_j[ 0 ] * b_x;
            *bvy_x += v_j[ 1 ] * b_x;
            *bvz_x += v_j[ 2 ] * b_x;
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_x = view_points_overlap[ i ];
         MfdMatrixType M_x = 0.f;
         MfdVectorType brho_x = 0.f;
         MfdVectorType bvx_x = 0.f;
         MfdVectorType bvy_x = 0.f;
         MfdVectorType bvz_x = 0.f;

         ParticlesType::NeighborsLoop::exec( i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bvx_x, &bvy_x, &bvz_x );

         if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold ) {
            const RealType rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            const RealType vx_x =  Matrices::solve( M_x, vx_x   )[ 0 ];
            const RealType vy_x =  Matrices::solve( M_x, vy_x   )[ 0 ];
            const RealType vz_x =  Matrices::solve( M_x, vx_x   )[ 0 ];

            view_rho_overlap[ i ] = rho_x;
            view_v_overlap[ i ] = { vx_x, vy_x, vz_x };
         }
         else if( M_x( 0, 0 ) > 0.f ) {
            const RealType rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            const RealType vx_x =  vx_x[ 0 ]   / M_x( 0, 0 );
            const RealType vy_x =  vy_x[ 0 ]   / M_x( 0, 0 );
            const RealType vz_x =  vx_x[ 0 ]   / M_x( 0, 0 );

            view_rho_overlap[ i ] = rho_x;
            view_v_overlap[ i ] = { vx_x, vy_x, vz_x };
         }
         else{
            // not sure what to do here
         }
      };
      this->getParticles()->forAll( particleLoop );
   }

   //shift the buffer particles

   // open bc logic
   // - fluid into buffer -> buffer
   // - buffer out of buffer -> remove
   // - buffer into fluid -> fluid
   void
   moveBufferParticles( const RealType dt )
   {
      auto view_r_buffer = this->getParticles()->getPoints().getView();
      const auto view_v_buffer = this->getVariables()->v.getConstView();

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const VectorType bufferPosition = this->bufferPosition;
      const VectorType bufferOrientation = this->bufferOrientation;
      const RealType bufferWidth = this->bufferWidth;

      auto retypeMarker_view = retypeMarker.getView();
      auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;
      };
      this->getParticles()->forAll( moveBufferParticles );

      // identify particles moving to fluid
      auto identifyParticlesToRetype = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;
         if( ( r_relative, bufferOrientation ) <= 0.f ){
            retypeMarker_view[ i ] = -1;
            return 1;
         }
         else
            return 0;

      };
      this->numOfPtcsToRetype = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, identifyParticlesToRetype );

      // identify particles moving out of the buffer
      auto identifyPtcsToRemove = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;
         if( ( r_relative, bufferOrientation ) > bufferWidth ) {
            retypeMarker_view[ i ] = 1;
            return 1;
         }
         else
            return 0;
      };
      this->numOfPtcsToRemove = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, identifyParticlesToRetype );
   }

   void
   sortBufferParticles()
   {
      const IndexType numberOfBufferPtcs = this->getNumberOfParticles();

      auto r_view = this->getParticles()->getPoints().getView();
      auto v_view = this->getVariables()->v.getView();
      auto rho_view = this->getVariables()->rho.getView();
      auto retypeMarker_view = this->retypeMarker.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;

      thrust::sort_by_key( thrustDevice,
                           retypeMarker_view.getArrayData(),
                           retypeMarker_view.getArrayData() + numberOfBufferPtcs,
                           thrust::make_zip_iterator( thrust::make_tuple( r_view.getArrayData(),
                                                                          v_view.getArrayData(),
                                                                          rho_view.getArrayData() ) ) );
   }

   template< typename FluidPointer >
   void
   convertBufferToFluid( FluidPointer& fluid )
   {
      const IndexType numberOfFluidPtcs = fluid->getParticles()->getNumberOfParticles();

      auto view_r_fluid = fluid->getParticles()->getPoints().getView();
      auto view_v_fluid = fluid->getVariables()->v.getView();
      auto view_rho_fluid = fluid->getVariables()->rho.getView();
      auto view_rho_old = fluid->getIntegratorVariables()->rho_old.getView();
      auto view_v_old = fluid->getIntegratorVariables()->v_old.getView();

      auto view_r_buffer = this->getParticles()->getPoints().getView();
      auto view_v_buffer = this->getVariables()->v.getView();
      auto view_rho_buffer = this->getVariables()->rho.getView();

      auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_fluid[ numberOfFluidPtcs + i ] = view_r_buffer[ i ];
         view_rho_fluid[ numberOfFluidPtcs + i ] = view_rho_buffer[ i ];
         view_v_fluid[ numberOfFluidPtcs + i ] = view_v_buffer[ i ];

         view_rho_old[ numberOfFluidPtcs + i ] = view_rho_buffer[ i ];
         view_v_old[ numberOfFluidPtcs + i ] = view_v_buffer[ i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numOfPtcsToRetype, createNewFluidParticles );


   }

   void
   removeBufferParticles()
   {
      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs - this->numOfPtcsToRemove );
   }

   template< typename ModelParams >
   void
   updateMassNodes( ModelParams& modelParams, const RealType dt )
   {
      const auto massFlux_view = massNodes->massFlux.getConstView();
      auto mass_view = massNodes->mass.getView();
      auto view_particlesToCreate = massNodes->particlesToCreate.getView();

      const RealType particleMass = modelParams.mass;
      const IndexType numberOfBufferPartices = this->getParticles()->getNumberOfParticles();

      // reset list with markers
      view_particlesToCreate = 0;

      auto identifyParticlesToCreate = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         mass_view[ i ] = dt * massFlux_view[ i ];
         if( mass_view[ i ] > particleMass )
         {
            mass_view[ i ] -= particleMass;
            view_particlesToCreate[ i ] = 1;
            return 1;
         }
         else
            return 0;
      };
      this->numberOfPtcsToCreate = Algorithms::reduce< DeviceType >( 0, numberOfBufferPartices, identifyParticlesToCreate );
   }

   template< typename ModelParams >
   void
   createBufferParticles( ModelParams& modelParams )
   {
      auto r_buffer_view = this->getParticles()->getPoints().getView();
      const auto r_massNodes_view = this->massNodes.points.getConstView();
      const auto normal_massNodes_view = this->massNodes.normal.getConstView();

      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();
      const RealType dp = modelParams.dp;

      auto createNewBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         const VectorType r_new = r_massNodes_view[ i ] + ( 0.5f * dp ) * normal_massNodes_view[ i ];
         r_buffer_view[ numberOfBufferPtcs + i ] = r_new;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfPtcsToCreate, createNewBufferParticles );

      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs + this->numberOfPtcsToCreate );
   }

   template< typename FluidPointer >
   void
   getFluidParticlesEneringTheBuffer( FluidPointer& fluid )
   {
      auto particlesToBuffer_view = this->particlesToBuffer.getView();
      particlesToBuffer_view = INT_MAX;

      const auto r_view = fluid->getParticles()->getPoints()->getConstView();
      const auto zoneParticleIndices_view = this->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = this->zone.getNumberOfParticles();

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const VectorType bufferPosition = this->bufferPosition;
      const VectorType bufferOrientation = this->bufferOrientation;

      auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         const VectorType r_relative = bufferPosition - r;

         if( ( r_relative, bufferOrientation ) > 0 ){
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         return 0;
      };
      this->numberOfPtcsToBuffer = Algorithms::reduce< DeviceType >(
            0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );

      // sort the indices
      const IndexType rangeToSort = ( numberOfZoneParticles > numberOfBufferParticles ) ? numberOfZoneParticles : numberOfBufferParticles;
      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::sort( thrustDevice,
                    particlesToBuffer_view.getArrayData(),
                    particlesToBuffer_view.getArrayData() + rangeToSort );

   }

   template< typename FluidPointer >
   void
   convertFluidToBuffer( FluidPointer& fluid )
   {
      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();

      //FIXME: I'm sure that it is enough to copy the positioin.
      auto r_fluid_view = fluid->getParticles()->getPoints().getConstView();
      const auto v_fluid_view = fluid->getVariables()->v.getConstView();
      const auto rho_fluid_view = fluid->getVariables()->rho.getConstView();

      auto r_buffer_view = this->getParticles()->getPoints().getView();
      auto v_buffer_view = this->getVariables()->v.getView();
      auto rho_buffer_view = this->getVariables()->rho.getView();
      const auto particlesToBuffer_view = this->particlesToBuffer.getConstView();

      auto retypeFluidToBuffer = [=] __cuda_callable__ ( int i ) mutable
      {
         const IndexType p = particlesToBuffer_view[ i ];

         r_buffer_view[ numberOfBufferPtcs + i ] = r_fluid_view[ p ];
         v_buffer_view[ numberOfBufferPtcs + i ] = v_fluid_view[ p ];
         rho_buffer_view[ numberOfBufferPtcs + i ] = rho_fluid_view[ p ];

         r_fluid_view[ p ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, fluidToBufferCount, retypeFluidToBuffer );

      // update the particles counts
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs + fluidToBufferCount );
      fluid->getParticles()->setNumberOfParticlesToRemove(
            fluid->getParticles()->getNumberOfParticlesToRemove() + fluidToBufferCount );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   updateInterfaceBuffer( FluidPointer& fluid, FluidPointer& fluid_neihgbor, ModelParams& modelParams, const RealType dt )
   {
      // update buffer: buffer -> fluid, buffer -> remove
      moveBufferParticles();
      sortBufferParticles();
      removeBufferParticles();
      convertBufferToFluid( fluid );

      // fluid -> buffer
      getFluidParticlesEneringTheBuffer( fluid );
      convertFluidToBuffer( fluid );

      // mass-nodes -> buffer
      updateMassNodes();
      massNodes->sort();
      createBufferParticles();

      // update the values
      accumulateMasses( fluid_neihgbor, modelParams, dt ); // ?
      interpolateVariables( fluid_neihgbor, modelParams );

      // reset temorary variables and fields
      numberOfPtcsToRemove = 0;
      numberOfPtcsToRetype = 0;
      numberOfPtcsToCreate = 0;
      numberOfPtcsToBuffer = 0;

   }

   int ownerSetID;
   int neighborSetID;

   MassNodes massNodes;

   // temporary constants
   IndexType numberOfPtcsToRemove = 0;
   IndexType numberOfPtcsToRetype = 0;
   IndexType numberOfPtcsToCreate = 0;
   IndexType numberOfPtcsToBuffer = 0;

   IndexType fluidToBufferCount;

   // temporary arrays
   IndexArrayType retypeMarker;


   IndexArrayType particlesToFluid;
   IndexArrayType particlesToRemove;
   IndexArrayType particlesToBuffer;

   // buffer referential position (TODO: Use openbc config?)
   VectorType bufferPosition;
   VectorType bufferOrientation; //FIXME: temp
   RealType bufferWidth; // FIXME: temp

   //
   ParticleZone zone;

};

} // SPH
} // TNL

