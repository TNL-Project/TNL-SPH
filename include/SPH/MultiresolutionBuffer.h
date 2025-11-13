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

//FIXME - debug use MDBC corrections
#include "Models/WCSPH_DBC/details.h"

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

   MassNodes() = default;

   void
   setSize( const IndexType& size )
   {
      points.setSize( size );
      normal.setSize( size );
      massFlux.setSize( size );
      mass.setSize( size );
      particlesToCreate.setSize( size );

      numberOfMassNodes = size;
      mass = 0.f;
      massFlux = 0.f;
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
   using BaseType = ParticleSet< ParticlesType, SPHCaseConfig, Variables, IntegratorVariables >;

   using DeviceType = typename SPHCaseConfig::DeviceType;
   using SolverTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using IndexType = typename SolverTraitsType::GlobalIndexType;
   using RealType = typename SolverTraitsType::RealType;
   using VectorType = typename SolverTraitsType::VectorType;
   using IndexArrayType = typename SolverTraitsType::IndexArrayType;

   using ParticleZone = TNL::ParticleSystem::ParticleZone< typename ParticlesType::Config, typename ParticlesType::DeviceType >;
   using MassNodes = MassNodes< SPHCaseConfig >;

   //TODO: Due to init
   using IndexVectorType = typename SolverTraitsType::IndexVectorType;

   // matrix and vector type for interpolation
   // FIXME: Hard-def ABFs

   //using ABFs = Interpolation::WendlandC2ABFs< 2, 2, SPHCaseConfig >;
   using KernelFunction = typename SPHDefs::KernelFunction;
   using MFD = Interpolation::MFD< 2, 1, RealType, Interpolation::WendlandC2ABFs, KernelFunction, SPHCaseConfig >;
   using ABFs = typename MFD::ABFs;

   using MfdVectorType = typename MFD::BaseVectorType;
   using MfdMatrixType = typename MFD::BaseMatrixType;

   MultiresolutionBoundary() = default;


   // -----------------------------------------------------------------------------------------------------------------------
   void
   initZones( const IndexVectorType zoneOriginIdx_left,
              const IndexVectorType zoneDimensions_left,
              const IndexVectorType zoneOriginIdx_right,
              const IndexVectorType zoneDimensions_right,
              const IndexVectorType gridDimensionsWithOverlap, //FIXME: Which of two grids is this one?
              const int subdomainIdx,
              const IndexType numberOfParticlesPerCell = 75 )
   {
      //dp get from arguments
      //const IndexVectorType subdomainGridDimension = this->getParticles()->getGridDimensions();
      //const IndexVectorType subdomainGridDimensionWithOverlap = this->getParticles()->getGridDimensionsWithOverlap();
      //const int numberOfOverlapLayers = this->getParticles()->getNumberOfOverlapLayers();
      //const IndexVectorType zoneOriginIdx_left = { 0, 0 };
      //const IndexVectorType zoneDimensions_left = { 2, subdomainGridDimension[ 1 ] };
      //const IndexVectorType zoneOriginIdx_right = { subdomainGridDimension[ 0 ] + numberOfOverlapLayers, 0 };
      //const IndexVectorType zoneDimensions_right = { 2, subdomainGridDimension[ 1 ] };

      //FIXME: Two zones
      //zone_left.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      //zone_left.assignCells( zoneOriginIdx_left, zoneDimensions_left, gridDimensionsWithOverlap );

      //zone_right.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      //zone_right.assignCells( zoneOriginIdx_right, zoneDimensions_right, gridDimensionsWithOverlap );
      const RealType searchRadius = this->getParticles()->getSearchRadius();

      // FIXME FIXME Requires dp, dp is in model params, what to do!
      // FIXME FIXME just put dp to multiresolution file
      const RealType dp = 0.002;


      //FIXME: With two subdomains, we can start with signle zone
      if( subdomainIdx == 0 ){
         zone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
         zone.assignCells( zoneOriginIdx_right, zoneDimensions_right, gridDimensionsWithOverlap );
         this->bufferOrientation = { -1.f, 0 };
         this->bufferPosition = this->getParticles()->getGridOrigin() + this->getParticles()->getGridDimensions() * searchRadius;
         this->bufferWidth = 4 * 0.002;
      }
      else if( subdomainIdx == 1 ){
         zone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
         zone.assignCells( zoneOriginIdx_left, zoneDimensions_left, gridDimensionsWithOverlap );
         this->bufferOrientation = { 1.f, 0 };
         this->bufferPosition = this->getParticles()->getGridOrigin();
         this->bufferWidth = 4 * 0.001; // dp times refinement factor
      }

      const IndexType numberOfAllocatedParticles = this->getNumberOfAllocatedParticles();
      particlesToFluid.setSize( numberOfAllocatedParticles );
      particlesToRemove.setSize( numberOfAllocatedParticles );
      particlesToBuffer.setSize( numberOfAllocatedParticles );      //TODO: Set sizes
      retypeMarker.setSize( numberOfAllocatedParticles ); //TODO: rename


      //FIXME: Set variables on boundary particles so the paraview doesn't ouput nonsense

   };
   // -----------------------------------------------------------------------------------------------------------------------
   template< typename ModelParams >
   void
   initMassNodes( ModelParams& modelParams, const int subdomainIdx, const RealType refinemnetFactor )
   {
      std::cout << "=== INIT MASS NODES =================== " << std::endl;
      const RealType searchRadius = this->getParticles()->getSearchRadius();
      std::cout << "Search radius: " << searchRadius << std::endl;
      std::cout << "Grid origin: " << this->getParticles()->getGridOrigin() << std::endl;
      std::cout << "Grid origin with overlap: " << this->getParticles()->getGridOriginWithOverlap() << std::endl;
      std::cout << "Grid dimensions: " << this->getParticles()->getGridDimensions() << std::endl;
      std::cout << "Grid dimensions with overlap: " << this->getParticles()->getGridDimensionsWithOverlap() << std::endl;
      std::cout << "Domain size: " << this->getParticles()->getGridOrigin() + this->getParticles()->getGridDimensions() * searchRadius << std::endl;
      std::cout << "Domain size with overlap: " << this->getParticles()->getGridOriginWithOverlap() + this->getParticles()->getGridDimensionsWithOverlap() * searchRadius << std::endl;
      std::cout << "======================================= " << std::endl;

      const RealType height = 1.f;
      const RealType dp = refinemnetFactor * modelParams.dp;
      const int numberOfMassNodes = height / dp;
      massNodes.setSize( numberOfMassNodes );

      RealType massNodes_xCoord;
      VectorType massNodes_normal;
      if( subdomainIdx == 0 ){
         massNodes_xCoord = this->getParticles()->getGridOriginWithOverlap()[ 0 ] +
                            this->getParticles()->getGridDimensionsWithOverlap()[ 0 ] * searchRadius;
         massNodes_normal = { -1.f, 0 };
      }
      else if( subdomainIdx == 1 ){
         massNodes_xCoord = this->getParticles()->getGridOriginWithOverlap()[ 0 ];
         massNodes_normal = { 1.f, 0 };
      }
      else{
         std::cerr << "Invalid buffer count." << std::endl;
      }

      auto points_massNodes_view = this->massNodes.points.getView();
      auto normals_massNodes_view = this->massNodes.normal.getView();

      auto generateMassNodesCoordinates = [=] __cuda_callable__ ( int i ) mutable
      {
         const VectorType r = { massNodes_xCoord, dp * ( i + 1 ) };
         points_massNodes_view[ i ] = r;
         normals_massNodes_view[ i ] = massNodes_normal;
      };
      Algorithms::parallelFor< DeviceType >( 1, numberOfMassNodes, generateMassNodesCoordinates ); //FIXME 0 is right, 1 for debug
   }

   // -----------------------------------------------------------------------------------------------------------------------



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
      const auto view_m_massPoints = massNodes.mass.getView(); //TODO: temp due to condition
      auto view_mFlux_massPoints = massNodes.massFlux.getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const unsigned int dim = SPHCaseConfig::spaceDimension;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius(); //TODO: Is it this one?
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      const RealType h = refinementFactor * modelParams.h;
      const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
      const RealType dx = refinementFactor * modelParams.dp; //FIXME
      const VectorType v_subdomain = 0.f; // velocity of moving subdomain
      //const RealType div_r_trashold = modelParams.div_r_trashold;
      const RealType div_r_trashold = 1.5f; //FIXME add to model params
      const RealType extrapolationDetTreshold = modelParams.mdbcExtrapolationDetTreshold; //TODO: rename -> introduce extrapolationDetTrashold

      //debug mdbc kernels
      std::cout << "Refinemet factor: " << refinementFactor << std::endl;

      auto interpolateFluid = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_x,
            MfdMatrixType* M_x,
            MfdVectorType* brho_x,
            MfdVectorType* bvx_x,
            MfdVectorType* bvy_x,
            //MfdVectorType* bvz_x,
            RealType* div_r_x,
            RealType* drs_min ) mutable
      {
         const VectorType r_j = view_points_fluid[ j ];
         const VectorType r_xj = r_x - r_j;
         const RealType drs = l2Norm( r_xj );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho_fluid[ j ];
            const VectorType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const MfdVectorType b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            *bvx_x += v_j[ 0 ] * b_x;
            *bvy_x += v_j[ 1 ] * b_x;
            //*bvz_x += v_j[ 2 ] * b_x;

            // div r
            const VectorType gradW = r_xj * KernelFunction::F( drs, h );
            *div_r_x += (-1.f) * ( r_xj, gradW ) * V_j; //TODO: Added - sign, should be there?

            // find nearest neighbor
            *drs_min = std::min( *drs_min, drs );
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
         //MfdVectorType bvz_x = 0.f;
         RealType div_r_x = 0.f;
         RealType drs_min = FLT_MAX;

         ParticlesType::NeighborsLoop::exec(
               //i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bvx_x, &bvy_x, &bvz_x, &div_r_x, &drs_min );
               i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bvx_x, &bvy_x, &div_r_x, &drs_min );

         RealType rho_x;
         RealType vx_x;
         RealType vy_x;
         //RealType vz_x;

         //if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold ) {
         if( M_x( 0, 0 ) > 0.05 ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            vx_x  = Matrices::solve( M_x, bvx_x  )[ 0 ];
            vy_x  = Matrices::solve( M_x, bvy_x  )[ 0 ];
            //vz_x  = Matrices::solve( M_x, bvx_x  )[ 0 ];
         }
         else if( M_x( 0, 0 ) > 0.f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            vx_x =  bvx_x[ 0 ]  / M_x( 0, 0 );
            vy_x =  bvy_x[ 0 ]  / M_x( 0, 0 );
            //vz_x =  bvx_x[ 0 ]  / M_x( 0, 0 );
         }
         else{
            // not sure what to do here
         }

         //const VectorType v_x = { vx_x, vy_x, vz_x };
         const VectorType v_x = { vx_x, vy_x };
         const RealType m_x_lessZero = TNL::max( 0,rho_x * ( ( v_x - v_subdomain ), normal_x ) * dx * dt ); // NOTE: Ricci uses (-1) * rho, but I think I have different normal orientation
         const RealType m_x_geqZero = rho_x * ( ( v_x - v_subdomain ), normal_x ) * dx * dt; // NOTE: Ricci uses (-1) * rho, but I think I have different normal orientation

         // sum up the flux with free surface corrections
         RealType m_flux_x = 0;

         //if( ( div_r_x > div_r_trashold || drs_min < dx ) && ( view_mFlux_massPoints[ i ] < 0.f  ) ) //TODO: Original implementation uses only fluxes
         if( ( div_r_x > div_r_trashold || drs_min < dx ) && ( view_m_massPoints[ i ] < 0.f  ) )
            m_flux_x = m_x_lessZero;
         else if( div_r_x > div_r_trashold || drs_min < dx )
            m_flux_x = m_x_geqZero;
         else
            m_flux_x = 0;

         view_mFlux_massPoints[ i ] += m_flux_x;

         // Debug playground
         bool firstCond = ( div_r_x > div_r_trashold || drs_min < dx ) && ( view_m_massPoints[ i ] < 0.f  );
         bool secondCond = ( div_r_x > div_r_trashold || drs_min < dx );
         //if( M_x( 0, 0 ) > 0 )
         //   printf( "!!!!! r_x: %.3f, %.3f, c1: %d, c2: %d, M00: %.3f, m: %f, brho[0]: %f, rho_x %.2f, v_x %.3f, %.3f massFlux: %f, div_r: %.2f,   ( nx: %.2f,  ny: %.2f , dx: %.3f, dt %f ), detM: %f \n",
         //         r_x[ 0 ], r_x[ 1 ], firstCond, secondCond, M_x( 0, 0 ), view_m_massPoints[ i ], brho_x[ 0 ], rho_x, vx_x, vy_x, m_flux_x, div_r_x, normal_x[ 0 ], normal_x[ 1 ], dx, dt, Matrices::determinant( M_x ) );

         //if( M_x( 0, 0 ) > 0 )
         //printf( "fluxik: %f \n ", view_mFlux_massPoints[ i ] );
      };
      //massNodes->getParticles()->forAll( particleLoop );
      Algorithms::parallelFor< DeviceType >( 0, massNodes.numberOfMassNodes, particleLoop );
   }

   //interpolate values to buffer particles
   template< typename FluidPointer, typename ModelParams >
   void
   interpolateVariables( FluidPointer& fluid_neihgbor, ModelParams& modelParams )
   {
      //auto searchInFluid = this->getParticles()->getSearchToken( fluid_neihgbor->getParticles() )
      typename ParticlesType::NeighborsLoopParams searchInFluid( fluid_neihgbor->getParticles() );
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();

      std::cout << "[iV] numberOfBufferPtcs: " << numberOfBufferParticles << std::endl;

      const auto view_points_overlap = this->getParticles()->getPoints().getConstView();
      auto view_rho_overlap = this->getVariables()->rho.getView();
      auto view_v_overlap = this->getVariables()->v.getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      //const RealType h = modelParams.h;
      //const RealType m = modelParams.mass;
      //const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius(); //TODO: Is it this one?
      //const RealType extrapolationDetTreshold = modelParams.mdbcExtrapolationDetTreshold; //TODO: rename -> introduce extrapolationDetTrashold

      const unsigned int dim = SPHCaseConfig::spaceDimension;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius(); //TODO: Is it this one?
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      const RealType h = refinementFactor * modelParams.h;
      const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
      const RealType dx = refinementFactor * modelParams.dp; //FIXME

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
            const VectorType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const MfdVectorType b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            *bvx_x += v_j[ 0 ] * b_x;
            *bvy_x += v_j[ 1 ] * b_x;
            //*bvz_x += v_j[ 2 ] * b_x;
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
         printf( "| i: %d buffer r: %f, %f | \n", i, r_x[ 0 ], r_x[ 0 ] );

         ParticlesType::NeighborsLoop::exec( i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bvx_x, &bvy_x, &bvz_x );

         //if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold ) {
         RealType rho_x;
         RealType vx_x;
         RealType vy_x;

         //TODO: the conditions needs to be different
         const RealType detM = Matrices::determinant( M_x );
         if( M_x( 0, 0 ) > 0.05  && detM > 0.001f ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            vx_x =  Matrices::solve( M_x, bvx_x   )[ 0 ];
            vy_x =  Matrices::solve( M_x, bvy_x   )[ 0 ];
            //const RealType vz_x =  Matrices::solve( M_x, vx_x   )[ 0 ];

            view_rho_overlap[ i ] = rho_x;
            //view_v_overlap[ i ] = { vx_x, vy_x, vz_x };
            view_v_overlap[ i ] = { vx_x, vy_x };
         }
         else if( M_x( 0, 0 ) > 0.05f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            vx_x =  bvx_x[ 0 ]  / M_x( 0, 0 );
            vy_x =  bvy_x[ 0 ]  / M_x( 0, 0 );
            //const RealType vz_x =  bvz_x[ 0 ]   / M_x( 0, 0 );

            view_rho_overlap[ i ] = rho_x;
            //view_v_overlap[ i ] = { vx_x, vy_x, vz_x };
            view_v_overlap[ i ] = { vx_x, vy_x };
         }
         else{
            // FIXME: not sure what to do here
         }


         // Debug playground
         if( M_x( 0, 0 ) > 0 )
            printf( "!! intp. !! r_x: %.5f, %.5f, M00: %.3f, brho[0]: %f, rho_x %.2f, v_x %.3f, %.3f, detM: %f \n",
                  r_x[ 0 ], r_x[ 1 ], M_x( 0, 0 ), brho_x[ 0 ], rho_x, vx_x, vy_x,  detM );

      };
      //this->getParticles()->forAll( particleLoop );
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, particleLoop );
   }

   //shift the buffer particles

   // open bc logic
   // - fluid into buffer -> buffer
   // - buffer out of buffer -> remove
   // - buffer into fluid -> fluid
   void
   moveBufferParticles( const RealType dt, const int subdomain )
   {
      auto view_r_buffer = this->getParticles()->getPoints().getView();
      const auto view_v_buffer = this->getVariables()->v.getConstView();

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const VectorType bufferPosition = this->bufferPosition;
      const VectorType bufferOrientation = this->bufferOrientation;
      const RealType bufferWidth = this->bufferWidth;
      std::cout << "[mBP]: subdomain: " << subdomain << " nbp: " << numberOfBufferParticles << " bufferPosition: " << bufferPosition << " bufferOrientation: " << bufferOrientation << " bufferWidth: " << bufferWidth << std::endl;

      auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;

         //const VectorType r = view_r_buffer[ i ];
         //const VectorType r_relative = bufferPosition - r;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles );

      // reset retype marker
      auto retypeMarker_view = retypeMarker.getView();
      retypeMarker_view = 0;

      // identify particles moving to fluid
      auto identifyParticlesToRetype = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;
         if( ( r_relative, bufferOrientation ) <= 0.f ){
            printf(" <<r_new: %f, %f>>\n", r[ 0 ], r[ 1 ]);
            retypeMarker_view[ i ] = -1;
            return 1;
         }
         else
            return 0;

      };
      this->numberOfPtcsToRetype = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, identifyParticlesToRetype );

      // identify particles moving out of the buffer
      auto identifyParticlesToRemove = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;
         //printf( "---> r: %f, r_r: %f vx: %f \n", r[ 0 ], r_relative[ 0 ], view_v_buffer[ i ][ 0 ] );
         if( ( r_relative, bufferOrientation ) > bufferWidth ) {
            retypeMarker_view[ i ] = 1;
            //printf( " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ RBP \n" );
            return 1;
         }
         else
            return 0;
      };
      this->numberOfPtcsToRemove = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, identifyParticlesToRemove );
      std::cout << "[mBP]: numberOfPtcsToRetype: " << numberOfPtcsToRetype << " numberOfPtcsToRemove: " << numberOfPtcsToRemove << std::endl;
   }

   void
   sortBufferParticles()
   {
      const IndexType numberOfBufferPtcs = this->getNumberOfParticles();
      std::cout << "[sBF]: number of particles: " << numberOfBufferPtcs << std::endl;

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

      std::cout << "[cBTF]: numberOfPtcsToRetype: " << numberOfPtcsToRetype << std::endl;

      auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_fluid[ numberOfFluidPtcs + i ] = view_r_buffer[ i ];
         view_rho_fluid[ numberOfFluidPtcs + i ] = view_rho_buffer[ i ];
         view_v_fluid[ numberOfFluidPtcs + i ] = view_v_buffer[ i ];

         view_rho_old[ numberOfFluidPtcs + i ] = view_rho_buffer[ i ];
         view_v_old[ numberOfFluidPtcs + i ] = view_v_buffer[ i ];

         //TODO: I need to remove the used particles, not sure if this is the right way
         view_r_buffer[ i ] = FLT_MAX;

         //debug: print new fluid particles
         printf(" [r_new: %f, %f]\n", view_r_fluid[ numberOfFluidPtcs + i ][ 0 ], view_r_fluid[ numberOfFluidPtcs + i ][ 1 ]);
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfPtcsToRetype, createNewFluidParticles );
      fluid->getParticles()->setNumberOfParticles( numberOfFluidPtcs + this->numberOfPtcsToRetype );
      this->getParticles()->setNumberOfParticlesToRemove(
            this->getParticles()->getNumberOfParticlesToRemove() + this->numberOfPtcsToRetype );

      std::cout << "[cBTF][end]: " <<
         " fluid n: " << fluid->getNumberOfParticles() <<
         " mr-buffer n: " << this->getNumberOfParticles() <<
         " numberOfPtcsToRetype: " << this->numberOfPtcsToRetype <<
         " mr-buffer to remove:" << this->getParticles()->getNumberOfParticlesToRemove() << std::endl;
   }

   void
   removeBufferParticles()
   {
      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();
      std::cout << "[rBP] numberOfParticlesToRemove: " << this->numberOfPtcsToRemove << std::endl;
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs - this->numberOfPtcsToRemove );
   }

   template< typename ModelParams >
   void
   updateMassNodes( ModelParams& modelParams, const RealType dt )
   {
      const auto massFlux_view = massNodes.massFlux.getConstView();
      auto mass_view = massNodes.mass.getView();
      auto view_particlesToCreate = massNodes.particlesToCreate.getView();

      //FIXME: I need refinemet factor here! Subdomain mass is different
      const RealType particleMass = 0.25f * modelParams.mass ;
      const IndexType numberOfMassNodes = this->massNodes.numberOfMassNodes;

      // debug
      //const auto massNodes_coords = this->massNodes.points.getConstView();

      // reset list with markers
      view_particlesToCreate = 0;

      auto identifyParticlesToCreate = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         //mass_view[ i ] += dt * massFlux_view[ i ];
         mass_view[ i ] += massFlux_view[ i ]; // dt is already added
         //*****
         //if( massFlux_view[ i ] != 0.f )
         //   printf( " ================================================================ m: %f, mass_view: %f, mass_flux: %f \n", particleMass, mass_view[ i ], massFlux_view[ i ] );
         //bool myCond = mass_view[ i ] > particleMass;
         //if( mass_view[ i ] != 0.f )
         //   printf( " ================================================================ m: %f, x: %.3f, y: %.3f , cond: %d, mass_view: %f, mass_flux: %f \n",
         //         particleMass, massNodes_coords[ i ][ 0 ], massNodes_coords[ i ][ 1 ], myCond, mass_view[ i ], massFlux_view[ i ] );
         //*****
         if( mass_view[ i ] > particleMass )
         {
            mass_view[ i ] -= particleMass;
            //view_particlesToCreate[ i ] = 1;
            view_particlesToCreate[ i ] = -1;
            return 1;
         }
         else
            return 0;
      };
      this->numberOfPtcsToCreate = Algorithms::reduce< DeviceType >( 0, numberOfMassNodes, identifyParticlesToCreate );
      std::cout << "[uMN] END - numberOfParticles: " << this->getNumberOfParticles() << " numberOfParticleToCreate: " << this->numberOfPtcsToCreate << " numberOfMassNodes: " << massNodes.numberOfMassNodes << std::endl;
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

      std::cout << "[cBP]: numberOfParticlesToCreate: " << numberOfPtcsToCreate << std::endl;

      auto createNewBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         const VectorType r_new = r_massNodes_view[ i ] + ( 0.5f * dp ) * normal_massNodes_view[ i ];
         r_buffer_view[ numberOfBufferPtcs + i ] = r_new;
         printf(" ** created r_new: %f, %f ** \n", r_new[ 0 ], r_new[ 1 ]);
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfPtcsToCreate, createNewBufferParticles );

      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs + this->numberOfPtcsToCreate );
      std::cout << "[cBP] END - number of buffer particles : " << this->getNumberOfParticles() << std::endl;
   }

   template< typename FluidPointer >
   void
   getFluidParticlesEneringTheBuffer( FluidPointer& fluid )
   {
      auto particlesToBuffer_view = this->particlesToBuffer.getView();
      particlesToBuffer_view = INT_MAX;

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto zoneParticleIndices_view = this->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = this->zone.getNumberOfParticles();
      std::cout << "[gFPETB] number of zone particles: " << numberOfZoneParticles << " n: " <<  bufferOrientation << " r: "<< bufferPosition << std::endl;

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const VectorType bufferPosition = this->bufferPosition;
      const VectorType bufferOrientation = this->bufferOrientation;

      auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         const VectorType r_relative = bufferPosition - r;

         if( ( r_relative, bufferOrientation ) > 0 ){
            printf(" @ bump: r_x: %f, r_x_buffer: %f \n", r[ 0 ], bufferPosition[ 0 ]);
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         return 0;
      };
      this->numberOfPtcsToBuffer = Algorithms::reduce< DeviceType >(
            0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );
      std::cout << "[gFPETB] Number of particles to buffer: " << numberOfPtcsToBuffer << std::endl;

      // sort the indices
      const IndexType rangeToSort = ( numberOfZoneParticles > numberOfBufferParticles ) ? numberOfZoneParticles : numberOfBufferParticles;
      std::cout << "[gFPETB] range to sort: " << rangeToSort << std::endl;
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
      auto r_fluid_view = fluid->getParticles()->getPoints().getView();
      const auto v_fluid_view = fluid->getVariables()->v.getConstView();
      const auto rho_fluid_view = fluid->getVariables()->rho.getConstView();

      auto r_buffer_view = this->getParticles()->getPoints().getView();
      auto v_buffer_view = this->getVariables()->v.getView();
      auto rho_buffer_view = this->getVariables()->rho.getView();
      const auto particlesToBuffer_view = this->particlesToBuffer.getConstView();

      std::cout << "[cFTB] Number of particles to buffer: " << numberOfPtcsToBuffer << std::endl;

      auto retypeFluidToBuffer = [=] __cuda_callable__ ( int i ) mutable
      {
         const IndexType p = particlesToBuffer_view[ i ];

         r_buffer_view[ numberOfBufferPtcs + i ] = r_fluid_view[ p ];
         v_buffer_view[ numberOfBufferPtcs + i ] = v_fluid_view[ p ];
         rho_buffer_view[ numberOfBufferPtcs + i ] = rho_fluid_view[ p ];

         r_fluid_view[ p ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfPtcsToBuffer, retypeFluidToBuffer );

      std::cout << "[cFTB] New number of buffer particles: " << numberOfBufferPtcs + numberOfPtcsToBuffer << std::endl;

      // update the particles counts
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs + numberOfPtcsToBuffer );
      fluid->getParticles()->setNumberOfParticlesToRemove(
            fluid->getParticles()->getNumberOfParticlesToRemove() + numberOfPtcsToBuffer );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   updateInterfaceBuffer( FluidPointer& fluid, FluidPointer& fluid_neihgbor, ModelParams& modelParams, const RealType dt, const int subdomainIdx )
   {
      // update buffer: buffer -> fluid, buffer -> remove
      std::cout << "Move buffer particles: " << std::endl;
      moveBufferParticles( dt, subdomainIdx ); //subdomain idx for debugging
      std::cout << "Sort buffer particles: " << std::endl;
      sortBufferParticles();

      /*
            //debug inlet
            if( ( subdomainIdx == 1 ) && ( numberOfPtcsToRetype > 0 ) )
               std::cout << "Ptcs: " << this->getParticles()->getPoints() << std::endl;
      */

      std::cout << "Remove buffer particles: " << std::endl;
      removeBufferParticles();
      std::cout << "Convert buffer to fluid: " << std::endl;
      convertBufferToFluid( fluid );

      // fluid -> buffer
      std::cout << "Update particles in zone: " << std::endl;
      zone.updateParticlesInZone( fluid->getParticles() );
      std::cout << "Get fluid particles entering the zone: " << std::endl;
      getFluidParticlesEneringTheBuffer( fluid );
      std::cout << "Conver fluid to buffer: " << std::endl;
      convertFluidToBuffer( fluid );

      // mass-nodes -> buffer
      std::cout << "Update mass nodes: " << std::endl;
      updateMassNodes( modelParams, dt );
      std::cout << "Sort mass nodes: " << std::endl;
      massNodes.sort();
      std::cout << "Create buffer particles: " << std::endl;
      createBufferParticles( modelParams );

      std::cout << "BUFFER UPDATED: buffer particles count: " << this->getNumberOfParticles() << " fluid particles count: " << fluid->getNumberOfParticles() << std::endl;

      //// update the values
      massNodes.massFlux = 0;
      std::cout << "Accumulate mass:" << std::endl;
      if( subdomainIdx == 1 )
      accumulateMasses( fluid_neihgbor, modelParams, dt ); // ?
      std::cout << "Interpolate variables: " << std::endl;
      if( subdomainIdx == 1 ) //LARGEST FIXME
      interpolateVariables( fluid_neihgbor, modelParams );

      std::cout << "# MR buffer updated." << std::endl;

      //sun search to rmeove the retyped particles
      if( numberOfPtcsToRetype > 0 )
         this->searchForNeighbors();

      // reset temorary variables and fields (i prefere names such as: fluidToBufferCount )
      numberOfPtcsToRemove = 0;
      numberOfPtcsToRetype = 0; //TODO: Move to fluid
      numberOfPtcsToCreate = 0;
      numberOfPtcsToBuffer = 0;
      //reset the flux


   }

   void
   writeProlog( TNL::Logger& logger, const int subdomainIdx )
   {
      logger.writeParameter( "Subdomain index:", subdomainIdx );
      BaseType::writeProlog( logger );
      logger.writeSeparator();
      logger.writeParameter( "Buffer position:", this->bufferPosition );
      logger.writeParameter( "Buffer orientation:", this->bufferOrientation );
      logger.writeParameter( "Buffer width:", this->bufferWidth );
      //if( subdomainIdx == 0 ){
      //   zone_right.writeProlog( logger );
      //}
      //if( subdomainIdx == 1 ){
      //   zone_left.writeProlog( logger );
      //}
      zone.writeProlog( logger );
   }

   int ownerSetID;
   int neighborSetID;

   MassNodes massNodes;

   // temporary constants
   IndexType numberOfPtcsToRemove = 0;
   IndexType numberOfPtcsToRetype = 0;
   IndexType numberOfPtcsToCreate = 0;
   IndexType numberOfPtcsToBuffer = 0;

   //IndexType fluidToBufferCount;

   // temporary arrays
   IndexArrayType retypeMarker;


   IndexArrayType particlesToFluid;
   IndexArrayType particlesToRemove;
   IndexArrayType particlesToBuffer;

   // buffer referential position (TODO: Use openbc config?)
   VectorType bufferPosition;
   VectorType bufferOrientation; //FIXME: temp
   RealType bufferWidth; // FIXME: temp

   // FIXME: There should be two zones. Or one big unified zone
   // ParticleZone zone_left;
   // ParticleZone zone_right;
   ParticleZone zone;

};

} // SPH
} // TNL

