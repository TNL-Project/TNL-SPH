#pragma once

#include <climits>
#include <string>
#include "ParticleSet.h"
#include <TNL/Particles/GhostZone.h>
#include "OpenBoundaryBuffers.h"
#include "SPHTraits.h"

//FIXME: Include interpolatiion as a single module
#include "shared/Interpolation.h"
#include "shared/WendlandC2ABFs.h"


namespace TNL {
namespace SPH {

template< typename VectorType, typename RealType, typename IndexType >
struct BufferSide
{
   VectorType orientation;    // inward unit normal for this face
   VectorType position;       // a point on the interface plane (face centre)
   RealType width;            // buffer width (physical, outward from position)
   IndexType massNodeOffset;  // start index into the shared mass node array
   IndexType massNodeCount;   // number of mass nodes on this face
};

//TODO: Consider MassNodes to be standard ParticleSet with custom variables
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
                           particlesToCreate.getData() + numberOfMassNodes,
                           thrust::make_zip_iterator( thrust::make_tuple(
                              points.getArrayData(), normal.getArrayData(), massFlux.getArrayData(), mass.getArrayData() ) ) );
   }

   IndexType numberOfMassNodes = 0;
   VectorArrayType points;
   VectorArrayType normal;
   ScalarArrayType massFlux;
   ScalarArrayType mass;
   IndexArrayType particlesToCreate;
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
   using IndexVectorType = typename SolverTraitsType::IndexVectorType;  //TODO: Due to init

   using ParticleZone = TNL::ParticleSystem::ParticleZone< typename ParticlesType::Config, typename ParticlesType::DeviceType >;
   using MassNodes = MassNodes< SPHCaseConfig >;

   using KernelFunction = typename SPHDefs::KernelFunction;
   using MFD = Interpolation::MFD< 2, 1, RealType, Interpolation::WendlandC2ABFs, KernelFunction, SPHCaseConfig >;
   using ABFs = typename MFD::ABFs;
   using MfdVectorType = typename MFD::BaseVectorType;
   using MfdMatrixType = typename MFD::BaseMatrixType;
   using MfdVectorPackType = Containers::StaticArray< VectorType::getSize(), MfdVectorType >; //TODO: Remove for pure 3D

   static constexpr RealType bufferWidthFactorConst = 1.5f;
   static constexpr int frameWidth = 2;

   MultiresolutionBoundary() = default;

   // helper functions
   __cuda_callable__
   bool
   isInsideBox( const VectorType& point; const VectorType& boxOrigin, const VectorType& boxSize )
   {
      bool isInside = true;
      for( int d = 0; d < VectorType::getSize(); d++ )
         if( ( point[ d ] < boxOrigin[ d ] ) || ( point[ d ] >= ( domainOrigin[ d ] + boxSize[ d ] ) ) )
            isInside = false;
      return isInside;
   }

   __cuda_callable__
   bool
   isInsideBox( const IndexVectorType& coords; const IndexArrayType& boxOriginCoords, const IndexVectorType& boxDims )
   {
      bool isInside = true;
      for( int i = 0; i < getParticlesDimension(); i++ )
         if( ( coords[ i ] < 0 ) || ( coords[ i ] >= gridDimensionsWithOverlap[ i ] ) )
            isInside = false;
      return isInside;
   }

   // needs to be called after the object is initialized
   template< typename ParticleSetPointer >
   void
   initZones( const ParticleSetPointer& ownParticles,
              const ParticleSetPointer& nbParticles,
              const int maxNumberOfPtcsPerCell = 75 )
/*
initZonesRectangular( const ParticleSetPointer& ownParticles,
                      const Containers::Array< BufferSideType, Devices::Host >& sides,
                      const VectorType& fineRegionOrigin,
                      const VectorType& fineRegionEnd,
                      const IndexType frameWidth = 2,
                      const int maxPtcsPerCell = 75 )
*/
   {
      sides_  = sides;

      const VectorType globalOrig = ownParticles->getGridReferentialOrigin();
      const VectorType ownOrigin = ownParticles->getGridOrigin();
      const VectorType nbOrigin = nbParticles->getGridOrigin();
      const IndexVectorType ownDims = ownParticles->getGridDimensions();
      const IndexVectorType nbDims = nbParticles->getGridDimensions();
      const IndexVectorType ownDimsWithOverlap = ownParticles->getGridDimensionsWithOverlap();
      const RealType own_sr = ownParticles->getSeachRadius();
      const RealType nb_sr = nbParticles->getSearchRadius();
      const IndexType overlapWidth = ownParticles->getOverlapWidth();
      const VectorType unitVect = 1;

      int resolutionFactor = own_sr / nb_sr;
      inner_overlap = isInsideBox( nbOrig, ownOrig, own_sr * ownDims );
      outer_overlap = isInsideBox( ownOrig, nbOrig, nb_sr * nbDims ); //TODO: Just use negation, right?

      assert( !(inner_overlap && outer_overlap) && "inner_overlap and outer_overlap cannot both be true!" );

      bufferWidth = bufferWidthFactorConst * searchRadius;
      // if zone is inner - compute from neighbors params
      if( inner_overlap ){
         frameOrigin = nbOrig;
         frameFrontOriginCoors = TNL::Floor( ( nbOrig - globalOrig ) / searchRadius );
         frameFrontDims = resolutionFactor * nbDims;
         frameFrontEnd = frameOrigin + frameDims;
         frameOrientation = -1;

         frameBackOrigin = frameFrontOrigin + bufferWidth * unitVect;
         frameBackSize = frameFrontSize - 2 * bufferWidth * unitVect
      }
      // if zone is outer - compute from local params
      else if( outer_overlap ) {
         frameFrontOrigin = ownOrig;
         frameFrontOriginCoords = 0;
         franeFrontDims = ownDims;
         frameFrontEnd = frameOrigin + frameDims;
         frameOrientation = 1;

         frameBackOrigin = frameFrontOrigin - bufferWidth * unitVect;
         frameBackSize = frameFrontSize + 2 * bufferWidth * unitVect
      }
      else {
         assert( false && "initZones: Invalid overlap state: neither inner_overlap nor outer_overlap is true!" );
      }

      zone.setNumberOfParticlesPerCell( maxPtcsPerCell );
      //TODO: Alternatively use the two concentric frameFront and frameBack
      zone.assignCellsFrame( frameFrontOrigin, frameDims, frameOrientation * frameWidth, ownDimesWithOverlap );
      bufferWidth = bufferWidthFactorConst * searchRadius;

      // Set size of multi-resoluton algorithm arrays (NOTE: Consider setSize function)
      const IndexType n_alloc = this->getNumberOfAllocatedParticles();
      particlesToFluid.setSize( n_alloc );
      particlesToRemove.setSize( n_alloc );
      particlesToBuffer.setSize( n_alloc );
      retypeMarker.setSize( n_alloc );
   }

   template< typename ModelParams >
   void
   initMassNodes( ModelParams& modelParams, const int subdomainIdx, const RealType refinemnetFactor )
   {
      const RealType searchRadius = this->getParticles()->getSearchRadius();
      const VectorType subdomainSize = searchRadius * this->getParticles()->getGridDimensions();
      const IndexType overlapWidth = this->getParticles()->getOverlapWidth();
      const RealType local_dp = refinemnetFactor * modelParams.dp;

      // Get the range to create mass nodes
      IndexVectorType begin = 0;
      IndexVectorType end = 0;
      IndexType n_massNodes = 1;

      for( int d = 0; d < VectorType::getSize(); d++ ){
         if( interfaceAxis[ d ] == 0 ){
            begin[ d ] = 1; //FIXME 0 is right, 1 is for debug
            end[ d ] = subdomainSize[ d ] / local_dp;
            n_massNodes *= end[ d ];
         }
         else{
            begin[ d ]= 0;
            end[ d ] = 1;
         }
      }
      massNodes.setSize( n_massNodes );

      // Initialize tha mass points coordinates
      const VectorType nodeNormal = bufferOrientation;
      const IndexVectorType unitVect = 1;
      const IndexVectorType perpAxis = unitVect - interfaceAxis;
      const VectorType nodeStartingPoint = bufferPosition + ( -1 ) * overlapWidth * searchRadius * bufferOrientation; //TODO: In for above?

      // Strides for linearisation over perpendicular axes only
      /*
         Strides for linearisation over perpendicular axes only.
         Interface axis stride stays 0 — idx[interfaceAxis] is always 0.
         In 2D with x as interface axis: stride = {0, 1} gives i = idx[1]
         In 2D with y as interface axis: stride = {1, 0} gives i = idx[0]
         In 3D with x as interface axis: stride = {0, ny, 1} gives i = idx[1]*nz + idx[2]
         General: stride[d] = product of end[k] for all perp k > d
       */
      IndexVectorType stride = 0;
      {
         IndexType running = 1;
         // Walk axes in reverse to build strides for perp axes only
         for( int d = VectorType::getSize() - 1; d >= 0; d-- ) {
            if( perpAxis[ d ] == 1 ) {
               stride[ d ] = running;
               running *= end[ d ];
            }
         }
      }

      auto points_nodes_view = this->massNodes.points.getView();
      auto normals_nodes_view = this->massNodes.normal.getView();
      auto generateMassNodesCoordinates = [ = ] __cuda_callable__( const IndexVectorType idx ) mutable
      {
         // Linearise: dot product of idx and stride (interface axis contributes 0)
         IndexType i = ( idx, stride );;
         const VectorType r = nodeStartingPoint + local_dp * ( idx + 1 ) * perpAxis;
         points_nodes_view[ i ] = r;
         normals_nodes_view[ i ] = nodeNormal;
      };
      Algorithms::parallelFor< DeviceType >( begin, end, generateMassNodesCoordinates );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   accumulateMasses( FluidPointer& fluid_neihgbor, ModelParams& modelParams, const RealType dt )
   {
      auto searchInFluid = fluid_neihgbor->getParticles()->getSearchToken( fluid_neihgbor->getParticles() );

      const auto view_points_massNodes = this->massNodes.points.getConstView();
      const auto view_normals_massNodes = this->massNodes.normal.getConstView();
      const auto view_m_massPoints = massNodes.mass.getView();  //TODO: temp due to condition
      auto view_mFlux_massPoints = massNodes.massFlux.getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const unsigned int dim = SPHCaseConfig::spaceDimension;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius();  //TODO: Is it this one?
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      const RealType h = refinementFactor * modelParams.h;
      const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
      const RealType dx = refinementFactor * modelParams.dp * 0.5f;  //FIXME I have no clue why there is 0.5 factor
      const VectorType v_subdomain = 0.f;  // velocity of moving subdomain
      const RealType div_r_trashold = 1.5f;  //FIXME add to model params, depends on dimension
      const RealType extrapolationDetTreshold = modelParams.mdbcExtrapolationDetTreshold;  //TODO: rename, remove mdbc

      auto interpolateFluid = [ = ] __cuda_callable__( IndexType i,
                                                       IndexType j,
                                                       VectorType & r_x,
                                                       MfdMatrixType * M_x,
                                                       MfdVectorType * brho_x,
                                                       MfdVectorPackType * bv_x,
                                                       RealType * div_r_x,
                                                       RealType * drs_min ) mutable
      {
         const VectorType r_j = view_points_fluid[ j ];
         const VectorType r_xj = r_x - r_j;
         const RealType drs = l2Norm( r_xj );
         if( drs <= searchRadius ) {
            const RealType rho_j = view_rho_fluid[ j ];
            const VectorType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const MfdVectorType b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            for( int d = 0; d < VectorType::getSize(); d++ )
               ( *bv_x )[ d ] += v_j[ d ] * b_x;

            // get div r
            const VectorType gradW = r_xj * KernelFunction::F( drs, h );
            *div_r_x += ( -1.f ) * ( r_xj, gradW ) * V_j;  //TODO: Added - sign, should be there?

            // find nearest neighbor
            *drs_min = std::min( *drs_min, drs );
         }
      };

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const VectorType r_x = view_points_massNodes[ i ];
         const VectorType normal_x = view_normals_massNodes[ i ];
         MfdMatrixType M_x = 0.f;
         MfdVectorType brho_x = 0.f;
         MfdVectorPackType bv_x;
         for( int d = 0; d < VectorType::getSize(); d++ )
            bv_x[ d ] = 0.f;
         RealType div_r_x = 0.f;
         RealType drs_min = FLT_MAX;

         ParticlesType::NeighborsLoop::exec(
               i,
               r_x,
               searchInFluid,
               interpolateFluid,
               &M_x,
               &brho_x,
               &bv_x,
               &div_r_x,
               &drs_min );

         RealType rho_x;
         VectorType v_x;

         //TODO: Use LU Decomposition so we can just reuse it different RHS
         //TODO: Proper condition should be if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold )
         if( M_x( 0, 0 ) > 0.05 ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = Matrices::solve( M_x, bv_x[ d ] )[ 0 ];
         }
         else if( M_x( 0, 0 ) > 0.f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = bv_x[ d ][ 0 ] / M_x( 0, 0 );
         }
         else {
            // TODO: not sure what to do here
         }

         // Ricci et al. uses (-1) * rho, but I think I have different normal orientation
         const RealType m_x_lessZero = TNL::max( 0, rho_x * ( v_x - v_subdomain, normal_x ) * dx * dt );
         const RealType m_x_geqZero = rho_x * ( v_x - v_subdomain, normal_x ) * dx * dt;

         // sum up the flux with free surface corrections
         RealType m_flux_x = 0;

         //TODO: Should I use ( view_mFlux_massPoints[ i ] < 0.f  ) here?
         if( ( div_r_x > div_r_trashold || drs_min < dx ) && ( view_m_massPoints[ i ] < 0.f ) )
            m_flux_x = m_x_lessZero;
         else if( div_r_x > div_r_trashold || drs_min < dx )
            m_flux_x = m_x_geqZero;
         else
            m_flux_x = 0;

         view_mFlux_massPoints[ i ] += m_flux_x;
      };
      Algorithms::parallelFor< DeviceType >( 0, massNodes.numberOfMassNodes, particleLoop );
   }

   //interpolate values to buffer particles
   template< typename FluidPointer, typename ModelParams >
   void
   interpolateVariables( FluidPointer& fluid_neihgbor, ModelParams& modelParams )
   {
      auto searchInFluid = this->getParticles()->getSearchToken( fluid_neihgbor->getParticles() );
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();

      auto view_points_overlap = this->getParticles()->getPoints().getView();
      auto view_rho_overlap = this->getVariables()->rho.getView();
      auto view_v_overlap = this->getVariables()->v.getView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const unsigned int dim = SPHCaseConfig::spaceDimension;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius();
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      const RealType h = refinementFactor * modelParams.h;
      const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
      //TODO: const RealType trasholdM00 = ..
      //TODO: const RealType extrapolationDetTreshold = modelParams.extrapolationDetTrashold;

      auto interpolateFluid = [ = ] __cuda_callable__( IndexType i,
                                                       IndexType j,
                                                       VectorType & r_x,
                                                       MfdMatrixType * M_x,
                                                       MfdVectorType * brho_x,
                                                       MfdVectorPackType * bv_x ) mutable
      {
         const VectorType r_j = view_points_fluid[ j ];
         const VectorType r_xj = r_x - r_j;
         const RealType drs = l2Norm( r_xj );
         if( drs <= searchRadius ) {
            const RealType rho_j = view_rho_fluid[ j ];
            const VectorType v_j = view_v_fluid[ j ];
            const RealType V_j = m / rho_j;

            *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
            const MfdVectorType b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
            *brho_x += rho_j * b_x;
            for( int d = 0; d < VectorType::getSize(); d++ )
               ( *bv_x )[ d ] += v_j[ d ] * b_x;
         }
      };

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const VectorType r_x = view_points_overlap[ i ];
         MfdMatrixType M_x = 0.f;
         MfdVectorType brho_x = 0.f;
         MfdVectorPackType bv_x;
         for( int d = 0; d < VectorType::getSize(); d++ )
            bv_x[ d ] = 0.f;

         ParticlesType::NeighborsLoop::exec( i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bv_x );

         RealType rho_x;
         VectorType v_x;

         //TODO: Use LU Decomposition so we can just reuse it different RHS
         //TODO: Proper condition should be if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold )
         const RealType detM = Matrices::determinant( M_x );
         if( M_x( 0, 0 ) > 0.05 && detM > 0.001f ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = Matrices::solve( M_x, bv_x[ d ] )[ 0 ];

            view_rho_overlap[ i ] = rho_x;
            view_v_overlap[ i ] = v_x;
            return 0;
         }
         else if( M_x( 0, 0 ) > 0.05f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = bv_x[ d ][ 0 ] / M_x( 0, 0 );

            view_rho_overlap[ i ] = rho_x;
            view_v_overlap[ i ] = v_x;
            return 0;
         }
         else {
            // FIXME: not sure what to do here, right now, I remove the particle
            view_points_overlap[ i ] = FLT_MAX;
            return 1;
         }

      };
      const IndexType numberOfInvalidBufferParticles =
         Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, particleLoop );
      this->getParticles()->setNumberOfParticlesToRemove( this->getParticles()->getNumberOfParticlesToRemove()
                                                          + numberOfInvalidBufferParticles );
   }

   /*
      template< typename FluidPointer, typename ModelParams >
      void
      shiftParticles( FluidPointer& fluid, ModelParams& modelParams )
      {

         auto interpolateFluid = [=] __cuda_callable__ (
               IndexType i,
               IndexType j,
               VectorType& r_i
               VectorType* gradC_i,
               RealType* gamma_i ) mutable
         {
            const VectorType r_j = view_points_fluid[ j ];
            const VectorType r_xj = r_x - r_j;
            const RealType drs = l2Norm( r_xj );
            if( drs <= searchRadius )
            {
               const RealType WV_j = KernelFunction::W( drs, h ) * m / rho_j;
               const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;
               *gradC_i += gradWV_j;
               *gamma_i += WV_j;
            }
         };


         auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
         {
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfParticles )
      }
   */

   void
   moveBufferParticles( const RealType dt )
   {
      auto view_r_buffer = this->getParticles()->getPoints().getView();
      const auto view_v_buffer = this->getVariables()->v.getConstView();
      const IndexType n_buffer = this->getParticles()->getNumberOfParticles();

      // move buffer particles
      auto moveBufferParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles );

      // reset retype marker
      auto retypeMarker_view = retypeMarker.getView();
      retypeMarker_view = 0;

      // Precompute conversion factors — captured by value into kernels.
      const RealType sr = this->getParticles()->getSearchRadius();
      const RealType inv_sr = 1.f / sr;
      //const VectorType frameFrontOrigin = frameFrontOrigin;
      //const IndexVectorType frameFronOriginCoords = frameFronOriginCoords;
      //const IndexVectorType frameFrontDims = frameFrontDimes;
      //const VectorType frameBackOrigin = frameBackOrigin;
      //const VectorType frameBackSize = frameBackSize;
      //bool inner_overlap = inner_overlap;
      //bool outer_overlap = outer_overlap;

      // Helper lambda: convert physical position to own grid integer coords.
      auto toGridCoord = [=] __cuda_callable__ ( const VectorType& r ) -> IndexVectorType
      {
         IndexVectorType gc;
         gc = TNL::Floor( ( r - frameFrontOrigin ) * inv_sr ) //TODO: Static cast?
         return gc;
      };

      // Retype to fluid if:
      // - outer zone + is inside
      // - inner zone + is outside
      auto identifyRetype = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const IndexVectorType gc = toGridCoord( view_r[ i ] );
         bool inside = isInsideBox( gc, frameFrontOriginCoords, frameFrontDims );

         if( inside && outer_overlap) {
            retypeMarker_view[ i ] = 1;
            return 1;
         }
         else if( !inside && inner_overlap  ) {
            retypeMarker_view[ i ] = 1;
            return 1;
         }
         return 0;
      };
      numberOfPtcsToRetype = Algorithms::reduce< DeviceType >( 0, n_buffer, identifyRetype );

      // Remove invalid buffe particles
      // - outer zone - is outside the ( subdomain box + zone width )
      // - inner zone - is insde the ( subdomain box - zone width )
      // TODO (?): If we would use buffer width corresponding to search radius, we can use grid to compare positions
      auto identifyRemove = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         bool inside = isInsideBox( view_r[ i ], frameBackOrigin, frameBackSize );
         if( inside && inner_overlap ) {
            retypeMarker_view[ i ] = 2;
            return 1;
         }
         else if( !inside && outer_overlap ) {
            retypeMarker_view[ i ] = 2;
            return 1;
         }
         return 0;
      };
      numberOfPtcsToRemove = Algorithms::reduce< DeviceType >( 0, n, identifyRemove );
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
                           thrust::make_zip_iterator(
                              thrust::make_tuple( r_view.getArrayData(), v_view.getArrayData(), rho_view.getArrayData() ) ) );
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

      const IndexType numberOfBufferParticles = this->getNumberOfParticles();
      const IndexType toFluidOffset = numberOfBufferParticles - this->numberOfPtcsToRetype;

      auto createNewFluidParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         view_r_fluid[ numberOfFluidPtcs + i ] = view_r_buffer[ toFluidOffset + i ];
         view_rho_fluid[ numberOfFluidPtcs + i ] = view_rho_buffer[ toFluidOffset + i ];
         view_v_fluid[ numberOfFluidPtcs + i ] = view_v_buffer[ toFluidOffset + i ];

         view_rho_old[ numberOfFluidPtcs + i ] = view_rho_buffer[ toFluidOffset + i ];
         view_v_old[ numberOfFluidPtcs + i ] = view_v_buffer[ toFluidOffset + i ];

         //TODO: I need to remove the used particles, not sure if this is the right way
         view_r_buffer[ toFluidOffset + i ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfPtcsToRetype, createNewFluidParticles );
      fluid->getParticles()->setNumberOfParticles( numberOfFluidPtcs + this->numberOfPtcsToRetype );
      this->getParticles()->setNumberOfParticles( numberOfBufferParticles - this->numberOfPtcsToRetype );
   }

   void
   removeBufferParticles()
   {
      const IndexType numberOfBufferPtcs = this->getNumberOfParticles();
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs - this->numberOfPtcsToRemove );
   }

   template< typename ModelParams >
   void
   updateMassNodes( ModelParams& modelParams, const RealType dt )
   {
      const auto massFlux_view = massNodes.massFlux.getConstView();
      auto mass_view = massNodes.mass.getView();
      auto view_particlesToCreate = massNodes.particlesToCreate.getView();

      //FIXME: I need refinemet factor here! Subdomain mass is different. But with refinement factor, it doesnt work!
      //const RealType refinementFactor = this->getParticles()->getSearchRadius() / ( 2.f * modelParams.h );
      //const RealType particleMass = std::pow( refinementFactor, ParticlesType::spaceDimension ) * modelParams.mass;
      const RealType particleMass = 0.25f * modelParams.mass;
      const IndexType numberOfMassNodes = this->massNodes.numberOfMassNodes;

      // reset list with markers
      view_particlesToCreate = 0;

      auto identifyParticlesToCreate = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         mass_view[ i ] += massFlux_view[ i ];  //NOTE: dt is already included!
         if( mass_view[ i ] > particleMass ) {
            mass_view[ i ] -= particleMass;
            view_particlesToCreate[ i ] = -1;
            return 1;
         }
         else
            return 0;
      };
      this->numberOfPtcsToCreate = Algorithms::reduce< DeviceType >( 0, numberOfMassNodes, identifyParticlesToCreate );
   }

   template< typename ModelParams >
   void
   createBufferParticles( ModelParams& modelParams )
   {
      auto r_buffer_view = this->getParticles()->getPoints().getView();
      const auto r_massNodes_view = this->massNodes.points.getConstView();
      const auto normal_massNodes_view = this->massNodes.normal.getConstView();

      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();
      const RealType refinementFactor = this->getParticles()->getSearchRadius() / ( 2.f * modelParams.h );
      const RealType dp = refinementFactor * modelParams.dp;

      auto createNewBufferParticles = [ = ] __cuda_callable__( int i ) mutable
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

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto zoneParticleIndices_view = this->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = this->zone.getNumberOfParticles();

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const VectorType bufferPosition = this->bufferPosition;
      const VectorType bufferOrientation = this->bufferOrientation;

      auto checkFluidParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         const VectorType r_relative = bufferPosition - r;

         if( ( r_relative, bufferOrientation ) > 0 ) {
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         return 0;
      };
      this->numberOfPtcsToBuffer =
         Algorithms::reduce< DeviceType >( 0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );

      // sort the indices
      const IndexType rangeToSort =
         ( numberOfZoneParticles > numberOfBufferParticles ) ? numberOfZoneParticles : numberOfBufferParticles;
      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::sort( thrustDevice, particlesToBuffer_view.getArrayData(), particlesToBuffer_view.getArrayData() + rangeToSort );
   }

   template< typename FluidPointer >
   void
   convertFluidToBuffer( FluidPointer& fluid )
   {
      const IndexType numberOfBufferPtcs = this->getParticles()->getNumberOfParticles();

      //TODO: I'm sure that it is enough to copy the positioin.
      auto r_fluid_view = fluid->getParticles()->getPoints().getView();
      const auto v_fluid_view = fluid->getVariables()->v.getConstView();
      const auto rho_fluid_view = fluid->getVariables()->rho.getConstView();

      auto r_buffer_view = this->getParticles()->getPoints().getView();
      auto v_buffer_view = this->getVariables()->v.getView();
      auto rho_buffer_view = this->getVariables()->rho.getView();
      const auto particlesToBuffer_view = this->particlesToBuffer.getConstView();

      auto retypeFluidToBuffer = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType p = particlesToBuffer_view[ i ];

         r_buffer_view[ numberOfBufferPtcs + i ] = r_fluid_view[ p ];
         v_buffer_view[ numberOfBufferPtcs + i ] = v_fluid_view[ p ];
         rho_buffer_view[ numberOfBufferPtcs + i ] = rho_fluid_view[ p ];
         r_fluid_view[ p ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfPtcsToBuffer, retypeFluidToBuffer );
      this->getParticles()->setNumberOfParticles( numberOfBufferPtcs + numberOfPtcsToBuffer );
      fluid->getParticles()->setNumberOfParticlesToRemove( fluid->getParticles()->getNumberOfParticlesToRemove()
                                                           + numberOfPtcsToBuffer );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   updateInterfaceBuffer( FluidPointer& fluid_own,
                          FluidPointer& fluid_neihgbor,
                          const ModelParams& modelParams,
                          const RealType dt,
                          const int subdomainIdx )
   {
      massNodes.massFlux = 0;

      accumulateMasses( fluid_neihgbor, modelParams, dt );
      moveBufferParticles( dt );
      sortBufferParticles();
      removeBufferParticles();
      convertBufferToFluid( fluid_own );
      zone.updateParticlesInZone( fluid_own->getParticles() );
      getFluidParticlesEneringTheBuffer( fluid_own );
      convertFluidToBuffer( fluid_own );
      updateMassNodes( modelParams, dt );
      massNodes.sort();
      createBufferParticles( modelParams );
      interpolateVariables( fluid_neihgbor, modelParams );
      this->searchForNeighbors();

      numberOfPtcsToRemove = 0;
      numberOfPtcsToRetype = 0;
      numberOfPtcsToCreate = 0;
      numberOfPtcsToBuffer = 0;
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
      zone.writeProlog( logger );
   }

protected:

   MassNodes massNodes;

   // temporary constants //TODO: Use names such as fluidToBufferCount
   IndexType numberOfPtcsToRemove = 0;
   IndexType numberOfPtcsToRetype = 0; //TODO: Rename to "moveToFluid"
   IndexType numberOfPtcsToCreate = 0;
   IndexType numberOfPtcsToBuffer = 0;

   // temporary arrays
   IndexArrayType retypeMarker; //TODO: rename

   IndexArrayType particlesToFluid;
   IndexArrayType particlesToRemove;
   IndexArrayType particlesToBuffer;

   // buffer referential specification (TODO: Consider to use some buffer config class)
   IndexVectorType interfaceAxis;
   VectorType bufferPosition;
   VectorType bufferOrientation;
   RealType bufferWidth;

   // UPDATE FOR RECTANGULAR ZONES //TODO: Consider to use topology instead of duplicitly stored adata
   IndexVectorType subdomainOriginCoords;
   IndexVectorType subdomainEndCoords;

   //renamed
   IndexVectorType frameFronOriginCoords;
   IndexVectorType frameFronDims;
   IndexVectorType frameFronEnd;

   VectorType frameBackOrigin;
   VectorType frameBackSize;

   int frameOrientation;

   bool inner_overlap;
   bool outer_overlap;

   //-----

   ParticleZone zone;
};

}  //namespace SPH
}  //namespace TNL

