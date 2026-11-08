// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <climits>
#include <functional>
#include <string>

#include <cuda/std/tuple>

#include <TNL/Particles/GhostZone.h>

#include <SPH/ParticleSet.h>
#include <SPH/MultiresolutionBufferHelpers.h>
#include <SPH/OpenBoundaryBuffers.h>
#include <SPH/SPHTraits.h>
//FIXME: Include interpolatiion as a single module
#include <SPH/shared/Interpolation.h>
#include <SPH/shared/WendlandC2ABFs.h>

namespace TNL::SPH {

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
      thrust::sort_by_key(
         thrustDevice,
         particlesToCreate.getData(),
         particlesToCreate.getData() + numberOfMassNodes,
         thrust::make_zip_iterator(
            cuda::std::make_tuple(
               points.getArrayData(), normal.getArrayData(), massFlux.getArrayData(), mass.getArrayData() ) ) );
   }

   IndexType numberOfMassNodes = 0;
   VectorArrayType points;
   VectorArrayType normal;
   ScalarArrayType massFlux;
   ScalarArrayType mass;
   IndexArrayType particlesToCreate;
};

template<
   typename ParticlesType,
   typename SPHCaseConfig,
   typename Variables,
   typename IntegratorVariables,
   typename OpenBoundaryConfig,
   typename SPHDefs >
// TODO: Derive from ParticleSet or OpenBoundaryBuffer?
/**
 * \brief Overlap-aware multiresolution buffer patch for uniform and local time stepping.
 *
 * The band width follows the patch's own ParticleSet overlap width: an intermediate
 * frame frameMid splits the band into layer 1 (frameMid XOR frameFront) and layer 2
 * (frameBack XOR frameMid - passive BC). With overlap width 1, layer 1 is empty and the
 * patch reduces to the single-timestep behavior. The LTS update is staged:
 * accumulateMasses() per substep, advectBufferParticles() per substep (both layers
 * drift ballistically), refreshVariablesFromNeighbor() to re-ground the band state
 * from the coarser neighbor at level-clock coincidences/midpoints, and
 * syncUpdateInterfaceBuffer() at sync points only (no advection, no accumulation,
 * no membership change elsewhere in the window).
 */
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
   using CoordinatesType = typename SolverTraitsType::CoordinatesType;

   using ParticleZone = TNL::Particles::ParticleZone< typename ParticlesType::Config, typename ParticlesType::DeviceType >;
   using MassNodes = MassNodes< SPHCaseConfig >;

   using KernelFunction = typename SPHDefs::KernelFunction;
   using MFD = Interpolation::
      MFD< ParticlesType::getParticlesDimension(), 1, RealType, Interpolation::WendlandC2ABFs, KernelFunction, SPHCaseConfig >;
   using ABFs = typename MFD::ABFs;
   using MfdVectorType = typename MFD::BaseVectorType;
   using MfdMatrixType = typename MFD::BaseMatrixType;
   using MfdVectorPackType = Containers::StaticArray< VectorType::getSize(), MfdVectorType >;  //TODO: Remove for pure 3D

   static constexpr RealType bufferWidthFactorConst = 1;  //FIXME IT CAN NOT BE LARGER THEN THE CELL WIDTH (I would like to have
                                                          //like this, but then, removal of particles needs to be changed)
   static constexpr int frameWidth = 2;

   MultiresolutionBoundary() = default;

   // helper functions
   // TODO: Used only for comparing origins, it can be removed
   __cuda_callable__
   static bool
   isInsideBox( const VectorType& point, const VectorType& boxOrigin, const VectorType& boxSize )
   {
      bool isInside = true;
      for( int d = 0; d < VectorType::getSize(); d++ )
         if( ( point[ d ] < boxOrigin[ d ] ) || ( point[ d ] > ( boxOrigin[ d ] + boxSize[ d ] ) ) )
            isInside = false;
      return isInside;
   }

   __cuda_callable__
   static bool
   isInsideBox( const CoordinatesType& coords, const CoordinatesType& boxDims )
   {
      bool isInside = true;
      for( int i = 0; i < ParticlesType::getParticlesDimension(); i++ )
         if( ( coords[ i ] < 0 ) || ( coords[ i ] >= boxDims[ i ] ) )
            isInside = false;
      return isInside;
   }

   // needs to be called after the object is initialized
   template< typename ParticleSetPointer >
   void
   initZones(
      const ParticleSetPointer& ownParticles,
      const ParticleSetPointer& nbParticles,
      const RealType refinementFraction,
      const int maxParticlesPerCell = 175 )
   {
      //const VectorType globalOrig = ownParticles->getReferentialOrigin();
      const VectorType ownOrig = ownParticles->getOrigin();
      const VectorType nbOrig = nbParticles->getOrigin();
      const CoordinatesType ownDims = ownParticles->getDimensions();
      const CoordinatesType nbDims = nbParticles->getDimensions();
      const CoordinatesType ownDimsWithOverlap = ownParticles->getDimensionsWithOverlap();
      const RealType own_sr = ownParticles->getSearchRadius();
      const RealType nb_sr = nbParticles->getSearchRadius();
      const IndexType overlapWidth = ownParticles->getOverlapWidth();
      //const VectorType unitVect = 1;

      const RealType resolutionFactor = nb_sr / own_sr;
      inner_overlap = isInsideBox( nbOrig, ownOrig, own_sr * ownDims );
      outer_overlap = isInsideBox( ownOrig, nbOrig, nb_sr * nbDims );  //TODO: Just use negation, right?

      // store the interface geometry for the init report (writeInitReport reads back the loaded values)
      refinementLevel = int( 1 / refinementFraction );
      this->resolutionFactor = resolutionFactor;
      ownParticleOrigin = ownOrig;
      nbParticleOrigin = nbOrig;
      ownGlobalOriginCoordinates = ownParticles->getGlobalOriginCoordinates();
      nbGlobalOriginCoordinates = nbParticles->getGlobalOriginCoordinates();

      assert( ! ( inner_overlap && outer_overlap ) && "inner_overlap and outer_overlap cannot both be true!" );

      bufferWidth = this->getParticles()->getOverlapWidth() * own_sr;
      // if zone is inner - compute from neighbors params
      if( inner_overlap ) {
         frameOriginGlobalCoordinates = nbParticles->getGlobalOriginCoordinates() * resolutionFactor;
         frameDimensions = resolutionFactor * nbDims;
         frameOrientation = -1;
      }
      // if zone is outer - compute from local params
      else if( outer_overlap ) {
         frameOriginGlobalCoordinates = ownParticles->getGlobalOriginCoordinates();
         frameDimensions = ownDims;
         frameOrientation = 1;
      }
      else {
         assert( false && "initZones: Invalid overlap state: neither inner_overlap nor outer_overlap is true!" );
      }

      zone.setNumberOfParticlesPerCell( maxParticlesPerCell );

      CoordinatesType zoneOrigin = 0;
      //TODO: I dont like the "overOffset"
      //NOTE: We increase the zone by one additional layers
      if( inner_overlap ) {
         // with-overlap shift (own fluid overlap width) + one-cell inset into the frame,
         // so the capture ring covers one cell inside the frame on every face for any
         // overlap width (base's literal +2 is the special case for overlap width 1)
         zoneOrigin = getFrameFrontOriginGlobalCoordinates() - ownParticles->getGlobalOriginCoordinates()
                    + ownParticles->getOverlapWidth() + 1;
         zoneDimensions = getFrameFrontDimensions() - 2;
      }
      if( outer_overlap ) {
         // zone box in with-overlap cell coords: shift by the own-fluid overlap width
         zoneOrigin = 2 + ownParticles->getOverlapWidth();
         zoneDimensions = getFrameFrontDimensions() - 4;
      }
      zone.assignCellsFrame( zoneOrigin, zoneDimensions, frameWidth + 1, ownDimsWithOverlap );
      //---------------
      const std::string outputFilename = "results/zone_L" + std::to_string( int( 1 / refinementFraction ) ) + ".vtk";
      zone.saveZoneToVTK(
         outputFilename,
         ownDimsWithOverlap,
         ownParticles->getOriginWithOverlap(),  //ownOrigin
         ownParticles->getSearchRadius() );

      // Set size of multi-resoluton algorithm arrays (NOTE: Consider setSize function)
      const IndexType n_alloc = this->getNumberOfAllocatedParticles();
      particlesToFluid.setSize( n_alloc );
      particlesToRemove.setSize( n_alloc );
      particlesToBuffer.setSize( zone.getNumberOfCells() * zone.getNumberOfParticlesPerCell() );
      retypeMarker.setSize( n_alloc );
   }

   /*
   template< typename ModelParams >
   void
   initMassNodes( ModelParams& modelParams, const int subdomainIdx, const RealType refinemnetFactor )
   {
      const RealType local_dp = refinemnetFactor * modelParams.dp;

      auto planeNodeCount = [&]( int faceAxis, int perpAxis ) -> IndexType {
         RealType extent = frameBackSize[ perpAxis ];
         // Each lower-priority face axis that is also perpendicular to perpAxis claims local_dp on each end
         for( int d = 0; d < faceAxis; d++ )
            if( d != perpAxis )
               extent -= 2.f * local_dp; //TODO: Why is there 2?
         return static_cast< IndexType >( TNL::max( 0.f, extent ) / local_dp );
      };

      auto perpOriginOffset = [&]( int faceAxis, int perpAxis ) -> RealType {
         RealType offset = 0.f;
         for( int d = 0; d < faceAxis; d++ )
            if( d != perpAxis )
               offset += local_dp;
         return offset;
      };

      // Count total nodes across all 2*dim faces
      IndexType n_massNodes = 0;
      for( int d = 0; d < VectorType::getSize(); d++ ) {
         IndexType faceNodes = 1;
         for( int pd = 0; pd < VectorType::getSize(); pd++ )
            if( pd != d )
               faceNodes *= planeNodeCount( d, pd );
         n_massNodes += 2 * faceNodes;   // min and max face
      }
      massNodes.setSize( n_massNodes );
      auto points_view = massNodes.points.getView();
      auto normals_view = massNodes.normal.getView();

      // Generate nodes face by face
      IndexType offset = 0;
      for( int d = 0; d < VectorType::getSize(); d++ ) {
         for( int sign : { -1, +1 } ) {

            // Inward normal
            VectorType normal = 0.f;
            normal[ d ] = ( sign < 0 ) ? 1.f : -1.f;
            if( inner_overlap )
               normal[ d ] *= -1.f;

            // Physical coordinate on the interface axis
            const RealType ifaceCoord = ( sign < 0 )
                  ? frameBackOrigin[ d ]
                  : frameBackOrigin[ d ] + frameBackSize[ d ];

            // parallelFor range: [0,1) on d, [0, count) on perp axes
            CoordinatesType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  end[ pd ] = planeNodeCount( d, pd );

            // Strides for linearisation (perp axes only, reverse order)
            CoordinatesType stride = 0;
            {
               IndexType running = 1;
               for( int pd = VectorType::getSize() - 1; pd >= 0; pd-- ) {
                  if( pd == d ) continue;
                  stride[ pd ] = running;
                  running *= end[ pd ];
               }
            }

            // Physical start of the node grid on each perp axis
            //VectorType perpStart = frameFrontOrigin;
            VectorType perpStart = frameBackOrigin;
            perpStart[ d ] = ifaceCoord;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  perpStart[ pd ] += perpOriginOffset( d, pd );

            //TODO: Can I just the variables present in the scope or do I need this?
            const IndexType faceOffset = offset;
            const int iAxis = d;
            const VectorType norm = normal;
            const VectorType pStart = perpStart;

            auto generate = [=] __cuda_callable__ ( const CoordinatesType idx ) mutable
            {
               // Linearise: sum over perp axes only
               IndexType i = faceOffset;
               for( int pd = 0; pd < VectorType::getSize(); pd++ )
                  i += idx[ pd ] * stride[ pd ];

               // Physical position
               VectorType r = pStart;
               for( int pd = 0; pd < VectorType::getSize(); pd++ )
                  if( pd != iAxis )
                     r[ pd ] += local_dp * ( idx[ pd ] + 1 );

               points_view[ i ] = r;
               normals_view[ i ] = norm;
            };
            Algorithms::parallelFor< DeviceType >( begin, end, generate );

            // Advance offset by number of nodes on this face
            IndexType faceNodes = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d ) faceNodes *= end[ pd ];
            offset += faceNodes;
         }
      }

      const std::string outputFileName = "results/massNodes_" + std::to_string( local_dp ) + ".vtk";
      writeMassNodesToVTK( outputFileName );

   }
   */

   template< typename ModelParams >
   void
   initMassNodes( ModelParams& modelParams, const int subdomainIdx, const RealType refinemnetFactor )
   {
      initMassNodes( modelParams, subdomainIdx, -1, refinemnetFactor );
   }

   template< typename ModelParams >
   void
   initMassNodes( ModelParams& modelParams, const int subdomainIdx, const int neighborIdx, const RealType refinemnetFactor )
   {
      //TODO: Make it general
      if constexpr( ParticlesType::getParticlesDimension() == 2 ) {
         TNL::Containers::Array< VectorType, TNL::Devices::Host > excluded( 2 );
         if( outer_overlap ) {
            excluded[ 0 ] = { -1.f, 0.f };  // exclude +x face
            excluded[ 1 ] = { 0.f, 1.f };   // exclude -y face
         }
         if( inner_overlap ) {
            excluded[ 0 ] = { 1.f, 0.f };   // exclude +x face
            excluded[ 1 ] = { 0.f, -1.f };  // exclude -y face
         }
         initMassNodesWithExcludedNormals( modelParams, subdomainIdx, neighborIdx, refinemnetFactor, excluded );
      }

      if constexpr( ParticlesType::getParticlesDimension() == 3 ) {
         TNL::Containers::Array< VectorType, TNL::Devices::Host > excluded( 1 );
         if( outer_overlap ) {
            excluded[ 0 ] = { 0.f, 0.f, 1.f };  // exclude +x face
         }
         if( inner_overlap ) {
            excluded[ 0 ] = { 0.f, 0.f, -1.f };  // exclude +x face
         }
         initMassNodesWithExcludedNormals( modelParams, subdomainIdx, neighborIdx, refinemnetFactor, excluded );
      }
   }

   template< typename ModelParams, typename SubdomainDescriptorT >
   void
   initMassNodesWithDomainExclusion(
      ModelParams& modelParams,
      const int subdomainIdx,
      const int neighborIdx,
      const RealType refinementFactor,
      const SubdomainDescriptorT& ownSd,
      const SubdomainDescriptorT& nbSd,
      const RealType globalCell,
      const VectorType& domMin,
      const VectorType& domMax )
   {
      const auto ctx = computeMassNodeFaceContext( ownSd, nbSd, globalCell );
      auto excludedFaces = computeExcludedDomainFaceNormals< SPHCaseConfig::spaceDimension >(
         ctx.faceMin, ctx.faceMax, domMin, domMax, ctx.tol, ctx.flipNormals );
      initMassNodesWithExcludedNormals( modelParams, subdomainIdx, neighborIdx, refinementFactor, excludedFaces );
   }

   template< typename ModelParams >
   void
   initMassNodesWithExcludedNormals(
      ModelParams& modelParams,
      const int subdomainIdx,
      const int neighborIdx,
      const RealType refinementFactor,
      const TNL::Containers::Array< VectorType, TNL::Devices::Host >& excludedNormals = {} )
   {
      const RealType local_dp = refinementFactor * modelParams.dp;
      const RealType eps = local_dp * 1e-4f;

      //added
      const RealType sr = this->getParticles()->getSearchRadius();
      const VectorType frameBackSize = getFrameBackDimensions() * sr;
      const VectorType frameBackOrigin =
         getFrameBackOriginGlobalCoordinates() * sr + this->getParticles()->getReferentialOrigin();

      // Returns true if this face should be skipped
      auto isExcluded = [ & ]( const VectorType& normal ) -> bool
      {
         for( int k = 0; k < excludedNormals.getSize(); k++ ) {
            bool match = true;
            for( int d = 0; d < VectorType::getSize(); d++ )
               if( std::fabs( normal[ d ] - excludedNormals[ k ][ d ] ) > eps ) {
                  match = false;
                  break;
               }
            if( match )
               return true;
         }
         return false;
      };

      // Build the normal for a given (axis, sign) so both lambdas use the same logic
      auto faceNormal = [ & ]( int d, int sign ) -> VectorType
      {
         VectorType normal = 0.f;
         normal[ d ] = ( sign < 0 ) ? 1.f : -1.f;
         if( inner_overlap )
            normal[ d ] *= -1.f;
         return normal;
      };
      //-------

      auto planeNodeCount = [ & ]( int faceAxis, int perpAxis ) -> IndexType
      {
         RealType extent = frameBackSize[ perpAxis ];
         for( int d = 0; d < faceAxis; d++ ) {
            if( d == perpAxis )
               continue;
            // Only active (non-excluded) faces claim corner strips
            if( ! isExcluded( faceNormal( d, -1 ) ) || ! isExcluded( faceNormal( d, +1 ) ) )
               extent -= 2.f * local_dp;
         }
         // One node row must be dropped so that no node lies exactly on the frameBack
         // boundary plane (frameBack is cell-aligned with the patch grid + overlap, so
         // particles created from that row would fall outside the searchable grid).
         const IndexType n = static_cast< IndexType >( TNL::max( 0.f, extent ) / local_dp );
         return ( n > 0 ) ? n - 1 : 0;
      };

      auto perpOriginOffset = [ & ]( int faceAxis, int perpAxis ) -> RealType
      {
         RealType offset = 0.f;
         for( int d = 0; d < faceAxis; d++ ) {
            if( d == perpAxis )
               continue;
            if( ! isExcluded( faceNormal( d, -1 ) ) || ! isExcluded( faceNormal( d, +1 ) ) )
               offset += local_dp;
         }
         return offset;
      };

      // Count total nodes across all 2*dim faces
      IndexType n_massNodes = 0;
      for( int d = 0; d < VectorType::getSize(); d++ ) {
         for( int sign : { -1, +1 } ) {
            if( isExcluded( faceNormal( d, sign ) ) )
               continue;  // ← only change
            IndexType faceNodes = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  faceNodes *= planeNodeCount( d, pd );
            n_massNodes += faceNodes;
         }
      }
      massNodes.setSize( n_massNodes );
      auto points_view = massNodes.points.getView();
      auto normals_view = massNodes.normal.getView();

      // Generate nodes face by face
      IndexType offset = 0;
      for( int d = 0; d < VectorType::getSize(); d++ ) {
         for( int sign : { -1, +1 } ) {
            if( isExcluded( faceNormal( d, sign ) ) )
               continue;

            // Inward normal
            VectorType normal = 0.f;
            normal[ d ] = ( sign < 0 ) ? 1.f : -1.f;
            if( inner_overlap )
               normal[ d ] *= -1.f;

            // Physical coordinate on the interface axis
            const RealType ifaceCoord = ( sign < 0 ) ? frameBackOrigin[ d ] : frameBackOrigin[ d ] + frameBackSize[ d ];

            // parallelFor range: [0,1) on d, [0, count) on perp axes
            CoordinatesType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  end[ pd ] = planeNodeCount( d, pd );

            // Strides for linearisation (perp axes only, reverse order)
            CoordinatesType stride = 0;
            {
               IndexType running = 1;
               for( int pd = VectorType::getSize() - 1; pd >= 0; pd-- ) {
                  if( pd == d )
                     continue;
                  stride[ pd ] = running;
                  running *= end[ pd ];
               }
            }

            // Physical start of the node grid on each perp axis
            //VectorType perpStart = frameFrontOrigin;
            VectorType perpStart = frameBackOrigin;
            perpStart[ d ] = ifaceCoord;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  perpStart[ pd ] += perpOriginOffset( d, pd );

            //TODO: Can I just the variables present in the scope or do I need this?
            const IndexType faceOffset = offset;
            const int iAxis = d;
            const VectorType norm = normal;
            const VectorType pStart = perpStart;

            auto generate = [ = ] __cuda_callable__( const CoordinatesType idx ) mutable
            {
               // Linearise: sum over perp axes only
               IndexType i = faceOffset;
               for( int pd = 0; pd < VectorType::getSize(); pd++ )
                  i += idx[ pd ] * stride[ pd ];

               // Physical position
               VectorType r = pStart;
               for( int pd = 0; pd < VectorType::getSize(); pd++ )
                  if( pd != iAxis )
                     r[ pd ] += local_dp * ( idx[ pd ] + 1 );
               //r[ pd ] += local_dp * ( idx[ pd ] ); //TODO: +1 should the +1 be there?

               points_view[ i ] = r;
               normals_view[ i ] = norm;
            };
            Algorithms::parallelFor< DeviceType >( begin, end, generate );

            // Advance offset by number of nodes on this face
            IndexType faceNodes = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  faceNodes *= end[ pd ];
            offset += faceNodes;
         }
      }

      const std::string outputFileName = "results/massNodes_subdomain" + std::to_string( subdomainIdx ) + "_to_"
                                       + std::to_string( neighborIdx ) + "_" + std::to_string( local_dp ) + ".vtk";
      writeMassNodesToVTK( outputFileName );
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
      //const RealType dx = refinementFactor * modelParams.dp * 0.5f;  //FIXME I have no clue why there is 0.5 factor
      //const RealType dx =
      //   //std::pow( refinementFactor * modelParams.dp * 0.5f, dim - 1 );  //FIXME I have no clue why there is 0.5 factor
      //   std::pow( refinementFactor * modelParams.dp, dim - 1 );  //FIXME I have no clue why there is 0.5 factor
      //----
      const RealType trueRefinementFactor = this->getParticles()->getSearchRadius() / ( 2.f * modelParams.h );
      const RealType dx = std::pow( trueRefinementFactor * modelParams.dp, dim - 1 );
      const RealType dx_nb = std::pow( refinementFactor * modelParams.dp, dim - 1 );
      const RealType dx_min = std::min( dx, dx_nb );
      //----
      const VectorType v_subdomain = 0.f;  // velocity of moving subdomain
      //const RealType div_r_trashold = 1.5f;  //FIXME add to model params, depends on dimension
      const RealType div_r_trashold = ( dim == 2 ) ? 1.f : 1.5f;
      const RealType extrapolationDetTreshold = modelParams.mdbcExtrapolationDetTreshold;  //TODO: rename, remove mdbc

      auto interpolateFluid = [ = ] __cuda_callable__(
                                 IndexType i,
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

         ParticlesType::NeighborsLoopAnotherSet::exec(
            i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bv_x, &div_r_x, &drs_min );

         RealType rho_x = 0.f;
         VectorType v_x = 0.f;

         //TODO: Use LU Decomposition so we can just reuse it different RHS
         // Same conditioning as in interpolateVariables: a sparse one-sided
         // support makes M singular and solve() then poisons the mass node
         // with Inf/NaN - the NaN mass would never spawn again
         const RealType detM = Matrices::determinant( M_x );
         if( M_x( 0, 0 ) > 0.05 && detM > 0.001f ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = Matrices::solve( M_x, bv_x[ d ] )[ 0 ];
         }
         else if( M_x( 0, 0 ) > 0.05f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = bv_x[ d ][ 0 ] / M_x( 0, 0 );
         }
         else {
            // no reliable support - keep rho_x = 0, v_x = 0 (zero flux)
         }

         // Ricci et al. uses (-1) * rho, but I think I have different normal orientation
         const RealType m_x_lessZero = TNL::max( 0, rho_x * ( v_x - v_subdomain, normal_x ) * dx * dt );
         const RealType m_x_geqZero = rho_x * ( v_x - v_subdomain, normal_x ) * dx * dt;

         // sum up the flux with free surface corrections
         RealType m_flux_x = 0;

         //TODO: Should I use ( view_mFlux_massPoints[ i ] < 0.f  ) here?
         if( ( div_r_x > div_r_trashold || drs_min < dx_min ) && ( view_m_massPoints[ i ] < 0.f ) )
            m_flux_x = m_x_lessZero;
         else if( div_r_x > div_r_trashold || drs_min < dx_min )
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

      auto interpolateFluid = [ = ] __cuda_callable__(
                                 IndexType i,
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

         ParticlesType::NeighborsLoopAnotherSet::exec( i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bv_x );

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
            // dry support: keep the last (captured / creation-initialized) values.
            // Deleting the buffer particle here absorbs the mass that crosses the
            // interface and the dry side of the face never seeds; the particle keeps
            // its values until the neighbor fluid shows up and refreshes them.
            return 1;
         }
      };
      const IndexType numberOfDryBufferParticles = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, particleLoop );
      (void) numberOfDryBufferParticles;
   }

   /* Band state refresh from the coarser neighbor (LTS schedule v2). Interpolates rho
      and v from the neighbor fluid - which has already advanced - and blends them with
      the past-endpoint anchors stored in rho_old/v_old (sync-time state, or the state
      re-stamped at the last clock coincidence):
         q <- theta * q_interp( now ) + ( 1 - theta ) * q_anchor.
      theta = 1 is the exact coincidence refresh. Positions are not changed - the band
      particles have advected to their current positions with the fine level. No
      particle is created, converted or removed here (membership changes happen only at
      sync points): particles without interpolation support keep their current values. */
   template< typename FluidPointer, typename ModelParams >
   void
   refreshVariablesFromNeighbor( FluidPointer& fluid_neihgbor, ModelParams& modelParams, const RealType theta )
   {
      auto searchInFluid = this->getParticles()->getSearchToken( fluid_neihgbor->getParticles() );
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();

      auto view_points_overlap = this->getParticles()->getPoints().getView();
      auto view_rho_overlap = this->getVariables()->rho.getView();
      auto view_v_overlap = this->getVariables()->v.getView();
      const auto view_rho_anchor = this->getIntegratorVariables()->rho_old.getConstView();
      const auto view_v_anchor = this->getIntegratorVariables()->v_old.getConstView();
      const auto view_points_fluid = fluid_neihgbor->getParticles()->getPoints().getConstView();
      const auto view_rho_fluid = fluid_neihgbor->getVariables()->rho.getConstView();
      const auto view_v_fluid = fluid_neihgbor->getVariables()->v.getConstView();

      const unsigned int dim = SPHCaseConfig::spaceDimension;
      const RealType searchRadius = fluid_neihgbor->getParticles()->getSearchRadius();
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      const RealType h = refinementFactor * modelParams.h;
      const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
      const RealType blend = theta;
      const RealType keep = 1.f - theta;

      auto interpolateFluid = [ = ] __cuda_callable__(
                                 IndexType i,
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

         ParticlesType::NeighborsLoopAnotherSet::exec( i, r_x, searchInFluid, interpolateFluid, &M_x, &brho_x, &bv_x );

         RealType rho_x;
         VectorType v_x;

         //TODO: Use LU Decomposition so we can just reuse it different RHS
         //TODO: Proper condition should be if( std::fabs( Matrices::determinant( M_x ) ) > extrapolationDetTreshold )
         const RealType detM = Matrices::determinant( M_x );
         if( M_x( 0, 0 ) > 0.05 && detM > 0.001f ) {
            rho_x = Matrices::solve( M_x, brho_x )[ 0 ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = Matrices::solve( M_x, bv_x[ d ] )[ 0 ];

            view_rho_overlap[ i ] = blend * rho_x + keep * view_rho_anchor[ i ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               view_v_overlap[ i ][ d ] = blend * v_x[ d ] + keep * view_v_anchor[ i ][ d ];
         }
         else if( M_x( 0, 0 ) > 0.05f ) {
            rho_x = brho_x[ 0 ] / M_x( 0, 0 );
            for( int d = 0; d < VectorType::getSize(); d++ )
               v_x[ d ] = bv_x[ d ][ 0 ] / M_x( 0, 0 );

            view_rho_overlap[ i ] = blend * rho_x + keep * view_rho_anchor[ i ];
            for( int d = 0; d < VectorType::getSize(); d++ )
               view_v_overlap[ i ][ d ] = blend * v_x[ d ] + keep * view_v_anchor[ i ][ d ];
         }
         // else: no interpolation support - keep the current state until the next sync
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, particleLoop );
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
      advectBufferParticles( dt );
      classifyParticlesByFrames();
   }

   /* Ballistic advection of the whole band (both layers drift with their current velocity). */
   void
   advectBufferParticles( const RealType dt )
   {
      auto view_r_buffer = this->getParticles()->getPoints().getView();
      const auto view_v_buffer = this->getVariables()->v.getConstView();
      const IndexType n_buffer = this->getParticles()->getNumberOfParticles();

      auto advectBufferParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      };
      Algorithms::parallelFor< DeviceType >( 0, n_buffer, advectBufferParticles );
   }

   /* Frame-based retype/remove classification on the overlap-width-aware frameBack box. */
   void
   classifyParticlesByFrames()
   {
      const IndexType n_buffer = this->getParticles()->getNumberOfParticles();
      auto view_r_buffer = this->getParticles()->getPoints().getView();

      // reset retype marker
      auto retypeMarker_view = retypeMarker.getView();
      retypeMarker_view = 0;

      const RealType sr = this->getParticles()->getSearchRadius();
      const RealType inv_sr = 1.f / sr;
      const VectorType refOrig = this->getParticles()->getReferentialOrigin();
      const CoordinatesType frameFrontOriginGlobalCoords = getFrameFrontOriginGlobalCoordinates();
      const CoordinatesType frameBackGlobCoords = getFrameBackOriginGlobalCoordinates();
      const CoordinatesType frameBackDims = getFrameBackDimensions();
      const CoordinatesType frameFrontDims = getFrameFrontDimensions();
      const bool inner_overlap = this->inner_overlap;
      const bool outer_overlap = this->outer_overlap;

      // Retype to fluid if:
      // - outer zone + is inside
      // - inner zone + is outside
      auto identifyRetype = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         //const CoordinatesType gc = TNL::floor( ( view_r_buffer[ i ] - frameFrontOrigin ) * inv_sr );
         const CoordinatesType giGlobCoords = TNL::floor( ( view_r_buffer[ i ] - refOrig ) * inv_sr );
         const CoordinatesType gc = giGlobCoords - frameFrontOriginGlobalCoords;
         bool inside = isInsideBox( gc, frameFrontDims );

         if( inside && outer_overlap ) {
            retypeMarker_view[ i ] = 1;
            return 1;
         }
         else if( ! inside && inner_overlap ) {
            retypeMarker_view[ i ] = 1;
            return 1;
         }
         return 0;
      };
      numberOfParticlesToRetype = Algorithms::reduce< DeviceType >( 0, n_buffer, identifyRetype );

      // Remove invalid buffe particles
      // - outer zone - is outside the ( subdomain box + zone width )
      // - inner zone - is insde the ( subdomain box - zone width )
      auto identifyRemove = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         //const CoordinatesType gc = TNL::floor( ( view_r_buffer[ i ] - frameBackOrigin ) * inv_sr );
         const CoordinatesType giGlobCoords = TNL::floor( ( view_r_buffer[ i ] - refOrig ) * inv_sr );
         const CoordinatesType gc = giGlobCoords - frameBackGlobCoords;
         bool inside = isInsideBox( gc, frameBackDims );
         if( inside && inner_overlap ) {
            retypeMarker_view[ i ] = 2;
            return 1;
         }
         else if( ! inside && outer_overlap ) {
            retypeMarker_view[ i ] = 2;
            return 1;
         }
         return 0;
      };
      numberOfParticlesToRemove = Algorithms::reduce< DeviceType >( 0, n_buffer, identifyRemove );
   }

   void
   sortBufferParticles()
   {
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();
      auto r_view = this->getParticles()->getPoints().getView();
      auto v_view = this->getVariables()->v.getView();
      auto rho_view = this->getVariables()->rho.getView();
      auto drho_view = this->getVariables()->drho.getView();
      auto p_view = this->getVariables()->p.getView();
      auto a_view = this->getVariables()->a.getView();
      auto gamma_view = this->getVariables()->gamma.getView();
      auto marker_view = this->getVariables()->marker.getView();
      auto referentialIdx_view = this->getVariables()->referentialIdx.getView();
      auto retypeMarker_view = this->retypeMarker.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;

      thrust::sort_by_key(
         thrustDevice,
         retypeMarker_view.getArrayData(),
         retypeMarker_view.getArrayData() + numberOfBufferParticles,
         thrust::make_zip_iterator( cuda::std::make_tuple(
            r_view.getArrayData(),
            v_view.getArrayData(),
            rho_view.getArrayData(),
            drho_view.getArrayData(),
            p_view.getArrayData(),
            a_view.getArrayData(),
            gamma_view.getArrayData(),
            marker_view.getArrayData(),
            referentialIdx_view.getArrayData() ) ) );
   }

   template< typename FluidPointer >
   void
   convertBufferToFluid( FluidPointer& fluid )
   {
      convertBufferToFluidImpl( fluid, false, 0.f );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   convertBufferToFluid( FluidPointer& fluid, ModelParams& modelParams )
   {
      convertBufferToFluidImpl( fluid, true, modelParams.rho0 );
   }

   void
   removeBufferParticles()
   {
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();
      this->getParticles()->setNumberOfParticles( numberOfBufferParticles - this->numberOfParticlesToRemove );
   }

   template< typename ModelParams >
   void
   updateMassNodes( ModelParams& modelParams, const RealType dt )
   {
      const auto massFlux_view = massNodes.massFlux.getConstView();
      auto mass_view = massNodes.mass.getView();
      auto view_particlesToCreate = massNodes.particlesToCreate.getView();

      //FIXME: I need refinemet factor here! Subdomain mass is different. But with refinement factor, it doesnt work!
      const RealType refinementFactor = this->getParticles()->getSearchRadius() / ( 2.f * modelParams.h );
      const RealType particleMass = std::pow( refinementFactor, ParticlesType::spaceDimension ) * modelParams.mass;
      //const RealType particleMass =
      //   ( ParticlesType::getParticlesDimension() == 2 ) ? 0.25f * modelParams.mass : 0.125f * modelParams.mass;
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
      this->numberOfParticlesToCreate = Algorithms::reduce< DeviceType >( 0, numberOfMassNodes, identifyParticlesToCreate );
   }

   template< typename ModelParams >
   void
   createBufferParticles( ModelParams& modelParams )
   {
      (void) ensureBufferCapacity( this->getNumberOfParticles() + this->numberOfParticlesToCreate );
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();

      auto r_buffer_view = this->getParticles()->getPoints().getView();
      auto v_buffer_view = this->getVariables()->v.getView();
      auto rho_buffer_view = this->getVariables()->rho.getView();
      auto drho_buffer_view = this->getVariables()->drho.getView();
      auto p_buffer_view = this->getVariables()->p.getView();
      auto a_buffer_view = this->getVariables()->a.getView();
      auto gamma_buffer_view = this->getVariables()->gamma.getView();
      auto marker_buffer_view = this->getVariables()->marker.getView();
      auto referentialIdx_buffer_view = this->getVariables()->referentialIdx.getView();
      auto rho_old_buffer_view = this->getIntegratorVariables()->rho_old.getView();
      auto v_old_buffer_view = this->getIntegratorVariables()->v_old.getView();
      const auto r_massNodes_view = this->massNodes.points.getConstView();
      const auto normal_massNodes_view = this->massNodes.normal.getConstView();
      const RealType refinementFactor = this->getParticles()->getSearchRadius() / ( 2.f * modelParams.h );
      const RealType dp = refinementFactor * modelParams.dp;
      const RealType rho0 = modelParams.rho0;
      const IndexType firstReferentialIdx = this->getVariables()->highestReferentialIdx;

      auto createNewBufferParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType target = numberOfBufferParticles + i;
         const VectorType r_new = r_massNodes_view[ i ] + ( 0.5f * dp ) * normal_massNodes_view[ i ];
         r_buffer_view[ target ] = r_new;
         // initialize values so a failed first interpolation keeps sane
         // last values instead of uninitialized memory
         rho_buffer_view[ target ] = rho0;
         v_buffer_view[ target ] = 0.f;
         drho_buffer_view[ target ] = 0.f;
         p_buffer_view[ target ] = 0.f;
         a_buffer_view[ target ] = 0.f;
         gamma_buffer_view[ target ] = 0.f;
         marker_buffer_view[ target ] = 0;
         referentialIdx_buffer_view[ target ] = firstReferentialIdx + i;
         rho_old_buffer_view[ target ] = rho0;
         v_old_buffer_view[ target ] = 0.f;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfParticlesToCreate, createNewBufferParticles );
      this->getVariables()->highestReferentialIdx = firstReferentialIdx + this->numberOfParticlesToCreate;
      this->getParticles()->setNumberOfParticles( numberOfBufferParticles + this->numberOfParticlesToCreate );
   }

   template< typename FluidPointer >
   void
   getFluidParticlesEneringTheBuffer( FluidPointer& fluid )
   {
      // reset the marker field
      auto particlesToBuffer_view = this->particlesToBuffer.getView();
      particlesToBuffer_view = INT_MAX;

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto zoneParticleIndices_view = this->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = this->zone.getNumberOfParticles();

      const CoordinatesType frameFrontDims = getFrameFrontDimensions();
      //const VectorType frameFrontOrigin = this->frameFrontOrigin;
      const bool inner_overlap = this->inner_overlap;
      const bool outer_overlap = this->outer_overlap;
      const RealType sr = this->getParticles()->getSearchRadius();
      const RealType inv_sr = 1.f / sr;
      const VectorType refOrig = this->getParticles()->getReferentialOrigin();
      //const CoordinatesType gridOriginGlobCoords = this->getParticles()->getGlobalOriginCoordinates();
      const CoordinatesType frameFrontOriginGlobCoords = getFrameFrontOriginGlobalCoordinates();
      //const CoordinatesType gridOriginGlobCoordsWithOverlap = this->particles->getGlobalOriginCoordinatesWithOverlap();

      // Retype fluid to buffer:
      // - outer zone and is outside the frame
      // - inner zone and is inside the frame
      auto checkFluidParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         // A degenerate (NaN/Inf) position must never be captured into the band:
         // its cell index is undefined and the raw copy would seed the buffer
         // with NaN, which conversions then recycle back into the fluid forever.
         // Skip the capture - removeParticlesOutOfGrid kills such particles.
         bool finite = true;
         for( int d = 0; d < VectorType::getSize(); d++ )
            if( ! ( fabsf( r[ d ] ) < 1e30f ) )
               finite = false;
         if( ! finite )
            return 0;
         const CoordinatesType giGlobCoords = TNL::floor( ( r - refOrig ) * inv_sr );
         const CoordinatesType gc = giGlobCoords - frameFrontOriginGlobCoords;
         const bool inside = isInsideBox( gc, frameFrontDims );

         if( inside && inner_overlap ) {
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         else if( ! inside && outer_overlap ) {
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         return 0;
      };
      this->numberOfParticlesToBuffer =
         Algorithms::reduce< DeviceType >( 0, numberOfZoneParticles, checkFluidParticles, TNL::Plus() );

      // sort the indices
      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::sort(
         thrustDevice, particlesToBuffer_view.getArrayData(), particlesToBuffer_view.getArrayData() + numberOfZoneParticles );
   }

   template< typename FluidPointer >
   void
   convertFluidToBuffer( FluidPointer& fluid )
   {
      convertFluidToBufferImpl( fluid, false, 0.f );
   }

   template< typename FluidPointer, typename ModelParams >
   void
   convertFluidToBuffer( FluidPointer& fluid, ModelParams& modelParams )
   {
      convertFluidToBufferImpl( fluid, true, modelParams.rho0 );
   }

   // Marks particles whose cell lies outside the own grid including overlap, so the
   // cell-list build (which has no range check) cannot touch out-of-bounds memory.
   template< typename ParticlesPointer >
   IndexType
   removeParticlesOutOfGrid( ParticlesPointer& particles )
   {
      return BaseType::removeParticlesOutOfGrid( particles );
   }

   /* New buffer particles (spawned or captured) inherit consistent old-state for the
      Verlet kick of the next cycle, same convention as convertBufferToFluid. */
   void
   initializeIntegratorVariables()
   {
      auto view_v = this->getVariables()->v.getView();
      auto view_rho = this->getVariables()->rho.getView();
      auto view_v_old = this->getIntegratorVariables()->v_old.getView();
      auto view_rho_old = this->getIntegratorVariables()->rho_old.getView();

      const IndexType n_buffer = this->getParticles()->getNumberOfParticles();
      auto initialize = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         view_rho_old[ i ] = view_rho[ i ];
         view_v_old[ i ] = view_v[ i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, n_buffer, initialize );
   }

   /* Sync-point lifecycle. Unlike the monolithic single-step update, no mass accumulation
      and no particle advection happens here: accumulateMasses is scheduled by the solver
      per substep (with the substep dt) and positions are already at the new sync time. */
   template< typename FluidPointer, typename ModelParams >
   void
   syncUpdateInterfaceBuffer(
      FluidPointer& fluid_own,
      FluidPointer& fluid_neihgbor,
      ModelParams& modelParams,
      const RealType dt )
   {
      removeParticlesOutOfGrid( fluid_neihgbor->getParticles() );
      fluid_neihgbor->searchForNeighbors();

      classifyParticlesByFrames();
      sortBufferParticles();
      removeBufferParticles();
      convertBufferToFluid( fluid_own, modelParams );
      removeParticlesOutOfGrid( fluid_own->getParticles() );
      fluid_own->searchForNeighbors();
      zone.updateParticlesInZone( fluid_own->getParticles() );
      getFluidParticlesEneringTheBuffer( fluid_own );
      convertFluidToBuffer( fluid_own, modelParams );
      updateMassNodes( modelParams, dt );
      massNodes.massFlux = 0;
      massNodes.sort();
      createBufferParticles( modelParams );
      interpolateVariables( fluid_neihgbor, modelParams );
      removeParticlesOutOfGrid( this->getParticles() );
      this->searchForNeighbors();
      initializeIntegratorVariables();

      numberOfParticlesToRemove = 0;
      numberOfParticlesToRetype = 0;
      numberOfParticlesToCreate = 0;
      numberOfParticlesToBuffer = 0;
   }

   template< typename FluidPointer, typename ModelParams >
   void
   updateInterfaceBuffer(
      FluidPointer& fluid_own,
      FluidPointer& fluid_neihgbor,
      const ModelParams& modelParams,
      const RealType dt,
      const int subdomainIdx )
   {
      removeParticlesOutOfGrid( fluid_neihgbor->getParticles() );
      fluid_neihgbor->searchForNeighbors();

      massNodes.massFlux = 0;
      accumulateMasses( fluid_neihgbor, modelParams, dt );
      moveBufferParticles( dt );
      sortBufferParticles();
      removeBufferParticles();
      convertBufferToFluid( fluid_own );
      removeParticlesOutOfGrid( fluid_own->getParticles() );
      fluid_own->searchForNeighbors();
      zone.updateParticlesInZone( fluid_own->getParticles() );
      getFluidParticlesEneringTheBuffer( fluid_own );
      convertFluidToBuffer( fluid_own );
      updateMassNodes( modelParams, dt );
      massNodes.sort();
      createBufferParticles( modelParams );
      interpolateVariables( fluid_neihgbor, modelParams );
      removeParticlesOutOfGrid( this->getParticles() );
      this->searchForNeighbors();

      numberOfParticlesToRemove = 0;
      numberOfParticlesToRetype = 0;
      numberOfParticlesToCreate = 0;
      numberOfParticlesToBuffer = 0;
   }

   // links
   [[nodiscard]] const CoordinatesType
   getFrameFrontOriginGlobalCoordinates() const
   {
      return frameOriginGlobalCoordinates;
   }

   [[nodiscard]] const CoordinatesType
   getFrameFrontDimensions() const
   {
      return frameDimensions;
   }

   [[nodiscard]] const CoordinatesType
   getFrameBackOriginGlobalCoordinates() const
   {
      return frameOriginGlobalCoordinates + ( -1 ) * this->getParticles()->getOverlapWidth() * frameOrientation;
   }

   [[nodiscard]] const CoordinatesType
   getFrameBackDimensions() const
   {
      return frameDimensions + 2 * this->getParticles()->getOverlapWidth() * frameOrientation;
   }

   /* Intermediate frame splitting the band into layer 1 (frameMid XOR frameFront) and
      layer 2 (frameBack XOR frameMid). Coincides with frameFront for overlap width 1. */
   [[nodiscard]] const CoordinatesType
   getFrameMidOriginGlobalCoordinates() const
   {
      return frameOriginGlobalCoordinates + ( -1 ) * ( this->getParticles()->getOverlapWidth() - 1 ) * frameOrientation;
   }

   [[nodiscard]] const CoordinatesType
   getFrameMidDimensions() const
   {
      return frameDimensions + 2 * ( this->getParticles()->getOverlapWidth() - 1 ) * frameOrientation;
   }

   void
   writeProlog( TNL::Logger& logger, const int subdomainIdx )
   {
      logger.writeParameter( "Subdomain index:", subdomainIdx );
      BaseType::writeProlog( logger );
      logger.writeSeparator();

      // Buffer configuration
      logger.writeParameter( "Buffer width:", this->bufferWidth );
      logger.writeSeparator();

      // Overlap type
      logger.writeParameter( "Overlap type:", inner_overlap ? "INNER" : outer_overlap ? "OUTER" : "INVALID" );
      logger.writeParameter( "Inner overlap:", inner_overlap );
      logger.writeParameter( "Outer overlap:", outer_overlap );
      logger.writeSeparator();

      // Frame front configuration
      logger.writeParameter(
         "Frame front origin:",
         this->getParticles()->getReferentialOrigin()
            + getFrameFrontOriginGlobalCoordinates() * this->getParticles()->getSearchRadius() );
      logger.writeParameter(
         "Frame front origin coords:",
         getFrameFrontOriginGlobalCoordinates() - this->getParticles()->getGlobalOriginCoordinates() );
      logger.writeParameter( "Frame front origin global coords:", getFrameFrontOriginGlobalCoordinates() );
      logger.writeParameter( "Frame front dims:", getFrameFrontDimensions() );
      logger.writeParameter(
         "Frame front end:",
         getFrameFrontOriginGlobalCoordinates() + getFrameFrontDimensions()
            - this->getParticles()->getGlobalOriginCoordinates() );
      logger.writeParameter( "Frame orientation:", frameOrientation );
      logger.writeSeparator();

      // Frame back configuration
      logger.writeParameter(
         "Frame back origin:",
         this->getParticles()->getReferentialOrigin()
            + getFrameBackOriginGlobalCoordinates() * this->getParticles()->getSearchRadius() );
      logger.writeParameter( "Frame back size:", getFrameBackDimensions() * this->getParticles()->getSearchRadius() );
      logger.writeParameter( "Frame back dims:", getFrameBackDimensions() );
      logger.writeSeparator();

      // Frame mid configuration (LTS band split between layer 1 and layer 2)
      logger.writeParameter( "Frame mid origin global coords:", getFrameMidOriginGlobalCoordinates() );
      logger.writeParameter( "Frame mid dims:", getFrameMidDimensions() );
      logger.writeSeparator();

      // Zone information
      zone.writeProlog( logger );
   }

   /**
    * \brief Collects the multiresolution buffer patch geometry into the
    *        initialization report.
    *
    * Extends the base ParticleSet report (the buffer band's own counts and
    * grid geometry) with the interface geometry as loaded by initZones() -
    * the values as actually computed and stored there, read back through
    * the class' own members and the same frame getters writeProlog prints.
    * Vector entries are three-component JSON arrays (z = 0 in 2D).
    *
    * \returns A JSON object with the ParticleSet base keys plus \e level,
    *          \e resolution_factor, \e own_particle_origin,
    *          \e nb_particle_origin, \e own_global_origin_coordinates,
    *          \e nb_global_origin_coordinates,
    *          \e frame_front_origin_global_coordinates,
    *          \e frame_back_origin_global_coordinates,
    *          \e frame_mid_origin_global_coordinates,
    *          \e frame_front_dimensions, \e frame_back_dimensions,
    *          \e frame_mid_dimensions, \e frame_width, \e zone_dimensions, \e zone_cell_count
    *          and \e number_of_mass_nodes.
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      nlohmann::json report = BaseType::writeInitReport();
      report[ "level" ] = refinementLevel;
      report[ "resolution_factor" ] = resolutionFactor;
      report[ "own_particle_origin" ] = this->vectorToInitReportArray( ownParticleOrigin );
      report[ "nb_particle_origin" ] = this->vectorToInitReportArray( nbParticleOrigin );
      report[ "own_global_origin_coordinates" ] = this->vectorToInitReportArray( ownGlobalOriginCoordinates );
      report[ "nb_global_origin_coordinates" ] = this->vectorToInitReportArray( nbGlobalOriginCoordinates );
      report[ "frame_front_origin_global_coordinates" ] =
         this->vectorToInitReportArray( getFrameFrontOriginGlobalCoordinates() );
      report[ "frame_back_origin_global_coordinates" ] = this->vectorToInitReportArray( getFrameBackOriginGlobalCoordinates() );
      report[ "frame_mid_origin_global_coordinates" ] = this->vectorToInitReportArray( getFrameMidOriginGlobalCoordinates() );
      report[ "frame_front_dimensions" ] = this->vectorToInitReportArray( getFrameFrontDimensions() );
      report[ "frame_back_dimensions" ] = this->vectorToInitReportArray( getFrameBackDimensions() );
      report[ "frame_mid_dimensions" ] = this->vectorToInitReportArray( getFrameMidDimensions() );
      report[ "frame_width" ] = frameWidth;
      report[ "zone_dimensions" ] = this->vectorToInitReportArray( zoneDimensions );
      report[ "zone_cell_count" ] = zone.getNumberOfCells();
      report[ "number_of_mass_nodes" ] = massNodes.numberOfMassNodes;
      return report;
   }

   void
   writeMassNodesToVTK( const std::string& outputFileName )
   {
      const int n_mn = massNodes.numberOfMassNodes;
      ParticlesType nodes;
      nodes.setSize( n_mn );
      nodes.setNumberOfParticles( n_mn );
      nodes.getPoints() = massNodes.points;

      using WriterType = TNL::Particles::Writers::VTKWriter< ParticlesType >;
      std::ofstream outputFileFluid( outputFileName, std::ofstream::out );
      WriterType writer( outputFileFluid );
      writer.writeParticles( nodes );
      writer.template writePointData< typename MassNodes::VectorArrayType >( massNodes.normal, "Normal", n_mn );
   }

   [[nodiscard]] const MassNodes&
   getMassNodes() const
   {
      return massNodes;
   }

public:
   // public: nvcc rejects extended __cuda_callable__ lambdas in private member functions
   template< typename FluidPointer >
   void
   convertBufferToFluidImpl( FluidPointer& fluid, const bool sanitizeVelocity, const RealType rho0 )
   {
      (void) fluid->ensureCapacity( fluid->getNumberOfParticles() + this->numberOfParticlesToRetype );
      const IndexType numberOfFluidParticles = fluid->getNumberOfParticles();
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();
      const IndexType toFluidOffset = numberOfBufferParticles - this->numberOfParticlesToRetype;

      auto view_r_fluid = fluid->getParticles()->getPoints().getView();
      auto view_v_fluid = fluid->getVariables()->v.getView();
      auto view_rho_fluid = fluid->getVariables()->rho.getView();
      auto view_drho_fluid = fluid->getVariables()->drho.getView();
      auto view_p_fluid = fluid->getVariables()->p.getView();
      auto view_a_fluid = fluid->getVariables()->a.getView();
      auto view_gamma_fluid = fluid->getVariables()->gamma.getView();
      auto view_marker_fluid = fluid->getVariables()->marker.getView();
      auto view_referentialIdx_fluid = fluid->getVariables()->referentialIdx.getView();
      auto view_rho_old = fluid->getIntegratorVariables()->rho_old.getView();
      auto view_v_old = fluid->getIntegratorVariables()->v_old.getView();

      auto view_r_buffer = this->getParticles()->getPoints().getView();
      const auto view_v_buffer = this->getVariables()->v.getConstView();
      const auto view_rho_buffer = this->getVariables()->rho.getConstView();
      const auto view_marker_buffer = this->getVariables()->marker.getConstView();
      const auto view_referentialIdx_buffer = this->getVariables()->referentialIdx.getConstView();

      auto createNewFluidParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType source = toFluidOffset + i;
         const IndexType target = numberOfFluidParticles + i;
         RealType rho = view_rho_buffer[ source ];
         VectorType v = view_v_buffer[ source ];
         if( sanitizeVelocity ) {
            // Only the ModelParams wrapper repairs NaN band velocity at hand-off.
            bool bad = false;
            for( int d = 0; d < VectorType::getSize(); d++ )
               if( v[ d ] != v[ d ] )
                  bad = true;
            if( bad ) {
               rho = rho0;
               v = 0.f;
            }
         }

         view_r_fluid[ target ] = view_r_buffer[ source ];
         view_rho_fluid[ target ] = rho;
         view_v_fluid[ target ] = v;
         view_drho_fluid[ target ] = 0.f;
         view_p_fluid[ target ] = 0.f;
         view_a_fluid[ target ] = 0.f;
         view_gamma_fluid[ target ] = 0.f;
         view_marker_fluid[ target ] = view_marker_buffer[ source ];
         view_referentialIdx_fluid[ target ] = view_referentialIdx_buffer[ source ];

         // A resolution hand-off restarts Verlet from the transferred current state.
         view_rho_old[ target ] = rho;
         view_v_old[ target ] = v;

         //TODO: I need to remove the used particles, not sure if this is the right way
         view_r_buffer[ source ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfParticlesToRetype, createNewFluidParticles );
      fluid->getVariables()->highestReferentialIdx = std::max(
         fluid->getVariables()->highestReferentialIdx, this->getVariables()->highestReferentialIdx );
      fluid->getParticles()->setNumberOfParticles( numberOfFluidParticles + this->numberOfParticlesToRetype );
      this->getParticles()->setNumberOfParticles( numberOfBufferParticles - this->numberOfParticlesToRetype );
   }

   template< typename FluidPointer >
   void
   convertFluidToBufferImpl( FluidPointer& fluid, const bool sanitizeVelocity, const RealType rho0 )
   {
      (void) ensureBufferCapacity( this->getNumberOfParticles() + this->numberOfParticlesToBuffer );
      const IndexType numberOfBufferParticles = this->getNumberOfParticles();

      auto r_fluid_view = fluid->getParticles()->getPoints().getView();
      const auto v_fluid_view = fluid->getVariables()->v.getConstView();
      const auto rho_fluid_view = fluid->getVariables()->rho.getConstView();
      const auto marker_fluid_view = fluid->getVariables()->marker.getConstView();
      const auto referentialIdx_fluid_view = fluid->getVariables()->referentialIdx.getConstView();

      auto r_buffer_view = this->getParticles()->getPoints().getView();
      auto v_buffer_view = this->getVariables()->v.getView();
      auto rho_buffer_view = this->getVariables()->rho.getView();
      auto drho_buffer_view = this->getVariables()->drho.getView();
      auto p_buffer_view = this->getVariables()->p.getView();
      auto a_buffer_view = this->getVariables()->a.getView();
      auto gamma_buffer_view = this->getVariables()->gamma.getView();
      auto marker_buffer_view = this->getVariables()->marker.getView();
      auto referentialIdx_buffer_view = this->getVariables()->referentialIdx.getView();
      auto rho_old_buffer_view = this->getIntegratorVariables()->rho_old.getView();
      auto v_old_buffer_view = this->getIntegratorVariables()->v_old.getView();
      const auto particlesToBuffer_view = this->particlesToBuffer.getConstView();

      auto retypeFluidToBuffer = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType p = particlesToBuffer_view[ i ];
         const IndexType target = numberOfBufferParticles + i;
         RealType rho = rho_fluid_view[ p ];
         VectorType v = v_fluid_view[ p ];
         if( sanitizeVelocity ) {
            // Only the ModelParams wrapper repairs NaN fluid velocity at capture.
            bool bad = false;
            for( int d = 0; d < VectorType::getSize(); d++ )
               if( v[ d ] != v[ d ] )
                  bad = true;
            if( bad ) {
               rho = rho0;
               v = 0.f;
            }
         }

         r_buffer_view[ target ] = r_fluid_view[ p ];
         v_buffer_view[ target ] = v;
         rho_buffer_view[ target ] = rho;
         drho_buffer_view[ target ] = 0.f;
         p_buffer_view[ target ] = 0.f;
         a_buffer_view[ target ] = 0.f;
         gamma_buffer_view[ target ] = 0.f;
         marker_buffer_view[ target ] = marker_fluid_view[ p ];
         referentialIdx_buffer_view[ target ] = referentialIdx_fluid_view[ p ];
         rho_old_buffer_view[ target ] = rho;
         v_old_buffer_view[ target ] = v;
         r_fluid_view[ p ] = FLT_MAX;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToBuffer, retypeFluidToBuffer );
      this->getVariables()->highestReferentialIdx = std::max(
         this->getVariables()->highestReferentialIdx, fluid->getVariables()->highestReferentialIdx );
      this->getParticles()->setNumberOfParticles( numberOfBufferParticles + numberOfParticlesToBuffer );
      fluid->getParticles()->setNumberOfParticlesToRemove(
         fluid->getParticles()->getNumberOfParticlesToRemove() + numberOfParticlesToBuffer );
   }

private:
   void
   ensureBufferCapacity( const IndexType required )
   {
      if( BaseType::ensureCapacity( required ) ) {
         const IndexType capacity = this->getNumberOfAllocatedParticles();
         particlesToFluid.resize( capacity );
         particlesToRemove.resize( capacity );
         // The growth search does not reorder MR scratch. The next classification
         // reinitializes retypeMarker before it is used again.
         retypeMarker.resize( capacity );
      }
   }

protected:
   MassNodes massNodes;

   // temporary constants //TODO: Use names such as fluidToBufferCount
   IndexType numberOfParticlesToRemove = 0;
   IndexType numberOfParticlesToRetype = 0;  //TODO: Rename to "moveToFluid"
   IndexType numberOfParticlesToCreate = 0;
   IndexType numberOfParticlesToBuffer = 0;

   // temporary arrays
   IndexArrayType retypeMarker;  //TODO: rename
   IndexArrayType particlesToFluid;
   IndexArrayType particlesToRemove;
   IndexArrayType particlesToBuffer;

   // Geometric definition of the interface
   bool inner_overlap;
   bool outer_overlap;

   int refinementLevel = 0;
   int frameOrientation;

   CoordinatesType frameOriginGlobalCoordinates;
   CoordinatesType frameDimensions;  // = inner subdomain dimensions
   RealType bufferWidth;             // useful for dual time stepping

   // interface geometry stored by initZones for the init report
   RealType resolutionFactor = 0.f;
   VectorType ownParticleOrigin = 0.f;
   VectorType nbParticleOrigin = 0.f;
   CoordinatesType ownGlobalOriginCoordinates = 0;
   CoordinatesType nbGlobalOriginCoordinates = 0;
   CoordinatesType zoneDimensions = 0;

   ParticleZone zone;
};

template<
   typename ParticlesType,
   typename SPHCaseConfig,
   typename Variables,
   typename IntegratorVariables,
   typename OpenBoundaryConfig,
   typename SPHDefs >
using MultiresolutionBoundaryLTS = MultiresolutionBoundary<
   ParticlesType,
   SPHCaseConfig,
   Variables,
   IntegratorVariables,
   OpenBoundaryConfig,
   SPHDefs >;

}  // namespace TNL::SPH
