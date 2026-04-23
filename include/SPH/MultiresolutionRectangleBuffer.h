#pragma once

#include <climits>
#include <functional>
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

   //static constexpr RealType bufferWidthFactorConst = 1.5f; //FIXME IT CAN NOT BE LARGER THEN THE CELL WIDTH (I would like to have like this, but then, removal of particles needs to be changed)
   static constexpr RealType bufferWidthFactorConst = 1; //FIXME IT CAN NOT BE LARGER THEN THE CELL WIDTH (I would like to have like this, but then, removal of particles needs to be changed)
   static constexpr int frameWidth = 2;

   MultiresolutionBoundary() = default;

   // helper functions
   __cuda_callable__
   static bool
   isInsideBox( const VectorType& point, const VectorType& boxOrigin, const VectorType& boxSize )
   {
      bool isInside = true;
      for( int d = 0; d < VectorType::getSize(); d++ )
         if( ( point[ d ] < boxOrigin[ d ] ) || ( point[ d ] > ( boxOrigin[ d ] + boxSize[ d ] ) ) ) //FIXME FIXME Disgusting trick to prevent PhysCoords instead of GirdCoords
            isInside = false;
      return isInside;
   }

   //FIXME: Do I need the boxOriginCoords or not? (I think i don't need them)
   __cuda_callable__
   static bool
   //isInsideBox( const IndexVectorType& coords, const IndexArrayType& boxOriginCoords, const IndexVectorType& boxDims )
   isInsideBox( const IndexVectorType& coords, const IndexVectorType& boxDims )
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
   initZones( const ParticleSetPointer& ownParticles,
              const ParticleSetPointer& nbParticles,
              const RealType refinementFraction,
              const int maxPtcsPerCell = 75 )
/*
initZonesRectangular( const ParticleSetPointer& ownParticles,
                      const Containers::Array< BufferSideType, Devices::Host >& sides,
                      const VectorType& fineRegionOrigin,
                      const VectorType& fineRegionEnd,
                      const IndexType frameWidth = 2,
                      const int maxPtcsPerCell = 75 )
*/
   {
      //sides_  = sides;

      std::cout << "========================================================================" << std::endl;
      const VectorType globalOrig = ownParticles->getGridReferentialOrigin();
      const VectorType ownOrig = ownParticles->getGridOrigin();
      const VectorType nbOrig = nbParticles->getGridOrigin();
      const IndexVectorType ownDims = ownParticles->getGridDimensions();
      const IndexVectorType nbDims = nbParticles->getGridDimensions();
      const IndexVectorType ownDimsWithOverlap = ownParticles->getGridDimensionsWithOverlap();
      const RealType own_sr = ownParticles->getSearchRadius();
      const RealType nb_sr = nbParticles->getSearchRadius();
      const IndexType overlapWidth = ownParticles->getOverlapWidth();
      const VectorType unitVect = 1;

      //const RealType resolutionFactor = own_sr / nb_sr;
      //int resolutionLevel = 1 / resolutionFactor; //FIXME, This is not correct, should be only 1 or 2, right?
      const RealType resolutionFactor = ( own_sr > nb_sr ) ? 0.5 : 2; //FIXME, This is not correct, should be only 1 or 2, right?
      inner_overlap = isInsideBox( nbOrig, ownOrig, own_sr * ownDims );
      outer_overlap = isInsideBox( ownOrig, nbOrig, nb_sr * nbDims ); //TODO: Just use negation, right?

      std::cout << "globalOrig: " << globalOrig << std::endl;
      std::cout << "ownOrig: " << ownOrig << std::endl;
      std::cout << "nbOrig: " << nbOrig << std::endl;

      std::cout << "ownDims: " << ownDims << std::endl;
      std::cout << "nbDims: " << nbDims << std::endl;
      std::cout << "ownDimsWithOverlap: " << ownDimsWithOverlap << std::endl;

      std::cout << "own_sr: " << own_sr << std::endl;
      std::cout << "nb_sr: " << nb_sr << std::endl;

      std::cout << "overlapWidth: " << overlapWidth << std::endl;
      std::cout << "unitVect: " << unitVect << std::endl;

      std::cout << "resolutionFactor: " << refinementFraction << std::endl;
      std::cout << "inner overlap: " << inner_overlap << std::endl;
      std::cout << "outer overlap: " << outer_overlap << std::endl;

      assert( !(inner_overlap && outer_overlap) && "inner_overlap and outer_overlap cannot both be true!" );

      bufferWidth = bufferWidthFactorConst * own_sr;
      // if zone is inner - compute from neighbors params
      if( inner_overlap ){
         frameFrontOrigin = nbOrig;
         frameFrontOriginCoords = TNL::floor( ( nbOrig - globalOrig ) / own_sr );
         frameFrontDims = resolutionFactor * nbDims;
         frameFrontEnd = frameFrontOrigin + frameFrontDims;
         frameOrientation = -1;

         frameBackOrigin = frameFrontOrigin + bufferWidth * unitVect;
         //frameBackSize = frameFrontDims * nb_sr - 2 * bufferWidth * unitVect;
         frameBackSize = frameFrontDims * own_sr - 2 * bufferWidth * unitVect;
         frameBackDims = frameFrontDims - 2;
      }
      // if zone is outer - compute from local params
      else if( outer_overlap ) {
         frameFrontOrigin = ownOrig;
         frameFrontOriginCoords = 0;
         frameFrontDims = ownDims;
         frameFrontEnd = frameFrontOrigin + frameFrontDims;
         frameOrientation = 1;

         frameBackOrigin = frameFrontOrigin - bufferWidth * unitVect;
         frameBackSize = frameFrontDims * own_sr + 2 * bufferWidth * unitVect;
         frameBackDims = frameFrontDims + 2;
      }
      else {
         assert( false && "initZones: Invalid overlap state: neither inner_overlap nor outer_overlap is true!" );
      }

      std::cout << "=== DEBUG ZONE INIT START ===\n";

      std::cout << "own_sr: " << own_sr << "\n";
      std::cout << "bufferWidthFactorConst: " << bufferWidthFactorConst << "\n";
      std::cout << "bufferWidth: " << bufferWidth << "\n";

      std::cout << "inner_overlap: " << inner_overlap << "\n";
      std::cout << "outer_overlap: " << outer_overlap << "\n";

      std::cout << "nbOrig: " << nbOrig << "\n";
      std::cout << "globalOrig: " << globalOrig << "\n";
      std::cout << "nbDims: " << nbDims << "\n";
      std::cout << "ownDims: " << ownDims << "\n";
      std::cout << "resolutionFactor: " << resolutionFactor << "\n";

      std::cout << "unitVect: " << unitVect << "\n";

      std::cout << "frameFrontOrigin: " << frameFrontOrigin << "\n";
      std::cout << "frameFrontOriginCoords: " << frameFrontOriginCoords << "\n";
      std::cout << "frameFrontDims: " << frameFrontDims << "\n";
      std::cout << "frameFrontEnd: " << frameFrontEnd << "\n";
      std::cout << "frameOrientation: " << frameOrientation << "\n";

      std::cout << "frameBackOrigin: " << frameBackOrigin << "\n";
      std::cout << "frameBackSize: " << frameBackSize << "\n";

      std::cout << "=== DEBUG ZONE INIT END ===\n";

      zone.setNumberOfParticlesPerCell( maxPtcsPerCell );
      //TODO: Alternatively use the two concentric frameFront and frameBack
      std::cout << "INIT: assignCellsFrame:\n"
          << "  frameFrontOrigin: " << frameFrontOrigin << "\n"
          << "  frameFrontDims:   " << frameFrontDims << "\n"
          << "  frameOrientation: " << frameOrientation << "\n"
          << "  frameWidth:       " << frameWidth << "\n"
          << "  computed orient*width: " << (frameOrientation * frameWidth) << "\n"
          << "  ownDimsWithOverlap: " << ownDimsWithOverlap << std::endl;

      //zone.assignCellsFrame( frameFrontOrigin, frameFrontDims, frameOrientation * frameWidth, ownDimsWithOverlap );
      zone.assignCellsFrame( frameFrontOriginCoords, frameFrontDims, frameOrientation * frameWidth, ownDimsWithOverlap );
      //zone.assignCellsFrame( frameFrontOriginCoords, (-1) * frameFrontDims, frameOrientation * frameWidth, ownDimsWithOverlap ); //FIXME: The zone is extruded in different dirrection then massnodes
      const std::string outpufilename = "results/zone_debug" + std::to_string(own_sr) + ".vtk";
      zone.saveZoneToVTK( outpufilename,
                    frameFrontDims,
                    frameFrontOrigin,
                    ownParticles->getSearchRadius() );

      // Set size of multi-resoluton algorithm arrays (NOTE: Consider setSize function)
      const IndexType n_alloc = this->getNumberOfAllocatedParticles();
      particlesToFluid.setSize( n_alloc );
      particlesToRemove.setSize( n_alloc );
      particlesToBuffer.setSize( n_alloc );
      retypeMarker.setSize( n_alloc );

      std::cout << "========================================================================" << std::endl;
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
            IndexVectorType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  end[ pd ] = planeNodeCount( d, pd );

            // Strides for linearisation (perp axes only, reverse order)
            IndexVectorType stride = 0;
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

            auto generate = [=] __cuda_callable__ ( const IndexVectorType idx ) mutable
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
      TNL::Containers::Array< VectorType, TNL::Devices::Host > excluded( 2 );
      if( outer_overlap ){
         excluded[ 0 ] = {  -1.f, 0.f };   // exclude +x face
         excluded[ 1 ] = {  0.f, 1.f };  // exclude -y face
      }
      if( inner_overlap ){
         excluded[ 0 ] = {  1.f, 0.f };   // exclude +x face
         excluded[ 1 ] = {  0.f, -1.f };  // exclude -y face
      }
      initMassNodesWithExcludedNormals( modelParams, refinemnetFactor, excluded );
   }

   template< typename ModelParams >
   void
   initMassNodesWithExcludedNormals( ModelParams& modelParams,
                  const RealType refinementFactor,
                  const TNL::Containers::Array< VectorType, TNL::Devices::Host >& excludedNormals = {} )
   {
      const RealType local_dp = refinementFactor * modelParams.dp;
      const RealType eps      = local_dp * 1e-4f;

      // Returns true if this face should be skipped
      auto isExcluded = [&]( const VectorType& normal ) -> bool {
         for( int k = 0; k < excludedNormals.getSize(); k++ ) {
            bool match = true;
            for( int d = 0; d < VectorType::getSize(); d++ )
               if( std::fabs( normal[ d ] - excludedNormals[ k ][ d ] ) > eps )
                  { match = false; break; }
            if( match ) return true;
         }
         return false;
      };

      // Build the normal for a given (axis, sign) so both lambdas use the same logic
      auto faceNormal = [&]( int d, int sign ) -> VectorType {
         VectorType normal = 0.f;
         normal[ d ] = ( sign < 0 ) ? 1.f : -1.f;
         if( inner_overlap ) normal[ d ] *= -1.f;
         return normal;
      };
      //-------

      auto planeNodeCount = [&]( int faceAxis, int perpAxis ) -> IndexType {
         RealType extent = frameBackSize[ perpAxis ];
         for( int d = 0; d < faceAxis; d++ ) {
            if( d == perpAxis ) continue;
            // Only active (non-excluded) faces claim corner strips
            if( !isExcluded( faceNormal( d, -1 ) ) || !isExcluded( faceNormal( d, +1 ) ) )
               extent -= 2.f * local_dp;
         }
         return static_cast< IndexType >( TNL::max( 0.f, extent ) / local_dp );
      };

      auto perpOriginOffset = [&]( int faceAxis, int perpAxis ) -> RealType {
         RealType offset = 0.f;
         for( int d = 0; d < faceAxis; d++ ) {
            if( d == perpAxis ) continue;
            if( !isExcluded( faceNormal( d, -1 ) ) || !isExcluded( faceNormal( d, +1 ) ) )
               offset += local_dp;
         }
         return offset;
      };

      // Count total nodes across all 2*dim faces
      IndexType n_massNodes = 0;
      for( int d = 0; d < VectorType::getSize(); d++ ) {
         for( int sign : { -1, +1 } ) {
            if( isExcluded( faceNormal( d, sign ) ) ) continue;  // ← only change
            IndexType faceNodes = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d ) faceNodes *= planeNodeCount( d, pd );
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
             if( isExcluded( faceNormal( d, sign ) ) ) continue;

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
            IndexVectorType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < VectorType::getSize(); pd++ )
               if( pd != d )
                  end[ pd ] = planeNodeCount( d, pd );

            // Strides for linearisation (perp axes only, reverse order)
            IndexVectorType stride = 0;
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

            auto generate = [=] __cuda_callable__ ( const IndexVectorType idx ) mutable
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
      Algorithms::parallelFor< DeviceType >( 0, n_buffer, moveBufferParticles );

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
      const VectorType frameFrontOrigin = this->frameFrontOrigin;
      const VectorType frameBackOrigin = this->frameBackOrigin;
      const VectorType frameBackSize = this->frameBackSize;
      const IndexVectorType frameBackDims = this->frameBackDims;
      const IndexVectorType frameFrontDims = this->frameFrontDims;
      const IndexVectorType frameFrontOriginCoords = this->frameFrontOriginCoords;
      const bool inner_overlap = this->inner_overlap;
      const bool outer_overlap = this->outer_overlap;
      std::cout << "===================================================" << std::endl;
      std::cout << "sr: " << sr << "\n";
      std::cout << "inv_sr: " << inv_sr << "\n";

      std::cout << "frameFrontOrigin: " << frameFrontOrigin << "\n";
      std::cout << "frameBackOrigin: " << frameBackOrigin << "\n";
      std::cout << "frameBackSize: " << frameBackSize << "\n";
      std::cout << "frameBackDims: " << frameBackDims << "\n";

      std::cout << "frameFrontDims: " << frameFrontDims << "\n";
      std::cout << "frameFrontOriginCoords: " << frameFrontOriginCoords << "\n";

      std::cout << "inner_overlap: " << inner_overlap << "\n";
      std::cout << "outer_overlap: " << outer_overlap << "\n";
      std::cout << "******** \n";
      std::cout << "particles_orig: " << this->getParticles()->getGridOrigin() << "\n";
      std::cout << "particles_orig_with_overlap: " << this->getParticles()->getGridOriginWithOverlap() << "\n";
      std::cout << "particles_ref_orig: " << this->getParticles()->getGridReferentialOrigin() << "\n";
      std::cout << "particles_dim: " << this->getParticles()->getGridDimensions() << "\n";
      std::cout << "particles_dim_with_overlap: " <<  this->getParticles()->getGridDimensionsWithOverlap() << "\n";
      std::cout << "===================================================" << std::endl;


      // Helper lambda: convert physical position to own grid integer coords.
      auto toGridCoord = [=] __cuda_callable__ ( const VectorType& r ) -> IndexVectorType
      {
         IndexVectorType gc;
         gc = TNL::floor( ( r - frameFrontOrigin ) * inv_sr ); //TODO: Static cast?
         return gc;
      };

      // Retype to fluid if:
      // - outer zone + is inside
      // - inner zone + is outside
      auto identifyRetype = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const IndexVectorType gc = toGridCoord( view_r_buffer[ i ] );
         //bool inside = isInsideBox( gc, frameFrontOriginCoords, frameFrontDims );
         bool inside = isInsideBox( gc, frameFrontDims );

         if( inside && outer_overlap ) {
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

      // // Remove invalid buffe particles
      // // - outer zone - is outside the ( subdomain box + zone width )
      // // - inner zone - is insde the ( subdomain box - zone width )
      // // TODO (?): If we would use buffer width corresponding to search radius, we can use grid to compare positions
      // auto identifyRemove = [=] __cuda_callable__ ( IndexType i ) mutable
      // {
      //    //const IndexVectorType gc =
      //    bool inside = isInsideBox( view_r_buffer[ i ], frameBackOrigin, frameBackSize );
      //    //DEBUG
      //    if( view_r_buffer[ i ][1] >= 0.3 )
      //       printf( "[ x:%f, y:%f,  inside: %d ] ", view_r_buffer[ i ][0], view_r_buffer[ i ][1], inside );
      //    //\DEBUG
      //    if( inside && inner_overlap ) {
      //       retypeMarker_view[ i ] = 2;
      //       return 1;
      //    }
      //    else if( !inside && outer_overlap ) {
      //       retypeMarker_view[ i ] = 2;
      //       return 1;
      //    }
      //    return 0;
      // };
      // numberOfPtcsToRemove = Algorithms::reduce< DeviceType >( 0, n_buffer, identifyRemove );

      // Remove invalid buffe particles
      // - outer zone - is outside the ( subdomain box + zone width )
      // - inner zone - is insde the ( subdomain box - zone width )
      // TODO (?): If we would use buffer width corresponding to search radius, we can use grid to compare positions
      auto identifyRemove = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const IndexVectorType gc = TNL::floor( ( view_r_buffer[ i ] - frameBackOrigin ) * inv_sr ); //for outer
         bool inside = isInsideBox( gc, frameBackDims );
         //DEBUG
         if( view_r_buffer[ i ][1] >= 0.3 )
            printf( "[ x:%f, y:%f,  inside: %d ] ", view_r_buffer[ i ][0], view_r_buffer[ i ][1], inside );
         //\DEBUG
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
      numberOfPtcsToRemove = Algorithms::reduce< DeviceType >( 0, n_buffer, identifyRemove );
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
      //const IndexVectorType ownOrigCoords = fluid->getParticles()->getGridOriginGlobalCoords();
      //const VectorType refOrig = fluid->getParticles()->getGridReferentialOrigin();
      //const VectorType ownOrig = fluid->getParticles()->getGridOrigin();

      const IndexType numberOfBufferParticles = this->getParticles()->getNumberOfParticles();
      const IndexVectorType frameFrontDims = this->frameFrontDims;
      const VectorType frameFrontOrigin = this->frameFrontOrigin;
      const IndexVectorType frameFrontOriginCoords = this->frameFrontOriginCoords;
      const bool inner_overlap = this->inner_overlap;
      const bool outer_overlap = this->outer_overlap;
      const RealType sr = this->getParticles()->getSearchRadius();
      const RealType inv_sr = 1.f / sr;

      std::cout << "___________________________________________________________" << std::endl;
      std::cout << "numberOfZoneParticles: " << numberOfZoneParticles << "\n";
      std::cout << "numberOfBufferParticles: " << numberOfBufferParticles << "\n";

      //std::cout << "refOrig: " << refOrigin << "\n";
      //std::cout << "ownOrig: " << refOrigin << "\n";
      std::cout << "frameFrontDims: " << frameFrontDims << "\n";
      std::cout << "frameFrontOrigin: " << frameFrontOrigin << "\n";
      std::cout << "frameFrontOriginCoords: " << frameFrontOriginCoords << "\n";

      std::cout << "inner_overlap: " << inner_overlap << "\n";
      std::cout << "outer_overlap: " << outer_overlap << "\n";

      std::cout << "sr: " << sr << "\n";
      std::cout << "inv_sr: " << inv_sr << "\n";
      std::cout << "___________________________________________________________" << std::endl;

      // Retype fluid to buffer:
      // - outer zone and is outside the frame
      // - inner zone and is inside the frame
      auto checkFluidParticles = [ = ] __cuda_callable__( int i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         const IndexVectorType gc = TNL::floor( ( r - frameFrontOrigin ) * inv_sr );
         //const bool inside = isInsideBox( gc, frameFrontOriginCoords, frameFrontDims );
         const bool inside = isInsideBox( gc, frameFrontDims );
         //DEBUG
         //if(i<5){
         //   printf( "[ x:%f, y:%f, gcx: %d, gcy: %d, inside: %d ] ", r[0], r[1], gc[0], gc[1], inside );
         //}
         //\DEBUG

         if( inside && inner_overlap ) {
            particlesToBuffer_view[ i ] = p;
            return 1;
         }
         else if( !inside && outer_overlap ) {
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

      /*
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
      */

      std::cout << "accumulateMasses start" << std::endl;
      accumulateMasses( fluid_neihgbor, modelParams, dt );
      std::cout << "moveBufferParticles start" << std::endl;
      moveBufferParticles( dt );
      std::cout << "sortBufferParticles start" << std::endl;
      sortBufferParticles();
      std::cout << "removeBufferParticles start" << std::endl;
      removeBufferParticles();
      std::cout << "convertBufferToFluid start" << std::endl;
      convertBufferToFluid( fluid_own );
      std::cout << "updateParticlesInZone start" << std::endl;
      zone.updateParticlesInZone( fluid_own->getParticles() );
      std::cout << "getFluidParticlesEneringTheBuffer start" << std::endl;
      getFluidParticlesEneringTheBuffer( fluid_own );
      std::cout << "convertFluidToBuffer start" << std::endl;
      convertFluidToBuffer( fluid_own );
      std::cout << "updateMassNodes start" << std::endl;
      updateMassNodes( modelParams, dt );
      std::cout << "massNodes.sort start" << std::endl;
      massNodes.sort();
      std::cout << "createBufferParticles start" << std::endl;
      createBufferParticles( modelParams );
      std::cout << "interpolateVariables start" << std::endl;
      interpolateVariables( fluid_neihgbor, modelParams );
      std::cout << "searchForNeighbors start" << std::endl;
      this->searchForNeighbors();
      std::cout << "Pipeline finished" << std::endl;


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

   void
   writeMassNodesToVTK( const std::string& outputFileName )
   {
      std::cout << "Writing mass nodes: " << outputFileName << std::endl;
      const int n_mn = massNodes.numberOfMassNodes;
      ParticlesType nodes;
      nodes.setSize( n_mn );
      nodes.setNumberOfParticles( n_mn );
      nodes.getPoints() = massNodes.points;

      using WriterType = TNL::ParticleSystem::Writers::VTKWriter< ParticlesType >;
      std::ofstream outputFileFluid ( outputFileName, std::ofstream::out );
      WriterType writer( outputFileFluid );
      writer.writeParticles( nodes );
      writer.template writeVector< typename MassNodes::VectorArrayType, RealType >( massNodes.normal, "Normal", n_mn, 0, 3 );
      std::cout << "Done." << std::endl;
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
   VectorType frameFrontOrigin;
   IndexVectorType frameFrontOriginCoords;
   IndexVectorType frameFrontDims;
   IndexVectorType frameFrontEnd;

   IndexVectorType frameBackDims;
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

