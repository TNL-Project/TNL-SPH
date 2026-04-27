#pragma once

#include <TNL/Algorithms/scan.h>
#include "GhostZone.h"

namespace TNL {
namespace ParticleSystem {

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCells( IndexVectorType firstPointIdx,
                                                         IndexVectorType zoneSizeInCells,
                                                         IndexVectorType gridSize )
{
   if constexpr( ParticleConfig::spaceDimension == 2 )
      this->numberOfCellsInZone = zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ];
   if constexpr( ParticleConfig::spaceDimension == 3 )
      this->numberOfCellsInZone = zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ] * zoneSizeInCells[ 2 ];

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();

   if constexpr( ParticleConfig::spaceDimension == 2 ) {
      auto init = [=] __cuda_callable__ ( const IndexVectorType i ) mutable
      {
         const GlobalIndexType idxLinearized = i[ 0 ] + i[ 1 ] * zoneSizeInCells[ 0 ];
         cellsInZone_view[ idxLinearized ] = CellIndexer::EvaluateCellIndex( firstPointIdx + i, gridSize );
      };
      const IndexVectorType begin = { 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, zoneSizeInCells, init );
   }

   if constexpr( ParticleConfig::spaceDimension == 3 ) {
      auto init = [=] __cuda_callable__ ( const IndexVectorType i ) mutable
      {
         const GlobalIndexType idxLinearized = i[ 0 ] + i[ 1 ] * zoneSizeInCells[ 0 ] + i[ 2 ] * zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ];
         cellsInZone_view[ idxLinearized ] = CellIndexer::EvaluateCellIndex( firstPointIdx + i, gridSize );

         //cellsInZone_view[ idxLinearized ] = i[ 2 ] * gridSize[ 0 ] * gridSize[ 1 ] + i[ 1 ] * gridSize[ 0 ] + i[ 0 ];
         //cellsInZone_view[ idxLinearized ] = i[ 2 ] * gridSize[ 0 ] * gridSize[ 1 ] + i[ 0 ] * gridSize[ 1 ] + i[ 1 ];
      };
      const IndexVectorType begin = { 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, zoneSizeInCells, init );
   }
}

//TODO: Merge both assign functions together
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCells( const PointType firstPoint,
                                                         const PointType secondPoint,
                                                         IndexVectorType gridSize,
                                                         PointType gridOrigin,
                                                         RealType searchRadius )
{
   const PointType zoneSize = TNL::abs( secondPoint - firstPoint );
   const IndexVectorType zoneSizeInCells = TNL::ceil( zoneSize / searchRadius );
   const IndexVectorType firstPointIdx = (firstPoint - gridOrigin ) / searchRadius;

   if constexpr( ParticleConfig::spaceDimension == 2 )
      this->numberOfCellsInZone = zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ];

   if constexpr( ParticleConfig::spaceDimension == 3 )
      this->numberOfCellsInZone = zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ] * zoneSizeInCells[ 2 ];

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();

   if constexpr( ParticleConfig::spaceDimension == 2 ) {
      auto init = [=] __cuda_callable__ ( const IndexVectorType i ) mutable
      {
         const GlobalIndexType idxLinearized = i[ 0 ] + i[ 1 ] * zoneSizeInCells[ 0 ];
         cellsInZone_view[ idxLinearized ] = CellIndexer::EvaluateCellIndex( firstPointIdx + i, gridSize );
      };
      const IndexVectorType begin = { 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, zoneSizeInCells, init );
   }
   if constexpr( ParticleConfig::spaceDimension == 3 ) {
      auto init = [=] __cuda_callable__ ( const IndexVectorType i ) mutable
      {
         const GlobalIndexType idxLinearized = i[ 0 ] + i[ 1 ] * zoneSizeInCells[ 0 ] + i[ 2 ] * zoneSizeInCells[ 0 ] * zoneSizeInCells[ 1 ];
         cellsInZone_view[ idxLinearized ] = CellIndexer::EvaluateCellIndex( firstPointIdx + i, gridSize );
      };
      const IndexVectorType begin = { 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, zoneSizeInCells, init );
   }
}

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCellsFrame(
   const IndexVectorType frameFrontOrigin,
   const IndexVectorType frameFrontDims,
   const int             frameWidth,
   const IndexVectorType gridSize )
{
   constexpr int dim = ParticleConfig::spaceDimension;
   const int absWidth = TNL::abs( frameWidth );

   const IndexVectorType outerBegin = frameFrontOrigin - frameWidth;
   const IndexVectorType outerEnd   = frameFrontOrigin + frameFrontDims + frameWidth;
   const IndexVectorType outerDims  = outerEnd - outerBegin;

   // -----------------------------------------------------------------------
   // Cell count: 2*dim face slabs, each absWidth thick in normal direction
   // and full outerDims in all perp directions. Corners/edges are counted
   // multiple times intentionally.
   // -----------------------------------------------------------------------
   this->numberOfCellsInZone = 0;
   for( int d = 0; d < dim; d++ ) {
      GlobalIndexType faceArea = absWidth;  // thickness in normal direction
      for( int pd = 0; pd < dim; pd++ )
         if( pd != d )
            faceArea *= outerDims[ pd ];
      this->numberOfCellsInZone += 2 * faceArea;  // lo and hi face
   }

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( this->numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();

   // -----------------------------------------------------------------------
   // One parallelFor per face slab, writing to known offset in cellsInZone
   // -----------------------------------------------------------------------
   GlobalIndexType writeOffset = 0;

   for( int d = 0; d < dim; d++ ) {
      for( int sign : { -1, +1 } ) {

         // Slab spans full outerDims in perp directions, absWidth in dim d
         IndexVectorType slabBegin = outerBegin;
         IndexVectorType slabEnd   = outerEnd;

         if( sign == -1 ) {
            // lo face: outerBegin[d] .. outerBegin[d] + absWidth
            slabEnd[ d ] = outerBegin[ d ] + absWidth;
         } else {
            // hi face: outerEnd[d] - absWidth .. outerEnd[d]
            slabBegin[ d ] = outerEnd[ d ] - absWidth;
         }

         const IndexVectorType slabDims   = slabEnd - slabBegin;
         const GlobalIndexType slabOffset = writeOffset;

         if constexpr( dim == 2 ) {
            auto fill = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
               const GlobalIndexType lin = i[ 0 ] + i[ 1 ] * slabDims[ 0 ];
               cellsInZone_view[ slabOffset + lin ] =
                  CellIndexer::EvaluateCellIndex( slabBegin + i, gridSize );
            };
            Algorithms::parallelFor< DeviceType >( IndexVectorType{ 0, 0 }, slabDims, fill );
         }
         if constexpr( dim == 3 ) {
            auto fill = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
               const GlobalIndexType lin = i[ 0 ]
                                         + i[ 1 ] * slabDims[ 0 ]
                                         + i[ 2 ] * slabDims[ 0 ] * slabDims[ 1 ];
               cellsInZone_view[ slabOffset + lin ] =
                  CellIndexer::EvaluateCellIndex( slabBegin + i, gridSize );
            };
            Algorithms::parallelFor< DeviceType >( IndexVectorType{ 0, 0, 0 }, slabDims, fill );
         }

         GlobalIndexType slabVolume = 1;
         for( int pd = 0; pd < dim; pd++ )
            slabVolume *= slabDims[ pd ];
         writeOffset += slabVolume;
      }
   }

   TNL_ASSERT_EQ( writeOffset, this->numberOfCellsInZone,
      "Shell slab decomposition wrote a different count than geometric formula predicted." );
}
*/

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCellsFrame(
   const IndexVectorType frameFrontOrigin,
   const IndexVectorType frameFrontDims,
   const int             frameWidth,
   const IndexVectorType gridSize )
{
   constexpr int dim = ParticleConfig::spaceDimension;
   const int absWidth = TNL::abs( frameWidth );

   const IndexVectorType outerBegin = frameFrontOrigin - frameWidth;
   const IndexVectorType outerEnd   = frameFrontOrigin + frameFrontDims + frameWidth;

   // -----------------------------------------------------------------------
   // Pre-compute clamped slab begin/end/dims and volumes on CPU,
   // accumulating the total cell count before any allocation.
   // -----------------------------------------------------------------------
   struct SlabInfo {
      IndexVectorType begin, dims;
      GlobalIndexType volume;
   };
   TNL::Containers::Array< SlabInfo, TNL::Devices::Host > slabs( 2 * dim );

   this->numberOfCellsInZone = 0;
   for( int d = 0; d < dim; d++ ) {
      for( int si = 0; si < 2; si++ ) {
         const int sign = ( si == 0 ) ? -1 : +1;

         IndexVectorType slabBegin = outerBegin;
         IndexVectorType slabEnd   = outerEnd;

         if( sign == -1 )
            slabEnd  [ d ] = outerBegin[ d ] + absWidth;
         else
            slabBegin[ d ] = outerEnd  [ d ] - absWidth;

         // Clamp to domain [0, gridSize)
         for( int pd = 0; pd < dim; pd++ ) {
            slabBegin[ pd ] = TNL::max( slabBegin[ pd ], 0 );
            slabEnd  [ pd ] = TNL::min( slabEnd  [ pd ], gridSize[ pd ] );
         }

         GlobalIndexType volume = 1;
         IndexVectorType slabDims;
         for( int pd = 0; pd < dim; pd++ ) {
            slabDims[ pd ] = TNL::max( 0, slabEnd[ pd ] - slabBegin[ pd ] );
            volume *= slabDims[ pd ];
         }

         slabs[ 2 * d + si ] = { slabBegin, slabDims, volume };
         this->numberOfCellsInZone += volume;
      }
   }

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( this->numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();

   // -----------------------------------------------------------------------
   // Launch one parallelFor per slab with known write offset
   // -----------------------------------------------------------------------
   GlobalIndexType writeOffset = 0;
   for( int s = 0; s < 2 * dim; s++ ) {
      const auto [ slabBegin, slabDims, slabVolume ] = slabs[ s ];

      if( slabVolume == 0 ) {
         continue;  // slab fully outside domain
      }

      const GlobalIndexType slabOffset = writeOffset;

      if constexpr( dim == 2 ) {
         auto fill = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
            const GlobalIndexType lin = i[ 0 ] + i[ 1 ] * slabDims[ 0 ];
            cellsInZone_view[ slabOffset + lin ] =
               CellIndexer::EvaluateCellIndex( slabBegin + i, gridSize );
         };
         Algorithms::parallelFor< DeviceType >( IndexVectorType{ 0, 0 }, slabDims, fill );
      }
      if constexpr( dim == 3 ) {
         auto fill = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
            const GlobalIndexType lin = i[ 0 ]
                                      + i[ 1 ] * slabDims[ 0 ]
                                      + i[ 2 ] * slabDims[ 0 ] * slabDims[ 1 ];
            cellsInZone_view[ slabOffset + lin ] =
               CellIndexer::EvaluateCellIndex( slabBegin + i, gridSize );
         };
         Algorithms::parallelFor< DeviceType >( IndexVectorType{ 0, 0, 0 }, slabDims, fill );
      }

      writeOffset += slabVolume;
   }

   TNL_ASSERT_EQ( writeOffset, this->numberOfCellsInZone,
      "Clamped shell slab decomposition wrote a different count than expected." );
}
*/

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCellsFrame(
   const IndexVectorType frameFrontOrigin,
   const IndexVectorType frameFrontDims,
   const int             frameWidth,
   const IndexVectorType gridSize )
{
   constexpr int dim = ParticleConfig::spaceDimension;
   const int w = TNL::abs( frameWidth );

   IndexVectorType outerBegin, outerEnd;
   IndexVectorType innerBegin, innerEnd;

   if( frameWidth >= 0 ) {
      for( int d = 0; d < dim; d++ ) {
         outerBegin[ d ] = frameFrontOrigin[ d ] - w;
         outerEnd  [ d ] = frameFrontOrigin[ d ] + frameFrontDims[ d ] + w;
         innerBegin[ d ] = frameFrontOrigin[ d ];
         innerEnd  [ d ] = frameFrontOrigin[ d ] + frameFrontDims[ d ];
      }
   } else {
      for( int d = 0; d < dim; d++ ) {
         outerBegin[ d ] = frameFrontOrigin[ d ];
         outerEnd  [ d ] = frameFrontOrigin[ d ] + frameFrontDims[ d ];
         innerBegin[ d ] = frameFrontOrigin[ d ] + w;
         innerEnd  [ d ] = frameFrontOrigin[ d ] + frameFrontDims[ d ] - w;
      }
   }

   // Clamp outer and inner boxes to [0, gridSize)
   for( int d = 0; d < dim; d++ ) {
      outerBegin[ d ] = TNL::max( outerBegin[ d ], 0 );
      outerEnd  [ d ] = TNL::min( outerEnd  [ d ], gridSize[ d ] );
      innerBegin[ d ] = TNL::max( innerBegin[ d ], 0 );
      innerEnd  [ d ] = TNL::min( innerEnd  [ d ], gridSize[ d ] );
   }

   // Outer box iteration dims and flat volume
   IndexVectorType outerDims;
   GlobalIndexType outerVolume = 1;
   for( int d = 0; d < dim; d++ ) {
      outerDims[ d ] = TNL::max( 0, outerEnd[ d ] - outerBegin[ d ] );
      outerVolume   *= outerDims[ d ];
   }

   // Inner box validity — if clamping collapsed it, no cell is excluded
   bool innerValid = true;
   for( int d = 0; d < dim; d++ )
      if( innerBegin[ d ] >= innerEnd[ d ] )
         innerValid = false;

   // isShell: 1 if this outer-box-local flat index is a shell cell, 0 otherwise
   TNL::Containers::Array< int, DeviceType, GlobalIndexType > isShell( outerVolume );

   auto isShell_view   = isShell.getView();
   const auto ob       = outerBegin;
   const auto od       = outerDims;
   const auto ib       = innerBegin;
   const auto ie       = innerEnd;
   const bool iv       = innerValid;

   // Pass 1: mark shell cells
   if constexpr( dim == 2 ) {
      auto mark = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
         const IndexVectorType c = ob + i;
         bool inner = iv;
         if( iv )
            for( int d = 0; d < dim; d++ )
               if( c[ d ] < ib[ d ] || c[ d ] >= ie[ d ] )
                  { inner = false; break; }
         const GlobalIndexType lin = i[ 0 ] + i[ 1 ] * od[ 0 ];
         isShell_view[ lin ] = inner ? 0 : 1;
      };
      Algorithms::parallelFor< DeviceType >(
            IndexVectorType{ 0, 0 }, outerDims, mark );
   } else {
      auto mark = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
         const IndexVectorType c = ob + i;
         bool inner = iv;
         if( iv )
            for( int d = 0; d < dim; d++ )
               if( c[ d ] < ib[ d ] || c[ d ] >= ie[ d ] )
                  { inner = false; break; }
         const GlobalIndexType lin = i[ 0 ]
                                   + i[ 1 ] * od[ 0 ]
                                   + i[ 2 ] * od[ 0 ] * od[ 1 ];
         isShell_view[ lin ] = inner ? 0 : 1;
      };
      Algorithms::parallelFor< DeviceType >(
            IndexVectorType{ 0, 0, 0 }, outerDims, mark );
   }

   // Pass 2: exclusive prefix sum → write indices
   TNL::Containers::Array< GlobalIndexType, DeviceType, GlobalIndexType >
         writeIdx( outerVolume );
   //TNL::Algorithms::exclusive_scan( isShell, writeIdx, TNL::Plus{}, 0 );
   TNL::Algorithms::inplaceExclusiveScan( isShell );

   // Total shell cells = last write index + last isShell value
   // Both are single elements — fetch to host
   this->numberOfCellsInZone =
         writeIdx.getElement( outerVolume - 1 )
       + isShell.getElement(  outerVolume - 1 );

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( this->numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();
   auto writeIdx_view    = writeIdx.getConstView();
   const auto gs         = gridSize;

   // Pass 3: scatter shell cell flat indices into cellsInZone
   if constexpr( dim == 2 ) {
      auto scatter = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
         const GlobalIndexType lin = i[ 0 ] + i[ 1 ] * od[ 0 ];
         if( isShell_view[ lin ] ) {
            const IndexVectorType c = ob + i;
            cellsInZone_view[ writeIdx_view[ lin ] ] =
               CellIndexer::EvaluateCellIndex( c, gs );
         }
      };
      Algorithms::parallelFor< DeviceType >(
            IndexVectorType{ 0, 0 }, outerDims, scatter );
   } else {
      auto scatter = [=] __cuda_callable__ ( const IndexVectorType i ) mutable {
         const GlobalIndexType lin = i[ 0 ]
                                   + i[ 1 ] * od[ 0 ]
                                   + i[ 2 ] * od[ 0 ] * od[ 1 ];
         if( isShell_view[ lin ] ) {
            const IndexVectorType c = ob + i;
            cellsInZone_view[ writeIdx_view[ lin ] ] =
               CellIndexer::EvaluateCellIndex( c, gs );
         }
      };
      Algorithms::parallelFor< DeviceType >(
            IndexVectorType{ 0, 0, 0 }, outerDims, scatter );
   }
}
*/

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCellsFrame(
   const IndexVectorType frameFrontOrigin,
   const IndexVectorType frameFrontDims,
   const int             frameWidth,
   const IndexVectorType gridSize )
{
   constexpr int dim = ParticleConfig::spaceDimension;
   const int w    = TNL::abs( frameWidth );
   const int sign = ( frameWidth >= 0 ) ? 1 : -1;  // +1 outward, -1 inward

   std::cout << "\n[assignCellsFrame] ========================================\n"
             << "  frameFrontOrigin : " << frameFrontOrigin << "\n"
             << "  frameFrontDims   : " << frameFrontDims   << "\n"
             << "  frameWidth       : " << frameWidth << "  (w=" << w << ", sign=" << sign << ")\n"
             << "  gridSize         : " << gridSize         << "\n"
             << "  spaceDimension   : " << dim              << "\n"
             << "============================================================\n";

   // -----------------------------------------------------------------------
   // Pass 1 — count total cells across all layers and faces on CPU.
   // One layer = one shell at offset k from the frame boundary.
   // -----------------------------------------------------------------------
   this->numberOfCellsInZone = 0;

   auto clampedFaceCount = [&]( int layer, int faceAxis, int sign, int perpAxis ) -> GlobalIndexType
   {
      // Physical slab extent on perpAxis for this face at this layer,
      // shrunk by lower-priority axes (corner ownership, same as initMassNodes).
      // Then additionally clamped to [0, gridSize[perpAxis]).
      //
      // The face at layer k is offset by k cells from the frame boundary:
      //   outward (sign=+1): frame expands by k on each side
      //   inward  (sign=-1): frame shrinks by k on each side
      const int expand = sign * layer;   // how many cells the frame has grown

      // ------ //NOTE: FIXME: My attampt to fix the corners
      // // Origin and end of the expanded/shrunk frame on perpAxis
      // GlobalIndexType faceOrigin = frameFrontOrigin[ perpAxis ] - expand;
      // GlobalIndexType faceEnd    = frameFrontOrigin[ perpAxis ] + frameFrontDims[ perpAxis ] + expand;

      // ---
      // For outward (sign=+1): perp extent expands by (layer+1) because the face
      // sits one cell AHEAD of the expansion, so perp must match that full extent.
      // For inward  (sign=-1): perp extent shrinks by layer (unchanged — works correctly).
      const int perpExpand = ( sign > 0 ) ? expand + 1 : expand;

      GlobalIndexType faceOrigin = frameFrontOrigin[ perpAxis ] - perpExpand;
      GlobalIndexType faceEnd    = frameFrontOrigin[ perpAxis ] + frameFrontDims[ perpAxis ] + perpExpand;


      // ----- NOTE: Attempt end

      // Corner ownership: lower-priority axes shrink by 1 on each end
      for( int d = 0; d < faceAxis; d++ )
         if( d != perpAxis ) {
            faceOrigin++;
            faceEnd--;
         }

      // Clamp to domain
      faceOrigin = TNL::max( faceOrigin, 0 );
      faceEnd    = TNL::min( faceEnd,    gridSize[ perpAxis ] );

      return static_cast< GlobalIndexType >( TNL::max( 0, faceEnd - faceOrigin ) );
   };

   for( int layer = 0; layer < w; layer++ ) {
      for( int d = 0; d < dim; d++ ) {
         for( int s : { -1, +1 } ) {
            GlobalIndexType faceNodes = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d )
                  faceNodes *= clampedFaceCount( layer, d, sign, pd );
            this->numberOfCellsInZone += faceNodes;
         }
      }
   }

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( this->numberOfCellsInZone * numberOfParticlesPerCell );

   auto cellsInZone_view = this->cellsInZone.getView();

   // -----------------------------------------------------------------------
   // Pass 2 — fill cells layer by layer, face by face.
   // Mirrors initMassNodes exactly: parallelFor over the face slab,
   // linearise to a write index, store the flat cell index.
   // -----------------------------------------------------------------------
   GlobalIndexType offset = 0;

   for( int layer = 0; layer < w; layer++ ) {
      const int expand = sign * layer;
      std::cout << " &&& expand: " << expand << std::endl;

      for( int d = 0; d < dim; d++ ) {
         for( int s : { -1, +1 } ) {
            std::cout << " &&&&& d: " << d << " s: " << s << std::endl;

            // Cell coordinate on the interface axis for this face+layer
            const GlobalIndexType ifaceCoord = ( s < 0 )
                  ? frameFrontOrigin[ d ] - expand - ( sign > 0 ? 1 : 0 )
                  : frameFrontOrigin[ d ] + frameFrontDims[ d ] + expand - ( sign > 0 ? 0 : 1 );

            // parallelFor range: 1 on interface axis, clamped count on perp axes
            IndexVectorType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d )
                  end[ pd ] = clampedFaceCount( layer, d, sign, pd );

            // Strides for linearisation (same as initMassNodes)
            IndexVectorType stride = 0;
            {
               GlobalIndexType running = 1;
               for( int pd = dim - 1; pd >= 0; pd-- ) {
                  if( pd == d ) continue;
                  stride[ pd ] = running;
                  running     *= end[ pd ];
               }
            }

            // Origin on each perp axis (corner ownership + clamping)
            IndexVectorType perpOrigin = 0;
            perpOrigin[ d ] = ifaceCoord;
            for( int pd = 0; pd < dim; pd++ ) {
               if( pd == d ) continue;
                  // //NOTE: FIXME: My attampt to fix the corners
                  // GlobalIndexType o = frameFrontOrigin[ pd ] - expand;
                  const int perpExpand = ( sign > 0 ) ? expand + 1 : expand;
                  GlobalIndexType o = frameFrontOrigin[ pd ] - perpExpand;
               // Corner ownership: lower-priority axes push origin inward
               for( int d2 = 0; d2 < d; d2++ )
                  if( d2 != pd ) o++;
               // Clamp
               perpOrigin[ pd ] = TNL::max( o, 0 );
            }

            // Clamp ifaceCoord — skip face entirely if outside domain
            if( ifaceCoord < 0 || ifaceCoord >= gridSize[ d ] ) { //FIXME:
            //if( ifaceCoord < -1 || ifaceCoord >= ( gridSize[ d ] - 1 ) ) { //FIXME:
               // Still advance offset by the face node count
               GlobalIndexType faceNodes = 1;
               for( int pd = 0; pd < dim; pd++ )
                  if( pd != d ) faceNodes *= end[ pd ];
               offset += faceNodes;
               continue;
            }

            const GlobalIndexType    faceOffset  = offset;
            const int          iAxis       = d;
            const GlobalIndexType    iCoord      = ifaceCoord;
            const IndexVectorType pOrigin  = perpOrigin;
            const IndexVectorType gs       = gridSize;

            // DEBUG PLAYGROUND
            std::cout
                << " ===> layer: " << layer
                << " ===> , faceOffset: " << faceOffset
                << " ===> , iAxis: " << iAxis
                << " ===> , iCoord: " << iCoord
                << " ===> , pOrigin: " << pOrigin
                << " ===> , end: " << end
                << " ===> , stride: " << stride
                << " ===> , gs: " << gs
                << std::endl;
            //\DEBUG PLAYGROUND

            auto fill = [=] __cuda_callable__ ( const IndexVectorType idx ) mutable
            {
               // Linearise over perp axes
               GlobalIndexType i = faceOffset;
               for( int pd = 0; pd < dim; pd++ )
                  i += idx[ pd ] * stride[ pd ];

               // Absolute cell coords
               IndexVectorType c = pOrigin;
               c[ iAxis ] = iCoord;
               for( int pd = 0; pd < dim; pd++ )
                  if( pd != iAxis )
                     c[ pd ] += idx[ pd ];

               // Skip if outside domain (clamped count should prevent this,
               // but guard against rounding at corners)
               for( int pd = 0; pd < dim; pd++ )
                  if( c[ pd ] < 0 || c[ pd ] >= gs[ pd ] ) return; //FIXME
                  //if( c[ pd ] < -1 || c[ pd ] >= ( gs[ pd ] - 1 ) ) return; //FIXME

               printf(" [ %d, %d ] ", c[0], c[1]);
               cellsInZone_view[ i ] = CellIndexer::EvaluateCellIndex( c, gs );
            };

            if constexpr( dim == 2 )
               Algorithms::parallelFor< DeviceType >( IndexVectorType{0,0}, end, fill );
            else
               Algorithms::parallelFor< DeviceType >( IndexVectorType{0,0,0}, end, fill );
            std::cout << std::endl;

            GlobalIndexType faceNodes = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d ) faceNodes *= end[ pd ];
            offset += faceNodes;
         }
      }
   }
}

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::assignCellsFrame(
   const IndexVectorType frameFrontOrigin,
   const IndexVectorType frameFrontDims,
   const int             frameWidth,
   const IndexVectorType gridSize )
{
   constexpr int dim = ParticleConfig::spaceDimension;
   const int w    = TNL::abs( frameWidth );
   const int sign = ( frameWidth >= 0 ) ? 1 : -1;  // +1 outward, -1 inward

   std::cout << "\n[assignCellsFrame] ========================================\n"
             << "  frameFrontOrigin : " << frameFrontOrigin << "\n"
             << "  frameFrontDims   : " << frameFrontDims   << "\n"
             << "  frameWidth       : " << frameWidth << "  (w=" << w << ", sign=" << sign << ")\n"
             << "  gridSize         : " << gridSize         << "\n"
             << "  spaceDimension   : " << dim              << "\n"
             << "============================================================\n";

   // -----------------------------------------------------------------------
   // Pass 1 — count total cells across all layers and faces on CPU.
   // One layer = one shell at offset k from the frame boundary.
   // -----------------------------------------------------------------------
   this->numberOfCellsInZone = 0;

   auto clampedFaceCount = [&]( int layer, int faceAxis, int sign, int perpAxis ) -> GlobalIndexType
   {
      const int expand = sign * layer;

      GlobalIndexType faceOrigin = frameFrontOrigin[ perpAxis ] - expand;
      GlobalIndexType faceEnd    = frameFrontOrigin[ perpAxis ] + frameFrontDims[ perpAxis ] + expand;

      for( int d = 0; d < faceAxis; d++ )
         if( d != perpAxis ) {
            faceOrigin++;
            faceEnd--;
         }

      faceOrigin = TNL::max( faceOrigin, 0 );
      faceEnd    = TNL::min( faceEnd,    gridSize[ perpAxis ] );

      return static_cast< GlobalIndexType >( TNL::max( 0, faceEnd - faceOrigin ) );
   };

   std::cout << "\n[Pass 1] Counting cells per layer/face\n";
   for( int layer = 0; layer < w; layer++ ) {
      const int expand = sign * layer;
      GlobalIndexType layerTotal = 0;
      std::cout << "  [layer " << layer << "] expand=" << expand << "\n";

      for( int d = 0; d < dim; d++ ) {
         std::cout << "    [axis d=" << d << "]\n";

         for( int s : { -1, +1 } ) {
            GlobalIndexType faceNodes = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d )
                  faceNodes *= clampedFaceCount( layer, d, sign, pd );
            this->numberOfCellsInZone += faceNodes;
            layerTotal += faceNodes;
            std::cout << "      [face s=" << std::setw(2) << s << "]  faceNodes=" << faceNodes
                      << "  (runningTotal=" << this->numberOfCellsInZone << ")\n";
         }
      }
      std::cout << "  [layer " << layer << "] subtotal=" << layerTotal << "\n";
   }
   std::cout << "[Pass 1] DONE — numberOfCellsInZone=" << this->numberOfCellsInZone << "\n";

   cellsInZone.resize( this->numberOfCellsInZone );
   numberOfParticlesInCell.resize( this->numberOfCellsInZone );
   particlesInZone.resize( this->numberOfCellsInZone * numberOfParticlesPerCell );

   std::cout << "\n[Resize] cellsInZone         = " << this->numberOfCellsInZone << "\n"
             << "[Resize] numberOfParticlesInCell = " << this->numberOfCellsInZone << "\n"
             << "[Resize] particlesInZone         = " << this->numberOfCellsInZone * numberOfParticlesPerCell
             << "  (perCell=" << numberOfParticlesPerCell << ")\n";

   auto cellsInZone_view = this->cellsInZone.getView();

   // -----------------------------------------------------------------------
   // Pass 2 — fill cells layer by layer, face by face.
   // -----------------------------------------------------------------------
   GlobalIndexType offset = 0;

   std::cout << "\n[Pass 2] Filling cells\n";
   for( int layer = 0; layer < w; layer++ ) {
      const int expand = sign * layer;
      std::cout << "\n  [layer " << layer << "] expand=" << expand
                << "  offset-at-layer-start=" << offset << "\n";

      for( int d = 0; d < dim; d++ ) {
         std::cout << "    [axis d=" << d << "]\n";

         for( int s : { -1, +1 } ) {
            const GlobalIndexType ifaceCoord = ( s < 0 )
                  ? frameFrontOrigin[ d ] - expand - ( sign > 0 ? 1 : 0 )
                  : frameFrontOrigin[ d ] + frameFrontDims[ d ] + expand - ( sign > 0 ? 0 : 1 );

            IndexVectorType begin = 0, end = 0;
            end[ d ] = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d )
                  end[ pd ] = clampedFaceCount( layer, d, sign, pd );

            IndexVectorType stride = 0;
            {
               GlobalIndexType running = 1;
               for( int pd = dim - 1; pd >= 0; pd-- ) {
                  if( pd == d ) continue;
                  stride[ pd ] = running;
                  running     *= end[ pd ];
               }
            }

            IndexVectorType perpOrigin = 0;
            perpOrigin[ d ] = ifaceCoord;
            for( int pd = 0; pd < dim; pd++ ) {
               if( pd == d ) continue;
               GlobalIndexType o = frameFrontOrigin[ pd ] - expand;
               for( int d2 = 0; d2 < d; d2++ )
                  if( d2 != pd ) o++;
               perpOrigin[ pd ] = TNL::max( o, 0 );
            }

            GlobalIndexType faceNodes = 1;
            for( int pd = 0; pd < dim; pd++ )
               if( pd != d ) faceNodes *= end[ pd ];

            const bool outOfDomain = ( ifaceCoord < 0 || ifaceCoord >= gridSize[ d ] );

            std::cout << "      [face s=" << std::setw(2) << s << "]"
                      << "  ifaceCoord=" << ifaceCoord
                      << "  perpOrigin=" << perpOrigin
                      << "  end="        << end
                      << "  stride="     << stride
                      << "  faceNodes="  << faceNodes
                      << "  offset="     << offset;

            if( outOfDomain ) {
               std::cout << "  => SKIPPED (ifaceCoord=" << ifaceCoord
                         << " outside [0," << gridSize[ d ] << "))\n";
               offset += faceNodes;
               continue;
            }

            std::cout << "  => FILLING cellsInZone[" << offset
                      << ".." << offset + faceNodes - 1 << "]\n";

            const GlobalIndexType    faceOffset  = offset;
            const int          iAxis       = d;
            const GlobalIndexType    iCoord      = ifaceCoord;
            const IndexVectorType pOrigin  = perpOrigin;
            const IndexVectorType gs       = gridSize;

            auto fill = [=] __cuda_callable__ ( const IndexVectorType idx ) mutable
            {
               GlobalIndexType i = faceOffset;
               for( int pd = 0; pd < dim; pd++ )
                  i += idx[ pd ] * stride[ pd ];

               IndexVectorType c = pOrigin;
               c[ iAxis ] = iCoord;
               for( int pd = 0; pd < dim; pd++ )
                  if( pd != iAxis )
                     c[ pd ] += idx[ pd ];

               for( int pd = 0; pd < dim; pd++ )
                  if( c[ pd ] < 0 || c[ pd ] >= gs[ pd ] ) return;

               cellsInZone_view[ i ] = CellIndexer::EvaluateCellIndex( c, gs );
            };

            if constexpr( dim == 2 )
               Algorithms::parallelFor< DeviceType >( IndexVectorType{0,0}, end, fill );
            else
               Algorithms::parallelFor< DeviceType >( IndexVectorType{0,0,0}, end, fill );

            offset += faceNodes;
         }
      }
      std::cout << "  [layer " << layer << "] offset-at-layer-end=" << offset << "\n";
   }

   std::cout << "\n[Pass 2] DONE — final offset=" << offset
             << "  (expected numberOfCellsInZone=" << this->numberOfCellsInZone << ")"
             << ( offset == this->numberOfCellsInZone ? "  [OK]\n" : "  [MISMATCH!]\n" );
   std::cout << "[assignCellsFrame] ==================== END ================\n\n";
}
*/


template< typename ParticleConfig, typename DeviceType >
template< typename Array >
void
ParticleZone< ParticleConfig, DeviceType >::assignCells( Array& inputCells )
{

}

template< typename ParticleConfig, typename DeviceType >
const typename ParticleZone< ParticleConfig, DeviceType >::IndexArrayType&
ParticleZone< ParticleConfig, DeviceType >::getCellsInZone() const
{
   return cellsInZone;
}

template< typename ParticleConfig, typename DeviceType >
typename ParticleZone< ParticleConfig, DeviceType >::IndexArrayType&
ParticleZone< ParticleConfig, DeviceType >::getCellsInZone()
{
   return cellsInZone;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticleZone< ParticleConfig, DeviceType >::IndexArrayType&
ParticleZone< ParticleConfig, DeviceType >::getParticlesInZone() const
{
   return particlesInZone;
}

template< typename ParticleConfig, typename DeviceType >
typename ParticleZone< ParticleConfig, DeviceType >::IndexArrayType&
ParticleZone< ParticleConfig, DeviceType >::getParticlesInZone()
{
   return particlesInZone;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticleZone< ParticleConfig, DeviceType >::GlobalIndexType
ParticleZone< ParticleConfig, DeviceType >::getNumberOfParticles() const
{
   return numberOfParticlesInZone;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticleZone< ParticleConfig, DeviceType >::GlobalIndexType
ParticleZone< ParticleConfig, DeviceType >::getNumberOfCells() const
{
   return numberOfCellsInZone;
}

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::setNumberOfParticlesPerCell( const GlobalIndexType numberOfParticlesPerCell )
{
   this->numberOfParticlesPerCell = numberOfParticlesPerCell;
}

template< typename ParticleConfig, typename DeviceType >
const typename ParticleZone< ParticleConfig, DeviceType >::GlobalIndexType
ParticleZone< ParticleConfig, DeviceType >::getNumberOfParticlesPerCell() const
{
   return this->numberOfParticlesPerCell;
}

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::resetParticles()
{
   numberOfParticlesInZone = 0;
   numberOfParticlesInCell = 0;
   particlesInZone = 0;
}

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::resetZoneCells()
{
   numberOfParticlesInZone = 0;
   particlesInZone = 0;
   cellsInZone = 0;
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig, DeviceType >::collectNumbersOfParticlesInCells( const ParticlesPointer& particles )
{
   const auto firstLastParticle_view = particles->getCellFirstLastParticleList().getConstView();
   const auto cellsInZone_view = this->cellsInZone.getConstView();
   auto numberOfParticlesInCell_view = this->numberOfParticlesInCell.getView();

   auto collectParticlesCounts = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType cell = cellsInZone_view[ i ];
      const PairIndexType firstAndLastParticleInCell = firstLastParticle_view[ cell ];

      if( firstAndLastParticleInCell[ 0 ] != INT_MAX )
      {
         numberOfParticlesInCell_view[ i ] = firstAndLastParticleInCell[ 1 ] - firstAndLastParticleInCell[ 0 ] + 1;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticlesCounts );
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig, DeviceType >::buildParticleList( const ParticlesPointer& particles )
{

   Algorithms::inplaceExclusiveScan( this->numberOfParticlesInCell );
   this->numberOfParticlesInZone = this->numberOfParticlesInCell.getElement( numberOfCellsInZone - 1 ); //without last cell!

   const auto firstLastCellParticle_view = particles->getCellFirstLastParticleList().getConstView();
   const auto cellsInZone_view = this->cellsInZone.getConstView();
   const auto numberOfParticlesInCell_view = this->numberOfParticlesInCell.getConstView();
   auto particlesInZone_view = this->particlesInZone.getView();

   auto collectParticles = [=] __cuda_callable__ ( int i ) mutable //TODO: This i is cell index, rename it
   {
      const GlobalIndexType cell = cellsInZone_view[ i ];
      const PairIndexType firstLastParticle = firstLastCellParticle_view[ cell ];

      if( firstLastParticle[ 0 ] != INT_MAX )
      {
         const GlobalIndexType numberOfPtcsInThisCell = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;
         const GlobalIndexType particleListCellPrefix = numberOfParticlesInCell_view[ i ];

         for( int j = 0; j < numberOfPtcsInThisCell; j++ )
         {
            particlesInZone_view[ particleListCellPrefix + j ] = firstLastParticle[ 0 ] + j;
         }
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfCellsInZone, collectParticles );

   //Add the particle from last cell
   const PairIndexType firstLastParticleLastCell = firstLastCellParticle_view.getElement(
         cellsInZone_view.getElement( this->numberOfCellsInZone - 1 ) );
   if( firstLastParticleLastCell[ 0 ] != INT_MAX )
   {
      this->numberOfParticlesInZone += ( firstLastParticleLastCell[ 1 ] - firstLastParticleLastCell[ 0 ] + 1 );
   }

}


template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointer >
void
ParticleZone< ParticleConfig, DeviceType >::updateParticlesInZone( const ParticlesPointer& particles )
{
   this->resetParticles();
   this->collectNumbersOfParticlesInCells( particles );
   this->buildParticleList( particles );
}

template< typename ParticleConfig, typename DeviceType >
template< typename ParticlesPointer, typename TimeMeasurement >
void
ParticleZone< ParticleConfig, DeviceType >::updateParticlesInZone( const ParticlesPointer& particles, TimeMeasurement& timeMeasurement )
{
   timeMeasurement.start( "zone-reset" );
   this->resetParticles();
   timeMeasurement.stop( "zone-reset" );
   timeMeasurement.start( "zone-collect" );
   this->collectNumbersOfParticlesInCells( particles );
   timeMeasurement.stop( "zone-collect" );
   timeMeasurement.start( "zone-build" );
   this->buildParticleList( particles );
   timeMeasurement.stop( "zone-build" );
}

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::writeProlog( TNL::Logger& logger ) const noexcept
{
   logger.writeParameter( "Particle zone information:", "" );
   logger.writeParameter( "Number of particles per cell:", numberOfParticlesPerCell, 1 );
   logger.writeParameter( "Number of cells in zone:", numberOfCellsInZone, 1 );
}

/*
template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::saveZoneToVTK(
   const std::string&    filename,
   const IndexVectorType gridSize,
   const PointType      gridOrigin,
   const RealType        searchRadius ) const
{
   constexpr int dim = ParticleConfig::spaceDimension;

   // Copy cell indices to host for writing
   TNL::Containers::Array< GlobalIndexType, TNL::Devices::Host, GlobalIndexType >
         cellsHost( this->numberOfCellsInZone );
   cellsHost = this->cellsInZone;  // device → host copy

   std::ofstream f( filename );
   if( !f.is_open() )
      throw std::runtime_error( "saveZoneToVTK: cannot open " + filename );

   const GlobalIndexType n = this->numberOfCellsInZone;

   // Each cell is written as a voxel (VTK_VOXEL = 11 in 3D) or
   // a pixel (VTK_PIXEL = 8 in 2D). A cell at grid coords (ix, iy)
   // has its corner at gridOrigin + {ix, iy} * searchRadius.

   const int pointsPerCell = ( dim == 2 ) ? 4 : 8;

   f << "# vtk DataFile Version 3.0\n";
   f << "ParticleZone frame\n";
   f << "ASCII\n";
   f << "DATASET UNSTRUCTURED_GRID\n";

   // --- Points ---
   f << "POINTS " << n * pointsPerCell << " float\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ ) {
      // Decode flat cell index back to grid coords
      GlobalIndexType idx = cellsHost[ ci ];
      IndexVectorType gc;
      if constexpr( dim == 2 ) {
         gc[ 0 ] = idx % gridSize[ 0 ];
         gc[ 1 ] = idx / gridSize[ 0 ];
      } else {
         gc[ 0 ] = idx % gridSize[ 0 ];
         gc[ 1 ] = ( idx / gridSize[ 0 ] ) % gridSize[ 1 ];
         gc[ 2 ] = idx / ( gridSize[ 0 ] * gridSize[ 1 ] );
      }

      // Cell corner in physical coords
      const float ox = gridOrigin[ 0 ] + gc[ 0 ] * searchRadius;
      const float oy = gridOrigin[ 1 ] + gc[ 1 ] * searchRadius;
      const float oz = ( dim == 3 ) ? gridOrigin[ 2 ] + gc[ 2 ] * searchRadius : 0.f;
      const float sr = searchRadius;

      if constexpr( dim == 2 ) {
         // 4 corners of the pixel (z = 0)
         f << ox      << " " << oy      << " 0\n";
         f << ox + sr << " " << oy      << " 0\n";
         f << ox      << " " << oy + sr << " 0\n";
         f << ox + sr << " " << oy + sr << " 0\n";
      } else {
         // 8 corners of the voxel
         f << ox      << " " << oy      << " " << oz      << "\n";
         f << ox + sr << " " << oy      << " " << oz      << "\n";
         f << ox      << " " << oy + sr << " " << oz      << "\n";
         f << ox + sr << " " << oy + sr << " " << oz      << "\n";
         f << ox      << " " << oy      << " " << oz + sr << "\n";
         f << ox + sr << " " << oy      << " " << oz + sr << "\n";
         f << ox      << " " << oy + sr << " " << oz + sr << "\n";
         f << ox + sr << " " << oy + sr << " " << oz + sr << "\n";
      }
   }

   // --- Cells ---
   f << "CELLS " << n << " " << n * ( pointsPerCell + 1 ) << "\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ ) {
      f << pointsPerCell;
      for( int p = 0; p < pointsPerCell; p++ )
         f << " " << ci * pointsPerCell + p;
      f << "\n";
   }

   // --- Cell types ---
   // VTK_PIXEL = 8 (2D),  VTK_VOXEL = 11 (3D)
   const int vtkCellType = ( dim == 2 ) ? 8 : 11;
   f << "CELL_TYPES " << n << "\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << vtkCellType << "\n";

   // --- Cell data: flat cell index for debugging ---
   f << "CELL_DATA " << n << "\n";
   f << "SCALARS cell_index int 1\n";
   f << "LOOKUP_TABLE default\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << cellsHost[ ci ] << "\n";

   f.close();
   std::cout << "saveZoneToVTK: wrote " << n << " cells to " << filename << std::endl;
}
*/

template< typename ParticleConfig, typename DeviceType >
void
ParticleZone< ParticleConfig, DeviceType >::saveZoneToVTK(
   const std::string&    filename,
   const IndexVectorType gridSize,
   const PointType       gridOrigin,
   const RealType        searchRadius ) const
{
   constexpr int dim = ParticleConfig::spaceDimension;

   // Copy to host
   TNL::Containers::Array< GlobalIndexType, TNL::Devices::Host, GlobalIndexType >
         cells( this->numberOfCellsInZone );
   cells = this->cellsInZone;

   const GlobalIndexType n          = this->numberOfCellsInZone;
   const int             cornersPerCell = ( dim == 2 ) ? 4 : 8;
   const int             vtkType        = ( dim == 2 ) ? 8 : 11;  // VTK_PIXEL / VTK_VOXEL

   std::ofstream f( filename );
   if( !f )
      throw std::runtime_error( "saveZoneToVTK: cannot open " + filename );

   // -----------------------------------------------------------------------
   // Header
   // -----------------------------------------------------------------------
   f << "# vtk DataFile Version 3.0\n"
     << "ParticleZone\n"
     << "ASCII\n"
     << "DATASET UNSTRUCTURED_GRID\n";

   // -----------------------------------------------------------------------
   // Points — one cell = cornersPerCell points
   // -----------------------------------------------------------------------
   f << "POINTS " << n * cornersPerCell << " float\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ ) {

      //// Decode flat index → integer grid coords
      //GlobalIndexType idx = cells[ ci ];
      //IndexVectorType gc  = 0;
      //if constexpr( dim == 2 ) {
      //   gc[ 0 ] = idx % gridSize[ 0 ];
      //   gc[ 1 ] = idx / gridSize[ 0 ];
      //} else {
      //   gc[ 0 ] = idx % gridSize[ 0 ];
      //   gc[ 1 ] = ( idx / gridSize[ 0 ] ) % gridSize[ 1 ];
      //   gc[ 2 ] = idx / ( gridSize[ 0 ] * gridSize[ 1 ] );
      //}

      // Physical origin of this cell
      const IndexVectorType gc = CellIndexer::GetCellCoordinates( cells[ ci ], gridSize );
      const float x0 = gridOrigin[ 0 ] + gc[ 0 ] * searchRadius;
      const float y0 = gridOrigin[ 1 ] + gc[ 1 ] * searchRadius;
      const float z0 = ( dim == 3 ) ? gridOrigin[ 2 ] + gc[ 2 ] * searchRadius : 0.f;
      const float sr = static_cast< float >( searchRadius );

      // Write corners in VTK_PIXEL / VTK_VOXEL order
      if constexpr( dim == 2 ) {
         f << x0      << " " << y0      << " 0\n";
         f << x0 + sr << " " << y0      << " 0\n";
         f << x0      << " " << y0 + sr << " 0\n";
         f << x0 + sr << " " << y0 + sr << " 0\n";
      } else {
         f << x0      << " " << y0      << " " << z0      << "\n";
         f << x0 + sr << " " << y0      << " " << z0      << "\n";
         f << x0      << " " << y0 + sr << " " << z0      << "\n";
         f << x0 + sr << " " << y0 + sr << " " << z0      << "\n";
         f << x0      << " " << y0      << " " << z0 + sr << "\n";
         f << x0 + sr << " " << y0      << " " << z0 + sr << "\n";
         f << x0      << " " << y0 + sr << " " << z0 + sr << "\n";
         f << x0 + sr << " " << y0 + sr << " " << z0 + sr << "\n";
      }
   }

   // -----------------------------------------------------------------------
   // Cells
   // -----------------------------------------------------------------------
   f << "CELLS " << n << " " << n * ( cornersPerCell + 1 ) << "\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ ) {
      f << cornersPerCell;
      for( int p = 0; p < cornersPerCell; p++ )
         f << " " << ci * cornersPerCell + p;
      f << "\n";
   }

   // -----------------------------------------------------------------------
   // Cell types
   // -----------------------------------------------------------------------
   f << "CELL_TYPES " << n << "\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << vtkType << "\n";

   // -----------------------------------------------------------------------
   // Cell data — flat index and grid coords for inspection in ParaView
   // -----------------------------------------------------------------------
   f << "CELL_DATA " << n << "\n";

   f << "SCALARS flat_index int 1\nLOOKUP_TABLE default\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << cells[ ci ] << "\n";

   f << "SCALARS gc_x int 1\nLOOKUP_TABLE default\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << cells[ ci ] % gridSize[ 0 ] << "\n";

   f << "SCALARS gc_y int 1\nLOOKUP_TABLE default\n";
   for( GlobalIndexType ci = 0; ci < n; ci++ )
      f << ( cells[ ci ] / gridSize[ 0 ] ) % gridSize[ 1 ] << "\n";

   if constexpr( dim == 3 ) {
      f << "SCALARS gc_z int 1\nLOOKUP_TABLE default\n";
      for( GlobalIndexType ci = 0; ci < n; ci++ )
         f << cells[ ci ] / ( gridSize[ 0 ] * gridSize[ 1 ] ) << "\n";
   }

   std::cout << "saveZoneToVTK: " << n << " cells → " << filename << "\n";
}

} // Particles
} // TNL
