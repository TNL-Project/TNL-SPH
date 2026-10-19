#include "dummyConfigLTS2D.h"
#include <TNL/Logger.h>
#include <cmath>
#include <fstream>
#include <string>

/* Dummy driver for the local-timestepping multiresolution solver: builds the dummy
 * multiresolution setup, writes the prolog, and runs host-side geometry assertions over
 * the LTS buffer patches (frame boxes at overlap 1/2, band-ring classification via the
 * frames, mass-node placement). Assertion results are piped into the same log for the
 * log-compare harness; the exit code is nonzero on any failure. */

using VectorType = typename Simulation::MultiresolutionBoundaryLTS::VectorType;
using CoordinatesType = typename Simulation::MultiresolutionBoundaryLTS::CoordinatesType;
using RealType = typename Simulation::MultiresolutionBoundaryLTS::RealType;

static bool
coordsEq( const CoordinatesType& a, const CoordinatesType& b )
{
   for( int d = 0; d < CoordinatesType::getSize(); d++ )
      if( a[ d ] != b[ d ] )
         return false;
   return true;
}

int main( int argc, char* argv[] )
{
   std::ofstream logFile( "results/dummyMultiresolutionLTSSimulation2D.log" );
   Simulation sph( logFile );
   sph.init( argc, argv );
   sph.writeProlog();

   TNL::Logger log( 100, logFile );
   using Patch = typename Simulation::MultiresolutionBoundaryLTS;

   int failures = 0;
   auto expect = [ & ]( const std::string& name, bool ok )
   {
      if( ! ok )
         failures++;
      log.writeParameter( name, ok ? "PASSED" : "FAILED" );
   };

   for( std::size_t p = 0; p < sph.multiresolutionBoundaryPatchesLTS.size(); p++ ) {
      auto& patch = *sph.multiresolutionBoundaryPatchesLTS[ p ];
      const auto& iface = sph.multiresolutionBoundaryPatchInterfaces[ p ];
      const std::string tag =
         "patch " + std::to_string( iface.ownIdx ) + "->" + std::to_string( iface.neighborIdx );

      const CoordinatesType frontO = patch.getFrameFrontOriginGlobalCoordinates();
      const CoordinatesType frontD = patch.getFrameFrontDimensions();
      const CoordinatesType midO = patch.getFrameMidOriginGlobalCoordinates();
      const CoordinatesType midD = patch.getFrameMidDimensions();
      const CoordinatesType backO = patch.getFrameBackOriginGlobalCoordinates();
      const CoordinatesType backD = patch.getFrameBackDimensions();

      const bool fineSide = sph.getLevelRefinementFactor( iface.ownIdx ) < sph.getLevelRefinementFactor( iface.neighborIdx );

      if( fineSide ) {
         expect( tag + ": frameMid = frameFront + 1 cell", coordsEq( midO, frontO - 1 ) && coordsEq( midD, frontD + 2 ) );
         expect( tag + ": frameBack = frameFront + 2 cells", coordsEq( backO, frontO - 2 ) && coordsEq( backD, frontD + 4 ) );
      }
      else {
         expect( tag + ": frameMid == frameFront (single layer)", coordsEq( midO, frontO ) && coordsEq( midD, frontD ) );
         expect( tag + ": frameBack = frameFront - 1 cell", coordsEq( backO, frontO + 1 ) && coordsEq( backD, frontD - 2 ) );
      }

       const RealType sr = patch.getParticles()->getSearchRadius();
       const RealType inv_sr = 1.f / sr;
       const VectorType refOrig = patch.getParticles()->getReferentialOrigin();

       /* Band-ring classification expressed directly through the front/mid/back frame
          boxes (the v2 patch keeps the two rings as band geometry: layer 1 = front XOR
          mid, layer 2 = mid XOR back). */
       const auto inFrame = [ & ]( const VectorType& point, const CoordinatesType& origin, const CoordinatesType& dims )
       {
          const CoordinatesType iGlob = TNL::floor( ( point - refOrig ) * inv_sr );
          return Patch::isInsideBox( iGlob - origin, dims );
       };
       const auto isInLayer1 = [ & ]( const VectorType& point )
       {
          return inFrame( point, frontO, frontD ) != inFrame( point, midO, midD );
       };
       const auto isInLayer2 = [ & ]( const VectorType& point )
       {
          return inFrame( point, midO, midD ) != inFrame( point, backO, backD );
       };

      const auto cellPoint = [ & ]( int cx, int cy ) -> VectorType
      {
         return refOrig + VectorType( 0.5f * sr ) + ( frontO + CoordinatesType{ cx, cy } ) * sr;
      };

      if( fineSide ) {
         expect( tag + ": first outer ring is layer 1",
                isInLayer1( cellPoint( frontD[ 0 ], 5 ) )
                   && isInLayer1( cellPoint( -1, 8 ) )
                   && ! isInLayer2( cellPoint( frontD[ 0 ], 5 ) ) );
         expect( tag + ": second outer ring is layer 2",
                isInLayer2( cellPoint( frontD[ 0 ] + 1, 5 ) )
                   && isInLayer2( cellPoint( -2, 8 ) )
                   && ! isInLayer1( cellPoint( frontD[ 0 ] + 1, 5 ) ) );
         expect( tag + ": third outer ring outside the band",
                ! isInLayer1( cellPoint( frontD[ 0 ] + 2, 5 ) )
                   && ! isInLayer2( cellPoint( frontD[ 0 ] + 2, 5 ) ) );
         expect( tag + ": interior cells in no layer",
                ! isInLayer1( cellPoint( 5, 5 ) ) && ! isInLayer2( cellPoint( 5, 5 ) ) );
         expect( tag + ": corner rings classified",
                isInLayer1( cellPoint( frontD[ 0 ], frontD[ 1 ] ) )
                   && isInLayer2( cellPoint( frontD[ 0 ] + 1, frontD[ 1 ] + 1 ) ) );
      }
      else {
         expect( tag + ": no layer-1 particles on single-layer side",
                ! isInLayer1( cellPoint( frontD[ 0 ], 3 ) )
                   && ! isInLayer1( cellPoint( -1, 3 ) )
                   && ! isInLayer1( cellPoint( 3, 3 ) ) );
         expect( tag + ": inner band ring is layer 2",
                isInLayer2( cellPoint( frontD[ 0 ] - 1, 3 ) )
                   && isInLayer2( cellPoint( 0, 3 ) ) );
         expect( tag + ": cells outside frameFront in no layer",
                ! isInLayer1( cellPoint( frontD[ 0 ], 3 ) )
                   && ! isInLayer2( cellPoint( frontD[ 0 ], 3 ) )
                   && ! isInLayer2( cellPoint( -1, 3 ) ) );
      }

      const auto& massNodes = patch.getMassNodes();
      expect( tag + ": mass nodes initialized", massNodes.numberOfMassNodes > 0 );

      TNL::Containers::Array< VectorType, TNL::Devices::Host > hostPoints;
      hostPoints = massNodes.points;

      bool allInsideBack = true;
      bool allOnBackFace = true;
      for( int i = 0; i < massNodes.numberOfMassNodes; i++ ) {
         const VectorType r = hostPoints[ i ];
         for( int d = 0; d < CoordinatesType::getSize(); d++ ) {
            const RealType lo = refOrig[ d ] + backO[ d ] * sr;
            const RealType hi = refOrig[ d ] + ( backO[ d ] + backD[ d ] ) * sr;
            if( r[ d ] < lo - 1e-3f * sr || r[ d ] > hi + 1e-3f * sr )
               allInsideBack = false;
         }
         bool onFace = false;
         for( int d = 0; d < CoordinatesType::getSize(); d++ ) {
            const RealType lo = refOrig[ d ] + backO[ d ] * sr;
            const RealType hi = lo + backD[ d ] * sr;
            if( std::fabs( r[ d ] - lo ) < 1e-3f * sr || std::fabs( r[ d ] - hi ) < 1e-3f * sr )
               onFace = true;
         }
         if( ! onFace )
            allOnBackFace = false;
      }
      expect( tag + ": mass nodes inside frameBack", allInsideBack );
      expect( tag + ": mass nodes on frameBack faces", allOnBackFace );
   }

   expect( "LTS geometry checks overall", failures == 0 );
   log.writeSeparator();

   return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
