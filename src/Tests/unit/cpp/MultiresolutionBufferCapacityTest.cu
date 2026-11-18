#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>

#include "MultiresolutionBufferCapacityTestSupport.h"

namespace {

using namespace MultiresolutionBufferCapacityTest;

template< typename Pointer, bool lts >
void
checkBufferToFluid( int& failures )
{
   const std::string label = lts ? "LTS buffer-to-fluid" : "uniform buffer-to-fluid";
   FluidPointer fluid = makeSet< FluidPointer >( 3, 4 );
   Pointer buffer = makeBuffer< Pointer >( 3, 3 );
   buffer->numberOfParticlesToRetype = 2;

   const Point sourcePoints[] = { buffer->getPoints().getElement( 1 ), buffer->getPoints().getElement( 2 ) };
   const Point sourceVelocities[] = { buffer->getVariables()->v.getArray().getElement( 1 ),
                                       buffer->getVariables()->v.getArray().getElement( 2 ) };
   const float sourceDensities[] = { buffer->getVariables()->rho.getArray().getElement( 1 ),
                                      buffer->getVariables()->rho.getArray().getElement( 2 ) };
   const int sourceIds[] = { buffer->getVariables()->referentialIdx.getElement( 1 ),
                             buffer->getVariables()->referentialIdx.getElement( 2 ) };
   const int sourceMarkers[] = { buffer->getVariables()->marker.getArray().getElement( 1 ),
                                 buffer->getVariables()->marker.getArray().getElement( 2 ) };
   SnapshotMap expected = snapshotOf( fluid );
   const int removedId = fluid->getVariables()->referentialIdx.getElement( 0 );
   expected.erase( removedId );
   fluid->getPoints().setElement( 0, Point{ FLT_MAX, FLT_MAX } );
   fluid->getParticles()->setNumberOfParticlesToRemove( 1 );

   ModelParams params;
   if constexpr( lts )
      buffer->convertBufferToFluid( fluid, params );
   else
      buffer->convertBufferToFluid( fluid );

   expect( failures, fluid->getNumberOfAllocatedParticles() == 5, label + ": fluid grew before append" );
   expect( failures, fluid->getNumberOfParticles() == 4, label + ": compacted count refreshed before append" );
   expect( failures, buffer->getNumberOfParticles() == 1, label + ": converted buffer count" );
   checkPreserved( failures, fluid, expected, label );
   checkAbsent( failures, fluid, removedId, label );
   for( int i = 0; i < 2; ++i ) {
      const int target = i + 2;
      expect( failures, fluid->getPoints().getElement( target ) == sourcePoints[ i ], label + ": appended position" );
      expect( failures, fluid->getVariables()->rho.getArray().getElement( target ) == sourceDensities[ i ],
              label + ": appended density" );
      expect( failures, fluid->getVariables()->v.getArray().getElement( target ) == sourceVelocities[ i ],
              label + ": appended velocity" );
      expect( failures, fluid->getVariables()->drho.getArray().getElement( target ) == 0.f,
              label + ": appended neutral density rate" );
      expect( failures, fluid->getVariables()->p.getArray().getElement( target ) == 0.f,
              label + ": appended neutral pressure" );
      expect( failures, fluid->getVariables()->a.getArray().getElement( target ) == ( Point{ 0.f, 0.f } ),
              label + ": appended neutral acceleration" );
      expect( failures, fluid->getVariables()->gamma.getArray().getElement( target ) == 0.f,
              label + ": appended neutral gamma" );
      expect( failures, fluid->getVariables()->referentialIdx.getElement( target ) == sourceIds[ i ],
              label + ": appended stable referential ID" );
      expect( failures, fluid->getVariables()->marker.getArray().getElement( target ) == sourceMarkers[ i ],
              label + ": appended marker" );
      expect( failures, fluid->getIntegratorVariables()->rho_old.getElement( target ) == sourceDensities[ i ],
              label + ": appended density history" );
      expect( failures, fluid->getIntegratorVariables()->v_old.getElement( target ) == sourceVelocities[ i ],
              label + ": appended velocity history" );
   }
}

template< typename Pointer, bool lts >
void
checkFluidToBuffer( int& failures )
{
   const std::string label = lts ? "LTS fluid-to-buffer" : "uniform fluid-to-buffer";
   FluidPointer fluid = makeSet< FluidPointer >( 2, 3 );
   Pointer buffer = makeBuffer< Pointer >( 3, 4 );
   SnapshotMap expected = snapshotOf( buffer );
   const int removedId = buffer->getVariables()->referentialIdx.getElement( 0 );
   expected.erase( removedId );
   buffer->getPoints().setElement( 0, Point{ FLT_MAX, FLT_MAX } );
   buffer->getParticles()->setNumberOfParticlesToRemove( 1 );
   buffer->numberOfParticlesToBuffer = 2;
   buffer->particlesToBuffer.setSize( 2 );
   buffer->particlesToBuffer.setElement( 0, 0 );
   buffer->particlesToBuffer.setElement( 1, 1 );
   buffer->massNodes.mass.setElement( 0, 88.f );
   fluid->getVariables()->referentialIdx.setElement( 0, 3000 );
   fluid->getVariables()->referentialIdx.setElement( 1, 3037 );

   const Point sourcePoints[] = { fluid->getPoints().getElement( 0 ), fluid->getPoints().getElement( 1 ) };
   const Point sourceVelocities[] = { fluid->getVariables()->v.getArray().getElement( 0 ),
                                       fluid->getVariables()->v.getArray().getElement( 1 ) };
   const float sourceDensities[] = { fluid->getVariables()->rho.getArray().getElement( 0 ),
                                      fluid->getVariables()->rho.getArray().getElement( 1 ) };
   const int sourceIds[] = { fluid->getVariables()->referentialIdx.getElement( 0 ),
                             fluid->getVariables()->referentialIdx.getElement( 1 ) };
   const int sourceMarkers[] = { fluid->getVariables()->marker.getArray().getElement( 0 ),
                                 fluid->getVariables()->marker.getArray().getElement( 1 ) };
   ModelParams params;
   if constexpr( lts )
      buffer->convertFluidToBuffer( fluid, params );
   else
      buffer->convertFluidToBuffer( fluid );

   expect( failures, buffer->getNumberOfAllocatedParticles() == 5, label + ": buffer grew before append" );
   expect( failures, buffer->getNumberOfParticles() == 4, label + ": compacted count refreshed before append" );
   checkPreserved( failures, buffer, expected, label );
   checkAbsent( failures, buffer, removedId, label );
   checkScratchGrowth( failures, buffer, 5, 2 );
   expect( failures, buffer->massNodes.mass.getElement( 0 ) == 88.f, label + ": mass nodes untouched" );
   for( int i = 0; i < 2; ++i ) {
      const int target = i + 2;
      expect( failures, buffer->getPoints().getElement( target ) == sourcePoints[ i ], label + ": captured position" );
      expect( failures, buffer->getVariables()->rho.getArray().getElement( target ) == sourceDensities[ i ],
              label + ": captured density" );
      expect( failures, buffer->getVariables()->v.getArray().getElement( target ) == sourceVelocities[ i ],
              label + ": captured velocity" );
      expect( failures, buffer->getVariables()->drho.getArray().getElement( target ) == 0.f,
              label + ": captured neutral density rate" );
      expect( failures, buffer->getVariables()->p.getArray().getElement( target ) == 0.f,
              label + ": captured neutral pressure" );
      expect( failures, buffer->getVariables()->a.getArray().getElement( target ) == ( Point{ 0.f, 0.f } ),
              label + ": captured neutral acceleration" );
      expect( failures, buffer->getVariables()->gamma.getArray().getElement( target ) == 0.f,
              label + ": captured neutral gamma" );
      expect( failures, buffer->getVariables()->referentialIdx.getElement( target ) == sourceIds[ i ],
              label + ": captured stable referential ID" );
      expect( failures, buffer->getVariables()->marker.getArray().getElement( target ) == sourceMarkers[ i ],
              label + ": captured marker" );
      expect( failures, buffer->getIntegratorVariables()->rho_old.getElement( target ) == sourceDensities[ i ],
              label + ": captured density history" );
      expect( failures, buffer->getIntegratorVariables()->v_old.getElement( target ) == sourceVelocities[ i ],
              label + ": captured velocity history" );
      expect( failures, fluid->getPoints().getElement( i )[ 0 ] == FLT_MAX, label + ": source marked removed" );
   }
   expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 2, label + ": source removal count" );
   if( buffer->getNumberOfAllocatedParticles() >= buffer->getNumberOfParticles() ) {
      buffer->moveBufferParticles( 0.f );
      expect( failures, buffer->numberOfParticlesToRetype == buffer->getNumberOfParticles(),
              label + ": lifecycle reclassifies stale markers" );
   }
}

template< typename Pointer, bool lts >
void
checkCreateBufferParticles( int& failures )
{
   const std::string label = lts ? "LTS create" : "uniform create";
   Pointer buffer = makeBuffer< Pointer >( 3, 4 );
   SnapshotMap expected = snapshotOf( buffer );
   const int removedId = buffer->getVariables()->referentialIdx.getElement( 0 );
   expected.erase( removedId );
   buffer->getPoints().setElement( 0, Point{ FLT_MAX, FLT_MAX } );
   buffer->getParticles()->setNumberOfParticlesToRemove( 1 );
   buffer->numberOfParticlesToCreate = 2;
   buffer->particlesToBuffer.setSize( 7 );
   buffer->particlesToBuffer.setElement( 0, 19 );
   buffer->massNodes.setSize( 2 );
   buffer->massNodes.points.setElement( 0, Point{ 0.5f, 3.5f } );
   buffer->massNodes.points.setElement( 1, Point{ 1.5f, 3.5f } );
   buffer->massNodes.normal.setElement( 0, Point{ 0.f, 0.f } );
   buffer->massNodes.normal.setElement( 1, Point{ 2.f, 0.f } );
   buffer->massNodes.mass.setElement( 0, 41.f );
   buffer->massNodes.mass.setElement( 1, 42.f );

   ModelParams params;
   buffer->createBufferParticles( params );
   expect( failures, buffer->getNumberOfAllocatedParticles() == 5, label + ": buffer grew before append" );
   expect( failures, buffer->getNumberOfParticles() == 4, label + ": compacted count refreshed before append" );
   checkPreserved( failures, buffer, expected, label, 4323 );
   checkAbsent( failures, buffer, removedId, label );
   checkScratchGrowth( failures, buffer, 5, 7 );
   expect( failures, buffer->particlesToBuffer.getElement( 0 ) == 19, label + ": zone scratch contents untouched" );
   expect( failures, buffer->massNodes.mass.getElement( 0 ) == 41.f && buffer->massNodes.mass.getElement( 1 ) == 42.f,
           label + ": mass-node accumulation untouched" );
   expect( failures, buffer->getPoints().getElement( 2 ) == Point( { 0.5f, 3.5f } ), label + ": first spawned position" );
   expect( failures, buffer->getPoints().getElement( 3 ) == Point( { 2.f, 3.5f } ), label + ": second spawned position" );
   for( int i = 2; i < 4; ++i ) {
      expect( failures, buffer->getVariables()->rho.getArray().getElement( i ) == params.rho0,
              label + ": spawned density" );
      expect( failures, buffer->getVariables()->v.getArray().getElement( i ) == Point{ 0.f, 0.f },
              label + ": spawned velocity" );
   }
}

template< typename Pointer, bool lts >
void
checkCaptureRange( int& failures, const int width, const int bufferCount )
{
   const std::string label = lts ? "LTS capture" : "uniform capture";
   FluidPointer fluid = makeSet< FluidPointer >( 4, 4 );
   Pointer buffer = makeBuffer< Pointer >( 0, bufferCount );
   buffer->getParticles()->setNumberOfParticles( bufferCount );
   buffer->inner_overlap = true;
   buffer->outer_overlap = false;
   buffer->frameOrientation = -1;
   buffer->frameDimensions = Coordinates{ width, width };
   buffer->zone.setNumberOfParticlesPerCell( 8 );
   buffer->zone.assignCells( Coordinates{ 0, 0 }, Coordinates{ width, width }, Coordinates{ 4, 4 } );
   buffer->zone.updateParticlesInZone( fluid->getParticles() );
   const int zoneCount = buffer->zone.getNumberOfParticles();
   std::vector< int > expected;
   for( int i = 0; i < zoneCount; ++i )
      expected.push_back( buffer->zone.getParticlesInZone().getElement( i ) );
   std::sort( expected.begin(), expected.end() );
   const int zoneCapacity = buffer->zone.getNumberOfCells() * buffer->zone.getNumberOfParticlesPerCell();
   buffer->particlesToBuffer.setSize( zoneCapacity );

   buffer->getFluidParticlesEneringTheBuffer( fluid );
   expect( failures, zoneCount == width * width, label + ": fixture zone count" );
   expect( failures, buffer->numberOfParticlesToBuffer == zoneCount, label + ": captured count" );
   expect( failures, buffer->particlesToBuffer.getSize() == zoneCapacity, label + ": zone storage unchanged" );
   for( int i = 0; i < zoneCount; ++i )
      expect( failures, buffer->particlesToBuffer.getElement( i ) == expected[ i ], label + ": sorted capture index" );
}

template< typename Pointer, bool lts >
void
checkAll( int& failures )
{
   checkBufferToFluid< Pointer, lts >( failures );
   checkFluidToBuffer< Pointer, lts >( failures );
   checkCreateBufferParticles< Pointer, lts >( failures );
   checkCaptureRange< Pointer, lts >( failures, 1, 24 );  // buffer count exceeds zone capacity
   checkCaptureRange< Pointer, lts >( failures, 2, 1 );   // zone count exceeds buffer count
}

}  // namespace

int
main()
{
   int failures = 0;
   checkAll< BufferPointer, false >( failures );
   checkAll< BufferLTSPointer, true >( failures );
   if( failures != 0 ) {
      std::cerr << failures << " multiresolution append checks failed\n";
      return 1;
   }
   std::cout << "PASS: both MR buffer variants grow preserving appends and bound capture sorting by the zone\n";
   return 0;
}
