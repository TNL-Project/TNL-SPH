#include <algorithm>
#include <cfloat>
#include <vector>

#include "ParticleSetCapacityTestSupport.h"

namespace {

using namespace ParticleSetCapacityTest;

std::vector< int >
currentIds( const FluidPointer& fluid )
{
   std::vector< int > ids;
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i )
      ids.push_back( fluid->getVariables()->referentialIdx.getElement( i ) );
   return ids;
}

void
checkNeighborSearch( int& failures, const FluidPointer& fluid )
{
   auto particles = fluid->getParticles();
   auto searchToken = particles->getSearchToken();
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i ) {
      Point query = fluid->getPoints().getElement( i );
      int visits = 0;
      auto countVisit = [] __cuda_callable__( int, int, Point&, int* count ) mutable
      {
         ++( *count );
      };
      Particles::NeighborsLoopAnotherSet::exec( i, query, searchToken, countVisit, &visits );
      expect( failures, visits > 0, "rebuilt cell list finds each active particle" );
   }
}

void
checkParticleSetCapacity( int& failures )
{
   FluidPointer fluid = makeFluid();
   const SnapshotMap initial = snapshot( fluid );
   SnapshotMap expected = initial;
   const std::vector< int > initialOrder = currentIds( fluid );
   auto* const initialPointStorage = fluid->getPoints().getData();
   auto* const initialReferentialStorage = fluid->getVariables()->referentialIdx.getData();

   expect( failures, ! fluid->ensureCapacity( 4 ), "smaller request is a no-op" );
   expect( failures, ! fluid->ensureCapacity( 6 ), "equal-capacity request is a no-op" );
   expect( failures, fluid->getNumberOfAllocatedParticles() == 6, "no-op requests never shrink" );
   expect( failures, fluid->getNumberOfParticles() == 5, "no-op requests preserve active count" );
   expect( failures, fluid->getPoints().getData() == initialPointStorage, "no-op requests retain point storage" );
   expect(
      failures,
      fluid->getVariables()->referentialIdx.getData() == initialReferentialStorage,
      "no-op requests retain referential storage" );
   expect( failures, currentIds( fluid ) == initialOrder, "no-op requests do not reorder particles" );
   checkSnapshot( failures, fluid, expected );

   expect( failures, fluid->ensureCapacity( 8 ), "larger request grows" );
   expect( failures, fluid->getNumberOfAllocatedParticles() == 8, "growth is exact" );
   expect( failures, fluid->getNumberOfParticles() == 5, "appended capacity remains inactive" );
   std::apply(
      [ & ]( auto&... fields )
      {
         ( expect( failures, fields.getArray().getSize() == 8, "every fluid field grows exactly" ), ... );
      },
      fluid->getVariables()->allFields() );
   expect( failures, fluid->getVariables()->referentialIdx.getSize() == 8, "referential IDs grow exactly" );
   expect( failures, fluid->getVariables()->referentialIdx_swap.getSize() == 8, "referential swap grows exactly" );
   expect( failures, fluid->getIntegratorVariables()->rho_old.getSize() == 8, "rho_old grows exactly" );
   expect( failures, fluid->getIntegratorVariables()->rho_old_swap.getSize() == 8, "rho_old swap grows exactly" );
   expect( failures, fluid->getIntegratorVariables()->v_old.getSize() == 8, "v_old grows exactly" );
   expect( failures, fluid->getIntegratorVariables()->v_old_swap.getSize() == 8, "v_old swap grows exactly" );
   checkSnapshot( failures, fluid, expected );
   checkNeighborSearch( failures, fluid );

   const int removedId = fluid->getVariables()->referentialIdx.getElement( 0 );
   fluid->getPoints().setElement( 0, Point{ FLT_MAX, FLT_MAX } );
   fluid->getParticles()->setNumberOfParticlesToRemove( 1 );
   expected.erase( removedId );
   expect( failures, fluid->ensureCapacity( 10 ), "second larger request grows" );
   expect( failures, fluid->getNumberOfAllocatedParticles() == 10, "second growth is exact" );
   expect( failures, fluid->getNumberOfParticles() == 4, "search compacts a marked removal" );
   checkSnapshot( failures, fluid, expected );
   const auto remainingIds = currentIds( fluid );
   expect(
      failures,
      std::find( remainingIds.begin(), remainingIds.end(), removedId ) == remainingIds.end(),
      "removed stable ID is absent after compaction" );
   expect( failures, ! fluid->ensureCapacity( 8 ), "post-growth smaller request remains a no-op" );
   expect( failures, fluid->getNumberOfAllocatedParticles() == 10, "post-growth capacity does not shrink" );
}

void
checkBoundaryAndOpenBoundaryStorage( int& failures )
{
   BoundaryVariables boundary;
   seedLeaf( boundary );
   boundary.n.getArray().setElement( 1, Point{ 4.f, 5.f } );
   boundary.elementSize.getArray().setElement( 1, 0.25f );
   checkLeafGrowth( failures, boundary );
   expect( failures, boundary.n.getArray().getElement( 1 ) == Point{ 4.f, 5.f }, "boundary normal preserved" );
   expect( failures, boundary.elementSize.getArray().getElement( 1 ) == 0.25f, "boundary element size preserved" );

   OpenBoundaryVariables openBoundary;
   seedLeaf( openBoundary );
   openBoundary.n.getArray().setElement( 1, Point{ 6.f, 7.f } );
   openBoundary.elementSize.getArray().setElement( 1, 0.5f );
   openBoundary.particleMark.getArray().setElement( 1, 31 );
   openBoundary.receivingParticleMark.getArray().setElement( 1, 47 );
   checkLeafGrowth( failures, openBoundary );
   expect( failures, openBoundary.n.getArray().getElement( 1 ) == Point{ 6.f, 7.f }, "open-boundary normal preserved" );
   expect( failures, openBoundary.elementSize.getArray().getElement( 1 ) == 0.5f, "open-boundary element preserved" );
   expect( failures, openBoundary.particleMark.getArray().getElement( 1 ) == 31, "particle mark preserved" );
   expect( failures, openBoundary.receivingParticleMark.getArray().getElement( 1 ) == 47, "receiving mark preserved" );
}

}  // namespace

int
main()
{
   int failures = 0;
   checkParticleSetCapacity( failures );
   checkBoundaryAndOpenBoundaryStorage( failures );
   if( failures != 0 ) {
      std::cerr << failures << " capacity-preservation checks failed\n";
      return 1;
   }
   std::cout << "PASS: grow-only capacity preserves fields, history, positions, referential IDs, and search data\n";
   return 0;
}
