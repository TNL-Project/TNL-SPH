#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "MultiresolutionGhostRefreshTestModel.h"

namespace {

template< typename T >
std::string
join( const std::vector< T >& values )
{
   std::ostringstream out;
   for( std::size_t i = 0; i < values.size(); i++ ) {
      if( i != 0 )
         out << ", ";
      out << values[ i ];
   }
   return out.str();
}

void
expectEvents( int& failures, const std::string& name, const Probe& probe, const std::vector< std::string >& expected )
{
   if( probe.events != expected ) {
      failures++;
      std::cerr << "FAIL: " << name << " runtime event order\n  actual:   " << join( probe.events )
                << "\n  expected: " << join( expected ) << '\n';
   }
}

void
expectGhost( int& failures, const std::string& name, const Probe& probe, int level, float expected )
{
   if( probe.ghostAtInteraction[ level ] != expected ) {
      failures++;
      std::cerr << "FAIL: " << name << " level " << level << " observed ghost value "
                << probe.ghostAtInteraction[ level ] << ", expected assigned source value " << expected << '\n';
   }
}

void
checkUniformInteraction( int& failures )
{
   using Solver = TNL::SPH::SolverMultiSetBlockMultiresolution< OrderingModel >;
   Solver solver;
   Probe probe;
   solver.model.probe = &probe;
   setupSolver( solver );

   solver.interact();

   expectEvents(
      failures,
      "uniform MR",
      probe,
      { "wall:0",
        "finalize-wall:0",
        "wall:1",
        "finalize-wall:1",
        "refresh:0<-1",
        "refresh:1<-0",
        "interact:0",
        "finalize:0",
        "interact:1",
        "finalize:1" } );
   expectGhost( failures, "uniform MR", probe, 0, sourceValue( 1 ) );
   expectGhost( failures, "uniform MR", probe, 1, sourceValue( 0 ) );
}

void
checkLTSInteractionAfterEarlyRefresh( int& failures )
{
   using Solver = TNL::SPH::SolverMultiSetBlockMultiresolutionLocalTimeStepping< OrderingModel >;
   Solver solver;
   Probe probe;
   solver.model.probe = &probe;
   setupSolver( solver );

   // Reproduce the old active driver's pre-window trigger. A correct per-level
   // interaction must restore its own ghosts after its wall update even if this
   // stale all-level refresh happened earlier.
   solver.refreshBoundaryGhostValues();
   solver.interact( 0 );
   solver.interact( 1 );

   expectEvents(
      failures,
      "LTS MR",
      probe,
      { "refresh:0<-1",
        "refresh:1<-0",
        "wall:0",
        "finalize-wall:0",
        "refresh:0<-1",
        "interact:0",
        "finalize:0",
        "wall:1",
        "finalize-wall:1",
        "refresh:1<-0",
        "interact:1",
        "finalize:1" } );
   expectGhost( failures, "LTS MR", probe, 0, sourceValue( 1 ) );
   expectGhost( failures, "LTS MR", probe, 1, sourceValue( 0 ) );
}

}  // namespace

int
main()
{
   int failures = 0;
   checkUniformInteraction( failures );
   checkLTSInteractionAfterEarlyRefresh( failures );

   if( failures != 0 ) {
      std::cerr << failures << " boundary-ghost ordering regression checks failed\n";
      return 1;
   }
   std::cout << "PASS: wall updates, own-level ghost refreshes, and interactions have correct runtime order\n";
   return 0;
}
