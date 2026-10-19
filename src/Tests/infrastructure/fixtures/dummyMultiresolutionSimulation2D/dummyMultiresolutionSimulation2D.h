#include "dummyConfig2D.h"

int main( int argc, char* argv[] )
{
   std::ofstream logFile( "results/dummyMultiresolutionSimulation2D.log" );
   Simulation sph( logFile );
   sph.init( argc, argv );
   sph.writeProlog();
}

