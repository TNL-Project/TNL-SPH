#include "dummyConfig3D.h"

int main( int argc, char* argv[] )
{
   std::ofstream logFile( "results/dummyMultiresolutionSimulation3D.log" );
   Simulation sph( logFile );
   sph.init( argc, argv );
   sph.writeProlog();
}
