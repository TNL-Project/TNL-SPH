#include "template/config.h"

int main( int argc, char* argv[] )
{
   TNL::MPI::ScopedInitializer mpi( argc, argv );

   std::ofstream logFile( "results/simulationLog_rank" + std::to_string( TNL::MPI::GetRank() ) );
   Simulation sph( logFile );
   sph.init( argc, argv );
   sph.writeProlog();
   TNL::MPI::Barrier( sph.communicator );

   while( sph.timeStepping.runTheSimulation() )
   {
      // keep me informed that everything is running
      if( ( sph.timeStepping.getStep() % 50  == 0 )  )
         sph.writeInfo();

      // search for neighbros
      sph.performNeighborSearch( true );
      TNL::MPI::Barrier( sph.communicator );

      // synchronize subdomains
      sph.synchronizeDistributedSimulation();
      TNL::MPI::Barrier( sph.communicator );

      // check timers and if load balancing should be preformed, it is performed
      sph.balanceSubdomains();

      // search for neighbros
      sph.performNeighborSearch( true );
      TNL::MPI::Barrier( sph.communicator );

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();
      TNL::MPI::Barrier( sph.communicator );

      // perform interaction with given model
      sph.interact();

      // reset subdomains overlaps
      sph.resetOverlaps();

      // make integration step with Verlet scheme
      sph.integrateVerletStep( SPHDefs::BCType::integrateInTime() );

      // check timers and if measurement or interpolation should be performed, it is performed
      sph.measure();

      // update time step - FIXME: Synchronize dt over all subdomains!
      sph.timeStepping.updateTimeStep();
      TNL::MPI::Barrier( sph.communicator );
   }

   sph.writeEpilog();
}

