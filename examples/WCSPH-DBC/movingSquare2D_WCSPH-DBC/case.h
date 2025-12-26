#include "template/config.h"
#include "template/userCodedFunctions.h"

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
{
   // FEATURE: motion and energy monitor
   userCodedFunctions::CustomMotion motion( "template/Motion_Body.dat" );
   EnergyMonitor energyMonitor( sph.fluid, true  );
   ForceMonitor forceMonitor( sph.boundary );

   while( sph.timeStepping.runTheSimulation() )
   {
      // FEATURE: asssing square motion
      motion.assignMotion( sph.boundary, sph.timeStepping );

      // search for neighbros
      sph.performNeighborSearch( true );

      // perform interaction with given model
      sph.interact();

      // FEATURE: monitor energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );
      // FEATURE: monitor forces
      forceMonitor.computeForces( sph.fluid, sph.boundary, sph.modelParams, 1 );

      //integrate
      sph.integrateVerletStep( SPHDefs::BCType::integrateInTime() );

      // FEATURE: shift particles
      PST::shift( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // output particle data
      sph.makeSnapshot();
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );
      forceMonitor.output( sph.outputDirectory + "/force.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.timeStepping.updateTimeStep();
   }
}

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

