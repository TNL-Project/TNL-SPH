#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

template< typename Simulation >
requires std::is_same_v<
   typename Simulation::ModelParams::IntegrationScheme,
   TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void
exec( Simulation& sph )
{
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // in case of variable time step, compute the step
      sph.computeTimeStep();

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      // make integration step with Verlet scheme
      sph.integrateVerletStep();

      // check timers and if output should be performed, it is performed
      sph.makeSnapshot();
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // check timers and if measurement or interpolation should be performed, it is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }
}

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename SPHDefs::SPHConfig >
>
void
exec( Simulation& sph )
{
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // compute new time step
      sph.computeTimeStep();
      sph.timeStepping.outputTimeStep( sph.outputDirectory + "/timeStep.dat" );

      //integrate predictor step
      sph.symplecticVerletPredictor();

      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      //integrate
      sph.symplecticVerletCorrector();

      // compute and outpute energy levels
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );

      // output particle data
      sph.makeSnapshot();
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }
}

template< typename Simulation >
requires std::is_same_v<
   typename Simulation::ModelParams::IntegrationScheme,
   TNL::SPH::IntegrationSchemes::MidpointScheme< typename Simulation::SPHConfig >
>
void
exec( Simulation& sph )
{
   EnergyFields energyMonitor;
   energyMonitor.init( sph.fluid, true );

   // search for neighbros
   sph.performNeighborSearch( true );

   while( sph.timeStepping.runTheSimulation() )
   {
      sph.midpointPredictor();

      while( sph.integrator->runMidpointSubiteration( sph.modelParams ) )
      {
         // update inner loop variables
         sph.midpointUpdateVariables();

         // search for neighbros
         sph.performNeighborSearch( true );

         // perform interaction with given model
         sph.interact();
         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

         // compute residuals and update relaxation factors
         sph.midpointResidualsAndRelaxationFactor();

         // relax
         sph.midpointRelax();

         // save information about midpoint iteration
         TNL::SPH::Info::midpointSchemeOutputLog(
               sph.integrator, sph.modelParams, sph.timeStepping, sph.outputDirectory + "/midpointInfo.dat" );
      }

      // compute and outpute energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // integration corrector
      sph.midpointCorrector();

      // post integration domain adjustemnt
      BoundaryCorrection::boundaryCorrectionPST( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // output particle data
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }
}

template< typename Simulation >
requires std::is_same_v<
   typename Simulation::ModelParams::IntegrationScheme,
   TNL::SPH::IntegrationSchemes::RK4Scheme< typename Simulation::SPHConfig >
>
void
exec( Simulation& sph )
{
   while( sph.timeStepping.runTheSimulation() ) {
      for( int predictorStep = 0; predictorStep < 4; predictorStep++ ) {

         // predictor step
         sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid, predictorStep );

         // search for neighbros
         sph.performNeighborSearch( true );

         // perform interaction with given model
         sph.interact();
         // custom: no-penetration bc
         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );
      }

      // predictor step
      sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );

      // output particle data
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }
}

int
main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

