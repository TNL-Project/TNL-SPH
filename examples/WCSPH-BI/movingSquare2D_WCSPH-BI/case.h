#include "template/config.h"
#include "template/userCodedFunctions.h"

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
{
   userCodedFunctions::CustomMotion motion( "template/Motion_Body.dat" );
   EnergyMonitor energyMonitor( sph.fluid, true  );
   ForceMonitor forceMonitor( sph.boundary );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch( true );

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // FEATURE: monitor energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );
      // FEATURE: monitor forces
      forceMonitor.computeForces( sph.fluid, sph.boundary, sph.modelParams, 1 );

      // asssing square motion
      motion.assignMotion( sph.boundary, sph.timeStepping );

      //integrate
      sph.integrateVerletStep();

      // FEATURE: shift particles
      PST::shift( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );
      // FEATURE: filter density
      const int filteringFrequency = sph.userParams.template getParameter< int >( "filtering-steps-interval" );
      if( sph.timeStepping.getStep() % filteringFrequency == 0 )
         DensityFilter::filterDensity( sph.fluid, sph.modelParams );

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

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
{
   userCodedFunctions::CustomMotion motion( "template/Motion_Body.dat" );
   EnergyMonitor energyMonitor( sph.fluid, true  );
   ForceMonitor forceMonitor( sph.boundary );

   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch( true );

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // // compute new time step
      // sph.computeTimeStep();
      // sph.timeStepping.outputTimeStep( sph.outputDirectory + "/timeStep.dat" );

      //integrate predictor step
      sph.symplecticVerletPredictor();

      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // FEATURE: monitor energy levels
      energyMonitor.computeEnergyDerivatives( sph.fluid, sph.modelParams );
      energyMonitor.integrate( sph.timeStepping.getTimeStep() );

      //integrate
      sph.symplecticVerletCorrector();

      // asssing square motion
      motion.assignMotion( sph.boundary, sph.timeStepping );

      // FEATURE: monitor energy levels
      energyMonitor.computeEnergyLevels( sph.fluid, sph.modelParams );
      // FEATURE: monitor forces
      forceMonitor.computeForces( sph.fluid, sph.boundary, sph.modelParams, 1 );

      // output particle data
      sph.makeSnapshot();
      energyMonitor.output( sph.outputDirectory + "/energy.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );
      forceMonitor.output( sph.outputDirectory + "/force.dat", sph.timeStepping.getStep(), sph.timeStepping.getTime() );

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }
}

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   // FEATURE: parse custom user defined params
   sph.initUserConfig( userCodedFunctions::userConfigSetup );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

