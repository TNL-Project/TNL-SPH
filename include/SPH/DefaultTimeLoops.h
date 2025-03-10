namespace TNL {
namespace SPH {

template< typename Simulation >
static void
exec( Simulation& sph, TNL::Logger log )
{
   while( sph.timeStepping.runTheSimulation() )
   {
      //integrate predictor step
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integratePredictorStep( sph.fluid, sph.boundary, sph.timeStepping );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate - predictor step...", "Done." );

      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      // perform interaction with given model
      //TODO: After the predictor step, there is apparently no reason to update BC
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      //integrate
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->integrateCorrectorStep( sph.fluid, sph.boundary, sph.timeStepping );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate - corrector step...", "Done." );

      // search for neighbros
      sph.timeMeasurement.start( "search" );
      sph.performNeighborSearch( log );
      sph.timeMeasurement.stop( "search" );
      sph.writeLog( log, "Search...", "Done." );

      // perform interaction with given model
      sph.timeMeasurement.start( "interact" );
      sph.interact();
      sph.timeMeasurement.stop( "interact" );
      sph.writeLog( log, "Interact...", "Done." );

      // output particle data
      if( sph.timeStepping.checkOutputTimer( "save_results" ) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * It's useful for output anyway.
          */
         sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
         sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

         sph.save( log );
      }

      // check timers and if measurement or interpolation should be performed, is performed
      sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

      // update time step
      sph.timeStepping.updateTimeStep();
   }
}

} // SPH
} // TNL

