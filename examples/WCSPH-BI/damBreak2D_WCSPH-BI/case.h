#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

template< typename Simulation,
          typename IntegrationScheme = typename Simulation::ModelParams::IntegrationScheme,
          typename std::enable_if_t< std::is_same_v< IntegrationScheme, TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename Simulation::SPHConfig > >, bool > Enabled = true >
void
exec( Simulation& sph, TNL::Logger& log )
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

      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // output particle data
      if( sph.timeStepping.checkOutputTimer( "save_results" ) ){
         // compute pressure from density
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

template< typename Simulation,
          typename IntegrationScheme = typename Simulation::ModelParams::IntegrationScheme,
          typename std::enable_if_t< std::is_same_v< IntegrationScheme, TNL::SPH::IntegrationSchemes::MidpointScheme< typename Simulation::SPHConfig > >, bool > Enabled = true >
void
exec( Simulation& sph, TNL::Logger& log )
{
   // search for neighbros
   sph.timeMeasurement.start( "search" );
   sph.performNeighborSearch( log, true );
   sph.timeMeasurement.stop( "search" );
   sph.writeLog( log, "Search...", "Done." );

   // perform interaction with given model
   sph.timeMeasurement.start( "interact" );
   sph.interact(); //TODO: What about BC conditions?
   sph.timeMeasurement.stop( "interact" );
   sph.writeLog( log, "Interact...", "Done." );

   while( sph.timeStepping.runTheSimulation() ){

      sph.timeMeasurement.start( "integrate" );
      sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate: predictor...", "Done." );

      int midpointIteration = 0;
      float residualPrevious = 0.f;
      float midpointRelaxCoef = sph.modelParams.midpointRelaxCoef;

      while( midpointIteration < sph.modelParams.midpointMaxInterations ) {

         // backup derivatives
         sph.fluid->integratorVariables->drhodt_in = sph.fluid->variables->drho;
         sph.fluid->integratorVariables->dvdt_in = sph.fluid->variables->a;

         // update inner loop variables
         sph.timeMeasurement.start( "integrate" );
         sph.integrator->midpointUpdateVariables( sph.timeStepping.getTimeStep(), sph.fluid );
         //sph.integrator->midpointUpdatePositions( sph.timeStepping.getTimeStep(), sph.fluid );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: midpoint update...", "Done." );

         // search for neighbros
         sph.timeMeasurement.start( "search" );
         sph.performNeighborSearch( log, true );
         sph.timeMeasurement.stop( "search" );
         sph.writeLog( log, "Search...", "Done." );

         // perform interaction with given model
         sph.timeMeasurement.start( "interact" );
         sph.interact();
         sph.timeMeasurement.stop( "interact" );
         sph.writeLog( log, "Interact...", "Done." );

         // custom: no-penetration bc
         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

         // compute residua
         sph.timeMeasurement.start( "integrate" );
         const float residual = sph.integrator->midpointResiduals( sph.fluid, sph.modelParams );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: compute residuals...", "Done." );

         // stop midpoint iterations
         if( residual < sph.modelParams.midpointResidualTolerance )
            midpointIteration = sph.modelParams.midpointMaxInterations;
         // constrol residua decay
         if( midpointIteration > 0 )
            if( residual / residualPrevious > sph.modelParams.midpointResidualMinimalDecay )
                midpointRelaxCoef = sph.modelParams.midpointRelaxCoefIncrement + ( 1.0 - sph.modelParams.midpointRelaxCoefIncrement ) * midpointRelaxCoef;
         // backup residua
         residualPrevious = residual;

         // relax
         sph.timeMeasurement.start( "integrate" );
         sph.integrator->relax( sph.fluid, sph.modelParams, midpointRelaxCoef, midpointIteration );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: relax...", "Done." );

         midpointIteration++;
      }
      std::cout  << "Midpoint iteractions residua: " << residualPrevious << " relax coef: " << midpointRelaxCoef << std::endl;

      sph.timeMeasurement.start( "integrate" );
      if( sph.timeStepping.getStep() == 0 )
         sph.integrator->corrector( 0, sph.fluid );
      else
         sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate: corrector...", "Done." );

      // output particle data
      if( sph.timeStepping.checkOutputTimer( "save_results" ) ){
         // compute pressure from density
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

template< typename Simulation,
          typename IntegrationScheme = typename Simulation::ModelParams::IntegrationScheme,
          typename std::enable_if_t< std::is_same_v< IntegrationScheme, TNL::SPH::IntegrationSchemes::RK45Scheme< typename Simulation::SPHConfig > >, bool > Enabled = true >
void
exec( Simulation& sph, TNL::Logger& log )
{
   while( sph.timeStepping.runTheSimulation() ){
      for( int predictorStep = 0; predictorStep < 4; predictorStep++ ){
         // predictor step
         sph.timeMeasurement.start( "integrate" );
         sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid, predictorStep );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: predictor...", "Done." );

         // search for neighbros
         sph.timeMeasurement.start( "search" );
         sph.performNeighborSearch( log );
         sph.timeMeasurement.stop( "search" );
         sph.writeLog( log, "Search...", "Done." );

         // perform interaction with given model
         sph.timeMeasurement.start( "interact" );
         sph.interact(); //TODO: What about BC conditions?
         sph.timeMeasurement.stop( "interact" );
         sph.writeLog( log, "Interact...", "Done." );

         // custom: no-penetration bc
         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );
      }

      // predictor step
      sph.timeMeasurement.start( "integrate" );
      sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate: predictor...", "Done." );

       // output particle data
       if( sph.timeStepping.checkOutputTimer( "save_results" ) ) {
          // compute pressure from density
          sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
          sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

          sph.save( log );
       }
      sph.timeStepping.updateTimeStep();
   }
}

int main( int argc, char* argv[] )
{
   // prepare client parameters
   TNL::Config::ParameterContainer cliParams;
   TNL::Config::ConfigDescription cliConfig;

   // prepare sph parameters
   TNL::Config::ParameterContainer parameters;
   TNL::Config::ConfigDescription config;

   try {
      TNL::SPH::template initialize< Simulation >( argc, argv, cliParams, cliConfig, parameters, config );
   }
   catch ( ... ) {
      return EXIT_FAILURE;
   }

   TNL::Logger log( 100, std::cout );
   Simulation sph;
   sph.init( parameters, log );
   sph.writeProlog( log );
   exec( sph, log );
   sph.writeEpilog( log );
}

