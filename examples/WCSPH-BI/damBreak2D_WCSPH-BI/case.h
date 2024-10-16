#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>

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

   // Solver model:
   //sph.init( parameters );
   //sph.writeProlog( parameters );
   //sph.exec();
   //sph.writeEpilog( parameters );

   // Library model:
   //:if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename SPHDefs::SPHConfig > > ){
   //:   while( sph.timeStepping.runTheSimulation() )
   //:   {
   //:      //integrate predictor step
   //:      sph.timeMeasurement.start( "integrate" );
   //:      sph.integrator->integratePredictorStep( sph.fluid, sph.boundary, sph.timeStepping );
   //:      sph.timeMeasurement.stop( "integrate" );
   //:      sph.writeLog( log, "Integrate - predictor step...", "Done." );

   //:      // search for neighbros
   //:      sph.timeMeasurement.start( "search" );
   //:      sph.performNeighborSearch( log );
   //:      sph.timeMeasurement.stop( "search" );
   //:      sph.writeLog( log, "Search...", "Done." );

   //:      // perform interaction with given model
   //:      sph.timeMeasurement.start( "interact" );
   //:      sph.interact(); //TODO: After the predictor step, there is apparently no reason to update BC
   //:      sph.timeMeasurement.stop( "interact" );
   //:      sph.writeLog( log, "Interact...", "Done." );

   //:      //integrate
   //:      sph.timeMeasurement.start( "integrate" );
   //:      sph.integrator->integrateCorrectorStep( sph.fluid, sph.boundary, sph.timeStepping );
   //:      sph.timeMeasurement.stop( "integrate" );
   //:      sph.writeLog( log, "Integrate - corrector step...", "Done." );

   //:      // search for neighbros
   //:      sph.timeMeasurement.start( "search" );
   //:      sph.performNeighborSearch( log );
   //:      sph.timeMeasurement.stop( "search" );
   //:      sph.writeLog( log, "Search...", "Done." );

   //:      // perform interaction with given model
   //:      sph.timeMeasurement.start( "interact" );
   //:      sph.interact();
   //:      sph.timeMeasurement.stop( "interact" );
   //:      sph.writeLog( log, "Interact...", "Done." );

   //:      // output particle data
   //:      if( sph.timeStepping.checkOutputTimer( "save_results" ) )
   //:      {
   //:         /**
   //:          * Compute pressure from density.
   //:          * This is not necessary since we do this localy, if pressure is needed.
   //:          * It's useful for output anyway.
   //:          */
   //:         sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
   //:         sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

   //:         sph.save( log );
   //:      }

   //:      // check timers and if measurement or interpolation should be performed, is performed
   //:      sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

   //:      // update time step
   //:      sph.timeStepping.updateTimeStep();
   //:   }
   //:} // SymplecticVerletScheme

   if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::MidpointScheme< typename SPHDefs::SPHConfig > > ){
      using RealType = typename Simulation::RealType;

      while( sph.timeStepping.runTheSimulation() )
      {
         int midpointIteration = 0;
         const int midpointMaxInterations = 10;
         const RealType residuaTolerance = 0.001f;
         RealType midpointRelaxCoef = 0;
         const RealType midpointRelaxCoef_0 = midpointRelaxCoef;
         const RealType residaMinimualDecay = 1.f / 5.f;
         RealType residuaPrevious = 1.f;
         const RealType midpointRelaxCoefIncrement = 0.25f;

         sph.timeMeasurement.start( "integrate" );
         sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: predictor...", "Done." );

         while( midpointIteration < midpointMaxInterations )
         {
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

            // update inner loop variables
            sph.timeMeasurement.start( "integrate" );
            sph.integrator->midpointUpdatePositions( sph.timeStepping.getTimeStep(), sph.fluid );
            sph.integrator->midpointUpdateVariables( sph.timeStepping.getTimeStep(), sph.fluid );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: predictor...", "Done." );

            // compute residua
            sph.timeMeasurement.start( "integrate" );
            sph.integrator->midpointResiduals( sph.fluid, sph.modelParams );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: predictor...", "Done." );

            const RealType maxResidua = sph.integrator->getMaxResidua( sph.fluid );
            if( maxResidua < residuaTolerance )
               midpointIteration = midpointMaxInterations;

            // constrol residua decay
            if( maxResidua / residuaPrevious > residaMinimualDecay )
                midpointRelaxCoef = midpointRelaxCoefIncrement + ( 1.0 - midpointRelaxCoefIncrement ) * midpointRelaxCoef; //TODO: Not sure here

            // relax
            sph.timeMeasurement.start( "integrate" );
            if( midpointIteration == 0 )
               sph.integrator->relax( sph.fluid, sph.modelParams, midpointRelaxCoef_0 );
            else if( midpointIteration == midpointMaxInterations )
               sph.integrator->relax( sph.fluid, sph.modelParams, 0.f );
            else
               sph.integrator->relax( sph.fluid, sph.modelParams, midpointRelaxCoef );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: predictor...", "Done." );

            midpointIteration++;
         }
         midpointIteration = 0;

         sph.timeMeasurement.start( "integrate" );
         sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: predictor...", "Done." );
      }
   }

   sph.writeEpilog( log );
}

