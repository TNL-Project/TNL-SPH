#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>
#include <SPH/DefaultTimeLoops.h>

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

   TNL::SPH::template exec( sph, log );

   //if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename SPHDefs::SPHConfig > > ) {
   //   while( sph.timeStepping.runTheSimulation() )
   //   {
   //      //integrate predictor step
   //      sph.timeMeasurement.start( "integrate" );
   //      sph.integrator->integratePredictorStep( sph.fluid, sph.boundary, sph.timeStepping );
   //      sph.timeMeasurement.stop( "integrate" );
   //      sph.writeLog( log, "Integrate - predictor step...", "Done." );

   //      // search for neighbros
   //      sph.timeMeasurement.start( "search" );
   //      sph.performNeighborSearch( log );
   //      sph.timeMeasurement.stop( "search" );
   //      sph.writeLog( log, "Search...", "Done." );

   //      // perform interaction with given model
   //      sph.timeMeasurement.start( "interact" );
   //      sph.interact(); //TODO: After the predictor step, there is apparently no reason to update BC
   //      sph.timeMeasurement.stop( "interact" );
   //      sph.writeLog( log, "Interact...", "Done." );
   //      std::cout << "A: " << sph.fluid->variables->a << std::endl;

   //      //integrate
   //      sph.timeMeasurement.start( "integrate" );
   //      sph.integrator->integrateCorrectorStep( sph.fluid, sph.boundary, sph.timeStepping );
   //      sph.timeMeasurement.stop( "integrate" );
   //      sph.writeLog( log, "Integrate - corrector step...", "Done." );

   //      // search for neighbros
   //      sph.timeMeasurement.start( "search" );
   //      sph.performNeighborSearch( log );
   //      sph.timeMeasurement.stop( "search" );
   //      sph.writeLog( log, "Search...", "Done." );

   //      // perform interaction with given model
   //      sph.timeMeasurement.start( "interact" );
   //      sph.interact();
   //      sph.timeMeasurement.stop( "interact" );
   //      sph.writeLog( log, "Interact...", "Done." );

   //      // output particle data
   //      if( sph.timeStepping.checkOutputTimer( "save_results" ) )
   //      {
   //         /**
   //          * Compute pressure from density.
   //          * This is not necessary since we do this localy, if pressure is needed.
   //          * It's useful for output anyway.
   //          */
   //         sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
   //         sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

   //         sph.save( log );
   //      }

   //      // check timers and if measurement or interpolation should be performed, is performed
   //      sph.template measure< SPHDefs::KernelFunction, SPHDefs::EOS >( log );

   //      // update time step
   //      sph.timeStepping.updateTimeStep();
   //   }
   //}
   //else{}// SymplecticVerletScheme

   //midpoint//sph.fluid->variables->a = 0.f;
   //midpoint//sph.fluid->variables->drho = 0.f;

   //midpoint//if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::MidpointScheme< typename SPHDefs::SPHConfig > > ){
   //midpoint//   using RealType = typename Simulation::RealType;

   //midpoint//   while( sph.timeStepping.runTheSimulation() )
   //midpoint//   {
   //midpoint//      std::cout << "Step: " << sph.timeStepping.getStep() << std::endl;

   //midpoint//      //int midpointIteration = 0;
   //midpoint//      //const int midpointMaxInterations = 50; //10 bdflt
   //midpoint//      //const RealType residuaTolerance = 1.e-6 * 264.0753900000003f;
   //midpoint//      //RealType midpointRelaxCoef = 0;
   //midpoint//      //const RealType midpointRelaxCoef_0 = midpointRelaxCoef;
   //midpoint//      //const RealType residaMinimualDecay = 1.f / 5.f;
   //midpoint//      //RealType residuaPrevious = 1.f;
   //midpoint//      //const RealType midpointRelaxCoefIncrement = 0.25f;

   //midpoint//      int midpointIteration = 0;
   //midpoint//      const int midpointMaxInterations = 3000; //10 bdflt
   //midpoint//      const RealType residuaTolerance = 1e-4f; // 1.e-6 * 264.0753900000003f
   //midpoint//      RealType midpointRelaxCoef = 0.;
   //midpoint//      const RealType midpointRelaxCoef_0 = 0;
   //midpoint//      const RealType residaMinimualDecay = 1.f / 2.f;
   //midpoint//      RealType residuaPrevious = 1.f;
   //midpoint//      //const RealType midpointRelaxCoefIncrement = 0.25;
   //midpoint//      //const RealType midpointRelaxCoefIncrement = 0.1;
   //midpoint//      const RealType midpointRelaxCoefIncrement = 0.001;

   //midpoint//      //// search for neighbros
   //midpoint//      //sph.timeMeasurement.start( "search" );
   //midpoint//      //sph.performNeighborSearch( log );
   //midpoint//      //sph.timeMeasurement.stop( "search" );
   //midpoint//      //sph.writeLog( log, "Search...", "Done." );

   //midpoint//      //// perform interaction with given model
   //midpoint//      //sph.timeMeasurement.start( "interact" );
   //midpoint//      //sph.interact(); //TODO: What about BC conditions?
   //midpoint//      //sph.timeMeasurement.stop( "interact" );
   //midpoint//      //sph.writeLog( log, "Interact...", "Done." );

   //midpoint//      //THIS IS TERRIBLY WRONG! THIS OVERWRITE THE LAST ITERAION FIELDS
   //midpoint//      //sph.fluid->variables->drho = sph.fluid->integratorVariables->drhodt_in;
   //midpoint//      //sph.fluid->variables->a = sph.fluid->integratorVariables->dvdt_in;

   //midpoint//      sph.timeMeasurement.start( "integrate" );
   //midpoint//      sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid );
   //midpoint//      sph.timeMeasurement.stop( "integrate" );
   //midpoint//      sph.writeLog( log, "Integrate: predictor...", "Done." );

   //midpoint//      //THIS IS IMHO NOT NECESSARY
   //midpoint//      sph.fluid->variables->drho = sph.fluid->integratorVariables->drhodt_in;
   //midpoint//      sph.fluid->variables->a = sph.fluid->integratorVariables->dvdt_in;

   //midpoint//      while( midpointIteration < midpointMaxInterations )
   //midpoint//      {
   //midpoint//         //
   //midpoint//         sph.fluid->integratorVariables->drhodt_in = sph.fluid->variables->drho;
   //midpoint//         sph.fluid->integratorVariables->dvdt_in = sph.fluid->variables->a;

   //midpoint//         // update inner loop variables
   //midpoint//         sph.timeMeasurement.start( "integrate" );
   //midpoint//         sph.integrator->midpointUpdateVariables( sph.timeStepping.getTimeStep(), sph.fluid );
   //midpoint//         //sph.integrator->midpointUpdatePositions( sph.timeStepping.getTimeStep(), sph.fluid );
   //midpoint//         sph.timeMeasurement.stop( "integrate" );
   //midpoint//         sph.writeLog( log, "Integrate: midpoint update...", "Done." );

   //midpoint//         // search for neighbros
   //midpoint//         sph.timeMeasurement.start( "search" );
   //midpoint//         sph.performNeighborSearch( log, true );
   //midpoint//         sph.timeMeasurement.stop( "search" );
   //midpoint//         sph.writeLog( log, "Search...", "Done." );

   //midpoint//         // perform interaction with given model
   //midpoint//         sph.timeMeasurement.start( "interact" );
   //midpoint//         sph.interact(); //TODO: What about BC conditions?
   //midpoint//         sph.timeMeasurement.stop( "interact" );
   //midpoint//         sph.writeLog( log, "Interact...", "Done." );

   //midpoint//         // custom: no-penetration bc
   //midpoint//         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

   //midpoint//         // compute residua
   //midpoint//         sph.timeMeasurement.start( "integrate" );
   //midpoint//         sph.integrator->midpointResiduals( sph.fluid, sph.modelParams );
   //midpoint//         sph.timeMeasurement.stop( "integrate" );
   //midpoint//         sph.writeLog( log, "Integrate: compute residuals...", "Done." );

   //midpoint//         // backup residua and get new maxima
   //midpoint//         sph.timeMeasurement.start( "integrate" );
   //midpoint//         const RealType maxResidua = sph.integrator->getMaxResidua( sph.fluid );
   //midpoint//         sph.timeMeasurement.stop( "integrate" );
   //midpoint//         sph.writeLog( log, "Integrate: reduce residuals...", "Done." );
   //midpoint//         if( maxResidua < residuaTolerance )
   //midpoint//            midpointIteration = midpointMaxInterations;

   //midpoint//         // constrol residua decay
   //midpoint//         if( midpointIteration > 0)
   //midpoint//            if( maxResidua / residuaPrevious > residaMinimualDecay )
   //midpoint//                midpointRelaxCoef = midpointRelaxCoefIncrement + ( 1.0 - midpointRelaxCoefIncrement ) * midpointRelaxCoef; //TODO: Not sure here

   //midpoint//         // backup residua
   //midpoint//         residuaPrevious = maxResidua;

   //midpoint//         // relax
   //midpoint//         sph.timeMeasurement.start( "integrate" );
   //midpoint//         if( midpointIteration == 0 )
   //midpoint//            sph.integrator->relax( sph.fluid, midpointRelaxCoef_0 );
   //midpoint//         else if( midpointIteration == midpointMaxInterations )
   //midpoint//            sph.integrator->relax( sph.fluid, 0.f );
   //midpoint//         else
   //midpoint//            sph.integrator->relax( sph.fluid, midpointRelaxCoef );
   //midpoint//         sph.timeMeasurement.stop( "integrate" );
   //midpoint//         sph.writeLog( log, "Integrate: relax...", "Done." );

   //midpoint//         std::cout  << "Midpoint iteractions: " << midpointIteration << " residua: " << maxResidua << " relax coef: " << midpointRelaxCoef << std::endl;
   //midpoint//         midpointIteration++;
   //midpoint//      }
   //midpoint//      //midpointIteration = 0;

   //midpoint//      sph.timeMeasurement.start( "integrate" );
   //midpoint//      if( sph.timeStepping.getStep() == 0 )
   //midpoint//         sph.integrator->corrector( 0, sph.fluid );
   //midpoint//      else
   //midpoint//         sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
   //midpoint//      sph.timeMeasurement.stop( "integrate" );
   //midpoint//      sph.writeLog( log, "Integrate: corrector...", "Done." );

   //midpoint//     // output particle data
   //midpoint//     if( sph.timeStepping.checkOutputTimer( "save_results" ) )
   //midpoint//     {
   //midpoint//        /**
   //midpoint//         * Compute pressure from density.
   //midpoint//         * This is not necessary since we do this localy, if pressure is needed.
   //midpoint//         * It's useful for output anyway.
   //midpoint//         */
   //midpoint//        sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
   //midpoint//        sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

   //midpoint//        sph.save( log );
   //midpoint//     }

   //midpoint//         //// search for neighbros
   //midpoint//         //sph.timeMeasurement.start( "search" );
   //midpoint//         //sph.performNeighborSearch( log );
   //midpoint//         //sph.timeMeasurement.stop( "search" );
   //midpoint//         //sph.writeLog( log, "Search...", "Done." );

   //midpoint//         //// perform interaction with given model
   //midpoint//         //sph.timeMeasurement.start( "interact" );
   //midpoint//         //sph.interact(); //TODO: What about BC conditions?
   //midpoint//         //sph.timeMeasurement.stop( "interact" );
   //midpoint//         //sph.writeLog( log, "Interact...", "Done." );

   //midpoint//      sph.timeStepping.updateTimeStep();
   //midpoint//   }

   //midpoint//}

   //if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::RK45Scheme< typename SPHDefs::SPHConfig > > ){
   //   while( sph.timeStepping.runTheSimulation() ){
   //      std::cout << "Time: " << sph.timeStepping.getTime() << " step: " << sph.timeStepping.getStep() << std::endl;

   //      for( int predictorStep = 0; predictorStep < 4; predictorStep++ ){
   //         // predictor step
   //         sph.timeMeasurement.start( "integrate" );
   //         sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid, predictorStep );
   //         sph.timeMeasurement.stop( "integrate" );
   //         sph.writeLog( log, "Integrate: predictor...", "Done." );

   //         // search for neighbros
   //         sph.timeMeasurement.start( "search" );
   //         sph.performNeighborSearch( log );
   //         sph.timeMeasurement.stop( "search" );
   //         sph.writeLog( log, "Search...", "Done." );

   //         // perform interaction with given model
   //         sph.timeMeasurement.start( "interact" );
   //         sph.interact(); //TODO: What about BC conditions?
   //         sph.timeMeasurement.stop( "interact" );
   //         sph.writeLog( log, "Interact...", "Done." );
   //         // custom: no-penetration bc
   //         //BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );
   //      }

   //      // predictor step
   //      sph.timeMeasurement.start( "integrate" );
   //      sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
   //      sph.timeMeasurement.stop( "integrate" );
   //      sph.writeLog( log, "Integrate: predictor...", "Done." );

   //       // output particle data
   //       if( sph.timeStepping.checkOutputTimer( "save_results" ) )
   //       {
   //          /**
   //           * Compute pressure from density.
   //           * This is not necessary since we do this localy, if pressure is needed.
   //           * It's useful for output anyway.
   //           */
   //          sph.model.computePressureFromDensity( sph.fluid, sph.modelParams );
   //          sph.model.computePressureFromDensity( sph.boundary, sph.modelParams );

   //          sph.save( log );
   //       }
   //      sph.timeStepping.updateTimeStep();
   //   }
   //}

   sph.writeEpilog( log );
}

