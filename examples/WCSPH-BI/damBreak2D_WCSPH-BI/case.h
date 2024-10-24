#include <TNL/Config/parseCommandLine.h>
#include <TNL/Config/parseINIConfigFile.h>
#include <TNL/Logger.h>

#include <SPH/configSetup.h>
#include <SPH/configInit.h>
#include "template/config.h"

#include <SPH/Models/WCSPH_BI/control.h>
#include <SPH/DefaultTimeLoops.h>

template< typename FluidPointer,
          typename BoudaryPointer,
          typename TimeScheme = typename SPHDefs::IntegrationScheme,
          std::enable_if_t< std::is_same_v< TimeScheme, TNL::SPH::IntegrationSchemes::SymplecticVerletScheme< typename SPHDefs::SPHConfig > >, bool > Enabled = true >
void
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

template< typename FluidPointer,
          typename BoudaryPointer,
          typename TimeScheme = typename SPHDefs::IntegrationScheme,
          std::enable_if_t< std::is_same_v< TimeScheme, TNL::SPH::IntegrationSchemes::MidpointScheme< typename SPHDefs::SPHConfig > >, bool > Enabled = true >
void
exec( Simulation& sph, TNL::Logger log )
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

      //int midpointIteration = 0;
      //const int midpointMaxInterations = 50; //10 bdflt
      //const RealType residuaTolerance = 1.e-6 * 264.0753900000003f;
      //RealType midpointRelaxCoef = 0;
      //const RealType midpointRelaxCoef_0 = midpointRelaxCoef;
      //const RealType residaMinimualDecay = 1.f / 5.f;
      //RealType residuaPrevious = 1.f;
      //const RealType midpointRelaxCoefIncrement = 0.25f;

      const int midpointMaxInterations = 3000; //10 bdflt
      const RealType residuaTolerance = 1e-4f; // 1.e-6 * 264.0753900000003f
      RealType midpointRelaxCoef = 0.;
      const RealType midpointRelaxCoef_0 = 0;
      const RealType residaMinimualDecay = 1.f / 2.f;
      RealType residuaPrevious = 1.f;
      //const RealType midpointRelaxCoefIncrement = 0.25;
      //const RealType midpointRelaxCoefIncrement = 0.1;
      const RealType midpointRelaxCoefIncrement = 0.001;

      sph.timeMeasurement.start( "integrate" );
      sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate: predictor...", "Done." );

      int midpointIteration = 0;
      float residualPrevious = 0.f;
      float midpointRelaxCoef = sph.modelParams.midpointRelaxCoef;

      while( midpointIteration < midpointMaxInterations ) {

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
         sph.interact(); //TODO: What about BC conditions?
         sph.timeMeasurement.stop( "interact" );
         sph.writeLog( log, "Interact...", "Done." );

         // custom: no-penetration bc
         BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

         // compute residua
         sph.timeMeasurement.start( "integrate" );
         const RealType residual = sph.integrator->midpointResiduals( sph.fluid, sph.modelParams );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: compute residuals...", "Done." );

         // stop midpoint iterations
         if( residual < residuaTolerance )
            midpointIteration = midpointMaxInterations;
         // constrol residua decay
         if( midpointIteration > 0 )
            if( residual / residualPrevious > sph.modelParams.residaMinimualDecay )
                midpointRelaxCoef = sph.modelParams.midpointRelaxCoefIncrement + ( 1.0 - sph.modelParams.midpointRelaxCoefIncrement ) * midpointRelaxCoef;
         // backup residua
         residualPrevious = residual;

         // relax
         sph.timeMeasurement.start( "integrate" );
         sph.integrator->relax( sph.fluid, sph.modelParams, midpointRelaxCoef, midpointIteration );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: relax...", "Done." );

         std::cout  << "Midpoint iteractions: " << midpointIteration << " residua: " << maxResidua << " relax coef: " << midpointRelaxCoef << std::endl;
         midpointIteration++;
      }

      sph.timeMeasurement.start( "integrate" );
      if( sph.timeStepping.getStep() == 0 )
         sph.integrator->corrector( 0, sph.fluid );
      else
         sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
      sph.timeMeasurement.stop( "integrate" );
      sph.writeLog( log, "Integrate: corrector...", "Done." );

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

         //// search for neighbros
         //sph.timeMeasurement.start( "search" );
         //sph.performNeighborSearch( log );
         //sph.timeMeasurement.stop( "search" );
         //sph.writeLog( log, "Search...", "Done." );

         //// perform interaction with given model
         //sph.timeMeasurement.start( "interact" );
         //sph.interact(); //TODO: What about BC conditions?
         //sph.timeMeasurement.stop( "interact" );
         //sph.writeLog( log, "Interact...", "Done." );

      sph.timeStepping.updateTimeStep();
   }

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


   sph.fluid->variables->a = 0.f;
   sph.fluid->variables->drho = 0.f;

   if constexpr ( std::is_same_v< SPHDefs::IntegrationScheme, TNL::SPH::IntegrationSchemes::MidpointScheme< typename SPHDefs::SPHConfig > > ){
      using RealType = typename Simulation::RealType;

      while( sph.timeStepping.runTheSimulation() )
      {
         std::cout << "Step: " << sph.timeStepping.getStep() << std::endl;

         //int midpointIteration = 0;
         //const int midpointMaxInterations = 50; //10 bdflt
         //const RealType residuaTolerance = 1.e-6 * 264.0753900000003f;
         //RealType midpointRelaxCoef = 0;
         //const RealType midpointRelaxCoef_0 = midpointRelaxCoef;
         //const RealType residaMinimualDecay = 1.f / 5.f;
         //RealType residuaPrevious = 1.f;
         //const RealType midpointRelaxCoefIncrement = 0.25f;

         int midpointIteration = 0;
         const int midpointMaxInterations = 3000; //10 bdflt
         const RealType residuaTolerance = 1e-4f; // 1.e-6 * 264.0753900000003f
         RealType midpointRelaxCoef = 0.;
         const RealType midpointRelaxCoef_0 = 0;
         const RealType residaMinimualDecay = 1.f / 2.f;
         RealType residuaPrevious = 1.f;
         //const RealType midpointRelaxCoefIncrement = 0.25;
         //const RealType midpointRelaxCoefIncrement = 0.1;
         const RealType midpointRelaxCoefIncrement = 0.001;

         //// search for neighbros
         //sph.timeMeasurement.start( "search" );
         //sph.performNeighborSearch( log );
         //sph.timeMeasurement.stop( "search" );
         //sph.writeLog( log, "Search...", "Done." );

         //// perform interaction with given model
         //sph.timeMeasurement.start( "interact" );
         //sph.interact(); //TODO: What about BC conditions?
         //sph.timeMeasurement.stop( "interact" );
         //sph.writeLog( log, "Interact...", "Done." );

         //THIS IS TERRIBLY WRONG! THIS OVERWRITE THE LAST ITERAION FIELDS
         //sph.fluid->variables->drho = sph.fluid->integratorVariables->drhodt_in;
         //sph.fluid->variables->a = sph.fluid->integratorVariables->dvdt_in;

         sph.timeMeasurement.start( "integrate" );
         sph.integrator->predictor( sph.timeStepping.getTimeStep(), sph.fluid );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: predictor...", "Done." );

         //THIS IS IMHO NOT NECESSARY
         sph.fluid->variables->drho = sph.fluid->integratorVariables->drhodt_in;
         sph.fluid->variables->a = sph.fluid->integratorVariables->dvdt_in;

         while( midpointIteration < midpointMaxInterations )
         {
            //
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
            sph.interact(); //TODO: What about BC conditions?
            sph.timeMeasurement.stop( "interact" );
            sph.writeLog( log, "Interact...", "Done." );

            // custom: no-penetration bc
            BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

            // compute residua
            sph.timeMeasurement.start( "integrate" );
            sph.integrator->midpointResiduals( sph.fluid, sph.modelParams );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: compute residuals...", "Done." );

            // backup residua and get new maxima
            sph.timeMeasurement.start( "integrate" );
            const RealType maxResidua = sph.integrator->getMaxResidua( sph.fluid );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: reduce residuals...", "Done." );
            if( maxResidua < residuaTolerance )
               midpointIteration = midpointMaxInterations;

            // constrol residua decay
            if( midpointIteration > 0)
               if( maxResidua / residuaPrevious > residaMinimualDecay )
                   midpointRelaxCoef = midpointRelaxCoefIncrement + ( 1.0 - midpointRelaxCoefIncrement ) * midpointRelaxCoef; //TODO: Not sure here

            // backup residua
            residuaPrevious = maxResidua;

            // relax
            sph.timeMeasurement.start( "integrate" );
            if( midpointIteration == 0 )
               sph.integrator->relax( sph.fluid, midpointRelaxCoef_0 );
            else if( midpointIteration == midpointMaxInterations )
               sph.integrator->relax( sph.fluid, 0.f );
            else
               sph.integrator->relax( sph.fluid, midpointRelaxCoef );
            sph.timeMeasurement.stop( "integrate" );
            sph.writeLog( log, "Integrate: relax...", "Done." );

            std::cout  << "Midpoint iteractions: " << midpointIteration << " residua: " << maxResidua << " relax coef: " << midpointRelaxCoef << std::endl;
            midpointIteration++;
         }
         //midpointIteration = 0;

         sph.timeMeasurement.start( "integrate" );
         if( sph.timeStepping.getStep() == 0 )
            sph.integrator->corrector( 0, sph.fluid );
         else
            sph.integrator->corrector( sph.timeStepping.getTimeStep(), sph.fluid );
         sph.timeMeasurement.stop( "integrate" );
         sph.writeLog( log, "Integrate: corrector...", "Done." );

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

            //// search for neighbros
            //sph.timeMeasurement.start( "search" );
            //sph.performNeighborSearch( log );
            //sph.timeMeasurement.stop( "search" );
            //sph.writeLog( log, "Search...", "Done." );

            //// perform interaction with given model
            //sph.timeMeasurement.start( "interact" );
            //sph.interact(); //TODO: What about BC conditions?
            //sph.timeMeasurement.stop( "interact" );
            //sph.writeLog( log, "Interact...", "Done." );

         sph.timeStepping.updateTimeStep();
      }

   }

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

