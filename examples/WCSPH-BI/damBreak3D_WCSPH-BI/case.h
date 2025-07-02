#include "template/config.h"

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
{
   while( sph.timeStepping.runTheSimulation() )
   {
      // search for neighbros
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      //integrate
      sph.integrateVerletStep();

      // output particle data
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.timeStepping.updateTimeStep();
   }
}

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::MidpointScheme< typename Simulation::SPHConfig >
>
void exec( Simulation& sph )
{
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
      }

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

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

