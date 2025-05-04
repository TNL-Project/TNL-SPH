#include "template/config.h"

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();

   while( sph.timeStepping.runTheSimulation() )
   {
      sph.performNeighborSearch(  );

      // perform interaction
      sph.interact();

      //integrate - update variables
      //sph.integrate( SPHDefs::IntegrationScheme::Stages::updateVariables );
      sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateVariables );

      // relax
      sph.model.relaxDistortion( sph.fluid, sph.timeStepping, sph.modelParams );

      //integrate
      //sph.integrate( SPHDefs::IntegrationScheme::Stages::moveParticles );
      sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::moveParticles );

      // output particle data
      sph.makeSnapshot();

      // check timers and if measurement or interpolation should be performed, is performed
      sph.measure();

      // update time step
      sph.updateTime();
   }

   sph.writeEpilog();
}

