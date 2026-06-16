#include "template/config.h"
#include <SPH/shared/removeParticlesOutOfDensityLimits.h>

template< typename Simulation >
requires std::is_same_v<
    typename Simulation::ModelParams::IntegrationScheme,
    TNL::SPH::IntegrationSchemes::VerletScheme< typename SPHDefs::SPHConfig >
>
void exec( Simulation& sph )
{
   while( sph.timeStepping.runTheSimulation() )
   {
      std::cout << "Step: " << sph.timeStepping.getStep() << " Time: " << sph.timeStepping.getTime()
         << " f0: "<< sph.fluidSets[ 0 ]->getNumberOfParticles()
         << " mb0: "<< sph.multiresolutionBoundaryPatches[ 0 ]->getNumberOfParticles()
         << " f1: "<< sph.fluidSets[ 1 ]->getNumberOfParticles()
         << " mb1: "<< sph.multiresolutionBoundaryPatches[ 1 ]->getNumberOfParticles()
         << std::endl;


      // search for neighbros
      //sph.removeParticlesOutOfDensityLimits();
      TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ 0 ], sph.modelParams );
      TNL::SPH::customFunctions::removeParticlesOutOfDensityLimits( sph.fluidSets[ 1 ], sph.modelParams );
      sph.removeParticlesOutOfDomain();
      sph.performNeighborSearch();

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 0 ], sph.boundarySets[ 0 ], sph.modelParams, sph.timeStepping.getTimeStep() );
      BoundaryCorrection::boundaryCorrection( sph.fluidSets[ 1 ], sph.boundarySets[ 1 ], sph.modelParams, sph.timeStepping.getTimeStep() );

      // integrate
      sph.integrate();

      // output particle data
      sph.makeSnapshot();

      ////DEBUG
      if( sph.timeStepping.getStep() == 12100) sph.save( false );

      // apply the multi-resoltion bc
      sph.multiresolutionUpdate();

      //// check timers and if measurement or interpolation should be performed, is performed
      //sph.measure();

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

