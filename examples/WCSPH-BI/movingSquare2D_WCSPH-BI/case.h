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

   while( sph.timeStepping.runTheSimulation() )
   {
      std::cout << "Step: " << sph.timeStepping.getStep() << " time: " << sph.timeStepping.getTime() << "." << std::endl;
      // search for neighbros
      sph.performNeighborSearch( true );

      // perform interaction with given model
      sph.interact();
      // custom: no-penetration bc
      BoundaryCorrection::boundaryCorrection( sph.fluid, sph.boundary, sph.modelParams, sph.timeStepping.getTimeStep() );

      // asssing square motion
      motion.assignMotion( sph.boundary, sph.timeStepping );

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

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();
   exec( sph );
   sph.writeEpilog();
}

