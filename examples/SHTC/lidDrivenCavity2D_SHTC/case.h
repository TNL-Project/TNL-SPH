#include "template/config.h"

template< typename ObjectPointer >
void setInitialDistortion( ObjectPointer object )
{
   auto A_view = object->getVariables()->A.getView();
   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      A_view[ i ] = { 1, 0, 0 , 1 };
      //A_view[ i ] = { 1, 0, 0,
      //                0, 1, 0,
      //                0, 0, 1 };
   };
   object->getParticles()->forAll( init );
}

int main( int argc, char* argv[] )
{
   Simulation sph;
   sph.init( argc, argv );
   sph.writeProlog();

   //FIXME: Custom initial condition
   setInitialDistortion( sph.fluid );
   setInitialDistortion( sph.boundary );

   while( sph.timeStepping.runTheSimulation() )
   //while( sph.timeStepping.getStep() < 4 )
   {
      std::cout << "\n-------- dt: " << sph.timeStepping.getTime() << " k: " << sph.timeStepping.getStep() << "--------------------------------------------\n" << std::endl;

      sph.performNeighborSearch();

      //not symplecic: // perform interaction
      //not symplecic: sph.interact();
      //not symplecic: //integrate - update variables
      //not symplecic: //sph.integrate( SPHDefs::IntegrationScheme::Stages::updateVariables );
      //not symplecic: sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateVariables );
      //not symplecic: sph.integrator->integrate( sph.boundary, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateVariables ); //FIXME: UGLY HACK

      sph.interact();
      sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateVelocity ); //FIXME: UGLY HACK
      sph.interact();
      sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateDensity ); //FIXME: UGLY HACK
      sph.interact();
      sph.integrator->integrate( sph.fluid, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateDistortion ); //FIXME: UGLY HACK
      sph.integrator->integrate( sph.boundary, sph.boundary, sph.timeStepping, SPHDefs::IntegrationScheme::Stages::updateDistortion ); //FIXME: UGLY HACK

      // relax
      sph.model.relaxDistortion( sph.fluid, sph.timeStepping, sph.modelParams );
      sph.model.relaxDistortion( sph.boundary, sph.timeStepping, sph.modelParams ); //FIMXE

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

