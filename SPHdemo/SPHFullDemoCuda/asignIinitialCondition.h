/* DEVICE */
//Copy particles

#include <new>
using PointArrayType = typename ParticleSystem::PointArrayType;
using PointArrayPointer = typename Pointers::SharedPointer< PointArrayType, Device >;

PointArrayPointer myArray( ParticlesConfig::numberOfParticles );


auto host_data_view = particlesToRead.getPoints().getView();
auto device_data_view = mySPHSimulation.particles->getPoints().getView();

Pointers::synchronizeSmartPointersOnDevice< Devices::Cuda >();

for( unsigned int p = 0; p < ParticlesConfig::numberOfParticles; p ++ )
{
	myArray->setElement( p, host_data_view[ p ] );
}

auto particleLoop = [=] __cuda_callable__ ( int i  ) mutable
{
	 //device_data_view[ i ] = host_data_view.getElement( i );
	 device_data_view[ i ] = myArray->getElement( i );
	 //host_data_view[ i ] = {1., 1.};
};
Algorithms::ParallelFor< Device >::exec( 0, ParticlesConfig::numberOfParticles, particleLoop );

for( unsigned int p = 0; p < ParticlesConfig::numberOfParticles; p ++ )
{
	//mySPHSimulation.particles->setPoint( p, particlesToRead.getPoint( p ));
	//std::cout << mySPHSimulation.particles->getPoint( p ) << std::endl;
}

//std::cout << mySPHSimulation.particles->getPoints() << std::endl;




/* HOST */
/*
for( unsigned int p = 0; p < mySPHSimulation.particles.getNumberOfParticles(); p ++ )
{

  if( ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] >  0.   - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] <  0.   + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] > -0.01 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] < -0.01 + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] > -0.02 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] < -0.02 + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 1 ] >  0.   - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 1 ] <  0.   + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 1 ] > -0.01 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 1 ] < -0.01 + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 1 ] > -0.02 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 1 ] < -0.02 + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] >  1.6  - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] <  1.6  + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] >  1.61 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] <  1.61 + eps ) ) ||
      ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] >  1.62 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] <  1.62 + eps ) ) )
  {
    mySPHSimulation.model.vars.type[ p ] = 1.;
  }
  else
  {
    mySPHSimulation.model.vars.type[ p ] = 0.;
  }

    mySPHSimulation.model.vars.rho[ p ] = 1000.;
    mySPHSimulation.model.vars.rho[ p ] = 1000.;
    mySPHSimulation.model.vars.v[ p ] = 0.;

    //fill in integrator arrays
    mySPHSimulation.model.integrator.rhoO[ p ] = 1000.;
    mySPHSimulation.model.integrator.rhoOO[ p ] = 1000.;

    mySPHSimulation.model.integrator.vO[ p ] = 0.;
    mySPHSimulation.model.integrator.vOO[ p ] = 0.;

}
*/

//auto device_data_view = mySPHSimulation.particles->getPoints().getView();
auto type_view = mySPHSimulation.model->getParticleType().getView();
auto rho_view = mySPHSimulation.model->getRho().getView();
auto p_view = mySPHSimulation.model->getPress().getView();
auto v_view = mySPHSimulation.model->getVel().getView();

auto initCond = [=] __cuda_callable__ ( int p  ) mutable
{

  if( ( ( device_data_view[ p ][ 0 ] >  0.   - eps ) && ( device_data_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.01 - eps ) && ( device_data_view[ p ][ 0 ] < -0.01 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.02 - eps ) && ( device_data_view[ p ][ 0 ] < -0.02 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] >  0.   - eps ) && ( device_data_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.01 - eps ) && ( device_data_view[ p ][ 1 ] < -0.01 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.02 - eps ) && ( device_data_view[ p ][ 1 ] < -0.02 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.6  - eps ) && ( device_data_view[ p ][ 0 ] <  1.6  + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.61 - eps ) && ( device_data_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.62 - eps ) && ( device_data_view[ p ][ 0 ] <  1.62 + eps ) ) )
  {
    type_view[ p ] = 1.;
  }
  else
  {
    type_view[ p ] = 0.;
  }

    rho_view[ p ] = 1000.;
    p_view[ p ] = 0.;
    v_view[ p ] = 0.;

    ////fill in integrator arrays
    //mySPHSimulation.model.integrator.rhoO[ p ] = 1000.;
    //mySPHSimulation.model.integrator.rhoOO[ p ] = 1000.;

    //mySPHSimulation.model.integrator.vO[ p ] = 0.;
    //mySPHSimulation.model.integrator.vOO[ p ] = 0.;

};
Algorithms::ParallelFor< Device >::exec( 0, ParticlesConfig::numberOfParticles, initCond );

//std::cout << mySPHSimulation.model->getParticleType() << std::endl;

//  for( unsigned int p = 0; p < mySPHSimulation.particles.getNumberOfParticles(); p ++ )
//  {
//
//    if( ( mySPHSimulation.particles.getPoint( p )[ 0 ] == 0. ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 0 ] == -0.01 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 0 ] == -0.02 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 0. ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == -0.01 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == -0.02 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 1.6 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 1.61 ) ||
//        ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 1.62 ) )
//    {
//      mySPHSimulation.model.vars.type[ p ] = 1.;
//    }
//    else
//    {
//      mySPHSimulation.model.vars.type[ p ] = 0.;
//    }
//
//      mySPHSimulation.model.vars.rho[ p ] = 1000.;
//      mySPHSimulation.model.vars.rho[ p ] = 1000.;
//      mySPHSimulation.model.vars.v[ p ] = 0.;
//
//      //fill in integrator arrays
//      mySPHSimulation.model.integrator.rhoO[ p ] = 1000.;
//      mySPHSimulation.model.integrator.rhoOO[ p ] = 1000.;
//
//      mySPHSimulation.model.integrator.vO[ p ] = 0.;
//      mySPHSimulation.model.integrator.vOO[ p ] = 0.;
//
//  }

