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





//auto device_data_view = mySPHSimulation.particles->getPoints().getView();
auto type_view = mySPHSimulation.model->getParticleType().getView();
auto rho_view = mySPHSimulation.model->getRho().getView();
auto p_view = mySPHSimulation.model->getPress().getView();
auto v_view = mySPHSimulation.model->getVel().getView();

auto initCond = [=] __cuda_callable__ ( int p  ) mutable
{

	/*
  if( ( ( device_data_view[ p ][ 0 ] >  0.   - eps ) && ( device_data_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.01 - eps ) && ( device_data_view[ p ][ 0 ] < -0.01 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.02 - eps ) && ( device_data_view[ p ][ 0 ] < -0.02 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] >  0.   - eps ) && ( device_data_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.01 - eps ) && ( device_data_view[ p ][ 1 ] < -0.01 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.02 - eps ) && ( device_data_view[ p ][ 1 ] < -0.02 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.6  - eps ) && ( device_data_view[ p ][ 0 ] <  1.6  + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.61 - eps ) && ( device_data_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.62 - eps ) && ( device_data_view[ p ][ 0 ] <  1.62 + eps ) ) )
	*/
  if( ( ( device_data_view[ p ][ 0 ] >  0.   - eps ) && ( device_data_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.002 - eps ) && ( device_data_view[ p ][ 0 ] < -0.002 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.004 - eps ) && ( device_data_view[ p ][ 0 ] < -0.004 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] >  0.   - eps ) && ( device_data_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.002 - eps ) && ( device_data_view[ p ][ 1 ] < -0.002 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.004 - eps ) && ( device_data_view[ p ][ 1 ] < -0.004 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.608  - eps ) && ( device_data_view[ p ][ 0 ] <  1.608  + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.61 - eps ) && ( device_data_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.612 - eps ) && ( device_data_view[ p ][ 0 ] <  1.612 + eps ) ) )
	/*
  if( ( ( device_data_view[ p ][ 0 ] >  0.   - eps ) && ( device_data_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.001 - eps ) && ( device_data_view[ p ][ 0 ] < -0.001 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] > -0.002 - eps ) && ( device_data_view[ p ][ 0 ] < -0.002 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] >  0.   - eps ) && ( device_data_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.001 - eps ) && ( device_data_view[ p ][ 1 ] < -0.001 + eps ) ) ||
      ( ( device_data_view[ p ][ 1 ] > -0.002 - eps ) && ( device_data_view[ p ][ 1 ] < -0.002 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.609  - eps ) && ( device_data_view[ p ][ 0 ] <  1.609  + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.61 - eps ) && ( device_data_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( device_data_view[ p ][ 0 ] >  1.611 - eps ) && ( device_data_view[ p ][ 0 ] <  1.611 + eps ) ) )
	*/
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

