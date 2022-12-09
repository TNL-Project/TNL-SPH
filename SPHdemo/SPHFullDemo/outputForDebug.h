
for( unsigned int p = 0; p < mySPHSimulation.particles->getNumberOfParticles(); p ++ )
{
  if( ( ( mySPHSimulation.particles->getPoint( p )[ 0 ] == 0. ) &&
      ( mySPHSimulation.particles->getPoint( p )[ 1 ] == 0. ) ) ||
      ( ( ( mySPHSimulation.particles->getPoint( p )[ 0 ] > 0.1 - eps ) && ( mySPHSimulation.particles->getPoint( p )[ 0 ] < 0.1 + eps ) ) &&
      ( ( mySPHSimulation.particles->getPoint( p )[ 1 ] > 0.09 - eps ) && ( mySPHSimulation.particles->getPoint( p )[ 1 ] < 0.09 + eps ) ) ) )
  {
    std::cout << std::endl << "  Particle found." << std::endl;
    std::cout << "  Coordinates: " << mySPHSimulation.particles->getPoint( p ) << " with type: " << mySPHSimulation.model->type[ p ]  << std::endl;
		std::cout << "  Velocity: " << mySPHSimulation.model->v[ p ] << ". Density: " << mySPHSimulation.model->rho[ p ]  <<  ". Pressure: " << mySPHSimulation.model->p[ p ] << std::endl;
    std::cout << "  Drho: " << mySPHSimulation.model->drho[ p ] << " acc:" << mySPHSimulation.model->a[ p ] << std::endl << std::endl;
  }
}

