
for( unsigned int p = 0; p < mySPHSimulation.particles.getNumberOfParticles(); p ++ )
{
  if( ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] == 0. ) &&
      ( mySPHSimulation.particles.getPoint( p )[ 1 ] == 0. ) ) ||
      ( ( ( mySPHSimulation.particles.getPoint( p )[ 0 ] > 0.1 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 0 ] < 0.1 + eps ) ) &&
      ( ( mySPHSimulation.particles.getPoint( p )[ 1 ] > 0.09 - eps ) && ( mySPHSimulation.particles.getPoint( p )[ 1 ] < 0.09 + eps ) ) ) )
  {
    std::cout << std::endl << "  Particle found." << std::endl;
    std::cout << "  Coordinates: " << mySPHSimulation.particles.getPoint( p ) << " with type: " << mySPHSimulation.model.vars.type[ p ]  << std::endl;
    std::cout << "  Velocity: " << mySPHSimulation.model.vars.v[ p ] << ". Density: " << mySPHSimulation.model.vars.rho[ p ]  <<  ". Pressure: " << mySPHSimulation.model.vars.p[ p ] << std::endl;
    std::cout << "  DrhoDv: " << mySPHSimulation.model.vars.DrhoDv[ p ] << "." << std::endl << std::endl;
  }
}
