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

