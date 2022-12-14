auto type_view = mySPHSimulation.model->getParticleType().getView();
auto rho_view = mySPHSimulation.model->getRho().getView();
auto p_view = mySPHSimulation.model->getPress().getView();
auto v_view = mySPHSimulation.model->getVel().getView();
auto r_view = mySPHSimulation.particles->getPoints().getView();

//integrator
//auto v_view = mySPHSimulation.model->getVel().getView();

for( unsigned int p = 0; p < mySPHSimulation.particles->getNumberOfParticles(); p ++ )
{

	/*
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
		*/

	/*
  if( ( ( r_view[ p ][ 0 ] >  0.   - eps ) && ( r_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.01 - eps ) && ( r_view[ p ][ 0 ] < -0.01 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.02 - eps ) && ( r_view[ p ][ 0 ] < -0.02 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] >  0.   - eps ) && ( r_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.01 - eps ) && ( r_view[ p ][ 1 ] < -0.01 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.02 - eps ) && ( r_view[ p ][ 1 ] < -0.02 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.6  - eps ) && ( r_view[ p ][ 0 ] <  1.6  + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.61 - eps ) && ( r_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.62 - eps ) && ( r_view[ p ][ 0 ] <  1.62 + eps ) ) )
	*/
  if( ( ( r_view[ p ][ 0 ] >  0.   - eps ) && (   r_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.002 - eps ) && (  r_view[ p ][ 0 ] < -0.002 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.004 - eps ) && (  r_view[ p ][ 0 ] < -0.004 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] >  0.   - eps ) && (   r_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.002 - eps ) && (  r_view[ p ][ 1 ] < -0.002 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.004 - eps ) && (  r_view[ p ][ 1 ] < -0.004 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.608  - eps ) && ( r_view[ p ][ 0 ] <  1.608  + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.61 - eps ) && (   r_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.612 - eps ) && (  r_view[ p ][ 0 ] <  1.612 + eps ) ) )
	/*
  if( ( ( r_view[ p ][ 0 ] >  0.   - eps ) && (   r_view[ p ][ 0 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.001 - eps ) && (  r_view[ p ][ 0 ] < -0.001 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] > -0.002 - eps ) && (  r_view[ p ][ 0 ] < -0.002 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] >  0.   - eps ) && (   r_view[ p ][ 1 ] <  0.   + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.001 - eps ) && (  r_view[ p ][ 1 ] < -0.001 + eps ) ) ||
      ( ( r_view[ p ][ 1 ] > -0.002 - eps ) && (  r_view[ p ][ 1 ] < -0.002 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.609  - eps ) && ( r_view[ p ][ 0 ] <  1.609  + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.61 - eps ) && (   r_view[ p ][ 0 ] <  1.61 + eps ) ) ||
      ( ( r_view[ p ][ 0 ] >  1.611 - eps ) && (  r_view[ p ][ 0 ] <  1.611 + eps ) ) )
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

