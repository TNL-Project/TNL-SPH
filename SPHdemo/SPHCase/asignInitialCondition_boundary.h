//using PointArrayType = typename ParticleSystem::PointArrayType;
//using ScalarArrayType = typename SPHModel::ScalarArrayType;
//using VectorArrayType = typename SPHModel::VectorArrayType;
//using ParticleTypeArrayType = typename SPHModel::ParticleTypeArrayType;

PointArrayType pointsLoaded_bound( ParticlesConfig_bound::numberOfParticles );
pointsLoaded_bound = particlesToRead_bound.getPoints();

ScalarArrayType pressureLoaded_bound( ParticlesConfig_bound::numberOfParticles );
pressureLoaded_bound = std::get< std::vector< float > >( myReader_bound.readPointData( "Pressure" ) );

ScalarArrayType densityLoaded_bound( ParticlesConfig_bound::numberOfParticles );
densityLoaded_bound = std::get< std::vector< float > >( myReader_bound.readPointData( "Density" ) );

VectorArrayType velocityLoaded_bound( ParticlesConfig_bound::numberOfParticles );
velocityLoaded_bound = std::get< std::vector< float > >( myReader_bound.readPointData( "Velocity" ) );


auto pointsLoaded_bound_view = pointsLoaded_bound.getView();
auto pLoaded_bound_view = pressureLoaded_bound.getView();
auto rhoLoaded_bound_view = densityLoaded_bound.getView();
auto vLoaded_bound_view = velocityLoaded_bound.getView();

auto points_bound_view = mySPHSimulation.particles_bound->getPoints().getView();
auto rho_bound_view = mySPHSimulation.model->getBoundaryVariables().rho.getView();
auto p_bound_view = mySPHSimulation.model->getBoundaryVariables().p.getView();
auto v_bound_view = mySPHSimulation.model->getBoundaryVariables().v.getView();

std::cout << std::endl << std::endl;
std::cout << "ParticlesConfig_bound::numberOfParticles: " << ParticlesConfig_bound::numberOfParticles << std::endl;
std::cout << "rhoLoaded_bound_view.getSize(): " << rhoLoaded_bound_view.getSize() << std::endl;
std::cout << std::endl << std::endl;

auto particleLoop_bound = [=] __cuda_callable__ ( int i  ) mutable
{
   points_bound_view[ i ] = pointsLoaded_bound_view.getElement( i );
   p_bound_view[ i ] = pLoaded_bound_view.getElement( i );
   rho_bound_view[ i ] = rhoLoaded_bound_view.getElement( i );
   v_bound_view[ i ] = vLoaded_bound_view.getElement( i );
};
Algorithms::ParallelFor< Device >::exec( 0, ParticlesConfig_bound::numberOfParticles, particleLoop_bound );

pressureLoaded_bound.reset();
densityLoaded_bound.reset();
velocityLoaded_bound.reset();
pointsLoaded_bound.reset();

