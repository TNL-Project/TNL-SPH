using PointArrayType = typename ParticleSystem::PointArrayType;
using ScalarArrayType = typename SPHModel::ScalarArrayType;
using VectorArrayType = typename SPHModel::VectorArrayType;
using ParticleTypeArrayType = typename SPHModel::ParticleTypeArrayType;

PointArrayType pointsLoaded( ParticlesConfig::numberOfParticles );
pointsLoaded = particlesToRead.getPoints();
ScalarArrayType pressureLoaded( ParticlesConfig::numberOfParticles );
pressureLoaded = std::get< std::vector< float > >( myReader.readPointData( "Pressure" ) );
ScalarArrayType densityLoaded( ParticlesConfig::numberOfParticles );
densityLoaded = std::get< std::vector< float > >( myReader.readPointData( "Density" ) );
VectorArrayType velocityLoaded( ParticlesConfig::numberOfParticles );
velocityLoaded = std::get< std::vector< float > >( myReader.readPointData( "Velocity" ) );

auto pointsLoaded_view = pointsLoaded.getView();
auto pLoaded_view = pressureLoaded.getView();
auto rhoLoaded_view = densityLoaded.getView();
auto vLoaded_view = velocityLoaded.getView();

auto points_view = mySPHSimulation.particles->getPoints().getView();
auto rho_view = mySPHSimulation.model->getFluidVariables().rho.getView();
auto p_view = mySPHSimulation.model->getFluidVariables().p.getView();
auto v_view = mySPHSimulation.model->getFluidVariables().v.getView();

auto particleLoop = [=] __cuda_callable__ ( int i  ) mutable
{
   points_view[ i ] = pointsLoaded_view.getElement( i );
   p_view[ i ] = pLoaded_view.getElement( i );
   rho_view[ i ] = rhoLoaded_view.getElement( i );
   v_view[ i ] = vLoaded_view.getElement( i );
};
Algorithms::ParallelFor< Device >::exec( 0, ParticlesConfig::numberOfParticles, particleLoop );

pressureLoaded.reset();
densityLoaded.reset();
velocityLoaded.reset();
pointsLoaded.reset();

