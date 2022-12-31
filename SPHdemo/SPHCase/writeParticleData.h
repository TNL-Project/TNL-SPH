std::ofstream outputFile;
outputFile.open( outputFileName );

auto of_r_view = mySPHSimulation.particles->getPoints().getView();
auto of_rho_view = mySPHSimulation.model->getFluidVariables().rho.getView();
auto of_p_view = mySPHSimulation.model->getFluidVariables().p.getView();
auto of_v_view = mySPHSimulation.model->getFluidVariables().v.getView();

for( unsigned int p = 0; p < ParticlesConfig::numberOfParticles; p++ )
  outputFile << \
  of_r_view.getElement( p )[ 0 ] << " " << \
  0 << " " << \
  of_r_view.getElement( p )[ 1 ]<< " " << \
  of_v_view.getElement( p )[ 0 ] << " " << \
  0 << " " << \
  of_v_view.getElement( p )[ 1 ] << " " << \
  of_rho_view.getElement( p ) << " " << \
  of_p_view.getElement( p ) << std::endl;

