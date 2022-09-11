
for(unsigned int i = 0; i < myParticleSystem.getNumberOfParticles(); i++)
{
  std::cout << "Particle id: " << i << " has number of nbs: ";
  std::cout << myParticleSystem.getNeighborsCount( i );
  std::cout << " with id: [ ";
  for(int j = 0; j < myParticleSystem.getNeighborsCount( i ); j++)
    std::cout << myParticleSystem.getNeighbor(i, j) << " ";
  std::cout << "]"<< std::endl;
}

//myParticleSystem.saveNeighborList( "nbs.txt" );

std::ofstream pointsFile;
pointsFile.open("points.txt");
for(unsigned int p = 0; p < myParticleSystem.getNumberOfParticles(); p++)
{
  pointsFile << myParticleSystem.getPoint( p )[0] << "," << myParticleSystem.getPoint( p )[1] << std::endl;
}
pointsFile.close();

std::ofstream neighborsFile;
neighborsFile.open("neighbors.txt");
for(unsigned int i = 0; i < myParticleSystem.getNumberOfParticles(); i++)
{
  neighborsFile << myParticleSystem.getNeighborsCount( i ) << " ";
  for(int j = 0; j < myParticleSystem.getNeighborsCount( i ); j++)
    neighborsFile << myParticleSystem.getNeighbor(i, j) << " ";
  neighborsFile << std::endl;
}
neighborsFile.close();


