
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


std::ofstream outputFile;
outputFile.open("points.txt");

for(unsigned int p = 0; p < myParticleSystem.getNumberOfParticles(); p++)
{
  outputFile << myParticleSystem.getPoint( p )[0] << "," << myParticleSystem.getPoint( p )[1] << std::endl;
}

outputFile.close();
