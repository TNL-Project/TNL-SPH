
for(int i = 0; i < myParticleSystem.getNumberOfParticles(); i++)
{
  std::cout << "Particle id: " << i << " has number of nbs: ";
  std::cout << myParticleSystem.getNeighborsCount( i );
  //for(int i = 0; i < this->)
  std::cout << std::endl;
}

myParticleSystem.saveNeighborList( "nbs.txt" );
