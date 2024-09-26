auto pressure_view = sph.fluid->variables->p.getView();
auto pressureBoundary_view = sph.boundary->variables->p.getView();
auto nbList_view = sph.fluid->getParticles()->getNeighborListStorage().getView();
auto nbListBoundary_view = sph.boundary->getParticles()->getNeighborListStorage().getView();
const int neighborsCountLimit = sph.fluid->getParticles()->getNeighborsCountLimit();
const int numberOfFluidParticles = sph.fluid->getNumberOfParticles();
const int numberOfBoundaryParticles = sph.boundary->getNumberOfParticles();

auto particleLoop = [=] __cuda_callable__ ( int i ) mutable
{
   const int fluidNeighborsCount = nbList_view[ i ];
   const int boundaryNeighborsCount = nbList_view[ numberOfFluidParticles * ( fluidNeighborsCount + 1 ) + i ];
   pressure_view[ i ] = fluidNeighborsCount + boundaryNeighborsCount;
};
sph.fluid->particles->forAll( particleLoop );

auto particleLoopBoundary = [=] __cuda_callable__ ( int i ) mutable
{
   const int fluidNeighborsCount = nbListBoundary_view[ i ];
   pressureBoundary_view[ i ] = fluidNeighborsCount;
};
sph.boundary->particles->forAll( particleLoopBoundary );

