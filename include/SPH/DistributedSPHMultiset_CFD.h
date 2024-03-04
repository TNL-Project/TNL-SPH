#pragma once

namespace TNL {
namespace SPH {

template< typename SimulationType >
class DistributedSPHMultiset_CFD
{
public:
   using ModelType = typename SimulationType::Model;
   using SPHConfig = typename SimulationType::SPHConfig;
   using ParticlesType = typename ModelType::ParticlesType;;

   using OpenBoundary = typename SimulationType::OpenBoundaryBuffers;
   using OpenBoundaryPointer = typename SimulationType::OpenBoundaryPointer;


void
init( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger )
{

}
/*
   Synchronizer needs to:
   0. transform particles which moved into another subdomains
   1. set local dimensions of subdomain
   2. set ghost zones for ghost patches

      a) if subdomain rebalanced -> synchronizer build the sendDimension, sendBegin and recvBegin and sendSizes.
         if not -> skip this point

         ghostPatch[ nbs ].setCells( sendBegin, sendBegin + sendDimension )

      b) update ghotstPatches

         for( int i = 0; i < this->getNeighborsCount(); i++ ) {

   3. send / move the particles
 */
void
setOverlaps()
{
   synchronizer.setDistributedGrid( distributedGrid, ghostPatchesBoundary );
}

void
synchronize()
{
   for( int i = 0; i < getNeighborsCount(); i++ ){
      ghostPatchesFluid[ i ]->zone.updateParticlesInZone( localSimulation.fluid->particles );
      ghostPatchesBoundary[ i ]->zone.updateParticlesInZone( localSimulation.boundary->particles );
   }

   synchronizer.synchronizeOverlapSizes();
   localSimulation.fluid->synchronize( synchronizer, ghostPatchesFluid );
   localSimulation.boundary->synchronize( synchronizer, ghostPatchesFluid );
}

void
moveParticles();

public:

   SimulationType localSimulation;
   DistributedGridType distributedGrid;

   SimulationSynchronzier synchronizer;

   std::vector< OpenBoundaryPointer > ghostPatchesFluid( getNeighborsCount() );
   std::vector< OpenBoundaryPointer > ghostPatchesBoundary( getNeighborsCount() );


};

} // SPH
} // TNL

