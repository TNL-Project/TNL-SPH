#include "DistributedSPH.h"
#include <TNL/Functional.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::setParticlesDecomposition( const IndexVectorType& domainDecomposition )
{
   this->domainDecomposition = domainDecomposition;
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::setCommunicator( const MPI::Comm& communicator )
{
   this->domainDecomposition = domainDecomposition;
}

template< typename SPHSimulation >
const typename DistributedSPHSimpleFluid< SPHSimulation >::IndexVectorType&
DistributedSPHSimpleFluid< SPHSimulation >::getParticleDimension() const
{
   return this->domainDecomposition;
}

//Load balancing
template< typename SPHSimulation >
DistributedSPHSimpleFluid< SPHSimulation >::RequestsVector
DistributedSPHSimpleFluid< SPHSimulation >::synchronizeSubdomainMetaDataArrayAsyncWorker( SimulationSubdomainInfo& subdomainInfo )
{
   //REAL FUNCTIONS
   const int rank = communicator.rank();
   const int nproc = communicator.size();

   // buffer for asynchronous communication requests
   RequestsVector requests;

   if( rank == 0 )
   {
      //Recieve
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInNextSubdomain, 1, 1, 0, communicator ) );
      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain, 1, 1, 0, communicator ) );
   }

   if( ( rank > 0 ) && ( rank < nproc - 1 ) )
   {
      //End - Recieve
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInNextSubdomain, 1, rank + 1, 0, communicator ) );
      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain, 1, rank + 1, 0, communicator ) );

      //Begin - Recieve
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInPreviousSubdomain, 1, rank - 1, 0, communicator ) );
      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain, 1, rank - 1, 0, communicator ) );
   }

   if( rank == nproc - 1 )
   {
      //Recieve
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInPreviousSubdomain, 1, nproc - 2, 0, communicator ) );
      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain, 1, nproc - 2, 0, communicator ) );
   }

   return requests;
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::updateSubdomainSize( SimulationSubdomainInfo& subdomainInfo,
                                                                 SimulationSubdomainInfo& subdomainInfo_boundary )
{

   //REAL FUNCTIONS
   const int rank = communicator.rank();
   const int nproc = communicator.size();

   // buffer for asynchronous communication requests
   RequestsVector requests;

   if( rank == 0 )
   {
      if( ( subdomainInfo.numberOfParticlesInNextSubdomain - subdomainInfo.numberOfParticlesInThisSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapEnd++;
         subdomainInfo.gridIdxEnd++;

         subdomainInfo_boundary.gridIdxOverlapEnd++;
         subdomainInfo_boundary.gridIdxEnd++;
      }
      //oposite if statement
      else if( ( subdomainInfo.numberOfParticlesInThisSubdomain - subdomainInfo.numberOfParticlesInNextSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapEnd--;
         subdomainInfo.gridIdxEnd--;

         subdomainInfo_boundary.gridIdxOverlapEnd--;
         subdomainInfo_boundary.gridIdxEnd--;
      }
   }

   if( rank == 1 )
   {
      //handle the begin
      if( ( subdomainInfo.numberOfParticlesInThisSubdomain - subdomainInfo.numberOfParticlesInPreviousSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapBegin++;
         subdomainInfo.gridIdxBegin++;

         subdomainInfo_boundary.gridIdxOverlapBegin++;
         subdomainInfo_boundary.gridIdxBegin++;
      }
      //oposite if statement
      else if( ( subdomainInfo.numberOfParticlesInPreviousSubdomain - subdomainInfo.numberOfParticlesInThisSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapBegin--;
         subdomainInfo.gridIdxBegin--;

         subdomainInfo_boundary.gridIdxOverlapBegin--;
         subdomainInfo_boundary.gridIdxBegin--;
      }

      //handle the end
      if( ( subdomainInfo.numberOfParticlesInNextSubdomain - subdomainInfo.numberOfParticlesInThisSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapEnd++;
         subdomainInfo.gridIdxEnd++;

         subdomainInfo_boundary.gridIdxOverlapEnd++;
         subdomainInfo_boundary.gridIdxEnd++;
      }
      //oposite if statement
      else if( (  subdomainInfo.numberOfParticlesInThisSubdomain - subdomainInfo.numberOfParticlesInNextSubdomain  ) > 500 )
      {
         subdomainInfo.gridIdxOverlapEnd--;
         subdomainInfo.gridIdxEnd--;

         subdomainInfo_boundary.gridIdxOverlapEnd--;
         subdomainInfo_boundary.gridIdxEnd--;
      }
   }

   if( rank == 2 )
   {
      if( ( subdomainInfo.numberOfParticlesInThisSubdomain - subdomainInfo.numberOfParticlesInPreviousSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapBegin++;
         subdomainInfo.gridIdxBegin++;

         subdomainInfo_boundary.gridIdxOverlapBegin++;
         subdomainInfo_boundary.gridIdxBegin++;
      }
      //oposite if statement
      else if( (  subdomainInfo.numberOfParticlesInPreviousSubdomain - subdomainInfo.numberOfParticlesInThisSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapBegin--;
         subdomainInfo.gridIdxBegin--;

         subdomainInfo_boundary.gridIdxOverlapBegin--;
         subdomainInfo_boundary.gridIdxBegin--;
      }
   }

}

//implement standard interaction functions aswell for distributed
template< typename SPHSimulation >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
void
DistributedSPHSimpleFluid< SPHSimulation >::interact( SPHState& sphState )
{
   localSimulation.template interact< SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >( sphState );
}

template< typename SPHSimulation >
template< typename Writer >
void
DistributedSPHSimpleFluid< SPHSimulation >::save( const std::string& outputFileName, const int step )
{
   const int rank = communicator.rank();
   std::string outputFileNameWithRank = outputFileName + "_rank" + std::to_string( rank ) + "_";
   localSimulation.template save< Writer >( outputFileNameWithRank, step );
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::performLoadBalancing()
{
   //Perform load balancing - the old version
   localSimulation.fluid->subdomainInfo.numberOfParticlesInThisSubdomain = localSimulation.fluid->particles->getNumberOfParticles();
   synchronizeSubdomainMetaData( localSimulation.fluid->subdomainInfo );
   TNL::MPI::Barrier( communicator );
   updateSubdomainSize( localSimulation.fluid->subdomainInfo, localSimulation.boundary->subdomainInfo );

   localSimulation.fluid->centerObjectArraysInMemory();
   localSimulation.boundary->centerObjectArraysInMemory();
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::updateLocalSubdomain()
{
   //Update information about subdomain - TODO: Not should be part of synchronizer
   synchronizer.updateLocalSimulationInfo( localSimulation.fluid );
   synchronizer.updateLocalSimulationInfo( localSimulation.boundary );
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::synchronize()
{
   localSimulation.fluid->synchronizeObject( synchronizer );
   localSimulation.boundary->synchronizeObject( synchronizer );

   localSimulation.fluid->completeSynchronization();
   localSimulation.boundary->completeSynchronization();
}


template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::writeProlog( TNL::Logger& logger ) const noexcept
{
   logger.writeParameter( "Number of fluid particles:",
                           this->localSimulation.fluid->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. fluid particles:",
                           this->localSimulation.fluid->particles->getNumberOfAllocatedParticles() );
   logger.writeParameter( "Number of boundary particles:",
                           this->localSimulation.boundary->particles->getNumberOfParticles() );
   logger.writeParameter( "Number of alloc. boundary particles:",
                           this->localSimulation.boundary->particles->getNumberOfAllocatedParticles() );

   //MOVE THIS INTO SYNCHRONIZER
   logger.writeParameter( "Grid start cell index: ",
                           this->localSimulation.fluid->subdomainInfo.gridIdxBegin );
   logger.writeParameter( "Grid end cell index: ",
                           this->localSimulation.fluid->subdomainInfo.gridIdxEnd );
   logger.writeParameter( "Grid real start cell index (start-overlap): ",
                           this->localSimulation.fluid->subdomainInfo.gridIdxOverlapBegin );
   logger.writeParameter( "Grid real end cell index (end-overlap): ",
                           this->localSimulation.fluid->subdomainInfo.gridIdxOverlapEnd );

   logger.writeParameter( "var: {firstParticleInFirstGridColumn} ",
                           this->localSimulation.fluid->subdomainInfo.firstParticleInFirstGridColumn );
   logger.writeParameter( "var: {lastParticleInFirstGridColumn} ",
                           this->localSimulation.fluid->subdomainInfo.lastParticleInFirstGridColumn );
   logger.writeParameter( "var: {firstParticleInLastGridColumn} ",
                           this->localSimulation.fluid->subdomainInfo.firstParticleInLastGridColumn );
   logger.writeParameter( "var: {lastParticleInLastGridColumn} ",
                           this->localSimulation.fluid->subdomainInfo.lastParticleInLastGridColumn );

   logger.writeParameter( "var: {recievedStart} ",
                           this->localSimulation.fluid->subdomainInfo.receivedBegin );
   logger.writeParameter( "var: {recievedEnd} ",
                           this->localSimulation.fluid->subdomainInfo.receivedEnd );

   //Boundary
   logger.writeParameter( "var: {firstParticleInFirstGridColumn-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.firstParticleInFirstGridColumn );
   logger.writeParameter( "var: {lastParticleInFirstGridColumn-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.lastParticleInFirstGridColumn );
   logger.writeParameter( "var: {firstParticleInLastGridColumn-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.firstParticleInLastGridColumn );
   logger.writeParameter( "var: {lastParticleInLastGridColumn-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.lastParticleInLastGridColumn );

   logger.writeParameter( "var: {recievedStart-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.receivedBegin );
   logger.writeParameter( "var: {recievedEnd-boundary} ",
                           this->localSimulation.boundary->subdomainInfo.receivedEnd );
}

//DEBUG
template< typename ParticleConfig >
void
SimulationSubdomainInfo< ParticleConfig >::writeProlog( TNL::Logger& logger ) const noexcept
{
   logger.writeParameter( "Grid start cell index: ",
                           this->gridIdxBegin );
   logger.writeParameter( "Grid end cell index: ",
                           this->gridIdxEnd );
   logger.writeParameter( "Grid real start cell index (start-overlap): ",
                           this->gridIdxOverlapBegin );
   logger.writeParameter( "Grid real end cell index (end-overlap): ",
                           this->gridIdxOverlapEnd );

   logger.writeParameter( "var: {firstParticleInFirstGridColumn} ",
                           this->firstParticleInFirstGridColumn );
   logger.writeParameter( "var: {lastParticleInFirstGridColumn} ",
                           this->lastParticleInFirstGridColumn );
   logger.writeParameter( "var: {firstParticleInLastGridColumn} ",
                           this->firstParticleInLastGridColumn );
   logger.writeParameter( "var: {lastParticleInLastGridColumn} ",
                           this->lastParticleInLastGridColumn );

   logger.writeParameter( "var: {numberOfParticlesToSendBegin} ",
                           this->numberOfParticlesToSendBegin );
   logger.writeParameter( "var: {numberOfParticlesToSendEnd} ",
                           this->numberOfParticlesToSendEnd );
   logger.writeParameter( "var: {recievedStart} ",
                           this->receivedBegin );
   logger.writeParameter( "var: {recievedEnd} ",
                           this->receivedEnd );
}

} // SPH
} // ParticleSystem
} // TNL

