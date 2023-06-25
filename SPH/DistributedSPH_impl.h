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

//Functions for updating the domain ranges
template< typename SPHSimulation >
template< typename SPHObjectPointer >
typename DistributedSPHSimpleFluid< SPHSimulation >::PairIndexType
DistributedSPHSimpleFluid< SPHSimulation >::getFirstLastParticleInColumnOfCells( const GlobalIndexType& gridColumn,
                                                                                 const SPHObjectPointer& sphObject )
{
   const IndexVectorType gridSize = sphObject->particles->getGridSize();
   const GlobalIndexType indexOfFirstColumnCell = ParticleSystem::CellIndexer::EvaluateCellIndex(
         gridColumn, 1, gridSize );
   const GlobalIndexType indexOfLastColumnCell = ParticleSystem::CellIndexer::EvaluateCellIndex(
         gridColumn, gridSize[ 1 ] - 1, gridSize );
   const auto view_firstLastCellParticle = sphObject->particles->getCellFirstLastParticleList().getConstView(
         indexOfFirstColumnCell, indexOfLastColumnCell );

   auto fetch_vect = [=] __cuda_callable__ ( int i ) -> PairIndexType  { return view_firstLastCellParticle[ i ]; };
   auto reduction_vect = [=] __cuda_callable__ ( const PairIndexType& a, const PairIndexType& b ) -> PairIndexType
   { return { min( a[ 0 ], b[ 0 ] ), max( a[ 1 ], ( b[ 1 ] < INT_MAX ) ? b[ 1 ] : -1 ) }; };

   PairIndexType identity = { INT_MAX , INT_MIN };
   PairIndexType firstLastParticle = Algorithms::reduce< Devices::Cuda >(
         0, view_firstLastCellParticle.getSize(), fetch_vect, reduction_vect, identity );

   return firstLastParticle;
}

template< typename SPHSimulation >
template< typename SPHObjectPointer >
void
DistributedSPHSimpleFluid< SPHSimulation >::updateLocalSimulationInfo( SimulationSubdomainInfo& subdomainInfo,
                                                                       SPHObjectPointer& sphObject )
{
   const int rank = communicator.rank();
   const int nproc = communicator.size();

   GlobalIndexType gridColumnBegin = 0;
   GlobalIndexType gridColumnEnd = 0;

   if( rank == 0 )
   {
      gridColumnBegin = subdomainInfo.gridIdxEnd;

      const PairIndexType firstLastParticle = getFirstLastParticleInColumnOfCells( gridColumnBegin,
                                                                                   sphObject );
      subdomainInfo.firstParticleInLastGridColumn = firstLastParticle[ 0 ];
      subdomainInfo.lastParticleInLastGridColumn = firstLastParticle[ 1 ];
      subdomainInfo.numberOfParticlesToSendEnd = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

      sphObject->firstActiveParticle = subdomainInfo.firstParticleInFirstGridColumn;
      sphObject->lastActiveParticle = subdomainInfo.lastParticleInLastGridColumn + 1; //TODO: FIX!
   }

   //if( rank in between )

   if( rank == nproc - 1 )
   {
      gridColumnEnd = subdomainInfo.gridIdxBegin;

      const PairIndexType firstLastParticle = getFirstLastParticleInColumnOfCells( gridColumnEnd,
                                                                                   sphObject );
      subdomainInfo.firstParticleInFirstGridColumn = firstLastParticle[ 0 ];
      subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];
      subdomainInfo.numberOfParticlesToSendBegin = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

      //is this safe? -in case that arrangeRecievedAndLocalData updates number of particles, then yes
      subdomainInfo.lastParticleInLastGridColumn = sphObject->particles->getNumberOfParticles() - 1;

      sphObject->firstActiveParticle = subdomainInfo.firstParticleInFirstGridColumn;
      sphObject->lastActiveParticle = subdomainInfo.lastParticleInLastGridColumn + 1; //TODO: FIX!
   }

   //For load balancing
   //subdomainInfo.numberOfParticlesInThisSubdomain = sphObject->particles->getNumberOfParticles();
}

//Functions for data transfer
template< typename SPHSimulation >
template< typename Array >
void
DistributedSPHSimpleFluid< SPHSimulation >::synchronizeArray( Array& arraySend, Array& arrayReceive, SimulationSubdomainInfo& subdomainInfo, int valuesPerElement )
{
   static_assert( std::is_same< typename Array::DeviceType, DeviceType >::value, "mismatched DeviceType of the array" );
   using ValueType = typename Array::ValueType;

   ByteArrayView viewSend;
   viewSend.bind( reinterpret_cast< std::uint8_t* >( arraySend.getData() ), sizeof( ValueType ) * arraySend.getSize() );
   ByteArrayView viewReceive;
   viewReceive.bind( reinterpret_cast< std::uint8_t* >( arrayReceive.getData() ), sizeof( ValueType ) * arrayReceive.getSize() );

   synchronizeByteArray( viewSend, viewReceive, subdomainInfo, sizeof( ValueType ) * valuesPerElement );
}


template< typename SPHSimulation >
DistributedSPHSimpleFluid< SPHSimulation >::RequestsVector
DistributedSPHSimpleFluid< SPHSimulation >::synchronizeByteArrayAsyncWorker( ByteArrayView arraySend, ByteArrayView arrayReceive, SimulationSubdomainInfo& subdomainInfo, int bytesPerValue )
{
   //REAL FUNCTIONS
   const int rank = communicator.rank();
   const int nproc = communicator.size();
   const int maxParticlesToSend = 2000; //TODO: Make this general.

   // buffer for asynchronous communication requests
   RequestsVector requests;

   if( rank == 0 )
   {
      //Recieve
      subdomainInfo.receivedBegin = 0;
      requests.push_back( MPI::Irecv( &subdomainInfo.receivedEnd,
                                      1, //count
                                      1, //destination
                                      0,
                                      communicator ) );


      const GlobalIndexType receiveToPosition = subdomainInfo.lastParticleInLastGridColumn + 1;
      requests.push_back( MPI::Irecv( arrayReceive.getData() + bytesPerValue * receiveToPosition,
                                      bytesPerValue * maxParticlesToSend,
                                      1,
                                      0,
                                      communicator ) );

      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendEnd,
                                      1, //count
                                      1, //denstination
                                      0,
                                      communicator ) );

      const GlobalIndexType sendFromPosition =  subdomainInfo.firstParticleInLastGridColumn;
      requests.push_back( MPI::Isend( arraySend.getData() +  bytesPerValue * sendFromPosition,
                                      bytesPerValue * subdomainInfo.numberOfParticlesToSendEnd,
                                      1,
                                      0,
                                      communicator ) );
   }

   if( rank == nproc - 1 )
   {
      //Recieve
      subdomainInfo.receivedEnd = 0;
      requests.push_back( MPI::Irecv( &subdomainInfo.receivedBegin,
                                      1, //count
                                      nproc - 2, //destination
                                      0,
                                      communicator ) );

      requests.push_back( MPI::Irecv( arrayReceive.getData(),
                                      bytesPerValue * maxParticlesToSend,
                                      nproc - 2,
                                      0,
                                      communicator ) );

      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendBegin,
                                      1, //count
                                      nproc - 2, //destination
                                      0,
                                      communicator ) );

      const GlobalIndexType sendFromPosition = subdomainInfo.firstParticleInFirstGridColumn;
      requests.push_back( MPI::Isend( arraySend.getData() + bytesPerValue * sendFromPosition,
                                      bytesPerValue * subdomainInfo.numberOfParticlesToSendBegin,
                                      nproc - 2,
                                      0,
                                      communicator ) );
   }

   return requests;
}

template< typename SPHSimulation >
template< typename Array, typename SPHObjectPointer >
void
DistributedSPHSimpleFluid< SPHSimulation >::arrangeRecievedAndLocalData( Array& arraySend,
                                                                                                 Array& arrayReceive,
                                                                                                 SPHObjectPointer& sphObject,
                                                                                                 SimulationSubdomainInfo& subdomainInfo,
                                                                                                 bool tempSetNumberOfPtcs )
{
   //REAL FUNCTIONS
   const int rank = communicator.rank();
   const int nproc = communicator.size();

   auto arraySendView = arraySend.getView();
   auto arrayReceiveView = arrayReceive.getView();

   auto copyToSwap = [ = ] __cuda_callable__( GlobalIndexType i,
                                              GlobalIndexType offsetSend,
                                              GlobalIndexType offsetReceive ) mutable
   {
      arrayReceiveView[ offsetReceive + i ] = arraySendView[ offsetSend + i ];
   };

   if( rank == 0 )
   {
      GlobalIndexType offsetSend = 0;
      GlobalIndexType offsetReceive = 0;
      GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + 1;
      GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + subdomainInfo.receivedEnd; //EXPERIMENT

      //TODO: Remove this ugly aberattion.
      if( tempSetNumberOfPtcs == true )
         sphObject->particles->setNumberOfParticles( numberOfParticlesToSet );

      Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
   }

   if( rank == nproc - 1 )
   {
      GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
      GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
      GlobalIndexType numberOfParticlesToCopy = sphObject->particles->getNumberOfParticles() - subdomainInfo.firstParticleInFirstGridColumn;
      GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;

      //TODO: Remove this ugly aberattion.
      if( tempSetNumberOfPtcs == true )
         sphObject->particles->setNumberOfParticles( numberOfParticlesToSet );

      Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
   }

   arraySend.swap( arrayReceive );
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
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInNextSubdomain,
                                      1, //count
                                      1, //destination
                                      0,
                                      communicator ) );


      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain,
                                      1, //count
                                      1, //denstination
                                      0,
                                      communicator ) );
   }

   if( rank == nproc - 1 )
   {
      //Recieve
      requests.push_back( MPI::Irecv( &subdomainInfo.numberOfParticlesInPreviousSubdomain,
                                      1, //count
                                      nproc - 2, //destination
                                      0,
                                      communicator ) );

      //Send
      requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesInThisSubdomain,
                                      1, //count
                                      nproc - 2, //destination
                                      0,
                                      communicator ) );
   }

   return requests;
}

template< typename SPHSimulation >
void
DistributedSPHSimpleFluid< SPHSimulation >::updateSubdomainSize( SimulationSubdomainInfo& subdomainInfo, SimulationSubdomainInfo& subdomainInfo_boundary )
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
   }

   if( rank == 1 )
   {
      if( ( subdomainInfo.numberOfParticlesInThisSubdomain - subdomainInfo.numberOfParticlesInPreviousSubdomain ) > 500 )
      {
         subdomainInfo.gridIdxOverlapBegin++;
         subdomainInfo.gridIdxBegin++;

         subdomainInfo_boundary.gridIdxOverlapBegin++;
         subdomainInfo_boundary.gridIdxBegin++;
      }
      //oposite if statement
   }

}

//implement standard interaction functions aswell for distributed
template< typename SPHSimulation >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
void
DistributedSPHSimpleFluid< SPHSimulation >::interact( SPHState& sphState )
{
   localSimulation.template Interact< SPHKernelFunction, DiffusiveTerm, ViscousTerm, EOS >( sphState );
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

   logger.writeParameter( "Grid start cell index: ",
                           this->localSimulationInfo.gridIdxBegin );
   logger.writeParameter( "Grid end cell index: ",
                           this->localSimulationInfo.gridIdxEnd );
   logger.writeParameter( "Grid real start cell index (start-overlap): ",
                           this->localSimulationInfo.gridIdxOverlapBegin );
   logger.writeParameter( "Grid real end cell index (end-overlap): ",
                           this->localSimulationInfo.gridIdxOverlapEnd );

   logger.writeParameter( "var: {firstParticleInFirstGridColumn} ",
                           this->localSimulationInfo.firstParticleInFirstGridColumn );
   logger.writeParameter( "var: {lastParticleInFirstGridColumn} ",
                           this->localSimulationInfo.lastParticleInFirstGridColumn );
   logger.writeParameter( "var: {firstParticleInLastGridColumn} ",
                           this->localSimulationInfo.firstParticleInLastGridColumn );
   logger.writeParameter( "var: {lastParticleInLastGridColumn} ",
                           this->localSimulationInfo.lastParticleInLastGridColumn );

   logger.writeParameter( "var: {recievedStart} ",
                           this->localSimulationInfo.receivedBegin );
   logger.writeParameter( "var: {recievedEnd} ",
                           this->localSimulationInfo.receivedEnd );

   //Boundary
   logger.writeParameter( "var: {firstParticleInFirstGridColumn-boundary} ",
                           this->localSimulationInfo_boundary.firstParticleInFirstGridColumn );
   logger.writeParameter( "var: {lastParticleInFirstGridColumn-boundary} ",
                           this->localSimulationInfo_boundary.lastParticleInFirstGridColumn );
   logger.writeParameter( "var: {firstParticleInLastGridColumn-boundary} ",
                           this->localSimulationInfo_boundary.firstParticleInLastGridColumn );
   logger.writeParameter( "var: {lastParticleInLastGridColumn-boundary} ",
                           this->localSimulationInfo_boundary.lastParticleInLastGridColumn );

   logger.writeParameter( "var: {recievedStart-boundary} ",
                           this->localSimulationInfo_boundary.receivedBegin );
   logger.writeParameter( "var: {recievedEnd-boundary} ",
                           this->localSimulationInfo_boundary.receivedEnd );
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

