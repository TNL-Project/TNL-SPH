#pragma once

#include "../Particles/ParticlesTraits.h"
#include "SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleConfig >
class DistributedPhysicalObjectInfo
{
   public:
   using ParticleTraitsType = ParticlesTraits< ParticleConfig, typename ParticleConfig::DeviceType >;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using RealType = typename ParticleTraitsType::RealType;
   using PointType = typename ParticleTraitsType::PointType;

   /**
    * Store the index of grid column (in 2D) or slince (in 3D), which
    * represents the origin and end of the current subdomain with respect
    * to the global grid.
    */
   GlobalIndexType gridIdxOverlapBegin;
   GlobalIndexType gridIdxBegin;
   GlobalIndexType gridIdxOverlapEnd;
   GlobalIndexType gridIdxEnd;

   /**
    * Store index of first and last particle of first and last
    * grid column (in 2D) or slice (in 3D) of the subdomain.
    */
   GlobalIndexType firstParticleInFirstGridColumn = 0;
   GlobalIndexType lastParticleInFirstGridColumn = 0;
   GlobalIndexType firstParticleInLastGridColumn = 0;
   GlobalIndexType lastParticleInLastGridColumn = 0;

   /**
    * Number of particles to send from current rank to its
    * left subdomain (numberOfParticlesToSet) and right
    * subdomain (numberOfParticlesToSendEnd).
    */
   GlobalIndexType numberOfParticlesToSendBegin = 0;
   GlobalIndexType numberOfParticlesToSendEnd = 0;

   /**
    * Number of particles received on current rank from its
    * left subdomain (receivedBegin) and right subdomain (receivedEnd).
    * This is necessary to incorporte the recieved data to
    * arrays that exits in the subdomain.
    */
   GlobalIndexType receivedBegin = 0;
   GlobalIndexType receivedEnd = 0;

   /**
    * Load balancing metadata - domainds are resized based on the
    * computation time on current, left and right subdomain and number of
    * particles in current, left and right subdomain.
    */
   GlobalIndexType numberOfParticlesInPreviousSubdomain;
   GlobalIndexType numberOfParticlesInNextSubdomain;
   GlobalIndexType numberOfParticlesInThisSubdomain;

   RealType solutionTimeInPreviousSubdomain;
   RealType solutionTimeInNextSubdomain;
   RealType solutionTimeInThisSubdomain;

   template< typename SubdomainInfoParams >
   void loadParameters( SubdomainInfoParams params )
   {
      this->gridIdxOverlapBegin = params.gridIdxOverlapStar;
      this->gridIdxBegin = params.gridIdxStart;

      this->gridIdxOverlapEnd = params.gridIdxOverlapEnd;
      this->gridIdxEnd = params.gridIdxEnd;

      //----- debug ------------------------------------------------------
      std::cout << "[ DistributedPhysicalObjectInfo, rank: " << TNL::MPI::GetRank() << " ] Variables loaded." << std::endl;
      //----- end-debug --------------------------------------------------
   }

   void
   debugWriteProlog( TNL::Logger& logger ) const noexcept
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
};

template< typename ParticleConfig, typename SPHConfig >
class DistributedSPHSimulationSynchronizer
{
public:

   //using SPHConfig = typename SPHSimulationType::ModelType::SPHConfig;
   using SPHTriatsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;
   using GlobalIndexType = typename SPHTriatsType::GlobalIndexType;
   using ByteArrayView = Containers::ArrayView< std::uint8_t, DeviceType, GlobalIndexType >;
   using RequestsVector = std::vector< MPI_Request >;

   using ParticleTraitsType = ParticlesTraits< ParticleConfig, DeviceType >;
   using PairIndexType = typename ParticleTraitsType::PairIndexType;
   using SimulationSubdomainInfo = DistributedPhysicalObjectInfo< ParticleConfig >; //RESOLVE

   template< typename SPHObjectPointer >
   void
   updateLocalSimulationInfo( SPHObjectPointer& sphObject )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();

      GlobalIndexType gridColumnBegin = 0;
      GlobalIndexType gridColumnEnd = 0;

      if( rank == 0 )
      {
         gridColumnEnd = sphObject->subdomainInfo.gridIdxEnd;

         PairIndexType firstLastParticle;
         if constexpr( SPHConfig::spaceDimension == 2 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInColumnOfCells( gridColumnEnd );
         else if constexpr( SPHConfig::spaceDimension == 3 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInBlockOfCells( gridColumnEnd );

         //added
         sphObject->subdomainInfo.firstParticleInFirstGridColumn = sphObject->getFirstActiveParticle();
         //sphObject->subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];

         sphObject->subdomainInfo.firstParticleInLastGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInLastGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendEnd = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         sphObject->firstActiveParticle = sphObject->subdomainInfo.firstParticleInFirstGridColumn;
         sphObject->lastActiveParticle = sphObject->subdomainInfo.lastParticleInLastGridColumn;
      }

      if( ( rank > 0 ) && ( rank < nproc - 1) )
      {
         //begin
         gridColumnBegin = sphObject->subdomainInfo.gridIdxBegin;

         PairIndexType firstLastParticle;
         if constexpr( SPHConfig::spaceDimension == 2 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInColumnOfCells( gridColumnBegin );
         else if constexpr( SPHConfig::spaceDimension == 3 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInBlockOfCells( gridColumnBegin );

         sphObject->subdomainInfo.firstParticleInFirstGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendBegin = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         sphObject->firstActiveParticle = sphObject->subdomainInfo.firstParticleInFirstGridColumn;

         //end
         gridColumnEnd = sphObject->subdomainInfo.gridIdxEnd;

         //turn off: PairIndexType firstLastParticle;
         if constexpr( SPHConfig::spaceDimension == 2 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInColumnOfCells( gridColumnEnd );
         else if constexpr( SPHConfig::spaceDimension == 3 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInBlockOfCells( gridColumnEnd );

         sphObject->subdomainInfo.firstParticleInLastGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInLastGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendEnd = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         sphObject->lastActiveParticle = sphObject->subdomainInfo.lastParticleInLastGridColumn;
      }

      if( rank == nproc - 1 )
      {
         gridColumnBegin = sphObject->subdomainInfo.gridIdxBegin;

         PairIndexType firstLastParticle;
         if constexpr( SPHConfig::spaceDimension == 2 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInColumnOfCells( gridColumnBegin );
         else if constexpr( SPHConfig::spaceDimension == 3 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInBlockOfCells( gridColumnBegin );

         sphObject->subdomainInfo.firstParticleInFirstGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendBegin = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         //is this safe? -in case that arrangeRecievedAndLocalData updates number of particles, then yes
         //TURN OFF: sphObject->subdomainInfo.lastParticleInLastGridColumn = sphObject->particles->getNumberOfParticles() - 1;
         sphObject->subdomainInfo.lastParticleInLastGridColumn = sphObject->getLastActiveParticle();

         sphObject->firstActiveParticle = sphObject->subdomainInfo.firstParticleInFirstGridColumn;
         sphObject->lastActiveParticle = sphObject->subdomainInfo.lastParticleInLastGridColumn;
      }

      //For load balancing
      //subdomainInfo.numberOfParticlesInThisSubdomain = sphObject->particles->getNumberOfParticles();
      synchronizeSubdomainInfo( sphObject->subdomainInfo );
   }

   template< typename Array >
   void
   synchronizeArray( Array& arraySend, Array& arrayReceive, SimulationSubdomainInfo& subdomainInfo, int valuesPerElement = 1 )
   {
      static_assert( std::is_same< typename Array::DeviceType, DeviceType >::value, "mismatched DeviceType of the array" );
      using ValueType = typename Array::ValueType;

      ByteArrayView view;
      view.bind( reinterpret_cast< std::uint8_t* >( arraySend.getData() ), sizeof( ValueType ) * arraySend.getSize() );
      synchronizeByteArray( view, subdomainInfo, sizeof( ValueType ) * valuesPerElement );
   }

   void
   synchronizeByteArray( ByteArrayView array,
                         SimulationSubdomainInfo& subdomainInfo,
                         int bytesPerValue )
   {
      auto requests = synchronizeByteArrayAsyncWorker( array, subdomainInfo, bytesPerValue );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView array,
                                    SimulationSubdomainInfo& subdomainInfo,
                                    int bytesPerValue )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();
      const int maxParticlesToSend = 15000; //TODO: Make this general.

      // buffer for asynchronous communication requests
      RequestsVector requests;

      if( rank == 0 )
      {
         //Recieve
         const GlobalIndexType receiveToPosition = subdomainInfo.lastParticleInLastGridColumn + 1;
         requests.push_back( MPI::Irecv( array.getData() + bytesPerValue * receiveToPosition,
                                         bytesPerValue * maxParticlesToSend,
                                         1,
                                         0,
                                         communicator ) );
         //Send
         const GlobalIndexType sendFromPosition =  subdomainInfo.firstParticleInLastGridColumn;
         requests.push_back( MPI::Isend( array.getData() +  bytesPerValue * sendFromPosition,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendEnd,
                                         1,
                                         0,
                                         communicator ) );
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         //End - Recieve
         const GlobalIndexType receiveToPositionEnd = subdomainInfo.lastParticleInLastGridColumn + 1;
         requests.push_back( MPI::Irecv( array.getData() + bytesPerValue * receiveToPositionEnd,
                                         //bytesPerValue * maxParticlesToSend,
                                         bytesPerValue * subdomainInfo.receivedEnd,
                                         rank + 1,
                                         0,
                                         communicator ) );


         //Begin - Recieve
         const GlobalIndexType receiveToPositionBegin = subdomainInfo.firstParticleInFirstGridColumn - subdomainInfo.receivedBegin;
         requests.push_back( MPI::Irecv( array.getData() + bytesPerValue * receiveToPositionBegin,
                                         //bytesPerValue * maxParticlesToSend,
                                         bytesPerValue * subdomainInfo.receivedBegin,
                                         rank - 1,
                                         0,
                                         communicator ) );

         //End - Send
         const GlobalIndexType sendFromPositionEnd =  subdomainInfo.firstParticleInLastGridColumn;
         requests.push_back( MPI::Isend( array.getData() +  bytesPerValue * sendFromPositionEnd,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendEnd,
                                         rank + 1,
                                         0,
                                         communicator ) );

         //Begin -Send
         const GlobalIndexType sendFromPositionBegin = subdomainInfo.firstParticleInFirstGridColumn;
         requests.push_back( MPI::Isend( array.getData() + bytesPerValue * sendFromPositionBegin,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendBegin,
                                         rank - 1,
                                         0,
                                         communicator ) );
      }

      if( rank == nproc - 1 )
      {
         //Recieve
         const GlobalIndexType receiveToPosition = subdomainInfo.firstParticleInFirstGridColumn - subdomainInfo.receivedBegin;
         requests.push_back( MPI::Irecv( array.getData() + bytesPerValue * receiveToPosition,
                                         //bytesPerValue * maxParticlesToSend,
                                         bytesPerValue * subdomainInfo.receivedBegin,
                                         nproc - 2,
                                         0,
                                         communicator ) );
         //Send
         const GlobalIndexType sendFromPosition = subdomainInfo.firstParticleInFirstGridColumn;
         requests.push_back( MPI::Isend( array.getData() + bytesPerValue * sendFromPosition,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendBegin,
                                         nproc - 2,
                                         0,
                                         communicator ) );
      }

      return requests;
   }

   void
   synchronizeSubdomainInfo( SimulationSubdomainInfo& subdomainInfo )
   {
      auto requests = synchronizeSubdomainInfoAsyncWorker( subdomainInfo );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeSubdomainInfoAsyncWorker( SimulationSubdomainInfo& subdomainInfo )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();

      // buffer for asynchronous communication requests
      RequestsVector requests;

      if( rank == 0 )
      {
         //Recieve
         subdomainInfo.receivedBegin = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedEnd, 1, 1, 0, communicator ) );
         //Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendEnd, 1, 1, 0, communicator ) );
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         //End - Recieve
         subdomainInfo.receivedEnd = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedEnd, 1, rank + 1, 0, communicator ) );
         //Begin - Recieve
         subdomainInfo.receivedBegin = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedBegin, 1, rank - 1, 0, communicator ) );
         //End - Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendEnd, 1, rank + 1, 0, communicator ) );
         //Begin - Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendBegin, 1, rank - 1, 0, communicator ) );
      }

      if( rank == nproc - 1 )
      {
         //Recieve
         subdomainInfo.receivedEnd = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedBegin, 1, nproc - 2, 0, communicator ) );
         //Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendBegin, 1, nproc - 2, 0, communicator ) );
      }

      return requests;
   }

   template< typename Array >
   void
   arrangeRecievedAndLocalData( Array& arraySend,
                                Array& arrayReceive,
                                SimulationSubdomainInfo& subdomainInfo )
   {
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
         const GlobalIndexType offsetSend = 0;
         const GlobalIndexType offsetReceive = 0;
         const GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + 1;

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         const GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         const GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         const GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn - \
            subdomainInfo.firstParticleInFirstGridColumn + 1 + subdomainInfo.receivedEnd;

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
      }

      if( rank == nproc - 1 )
      {
         const GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         const GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         const GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + \
            1 - subdomainInfo.firstParticleInFirstGridColumn;

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
      }

      arraySend.swap( arrayReceive );
   }

   template< typename SimulationSubdomainInfo >
   GlobalIndexType
   getNewParticleCount( SimulationSubdomainInfo& subdomainInfo, const GlobalIndexType currentNumberOfParticles )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();

      GlobalIndexType numberOfParticlesToSet;

      if( rank == 0 )
      {
         numberOfParticlesToSet = subdomainInfo.lastParticleInLastGridColumn + subdomainInfo.receivedEnd + 1;
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         numberOfParticlesToSet = subdomainInfo.lastParticleInLastGridColumn - subdomainInfo.firstParticleInFirstGridColumn +
                                  subdomainInfo.receivedEnd + subdomainInfo.receivedBegin + 1;
      }

      if( rank == nproc - 1 )
      {
         numberOfParticlesToSet = subdomainInfo.lastParticleInLastGridColumn - subdomainInfo.firstParticleInFirstGridColumn +
                                  subdomainInfo.receivedBegin + 1;
      }

      return numberOfParticlesToSet;
   }

   MPI::Comm communicator = MPI_COMM_WORLD;
   //const int maxParticlesToSend = 15000;

};




} // SPH
} // ParticleSystem
} // TNL

