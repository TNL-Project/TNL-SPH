#pragma once

#include "../Particles/ParticlesTraits.h"
#include "SPHTraits.h"

namespace TNL {
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
   synchronizeByteArray( ByteArrayView array, SimulationSubdomainInfo& subdomainInfo, int bytesPerValue )
   {
      auto requests = synchronizeByteArrayAsyncWorker( array, subdomainInfo, bytesPerValue );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView array, SimulationSubdomainInfo& subdomainInfo, int bytesPerValue )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();

      //buffer for asynchronous communication requests
      RequestsVector requests;

      if( rank == 0 )
      {
         //Recieve
         const GlobalIndexType receiveToPosition = subdomainInfo.lastParticleInLastGridColumn + 1;
         requests.push_back( MPI::Irecv( array.getData() + bytesPerValue * receiveToPosition,
                                         //bytesPerValue * maxParticlesToSend,
                                         bytesPerValue * subdomainInfo.receivedEnd,
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

      //buffer for asynchronous communication requests
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

   //Load balancing
   void
   synchronizeSubdomainMetaData( SimulationSubdomainInfo& subdomainInfo )
   {
      auto requests = synchronizeSubdomainMetaDataArrayAsyncWorker( subdomainInfo );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeSubdomainMetaDataArrayAsyncWorker( SimulationSubdomainInfo& subdomainInfo )
   {
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

//protected:

   MPI::Comm communicator = MPI_COMM_WORLD;

};

} // SPH
} // TNL

