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
   GlobalIndexType firstParticleInFirstGridColumn;
   GlobalIndexType lastParticleInFirstGridColumn;
   GlobalIndexType firstParticleInLastGridColumn;
   GlobalIndexType lastParticleInLastGridColumn;

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

      //Local simulation info:
      this->firstParticleInFirstGridColumn = 0;
      this->lastParticleInFirstGridColumn = 0;
      this->firstParticleInLastGridColumn = 0;
      this->lastParticleInLastGridColumn = 0;
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept;
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

         std::cout << "[ updateLocalSimulationInfo ] first and last particle for rank: " << rank << " in column-end: " << gridColumnEnd << " is: " <<  firstLastParticle << std::endl;

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

         std::cout << "[ updateLocalSimulationInfo ] first and last particle for rank: " << rank <<  " in column-begin: " << gridColumnBegin << " is: " <<  firstLastParticle << std::endl;

         sphObject->subdomainInfo.firstParticleInFirstGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendBegin = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         //is this safe? -in case that arrangeRecievedAndLocalData updates number of particles, then yes
         //turn off: subdomainInfo.lastParticleInLastGridColumn = sphObject->particles->getNumberOfParticles() - 1;

         sphObject->firstActiveParticle = sphObject->subdomainInfo.firstParticleInFirstGridColumn;
         //turn off: sphObject->lastActiveParticle = subdomainInfo.lastParticleInLastGridColumn;

         //end
         gridColumnEnd = sphObject->subdomainInfo.gridIdxEnd;

         //turn off: PairIndexType firstLastParticle;
         if constexpr( SPHConfig::spaceDimension == 2 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInColumnOfCells( gridColumnEnd );
         else if constexpr( SPHConfig::spaceDimension == 3 )
            firstLastParticle = sphObject->particles->getFirstLastParticleInBlockOfCells( gridColumnEnd );

         std::cout << "[ updateLocalSimulationInfo ] first and last particle for rank: " << rank << " in column-end: " << gridColumnEnd <<  " is: " <<  firstLastParticle << std::endl;

         sphObject->subdomainInfo.firstParticleInLastGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInLastGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendEnd = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         //turn off: sphObject->firstActiveParticle = subdomainInfo.firstParticleInFirstGridColumn;
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

         std::cout << "[ updateLocalSimulationInfo ] first and last particle for rank: " << rank <<  " in column-begin: " << gridColumnBegin <<   " is: " <<  firstLastParticle << std::endl;

         sphObject->subdomainInfo.firstParticleInFirstGridColumn = firstLastParticle[ 0 ];
         sphObject->subdomainInfo.lastParticleInFirstGridColumn = firstLastParticle[ 1 ];
         sphObject->subdomainInfo.numberOfParticlesToSendBegin = firstLastParticle[ 1 ] - firstLastParticle[ 0 ] + 1;

         //is this safe? -in case that arrangeRecievedAndLocalData updates number of particles, then yes
         sphObject->subdomainInfo.lastParticleInLastGridColumn = sphObject->particles->getNumberOfParticles() - 1;

         sphObject->firstActiveParticle = sphObject->subdomainInfo.firstParticleInFirstGridColumn;
         sphObject->lastActiveParticle = sphObject->subdomainInfo.lastParticleInLastGridColumn;
      }

      //For load balancing
      //subdomainInfo.numberOfParticlesInThisSubdomain = sphObject->particles->getNumberOfParticles();
   }

   template< typename Array >
   void
   synchronizeArray( Array& arraySend, Array& arrayReceive, SimulationSubdomainInfo& subdomainInfo, int valuesPerElement = 1 )
   {
      static_assert( std::is_same< typename Array::DeviceType, DeviceType >::value, "mismatched DeviceType of the array" );
      using ValueType = typename Array::ValueType;

      ByteArrayView viewSend;
      viewSend.bind( reinterpret_cast< std::uint8_t* >( arraySend.getData() ), sizeof( ValueType ) * arraySend.getSize() );
      ByteArrayView viewReceive;
      viewReceive.bind( reinterpret_cast< std::uint8_t* >( arrayReceive.getData() ), sizeof( ValueType ) * arrayReceive.getSize() );

      synchronizeByteArray( viewSend, viewReceive, subdomainInfo, sizeof( ValueType ) * valuesPerElement );
   }

   void
   synchronizeByteArray( ByteArrayView arraySend,
                         ByteArrayView arrayReceive,
                         SimulationSubdomainInfo& subdomainInfo,
                         int bytesPerValue )
   {
      auto requests = synchronizeByteArrayAsyncWorker( arraySend, arrayReceive, subdomainInfo, bytesPerValue );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView arraySend,
                                    ByteArrayView arrayReceive,
                                    SimulationSubdomainInfo& subdomainInfo,
                                    int bytesPerValue )
   {
      //REAL FUNCTIONS
      const int rank = communicator.rank();
      const int nproc = communicator.size();
      const int maxParticlesToSend = 15000; //TODO: Make this general.

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

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         //End
         //Recieve
         //turn off: subdomainInfo.receivedBegin = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedEnd,
                                         1, //count
                                         rank + 1, //destination
                                         0,
                                         communicator ) );


         //FIXME: ugly ugly workaround
         //const GlobalIndexType receiveToPositionWorkaround = subdomainInfo.lastParticleInLastGridColumn + 1;
         const GlobalIndexType receiveToPosition = subdomainInfo.lastParticleInLastGridColumn + 1;
         //requests.push_back( MPI::Irecv( arrayReceive.getData() + bytesPerValue * receiveToPosition, //FIXME: Workaround arrayReceiveReplaced
         requests.push_back( MPI::Irecv( arraySend.getData() + bytesPerValue * receiveToPosition, //FIXME: Workaround arrayReceiveReplaced
                                         bytesPerValue * maxParticlesToSend,
                                         rank + 1,
                                         0,
                                         communicator ) );


         //Begin
         //Recieve
         //turn off: subdomainInfo.receivedEnd = 0;
         requests.push_back( MPI::Irecv( &subdomainInfo.receivedBegin,
                                         1, //count
                                         rank - 1, //destination
                                         0,
                                         communicator ) );

         requests.push_back( MPI::Irecv( arrayReceive.getData(),
                                         bytesPerValue * maxParticlesToSend,
                                         rank - 1,
                                         0,
                                         communicator ) );

         //End
         //Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendEnd,
                                         1, //count
                                         rank + 1, //denstination
                                         0,
                                         communicator ) );

         const GlobalIndexType sendFromPositionEnd =  subdomainInfo.firstParticleInLastGridColumn;
         requests.push_back( MPI::Isend( arraySend.getData() +  bytesPerValue * sendFromPositionEnd,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendEnd,
                                         rank + 1,
                                         0,
                                         communicator ) );

         //Begin
         //Send
         requests.push_back( MPI::Isend( &subdomainInfo.numberOfParticlesToSendBegin,
                                         1, //count
                                         rank - 1, //destination
                                         0,
                                         communicator ) );

         const GlobalIndexType sendFromPositionBegin = subdomainInfo.firstParticleInFirstGridColumn;
         requests.push_back( MPI::Isend( arraySend.getData() + bytesPerValue * sendFromPositionBegin,
                                         bytesPerValue * subdomainInfo.numberOfParticlesToSendBegin,
                                         rank - 1,
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
         GlobalIndexType offsetSend = 0;
         GlobalIndexType offsetReceive = 0;
         GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + 1;
         GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + subdomainInfo.receivedEnd;

         //: //TODO: Remove this ugly aberattion.
         //: if( tempSetNumberOfPtcs == true ){
         //:    sphObject->particles->setNumberOfParticles( numberOfParticlesToSet );
         //:    sphObject->particles->setLastActiveParticle( numberOfParticlesToSet - 1 );
         //: }

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {


         //Begin
         GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn - subdomainInfo.firstParticleInFirstGridColumn + 1 + subdomainInfo.receivedEnd;
         GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;

         //: //TODO: Remove this ugly aberattion.
         //: if( tempSetNumberOfPtcs == true ){
         //:    sphObject->particles->setNumberOfParticles( numberOfParticlesToSet );
         //:    sphObject->particles->setLastActiveParticle( numberOfParticlesToSet - 1 );
         //: }

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );

      }

      if( rank == nproc - 1 )
      {
         GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         //GlobalIndexType numberOfParticlesToCopy = sphObject->particles->getNumberOfParticles() - subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + 1 - subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;

         //: //TODO: Remove this ugly aberattion.
         //: if( tempSetNumberOfPtcs == true ){
         //:    sphObject->particles->setNumberOfParticles( numberOfParticlesToSet );
         //:    sphObject->particles->setLastActiveParticle( numberOfParticlesToSet - 1 );
         //: }

         Algorithms::parallelFor< DeviceType >( 0, numberOfParticlesToCopy, copyToSwap, offsetSend, offsetReceive );
      }

      arraySend.swap( arrayReceive );
   }

   template< typename SimulationSubdomainInfo >
   GlobalIndexType
   getNumberOfParticlesAfterSynchronization( SimulationSubdomainInfo& subdomainInfo, const GlobalIndexType currentNumberOfParticles )
   {
      const int rank = communicator.rank();
      const int nproc = communicator.size();

      GlobalIndexType numberOfParticlesToSet;

      if( rank == 0 )
      {
         GlobalIndexType offsetSend = 0;
         GlobalIndexType offsetReceive = 0;
         GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn + 1;
         numberOfParticlesToSet = numberOfParticlesToCopy + subdomainInfo.receivedEnd;
         //GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + subdomainInfo.receivedEnd;
         //return numberOfParticlesToSet;
      }

      if( ( rank > 0 ) && ( rank < nproc - 1 ) )
      {
         GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         GlobalIndexType numberOfParticlesToCopy = subdomainInfo.lastParticleInLastGridColumn - subdomainInfo.firstParticleInFirstGridColumn + 1 + subdomainInfo.receivedEnd;
         numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;
         //GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;
         //return numberOfParticlesToSet;
      }

      if( rank == nproc - 1 )
      {
         GlobalIndexType offsetSend = subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType offsetReceive = subdomainInfo.receivedBegin;
         //GlobalIndexType numberOfParticlesToCopy = sphObject->particles->getNumberOfParticles() - subdomainInfo.firstParticleInFirstGridColumn;
         GlobalIndexType numberOfParticlesToCopy = currentNumberOfParticles - subdomainInfo.firstParticleInFirstGridColumn;
         numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;
         //GlobalIndexType numberOfParticlesToSet = numberOfParticlesToCopy + offsetReceive;
         //return numberOfParticlesToSet;
      }

      return numberOfParticlesToSet;
   }

      //DistributedPhysicalObjectInfo subdomainInfo;
   MPI::Comm communicator = MPI_COMM_WORLD;

};


} // SPH
} // ParticleSystem
} // TNL

