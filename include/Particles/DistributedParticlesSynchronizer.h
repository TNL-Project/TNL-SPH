#pragma once

#include "ParticlesTraits.h"
#include "GhostZone.h"
#include <vector>

#include <TNL/Containers/ByteArraySynchronizer.h>
#include <TNL/MPI/Comm.h>
#include <TNL/MPI/Wrappers.h>
#include <TNL/Algorithms/parallelFor.h>

namespace TNL {
namespace ParticleSystem {

template< typename DistributedParticlesType >
class DistributedParticlesSynchronizer
//: public Containers::ByteArraySynchronizer< typename Particles::DeviceType, typename Particles::GlobalIndexType >
{
public:

   using Particles = typename DistributedParticlesType::ParticleSystemType;
   using DistributedGridType = typename DistributedParticlesType::DistributedGridType;

   using RealType = typename Particles::RealType;
   using GlobalIndexType = typename Particles::GlobalIndexType;
   //using Index = typename DistributedGridType::Index;
   using Index = typename Particles::GlobalIndexType; //FIXME: Back to grid
   using ParticleZone = ParticleZone< typename Particles::Config >;
   using Device = typename Particles::DeviceType;
   using DeviceType = typename Particles::DeviceType;

   using CoordinatesType = typename DistributedGridType::CoordinatesType;
   using SubdomainOverlapsType = typename DistributedGridType::SubdomainOverlapsType;

   using ByteArrayView = Containers::ArrayView< std::uint8_t, Device, GlobalIndexType >;

   //using Base = Containers::ByteArraySynchronizer< typename Particles::DeviceType, typename Particles::GlobalIndexType >;
   //using RequestsVector = typename Base::RequestsVector;
   //using ByteArrayView = Containers::ArrayView< std::uint8_t, Device, GlobalIndexType >;
   using RequestsVector = std::vector< MPI_Request >;

   using DistributedParticlesPointerType = typename Pointers::SharedPointer< DistributedParticlesType >;

   //done
   [[nodiscard]] static constexpr int
   getMeshDimension()
   {
      return DistributedGridType::getMeshDimension();
   }

   //done
   [[nodiscard]] static constexpr int
   getNeighborsCount()
   {
      return DistributedGridType::getNeighborsCount();
   }

   //done
   DistributedParticlesSynchronizer()
   {
      isSet = false;
   }

   //done
   DistributedParticlesSynchronizer( const DistributedGridType* distributedGrid )
   {
      isSet = false;
      setDistributedGrid( distributedGrid );
   }

   void
   initialize( const DistributedParticlesPointerType& distributedParticles )
   {
      //this->distributedParticles = std::make_shared< DistributedParticlesPointerType >( distributedParticles );

      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();
      const auto innerOverlapsView = distributedParticles->getInnerOverlaps().getConstView();
      int bufferSize = 0;
      for( int i = 0; i < this->getNeighborsCount(); i++ )
         if( neighbors[ i ] != -1 )
             bufferSize += innerOverlapsView[ i ].getNumberOfCells() * innerOverlapsView[ i ].getNumberOfParticlesPerCell();

      this->sendBuffers.setSize( bufferSize * Particles::getParticlesDimension() );
   }


   //TODO: Temp, test
   void
   setCommunicator( MPI::Comm& comm )
   {
      this->communicator = comm;
   }

   template< typename DistributedParticlesPointer, typename ParticlesPointerType  >
   void
   synchronizeOverlapSizes( const DistributedParticlesPointer& distributedParticles,
                            const ParticlesPointerType& particles )
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedParticles->getDistributedGrid().isDistributed() )
         return;

      // resent send and recv sizes
      sendSizes = 0;
      recvSizes = 0;
      sendNeighborOffsets = 0;
      recvNeighborOffsets = 0;
      //if( TNL::MPI::GetRank() == 2 )
      //std::cout << "RANK: " << communicator.rank() << " EMPTY: [ sendSize ] = " << sendSizes << " [ sendNeighborOffsets ] = " << sendNeighborOffsets << " [ recvSizes ] = " << recvSizes << " [ recvNeighborOffsets ] = " << recvNeighborOffsets << std::endl;

      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();

      // fill send sizes
      const auto innerOverlaps_view = distributedParticles->getInnerOverlaps().getConstView();
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 )
            sendSizes[ i ] = innerOverlaps_view[ i ].getNumberOfParticles();
      }

      // async send and receive
      std::unique_ptr< MPI_Request[] > requests{ new MPI_Request[ 2 * this->getNeighborsCount() ] };
      //const MPI::Comm& communicator = distributedParticles->getCommunicator();
      int requestsCount = 0;

      // send everything, recieve everything
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            requests[ requestsCount++ ] = MPI::Isend( &sendSizes[ i ], 1, neighbors[ i ], 0, communicator );
            requests[ requestsCount++ ] = MPI::Irecv( &recvSizes[ i ], 1, neighbors[ i ], 0, communicator );
         }
         //TODO: There should be part to include periodicity
         //else if( sendSizes[ i ] != 0 ) {
         //   requests[ requestsCount++ ] = MPI::Isend( &sendSizes[ i ], 1, neighbors[ i ], 1, communicator );
         //   requests[ requestsCount++ ] = MPI::Irecv( &recvSizes[ i ], 1, neighbors[ i ], 1, communicator );
         //}
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );

      // set offsets
      //FIXME: It is not possible to perfrom scan on static vector. I should define sequential scan manually
      for( int i = 0; i < this->getNeighborsCount(); i++ ){
         sendNeighborOffsets[ i ] = 0;
         //NOTE: This should be splited into two synchronizer approaches
         //since we decided to put overlaps to the same sets, we should include index of last active ptcs
         recvNeighborOffsets[ i ] = particles->getNumberOfParticles();
         for( int j = 0; j < i; j++ ) {
            sendNeighborOffsets[ i ] += sendSizes[ j ];
            recvNeighborOffsets[ i ] += recvSizes[ j ];
         }
      }
      //if( TNL::MPI::GetRank() == 2 )
      //std::cout << "RANKEND " << communicator.rank() << " FILLED: [ sendSize ] = " << sendSizes << " [ sendNeighborOffsets ] = " << sendNeighborOffsets << " [ recvSizes ] = " << recvSizes << " [ recvNeighborOffsets ] = " << recvNeighborOffsets << std::endl;
   }

   template< typename DistributedParticlesPointer >
   void
   synchronizeBalancingMeasures( DistributedParticlesPointer& distributedParticles )
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedParticles->getDistributedGrid().isDistributed() )
         return;

      // resent send and recv sizes
      sendParticlesCount = 0;
      recvParticlesCount = 0;
      sendCompTime = 0;
      recvCompTime = 0;
      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();

      // fill send sizes
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ){
            sendParticlesCount[ i ] = distributedParticles->getNumberOfParticlesForLoadBalancing();
            sendCompTime[ i ] = distributedParticles->getCompTime();
         }
      }

      // async send and receive
      std::unique_ptr< MPI_Request[] > requests{ new MPI_Request[ 2 * this->getNeighborsCount() ] };
      //const MPI::Comm& communicator = distributedParticles->getCommunicator();
      int requestsCount = 0;

      // send everything, recieve everything
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            requests[ requestsCount++ ] = MPI::Isend( &sendParticlesCount[ i ], 1, neighbors[ i ], 0, communicator );
            requests[ requestsCount++ ] = MPI::Isend( &sendCompTime[ i ], 1, neighbors[ i ], 0, communicator );
            requests[ requestsCount++ ] = MPI::Irecv( &recvParticlesCount[ i ], 1, neighbors[ i ], 0, communicator );
            requests[ requestsCount++ ] = MPI::Irecv( &recvCompTime[ i ], 1, neighbors[ i ], 0, communicator );
         }
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );

      // store recv values
      //auto subdomainsParticlesCount_view = distributedParticles->getSubdomainsParticlesCount().getView();
      //auto subdomainsCompTime_view = distributedParticles->getSubdomainsCompTime().getView();
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ){
            //subdomainsParticlesCount_view[ i ] = recvParticlesCount[ i ];
            //subdomainsCompTime_view[ i ] = recvCompTime[ i ];
            distributedParticles->getSubdomainsParticlesCount()[ i ] = recvParticlesCount[ i ];
            distributedParticles->getSubdomainsCompTime()[ i ] = recvCompTime[ i ];
         }
      }
   }

   GlobalIndexType
   getNumberOfRecvParticles()
   {
      const GlobalIndexType ret = recvNeighborOffsets[ getNeighborsCount() - 1 ] - recvNeighborOffsets[ 0 ];
      return ret;
   }


   void
   updateOffsets( DistributedParticlesPointerType& distributedParticles )
   {
      //TODO: collect particles in overlaps HERE
      //distributedParticles->collectParticlesInInnerOverlaps( particles )

      //TODO: fill
      //ghostNeighbosOffsets
      //ghostNeihbors

      //TODO: Neighbors needs to be initialized

   }

   template< typename Array, typename DistributedParticlesPointer >
   void
   synchronize( Array& array, DistributedParticlesPointer& distributedParticles )
   {
      // wrapped only because nvcc is fucked up and does not like __cuda_callable__ lambdas in enable_if methods
      synchronizeArray( array, distributedParticles );
   }

   template< typename Array, typename DistributedParticlesPointer >
   void
   synchronizeArray( Array& array,
                     DistributedParticlesPointer& distributedParticles,
                     int valuesPerElement = 1 )
   {
      static_assert( std::is_same_v< typename Array::DeviceType, DeviceType >, "mismatched DeviceType of the array" );
      using ValueType = typename Array::ValueType;

      ByteArrayView arrayView;
      arrayView.bind( reinterpret_cast< std::uint8_t* >( array.getData() ), sizeof( ValueType ) * array.getSize() );
      synchronizeByteArray( arrayView, sizeof( ValueType ) * valuesPerElement, distributedParticles );
   }

   template< typename DistributedParticlesPointer >
   void
   synchronizeByteArray( ByteArrayView sendArray,
                         int bytesPerValue,
                         DistributedParticlesPointer& distributedParticles )
   {
      auto requests = synchronizeByteArrayAsyncWorker( sendArray, bytesPerValue, distributedParticles );
      MPI::Waitall( requests.data(), requests.size() );
   }

   template< typename DistributedParticlesPointer >
   [[nodiscard]] RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView sendArray,
                                    int bytesPerValue,
                                    DistributedParticlesPointer& distributedParticles )
   {
      if( sendArray.getSize() < bytesPerValue * this->sendBufferElementSize )
         throw std::invalid_argument( "synchronizeByteArrayAsyncWorker: the array does not have the expected size" );

      const int rank = communicator.rank();
      const int nproc = communicator.size();

      //TODO: Do this in preivous step, set max size ptcs per zone cells.
      // allocate send buffers (setSize does nothing if the array size is already correct)
      // sendBuffers.setSize( bytesPerValue * ghostNeighborOffsets[ nproc ] );

      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();

      // buffer for asynchronous communication requests
      RequestsVector requests;

      // issue all receive async operations
      for( int j = 0; j < this->getNeighborsCount(); j++ ) {
         if( neighbors[ j ] != -1 ) {
            //NOTE: This should be splited into two synchronizer approaches.
            //      For now, we use single field with mark, but in case of overlap, recvArray should be placed here.
            requests.push_back( MPI::Irecv( sendArray.getData() + bytesPerValue * recvNeighborOffsets[ j ],
                                            bytesPerValue * recvSizes[ j ],
                                            neighbors[ j ],
                                            0,
                                            communicator ) );
         }
      }

      //TODO: Here, I modified the size of the binded ByteArrayView
      ByteArrayView sendBuffersView;
      sendBuffersView.bind( sendBuffers.getData(), bytesPerValue * sendBuffers.getSize() );

      const auto innerOverlapsView = distributedParticles->getInnerOverlaps().getConstView();
      const auto sendArrayView = sendArray.getConstView();
      //TODO: The kernel could be probabaly defined here as in Distributed Mesh Synchronizer

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            const GlobalIndexType offset = sendNeighborOffsets[ i ];

            const auto ghostZoneView = innerOverlapsView[ i ].getParticlesInZone().getConstView();
            auto copy_kernel = [ sendBuffersView, sendArrayView, ghostZoneView, bytesPerValue ] __cuda_callable__(
                                  GlobalIndexType k, GlobalIndexType offset ) mutable
            {
               for( int i = 0; i < bytesPerValue; i++ )
                  sendBuffersView[ i + bytesPerValue * ( offset + k ) ] =
                     sendArrayView[ i + bytesPerValue * ghostZoneView[ k ] ];
            };

            // copy data to send buffers
            Algorithms::parallelFor< DeviceType >( 0, sendSizes[ i ], copy_kernel, offset );

            // issue async send operation
            requests.push_back( MPI::Isend( sendBuffersView.getData() + bytesPerValue * sendNeighborOffsets[ i ],
                                            bytesPerValue * sendSizes[ i ],
                                            neighbors[ i ],
                                            0,
                                            communicator ) );
         }
      }

      return requests;
   }


private:


   std::vector< ParticleZone > innerGhostZones;


   //Containers::Array< std::uint8_t, Device, Index > sendBuffers[ getNeighborsCount() ];
   Containers::Array< std::uint8_t, Device, Index > recieveBuffers[ getNeighborsCount() ];

   //PeriodicBoundariesCopyDirection periodicBoundariesCopyDirection = BoundaryToOverlap; //TODO: This will be necessary

   CoordinatesType sendDimensions[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType recieveDimensions[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType sendBegin[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType recieveBegin[ getNeighborsCount() ]; //TODO: not necessary

   //const DistributedSimulationType* distributedSimulation;
   bool isSet;

   MPI::Comm communicator;

   /**
    * Send buffers: array for buffering the mesh function values which will be
    * sent to other ranks. We use std::uint8_t as the value type to make this
    * class independent of the mesh function's real type. When cast to the real
    * type values, the send buffer for the i-th rank is the part of the array
    * starting at index ghostNeighborOffsets[i] (inclusive) and ending at index
    * ghostNeighborOffsets[i+1] (exclusive).
    */
   GlobalIndexType sendBufferElementSize = 0;
   Containers::Array< std::uint8_t, DeviceType, GlobalIndexType > sendBuffers;

   Containers::StaticArray< getNeighborsCount(), int > sendSizes;
   Containers::StaticArray< getNeighborsCount(), int > recvSizes;

   // load balancing props
   Containers::StaticArray< getNeighborsCount(), int > sendParticlesCount;
   Containers::StaticArray< getNeighborsCount(), int > recvParticlesCount;

   Containers::StaticArray< getNeighborsCount(), float > sendCompTime;
   Containers::StaticArray< getNeighborsCount(), float > recvCompTime;

   Containers::StaticArray< getNeighborsCount(), int > sendNeighborOffsets;
   Containers::StaticArray< getNeighborsCount(), int > recvNeighborOffsets;

   //Containers::Array< int, TNL::Devices::Host > sendNeighborOffsets;
   //Containers::Array< int, TNL::Devices::Host > recvNeighborOffsets;

   //const DistributedParticlesPointerType distributedParticles;


};

} // ParticleSystem
} // TNL

