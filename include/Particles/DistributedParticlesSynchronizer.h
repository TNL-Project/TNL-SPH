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

   //TODO: Temp, test
   void
   setCommunicator( MPI::Comm& comm )
   {
      this->communicator = comm;
   }

   void
   collectParticlesToSend()
   {
      // --> update cells in zone


      //1. get local grid -> from the local grid, extract neighbors
      //2. separate cells which are crutical for given neighbor

      // --> collect particles in zone cells

   }

   //template< typename GhostBoundaryPatches >
   //void
   //setDistributedGrid( const DistributedGridType* distributedGrid,
   //                    GhostBoundaryPatches& ghostBoundaryPatches )
   //{
   //   isSet = true;

   //   this->distributedGrid = distributedGrid;

   //   const SubdomainOverlapsType& lowerOverlap = this->distributedGrid->getLowerOverlap();
   //   const SubdomainOverlapsType& upperOverlap = this->distributedGrid->getUpperOverlap();

   //   const CoordinatesType& localBegin = this->distributedGrid->getLocalMesh().getLocalBegin();
   //   const CoordinatesType& localSize = this->distributedGrid->getLocalSize();

   //   const int* neighbors = distributedGrid->getNeighbors();

   //   for( int i = 0; i < this->getNeighborsCount(); i++ ) {
   //      Index sendSize = 1;  // send and receive  areas have the same size

   //      auto directions = Directions::template getXYZ< getMeshDimension() >( i );

   //      sendDimensions[ i ] = localSize;  // send and receive areas have the same dimensions
   //      sendBegin[ i ] = localBegin;
   //      recieveBegin[ i ] = localBegin;

   //      for( int j = 0; j < this->getMeshDimension(); j++ ) {
   //         if( directions[ j ] == -1 ) {
   //            sendDimensions[ i ][ j ] = lowerOverlap[ j ];
   //            recieveBegin[ i ][ j ] = 0;
   //         }

   //         if( directions[ j ] == 1 ) {
   //            sendDimensions[ i ][ j ] = upperOverlap[ j ];
   //            sendBegin[ i ][ j ] = localBegin[ j ] + localSize[ j ] - upperOverlap[ j ];
   //            recieveBegin[ i ][ j ] = localBegin[ j ] + localSize[ j ];
   //         }

   //         sendSize *= sendDimensions[ i ][ j ];
   //      }

   //      sendSizes[ i ] = sendSize;
   //   }

   //   for( int i = 0; i < this->getNeighborsCount(); i++ ) {
   //      innerGhostZones[ i ].updateCells( sendBegin[ i ], sendDimensions[ i ] );
   //   }
   //}

   template< typename DistributedParticlesPointer >
   void
   synchronizeOverlapSizes( const DistributedParticlesPointer& distributedParticles )
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedParticles->getDistributedGrid().isDistributed() )
         return;

      // resent send and recv sizes
      sendSizes = 0;
      recvSizes = 0;
      sendNeighborOffsets = 0;
      recvNeighborOffsets = 0;
      std::cout << "EMPTY: [ sendSize ] = " << sendSizes << " [ sendNeighborOffsets ] = " << sendNeighborOffsets << " [ recvSizes ] = " << recvSizes << " [ recvNeighborOffsets ] = " << recvNeighborOffsets << std::endl;

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
         recvNeighborOffsets[ i ] = 0;
         for( int j = 0; j < i; j++ ) {
            sendNeighborOffsets[ i ] += sendSizes[ j ];
            recvNeighborOffsets[ i ] += recvSizes[ j ];
         }
      }
      std::cout << "FILLED: [ sendSize ] = " << sendSizes << " [ sendNeighborOffsets ] = " << sendNeighborOffsets << " [ recvSizes ] = " << recvSizes << " [ recvNeighborOffsets ] = " << recvNeighborOffsets << std::endl;
   }


   template< typename DistributedParticlesPointer, typename ArrayType >
   void
   synchronize( ArrayType& sendFunction,
                ArrayType& receiveFunction,
                DistributedParticlesPointer& distributedParticles)
   {
      using ValueType = typename ArrayType::ValueType;

      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedParticles->getDistributedGrid().isDistributed() )
         return;

      // allocate buffers (setSize does nothing if the array size is already correct)
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         sendBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
         recieveBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
      }

      //TODO: Synchronizer should probabaly store pointer to distributed parricles
      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();
      const auto innerOverlaps_view = distributedParticles->getInnerOverlaps().getConstView();

      // fill send buffers
      auto array_view = sendFunction.getView();

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         //const auto zoneParticleIndices_view = ghostBoundaryPatches[ i ]->zone.getParticlesInZone().getConstView();
         //const GlobalIndexType numberOfZoneParticles = ghostBoundaryPatches[ i ]->zone.getNumberOfParticles();
         const auto zoneParticleIndices_view = innerOverlaps_view[ i ].getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = innerOverlaps_view[ i ].getNumberOfParticles();
         auto sendBufferView = this->sendBuffers[ i ].getView();

         if( numberOfZoneParticles == 0 )
            continue;

         //FIXME: Use receive array for linearization of sending field
         //auto recieveArray_view = receiveFunction.getView();
         //auto copyVariables = [=] __cuda_callable__ ( int j ) mutable
         //{
         //   const Index p = zoneParticleIndices_view[ j ];
         //   recieveArray_view[ j ] = array_view[ p ];
         //};
         //Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, copyVariables );
         //sendBufferView.bind( reinterpret_cast< std::uint8_t* >( recieveArray_view.getData() ), sizeof( ValueType ) * numberOfZoneParticles );

         //auto copyVariables = [=] __cuda_callable__ ( int j ) mutable
         //{
         //   const Index p = zoneParticleIndices_view[ j ];
         //   recieveArray_view[ j ] = array_view[ p ];
         //};
         //Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, copyVariables );

         //auto moveParticles = [ = ] __cuda_callable__( int k ) mutable
         //{
         //   const Index p = zoneParticleIndices_view[ k ];
         //   //sendBufferView[ k ] = static_cast< uint8_t >( array_view[ p ].getData() );
         //   //sendBufferView[ k ] = reinterpret_cast< uint8_t >( array_view[ p ].getData() );
         //};
         //Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, moveParticles );


      }



      /*
      // async send and receive
      std::unique_ptr< MPI_Request[] > requests{ new MPI_Request[ 4 * this->getNeighborsCount() ] };
      const MPI::Comm& communicator = distributedParticles->getCommunicator();
      int requestsCount = 0;

      // send everything, recieve everything
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            requests[ requestsCount++ ] = MPI::Isend( reinterpret_cast< RealType* >( sendBuffers[ i ].getData() ),
                                                      sendSizes[ i ],
                                                      neighbors[ i ],
                                                      0,
                                                      communicator );
            requests[ requestsCount++ ] = MPI::Irecv( reinterpret_cast< RealType* >( recieveBuffers[ i ].getData() ),
                                                      sendSizes[ i ],
                                                      neighbors[ i ],
                                                      0,
                                                      communicator );
         }
         //TODO: Periodic boundaries should be resolved
         //else if( sendSizes[ i ] != 0 ) {
         //   requests[ requestsCount++ ] = MPI::Isend( reinterpret_cast< RealType* >( sendBuffers[ i ].getData() ),
         //                                             sendSizes[ i ],
         //                                             periodicNeighbors[ i ],
         //                                             1,
         //                                             communicator );
         //   requests[ requestsCount++ ] = MPI::Irecv( reinterpret_cast< RealType* >( recieveBuffers[ i ].getData() ),
         //                                             sendSizes[ i ],
         //                                             periodicNeighbors[ i ],
         //                                             1,
         //                                             communicator );
         //}
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );
      */

      /*
      // copy data from receive buffers
      auto receiveArray_view = receiveFunction.getView();
      int offset = 0; //FIXME: RESOLVE OFFSETS

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         const auto recieveBuffers_view = recieveBuffers[ i ].getConstView();

         auto copy = [=]  __cuda_callable__ ( Index j )
         {
             receiveArray_view[ j + offset ] = recieveBuffers_view[ j ];
         };
         Algorithms::parallelFor< Device >( 0, recvSizes[ i ], copy );
      }
      */

   }

   template< typename Array >
   void
   synchronize( Array& sendArray, Array& recvArray )
   {
      // wrapped only because nvcc is fucked up and does not like __cuda_callable__ lambdas in enable_if methods
      synchronizeArray( sendArray, recvArray );
   }

   template< typename Array >
   void
   synchronizeArray( Array& sendArray, Array& recvArray, int valuesPerElement = 1 )
   {
      static_assert( std::is_same_v< typename Array::DeviceType, DeviceType >, "mismatched DeviceType of the array" );
      using ValueType = typename Array::ValueType;

      ByteArrayView sendView;
      ByteArrayView recvView;
      sendView.bind( reinterpret_cast< std::uint8_t* >( sendArray.getData() ), sizeof( ValueType ) * sendArray.getSize() );
      recvView.bind( reinterpret_cast< std::uint8_t* >( recvArray.getData() ), sizeof( ValueType ) * recvArray.getSize() );
      synchronizeByteArray( sendView, recvView, sizeof( ValueType ) * valuesPerElement );
   }

   void
   synchronizeByteArray( ByteArrayView sendArray, ByteArrayView recvArray, int bytesPerValue )
   {
      auto requests = synchronizeByteArrayAsyncWorker( sendArray, recvArray, bytesPerValue );
      MPI::Waitall( requests.data(), requests.size() );
   }

   [[nodiscard]] RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView sendArray, ByteArrayView recvArray, int bytesPerValue )
   {
      //TODO: Bring back this test.
      //if( sendArray.getSize() != bytesPerValue * ghostOffsets[ ghostOffsets.getSize() - 1 ] )
      //   throw std::invalid_argument( "synchronizeByteArrayAsyncWorker: the array does not have the expected size" );

      const int rank = communicator.rank();
      const int nproc = communicator.size();

      //TODO: Do this in preivous step, set max size ptcs per zone cells.
      // allocate send buffers (setSize does nothing if the array size is already correct)
      // sendBuffers.setSize( bytesPerValue * ghostNeighborOffsets[ nproc ] );

      const int* neighbors = this->distributedParticles->getDistributedGrid().getNeighbors();

      // buffer for asynchronous communication requests
      RequestsVector requests;

      // issue all receive async operations
      for( int j = 0; j < this->getNeighborsCount(); j++ ) {
         if( neighbors[ j ] != -1 ) {
            requests.push_back( MPI::Irecv( recvArray.getData() + bytesPerValue * recvNeighborOffsets[ j ],
                                            bytesPerValue * recvSizes[ j ],
                                            j,
                                            0,
                                            communicator ) );
         }
      }

      //TODO: Here, I modified the size of the binded ByteArrayView
      ByteArrayView sendBuffersView;
      //sendBuffersView.bind( sendBuffers.getData(), bytesPerValue * ghostNeighborOffsets[ nproc ] );
      sendBuffersView.bind( sendBuffers.getData(), bytesPerValue * sendBuffers.getSize() );

      const auto sendArrayView = sendArray.getConstView();
      //TODO: The kernel could be probabaly defined here as in Distributed Mesh Synchronizer


      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            const GlobalIndexType offset = sendNeighborOffsets[ i ];

            const auto innerOverlapsView = distributedParticles->getInnerOverlaps()[ i ].getParticlesInZone().getConstView();
            auto copy_kernel = [ sendBuffersView, sendArrayView, innerOverlapsView, bytesPerValue ] __cuda_callable__(
                                  GlobalIndexType k, GlobalIndexType offset ) mutable
            {
               for( int i = 0; i < bytesPerValue; i++ )
                  sendBuffersView[ i + bytesPerValue * ( offset + k ) ] =
                     sendArrayView[ i + bytesPerValue * innerOverlapsView[ k ] ];
            };

            // copy data to send buffers
            Algorithms::parallelFor< DeviceType >( 0, sendSizes[ i ], copy_kernel, offset );

            // issue async send operation
            requests.push_back( MPI::Isend( sendBuffersView.getData() + bytesPerValue * sendNeighborOffsets[ i ],
                                            bytesPerValue * sendSizes[ i ],
                                            i,
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
   Containers::Array< std::uint8_t, DeviceType, GlobalIndexType > sendBuffers;

   Containers::Array< std::uint8_t, DeviceType, GlobalIndexType > recvBuffers;

   Containers::StaticArray< getNeighborsCount(), int > sendSizes;
   Containers::StaticArray< getNeighborsCount(), int > recvSizes;

   Containers::StaticArray< getNeighborsCount(), int > sendNeighborOffsets;
   Containers::StaticArray< getNeighborsCount(), int > recvNeighborOffsets;

   //Containers::Array< int, TNL::Devices::Host > sendNeighborOffsets;
   //Containers::Array< int, TNL::Devices::Host > recvNeighborOffsets;

   const DistributedParticlesPointerType distributedParticles;


};

} // ParticleSystem
} // TNL

