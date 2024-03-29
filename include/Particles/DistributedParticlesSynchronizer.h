#pragma once

#include "ParticlesTraits.h"
#include "GhostZone.h"
#include <vector>

#include <TNL/Containers/ByteArraySynchronizer.h>
#include <TNL/MPI/Comm.h>
#include <TNL/MPI/Wrappers.h>

namespace TNL {
namespace ParticleSystem {

template< typename Particles, typename DistributedGridType >
class DistributedParticlesSynchronizer
//: public Containers::ByteArraySynchronizer< typename Particles::DeviceType, typename Particles::GlobalIndexType >
{
public:


   using RealType = typename Particles::RealType;
   using GlobalIndexType = typename Particles::RealType;
   //using Index = typename DistributedGridType::Index;
   using Index = typename Particles::GlobalIndexType; //FIXME: Back to grid
   using ParticleZone = ParticleZone< typename Particles::Config >;
   using Device = typename Particles::DeviceType;

   using CoordinatesType = typename DistributedGridType::CoordinatesType;
   using SubdomainOverlapsType = typename DistributedGridType::SubdomainOverlapsType;

   //using Base = Containers::ByteArraySynchronizer< typename Particles::DeviceType, typename Particles::GlobalIndexType >;
   //using RequestsVector = typename Base::RequestsVector;
   using ByteArrayView = Containers::ArrayView< std::uint8_t, Device, GlobalIndexType >;
   using RequestsVector = std::vector< MPI_Request >;

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

      const int* neighbors = distributedParticles->getDistributedGrid().getNeighbors();

      // fill send sizes
      const auto innerOverlaps_view = distributedParticles->getInnerOverlaps().getConstView();
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
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
            requests[ requestsCount++ ] = MPI::Irecv( &receivedSizes[ i ], 1, neighbors[ i ], 0, communicator );
         }
         //TODO: What is purpose of this part? Periodicity?
         //else if( sendSizes[ i ] != 0 ) {
         //   requests[ requestsCount++ ] = MPI::Isend( &sendSizes[ i ], 1, neighbors[ i ], 1, communicator );
         //   requests[ requestsCount++ ] = MPI::Irecv( &receivedSizes[ i ], 1, neighbors[ i ], 1, communicator );
         //}
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );
   }

   template< typename DistributedParticlesPointer, typename ArrayType, typename GhostBoundaryPatches >
   void
   synchronize( ArrayType& sendFunction,
                ArrayType& receiveFunction,
                GhostBoundaryPatches& ghostBoundaryPatches,
                DistributedParticlesPointer& distributedParticles)
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedParticles->getDistributedGrid()->isDistributed() )
         return;

      // allocate buffers (setSize does nothing if the array size is already correct)
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         sendBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
         recieveBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
      }

      //TODO: Synchronizer should probabaly store pointer to distributed parricles
      const int* neighbors = distributedParticles()->getDistributedGrid()->getNeighbors();

      // fill send buffers
      const auto array_view = sendFunction.getView();

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         const auto zoneParticleIndices_view = ghostBoundaryPatches[ i ]->zone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = ghostBoundaryPatches[ i ]->zone.getNumberOfParticles();
         auto sendBufferView = sendBuffers[ i ].getView();

         auto copy = [=] __cuda_callable__ ( Index j )
         {
            const Index p = zoneParticleIndices_view[ i ];
            sendBufferView[ j ] = array_view[ p ];
         };
         Algorithms::parallelFor< Device >( 0, numberOfZoneParticles, copy );
      }

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

      // copy data from receive buffers
      auto receiveArray_view = receiveFunction.getView();
      int offset = 0; //FIXME: RESOLVE OFFSETS

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         const auto recieveBuffers_view = recieveBuffers[ i ].getConstView();

         auto copy = [=]  __cuda_callable__ ( Index j )
         {
             receiveArray_view[ j + offset ] = recieveBuffers_view[ j ];
         };
         Algorithms::parallelFor< Device >( 0, receivedSizes[ i ], copy );
      }

   }

private:

   MPI::Comm communicator;

   std::vector< ParticleZone > innerGhostZones;

   Containers::StaticArray< getNeighborsCount(), int > sendSizes;
   Containers::StaticArray< getNeighborsCount(), int > receivedSizes;

   Containers::Array< std::uint8_t, Device, Index > sendBuffers[ getNeighborsCount() ];
   Containers::Array< std::uint8_t, Device, Index > recieveBuffers[ getNeighborsCount() ];

   //PeriodicBoundariesCopyDirection periodicBoundariesCopyDirection = BoundaryToOverlap; //TODO: This will be necessary

   CoordinatesType sendDimensions[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType recieveDimensions[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType sendBegin[ getNeighborsCount() ]; //TODO: not necessary
   CoordinatesType recieveBegin[ getNeighborsCount() ]; //TODO: not necessary

   //const DistributedSimulationType* distributedSimulation;
   bool isSet;


};

} // ParticleSystem
} // TNL

