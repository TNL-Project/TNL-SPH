#pragma once

#include "../Particles/ParticlesTraits.h"
#include "SPHTraits.h"
#include "../Particles/GhostZone.h"
#include <vector>

namespace TNL {
namespace SPH {

template< typename Particles, typename DistributedGridType >
class DistributedSimulationSynchronizer
{
public:


   using RealType = typename Particles::RealType;
   using GlobalIndexType = typename Particles::RealType;
   using Index = typename DistributedGridType::Index;
   using ParticleZone = ParticleZone< typename Particles::Config >;
   using Device = typename Particles::DeviceType;

   using CoordinatesType = typename DistributedGridType::CoordinatesType;
   using SubdomainOverlapsType = typename DistributedGridType::SubdomainOverlapsType;

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
   DistributedSimulationSynchronizer()
   {
      isSet = false;
   }

   //done
   DistributedSimulationSynchronizer( const DistributedGridType* distributedGrid )
   {
      isSet = false;
      setDistributedGrid( distributedGrid );
   }

   void
   synchronizeOverlapSizes()
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedSimulation->distributedGrid->isDistributed() )
         return;

      const int* neighbors = distributedSimulation->distributedGrid->getNeighbors();

      // async send and receive
      std::unique_ptr< MPI_Request[] > requests{ new MPI_Request[ 2 * this->getNeighborsCount() ] };
      const MPI::Comm& communicator = distributedGrid->getCommunicator();
      int requestsCount( 0 );

      // send everything, recieve everything
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         if( neighbors[ i ] != -1 ) {
            requests[ requestsCount++ ] = MPI::Isend( sendSizes[ i ], 1, neighbors[ i ], 0, communicator );
            requests[ requestsCount++ ] = MPI::Rrecv( recievedSizes[ i ], 1, neighbors[ i ], 0, communicator );
         }
         else if( sendSizes[ i ] != 0 ) {
            requests[ requestsCount++ ] = MPI::Isend( sendSizes[ i ], 1, neighbors[ i ], 1, communicator );
            requests[ requestsCount++ ] = MPI::Rrecv( recievedSizes[ i ], 1, neighbors[ i ], 1, communicator );
         }
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );
   }

   template< typename ArrayType, typename GhostBoundaryPatches >
   void
   synchronize( ArrayType& sendFunction,
                ArrayType& receiveFunction,
                GhostBoundaryPatches& ghostBoundaryPatches )
   {
      TNL_ASSERT_TRUE( isSet, "Synchronizer is not set, but used to synchronize" );
      if( ! distributedSimulation->distributedGrid->isDistributed() )
         return;

      // allocate buffers (setSize does nothing if the array size is already correct)
      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         sendBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
         recieveBuffers[ i ].setSize( sendSizes[ i ] * sizeof( RealType ) );
      }

      const int* neighbors = distributedSimulation->distributedGrid->getNeighbors();

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
      const MPI::Comm& communicator = distributedGrid->getCommunicator();
      int requestsCount( 0 );

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
         else if( sendSizes[ i ] != 0 ) {
            requests[ requestsCount++ ] = MPI::Isend( reinterpret_cast< RealType* >( sendBuffers[ i ].getData() ),
                                                      sendSizes[ i ],
                                                      periodicNeighbors[ i ],
                                                      1,
                                                      communicator );
            requests[ requestsCount++ ] = MPI::Irecv( reinterpret_cast< RealType* >( recieveBuffers[ i ].getData() ),
                                                      sendSizes[ i ],
                                                      periodicNeighbors[ i ],
                                                      1,
                                                      communicator );
         }
      }

      // wait until send is done
      MPI::Waitall( requests.get(), requestsCount );

      // copy data from receive buffers
      auto receiveArray_view = receiveFunction.getView();

      for( int i = 0; i < this->getNeighborsCount(); i++ ) {
         const auto recieveBuffers_view = sendBuffers[ i ].getConstView();

         auto copy = [=]  __cuda_callable__ ( Index j )
         {
            const Index p = zoneParticleIndices_view[ i ];
             sendBufferView[ j ] = array_view[ p ];
         };
         Algorithms::parallelFor< Device >( 0, recievedSizes[ i ], copy );
      }

   }

private:
   std::vector< ParticleZone > innerGhostZones;

   Containers::StaticArray< getNeighborsCount(), int > sendSizes;
   Containers::StaticArray< getNeighborsCount(), int > recievedSizes;

   Containers::Array< std::uint8_t, Device, Index > sendBuffers[ getNeighborsCount() ];
   Containers::Array< std::uint8_t, Device, Index > recieveBuffers[ getNeighborsCount() ];

   PeriodicBoundariesCopyDirection periodicBoundariesCopyDirection = BoundaryToOverlap;

   CoordinatesType sendDimensions[ getNeighborsCount() ]; //TODO: Not necessary
   CoordinatesType recieveDimensions[ getNeighborsCount() ]; //TODO: Not necessary
   CoordinatesType sendBegin[ getNeighborsCount() ]; //TODO: Not necessary
   CoordinatesType recieveBegin[ getNeighborsCount() ]; //TODO: Not necessary

   const DistributedSimulationType* distributedSimulation;
   bool isSet;


};

} // SPH
} // TNL

