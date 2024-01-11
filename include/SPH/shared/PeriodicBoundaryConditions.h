#include <TNL/Algorithms/parallelFor.h>

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename Model >
class PeriodicBoundaryConditionsShared
{
public:

   using SPHConfig = typename Model::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;
   using ParticlesType = typename Model::ParticlesType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using PairIndexType = typename ParticlesType::PairIndexType;

   template< typename Array >
   void
   copyPartOfArray( Array& array, GlobalIndexType fromPosition, GlobalIndexType toPosition, GlobalIndexType size )
   {
      auto view = array.getView();
      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename Array >
   void
   copyPartOfPoints( Array& array,
                     GlobalIndexType fromPosition,
                     GlobalIndexType toPosition,
                     GlobalIndexType size,
                     VectorType coordinatesDifference )
   {
      auto view = array.getView();
      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ] + coordinatesDifference;
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename Array >
   void
   arrangeArray( Array& arrayCopy,
                 Array& arrayPaste,
                 GlobalIndexType fromPosition,
                 GlobalIndexType toPosition,
                 GlobalIndexType size )
   {
      const auto viewArrayCopy = arrayCopy.getConstView();
      auto viewPaste = arrayPaste.getView();
      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         viewPaste[ toPosition + i ] = viewArrayCopy[ fromPosition + i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
      arrayCopy.swap( arrayPaste );
   }

   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   void
   initializePeriodicBoundaryTransfer( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const VectorType coordinatesDifference = particlesParams.periodicBoundaryDistance;

      firstAndLastParticleInFirstBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
      particlesParams.indexOfColumnWithLeftPeriodicity );
      firstAndLastParticleInLastBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
      particlesParams.indexOfColumnWithRightPeriodicity );

      firstActiveParticle = physicalObject->getFirstActiveParticle();
      lastActiveParticle = physicalObject->getLastActiveParticle();

      sizeToCopyFirstBlock = firstAndLastParticleInFirstBlock[ 1 ] - firstAndLastParticleInFirstBlock[ 0 ] + 1;
      sizeToCopyLastBlock = firstAndLastParticleInLastBlock[ 1 ] - firstAndLastParticleInLastBlock[ 0 ] + 1;

      if( firstAndLastParticleInFirstBlock[ 0 ] != firstActiveParticle ){
         std::cerr << "Periodic boundary error: first particle in first cell doesn't match the first particle of particle system. FirstAndLastParticleInFirstBlock[ 0 ]:"  <<
                       firstAndLastParticleInFirstBlock[ 0 ] << " and firstActiveParticle:" << firstActiveParticle << std::endl;
      }

      if(firstAndLastParticleInLastBlock[ 1 ] != lastActiveParticle ){
         std::cerr << "Periodic boundary error: last particle in last cell doesn't match the first particle of particle system. FirstAndLastParticleInLastBlock[ 1 ]:"  <<
                       firstAndLastParticleInLastBlock[ 1 ] << " and firstActiveParticle:" << lastActiveParticle << std::endl;
      }

      //Copy data from periodic patch A to periodic patch B:
      copyPartOfPoints( physicalObject->getPoints(),
                        firstAndLastParticleInFirstBlock[ 0 ],
                        lastActiveParticle + 1,
                        sizeToCopyFirstBlock,
                        coordinatesDifference );

      //Copy data from periodic patch B to periodic patch A:
      copyPartOfPoints( physicalObject->getPoints(),
                        firstAndLastParticleInLastBlock[ 0 ],
                        firstActiveParticle - sizeToCopyLastBlock,
                        sizeToCopyLastBlock,
                        ( -1.f ) * coordinatesDifference );

      physicalObject->particles->setFirstActiveParticle( firstActiveParticle - sizeToCopyLastBlock );
      physicalObject->particles->setLastActiveParticle( lastActiveParticle + sizeToCopyFirstBlock );
      const GlobalIndexType newNumberOfParticles = ( lastActiveParticle + sizeToCopyFirstBlock ) - ( firstActiveParticle - sizeToCopyLastBlock ) + 1;
      physicalObject->particles->setNumberOfParticles( newNumberOfParticles );

      this->synchronizingPeriodic = true;

   }

   template< typename Array >
   void
   periodicUpdateOfParticleField( Array& array )
   {
      if( synchronizingPeriodic == false ) return;

      //Copy data from periodic patch A to periodic patch B:
      copyPartOfArray( array,
                       firstAndLastParticleInFirstBlock[ 0 ],
                       lastActiveParticle + 1,
                       sizeToCopyFirstBlock );

      //Copy data from periodic patch B to periodic patch A:
      copyPartOfArray( array,
                       firstAndLastParticleInLastBlock[ 0 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyLastBlock );
   }

   void
   finalizePeriodicBoundaryTransfer()
   {
      this->synchronizingPeriodic = false;
   }

   //TODO: Remove this.
   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   void
   initialize( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const GlobalIndexType numberOfParticles = physicalObject->getNumberOfParticles();
      const GlobalIndexType numberOfAllocatedParticles = physicalObject->particles->getNumberOfAllocatedParticles();
      const GlobalIndexType shiftInMemory = static_cast< int >( ( numberOfAllocatedParticles - numberOfParticles ) / 2 );

      arrangeArray( physicalObject->getParticles()->getPoints(),
                    physicalObject->getParticles()->getPointsSwap(),
                    0,
                    shiftInMemory,
                    numberOfParticles );

      physicalObject->setFirstActiveParticle( shiftInMemory );
      physicalObject->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 );
      physicalObject->particles->setFirstActiveParticle( shiftInMemory );
      physicalObject->particles->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 );
   }

   //TODO: Remove this.
   template< typename PhysicalObjectPointer, typename Array >
   void
   initializeParticleField( PhysicalObjectPointer& physicalObject, Array& array, Array& arraySwap )
   {
      const GlobalIndexType numberOfParticles = physicalObject->getNumberOfParticles();
      const GlobalIndexType numberOfAllocatedParticles = physicalObject->particles->getNumberOfAllocatedParticles();
      const GlobalIndexType shiftInMemory = static_cast< int >( ( numberOfAllocatedParticles - numberOfParticles ) / 2 );

      arrangeArray( array,
                    arraySwap,
                    0,
                    shiftInMemory,
                    numberOfParticles );

   }

   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   void
   applyPeriodicBoundaryConditionPostIntegration( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const VectorType coordinatesDifference = particlesParams.periodicBoundaryDistance;

      auto points_view = physicalObject->getPoints().getView();

      const PairIndexType firstAndLastParticleInFirstBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithLeftPeriodicity );
      const PairIndexType firstAndLastParticleInLastBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithRightPeriodicity );

      auto last = [=] __cuda_callable__ ( int i ) mutable
      {
         VectorType r = points_view[ i ];
         if(  r[ 0 ] > 0.08f )
            points_view[ i ] = r - coordinatesDifference;

      };
      Algorithms::parallelFor< DeviceType >( firstAndLastParticleInLastBlock[ 0 ], firstAndLastParticleInLastBlock[ 1 ] + 1, last );

      auto first = [=] __cuda_callable__ ( int i ) mutable
      {
         VectorType r = points_view[ i ];
         if ( r[ 0 ] < -0.08f )
            points_view[ i ] = r + coordinatesDifference;
      };
      Algorithms::parallelFor< DeviceType >( firstAndLastParticleInFirstBlock[ 0 ], firstAndLastParticleInFirstBlock[ 1 ] + 1, first );
   }

protected:

   //Track the process of periodic BC synchronization.
   bool synchronizingPeriodic = false;

   //Store first and last particle in first and last particle block.
   PairIndexType firstAndLastParticleInFirstBlock = 0;
   PairIndexType firstAndLastParticleInLastBlock = 0;

   //Store first and last active particle for the synchronization process.
   GlobalIndexType firstActiveParticle = 0;
   GlobalIndexType lastActiveParticle = 0;

   //Store number of particles to send between patches during the synchronization.
   GlobalIndexType sizeToCopyFirstBlock = 0;
   GlobalIndexType sizeToCopyLastBlock = 0;

};


} // SPH
} // TNL

