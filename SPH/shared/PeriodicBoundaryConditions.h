#include <TNL/Algorithms/parallelFor.h>

#include "../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
class PeriodicBoundaryConditions
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
   static void
   copyPartOfArray( Array& array, GlobalIndexType fromPosition, GlobalIndexType toPosition, GlobalIndexType size )
   {
      //static assert
      auto view = array.getView();

      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename Array >
   static void
   copyPartOfPoints( Array& array,
                     GlobalIndexType fromPosition,
                     GlobalIndexType toPosition,
                     GlobalIndexType size,
                     VectorType coordinatesDifference )
   {
      //static assert
      auto view = array.getView();

      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ] + coordinatesDifference;
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename Array >
   static void
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
   static void
   applyPeriodicBoundaryCondition( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const VectorType coordinatesDifference = particlesParams.periodicBoundaryDistance;

      //find the limits at start and at the end
      const PairIndexType firstAndLastParticleInFirstBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithLeftPeriodicity );
      const PairIndexType firstAndLastParticleInLastBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithRightPeriodicity );

      const GlobalIndexType firstActiveParticle = physicalObject->getFirstActiveParticle();
      const GlobalIndexType lastActiveParticle = physicalObject->getLastActiveParticle();

      //workaround
      //const GlobalIndexType sizeToCopyFirstBlock = firstAndLastParticleInFirstBlock[ 1 ] - firstAndLastParticleInFirstBlock[ 0 ] + 1;
      const GlobalIndexType sizeToCopyFirstBlock = firstAndLastParticleInFirstBlock[ 1 ] - firstActiveParticle + 1;
      //const GlobalIndexType sizeToCopyLastBlock = firstAndLastParticleInLastBlock[ 1 ] - firstAndLastParticleInLastBlock[ 0 ] + 1;
      const GlobalIndexType sizeToCopyLastBlock = lastActiveParticle - firstAndLastParticleInLastBlock[ 0 ] + 1;

      std::cout << "coordinatesDifference: " << coordinatesDifference << std::endl;
      std::cout << "flpfb: " << firstAndLastParticleInFirstBlock << std::endl;
      std::cout << "flplb: " << firstAndLastParticleInLastBlock << std::endl;
      std::cout << "firstActiveParticle: " << firstActiveParticle << std::endl;
      std::cout << "lastActiveParticle: " << lastActiveParticle << std::endl;
      std::cout << "sizeToCopyFirstBlock: " << sizeToCopyFirstBlock << std::endl;
      std::cout << "sizeToCopyLastBlock: " << sizeToCopyLastBlock << std::endl;

      //check array capacitis
      //return;

      //Copy data from periodic patch A to periodic patch B:
      copyPartOfPoints( physicalObject->getPoints(),
                        firstAndLastParticleInFirstBlock[ 0 ],
                        lastActiveParticle + 1,
                        sizeToCopyFirstBlock,
                        coordinatesDifference );

      copyPartOfArray( physicalObject->getVariables()->rho,
                       firstAndLastParticleInFirstBlock[ 0 ],
                       lastActiveParticle + 1,
                       sizeToCopyFirstBlock );

      copyPartOfArray( physicalObject->getVariables()->v,
                       firstAndLastParticleInFirstBlock[ 0 ],
                       lastActiveParticle + 1,
                       sizeToCopyFirstBlock );

      //Copy data from periodic patch B to periodic patch A:
      copyPartOfPoints( physicalObject->getPoints(),
                        firstAndLastParticleInLastBlock[ 0 ],
                        firstActiveParticle - sizeToCopyLastBlock,
                        sizeToCopyLastBlock,
                        ( -1.f ) * coordinatesDifference );

      copyPartOfArray( physicalObject->getVariables()->rho,
                       firstAndLastParticleInLastBlock[ 0 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyLastBlock );

      copyPartOfArray( physicalObject->getVariables()->v,
                       firstAndLastParticleInLastBlock[ 0 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyLastBlock );

      //physicalObject->setFirstActiveParticle( shiftInMemory );
      //physicalObject->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 ); //FIXME!!!

      //std::cout << "setting the start of particle range to size: " << firstActiveParticle - sizeToCopyLastBlock << std::endl;
      //std::cout << "setting the end of particle range to size: " << lastActiveParticle + sizeToCopyLastBlock - 1 << std::endl;
      //physicalObject->particles->setFirstActiveParticle( firstActiveParticle - sizeToCopyLastBlock );
      //physicalObject->particles->setLastActiveParticle( lastActiveParticle + sizeToCopyLastBlock - 1 ); //FIXME!!!
      //const GlobalIndexType newNumberOfParticles = ( lastActiveParticle + sizeToCopyLastBlock - 1 ) - ( firstActiveParticle - sizeToCopyLastBlock ) + 1;
      //physicalObject->particles->setNumberOfParticles( newNumberOfParticles );

      std::cout << "setting the start of particle range to size: " << firstActiveParticle - sizeToCopyLastBlock << std::endl;
      std::cout << "setting the end of particle range to size: " << lastActiveParticle + sizeToCopyFirstBlock - 1 << std::endl;
      physicalObject->particles->setFirstActiveParticle( firstActiveParticle - sizeToCopyLastBlock );
      physicalObject->particles->setLastActiveParticle( lastActiveParticle + sizeToCopyFirstBlock - 1 ); //FIXME!!!
      const GlobalIndexType newNumberOfParticles = ( lastActiveParticle + sizeToCopyFirstBlock - 1 ) - ( firstActiveParticle - sizeToCopyLastBlock ) + 1;
      physicalObject->particles->setNumberOfParticles( newNumberOfParticles );

      //std::cout << "setting the start of particle range to size: " << firstActiveParticle - sizeToCopyLastBlock << std::endl;
      //std::cout << "setting the end of particle range to size: " << lastActiveParticle + sizeToCopyLastBlock - 2 << std::endl;
      //physicalObject->particles->setFirstActiveParticle( firstActiveParticle - sizeToCopyLastBlock );
      //physicalObject->particles->setLastActiveParticle( lastActiveParticle + sizeToCopyLastBlock - 2 ); //FIXME!!!
      //const GlobalIndexType newNumberOfParticles = ( lastActiveParticle + sizeToCopyLastBlock - 2 ) - ( firstActiveParticle - sizeToCopyLastBlock ) + 1;
      //physicalObject->particles->setNumberOfParticles( newNumberOfParticles );
   }

   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   static void
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

      arrangeArray( physicalObject->getVariables()->rho,
                    physicalObject->getVariables()->rho_swap,
                    0,
                    shiftInMemory,
                    numberOfParticles );

      arrangeArray( physicalObject->getVariables()->v,
                    physicalObject->getVariables()->v_swap,
                    0,
                    shiftInMemory,
                    numberOfParticles );

      physicalObject->setFirstActiveParticle( shiftInMemory );
      physicalObject->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 ); //FIXME!!!
      physicalObject->particles->setFirstActiveParticle( shiftInMemory );
      physicalObject->particles->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 ); //FIXME!!!
   }

   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   static void
   applyPeriodicBoundaryConditionPostIntegration( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const VectorType coordinatesDifference = particlesParams.periodicBoundaryDistance;
      auto points_view = physicalObject->getPoints().getView();

      const PairIndexType firstAndLastParticleInFirstBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithLeftPeriodicity );
      const PairIndexType firstAndLastParticleInLastBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithRightPeriodicity );

      //workaround
      const GlobalIndexType lastParticle = physicalObject->getLastActiveParticle();

      std::cout << "RETYPER: firstActiveParticle: " << firstAndLastParticleInFirstBlock << std::endl;
      std::cout << "RETYPER: lastActiveParticle: " << firstAndLastParticleInLastBlock << std::endl;

      auto last = [=] __cuda_callable__ ( int i ) mutable
      {
         VectorType r = points_view[ i ];
         if(  r[ 0 ] >= 0.3f ){
            points_view[ i ] = r - coordinatesDifference;
            printf(" Retyping i = %d, from [ %f, %f ] to [ %f, %f ]. \n", i, r[ 0 ], r[ 1 ] ,(r - coordinatesDifference)[ 0 ], (r - coordinatesDifference)[ 1 ]);
         }
      };
      //Algorithms::parallelFor< DeviceType >( firstAndLastParticleInLastBlock[ 0 ], firstAndLastParticleInLastBlock[ 1 ] + 1, last );
      Algorithms::parallelFor< DeviceType >( firstAndLastParticleInLastBlock[ 0 ], lastParticle, last );

      auto first = [=] __cuda_callable__ ( int i ) mutable
      {
         VectorType r = points_view[ i ];
         if ( r[ 0 ] <= 0.f )
            points_view[ i ] = r + coordinatesDifference;
      };
      Algorithms::parallelFor< DeviceType >( firstAndLastParticleInFirstBlock[ 0 ], firstAndLastParticleInFirstBlock[ 1 ] + 1, first );
   }

protected:

   //PairIndexType firstAndLastParticleInLastBlock;
   //PairIndexType lastAndLastParticleInLastBlock;

};


} // SPH
} // ParticleSystem
} // TNL

