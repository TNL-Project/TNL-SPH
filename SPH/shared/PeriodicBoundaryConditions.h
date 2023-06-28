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

   using ParirIndexType = typename ParticlesType::ParirIndexType;

   template< typename ArrayView >
   void
   copyPartOfArray( ArrayView& view, GlobalIndexType fromPosition, GlobalIndexType toPosition, GlobalIndexType size )
   {
      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ];
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename ArrayView >
   void
   copyPartOfPoints( ArrayView& view,
                     GlobalIndexType fromPosition,
                     GlobalIndexType toPosition,
                     GlobalIndexType size,
                     VectorType coordinatesDifference )
   {
      //static assert
      auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
      {
         view[ toPosition + i ] = view[ fromPosition + i ] + coordinatesDifference;
      };
      Algorithms::parallelFor< DeviceType >( 0, size, copyToSwap );
   }

   template< typename Array >
   void
   arrangeArray( Array& arrayCopy, Array& arrayPaste,
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
   applyPeriodicBoundaryCondition( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const VectorType coordinatesDifference = particlesParams.periodicBoundaryDistance;

      //find the limits at start and at the end
      const ParirIndexType firstAndLastParticleInFirstBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithLeftPeriodicity );
      const ParirIndexType firstAndLastParticleInLastBlock = physicalObject->particles->getFirstLastParticleInColumnOfCells(
            particlesParams.indexOfColumnWithRightPeriodicity );

      const GlobalIndexType firstActiveParticle = physicalObject->getFirstActiveParticle();
      const GlobalIndexType lastActiveParticle = physicalObject->getLastActiveParticle();

      const GlobalIndexType sizeToCopyFirstBlock = firstAndLastParticleInFirstBlock[ 1 ] - firstAndLastParticleInFirstBlock[ 0 ];
      const GlobalIndexType sizeToCopyLastBlock = firstAndLastParticleInLastBlock[ 1 ] - firstAndLastParticleInLastBlock[ 0 ];

      //check array capacitis

      //Copy data from periodic patch A to periodic patch B:
      copyPartOfPoints( physicalObject.getPoints().getView(),
                       firstAndLastParticleInFirstBlock[ 0 ],
                       firstActiveParticle,
                       sizeToCopyFirstBlock,
                       coordinatesDifference );

      copyPartOfArray( physicalObject.getVariables()->rho.getView(),
                       firstAndLastParticleInFirstBlock[ 0 ],
                       firstActiveParticle,
                       sizeToCopyFirstBlock );

      copyPartOfArray( physicalObject.getVariables()->v.getView(),
                       firstAndLastParticleInFirstBlock[ 0 ],
                       firstActiveParticle,
                       sizeToCopyFirstBlock );

      //Copy data from periodic patch B to periodic patch A:
      copyPartOfPoints( physicalObject.getPoints().getView(),
                       firstAndLastParticleInFirstBlock[ 1 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyFirstBlock,
                       ( -1.f ) * coordinatesDifference );

      copyPartOfArray( physicalObject.getVariables()->rho.getView(),
                       firstAndLastParticleInFirstBlock[ 1 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyFirstBlock );

      copyPartOfArray( physicalObject.getVariables()->v.getView(),
                       firstAndLastParticleInFirstBlock[ 1 ],
                       firstActiveParticle - sizeToCopyLastBlock,
                       sizeToCopyFirstBlock );

   }

   template< typename PhysicalObjectPointer, typename ParticlesConfig >
   void
   initialize( PhysicalObjectPointer& physicalObject, ParticlesConfig& particlesParams )
   {
      const GlobalIndexType numberOfParticles = physicalObject->getNumberOfParticles();
      const GlobalIndexType numberOfAllocatedParticles = physicalObject->particles->getNumberOfAllocatedParticles();
      const GlobalIndexType shiftInMemory = static_cast< int >( ( numberOfAllocatedParticles - numberOfParticles ) / 2 );

      arrangeArray( physicalObject.getParticles()->getPoints().getView(),
                    physicalObject.getParticles()->getPoints().getView(),
                    0,
                    shiftInMemory,
                    numberOfParticles );

      arrangeArray( physicalObject.getVariables()->rho.getView(),
                    physicalObject.getVariables()->rho_swap.getView(),
                    0,
                    shiftInMemory,
                    numberOfParticles );

      arrangeArray( physicalObject.getVariables()->v.getView(),
                    physicalObject.getVariables()->v_swap.getView(),
                    0,
                    shiftInMemory,
                    numberOfParticles );
   }


};


} // SPH
} // ParticleSystem
} // TNL

