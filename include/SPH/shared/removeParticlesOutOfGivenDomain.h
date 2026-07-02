#pragma once
#include <TNL/Algorithms/reduce.h>

namespace TNL {
namespace SPH {
namespace customFunctions {

//TODO: During refactoring, this could move to particle class
template< typename FluidPointer >
int
removeParticlesOutOfGivenBox(
   FluidPointer& fluid,
   const typename FluidPointer::ObjectType::IndexVectorType& boxOriginGlobalCoords,
   const typename FluidPointer::ObjectType::IndexVectorType& boxDimensions )
{
   using ParticleSetType = typename FluidPointer::ObjectType;
   using DeviceType = typename ParticleSetType::DeviceType;
   using IndexVectorType = typename ParticleSetType::IndexVectorType;
   using VectorType = typename ParticleSetType::VectorType;
   using RealType = typename ParticleSetType::RealType;

   const VectorType gridRefOrigin = fluid->getParticles()->getGridReferentialOrigin();
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType invSearchRadius = 1.0f / searchRadius;

   auto r_view = fluid->getPoints().getView();

   auto checkParticlePosition = [ = ] __cuda_callable__( int i ) mutable
   {
      const VectorType point = r_view[ i ];

      // if the particle is already removed, skip
      if( point[ 0 ] == FLT_MAX )
         return 0;

      const IndexVectorType cellGlobalCoords = TNL::floor( ( point - gridRefOrigin ) * invSearchRadius );
      const IndexVectorType cellCoords = cellGlobalCoords - boxOriginGlobalCoords;

      // Keep particles inside the extended box [0, boxDimensions) — including
      // overlap cells. Remove only particles truly outside.
      // (Unlike Particles::isInsideDomain which uses <= 0 and >= dim-1,
      // stripping the interior boundary ring as well.)
      bool isInside = true;
      for( int d = 0; d < IndexVectorType::getSize(); d++ )
         if( cellCoords[ d ] < 0 || cellCoords[ d ] >= boxDimensions[ d ] )
            isInside = false;

      if( isInside )
         return 0;
      else {
         r_view[ i ] = FLT_MAX;
         return 1;
      }
   };

   const int numberOfParticlesToRemove =
      TNL::Algorithms::reduce< DeviceType >( 0, fluid->getNumberOfParticles(), checkParticlePosition, TNL::Plus() );
   fluid->getParticles()->setNumberOfParticlesToRemove(
         fluid->getParticles()->getNumberOfParticlesToRemove() + numberOfParticlesToRemove );

   return numberOfParticlesToRemove;
}

}  //namespace customFunctions
}  //namespace SPH
}  //namespace TNL

