#pragma once

#include <string>
#include "ParticleSet.h"
#include "../Particles/GhostZone.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables >
class OpenBoundary : public ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;

   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using VectorType = typename SPHTraitsType::VectorType;

   using ParticleZone = ParticleZone< typename ParticleSystem::Config >;

   //remove
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;

   OpenBoundary() = default;

   OpenBoundary( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >( size, sizeAllocated, h, numberOfCells ) {};

   void
   initialize( int numberOfParticles,
               int numberOfAllocatedParticles,
               RealType searchRadius,
               IndexVectorType gridSize,
               VectorType gridOrigin,
               VectorType zoneBoxFirstPoint,
               VectorType zoneBoxSecondPoit,
               GlobalIndexType numberOfParticlesPerCell = 15 )
   {
      BaseType::initialize( numberOfParticles, numberOfAllocatedParticles, searchRadius, gridSize, gridOrigin );

      //initialize the zone
      zone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      zone.assignCells( zoneBoxFirstPoint, zoneBoxSecondPoit, gridSize, gridOrigin, searchRadius );
   }


   VariablesPointerType&
   getOpenBoundaryVariables()
   {
      return this->variables;
   }

   const VariablesPointerType&
   getOpenBoundaryVariables() const
   {
      return this->variables;
   }

   template< typename Params >
   void
   readOpenBoundaryParameters( Params config )
   {
      parameters.identifier = config.identifier;
      parameters.position = config.position;
      parameters.orientation = config.orientation;
      parameters.bufferWidth = config.bufferWidth;
   }

   template< typename Params, typename SPHParams >
   void
   readOpenBoundaryParameters( Params config, SPHParams params )
   {
      parameters.identifier = config.identifier;
      parameters.position = config.position;
      parameters.orientation = config.orientation;
      parameters.bufferWidth = config.bufferWidth;

      //initialize adjecent zone
      zone.setNumberOfParticlesPerCell( config.numberOfParticlesPerCell );
      VectorType zoneBoxFirstPoit = config.periodicityFirstPoint;
      //VectorType zoneBoxSecondPoit = config.periodicitySecondPoint + config.bufferWidth * config.orientation;
      VectorType zoneBoxSecondPoit = config.periodicitySecondPoint;
      zone.template assignCells< typename ParticleSystem::CellIndexer >(
            zoneBoxFirstPoit, zoneBoxSecondPoit, params.gridSize, params.gridOrigin, params.searchRadius  );

   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      BaseType::writeProlog( logger );
      zone.writeProlog( logger );
   }

   struct OpenBoundaryParameters
   {
      std::string identifier;
      VectorType position;
      VectorType orientation;
      VectorType bufferWidth;
   };
   OpenBoundaryParameters parameters;

   GlobalIndexType numberOfFluidParticlesToRemove = 0;

   //zone grid
   ParticleZone zone;

};

}
}

