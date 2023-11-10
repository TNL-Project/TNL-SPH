#pragma once

#include <string>
#include "PhysicalObject.h"
#include "../Particles/GhostZone.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables >
class OpenBoundary : public PhysicalObject< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = PhysicalObject< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;

   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using VectorType = typename SPHTraitsType::VectorType;

   using ParticleZone = ParticleZone< typename ParticleSystem::Config >;

   OpenBoundary( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : PhysicalObject< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >( size, sizeAllocated, h, numberOfCells ) {};


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
}

