#pragma once

#include <string>
#include "ParticleSet.h"
#include <TNL/Particles/GhostZone.h>
#include "OpenBoundaryConfig.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables,
          typename OpenBoundaryConfig >
class OpenBoundary : public ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;
   using OpenBoundaryConfigType = OpenBoundaryConfig;

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
               GlobalIndexType numberOfParticlesPerCell = 75 )
   {
      BaseType::initialize( numberOfParticles, numberOfAllocatedParticles, searchRadius, gridSize, gridOrigin );

      //initialize the zone
      zone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      zone.assignCells( config.zoneFirstPoint, config.zoneSecondPoint, gridSize, gridOrigin, searchRadius );

      //TODO: this is ungly and has to be adjusted somehow
      parameters.identifier = config.identifier;
      parameters.position = config.position;
      parameters.orientation = config.orientation;
      parameters.bufferWidth = config.bufferWidth;
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

   //open boundary config
   OpenBoundaryConfigType config;
   //zone grid
   ParticleZone zone;


};

} // SPH
} // TNL

