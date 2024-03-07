#pragma once

#include <string>
#include "ParticleSet.h"
#include "../Particles/GhostZone.h"
#include "OpenBoundaryConfig.h"

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

template< typename ParticlesType,
          typename SPHCaseConfig,
          typename FluidVariables,
          typename BoundaryVariables,
          typename IntegratorVariables >
class PeriodicBoundary
{
   using DeviceType = typename SPHCaseConfig::DeviceType;
   using RealType = typename ParticlesType::RealType;
   using VectorType = typename ParticlesType::PointType;
   using IndexVectorType = typename ParticlesType::PointType;

   using OpenBoundaryFluidType = OpenBoundary< ParticlesType, SPHCaseConfig, FluidVariables, IntegratorVariables >;
   using OpenBoundaryFluidPointer = Pointers::SharedPointer< OpenBoundaryFluidType, DeviceType >;
   using OpenBoundaryBoundaryType = OpenBoundary< ParticlesType, SPHCaseConfig, BoundaryVariables, IntegratorVariables >;
   using OpenBoundaryBoundaryPointer = Pointers::SharedPointer< OpenBoundaryBoundaryType, DeviceType >;
   using OpenBoundaryConfigType = OpenBoundaryConfig< SPHCaseConfig >;

   void
   initialize( TNL::Config::ParameterContainer& parameters, std::string prefix, TNL::Logger& logger )
   {
      //compute domain properetis
      const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
      const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
      const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
      const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

      fluidPeriodicPatch->initialize( parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                      parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                      searchRadius,
                                      gridSize,
                                      domainOrigin,
                                      config.zoneFirstPoint,
                                      config.zoneSecondPoint );

      boundaryPeriodicPatch->initialize( parameters.getParameter< int >( prefix + "numberOfParticles" ),
                                         parameters.getParameter< int >( prefix + "numberOfAllocatedParticles" ),
                                         searchRadius,
                                         gridSize,
                                         domainOrigin,
                                         config.zoneFirstPoint,
                                         config.zoneSecondPoint );
   }

   OpenBoundaryConfigType config;
   OpenBoundaryFluidPointer fluidPeriodicPatch;
   OpenBoundaryBoundaryPointer boundaryPeriodicPatch;
};


} // SPH
} // TNL

