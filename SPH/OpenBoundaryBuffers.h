#pragma once

#include <string>
#include "PhysicalObject.h"

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

      parameters.density = config.density;
      parameters.velocity = config.velocity;
   }

   struct OpenBoundaryParameters
   {
      std::string identifier;
      VectorType position;
      VectorType orientation;
      VectorType bufferWidth;

      RealType density;
      VectorType velocity;

   };
   OpenBoundaryParameters parameters;

   GlobalIndexType numberOfFluidParticlesToRemove = 0;

   //zone grid
};

}
}
}

