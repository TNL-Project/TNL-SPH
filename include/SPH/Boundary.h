#pragma once

#include "ParticleSet.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables >
class Boundary : public ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;

   Boundary() = default;

   Boundary( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >( size, sizeAllocated, h, numberOfCells ) {};

   VariablesPointerType&
   getBoundaryVariables()
   {
      return this->variables;
   }

   const VariablesPointerType&
   getBoundaryVariables() const
   {
      return this->variables;
   }

   /**
    * \brief Collects this boundary set's geometry into the initialization report.
    *
    * The boundary set stores no set-specific metadata besides the base
    * particle-set state (in particular, the number of boundary layers is a
    * case-setup input, not a stored member), so this intentionally matches
    * the base report - the same call interface as the other particle-set
    * types.
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      return BaseType::writeInitReport();
   }
};

}
}

