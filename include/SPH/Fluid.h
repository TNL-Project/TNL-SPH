#pragma once

#include "ParticleSet.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem,
          typename SPHCaseConfig,
          typename Variables,
          typename IntegratorVariables >
class Fluid : public ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >
{
   public:
   using BaseType = ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;
   using RealType = typename BaseType::RealType;
   using VariablesPointerType = typename BaseType::VariablesPointerType;

   Fluid() = default;

   Fluid( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : ParticleSet< ParticleSystem, SPHCaseConfig, Variables, IntegratorVariables >( size, sizeAllocated, h, numberOfCells ) {};

   VariablesPointerType&
   getFluidVariables()
   {
      return this->variables;
   }

   const VariablesPointerType&
   getFluidVariables() const
   {
      return this->variables;
   }

   /**
    * \brief Collects this fluid set's geometry into the initialization report.
    *
    * There is no fluid-specific metadata stored on this class beyond the
    * base particle-set state, so this is intentionally identical to the base
    * report (kept as a separate method to have the same call interface as
    * the other particle-set types).
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      return BaseType::writeInitReport();
   }
};

}
}

