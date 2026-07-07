#pragma once

#include "../../ParticleField.h"
#include "../../VariablesBase.h"
#include "../../SPHTraits.h"

#include <tuple>

#include "BoundaryConditionsTypes.h"

namespace TNL {
namespace SPH {

/**
 * \brief Shared field holder for fluid-like particle sets.
 *
 *  field  array type      swap   read   write  notes
 *  ------ --------------- -----  -----  -----  -----------------------------
 *  rho    ScalarArrayType Yes    true   true   state, reordered + I/O + sync
 *  drho   ScalarArrayType No     false  false  transient, recomputed
 *  p      ScalarArrayType No     false  true   output only
 *  v      VectorArrayType Yes    true   true   state, reordered + I/O + sync
 *  a      VectorArrayType No     false  false  transient, recomputed
 */
template< typename SPHState >
class FluidVariablesBase
{
public:
   using SPHTraitsType = SPHFluidTraits< typename SPHState::SPHConfig >;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   ParticleField< ScalarArrayType, Swap::Yes, true, true > rho{ "Density" };
   ParticleField< ScalarArrayType, Swap::No, false, false > drho{ "Drho" };
   ParticleField< ScalarArrayType, Swap::No, false, true > p{ "Pressure" };
   ParticleField< VectorArrayType, Swap::Yes, true, true > v{ "Velocity" };
   ParticleField< VectorArrayType, Swap::No, false, false > a{ "Accel" };

   auto
   allFields()
   {
      return std::tie( rho, drho, p, v, a );
   }
};

/**
 * \brief Fluid variables — leaf class.
 *
 * Inherits the fields from \ref FluidVariablesBase and the lifecycle methods
 * from \ref VariablesBase (CRTP).
 */
template< typename SPHState >
class FluidVariables : public FluidVariablesBase< SPHState >, public VariablesBase< FluidVariables< SPHState > >
{};

/**
 * \brief Open boundary variables - fluid fields plus two index marks.
 *
 * The marks are sized but neither read nor written.
 */
template< typename SPHState >
class OpenBoundaryVariables : public FluidVariablesBase< SPHState >, public VariablesBase< OpenBoundaryVariables< SPHState > >
{
public:
   using FieldBase = FluidVariablesBase< SPHState >;
   using SPHTraitsType = typename FieldBase::SPHTraitsType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   ParticleField< IndexArrayType, Swap::No, false, false > particleMark{ "ParticleMark" };
   ParticleField< IndexArrayType, Swap::No, false, false > receivingParticleMark{ "ReceivingParticleMark" };

   auto
   allFields()
   {
      return std::tuple_cat( FieldBase::allFields(), std::tie( particleMark, receivingParticleMark ) );
   }
};

/**
 * \brief Primary template — never instantiated; only the DBC and MDBC
 * specializations below are used.
 */
template< typename SPHState >
class BoundaryVariables : public FluidVariables< SPHState >
{};

/**
 * \brief Boundary variables for DBC.
 *
 * \c marker is \c Optional: "Ptype" may be absent from some VTK input files,
 * so \ref ParticleField::read wraps the read in try/catch automatically.
 */
template< typename SPHState >
requires std::same_as< typename SPHState::BCType, WCSPH_BCTypes::DBC > class BoundaryVariables< SPHState >
: public FluidVariablesBase< SPHState >, public VariablesBase< BoundaryVariables< SPHState > >
{
public:
   using FieldBase = FluidVariablesBase< SPHState >;
   using SPHTraitsType = typename FieldBase::SPHTraitsType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   ParticleField< IndexArrayType, Swap::Yes, true, true, true > marker{ "Ptype" };

   auto
   allFields()
   {
      return std::tuple_cat( FieldBase::allFields(), std::tie( marker ) );
   }
};

/**
 * \brief Boundary variables for MDBC — ghost-node boundary.
 *
 * \c ghostNodes and \c n are read and reordered. \c rhoGradRho_gn and
 * \c cMatrix_gn are computation buffers only: sized with the rest, but never
 * read, written, reordered, or synchronized.
 */
template< typename SPHState >
requires std::same_as< typename SPHState::BCType, WCSPH_BCTypes::MDBC > class BoundaryVariables< SPHState >
: public FluidVariablesBase< SPHState >, public VariablesBase< BoundaryVariables< SPHState > >
{
public:
   using FieldBase = FluidVariablesBase< SPHState >;
   using SPHTraitsType = typename FieldBase::SPHTraitsType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using VectorExtendedArrayType = typename SPHTraitsType::VectorExtendedArrayType;
   using MatrixExtendedArrayType = typename SPHTraitsType::MatrixExtendedArrayType;

   ParticleField< VectorArrayType, Swap::Yes, true, false > ghostNodes{ "GhostNodes" };
   ParticleField< VectorArrayType, Swap::Yes, true, false > n{ "Normals" };
   ParticleField< VectorExtendedArrayType, Swap::No, false, false > rhoGradRho_gn;
   ParticleField< MatrixExtendedArrayType, Swap::No, false, false > cMatrix_gn;

   auto
   allFields()
   {
      return std::tuple_cat( FieldBase::allFields(), std::tie( ghostNodes, n, rhoGradRho_gn, cMatrix_gn ) );
   }
};

}  // namespace SPH
}  // namespace TNL

