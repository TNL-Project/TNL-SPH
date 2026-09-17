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
 *  gamma  ScalarArrayType No     false  true   output only
 *  marker MarkerArrayType Yes    true   false  state, reordered, input only
 *
 * referentialIdx and its swap buffer are index bookkeeping, not physical
 * state: referentialIdx is identity-initialized in setSize and participates
 * only in sorting and output, so it stays outside allFields() and is handled
 * explicitly by the leaf classes (this also keeps initGhostRangeFrom() from
 * ever overwriting it).
 */
template< typename SPHState >
class FluidVariablesBase
{
public:
   using SPHTraitsType = SPHFluidTraits< typename SPHState::SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using MarkerArrayType = typename SPHTraitsType::MarkerArrayType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   ParticleField< ScalarArrayType, Swap::Yes, true, true > rho{ "Density" };
   ParticleField< ScalarArrayType, Swap::No, false, false > drho{ "Drho" };
   ParticleField< ScalarArrayType, Swap::No, false, true > p{ "Pressure" };
   ParticleField< VectorArrayType, Swap::Yes, true, true > v{ "Velocity" };
   ParticleField< VectorArrayType, Swap::No, false, false > a{ "Accel" };
   ParticleField< ScalarArrayType, Swap::No, false, true > gamma{ "Gamma" };
   ParticleField< MarkerArrayType, Swap::Yes, true, false > marker{ "Ptype" };

   GlobalIndexType highestReferentialIdx = 0;
   IndexArrayType referentialIdx;
   IndexArrayType referentialIdx_swap;

   auto
   allFields()
   {
      return std::tie( rho, drho, p, v, a, gamma, marker );
   }

   // public: nvcc rejects extended __cuda_callable__ lambdas in protected methods
   void
   setReferentialSize( const GlobalIndexType& size )
   {
      referentialIdx.setSize( size );
      referentialIdx_swap.setSize( size );
      referentialIdx.forAllElements(
         [] __cuda_callable__( GlobalIndexType i, GlobalIndexType& value )
         {
            value = i;
         } );
      highestReferentialIdx = size;
   }

   template< typename ParticlesPointer >
   void
   sortReferential( ParticlesPointer& particles )
   {
      particles->reorderArray( referentialIdx, referentialIdx_swap );
   }

   template< typename WriterType >
   void
   writeReferential( WriterType& writer )
   {
      writer.template writePointData< IndexArrayType >( referentialIdx, "ReferentialIndex" );
   }
};

/**
 * \brief Fluid variables - leaf class.
 *
 * Inherits the fields from \ref FluidVariablesBase and the lifecycle methods
 * from \ref VariablesBase (CRTP); referential indices are sized, sorted and
 * written explicitly.
 */
template< typename SPHState >
class FluidVariables : public FluidVariablesBase< SPHState >, public VariablesBase< FluidVariables< SPHState > >
{
public:
   using GlobalIndexType = typename FluidVariablesBase< SPHState >::GlobalIndexType;

   void
   setSize( const GlobalIndexType& size )
   {
      VariablesBase< FluidVariables >::setSize( size );
      this->setReferentialSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      VariablesBase< FluidVariables >::sortVariables( particles );
      this->sortReferential( particles );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      VariablesBase< FluidVariables >::writeVariables( writer );
      this->writeReferential( writer );
   }
};

/**
 * \brief Shared field holder for boundary-like particle sets - fluid fields
 * plus wall normal and element size, both reordered and file-carried.
 */
template< typename SPHState >
class BoundaryVariablesBase : public FluidVariablesBase< SPHState >
{
public:
   using FieldBase = FluidVariablesBase< SPHState >;
   using ScalarArrayType = typename FieldBase::ScalarArrayType;
   using VectorArrayType = typename FieldBase::VectorArrayType;

   ParticleField< VectorArrayType, Swap::Yes, true, true > n{ "Normals" };
   ParticleField< ScalarArrayType, Swap::Yes, true, false > elementSize{ "ElementSize" };

   auto
   allFields()
   {
      return std::tuple_cat( FieldBase::allFields(), std::tie( n, elementSize ) );
   }
};

/**
 * \brief Boundary variables - leaf class.
 */
template< typename SPHState >
class BoundaryVariables : public BoundaryVariablesBase< SPHState >, public VariablesBase< BoundaryVariables< SPHState > >
{
public:
   using GlobalIndexType = typename BoundaryVariablesBase< SPHState >::GlobalIndexType;

   void
   setSize( const GlobalIndexType& size )
   {
      VariablesBase< BoundaryVariables >::setSize( size );
      this->setReferentialSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      VariablesBase< BoundaryVariables >::sortVariables( particles );
      this->sortReferential( particles );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      VariablesBase< BoundaryVariables >::writeVariables( writer );
      this->writeReferential( writer );
   }
};

/**
 * \brief Open boundary variables - boundary fields plus two index marks,
 * both sized but neither read, written nor reordered.
 */
template< typename SPHState >
class OpenBoundaryVariables
: public BoundaryVariablesBase< SPHState >, public VariablesBase< OpenBoundaryVariables< SPHState > >
{
public:
   using FieldBase = BoundaryVariablesBase< SPHState >;
   using SPHTraitsType = typename FieldBase::SPHTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   ParticleField< IndexArrayType, Swap::No, false, false > particleMark{ "ParticleMark" };
   ParticleField< IndexArrayType, Swap::No, false, false > receivingParticleMark{ "ReceivingParticleMark" };

   auto
   allFields()
   {
      return std::tuple_cat( FieldBase::allFields(), std::tie( particleMark, receivingParticleMark ) );
   }

   void
   setSize( const GlobalIndexType& size )
   {
      VariablesBase< OpenBoundaryVariables >::setSize( size );
      this->setReferentialSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      VariablesBase< OpenBoundaryVariables >::sortVariables( particles );
      this->sortReferential( particles );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      VariablesBase< OpenBoundaryVariables >::writeVariables( writer );
      this->writeReferential( writer );
   }
};

}  //namespace SPH
}  //namespace TNL
