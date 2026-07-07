#pragma once

#include "SPHTraits.h"

#include <exception>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>

namespace TNL {
namespace SPH {

enum class Swap
{
   No,
   Yes
};

/**
 * \brief One particle variable together with its lifecycle metadata.
 *
 * \tparam ArrayType the array type used to store the variable (e.g.
 *                   \ref SPHFluidTraits::ScalarArrayType,
 *                   \ref SPHFluidTraits::VectorArrayType,
 *                   \ref SPHFluidTraits::IndexArrayType).
 * \tparam UseSwap   Swap::Yes if the field owns a swap buffer used by
 *                   sortVariables(); Swap::No otherwise. Swapped fields are
 *                   also the only ones synchronizeVariables() exchanges.
 * \tparam Readable  if true, readVariables() reads this field from input.
 * \tparam Writable  if true, writeVariables() writes this field to output.
 * \tparam Optional  if true, a missing field in input is warned, not fatal.
 *
 * The reader and writer each expose a single unified call for both scalar and
 * vector arrays, so ParticleField does not carry a scalar/vector kind: read()
 * and write() issue exactly one call.
 */
template< typename ArrayType, Swap UseSwap = Swap::Yes, bool Readable = true, bool Writable = true, bool Optional = false >
class ParticleField
{
public:
   using ValueType = typename ArrayType::ValueType;
   using GlobalIndexType = typename ArrayType::IndexType;

   ParticleField() = default;

   explicit ParticleField( std::string name )
   : name_( std::move( name ) )
   {}

   void
   setName( std::string name )
   {
      name_ = std::move( name );
   }

   const std::string&
   name() const
   {
      return name_;
   }

   void
   setSize( GlobalIndexType n )
   {
      field_.setSize( n );
      if constexpr( has_swap )
         field_swap_.setSize( n );
   }

   auto
   getView()
   {
      return field_.getView();
   }

   auto
   getConstView() const
   {
      return field_.getConstView();
   }

   ArrayType&
   getArray()
   {
      return field_;
   }

   const ArrayType&
   getArray() const
   {
      return field_;
   }

   // implicit conversion to the underlying array — keeps legacy call sites
   // (passing the field directly to a writer/gather/etc.) working unchanged.
   operator ArrayType&()
   {
      return field_;
   }

   operator const ArrayType&() const
   {
      return field_;
   }

   void
   swap( ArrayType& other )
   {
      field_.swap( other );
   }

   ParticleField&
   operator=( const ValueType& v )
   {
      field_ = v;
      return *this;
   }

   void
   fill( const ValueType& v )
   {
      field_ = v;
      if constexpr( has_swap )
         field_swap_ = v;
   }

   template< typename ParticlesPointer >
   void
   reorder( ParticlesPointer& particles )
   {
      if constexpr( has_swap )
         particles->reorderArray( field_, field_swap_ );
   }

   template< typename Synchronizer, typename DistributedParticlesPointer >
   void
   synchronize( Synchronizer& synchronizer, DistributedParticlesPointer& distributedParticles )
   {
      if constexpr( has_swap )
         synchronizer.synchronize( field_, distributedParticles );
   }

   template< typename ReaderType >
   void
   read( ReaderType& reader )
   {
      if constexpr( Readable ) {
         if constexpr( Optional ) {
            try {
               reader.template readParticleVariable< ArrayType >( field_, name_ );
            }
            catch( const std::exception& e ) {
               std::cout << "Warning: Unable to read particle variable '" << name_ << "': " << e.what() << std::endl;
            }
         }
         else {
            reader.template readParticleVariable< ArrayType >( field_, name_ );
         }
      }
   }

   template< typename WriterType >
   void
   write( WriterType& writer )
   {
      if constexpr( Writable )
         writer.template writePointData< ArrayType >( field_, name_ );
   }

private:
   static constexpr bool has_swap = ( UseSwap == Swap::Yes );

   // The swap buffer is only stored when UseSwap == Yes; otherwise an empty
   // placeholder takes no space thanks to [[no_unique_address]].
   struct no_swap_t
   {};
   using swap_storage_t = std::conditional_t< has_swap, ArrayType, no_swap_t >;

   ArrayType field_;
   std::string name_;
   [[no_unique_address]] swap_storage_t field_swap_;
};

}  // namespace SPH
}  // namespace TNL
