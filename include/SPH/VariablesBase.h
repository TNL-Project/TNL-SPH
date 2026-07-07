#pragma once

#include "ParticleField.h"

#include <tuple>

namespace TNL {
namespace SPH {

/**
 * \brief CRTP base providing the five lifecycle methods (setSize,
 *        sortVariables, readVariables, writeVariables, synchronizeVariables)
 *        by iterating over Derived::allFields() — a tuple of references to
 *        the Derived's ParticleField members.
 */
template< typename Derived >
class VariablesBase
{
public:
   template< typename Index >
   void
   setSize( const Index& size )
   {
      for_each_field(
         [ & ]( auto& field )
         {
            field.setSize( size );
         } );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      for_each_field(
         [ & ]( auto& field )
         {
            field.reorder( particles );
         } );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      for_each_field(
         [ & ]( auto& field )
         {
            field.read( reader );
         } );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      for_each_field(
         [ & ]( auto& field )
         {
            field.write( writer );
         } );
   }

#ifdef HAVE_MPI
   template< typename Synchronizer, typename DistributedParticlesPointer >
   void
   synchronizeVariables( Synchronizer& synchronizer, DistributedParticlesPointer& distributedParticles )
   {
      for_each_field(
         [ & ]( auto& field )
         {
            field.synchronize( synchronizer, distributedParticles );
         } );
   }
#endif

protected:
   template< typename Functor >
   void
   for_each_field( Functor&& fn )
   {
      std::apply(
         [ & ]( auto&... fields )
         {
            ( fn( fields ), ... );
         },
         static_cast< Derived* >( this )->allFields() );
   }
};

}  // namespace SPH
}  // namespace TNL
