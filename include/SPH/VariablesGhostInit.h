#pragma once

#include "VariablesBase.h"

#include <TNL/Algorithms/parallelFor.h>

#include <tuple>
#include <utility>

namespace TNL {
namespace SPH {
namespace detail {

template< typename Field, typename Index >
void
initGhostField( Field& dst, Field& src, const Index& dstBegin, const Index& count )
{
   auto dstView = dst.getView();
   auto srcView = src.getConstView();
   using DeviceType = typename std::decay_t< decltype( dstView ) >::DeviceType;

   // if constexpr must stay outside the extended lambdas - nvcc rejects it inside
   if constexpr( Field::readable ) {
      Algorithms::parallelFor< DeviceType >(
         0,
         count,
         [ = ] __cuda_callable__( Index k ) mutable
         {
            dstView[ dstBegin + k ] = srcView[ k ];
         } );
   }
   else {
      Algorithms::parallelFor< DeviceType >(
         0,
         count,
         [ = ] __cuda_callable__( Index k ) mutable
         {
            dstView[ dstBegin + k ] = 0;
         } );
   }
}

template< typename DstTuple, typename SrcTuple, typename Index, std::size_t... I >
void
initGhostFields( DstTuple& dst, SrcTuple& src, const Index& dstBegin, const Index& count, std::index_sequence< I... > )
{
   ( initGhostField( std::get< I >( dst ), std::get< I >( src ), dstBegin, count ), ... );
}

}  // namespace detail

/**
 * \brief Initialize the variable range [dstBegin, dstBegin + count) of `dst` from the
 * leading range of a compatible, freshly-read variable set `src`.
 *
 * Fields flagged ParticleField::readable (carried by the input file, exactly the
 * fields readVariables() reads) are copied; all other fields start at zero instead
 * of uninitialized memory. Members outside allFields() (e.g. referential indices,
 * identity-initialized by setSize) are never touched. Positions are not part of
 * the variable sets and must be copied by the caller separately.
 */
template< typename DstVariables, typename SrcVariables, typename Index >
void
initGhostRangeFrom( DstVariables& dst, SrcVariables& src, const Index& dstBegin, const Index& count )
{
   auto dstFields = dst.allFields();
   auto srcFields = src.allFields();
   static_assert( std::tuple_size_v< decltype( dstFields ) > == std::tuple_size_v< decltype( srcFields ) >,
                  "Source and destination variable lists must have the same number of fields." );
   detail::initGhostFields(
      dstFields, srcFields, dstBegin, count, std::make_index_sequence< std::tuple_size_v< decltype( dstFields ) > >{} );
}

}  // namespace SPH
}  // namespace TNL
