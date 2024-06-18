// SPDX-FileComment: This file is part of TNL - Template Numerical Library (https://tnl-project.org/)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/TypeTraits.h>
#include "customParallelForDetail_ParallelFor1D.h"

/**
 * \brief Namespace for fundamental TNL algorithms
 *
 * It contains algorithms like for-loops, memory operations, (parallel) reduction,
 * multireduction, scan etc.
 */
namespace TNL::ModifiedAlgorithms {

/**
 * \brief Parallel for-loop function for 1D range specified with integral values.
 *
 * \tparam Device is a type of the device where the reduction will be performed.
 * \tparam Begin must be an \e integral type.
 * \tparam End must be an \e integral type.
 *
 * \param begin is the left bound of the iteration range `[begin, end)`.
 * \param end is the right bound of the iteration range `[begin, end)`.
 * \param f is the function to be called in each iteration. Arguments of the
 *          function are the iteration index and arguments from the `args...`
 *          variadic pack.
 * \param launch_config specifies kernel launch parameters.
 * \param args are additional parameters to be passed to the function f.
 *
 * \par Example
 * \include Algorithms/parallelForScalarExample.cpp
 * \par Output
 * \include parallelForScalarExample.out
 */
template< typename Device, typename Begin, typename End, typename Function, typename... FunctionArgs >
std::enable_if_t< std::is_integral_v< Begin > && std::is_integral_v< End > >
parallelFor( const Begin& begin,
             const End& end,
             typename Device::LaunchConfiguration launch_config,
             Function f,
             FunctionArgs... args )
{
   using Index = std::common_type_t< Begin, End >;
   detail::ParallelFor1D< Device >::exec(
      static_cast< Index >( begin ), static_cast< Index >( end ), launch_config, f, args... );
}

/**
 * \brief Parallel for-loop function for 1D range specified with integral values with default launch configuration.
 */
template< typename Device, typename Begin, typename End, typename Function, typename... FunctionArgs >
std::enable_if_t< std::is_integral_v< Begin > && std::is_integral_v< End > >
parallelFor( const Begin& begin, const End& end, Function f, FunctionArgs... args )
{
   typename Device::LaunchConfiguration launch_config;
   parallelFor< Device >( begin, end, launch_config, f, args... );
}

/**
 * \brief Parallel for-loop function for range specified with multi-index values.
 *
 * \tparam Device is a type of the device where the reduction will be performed.
 * \tparam Begin must satisfy the constraints checked by the \ref TNL::IsStaticArrayType type trait.
 * \tparam End must satisfy the constraints checked by the \ref TNL::IsStaticArrayType type trait.
 *
 * \param begin is the left bound of the iteration range `[begin, end)`.
 * \param end is the right bound of the iteration range `[begin, end)`.
 * \param f is the function to be called in each iteration. Arguments of the
 *          function are the iteration multi-index, which is an instance of the
 *          `End` type, and arguments from the `args...` variadic pack.
 * \param launch_config specifies kernel launch parameters.
 * \param args are additional parameters to be passed to the function f.
 *
 * \par Example
 * \include Algorithms/parallelForMultiIndexExample.cpp
 * \par Output
 * \include parallelForMultiIndexExample.out
 */
template< typename Device, typename Begin, typename End, typename Function, typename... FunctionArgs >
std::enable_if_t< IsStaticArrayType< Begin >::value && IsStaticArrayType< End >::value >
parallelFor( const Begin& begin,
             const End& end,
             typename Device::LaunchConfiguration launch_config,
             Function f,
             FunctionArgs... args )
{
   static_assert( std::is_integral_v< typename Begin::ValueType >,
                  "the ValueType of the Begin multi-index must be an integral type" );
   static_assert( std::is_integral_v< typename End::ValueType >,
                  "the ValueType of the End multi-index must be an integral type" );
   static_assert( Begin::getSize() == End::getSize(), "the Begin multi-index must have the same size as the End multi-index" );

   if constexpr( Begin::getSize() == 1 ) {
      parallelFor< Device >( begin.x(), end.x(), launch_config, f, args... );
   }
}

/**
 * \brief Parallel for-loop function for range specified with multi-index values with default launch configuration.
 */
template< typename Device, typename Begin, typename End, typename Function, typename... FunctionArgs >
std::enable_if_t< IsStaticArrayType< Begin >::value && IsStaticArrayType< End >::value >
parallelFor( const Begin& begin, const End& end, Function f, FunctionArgs... args )
{
   typename Device::LaunchConfiguration launch_config;
   parallelFor< Device >( begin, end, launch_config, f, args... );
}

}  // namespace TNL::ModifiedAlgorithms

