// SPDX-FileComment: This file is part of TNL - Template Numerical Library (https://tnl-project.org/)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/Backend.h>
#include <TNL/Devices/Sequential.h>
#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>
#include <TNL/Math.h>

namespace TNL::ModifiedAlgorithms::detail {

template< typename Device = Devices::Sequential >
struct ParallelFor1D
{
   template< typename Index, typename Function, typename... FunctionArgs >
   static void
   exec( Index begin, Index end, typename Device::LaunchConfiguration launch_config, Function f, FunctionArgs... args )
   {
      for( Index i = begin; i < end; i++ )
         f( i, args... );
   }
};

template<>
struct ParallelFor1D< Devices::Host >
{
   template< typename Index, typename Function, typename... FunctionArgs >
   static void
   exec( Index begin, Index end, Devices::Host::LaunchConfiguration launch_config, Function f, FunctionArgs... args )
   {
#ifdef HAVE_OPENMP
      // Benchmarks show that this is significantly faster compared
      // to '#pragma omp parallel for if( Devices::Host::isOMPEnabled() && end - begin > 512 )'
      if( Devices::Host::isOMPEnabled() && end - begin > 512 ) {
         #pragma omp parallel for
         for( Index i = begin; i < end; i++ )
            f( i, args... );
      }
      else {
         Devices::Sequential::LaunchConfiguration sequential_config;
         ParallelFor1D< Devices::Sequential >::exec( begin, end, sequential_config, f, args... );
      }
#else
      Devices::Sequential::LaunchConfiguration sequential_config;
      ParallelFor1D< Devices::Sequential >::exec( begin, end, sequential_config, f, args... );
#endif
   }
};

template< bool gridStride = true, typename Index, typename Function, typename... FunctionArgs >
__global__
void
ParallelFor1DKernel( Index begin, Index end, Function f, FunctionArgs... args )
{
#if defined( __CUDACC__ ) || defined( __HIP__ )
   Index i = begin + blockIdx.x * blockDim.x + threadIdx.x;
   while( i < end ) {
      f( i, args... );
      if( gridStride )
         i += blockDim.x * gridDim.x;
      else
         break;
   }
#endif
}

template<>
struct ParallelFor1D< Devices::Cuda >
{
   // NOTE: launch_config must be passed by value so that the modifications of
   // blockSize and gridSize do not propagate to the caller
   template< typename Index, typename Function, typename... FunctionArgs >
   static void
   exec( Index begin, Index end, Devices::Cuda::LaunchConfiguration launch_config, Function f, FunctionArgs... args )
   {
      if( end <= begin )
         return;

      launch_config.blockSize.x = 128;
      launch_config.blockSize.y = 1;
      launch_config.blockSize.z = 1;
      launch_config.gridSize.x =
         TNL::min( Backend::getMaxGridXSize(), Backend::getNumberOfBlocks( end - begin, launch_config.blockSize.x ) );
      launch_config.gridSize.y = 1;
      launch_config.gridSize.z = 1;

      if( (std::size_t) launch_config.blockSize.x * launch_config.gridSize.x >= (std::size_t) end - begin ) {
         constexpr auto kernel = ParallelFor1DKernel< false, Index, Function, FunctionArgs... >;
         Backend::launchKernel( kernel, launch_config, begin, end, f, args... );
      }
      else {
         // decrease the grid size and align to the number of multiprocessors
         const int desGridSize = 32 * Backend::getDeviceMultiprocessors( Backend::getDevice() );
         launch_config.gridSize.x =
            TNL::min( desGridSize, Backend::getNumberOfBlocks( end - begin, launch_config.blockSize.x ) );
         constexpr auto kernel = ParallelFor1DKernel< true, Index, Function, FunctionArgs... >;
         Backend::launchKernel( kernel, launch_config, begin, end, f, args... );
      }
   }
};

}  // namespace TNL::ModifiedAlgorithms::detail

