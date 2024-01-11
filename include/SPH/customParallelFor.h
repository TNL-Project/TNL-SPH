#pragma once

#include <TNL/Devices/Cuda.h>
#include <TNL/Cuda/CheckDevice.h>
#include <TNL/Cuda/DeviceInfo.h>
#include <TNL/Cuda/LaunchHelpers.h>
#include <TNL/Cuda/KernelLaunch.h>
#include <TNL/Math.h>

namespace TNL {
namespace SPH {

template< bool gridStrideX = true, typename Index, typename Function, typename... FunctionArgs >
__global__
void
SPHParallelForKernel( Index start, Index end, Function f, FunctionArgs... args )
{
#ifdef HAVE_CUDA
   Index i = start + blockIdx.x * blockDim.x + threadIdx.x;
   while( i < end ) {
      f( i, args... );
      if( gridStrideX )
         i += blockDim.x * gridDim.x;
      else
         break;
   }
#endif
};

//template< >
//struct ParallelFor< Devices::Cuda >
struct SPHParallelFor
{
   template< typename Index, typename Function, typename... FunctionArgs >
   static void
   exec( Index start, Index end, Function f, FunctionArgs... args )
   {
      Devices::Cuda::LaunchConfiguration launch_config;
      exec( start, end, launch_config, f, args... );
   }

   // NOTE: launch_config must be passed by value so that the modifications of
   // blockSize and gridSize do not propagate to the caller
   template< typename Index, typename Function, typename... FunctionArgs >
   static void
   exec( Index start, Index end, Devices::Cuda::LaunchConfiguration launch_config, Function f, FunctionArgs... args )
   {
      if( end <= start )
         return;

      launch_config.blockSize.x = 128;
      launch_config.blockSize.y = 1;
      launch_config.blockSize.z = 1;
      launch_config.gridSize.x =
         TNL::min( Cuda::getMaxGridXSize(), Cuda::getNumberOfBlocks( end - start, launch_config.blockSize.x ) );

      launch_config.gridSize.y = 1;
      launch_config.gridSize.z = 1;

      if( (std::size_t) launch_config.blockSize.x * launch_config.gridSize.x >= (std::size_t) end - start ) {
         constexpr auto kernel = SPHParallelForKernel< false, Index, Function, FunctionArgs... >;
         Cuda::launchKernel( kernel, launch_config, start, end, f, args... );
      }
      else {
         // decrease the grid size and align to the number of multiprocessors
         const int desGridSize = 32 * Cuda::DeviceInfo::getCudaMultiprocessors( Cuda::DeviceInfo::getActiveDevice() );
         launch_config.gridSize.x = TNL::min( desGridSize, Cuda::getNumberOfBlocks( end - start, launch_config.blockSize.x ) );
         constexpr auto kernel = SPHParallelForKernel< true, Index, Function, FunctionArgs... >;
         Cuda::launchKernel( kernel, launch_config, start, end, f, args... );
      }
   }
};

} // SPH
} // TNL

