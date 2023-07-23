#pragma once

#include <TNL/Algorithms/parallelFor.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace utils {

template< typename Array, typename GlobalIndexType >
static void
shiftArray( Array& arrayCopy,
            Array& arrayPaste,
            GlobalIndexType fromPosition,
            GlobalIndexType toPosition,
            GlobalIndexType size )
{
   const auto viewArrayCopy = arrayCopy.getConstView();
   auto viewPaste = arrayPaste.getView();

   auto copyToSwap = [ = ] __cuda_callable__ ( GlobalIndexType i ) mutable
   {
      viewPaste[ toPosition + i ] = viewArrayCopy[ fromPosition + i ];
   };
   Algorithms::parallelFor< typename Array::DeviceType >( 0, size, copyToSwap );

   arrayCopy.swap( arrayPaste );
}

} // utils
} // SPH
} // ParticleSystem
} // TNL

