#pragma once

#include <thrust/execution_policy.h>
#include <TNL/Devices/Cuda.h>
#include <TNL/Devices/Host.h>

namespace TNL {
namespace Thrust {

template< typename DeviceType >
struct ThrustExecutionPolicySelector;

template<>
struct ThrustExecutionPolicySelector< TNL::Devices::Host >
{
   using type = typename thrust::detail::host_t;
};

template<>
struct ThrustExecutionPolicySelector< TNL::Devices::Cuda >
{
   using type = typename thrust::detail::device_t;
};

template< typename DeviceType >
using ThrustExecutionPolicy = typename ThrustExecutionPolicySelector< DeviceType >::type;

} // Thrust
} // TNL

