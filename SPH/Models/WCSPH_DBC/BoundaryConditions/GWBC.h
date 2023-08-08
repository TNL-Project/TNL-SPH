#pragma once

#include "../Interactions.h"
#include "../../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleSystem, typename SPHFluidConfig >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename SPHKernelFunction,
          typename DiffusiveTerm,
          typename ViscousTerm,
          typename EOS,
          typename SPHState,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::GWBC >, bool > Enabled >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                                  BoudaryPointer& boundary,
                                                                  SPHState& sphState )
{

}

} // SPH
} // ParticleSystem
} // TNL

