#pragma once

#include "../Interactions.h"
#include "../../../customParallelFor.h"

namespace TNL {
namespace SPH {

template< typename ParticleSystem, typename SPHFluidConfig >
template< typename FluidPointer,
          typename BoudaryPointer,
          typename ModelConfig,
          typename BCType,
          typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::GWBC >, bool > Enabled >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig >::updateSolidBoundary( FluidPointer& fluid,
                                                                  BoudaryPointer& boundary,
                                                                  ModelParams& sphState )
{

}

} // SPH
} // TNL

