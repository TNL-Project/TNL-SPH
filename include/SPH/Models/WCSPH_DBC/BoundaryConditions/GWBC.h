// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <SPH/Models/WCSPH_DBC/Interactions.h>

namespace TNL::SPH {

template< typename ParticleSystem, typename SPHFluidConfig >
template<
   typename FluidPointer,
   typename BoundaryPointer,
   typename ModelConfig,
   typename BCType,
   typename std::enable_if_t< std::is_same_v< BCType, WCSPH_BCTypes::GWBC >, bool > Enabled >
void
WCSPH_DBC< ParticleSystem, SPHFluidConfig >::updateSolidBoundary(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& sphState )
{}

}  // namespace TNL::SPH
