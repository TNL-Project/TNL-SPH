// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

namespace TNL::SPH {
namespace WCSPH_BCTypes {

struct BI_numeric
{};

struct BIConsistent_numeric
{
   static constexpr bool renormalize = true;
};

struct BIConsistent_numeric_interpolated
{
   static constexpr bool renormalize = true;
};

struct BIConservative_numeric
{
   static constexpr bool renormalize = false;
};

struct None
{};

enum class OpenBoundaryConditionsType
: std::uint8_t
{
    PeriodicBoundary = 0,

    Inlet = 1,

    Outlet = 2
};

} // BCTypes
}  // namespace TNL::SPH

