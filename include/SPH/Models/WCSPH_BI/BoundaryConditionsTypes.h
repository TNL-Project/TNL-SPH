#pragma once

namespace TNL {
namespace SPH {
namespace WCSPH_BCTypes {

struct BI_numeric
{};

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
} // SPH
} // TNL

