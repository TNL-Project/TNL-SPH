#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace WCSPH_BCTypes {

enum class OpenBoundaryConditionsType
: std::uint8_t
{
    PeriodicBoundary = 0,

    Inlet = 1,

    Outlet = 2
};

} // WCSPHSolidWallBCTypes
} // SPH
} // ParticleSystem
} // TNL

