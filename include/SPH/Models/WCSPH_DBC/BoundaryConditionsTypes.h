#pragma once

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace WCSPH_BCTypes {

struct DBC
{};

struct MDBC
{};

struct MDBCVelocity //not implement yet
{};

struct GWBC //not implemented
{};

struct OpenBoundary
{};

struct Inlet
{};

struct Outlet
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

} // WCSPHSolidWallBCTypes
} // SPH
} // ParticleSystem
} // TNL

