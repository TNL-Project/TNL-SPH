#pragma once

namespace TNL {
namespace SPH {
namespace WCSPH_BCTypes {

struct DBC
{
   static constexpr bool
   integrateInTime()
   {
      return true;
   }
};

struct MDBC
{
   static constexpr bool
   integrateInTime()
   {
      return false;
   }
};

struct MDBCVelocity //not implement yet
{
   static constexpr bool
   integrateInTime()
   {
      return false;
   }
};

struct GWBC //not implemented
{
   static constexpr bool
   integrateInTime()
   {
      return false;
   }
};

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

} // BCTypes
} // SPH
} // TNL

