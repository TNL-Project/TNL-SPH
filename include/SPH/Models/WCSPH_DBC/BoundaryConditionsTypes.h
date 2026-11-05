// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

namespace TNL::SPH {
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

struct MDBCVelocity  //not implement yet
{
   static constexpr bool
   integrateInTime()
   {
      return false;
   }
};

struct GWBC  //not implemented
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

enum class OpenBoundaryConditionsType : std::uint8_t
{
   PeriodicBoundary = 0,

   Inlet = 1,

   Outlet = 2
};

}  // namespace WCSPH_BCTypes
}  // namespace TNL::SPH
