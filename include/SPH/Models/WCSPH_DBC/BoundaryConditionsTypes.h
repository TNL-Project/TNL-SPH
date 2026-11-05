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

/**
 * \todo DRAFT - TNL-SPH draft boundary condition, still unimplemented.
 * Kept intentionally as a documented placeholder (user decision); tracked as follow-up
 * work in Documentation/refactor-followups.md (refactor-followups registry).
 */
struct MDBCVelocity
{
   static constexpr bool
   integrateInTime()
   {
      return false;
   }
};

/**
 * \todo DRAFT - TNL-SPH draft boundary condition, still unimplemented.
 * Kept intentionally as a documented placeholder (user decision); tracked as follow-up
 * work in Documentation/refactor-followups.md (refactor-followups registry).
 */
struct GWBC
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
