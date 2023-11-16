#pragma once
#include "../../OpenBoundaryConfig.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename OpenBoundaryType, typename SPHConfig >
class DBCOpenBoundaryConfig : public OpenBoundaryConfig< OpenBoundaryType, SPHConfig >
{
   public:
   using Base = OpenBoundaryConfig< OpenBoundaryType, SPHConfig >;
   using RealType = typename OpenBoundaryType::RealType;
   using VectorType = typename OpenBoundaryType::VectorType;

   /**
    * Define value handling on the buffer. Three options are possible:
    *
    * - "fixed" sets constat value for given variable. Together with
    *   this, values corresponding with given model needs to be specified.
    *
    * - "profile" sets profile for given variable. Together with this
    *   lambda function specifing the profil along the buffer needs to be specified.
    *
    * - "extrapolated" extrapolates given variable from fluid to boundary
    *   buffer. Corresponding values are set based on initial condition.
    *   Together with this, treshold for the interpolation matrix needs to be specify:
    *   use 1e-3 first order, 1e3 zero order.
    */
   std::string rho_bc = "undefine";
   RealType density = 0.f;
   RealType densityProfile = 0.f;

   std::string v_bc = "undefine";
   VectorType velocity = 0.f;
   VectorType velocityProfile = 0.f;

   RealType extrapolationDetTreshold = 0.f;
};

} // SPH
} // ParticleSystem
} // TNL

