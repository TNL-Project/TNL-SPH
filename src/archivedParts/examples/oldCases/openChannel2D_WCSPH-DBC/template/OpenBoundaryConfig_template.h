#include <string>
#include "../../../../SPH/SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

/**
 * PARAMETERS OF OPEN BOUNDARY PATCH
 *
 * This class is used to store core parameters for inlet boundary patch
 * i.e. inlet or outlet. The values are used only for initialization.
 */
template< typename SPHConfig >
class InletBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   /**
    * Define identifier of the open boundary buffer patch.
    */
   std::string identifier = "inlet";

   /**
    * Define geometrical identifier of buffer.
    * - orientation - unit normal buffer orientation [-]
    * - position - referential point of buffer (corner or centre ) [m]
    * - bufferWidth - with of buffer (depends on number of boundary layers) [m]
    */
   VectorType orientation = { placeholderInletOrientation_xf, placeholderInletOrientation_yf };
   VectorType position = { placeholderInletPosition_xf, placeholderInletPosition_yf };
   VectorType bufferWidth = { placeholderInletWidth_xf, placeholderInletWidth_yf };

   /**
    * Define value handling on the buffer. Three options are possible:
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
   std::string rho_bc = "fixed";
   RealType density = placeholderInletDensityf;

   std::string v_bc = "fiexd";
   VectorType velocity = { placeholderInletVelocity_xf, placeholderInletVelocity_yf };
};

/**
 * PARAMETERS OF OPEN BOUNDARY PATCH
 *
 * This class is used to store core parameters for inlet boundary patch
 * i.e. inlet or outlet. The values are used only for initialization.
 */
template< typename SPHConfig >
class OutletBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   /**
    * Define identifier of the open boundary buffer patch.
    */
   std::string identifier = "outlet";

   /**
    * Define geometrical identifier of buffer.
    * - orientation - unit normal buffer orientation [-]
    * - position - referential point of buffer (corner or centre ) [m]
    * - bufferWidth - with of buffer (depends on number of boundary layers) [m]
    */
   VectorType orientation = { placeholderOutletOrientation_xf, placeholderOutletOrientation_yf };
   VectorType position = { placeholderOutletPosition_xf, placeholderOutletPosition_yf };
   VectorType bufferWidth = { placeholderOutletWidth_xf, placeholderOutletWidth_yf };

   /**
    * Define value handling on the buffer. Three options are possible:
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
   std::string rho_bc = "extrapolated";
   RealType density = placeholderOutletDensityf;

   std::string v_bc = "extrapolated";
   VectorType velocity = { placeholderOutletVelocity_xf, placeholderOutletVelocity_yf };

   float extrapolationDetTreshold = 1000.f;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

