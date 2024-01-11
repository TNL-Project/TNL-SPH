#include <string>
#include <SPH/SPHTraits.h>

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
   VectorType orientation = { 1.0f, 0.0f };
   VectorType position = { 0.101f, 0.002f };
   VectorType bufferWidth = { 0.008f, 0.0f };

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
   RealType density = 1000.0f;

   std::string v_bc = "fiexd";
   VectorType velocity = { 1.0f, 0.0f };
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
   VectorType orientation = { -1.0f, 0.0f };
   VectorType position = { 0.701f, 0.002f };
   VectorType bufferWidth = { 0.008f, 0.0f };

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
   RealType density = 1000.0f;

   std::string v_bc = "extrapolated";
   VectorType velocity = { 1.5f, 0.0f };

   float extrapolationDetTreshold = 1000.f;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

