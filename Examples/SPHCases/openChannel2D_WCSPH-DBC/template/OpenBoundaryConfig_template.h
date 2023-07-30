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
 *
 * It is necessary to enter:
 *
 * - orientation_x - x component of normal buffer vector
 * - orientation_y - y component of normal buffer vector
 * - velocity_x - initial x component of open boundary patch velocity
 * - velocity_y - initial x component of open boundary patch velocity
 * - position_x - referential position of open boundary buffer TODO: Move to centre.
 * - inlet_density - referential position of open boundary buffer TODO: Move to centre.
 * - bufferWidth_x - width of buffer - dependent on number of layers
 * - bufferWidth_y - width of buffer - dependent on number of layers
 *
 */
template< typename SPHConfig >
class InletBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   std::string identifier = "inlet";
   VectorType orientation = { placeholderInletOrientation_xf, placeholderInletOrientation_yf };
   VectorType position = { placeholderInletPosition_xf, placeholderInletPosition_yf };
   VectorType bufferWidth = { placeholderInletWidth_xf, placeholderInletWidth_yf };
   RealType bufferEdge = placeholderInletBufferEdgef; //TODO: Remove, deprecated

   VectorType velocity = { placeholderInletVelocity_xf, placeholderInletVelocity_yf };
   RealType density = placeholderInletDensityf;

   std::string waterLevelHandling = "defined";
   RealType waterLevel = 0.15;
};

template< typename SPHConfig >
class OutletBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   std::string identifier = "outlet";
   VectorType orientation = { placeholderOutletOrientation_xf, placeholderOutletOrientation_yf };
   VectorType position = { placeholderOutletPosition_xf, placeholderOutletPosition_yf };
   VectorType bufferWidth = { placeholderOutletWidth_xf, placeholderOutletWidth_yf };
   RealType bufferEdge = placeholderOutletBufferEdgef; //TODO: Remove, deprecated

   VectorType velocity = { placeholderOutletVelocity_xf, placeholderOutletVelocity_yf };
   RealType density = placeholderOutletDensityf;

   std::string waterLevelHandling = "defined";
   RealType waterLevel = 0.15;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

