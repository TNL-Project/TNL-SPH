#include "../../../SPH/SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace MeasuretoolConfiguration {

/**
 * GRID INTERPOLATION (optional-add)
 *
 * This class stores data for data interpolation from particle representation
 * back to continuum, specifically on a uniform structured grid.
 *
 * It is necessary to enter:
 * - VectorType gridOrigin - origin of the interpolation grid.
 * - IndexVectorType gridSize - size of the interpolation grid in steps
 * - VectorType gridSize - size of the grid step.
 */
template< typename SPHConfig >
class GridInterpolationConfig
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using VectorType = typename SPHTraitsType::VectorType;

   VectorType gridOrigin = { 0.f, 0.f };
   IndexVectorType gridSize = { 150, 70 };
   VectorType gridStep = { SPHConfig::h, SPHConfig::h };
};

} //Measuretool configuration
} //SPH
} //ParticleSystem
} //TNL

