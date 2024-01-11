#include "../../../../SPH/SPHTraits.h"

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

/**
 * SENSORS FOR VALUE MEASUREMENT (optional-add)
 *
 * This class is used to initialize sensors for measuring a variable.
 *
 * It is necessary to enter:
 * - const float outputTime - frequency for evaluating a quantity
 * - std::vector< VectorType > - list of sensors coodinates
 */
template< typename SPHConfig >
class MeasuretoolConfigForPressure
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   const float outputTime = 2e-3f;

   std::vector< VectorType > points{
   { 1.60f - SPHConfig::h, 0.003f },
   { 1.60f - SPHConfig::h, 0.015f },
   { 1.60f - SPHConfig::h, 0.03f  },
   { 1.60f - SPHConfig::h, 0.08f  } };
};

/**
 * SENSORS FOR WATER LEVEL MEASUREMENT (optional-add)
 *
 * This class is used to initialize sensors for measuring a water level.
 *
 * It is necessary to enter:
 * - const float outputTime - frequency for evaluating a quantity
 * - std::vector< VectorType > points - list of sensors coodinates
 * - VectorType direction - direction, in which is the water level measured
 * - const RealType startMeasureAtLevel - level, at which the measurement start
 * - const RealType stopMeasureAtLevel - level, at which the measurement stops
 */
template< typename SPHConfig >
class MeasuretoolConfigForWaterLevel
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   const float outputTime = 2e-3f;

   std::vector< VectorType > points{
   { 0.3f,    0.f + SPHConfig::h },
   { 0.865f,  0.f + SPHConfig::h },
   { 1.114f,  0.f + SPHConfig::h },
   { 1.3625f, 0.f + SPHConfig::h } };

   const VectorType direction = { 0.f, 1.f };
   const RealType startMeasureAtLevel = 0.f;
   const RealType stopMeasureAtLevel = 0.4f;
};

} //Measuretool configuration
} //SPH
} //ParticleSystem
} //TNL

