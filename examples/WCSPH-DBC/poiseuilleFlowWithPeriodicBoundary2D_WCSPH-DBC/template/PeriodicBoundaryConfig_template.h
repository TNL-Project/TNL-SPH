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
class PeriodicityLeftBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   /**
    * Define identifier of the open boundary buffer patch.
    */
   std::string identifier = "periodicity-left";

   /**
    * Define geometrical identifier of buffer.
    * - orientation - unit normal buffer orientation [-]
    * - position - referential point of buffer (corner or centre ) [m]
    * - bufferWidth - with of buffer (depends on number of boundary layers) [m]
    */
   VectorType orientation = { placeholderPeriodicityLeftOrientation_xf, placeholderPeriodicityLeftOrientation_yf };
   VectorType position = { placeholderPeriodicityLeftPosition_xf, placeholderPeriodicityLeftPosition_yf };
   VectorType bufferWidth = { placeholderPeriodicityLeftWidth_xf, placeholderPeriodicityLeftWidth_yf };
   VectorType bufferHeight = { placeholderPeriodicityLeftHeigth_xf, placeholderPeriodicityLeftHeigth_yf }; //TODO: Merge with bufferWidth

   /**
    * Max number of particles per cell (required to create buffer zones)
    */
   int numberOfParticlesPerCell = placeholderNumberOfParticlesPerCell;

   /**
    * Shift of particle during the periodicity transfer.
    */
   VectorType shift = { placeholderPeriodicityLeftShiftVector_xf , placeholderPeriodicityLeftShiftVector_yf };

   /**
    * Coordinates of adjecent particle zone.
    */
   VectorType periodicityFirstPoint = { placeholderPeriodicityLeftFirstPoint_xf, placeholderPeriodicityLeftFirstPoint_yf };
   VectorType periodicitySecondPoint = { placeholderPeriodicityLeftSecondPoint_xf, placeholderPeriodicityLeftSecondPoint_yf } ;

};

/**
 * PARAMETERS OF OPEN BOUNDARY PATCH
 *
 * This class is used to store core parameters for inlet boundary patch
 * i.e. inlet or outlet. The values are used only for initialization.
 */
template< typename SPHConfig >
class PeriodicityRightBuffer
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   /**
    * Define identifier of the open boundary buffer patch.
    */
   std::string identifier = "periodicity-right";

   /**
    * Define geometrical identifier of buffer.
    * - orientation - unit normal buffer orientation [-]
    * - position - referential point of buffer (corner or centre ) [m]
    * - bufferWidth - with of buffer (depends on number of boundary layers) [m]
    */
   VectorType orientation = { placeholderPeriodicityRightOrientation_xf, placeholderPeriodicityRightOrientation_yf };
   VectorType position = { placeholderPeriodicityRightPosition_xf, placeholderPeriodicityRightPosition_yf };
   VectorType bufferWidth = { placeholderPeriodicityRightWidth_xf, placeholderPeriodicityRightWidth_yf };
   VectorType bufferHeight = { placeholderPeriodicityRightHeigth_xf, placeholderPeriodicityRightHeigth_yf }; //TODO: Merge with bufferWidth

   /**
    * Max number of particles per cell (required to create buffer zones)
    */
   int numberOfParticlesPerCell = placeholderNumberOfParticlesPerCell;

   /**
    * Shift of particle during the periodicity transfer.
    */
   VectorType shift = { placeholderPeriodicityRightShiftVector_xf , placeholderPeriodicityRightShiftVector_yf };

   /**
    * Coordinates of adjecent particle zone.
    */
   VectorType periodicityFirstPoint = { placeholderPeriodicityRightFirstPoint_xf, placeholderPeriodicityRightFirstPoint_yf };
   VectorType periodicitySecondPoint = { placeholderPeriodicityRightSecondPoint_xf, placeholderPeriodicityRightSecondPoint_yf } ;
};

} //SPH
} //namespace ParticleSystem
} //namespace TNL

