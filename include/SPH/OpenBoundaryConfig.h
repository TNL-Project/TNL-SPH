#pragma once

namespace TNL {
namespace SPH {

template< typename SPHConfig >
class OpenBoundaryConfig
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   OpenBoundaryConfig() = default;

   /**
    * Define identifier of the open boundary buffer patch.
    */
   std::string identifier = "empty";

   /**
    * Define geometrical identifier of buffer.
    * - orientation - unit normal buffer orientation [-]
    * - position - referential point of buffer (corner or centre ) [m]
    * - bufferWidth - with of buffer (depends on number of boundary layers) [m]
    */
   VectorType orientation = 0.f;
   VectorType position = 0.f;
   VectorType bufferWidth = 0.f;
   VectorType bufferHeight = 0.f; //TODO: Merge with bufferWidth

   /**
    * Shift of particle during the periodicity transfer.
    */
   int pairedPeriodicBuffer = -1;
   VectorType shift = 0.f;

   /**
    * Max number of particles per cell (required to create buffer zones)
    */
   int numberOfParticlesPerCell = 0;

   /**
    * Coordinates of adjecent particle zone.
    */
   VectorType zoneFirstPoint = 0.;
   VectorType zoneSecondPoint = 0. ;

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      logger.writeParameter( "Patch orientation:", orientation );
      logger.writeParameter( "Patch referential position:", position );
      logger.writeParameter( "Max. particle cound per zone cell:", numberOfParticlesPerCell );
      logger.writeParameter( "Zone first point:", zoneFirstPoint );
      logger.writeParameter( "Zone second point:", zoneSecondPoint );
   }



};


}
}

