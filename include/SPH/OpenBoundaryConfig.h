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
   init( TNL::Config::ParameterContainer& parameters, std::string prefix )
   {
      identifier = parameters.getParameter< std::string >( prefix + "identifier" );
      orientation = parameters.getXyz< VectorType >( prefix + "orientation" );
      position = parameters.getXyz< VectorType >( prefix + "position-1" );
      bufferWidth = parameters.getXyz< VectorType >( prefix + "bufferWidth" );
      //bufferHeight =
      pairedPeriodicBuffer = parameters.getParameter< int >( prefix + "paired-periodic-buffer" );
      shift = parameters.getXyz< VectorType >( prefix + "shiftVector" );
      numberOfParticlesPerCell = parameters.getParameter< int >( prefix + "numberOfParticlesPerCell" );

      //TODO: Check what is the shape of the buffer, following lines are valid for planar buffer
      computeZonePoints( parameters, prefix );
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      logger.writeParameter( "Patch orientation:", orientation );
      logger.writeParameter( "Patch referential position:", position );
      logger.writeParameter( "Max. particle cound per zone cell:", numberOfParticlesPerCell );
      logger.writeParameter( "Zone first point:", zoneFirstPoint );
      logger.writeParameter( "Zone second point:", zoneSecondPoint );
   }

   private:

   void
   computeZonePoints( TNL::Config::ParameterContainer& parameters, std::string prefix )
   {
      const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
      const VectorType firstPointOfBufferArea = parameters.getXyz< VectorType >( prefix + "position-1" );
      const VectorType secondPointOfBufferArea = parameters.getXyz< VectorType >( prefix + "position-2" );
      const VectorType bufferAreaDiagonal = secondPointOfBufferArea - firstPointOfBufferArea;
      const VectorType bufferUnitDiagonal = bufferAreaDiagonal / l2Norm( bufferAreaDiagonal );
      //TODO: Ugly, ugly code:
      if( orientation[ 0 ] != 0 ) {
         if( orientation[ 0 ] >= 0. ){
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
         }
         if( orientation[ 0 ] <= 0. ){
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal;
         }
      }
      if( orientation[ 0 ] == 0 ) {
         if( orientation[ 1 ] >= 0. ){
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
         }
         if( orientation[ 1 ] <= 0. ){
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal;
         }
      }
   }


};


}
}

