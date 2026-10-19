#pragma once

#include <nlohmann/json.hpp>

namespace TNL {
namespace SPH {

template< typename SPHConfig >
void
conigSetupOpenBoundaryPatch( TNL::Config::ConfigDescription& config, std::string prefix )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   config.addDelimiter( "Open boundary patch delimiter." );
   config.addEntry< int >( prefix + "numberOfParticles", "The initial number of fluid particles.", 0. );
   config.addEntry< int >( prefix + "numberOfAllocatedParticles", "The allocated number of fluid particles.", 0. );
   config.addEntry< std::string >( prefix + "particles", "Input fluid particles file path." );
   config.addEntry< std::string >( prefix + "identifier", "Identifier of the open boundary patch.", "empty" );
   config.addEntry< RealType >( prefix + "orientation-x", "Orientation of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "orientation-y", "Orientation of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "orientation-z", "Orientation of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-1-x", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-1-y", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-1-z", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-2-x", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-2-y", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "position-2-z", "Position of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferWidth-x", "Width of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferWidth-y", "Width of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferWidth-z", "Width of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferHeight-x", "Height of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferHeight-y", "Height of the open boundary buffer.", 0 );
   config.addEntry< RealType >( prefix + "bufferHeight-z", "Height of the open boundary buffer.", 0 );
   config.addEntry< int >( prefix + "numberOfParticlesPerCell", "Maximum allowed number of particles per grid cell.", 0 );

   // parameters for periodic boundary conditions
   config.addEntry< int >( prefix + "paired-periodic-buffer", "Index of paired periodic boundary buffer.", -1 );
   config.addEntry< RealType >( prefix + "shiftVector-x", "Shift vector within the periodic boundary conditions.", 0 );
   config.addEntry< RealType >( prefix + "shiftVector-y", "Shift vector within the periodic boundary conditions.", 0 );
   config.addEntry< RealType >( prefix + "shiftVector-z", "Shift vector within the periodic boundary conditions.", 0 );
}

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

   /**
    * Second referential point of the buffer area, kept exactly as loaded from
    * the jsonc `position-2-*` entries (needed verbatim by the initialization
    * report - the zone points below are derived from it, so it cannot be
    * recovered from them afterwards).
    */
   VectorType position2 = 0.f;
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
   init( TNL::Config::ParameterContainer& parameters,
         TNL::Config::ParameterContainer& parametersOpenBoundary,
         std::string prefix )
   {
      identifier = parametersOpenBoundary.getParameter< std::string >( prefix + "identifier" );
      orientation = parametersOpenBoundary.getXyz< VectorType >( prefix + "orientation" );
      position = parametersOpenBoundary.getXyz< VectorType >( prefix + "position-1" );
      position2 = parametersOpenBoundary.getXyz< VectorType >( prefix + "position-2" );
      bufferWidth = parametersOpenBoundary.getXyz< VectorType >( prefix + "bufferWidth" );
      //bufferHeight =
      pairedPeriodicBuffer = parametersOpenBoundary.getParameter< int >( prefix + "paired-periodic-buffer" );
      shift = parametersOpenBoundary.getXyz< VectorType >( prefix + "shiftVector" );
      numberOfParticlesPerCell = parametersOpenBoundary.getParameter< int >( prefix + "numberOfParticlesPerCell" );

      //TODO: Check what is the shape of the buffer, following lines are valid for planar buffer
      computeZonePoints( parameters, parametersOpenBoundary, prefix );
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      logger.writeParameter( "Patch orientation:", orientation );
      logger.writeParameter( "Patch referential position:", position );
      logger.writeParameter( "Max. particle count per zone cell:", numberOfParticlesPerCell );
      logger.writeParameter( "Zone first point:", zoneFirstPoint );
      logger.writeParameter( "Zone second point:", zoneSecondPoint );
      //temp: Logs to periodic bounday buffers
      logger.writeParameter( "Paired periodic boundary patch index: ", pairedPeriodicBuffer );
      logger.writeParameter( "Periodic boundary shift vector: ", shift );
   }

   /**
    * \brief Collects the boundary-patch geometry into the initialization report.
    *
    * Emits the values exactly as loaded from the open-boundary JSONC
    * configuration, with the entry names snake-cased: \e identifier,
    * \e orientation_x/y/z, \e position_1_x/y/z, \e position_2_x/y/z,
    * \e bufferWidth_x/y/z plus the periodic-patch entries
    * \e paired_periodic_buffer and \e shiftVector_x/y/z (both present in the
    * config description for every patch; their defaults report as loaded).
    * Every vector entry carries the three axis components (z = 0 for 2D
    * configurations), the same convention as the \e config section of
    * InitReport.hpp. Model-specific entries (\e type, \e velocity_*, ...)
    * are appended by the derived model configurations.
    *
    * \returns A JSON object with the patch entries as loaded.
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      nlohmann::json report = nlohmann::json::object();
      report[ "identifier" ] = identifier;
      vectorToInitReportEntries( report, "orientation", orientation );
      vectorToInitReportEntries( report, "position_1", position );
      vectorToInitReportEntries( report, "position_2", position2 );
      vectorToInitReportEntries( report, "bufferWidth", bufferWidth );
      report[ "paired_periodic_buffer" ] = pairedPeriodicBuffer;
      vectorToInitReportEntries( report, "shiftVector", shift );
      return report;
   }

   protected:
   /**
    * \brief Serializes a space-dimension-sized static vector as per-axis
    *        entries \<prefix\>_x/_y/_z of a report object.
    *
    * The z entry is emitted as 0 in 2D, so the reported entries always carry
    * the three axis components - mirroring the z-padded vector arrays of the
    * config section and particle-set reports.
    */
   template< typename Vector >
   static void
   vectorToInitReportEntries( nlohmann::json& report, const std::string& prefix, const Vector& v )
   {
      report[ prefix + "_x" ] = v[ 0 ];
      report[ prefix + "_y" ] = v[ 1 ];
      if constexpr( Vector::getSize() > 2 )
         report[ prefix + "_z" ] = v[ 2 ];
      else
         report[ prefix + "_z" ] = RealType( 0 );
   }

   private:

   void
   computeZonePoints( TNL::Config::ParameterContainer& parameters,
                      TNL::Config::ParameterContainer& parametersOpenBoundary,
                      std::string prefix )
   {
      const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
      const VectorType firstPointOfBufferArea = parametersOpenBoundary.getXyz< VectorType >( prefix + "position-1" );
      const VectorType secondPointOfBufferArea = parametersOpenBoundary.getXyz< VectorType >( prefix + "position-2" );
      const VectorType bufferAreaDiagonal = secondPointOfBufferArea - firstPointOfBufferArea;
      const VectorType bufferUnitDiagonal = bufferAreaDiagonal / l2Norm( bufferAreaDiagonal );
      //TODO: Ugly, ugly code:
      if( orientation[ 0 ] != 0 ) {
         if( orientation[ 0 ] >= 0. ){
            //zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal;
            //zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal - bufferWidth * orientation; //this could be more narrow
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + searchRadius * orientation;
         }
         if( orientation[ 0 ] <= 0. ){
            //zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
            //zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal;
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal + searchRadius * orientation;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal - bufferWidth * orientation; //this could be more narrow
         }
      }
      if( orientation[ 0 ] == 0 ) {
         if( orientation[ 1 ] >= 0. ){
            //zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal;
            //zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + bufferWidth * orientation;
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal - bufferWidth[ 0 ] * orientation; //this could be more narrow
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal + searchRadius * orientation;
         }
         if( orientation[ 1 ] <= 0. ){
            zoneFirstPoint = firstPointOfBufferArea - searchRadius * bufferUnitDiagonal + searchRadius * orientation;
            zoneSecondPoint = secondPointOfBufferArea + searchRadius * bufferUnitDiagonal - bufferWidth[ 0 ] * orientation; //this could be more narrow
         }
      }
   }
};


}
}

