#pragma once

#include <TNL/Particles/GhostZone.h>

#include <nlohmann/json.hpp>

namespace TNL {
namespace SPH {

template< typename ParticlesType, typename OpenBoundaryConfigType >
class PeriodicBoundary
{
   public:
   using DeviceType = typename ParticlesType::DeviceType;
   using RealType = typename ParticlesType::RealType;
   using VectorType = typename ParticlesType::PointType;
   using GlobalIndexType = typename ParticlesType::GlobalIndexType;
   using CoordinatesType = typename ParticlesType::PointType;

   using ParticleZone = TNL::Particles::ParticleZone< typename ParticlesType::Config, DeviceType >;

   void
   initialize( RealType searchRadius,
               CoordinatesType dimensions,
               VectorType origin,
               GlobalIndexType numberOfParticlesPerCell = 75 ) //TODO: Move this to config params.
   {
      //initialize the zone
      particleZone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      particleZone.assignCells( config.zoneFirstPoint, config.zoneSecondPoint, dimensions, origin, searchRadius );
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      config.writeProlog( logger );
      particleZone.writeProlog( logger );
   }

   /**
    * \brief Collects this periodic boundary patch into the initialization
    *        report.
    *
    * Emits the patch configuration exactly as loaded (identifier, geometry,
    * paired buffer index and shift vector) under \e config together with the
    * current zone state \e zone_cell_count / \e zone_particle_count - the
    * zone particle list is filled during the neighbor search, so the
    * particle count is 0 right after initialization, mirroring the two
    * \e writeProlog parameter lines of the particle zone.
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      nlohmann::json report = nlohmann::json::object();
      report[ "config" ] = config.writeInitReport();
      report[ "zone_cell_count" ] = particleZone.getNumberOfCells();
      report[ "zone_particle_count" ] = particleZone.getNumberOfParticles();
      return report;
   }

   OpenBoundaryConfigType config;
   ParticleZone particleZone;
};

} // SPH
} // TNL

