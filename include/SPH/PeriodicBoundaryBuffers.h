#include <TNL/Particles/GhostZone.h>

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

   OpenBoundaryConfigType config;
   ParticleZone particleZone;
};

} // SPH
} // TNL

