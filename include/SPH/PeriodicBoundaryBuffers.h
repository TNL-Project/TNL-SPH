#include "../Particles/GhostZone.h"

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
   using IndexVectorType = typename ParticlesType::PointType;

   using ParticleZoneType = ParticleZone< typename ParticlesType::Config >;

   void
   initialize( TNL::Config::ParameterContainer& parameters,
               std::string prefix,
               RealType searchRadius,
               IndexVectorType gridSize,
               VectorType domainOrigin,
               GlobalIndexType numberOfParticlesPerCell = 15 ) //TODO: Move this to config params.
   {
      config.init( parameters, prefix );

      //initialize the zone
      particleZone.setNumberOfParticlesPerCell( numberOfParticlesPerCell );
      particleZone.assignCells( config.zoneFirstPoint, config.zoneSecondPoint, gridSize, domainOrigin, searchRadius );
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      config.writeProlog( logger );
      particleZone.writeProlog( logger );
   }

   OpenBoundaryConfigType config;
   ParticleZoneType particleZone;
};

} // SPH
} // TNL

