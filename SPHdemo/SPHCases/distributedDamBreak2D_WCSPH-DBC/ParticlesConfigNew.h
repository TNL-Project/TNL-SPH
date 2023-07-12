#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

//#include "../../../Particles/GenerateCellIndex.h"
#include "../../../Particles/GenerateCellIndexFloating.h"
#include "../../../Particles/ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace ParticleSystemConfig {

/**
 * TYPES OF PARTICLE SYSTEM AND NEIGHBOR SEARCH (necessary)
 *
 * This class is used to store parameters necessary for particle system,
 * i.e. data types for quantities and indices. It also defines dimension
 * and attributes of particle system.
 *
 * It is necessary to enter TYPES for:
 * - GlobalIndexType
 * - LocalIndexType
 * - CellIndexType
 * - RealType
 *
 * - CoordinatesType
 * - CellIndexerType - algorithm to index particle for neighborSearch
 * - NeighborListType - NOT USED ATM
 */
template< typename Device >
class ParticleSystemConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;

   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig, std::index_sequence< 1, 0 > >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >; //deprecated
};

/**
 * PARAMETERS OF PARTICLE SYSTEM AND NEIGHBOR SEARCH (necessary)
 *
 * This class is used to store parameters necessary for particle system,
 * neighbor search and SPH simulation initialization. This is basically
 * just placeholder to load the data.
 *
 * It is necessary to enter:
 * - Number of loaded particles and nuber of allocated particles for
 *   all "objects" (i.e. fluid, boundary, inlet buffer, etc.) of the
 *   SPH simulation.
 * -- numberOfParticles - number of loaded fluid particles
 * -- numberOfAllocatedParticles - number of initially allocated fluid
 * -- numberOfBoundaryParticles - number of loaded bound particles
 * -- numberOfAllocatedBoundaryParticles - number of initially allocated bound
 *
 * - searchRadius - define the size of grid cell
 * - gridXsize - define the size of grid in dim X
 * - gridYsize - define the size of grid in dim Y
 * - gridXorigin - define the origin of grid in dim X
 * - gridYorigin - define the origin of grid in dim Y
 */
template< typename ParticleSystemConfig >
class ParticleInitialSetup
{
   public:
   using ParticlesTraitsType = ParticlesTraits< ParticleSystemConfig, typename ParticleSystemConfig::DeviceType >;
   using IndexVectorType = typename ParticlesTraitsType::IndexVectorType;
   using PointType = typename ParticlesTraitsType::PointType;

   int numberOfParticles;
   int numberOfAllocatedParticles;
   int numberOfBoundaryParticles;
   int numberOfAllocatedBoundaryParticles;

   float searchRadius;
   int gridXsize;
   int gridYsize;
   PointType gridOrigin;

   IndexVectorType gridSize;
   int numberOfGridCells;
};

template< typename ParticleSystemConfig >
class SubdomainInitialSetup
{
   public:
   public:
   using ParticlesTraitsType = ParticlesTraits< ParticleSystemConfig, typename ParticleSystemConfig::DeviceType >;

   int particleIdxStart;
   int particleIdxRealStart;
   int particleIdxEnd;
   int particleIdxRealEnd;

   int gridIdxOverlapStar;
   int gridIdxStart;
   int gridIdxOverlapEnd;
   int gridIdxEnd;
};

template< typename ParticleSystemConfig >
class DistributedInitialParticleSetup
{
   public:
   using ParticleInitialSetupType = ParticleInitialSetup< ParticleSystemConfig >;
   using SubdomainInitialSetupType = SubdomainInitialSetup< ParticleSystemConfig >;

   static const int numberOfSubdomains = 2;
   ParticleInitialSetupType particlesParams[ 2 ];
   SubdomainInitialSetupType subdomainParams[ 2 ];

   //Subdomain 1 - Particles Init info
   DistributedInitialParticleSetup()
   {
      particlesParams[ 0 ].numberOfParticles = 22500;
      particlesParams[ 0 ].numberOfAllocatedParticles = 45000;
      particlesParams[ 0 ].numberOfBoundaryParticles = 1656;
      particlesParams[ 0 ].numberOfAllocatedBoundaryParticles = 4824;

      particlesParams[ 0 ].searchRadius = 0.0056568 * 1.001;
      particlesParams[ 0 ].gridXsize = 290;
      particlesParams[ 0 ].gridYsize = 145;
      particlesParams[ 0 ].gridOrigin = { -0.0096768f, -0.0096768f };

      particlesParams[ 0 ].gridSize = { particlesParams[ 0 ].gridXsize, particlesParams[ 0 ].gridYsize };
      particlesParams[ 0 ].numberOfGridCells = particlesParams[ 0 ].gridXsize * particlesParams[ 0 ].gridYsize;

      //Subdomain 1 - Subdomain info
      subdomainParams[ 0 ].particleIdxStart = 0;
      subdomainParams[ 0 ].particleIdxRealStart = 0;
      subdomainParams[ 0 ].particleIdxEnd = 22049;
      subdomainParams[ 0 ].particleIdxRealEnd = 22499;

      subdomainParams[ 0 ].gridIdxOverlapStar = 0;
      subdomainParams[ 0 ].gridIdxStart = 0;
      subdomainParams[ 0 ].gridIdxOverlapEnd = 54;
      subdomainParams[ 0 ].gridIdxEnd = 53;

      //Subdomain 2
      particlesParams[ 1 ].numberOfParticles = 23250;
      particlesParams[ 1 ].numberOfAllocatedParticles = 45000;
      particlesParams[ 1 ].numberOfBoundaryParticles = 3183;
      particlesParams[ 1 ].numberOfAllocatedBoundaryParticles = 4824;

      particlesParams[ 1 ].searchRadius = 0.0056568 * 1.001;
      particlesParams[ 1 ].gridXsize = 290;
      particlesParams[ 1 ].gridYsize = 145;
      particlesParams[ 1 ].gridOrigin = { -0.0096768f, -0.0096768f };

      particlesParams[ 1 ].gridSize = { particlesParams[ 1 ].gridXsize, particlesParams[ 1 ].gridYsize };
      particlesParams[ 1 ].numberOfGridCells = particlesParams[ 1 ].gridXsize * particlesParams[ 1 ].gridYsize;

      //Subdomain 2 - Subdomain info
      subdomainParams[ 1 ].particleIdxStart = 0;
      subdomainParams[ 1 ].particleIdxRealStart = 0;
      subdomainParams[ 1 ].particleIdxEnd = 22049;
      subdomainParams[ 1 ].particleIdxRealEnd = 22499;

      subdomainParams[ 1 ].gridIdxOverlapStar = 53;
      subdomainParams[ 1 ].gridIdxStart = 54;
      subdomainParams[ 1 ].gridIdxOverlapEnd = 290;
      subdomainParams[ 1 ].gridIdxEnd = 290;
   }

};

} //namespace ParticleSystemConfig
} //namespace ParticleSystem
} //namespace TNL

