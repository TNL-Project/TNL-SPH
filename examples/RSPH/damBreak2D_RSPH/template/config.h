#include <TNL/Devices/Cuda.h>
using Device = TNL::Devices::Cuda;

#include <TNL/Containers/StaticVector.h>
#include <TNL/Particles/CellIndexer.h>
#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>

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

   using UseWithDomainDecomposition = std::false_type;
   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, ParticleSystemConfig, std::index_sequence< 0, 1 > >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

template< typename Device >
class SPHConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;
   static constexpr int numberOfBoundaryBuffers = 0;
   static constexpr int numberOfPeriodicBuffers = 0;
};

#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/RiemannSolvers.h>
#include <SPH/Kernels.h>
#include <SPH/Models/RSPH/IntegrationSchemes/VerletScheme.h>
#include <SPH/TimeStep.h>

/**
 * Particle system reader.
 */
#include <TNL/Particles/Readers/VTKReader.h>
#include <TNL/Particles/Writers/VTKWriter.h>
#include <TNL/Particles/Readers/readSPHSimulation.h>

template< typename Device >
class SPHParams
{
public:
   using SPHConfig = SPHConfig< Device >;

   using KernelFunction = TNL::SPH::KernelFunctions::WendlandKernel< SPHConfig >;
   using RiemannSolver = TNL::SPH::RiemannSolvers::RoeLinearized< SPHConfig >;
   using EOS = TNL::SPH::EquationsOfState::TaitWeaklyCompressibleEOS< SPHConfig >;
   using TimeStepping = TNL::SPH::ConstantTimeStep< SPHConfig >;
   using IntegrationScheme = TNL::SPH::IntegrationSchemes::VerletScheme< SPHConfig >;
};

using SPHDefs = SPHParams< Device >;
using ParticlesConfig = ParticleSystemConfig< Device >;

// particle system
#include <TNL/Particles/ParticlesLinkedList.h>
using ParticleSystemType = TNL::ParticleSystem::ParticlesLinkedList< ParticlesConfig, Device >;

// SPH model
#include <SPH/Models/RSPH/Interactions.h>
using Model = TNL::SPH::RSPH< ParticleSystemType, SPHDefs >;

// SPH simulation type
#include <SPH/SPHMultiset_CFD.h>
using Simulation = TNL::SPH::SPHMultiset_CFD< Model >;

