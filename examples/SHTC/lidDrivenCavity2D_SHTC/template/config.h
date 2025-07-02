#include <TNL/Devices/Cuda.h>
using Device = TNL::Devices::Cuda;

#include <TNL/Algorithms/Segments/CSR.h>
#include <TNL/Algorithms/Segments/Ellpack.h>
#include <TNL/Particles/CellIndexer.h>

class ParticleSystemConfig
{
   public:
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;

   using UseWithDomainDecomposition = std::false_type;
   using CellIndexerType = TNL::ParticleSystem::SimpleCellIndex< spaceDimension, std::index_sequence< 0, 1 > >;
   using NeighborListType = TNL::Algorithms::Segments::Ellpack< Device, int >;
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
   static constexpr int numberOfPeriodicBuffers = 0;
};

#include <SPH/Kernels.h>
#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Models/SHTC/IntegrationSchemes/ExplicitEulerScheme.h>
#include <SPH/Models/SHTC/stress.h>
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
   using DiffusiveTerm = TNL::SPH::DiffusiveTerms::MolteniDiffusiveTerm< SPHConfig >;
   using Stress = TNL::SPH::Stress::FluidStress< SPHConfig >;
   using TimeStepping = TNL::SPH::ConstantTimeStep< SPHConfig >;
   using IntegrationScheme = TNL::SPH::IntegrationSchemes::ExplicitEulerScheme< SPHConfig >;
};

using SPHDefs = SPHParams< Device >;
using ParticlesConfig = ParticleSystemConfig;

// particle system
#include <TNL/Particles/ParticlesLinkedList.h>
using ParticleSystemType = TNL::ParticleSystem::ParticlesLinkedList< ParticlesConfig, Device >;

// SPH model
#include <SPH/Models/SHTC/SHTC.h>
using Model = TNL::SPH::SHTC< ParticleSystemType, SPHDefs >;

// SPH simulation type
#include <SPH/SPHMultiset_CFD.h>
using Simulation = TNL::SPH::SPHMultiset_CFD< Model >;

