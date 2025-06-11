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

#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/DiffusiveTerms.h>
#include <SPH/Models/VisousTerms.h>
#include <SPH/Models/BoundaryViscousTerms.h>
#include <SPH/Models/DensityFilters.h>
#include <SPH/Kernels.h>
#include <SPH/Models/WCSPH_BI/BoundaryConditionsTypes.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/VerletScheme.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/SymplecticVerletScheme.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/MidpointScheme.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/RK45Scheme.h>
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
   using ViscousTerm = TNL::SPH::BIViscousTerms::PhysicalViscosity_MGVT< SPHConfig >;
   using BoundaryViscousTerm = TNL::SPH::BoundaryViscousTerms::None< SPHConfig >;
   using EOS = TNL::SPH::EquationsOfState::TaitLinearizedWeaklyCompressibleEOS< SPHConfig >;
   using BCType = TNL::SPH::WCSPH_BCTypes::BIConservative_numeric;
   using TimeStepping = TNL::SPH::ConstantTimeStep< SPHConfig >;
   using IntegrationScheme = TNL::SPH::IntegrationSchemes::MidpointScheme< SPHConfig >;
   using DensityFilter = TNL::SPH::DensityFilters::None;
   //using DensityFilter = TNL::SPH::DensityFilters::ShepardFilter< SPHConfig, KernelFunction >;
};

using SPHDefs = SPHParams< Device >;
using ParticlesConfig = ParticleSystemConfig;

/**
 * Include type of particle system.
 */
#include <TNL/Particles/ParticlesLinkedList.h>
using ParticlesSys = TNL::ParticleSystem::ParticlesLinkedList< ParticlesConfig, Device >;

/**
 * Include particular formulation of SPH method.
 */
#include <SPH/Models/WCSPH_BI/Interactions.h>
using Model = TNL::SPH::WCSPH_BI< ParticlesSys, SPHParams< Device > >;

#include <SPH/shared/ElasticBounce.h>
using BoundaryCorrection = TNL::SPH::ElasticBounceLight< ParticlesSys, SPHDefs::SPHConfig >;

/**
 * Include type of SPH simulation.
 */
#include <SPH/SPHMultiset_CFD.h>
using Simulation = TNL::SPH::SPHMultiset_CFD< Model >;

// Custom post processing tools
#include <SPH/shared/energyEvaluation/energyFields.h>
using EnergyFields = TNL::SPH::WCSPHEnergyFields< SPHDefs >;

