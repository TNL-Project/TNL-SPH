#pragma once

#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TNL/Algorithms/Segments/Ellpack.h>
#include <TNL/Devices/Host.h>
#include <TNL/Particles/CellIndexer.h>
#include <TNL/Particles/ParticlesCellList.h>
#include <TNL/Particles/Readers/VTKReader.h>
#include <TNL/Particles/Writers/VTKWriter.h>

#include <SPH/shared/thrustExecPolicySelector.h>
#include <SPH/Kernels.h>
#include <SPH/solvers/SolverMultiSetBlockMultiresolutionLocalTimeStepping.h>

namespace {

using Device = TNL::Devices::Host;

class ParticlesConfig
{
public:
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;
   using UseWithDomainDecomposition = std::true_type;
   using CellIndexerType = TNL::Particles::SimpleCellIndex< spaceDimension, std::index_sequence< 0, 1 > >;
   using NeighborListType = TNL::Algorithms::Segments::Ellpack< Device, int >;
};

class CaseConfig
{
public:
   using DeviceType = Device;
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = 2;
};

using Traits = TNL::SPH::SPHFluidTraits< CaseConfig >;

class Variables
{
public:
   void
   setSize( const int size )
   {
      rho.setSize( size );
      p.setSize( size );
      referentialIdx.setSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& )
   {}

   template< typename Writer >
   void
   writeVariables( Writer& )
   {}

   Traits::ScalarArrayType rho;
   Traits::ScalarArrayType p;
   Traits::IndexArrayType referentialIdx;
   int level = -1;
};

struct Probe
{
   std::vector< std::string > events;
   std::array< float, 2 > ghostAtInteraction{ -999.f, -999.f };
};

struct ModelParams
{
   using SPHConfig = CaseConfig;
   using KernelFunction = TNL::SPH::KernelFunctions::WendlandKernel< CaseConfig >;
   using TimeStepping = TNL::SPH::ConstantTimeStep< CaseConfig >;
};

struct IntegrationScheme
{};

struct OpenBoundaryModel
{};

using Particles = TNL::Particles::ParticlesCellList< ParticlesConfig, Device >;
using OpenBoundaryConfig = TNL::SPH::OpenBoundaryConfig< CaseConfig >;

float
sourceValue( const int level )
{
   return 101.f * ( level + 1 );
}

class OrderingModel
{
public:
   using DeviceType = Device;
   using ModelConfigType = ::ModelParams;
   using ModelParams = ::ModelParams;
   using ParticlesType = Particles;
   using IntegrationSchemeType = IntegrationScheme;
   using IntegrationSchemeVariables = Variables;
   using SPHConfig = CaseConfig;
   using FluidVariables = Variables;
   using BoundaryVariables = Variables;
   using OpenBoundaryVariables = Variables;
   using OpenBoundaryConfig = ::OpenBoundaryConfig;
   using OpenBoundaryModel = ::OpenBoundaryModel;

   template< typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, ModelParams& )
   {
      auto& rho = physicalObject->getVariables()->rho;
      auto& p = physicalObject->getVariables()->p;
      for( int i = 0; i < rho.getSize(); i++ )
         p.setElement( i, rho.getElement( i ) );
   }

   template< typename FluidPointer, typename BoundaryPointer >
   void
   updateSolidBoundary( FluidPointer&, BoundaryPointer& boundary, ModelParams& )
   {
      const int level = boundary->getVariables()->level;
      probe->events.push_back( "wall:" + std::to_string( level ) );
      auto& rho = boundary->getVariables()->rho;
      for( int i = 0; i < rho.getSize(); i++ )
         rho.setElement( i, -10.f * ( level + 1 ) );
   }

   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeBoundaryInteraction( FluidPointer&, BoundaryPointer& boundary, ModelParams& )
   {
      probe->events.push_back( "finalize-wall:" + std::to_string( boundary->getVariables()->level ) );
   }

   template< typename BoundaryPointer, typename SourceBoundaryPointer, typename IndexArray >
   void
   updateGhostBoundaryInterpolated(
      BoundaryPointer& ownBoundary,
      SourceBoundaryPointer& srcBoundary,
      const IndexArray& ghostIndices,
      ModelParams& )
   {
      refresh( ownBoundary, srcBoundary, ghostIndices );
   }

   template< typename BoundaryPointer, typename FluidPointer, typename IndexArray >
   void
   updateGhostBoundaryDirectFromSource(
      BoundaryPointer& ownBoundary,
      FluidPointer& srcFluid,
      const IndexArray& ghostIndices,
      ModelParams& )
   {
      refresh( ownBoundary, srcFluid, ghostIndices );
   }

   template< typename FluidPointer, typename BoundaryPointer >
   void
   interaction( FluidPointer&, BoundaryPointer& boundary, ModelParams& )
   {
      const int level = boundary->getVariables()->level;
      probe->events.push_back( "interact:" + std::to_string( level ) );
      probe->ghostAtInteraction[ level ] = boundary->getVariables()->rho.getElement( 1 );
   }

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   interactionWithOpenBoundary( FluidPointer&, OpenBoundaryPointer&, ModelParams& )
   {
      probe->events.push_back( "MR-interaction" );
   }

   template< typename FluidPointer, typename BoundaryPointer >
   void
   finalizeInteraction( FluidPointer&, BoundaryPointer& boundary, ModelParams& )
   {
      probe->events.push_back( "finalize:" + std::to_string( boundary->getVariables()->level ) );
   }

   Probe* probe = nullptr;

private:
   template< typename BoundaryPointer, typename SourcePointer, typename IndexArray >
   void
   refresh( BoundaryPointer& ownBoundary, SourcePointer& source, const IndexArray& ghostIndices )
   {
      const int ownLevel = ownBoundary->getVariables()->level;
      const int sourceLevel = source->getVariables()->level;
      probe->events.push_back(
         "refresh:" + std::to_string( ownLevel ) + "<-" + std::to_string( sourceLevel ) );
      for( int k = 0; k < ghostIndices.getSize(); k++ )
         ownBoundary->getVariables()->rho.setElement(
            ghostIndices.getElement( k ), source->getVariables()->rho.getElement( 0 ) );
   }
};

template< typename Solver >
void
setupSolver( Solver& solver )
{
   solver.numberOfSubsets = 2;
   solver.fluidSets.resize( 2 );
   solver.boundarySets.resize( 2 );
   solver.boundaryGhostUpdate = TNL::SPH::BoundaryGhostUpdate::DirectFromSource;

   for( int level = 0; level < 2; level++ ) {
      auto& fluid = solver.fluidSets[ level ];
      fluid->getParticles()->setSize( 1 );
      fluid->getParticles()->setNumberOfParticles( 1 );
      fluid->getVariables()->setSize( 1 );
      fluid->getVariables()->level = level;
      fluid->getVariables()->rho.setElement( 0, sourceValue( level ) );

      auto& boundary = solver.boundarySets[ level ];
      boundary->getParticles()->setSize( 2 );
      boundary->getParticles()->setNumberOfParticles( 2 );
      boundary->getVariables()->setSize( 2 );
      boundary->getVariables()->level = level;
      boundary->getVariables()->referentialIdx.setElement( 0, 0 );
      boundary->getVariables()->referentialIdx.setElement( 1, 1 );
      boundary->getVariables()->rho.setElement( 0, sourceValue( level ) );
      boundary->getVariables()->rho.setElement( 1, -999.f );
   }

   for( int ownLevel = 0; ownLevel < 2; ownLevel++ ) {
      TNL::SPH::BoundaryGhostParticles< CaseConfig > interface;
      interface.ownIdx = ownLevel;
      interface.neighborIdx = 1 - ownLevel;
      interface.ghostBegin = 1;
      interface.numberOfGhostParticles = 1;
      interface.ghostIndices.setSize( 1 );
      interface.ghostIndices.setElement( 0, 1 );
      solver.boundaryGhostInterfaces.push_back( interface );
   }
}

}  // namespace
