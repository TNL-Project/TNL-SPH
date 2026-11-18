#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include <TNL/Devices/Cuda.h>

#include "MultiresolutionBufferCapacityTestSupport.h"

namespace {

using namespace MultiresolutionBufferCapacityTest;

template< typename SetPointer >
void
checkTransferredState( int& failures,
                       const SetPointer& set,
                       const int index,
                       const int stableId,
                       const Snapshot& expected,
                       const std::string& label )
{
   const auto& variables = *set->getVariables();
   const auto& history = *set->getIntegratorVariables();
   expect( failures, variables.referentialIdx.getElement( index ) == stableId, label + ": stable referential ID" );
   expect( failures, set->getPoints().getElement( index ) == expected.point, label + ": position" );
   expect( failures, variables.rho.getArray().getElement( index ) == expected.density, label + ": density" );
   expect( failures, variables.drho.getArray().getElement( index ) == 0.f, label + ": neutral density rate" );
   expect( failures, variables.p.getArray().getElement( index ) == 0.f, label + ": neutral pressure" );
   expect( failures, variables.v.getArray().getElement( index ) == expected.velocity, label + ": velocity" );
   expect( failures, variables.a.getArray().getElement( index ) == Point{ 0.f, 0.f }, label + ": neutral acceleration" );
   expect( failures, variables.gamma.getArray().getElement( index ) == 0.f, label + ": neutral gamma" );
   expect( failures, variables.marker.getArray().getElement( index ) == expected.marker, label + ": marker" );
   expect( failures, history.rho_old.getElement( index ) == expected.density, label + ": density hand-off history" );
   expect( failures, history.v_old.getElement( index ) == expected.velocity, label + ": velocity hand-off history" );
}

template< typename Pointer, bool lts >
void
convertToBuffer( Pointer& buffer, FluidPointer& fluid, ModelParams& params )
{
   if constexpr( lts )
      buffer->convertFluidToBuffer( fluid, params );
   else
      buffer->convertFluidToBuffer( fluid );
}

template< typename Pointer, bool lts >
void
convertToFluid( Pointer& buffer, FluidPointer& fluid, ModelParams& params )
{
   if constexpr( lts )
      buffer->convertBufferToFluid( fluid, params );
   else
      buffer->convertBufferToFluid( fluid );
}

template< typename SetPointer >
void
checkNaNTransferredState( int& failures, const SetPointer& set, const bool sanitized, const std::string& label )
{
   const auto& variables = *set->getVariables();
   const auto& history = *set->getIntegratorVariables();
   const float expectedDensity = sanitized ? 1000.f : 1001.f;
   expect( failures, set->getPoints().getElement( 0 ) == ( Point{ 0.5f, 0.5f } ), label + ": position" );
   expect( failures, variables.referentialIdx.getElement( 0 ) == 1000, label + ": stable referential ID" );
   expect( failures, variables.marker.getArray().getElement( 0 ) == 100, label + ": marker" );
   expect( failures, variables.rho.getArray().getElement( 0 ) == expectedDensity, label + ": density policy" );
   expect( failures, variables.drho.getArray().getElement( 0 ) == 0.f, label + ": neutral density rate" );
   expect( failures, variables.p.getArray().getElement( 0 ) == 0.f, label + ": neutral pressure" );
   expect( failures, variables.a.getArray().getElement( 0 ) == ( Point{ 0.f, 0.f } ), label + ": neutral acceleration" );
   expect( failures, variables.gamma.getArray().getElement( 0 ) == 0.f, label + ": neutral gamma" );
   expect( failures, history.rho_old.getElement( 0 ) == expectedDensity, label + ": density hand-off history" );
   if( sanitized ) {
      expect( failures, variables.v.getArray().getElement( 0 ) == ( Point{ 0.f, 0.f } ), label + ": sanitized velocity" );
      expect( failures, history.v_old.getElement( 0 ) == ( Point{ 0.f, 0.f } ), label + ": sanitized velocity history" );
   }
   else {
      expect( failures, std::isnan( variables.v.getArray().getElement( 0 )[ 0 ] ), label + ": copied NaN velocity" );
      expect( failures, variables.v.getArray().getElement( 0 )[ 1 ] == -1.f, label + ": copied finite velocity component" );
      expect( failures, std::isnan( history.v_old.getElement( 0 )[ 0 ] ), label + ": copied NaN velocity history" );
      expect( failures, history.v_old.getElement( 0 )[ 1 ] == -1.f, label + ": copied finite velocity history component" );
   }
}

template< typename Pointer, bool sanitized >
void
checkNaNVelocityPolicies( int& failures )
{
   const std::string label = sanitized ? "ModelParams" : "plain";
   ModelParams params;
   {
      FluidPointer fluid = makeSet< FluidPointer >( 1, 2 );
      Pointer buffer = makeBuffer< Pointer >( 0, 1 );
      fluid->getVariables()->v.getArray().setElement(
         0, Point{ std::numeric_limits< float >::quiet_NaN(), -1.f } );
      buffer->numberOfParticlesToBuffer = 1;
      buffer->particlesToBuffer.setSize( 1 );
      buffer->particlesToBuffer.setElement( 0, 0 );
      if constexpr( sanitized )
         buffer->convertFluidToBuffer( fluid, params );
      else
         buffer->convertFluidToBuffer( fluid );
      checkNaNTransferredState( failures, buffer, sanitized, label + " capture" );
      expect( failures, fluid->getPoints().getElement( 0 )[ 0 ] == FLT_MAX, label + " capture: source quarantined" );
      expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 1, label + " capture: removal count" );
   }
   {
      FluidPointer fluid = makeSet< FluidPointer >( 0, 0 );
      Pointer buffer = makeBuffer< Pointer >( 1, 2 );
      buffer->getVariables()->v.getArray().setElement(
         0, Point{ std::numeric_limits< float >::quiet_NaN(), -1.f } );
      buffer->numberOfParticlesToRetype = 1;
      if constexpr( sanitized )
         buffer->convertBufferToFluid( fluid, params );
      else
         buffer->convertBufferToFluid( fluid );
      checkNaNTransferredState( failures, fluid, sanitized, label + " return" );
      expect( failures, buffer->getPoints().getElement( 0 )[ 0 ] == FLT_MAX, label + " return: source quarantined" );
      expect( failures, buffer->getNumberOfParticles() == 0, label + " return: active count" );
   }
}

template< typename Pointer, bool lts >
void
checkRoundTripState( int& failures )
{
   const std::string label = lts ? "LTS round trip" : "uniform round trip";
   FluidPointer fluid = makeSet< FluidPointer >( 1, 1 );
   Pointer buffer = makeBuffer< Pointer >( 0, 0 );
   const SnapshotMap source = snapshotOf( fluid );
   constexpr int stableId = 1000;
   ModelParams params;

   buffer->numberOfParticlesToBuffer = 1;
   buffer->particlesToBuffer.setSize( 1 );
   buffer->particlesToBuffer.setElement( 0, 0 );
   convertToBuffer< Pointer, lts >( buffer, fluid, params );

   expect( failures, buffer->getNumberOfAllocatedParticles() == 1, label + ": exact capture capacity" );
   expect( failures, buffer->getNumberOfParticles() == 1, label + ": capture active count" );
   expect( failures, fluid->getPoints().getElement( 0 )[ 0 ] == FLT_MAX, label + ": captured source quarantined" );
   expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 1, label + ": capture removal count" );
   checkTransferredState( failures, buffer, 0, stableId, source.at( stableId ), label + " captured state" );
   expect( failures, buffer->getVariables()->highestReferentialIdx == 4321, label + ": capture next-ID watermark" );

   buffer->numberOfParticlesToRetype = 1;
   buffer->retypeMarker.setElement( 0, 1 );
   buffer->sortBufferParticles();
   convertToFluid< Pointer, lts >( buffer, fluid, params );

   expect( failures, fluid->getNumberOfAllocatedParticles() == 2, label + ": return grows exactly after compaction" );
   expect( failures, fluid->getNumberOfParticles() == 1, label + ": return active count" );
   expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 0, label + ": return clears removal count" );
   checkTransferredState( failures, fluid, 0, stableId, source.at( stableId ), label + " returned state" );
   expect( failures, fluid->getVariables()->highestReferentialIdx == 4321, label + ": return next-ID watermark" );
}

template< typename Pointer, bool lts >
void
checkSortedTransferState( int& failures )
{
   const std::string label = lts ? "LTS sorted transfer" : "uniform sorted transfer";
   FluidPointer fluid = makeSet< FluidPointer >( 0, 0 );
   Pointer buffer = makeBuffer< Pointer >( 3, 3 );
   const SnapshotMap source = snapshotOf( buffer );
   const int transferredId = buffer->getVariables()->referentialIdx.getElement( 2 );
   ModelParams params;

   // Unique keys force the order (old 1, old 2, old 0). After removing key 2, old 2 is retyped.
   buffer->retypeMarker.setElement( 0, 2 );
   buffer->retypeMarker.setElement( 1, 0 );
   buffer->retypeMarker.setElement( 2, 1 );
   buffer->numberOfParticlesToRemove = 1;
   buffer->numberOfParticlesToRetype = 1;
   buffer->sortBufferParticles();
   buffer->removeBufferParticles();
   convertToFluid< Pointer, lts >( buffer, fluid, params );

   expect( failures, fluid->getNumberOfParticles() == 1, label + ": one converted particle" );
   expect( failures, buffer->getNumberOfParticles() == 1, label + ": removal and conversion counts" );
   checkTransferredState( failures, fluid, 0, transferredId, source.at( transferredId ), label );
}

template< typename Pointer >
void
checkSpawnedState( int& failures, const std::string& label )
{
   Pointer buffer = makeBuffer< Pointer >( 0, 0 );
   buffer->massNodes.setSize( 2 );
   buffer->massNodes.points.setElement( 0, Point{ 0.5f, 3.5f } );
   buffer->massNodes.points.setElement( 1, Point{ 1.5f, 3.5f } );
   buffer->massNodes.normal.setElement( 0, Point{ 0.f, 0.f } );
   buffer->massNodes.normal.setElement( 1, Point{ 2.f, 0.f } );
   buffer->numberOfParticlesToCreate = 2;
   ModelParams params;
   buffer->createBufferParticles( params );

   expect( failures, buffer->getNumberOfParticles() == 2, label + ": spawned active count" );
   expect( failures, buffer->getVariables()->highestReferentialIdx == 4323, label + ": spawned IDs advance watermark" );
   for( int i = 0; i < 2; ++i ) {
      const std::string at = label + " slot " + std::to_string( i );
      const auto& variables = *buffer->getVariables();
      const auto& history = *buffer->getIntegratorVariables();
      expect( failures, variables.referentialIdx.getElement( i ) == 4321 + i, at + ": fresh sequential ID" );
      expect( failures, variables.rho.getArray().getElement( i ) == params.rho0, at + ": initial density" );
      expect( failures, variables.v.getArray().getElement( i ) == Point{ 0.f, 0.f }, at + ": initial velocity" );
      expect( failures, variables.drho.getArray().getElement( i ) == 0.f, at + ": neutral density rate" );
      expect( failures, variables.p.getArray().getElement( i ) == 0.f, at + ": neutral pressure" );
      expect( failures, variables.a.getArray().getElement( i ) == Point{ 0.f, 0.f }, at + ": neutral acceleration" );
      expect( failures, variables.gamma.getArray().getElement( i ) == 0.f, at + ": neutral gamma" );
      expect( failures, variables.marker.getArray().getElement( i ) == 0, at + ": ordinary fluid marker" );
      expect( failures, history.rho_old.getElement( i ) == params.rho0, at + ": density hand-off history" );
      expect( failures, history.v_old.getElement( i ) == Point{ 0.f, 0.f }, at + ": velocity hand-off history" );
   }
}

void
checkHostTransfers( int& failures )
{
   checkRoundTripState< BufferPointer, false >( failures );
   checkRoundTripState< BufferLTSPointer, true >( failures );
   checkSortedTransferState< BufferPointer, false >( failures );
   checkSortedTransferState< BufferLTSPointer, true >( failures );
   checkNaNVelocityPolicies< BufferPointer, false >( failures );
   checkNaNVelocityPolicies< BufferLTSPointer, true >( failures );
   checkSpawnedState< BufferPointer >( failures, "uniform spawn" );
   checkSpawnedState< BufferLTSPointer >( failures, "LTS spawn" );
}

namespace CudaGrowth {

using Device = TNL::Devices::Cuda;

class ParticlesConfig
{
public:
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;
   static constexpr int spaceDimension = 2;
   using UseWithDomainDecomposition = std::false_type;
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

struct State
{
   using SPHConfig = CaseConfig;
};

using Traits = TNL::SPH::SPHFluidTraits< CaseConfig >;
using Point = Traits::VectorType;
using Coordinates = Traits::CoordinatesType;
using Particles = TNL::Particles::ParticlesCellList< ParticlesConfig, Device >;
using Variables = TNL::SPH::FluidVariables< State >;
using IntegratorVariables = TNL::SPH::IntegrationSchemes::IntegrationSchemeVariables< CaseConfig >;
using Fluid = TNL::SPH::Fluid< Particles, CaseConfig, Variables, IntegratorVariables >;
using FluidPointer = TNL::Pointers::SharedPointer< Fluid, Device >;

FluidPointer
makeExactlyFullInvalidSet()
{
   FluidPointer fluid;
   auto particles = fluid->getParticles();
   particles->setSize( 5 );
   particles->setSearchRadius( 1.f );
   particles->setDimensions( Coordinates{ 4, 2 } );
   particles->setOrigin( Point{ 0.f, 0.f } );
   particles->setReferentialOrigin( Point{ 0.f, 0.f } );
   particles->setGlobalOriginCoordinates( Coordinates{ 0, 0 } );
   particles->setNumberOfParticles( 5 );
   fluid->getVariables()->setSize( 5 );
   fluid->getIntegratorVariables()->setSize( 5 );

   const Point points[] = {
      { 0.5f, 0.5f },
      { -0.25f, 0.5f },
      { 20.f, 0.5f },
      { std::numeric_limits< float >::quiet_NaN(), 0.5f },
      { std::numeric_limits< float >::infinity(), 0.5f },
   };
   auto& variables = *fluid->getVariables();
   auto& history = *fluid->getIntegratorVariables();
   for( int i = 0; i < 5; ++i ) {
      particles->getPoints().setElement( i, points[ i ] );
      variables.referentialIdx.setElement( i, 700 + i );
      variables.rho.getArray().setElement( i, 1000.f + i );
      variables.drho.getArray().setElement( i, 10.f + i );
      variables.p.getArray().setElement( i, 100.f + i );
      variables.v.getArray().setElement( i, Point{ 1.f + i, -1.f - i } );
      variables.a.getArray().setElement( i, Point{ 2.f + i, 3.f + i } );
      variables.gamma.getArray().setElement( i, 0.5f + i );
      variables.marker.getArray().setElement( i, 100 + i );
      history.rho_old.setElement( i, 900.f + i );
      history.v_old.setElement( i, Point{ -2.f - i, 0.25f + i } );
   }
   variables.highestReferentialIdx = 705;
   return fluid;
}

}  // namespace CudaGrowth

void
checkCudaInvalidGrowth( int& failures )
{
   auto fluid = CudaGrowth::makeExactlyFullInvalidSet();
   expect( failures, fluid->ensureCapacity( 10 ), "CUDA invalid growth: exact capacity grows" );
   expect( failures, fluid->getNumberOfAllocatedParticles() == 10, "CUDA invalid growth: exact requested capacity" );
   expect( failures, fluid->getNumberOfParticles() == 1, "CUDA invalid growth: all invalid positions compacted" );
   expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 0,
           "CUDA invalid growth: pending removals cleared" );
   if( fluid->getNumberOfParticles() != 1 )
      return;

   const auto& variables = *fluid->getVariables();
   const auto& history = *fluid->getIntegratorVariables();
   expect( failures, variables.referentialIdx.getElement( 0 ) == 700, "CUDA invalid growth: valid ID survives" );
   expect( failures, fluid->getPoints().getElement( 0 ) == CudaGrowth::Point{ 0.5f, 0.5f },
           "CUDA invalid growth: valid position survives" );
   expect( failures, variables.rho.getArray().getElement( 0 ) == 1000.f, "CUDA invalid growth: density survives" );
   expect( failures, variables.drho.getArray().getElement( 0 ) == 10.f, "CUDA invalid growth: density rate survives" );
   expect( failures, variables.p.getArray().getElement( 0 ) == 100.f, "CUDA invalid growth: pressure survives" );
   expect( failures, variables.v.getArray().getElement( 0 ) == CudaGrowth::Point{ 1.f, -1.f },
           "CUDA invalid growth: velocity survives" );
   expect( failures, variables.a.getArray().getElement( 0 ) == CudaGrowth::Point{ 2.f, 3.f },
           "CUDA invalid growth: acceleration survives" );
   expect( failures, variables.gamma.getArray().getElement( 0 ) == 0.5f, "CUDA invalid growth: gamma survives" );
   expect( failures, variables.marker.getArray().getElement( 0 ) == 100, "CUDA invalid growth: marker survives" );
   expect( failures, history.rho_old.getElement( 0 ) == 900.f, "CUDA invalid growth: density history survives" );
   expect( failures, history.v_old.getElement( 0 ) == CudaGrowth::Point{ -2.f, 0.25f },
           "CUDA invalid growth: velocity history survives" );
   expect( failures, variables.highestReferentialIdx == 705, "CUDA invalid growth: next-ID watermark survives" );
}

}  // namespace

int
main( int argc, char** argv )
{
   const std::string mode = argc > 1 ? argv[ 1 ] : "";
   if( ! mode.empty() && mode != "--host-only" && mode != "--cuda-only" ) {
      std::cerr << "Usage: MultiresolutionGrowthSafetyTest [--host-only|--cuda-only]\n";
      return 2;
   }

   int failures = 0;
   if( mode != "--cuda-only" )
      checkHostTransfers( failures );
   if( mode != "--host-only" )
      checkCudaInvalidGrowth( failures );
   if( failures != 0 ) {
      std::cerr << failures << " multiresolution growth-safety checks failed\n";
      return 1;
   }
   std::cout << "PASS: invalid growth is quarantined and transferred state is complete and identity-preserving\n";
   return 0;
}
