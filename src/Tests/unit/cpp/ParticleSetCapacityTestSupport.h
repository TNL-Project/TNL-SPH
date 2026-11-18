#pragma once

#include <iostream>
#include <map>
#include <string>
#include <tuple>

#include <TNL/Algorithms/Segments/Ellpack.h>
#include <TNL/Devices/Host.h>
#include <TNL/Particles/CellIndexer.h>
#include <TNL/Particles/ParticlesCellList.h>
#include <TNL/Pointers/SharedPointer.h>

#include <SPH/Fluid.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/VerletScheme.h>
#include <SPH/Models/WCSPH_BI/Variables.h>

namespace ParticleSetCapacityTest {

using Device = TNL::Devices::Host;

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
using FluidVariables = TNL::SPH::FluidVariables< State >;
using BoundaryVariables = TNL::SPH::BoundaryVariables< State >;
using OpenBoundaryVariables = TNL::SPH::OpenBoundaryVariables< State >;
using IntegratorVariables = TNL::SPH::IntegrationSchemes::IntegrationSchemeVariables< CaseConfig >;
using Fluid = TNL::SPH::Fluid< Particles, CaseConfig, FluidVariables, IntegratorVariables >;
using FluidPointer = TNL::Pointers::SharedPointer< Fluid, Device >;

struct Snapshot
{
   Point point;
   Point velocity;
   Point acceleration;
   Point oldVelocity;
   float density;
   float densityRate;
   float pressure;
   float gamma;
   float oldDensity;
   int marker;
};

using SnapshotMap = std::map< int, Snapshot >;

inline void
expect( int& failures, bool condition, const std::string& message )
{
   if( ! condition ) {
      ++failures;
      std::cerr << "FAIL: " << message << '\n';
   }
}

inline FluidPointer
makeFluid( int activeCount = 5, int capacity = 6 )
{
   FluidPointer fluid;
   auto particles = fluid->getParticles();
   particles->setSize( capacity );
   particles->setSearchRadius( 1.f );
   particles->setDimensions( Coordinates{ 4, 2 } );
   particles->setOrigin( Point{ 0.f, 0.f } );
   particles->setGlobalOriginCoordinates( Coordinates{ 0, 0 } );
   particles->setNumberOfParticles( activeCount );
   fluid->getVariables()->setSize( capacity );
   fluid->getIntegratorVariables()->setSize( capacity );

   const Point points[] = { { 3.5f, 0.5f }, { 0.5f, 1.5f }, { 2.5f, 0.5f }, { 1.5f, 1.5f }, { 0.2f, 0.4f } };
   auto& variables = *fluid->getVariables();
   auto& history = *fluid->getIntegratorVariables();
   for( int i = 0; i < activeCount; ++i ) {
      const float x = static_cast< float >( i + 1 );
      particles->getPoints().setElement( i, points[ i ] );
      variables.referentialIdx.setElement( i, 1000 + 37 * i );
      variables.rho.getArray().setElement( i, 1000.f + x );
      variables.drho.getArray().setElement( i, 10.f + x );
      variables.p.getArray().setElement( i, 100.f + x );
      variables.v.getArray().setElement( i, Point{ x, -x } );
      variables.a.getArray().setElement( i, Point{ 2.f * x, 3.f * x } );
      variables.gamma.getArray().setElement( i, 0.5f + x );
      variables.marker.getArray().setElement( i, 100 + i );
      history.rho_old.setElement( i, 900.f + x );
      history.v_old.setElement( i, Point{ -2.f * x, 0.25f * x } );
   }
   variables.highestReferentialIdx = 4321;
   return fluid;
}

inline SnapshotMap
snapshot( const FluidPointer& fluid )
{
   SnapshotMap result;
   const auto& particles = *fluid->getParticles();
   const auto& variables = *fluid->getVariables();
   const auto& history = *fluid->getIntegratorVariables();
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i ) {
      const int id = variables.referentialIdx.getElement( i );
      result[ id ] = { particles.getPoints().getElement( i ),
                       variables.v.getArray().getElement( i ),
                       variables.a.getArray().getElement( i ),
                       history.v_old.getElement( i ),
                       variables.rho.getArray().getElement( i ),
                       variables.drho.getArray().getElement( i ),
                       variables.p.getArray().getElement( i ),
                       variables.gamma.getArray().getElement( i ),
                       history.rho_old.getElement( i ),
                       variables.marker.getArray().getElement( i ) };
   }
   return result;
}

inline void
checkSnapshot( int& failures, const FluidPointer& fluid, const SnapshotMap& expected )
{
   const auto& variables = *fluid->getVariables();
   const auto& history = *fluid->getIntegratorVariables();
   expect( failures, expected.size() == static_cast< std::size_t >( fluid->getNumberOfParticles() ), "active count" );
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i ) {
      const int id = variables.referentialIdx.getElement( i );
      const auto found = expected.find( id );
      expect( failures, found != expected.end(), "referential ID preserved at index " + std::to_string( i ) );
      if( found == expected.end() )
         continue;
      const auto& old = found->second;
      const std::string at = " at stable ID " + std::to_string( id );
      expect( failures, fluid->getPoints().getElement( i ) == old.point, "position" + at );
      expect( failures, variables.rho.getArray().getElement( i ) == old.density, "rho" + at );
      expect( failures, variables.drho.getArray().getElement( i ) == old.densityRate, "drho" + at );
      expect( failures, variables.p.getArray().getElement( i ) == old.pressure, "pressure" + at );
      expect( failures, variables.v.getArray().getElement( i ) == old.velocity, "velocity" + at );
      expect( failures, variables.a.getArray().getElement( i ) == old.acceleration, "acceleration" + at );
      expect( failures, variables.gamma.getArray().getElement( i ) == old.gamma, "gamma" + at );
      expect( failures, variables.marker.getArray().getElement( i ) == old.marker, "marker" + at );
      expect( failures, history.rho_old.getElement( i ) == old.oldDensity, "rho_old" + at );
      expect( failures, history.v_old.getElement( i ) == old.oldVelocity, "v_old" + at );
   }
   expect( failures, variables.highestReferentialIdx == 4321, "highest referential ID preserved" );
}

template< typename Variables >
void
seedLeaf( Variables& variables )
{
   variables.setSize( 3 );
   variables.rho.getArray().setElement( 1, 1002.f );
   variables.marker.getArray().setElement( 1, 77 );
   variables.referentialIdx.setElement( 1, 8080 );
   variables.highestReferentialIdx = 9090;
}

template< typename Variables >
void
checkLeafGrowth( int& failures, Variables& variables )
{
   variables.resize( 8 );
   std::apply(
      [ & ]( auto&... fields )
      {
         ( expect( failures, fields.getArray().getSize() == 8, "every leaf field has required capacity" ), ... );
      },
      variables.allFields() );
   expect( failures, variables.referentialIdx.getSize() == 8, "referential capacity" );
   expect( failures, variables.referentialIdx_swap.getSize() == 8, "referential swap capacity" );
   expect( failures, variables.rho.getArray().getElement( 1 ) == 1002.f, "leaf density preserved" );
   expect( failures, variables.marker.getArray().getElement( 1 ) == 77, "leaf marker preserved" );
   expect( failures, variables.referentialIdx.getElement( 1 ) == 8080, "leaf referential ID not reset" );
   expect( failures, variables.highestReferentialIdx == 9090, "leaf highest ID preserved" );
}

}  // namespace ParticleSetCapacityTest
