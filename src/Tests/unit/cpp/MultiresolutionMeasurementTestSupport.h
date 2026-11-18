#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TNL/Algorithms/Segments/Ellpack.h>
#include <TNL/Config/ParameterContainer.h>
#include <TNL/Devices/Host.h>
#include <TNL/Particles/CellIndexer.h>
#include <TNL/Particles/ParticlesCellList.h>
#include <TNL/Particles/Writers/VTKWriter.h>
#include <TNL/Pointers/SharedPointer.h>

#include <SPH/shared/thrustExecPolicySelector.h>
#include <SPH/Boundary.h>
#include <SPH/Fluid.h>
#include <SPH/Kernels.h>
#include <SPH/Models/EquationOfState.h>
#include <SPH/Models/WCSPH_BI/IntegrationSchemes/VerletScheme.h>
#include <SPH/Models/WCSPH_BI/Variables.h>
#include <SPH/MultiresolutionRectangleBuffer.h>
#include <SPH/shared/measurement/MeasureTool.h>
#include <SPH/shared/updateInletHydrostaticProfile.h>

namespace MultiresolutionMeasurementTest {

using Device = TNL::Devices::Host;

struct MeasurementTestKernel
{
   __cuda_callable__
   static float
   W( const float, const float h )
   {
      return 1.f + h;
   }
};

template< int Dimension >
class ParticleSystemConfig
{
public:
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;
   static constexpr int spaceDimension = Dimension;
   using UseWithDomainDecomposition = std::false_type;
   using CellIndexerType = TNL::Particles::SimpleCellIndex< Dimension >;
   using NeighborListType = TNL::Algorithms::Segments::Ellpack< Device, int >;
};

template< int Dimension >
class CaseConfig
{
public:
   using DeviceType = Device;
   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;
   static constexpr int spaceDimension = Dimension;
   static constexpr int numberOfPeriodicBuffers = 0;
};

template< int Dimension >
struct State
{
   using SPHConfig = CaseConfig< Dimension >;
};

template< int Dimension >
struct ModelParams
{
   using SPHConfig = CaseConfig< Dimension >;
   using SPHTraitsType = TNL::SPH::SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using KernelFunction = MeasurementTestKernel;
   using EOS = TNL::SPH::EquationsOfState::TaitLinearizedWeaklyCompressibleEOS< SPHConfig >;

   RealType h = 0.5f;
   RealType mass = 4.f;
   RealType dp = 0.25f;
   RealType rho0 = 0.f;
   RealType speedOfSound = 1.f;
   RealType p0 = 0.f;
   VectorType gravity = 0.f;
};

template< int Dimension >
struct Simulation
{
   using SPHConfig = CaseConfig< Dimension >;
   using SPHTraitsType = TNL::SPH::SPHFluidTraits< SPHConfig >;
   using ParticlesType = TNL::Particles::ParticlesCellList< ParticleSystemConfig< Dimension >, Device >;
   using FluidVariables = TNL::SPH::FluidVariables< State< Dimension > >;
   using BoundaryVariables = TNL::SPH::BoundaryVariables< State< Dimension > >;
   using IntegratorVariables = TNL::SPH::IntegrationSchemes::IntegrationSchemeVariables< SPHConfig >;
   using Fluid = TNL::SPH::Fluid< ParticlesType, SPHConfig, FluidVariables, IntegratorVariables >;
   using Boundary = TNL::SPH::Boundary< ParticlesType, SPHConfig, BoundaryVariables, IntegratorVariables >;
   using FluidPointer = TNL::Pointers::SharedPointer< Fluid, Device >;
   using BoundaryPointer = TNL::Pointers::SharedPointer< Boundary, Device >;
   using OpenBoundaryConfig = TNL::SPH::OpenBoundaryConfig< SPHConfig >;
   using MultiresolutionBoundary = TNL::SPH::MultiresolutionBoundary<
      ParticlesType,
      SPHConfig,
      FluidVariables,
      IntegratorVariables,
      OpenBoundaryConfig,
      ModelParams< Dimension > >;
   using MultiresolutionBoundaryPointer = TNL::Pointers::SharedPointer< MultiresolutionBoundary, Device >;
};

inline void
expect( int& failures, bool condition, const std::string& message )
{
   if( ! condition ) {
      ++failures;
      std::cerr << "FAIL: " << message << '\n';
   }
}

inline void
expectNear( int& failures, float actual, float expected, float tolerance, const std::string& message )
{
   if( ! std::isfinite( actual ) || std::abs( actual - expected ) > tolerance ) {
      ++failures;
      std::cerr << "FAIL: " << message << ": expected " << expected << ", got " << actual << '\n';
   }
}

template< int Dimension, typename Values >
void
addVectorParameters( TNL::Config::ParameterContainer& parameters, const std::string& prefix, const Values& values )
{
   static const char* axes[] = { "X", "Y", "Z" };
   for( int d = 0; d < Dimension; ++d )
      parameters.addParameter< typename Values::ValueType >( prefix + axes[ d ], values[ d ] );
}

template< int Dimension, typename Grid >
void
initializeCellGrid(
   Grid& grid,
   const typename Grid::VectorType& origin,
   const typename Grid::CoordinatesType& dimensions = typename Grid::CoordinatesType( 1 ) )
{
   using Point = typename Grid::VectorType;
   TNL::Config::ParameterContainer parameters;
   addVectorParameters< Dimension >( parameters, "GridOrigin", origin );
   addVectorParameters< Dimension >( parameters, "GridSize", dimensions );
   addVectorParameters< Dimension >( parameters, "GridStep", Point( 1.f ) );
   grid.init( parameters, "", Dimension );
}

template< typename SetPointer >
SetPointer
makeSet( const std::vector< typename SetPointer::ObjectType::VectorType >& points,
         const std::vector< float >& densities,
         const std::vector< typename SetPointer::ObjectType::VectorType >& velocities,
         float searchRadius )
{
   using Set = typename SetPointer::ObjectType;
   using Point = typename Set::VectorType;
   using Coordinates = typename Set::CoordinatesType;
   const int count = static_cast< int >( points.size() );
   const int capacity = std::max( 1, count );

   SetPointer set;
   auto particles = set->getParticles();
   particles->setSize( capacity );
   particles->setSearchRadius( searchRadius );
   const Point origin = -16.f;
   Coordinates cellDimensions( 3 );
   for( const auto& point : points )
      for( int d = 0; d < Point::getSize(); ++d ) {
         const int cell = static_cast< int >( std::floor( ( point[ d ] - origin[ d ] ) / searchRadius ) );
         cellDimensions[ d ] = std::max( cellDimensions[ d ], cell + 3 );
      }
   particles->setDimensions( cellDimensions );
   particles->setOrigin( origin );
   particles->setReferentialOrigin( origin );
   particles->setGlobalOriginCoordinates( Coordinates( 0 ) );
   particles->setNumberOfParticles( count );
   set->getVariables()->setSize( capacity );
   set->getIntegratorVariables()->setSize( capacity );
   for( int i = 0; i < count; ++i ) {
      particles->getPoints().setElement( i, points[ i ] );
      set->getVariables()->rho.getArray().setElement( i, densities.empty() ? 1.f : densities[ i ] );
      set->getVariables()->v.getArray().setElement( i, velocities.empty() ? Point( 0.f ) : velocities[ i ] );
   }
   set->searchForNeighbors();
   return set;
}

template< typename SetPointer >
SetPointer
makeSet( const std::vector< typename SetPointer::ObjectType::VectorType >& points,
         const std::vector< float >& densities,
         float searchRadius )
{
   using Point = typename SetPointer::ObjectType::VectorType;
   return makeSet< SetPointer >( points, densities, std::vector< Point >{}, searchRadius );
}

template< int Dimension >
class TestGrid : public TNL::SPH::InterpolateToGrid< CaseConfig< Dimension >, Simulation< Dimension > >
{
public:
   using Base = TNL::SPH::InterpolateToGrid< CaseConfig< Dimension >, Simulation< Dimension > >;
   using Base::interpolate;

   auto&
   outputVariables()
   {
      return *this->variables;
   }

   auto&
   outputGamma()
   {
      return this->gamma;
   }
};

template< int Dimension >
class TestPressureSensor : public TNL::SPH::SensorInterpolation< CaseConfig< Dimension >, Simulation< Dimension > >
{
public:
   using Base = TNL::SPH::SensorInterpolation< CaseConfig< Dimension >, Simulation< Dimension > >;
   using Base::interpolate;

   const typename Base::SensorsDataArray&
   data() const
   {
      return this->sensors;
   }

   int
   historySize() const
   {
      return this->sensorIndexer;
   }
};

template< int Dimension >
class TestWaterSensor : public TNL::SPH::SensorWaterLevel< CaseConfig< Dimension >, Simulation< Dimension > >
{
public:
   using Base = TNL::SPH::SensorWaterLevel< CaseConfig< Dimension >, Simulation< Dimension > >;
   using Base::interpolate;

   int
   historySize() const
   {
      return this->sensorIndexer;
   }
};

inline void
checkLegacyGrid( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 2.f, 0.f } }, { 2.f }, { Point{ 3.f, 5.f } }, 4.f );
   auto boundary =
      makeSet< Simulation2::BoundaryPointer >( { Point{ 2.f, 0.f } }, { 100.f }, { Point{ 100.f, 100.f } }, 4.f );
   TestGrid< 2 > grid;
   initializeCellGrid< 2 >( grid, Point{ 1.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 2, 1 } );
   grid.interpolate( fluid, boundary, params );

   // Hand calculation with base h/mass: W=1+h=1.5, V=m/rho=4/2=2,
   // hence gamma=3, normalized density=2, and normalized velocity=(3,5).
   // searchRadius=4 would instead give rf=4 if grid incorrectly derived MR values.
   expectNear( failures, grid.outputGamma().getElement( 1 ), 3.f, 1.e-6f, "legacy grid uses base h and mass" );
   expectNear( failures, grid.outputVariables().rho.getArray().getElement( 1 ), 2.f, 1.e-6f, "legacy grid density" );
   expectNear(
      failures, grid.outputVariables().v.getArray().getElement( 1 )[ 0 ], 3.f, 1.e-6f, "legacy grid velocity x" );
   expectNear(
      failures, grid.outputVariables().v.getArray().getElement( 1 )[ 1 ], 5.f, 1.e-6f, "legacy grid velocity y" );
}

inline void
checkLegacyPressure( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   params.mass = 8.f;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 4.f }, 2.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 8.f }, 2.f );
   std::vector< Point > positions = { Point{ 1.f, 0.f }, Point{ 20.f, 20.f } };

   // Linear test EOS is p=rho. Direct pressure uses only base h=0.5, mass=8, W=1.5
   // with no refinement-derived scaling. Fluid-only: gamma=3 and density numerator=12
   // give p=4; adding the same-set boundary raises gamma to 4.5 and numerator to 24,
   // preserving p=24/4.5=16/3.
   TestPressureSensor< 2 > fluidOnly;
   fluidOnly.init( positions, 2, 1, false );
   fluidOnly.interpolate( fluid, boundary, params );
   expectNear( failures, fluidOnly.data().getElement( 0, 0 ), 4.f, 1.e-5f, "legacy pressure fluid only" );
   expectNear( failures, fluidOnly.data().getElement( 0, 1 ), 0.f, 1.e-6f, "legacy uncovered pressure is zero" );
   fluidOnly.interpolate( fluid, boundary, params );
   expect( failures, fluidOnly.historySize() == 2, "legacy pressure advances once per direct call" );
   expectNear( failures, fluidOnly.data().getElement( 1, 0 ), 4.f, 1.e-5f, "legacy pressure second history row" );

   TestPressureSensor< 2 > withBoundary;
   withBoundary.init( positions, 2, 1, true );
   withBoundary.interpolate( fluid, boundary, params );
   expectNear( failures, withBoundary.data().getElement( 0, 0 ), 16.f / 3.f, 1.e-5f, "legacy pressure adds boundary" );
   expectNear(
      failures, withBoundary.data().getElement( 0, 1 ), 0.f, 1.e-6f, "legacy boundary leaves uncovered pressure dry" );
}

inline void
initializeWaterSensor( TestWaterSensor< 2 >& sensor,
                       const std::vector< Simulation< 2 >::SPHTraitsType::VectorType >& positions,
                       float endLevel = 2.f )
{
   using Point = Simulation< 2 >::SPHTraitsType::VectorType;
   std::vector< Point > mutablePositions = positions;
   sensor.init( mutablePositions, 2, 1.f, Point{ 0.f, 1.f }, 0.f, endLevel );
}

inline void
checkLegacyWaterLevel( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   const std::vector< Point > positions = { Point{ 0.f, 0.f } };

   ModelParams< 2 > lowGammaParams;
   lowGammaParams.mass = 0.1f;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 0.f, 0.f } }, { 2.f }, 4.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 4.f );
   TestWaterSensor< 2 > lowGamma;
   initializeWaterSensor( lowGamma, positions );
   lowGamma.interpolate( fluid, boundary, lowGammaParams );

   // Legacy base arithmetic gives gamma=(1+0.5)*(0.1/2)=0.075 at both samples,
   // so water level is zero. Derived rf=4 arithmetic would exceed 0.5 and is forbidden here.
   expectNear( failures, lowGamma.getSensorData().getElement( 0, 0 ), 0.f, 1.e-6f, "legacy water uses base h/mass" );
   expect( failures, lowGamma.historySize() == 1, "legacy water advances once per direct call" );

   ModelParams< 2 > params;
   auto distantFluid = makeSet< Simulation2::FluidPointer >( { Point{ 10.f, 10.f } }, { 2.f }, 1.f );
   auto wetBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 0.f, 0.f } }, { 2.f }, 1.f );
   TestWaterSensor< 2 > fluidOnly;
   initializeWaterSensor( fluidOnly, positions );
   fluidOnly.interpolate( distantFluid, wetBoundary, params );
   expectNear( failures, fluidOnly.getSensorData().getElement( 0, 0 ), 0.f, 1.e-6f, "legacy water excludes boundary" );
}

inline void
checkLegacyOpenBoundaryCaller( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   params.rho0 = 1000.f;
   params.speedOfSound = 2.f;
   params.p0 = 5.f;
   params.gravity = Point{ 0.f, -2.f };
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 0.f, 0.f } }, { 2.f }, 1.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 1.f );
   auto openBoundary = makeSet< Simulation2::FluidPointer >( { Point{ 0.f, -1.f } }, { -7.f }, 1.f );

   // Two wet samples at y=0 and y=0.25 give waterLevel=0.5. At open-boundary
   // depth 1.5, hydrostatic p=3000, hydrostatic rho=1748.75, and kinetic
   // density correction is -0.25, giving independently expected rho=1748.5.
   TNL::SPH::features::updateOpenBCHydrostaticProfile< Simulation2 >(
      fluid, boundary, openBoundary, params, Point{ 0.f, 0.f }, Point{ 0.f, 1.f }, 0.25f, 4.f );
   expectNear( failures,
               openBoundary->getVariables()->rho.getArray().getElement( 0 ),
               1748.5f,
               1.e-4f,
               "legacy updateOpenBCHydrostaticProfile direct caller" );
}

inline void
runLegacyMeasurementTests( int& failures )
{
   checkLegacyGrid( failures );
   checkLegacyPressure( failures );
   checkLegacyWaterLevel( failures );
   checkLegacyOpenBoundaryCaller( failures );
}

}  // namespace MultiresolutionMeasurementTest
