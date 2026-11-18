#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>

#include <unistd.h>

#include "MultiresolutionMeasurementTestSupport.h"
#include <TNL/Particles/Readers/VTKReader.h>

#include <SPH/shared/measurement/MeasurementLevelView.h>
#include <SPH/TimeStep.h>
#include <SPH/solvers/SolverMultiSetBlockMultiresolution.h>

namespace {

using namespace MultiresolutionMeasurementTest;

template< int Dimension >
using OwnerBox = TNL::SPH::MeasurementOwnerBox< CaseConfig< Dimension > >;

template< int Dimension >
using LevelView = TNL::SPH::MeasurementLevelView<
   CaseConfig< Dimension >,
   Simulation< Dimension >,
   typename Simulation< Dimension >::MultiresolutionBoundaryPointer >;

template< int Dimension >
OwnerBox< Dimension >
makeOwnerBox( const typename CaseConfig< Dimension >::RealType& minX,
              const typename CaseConfig< Dimension >::RealType& minY,
              const typename CaseConfig< Dimension >::RealType& maxX,
              const typename CaseConfig< Dimension >::RealType& maxY,
              const typename CaseConfig< Dimension >::RealType& refinementFactor )
requires( Dimension == 2 )
{
   using Traits = TNL::SPH::SPHFluidTraits< CaseConfig< Dimension > >;
   using Point = typename Traits::VectorType;
   OwnerBox< Dimension > box;
   box.minimum = Point{ minX, minY };
   box.maximum = Point{ maxX, maxY };
   box.refinementFactor = refinementFactor;
   return box;
}

template< int Dimension >
OwnerBox< Dimension >
makeOwnerBox( const typename Simulation< Dimension >::SPHTraitsType::VectorType& minimum,
              const typename Simulation< Dimension >::SPHTraitsType::VectorType& maximum,
              float refinementFactor )
{
   OwnerBox< Dimension > box;
   box.minimum = minimum;
   box.maximum = maximum;
   box.refinementFactor = refinementFactor;
   return box;
}

template< int Dimension >
LevelView< Dimension >
makeLevel( typename Simulation< Dimension >::FluidPointer& fluid,
           typename Simulation< Dimension >::BoundaryPointer& boundary,
           const std::vector< typename Simulation< Dimension >::MultiresolutionBoundaryPointer >& allPatches,
           const std::vector< InterfaceDescriptor >& interfaces,
           int ownIdx,
           const OwnerBox< Dimension >& ownerBox,
           ModelParams< Dimension >& modelParams )
{
   // Approved compile seam: the view is non-owning, selects all and only the
   // directed ownIdx patches, and derives h/mass from the own-fluid search
   // radius using rf=searchRadius/(2*modelParams.h).
   return LevelView< Dimension >( fluid, boundary, allPatches, interfaces, ownIdx, ownerBox, modelParams );
}

struct EmptyRoutingScheme
{};

struct EmptyRoutingOpenBoundaryModel
{};

struct RoutingParams
{
   using SPHConfig = CaseConfig< 2 >;
   using SPHTraitsType = TNL::SPH::SPHFluidTraits< SPHConfig >;
   using VectorType = typename SPHTraitsType::VectorType;
   using KernelFunction = MeasurementTestKernel;
   using EOS = ModelParams< 2 >::EOS;
   using TimeStepping = TNL::SPH::ConstantTimeStep< SPHConfig >;

   float h = 0.5f;
   float mass = 4.f;
   float dp = 0.25f;
   float rho0 = 0.f;
   float speedOfSound = 1.f;
   float p0 = 0.f;
   VectorType gravity = 0.f;
};

struct RoutingModel
{
   using DeviceType = Device;
   using ModelConfigType = RoutingParams;
   using ModelParams = RoutingParams;
   using ParticlesType = typename Simulation< 2 >::ParticlesType;
   using IntegrationSchemeType = EmptyRoutingScheme;
   using IntegrationSchemeVariables = typename Simulation< 2 >::IntegratorVariables;
   using SPHConfig = CaseConfig< 2 >;
   using FluidVariables = typename Simulation< 2 >::FluidVariables;
   using BoundaryVariables = typename Simulation< 2 >::BoundaryVariables;
   using OpenBoundaryVariables = TNL::SPH::OpenBoundaryVariables< State< 2 > >;
   using OpenBoundaryConfig = typename Simulation< 2 >::OpenBoundaryConfig;
   using OpenBoundaryModel = EmptyRoutingOpenBoundaryModel;

   template< typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, RoutingParams& modelParams )
   {
      auto view_rho = physicalObject->getVariables()->rho.getView();
      auto view_p = physicalObject->getVariables()->p.getView();
      const float rho0 = modelParams.rho0;
      physicalObject->getParticles()->forAll(
         [ = ] __cuda_callable__( int i ) mutable
         {
            view_p[ i ] = view_rho[ i ] - rho0;
         } );
   }
};

using RoutingSolver = TNL::SPH::SolverMultiSetBlockMultiresolution< RoutingModel >;

static_assert(
   std::is_same_v< RoutingSolver::MeasurementLevel,
                   TNL::SPH::MeasurementLevelView<
                      CaseConfig< 2 >,
                      RoutingSolver::BaseType,
                      typename RoutingSolver::MultiresolutionBoundaryPointer > >,
   "solver measurement views use the monitor's simulation type" );
static_assert( std::is_same_v< RoutingSolver::MeasurementLevels, std::vector< RoutingSolver::MeasurementLevel > >,
               "solver exposes a vector of measurement levels" );

template< typename Solver >
concept ExposesMultiresolutionBoundary = requires { typename Solver::MultiresolutionBoundary; };
template< typename Solver >
concept ExposesMultiresolutionBoundaryPointer = requires { typename Solver::MultiresolutionBoundaryPointer; };
using RoutingBase = RoutingSolver::BaseType;
static_assert( ! ExposesMultiresolutionBoundary< RoutingBase > );
static_assert( ! ExposesMultiresolutionBoundaryPointer< RoutingBase > );
static_assert( ExposesMultiresolutionBoundary< RoutingSolver > );
static_assert( ExposesMultiresolutionBoundaryPointer< RoutingSolver > );

void
initializeRoutingSolver( RoutingSolver& solver )
{
   using Point = Simulation< 2 >::SPHTraitsType::VectorType;
   using Coordinates = Simulation< 2 >::SPHTraitsType::CoordinatesType;
   TNL::Config::ParameterContainer parameters;
   TNL::Config::ParameterContainer subdomains;
   parameters.addParameter< float >( "searchRadius", 1.f );
   addVectorParameters< 2 >( parameters, "domainOrigin", Point{ -2.f, -2.f } );
   addVectorParameters< 2 >( parameters, "domainSize", Point{ 12.f, 10.f } );
   parameters.addParameter< int >( "numberOfSubdomains", 2 );

   subdomains.addParameter< float >( "subdomain0RefinementFactor", 1.f );
   addVectorParameters< 2 >( subdomains, "subdomain0GridDimensions", Coordinates{ 4, 4 } );
   addVectorParameters< 2 >( subdomains, "subdomain0Origin", Point{ 0.f, -1.f } );
   addVectorParameters< 2 >( subdomains, "subdomain0OriginGlobalCoords", Coordinates{ 0, 0 } );
   subdomains.addParameter< float >( "subdomain1RefinementFactor", 0.5f );
   addVectorParameters< 2 >( subdomains, "subdomain1GridDimensions", Coordinates{ 4, 4 } );
   addVectorParameters< 2 >( subdomains, "subdomain1Origin", Point{ 1.f, -1.f } );
   addVectorParameters< 2 >( subdomains, "subdomain1OriginGlobalCoords", Coordinates{ 1, 0 } );
   solver.topology.loadFromConfig( parameters, subdomains );

   solver.numberOfSubsets = 2;
   solver.fluidSets.resize( 2 );
   solver.boundarySets.resize( 2 );
   solver.fluidSets[ 0 ] = makeSet< RoutingSolver::FluidPointer >( { Point{ 1.5f, 0.f } }, { 1.f }, 1.f );
   solver.fluidSets[ 1 ] = makeSet< RoutingSolver::FluidPointer >( { Point{ 1.5f, 0.f } }, { 0.25f }, 0.5f );
   solver.boundarySets[ 0 ] = makeSet< RoutingSolver::BoundaryPointer >( { Point{ 12.f, 12.f } }, { 1.f }, 1.f );
   solver.boundarySets[ 1 ] = makeSet< RoutingSolver::BoundaryPointer >( { Point{ 12.f, 12.f } }, { 1.f }, 0.5f );

   using PatchPointer = RoutingSolver::MultiresolutionBoundaryPointer;
   solver.multiresolutionBoundaryPatches.resize( 5 );
   for( int p = 0; p < 5; ++p )
      solver.multiresolutionBoundaryPatches[ p ] =
         makeSet< PatchPointer >( { Point{ 1.5f, 0.f } }, { 0.25f }, p == 2 ? 0.5f : 1.f );
   // Only ownIdx selects sources. The neighbor indices are deliberately absurd
   // sentinels to expose any accidental reverse/other-level interpretation.
   solver.multiresolutionBoundaryPatchInterfaces = { { 0, 99 }, { 0, -99 }, { 1, 77 }, { 0, 77 }, { 2, 0 } };
}

void
checkSolverMeasurementBuilder( int& failures )
{
   using Point = Simulation< 2 >::SPHTraitsType::VectorType;
   RoutingSolver solver;
   initializeRoutingSolver( solver );
   const float baseH = solver.modelParams.h;
   const float baseMass = solver.modelParams.mass;
   const auto levels = solver.buildMeasurementLevels();

   expect( failures, levels.size() == 2, "actual solver builds one view per topology level" );
   if( levels.size() != 2 )
      return;
   expect( failures, levels[ 0 ].levelId == 0 && levels[ 1 ].levelId == 1, "builder preserves level ids" );
   expect( failures, levels[ 0 ].fluid == &solver.fluidSets[ 0 ] && levels[ 1 ].fluid == &solver.fluidSets[ 1 ],
           "builder selects each level's own fluid" );
   expect( failures, levels[ 0 ].boundary == &solver.boundarySets[ 0 ] && levels[ 1 ].boundary == &solver.boundarySets[ 1 ],
           "builder selects each level's own boundary" );
   expect( failures, levels[ 0 ].ownerBox.minimum == ( Point{ 0.f, -1.f } ), "coarse nominal owner minimum" );
   expect( failures, levels[ 0 ].ownerBox.maximum == ( Point{ 4.f, 3.f } ), "coarse nominal owner maximum from grid" );
   expect( failures, levels[ 1 ].ownerBox.minimum == ( Point{ 1.f, -1.f } ), "fine nominal owner minimum" );
   expect( failures, levels[ 1 ].ownerBox.maximum == ( Point{ 3.f, 1.f } ), "fine nominal owner maximum from scaled grid" );
   expectNear( failures, levels[ 0 ].ownerBox.refinementFactor, 1.f, 1.e-7f, "coarse owner refinement factor" );
   expectNear( failures, levels[ 1 ].ownerBox.refinementFactor, 0.5f, 1.e-7f, "fine owner refinement factor" );
   expectNear( failures, levels[ 0 ].h, 0.5f, 1.e-7f, "coarse builder h" );
   expectNear( failures, levels[ 0 ].mass, 4.f, 1.e-7f, "coarse builder mass" );
   expectNear( failures, levels[ 1 ].h, 0.25f, 1.e-7f, "fine builder h" );
   expectNear( failures, levels[ 1 ].mass, 1.f, 1.e-7f, "fine builder mass scales by rf squared" );

   expect( failures, levels[ 0 ].ownPatches.size() == 3, "builder includes all and only level-zero ownIdx patches" );
   expect( failures, levels[ 1 ].ownPatches.size() == 1, "builder includes all and only level-one ownIdx patches" );
   const std::size_t expectedCoarse[] = { 0, 1, 3 };
   for( std::size_t i = 0; i < 3; ++i )
      expect( failures,
              levels[ 0 ].ownPatches[ i ] == &solver.multiresolutionBoundaryPatches[ expectedCoarse[ i ] ],
              "builder preserves ownIdx patch vector order " + std::to_string( i ) );
   expect( failures,
           levels[ 1 ].ownPatches[ 0 ] == &solver.multiresolutionBoundaryPatches[ 2 ],
           "builder selects the exact level-one patch" );
   expectNear( failures, solver.modelParams.h, baseH, 0.f, "builder does not mutate base h" );
   expectNear( failures, solver.modelParams.mass, baseMass, 0.f, "builder does not mutate base mass" );

   const auto fresh = solver.buildMeasurementLevels();
   expect( failures, fresh.size() == 2 && fresh[ 1 ].ownPatches[ 0 ] == levels[ 1 ].ownPatches[ 0 ],
           "builder reconstructs fresh views from canonical patch storage" );
   expect( failures, TNL::SPH::validateMeasurementLevels( fresh ), "builder produces valid nested owner boxes" );
}

std::filesystem::path
configureRoutingMonitor( RoutingSolver& solver )
{
   static int monitorSequence = 0;
   const std::filesystem::path outputDirectory = std::filesystem::path( "/tmp/opencode" )
      / ( "mr-measurement-routing-" + std::to_string( getpid() ) + "-" + std::to_string( monitorSequence++ ) );
   std::filesystem::remove_all( outputDirectory );
   std::filesystem::create_directories( outputDirectory );
   solver.timeStepping.setTimeStep( 1.f );
   solver.timeStepping.setEndTime( 2.f );

   const std::string measurement = R"json({
      "plane1Identifier": "interpolation",
      "plane1OutputTime": 1.0,
      "plane1GridOriginX": 1.0,
      "plane1GridOriginY": 0.0,
      "plane1GridSizeX": 2,
      "plane1GridSizeY": 1,
      "plane1GridStepX": 1.0,
      "plane1GridStepY": 1.0,
      "sensorPOutputTime": 1.0,
      "sensorPIncludeBoundary": false,
      "sensorPPoint1X": 1.5,
      "sensorPPoint1Y": 0.0,
      "sensorWlOutputTime": 1.0,
      "sensorWlStartLevel": 0.0,
      "sensorWlEndLevel": 2.0,
      "sensorWlLevelIncrement": 1.0,
      "sensorWlDirectionX": 0.0,
      "sensorWlDirectionY": 1.0,
      "sensorWlPoint1X": 1.5,
      "sensorWlPoint1Y": 0.0
   })json";
   TNL::Config::ParameterContainer parameters;
   parameters.addParameter< std::string >( "outputDirectory", outputDirectory.string() );
   parameters.addParameter< float >( "finalTime", 2.f );
   parameters.addParameter< float >( "snapshotPeriod", 1.f );
   parameters.addParameter< int >( "interpolationPlanesCount", 1 );
   parameters.addParameter< int >( "pressureSensorsCount", 1 );
   parameters.addParameter< int >( "waterLevelSensorsCount", 1 );
   parameters.addParameter< int >( "volumetricFlowRatePlanesCount", 0 );
   parameters.addParameter< std::string >( "measuretoolConfig", "" );
   parameters.addParameter< std::string >( "measuretool", measurement );
   solver.simulationMonitor.init( parameters, solver.timeStepping, solver.logger );
   return outputDirectory;
}

std::vector< float >
readSensorValues( int& failures, const std::filesystem::path& path, std::vector< int >& indices )
{
   std::vector< float > values;
   std::ifstream input( path );
   expect( failures, input.good(), "sensor history exists: " + path.string() );
   if( ! input.good() )
      return values;
   std::string line;
   while( std::getline( input, line ) ) {
      if( line.empty() || line[ 0 ] == '#' )
         continue;
      std::istringstream row( line );
      int index = -1;
      float value = 0.f;
      if( !( row >> index >> value ) || ! std::isfinite( value ) ) {
         expect( failures, false, "finite sensor history row in " + path.string() );
         continue;
      }
      indices.push_back( index );
      values.push_back( value );
   }
   return values;
}

void
checkSolverMeasurementDispatch( int& failures )
{
   RoutingSolver solver;
   initializeRoutingSolver( solver );
   const std::filesystem::path outputDirectory = configureRoutingMonitor( solver );
   solver.measure();
   solver.timeStepping.updateTimeStep();
   solver.measure();
   solver.simulationMonitor.save( solver.logger );

   std::vector< int > pressureIndices;
   std::vector< int > waterIndices;
   const auto pressure = readSensorValues( failures, outputDirectory / "sensorsPressure.dat", pressureIndices );
   const auto water = readSensorValues( failures, outputDirectory / "sensorsWaterLevel.dat", waterIndices );
   expect( failures, pressureIndices == std::vector< int >( { 0, 1 } ), "multi-view pressure has one row per monitor event" );
   expect( failures, waterIndices == std::vector< int >( { 0, 1 } ), "multi-view water has one row per monitor event" );
   for( const float value : pressure )
      expectNear( failures, value, 0.25f, 1.e-5f, "actual solver pressure routes through measured levels" );
   for( const float value : water )
      expectNear( failures, value, 1.f, 1.e-5f, "actual solver water routes through measured levels" );

   int gridFiles = 0;
   for( const auto& entry : std::filesystem::directory_iterator( outputDirectory ) )
      if( entry.path().extension() == ".vtk" && entry.path().filename().string().find( "interpol" ) != std::string::npos )
         ++gridFiles;
   expect( failures, gridFiles == 2, "monitor saves one grid per plane and event, not per level" );
   std::filesystem::remove_all( outputDirectory );
}

void
checkHalfOpenOwnership( int& failures )
{
   using Config2 = CaseConfig< 2 >;
   using Point = Simulation< 2 >::SPHTraitsType::VectorType;
   using OwnerArray = TNL::Containers::Array< OwnerBox< 2 >, Device >;
   using PointArray = TNL::Containers::Array< Point, Device >;
   using IndexArray = TNL::Containers::Array< int, Device >;

   OwnerArray boxes;
   boxes.setSize( 1 );
   boxes.setElement( 0, makeOwnerBox< 2 >( 0.f, -1.f, 2.f, 1.f, 1.f ) );
   const Point queries[] = {
      Point{ 0.f, 0.f },  // min x is excluded
      Point{ 2.f, 1.f },  // max x and max y are included
      Point{ 0.f, 1.f },  // min x remains excluded on the max-y seam
      Point{ 2.f, -1.f }, // min y is excluded on the max-x seam
      Point{ 1.f, 0.f },
      Point{ 3.f, 0.f },
   };
   const int expected[] = { -1, 0, -1, -1, 0, -1 };
   PointArray points;
   IndexArray owners;
   points.setSize( 6 );
   owners.setSize( 6 );
   for( int i = 0; i < 6; ++i )
      points.setElement( i, queries[ i ] );
   owners = -77;

   TNL::SPH::computeMeasurementOwners( points, boxes, owners );
   for( int i = 0; i < 6; ++i ) {
      expect( failures, owners.getElement( i ) == expected[ i ], "(min,max] per-axis ownership seam " + std::to_string( i ) );
      expect( failures,
              TNL::SPH::findOwningMeasurementLevel( queries[ i ], boxes.getConstView() ) == expected[ i ],
              "single-query ownership agrees with batched ownership " + std::to_string( i ) );
   }
}

void
checkNestedDisjointAndReversedOwnership( int& failures )
{
   using Config2 = CaseConfig< 2 >;
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   using FluidPointer = Simulation2::FluidPointer;
   using BoundaryPointer = Simulation2::BoundaryPointer;
   using PatchPointer = Simulation2::MultiresolutionBoundaryPointer;

   auto fluidA = makeSet< FluidPointer >( { Point{ 1.5f, 0.f } }, { 2.f }, 1.f );
   auto fineFluid = makeSet< FluidPointer >( { Point{ 1.5f, 0.f } }, { 2.f }, 0.5f );
   auto fluidB = makeSet< FluidPointer >( { Point{ 5.5f, 0.f } }, { 2.f }, 1.f );
   auto boundaryA = makeSet< BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 1.f );
   auto fineBoundary = makeSet< BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 0.5f );
   auto boundaryB = makeSet< BoundaryPointer >( { Point{ 14.f, 14.f } }, { 2.f }, 1.f );
   const std::vector< PatchPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   ModelParams< 2 > params;

   LevelView< 2 > coarse = makeLevel< 2 >(
      fluidA, boundaryA, noPatches, noInterfaces, 0, makeOwnerBox< 2 >( 0.f, -1.f, 4.f, 1.f, 1.f ), params );
   LevelView< 2 > fine = makeLevel< 2 >(
      fineFluid, fineBoundary, noPatches, noInterfaces, 1, makeOwnerBox< 2 >( 1.f, -1.f, 2.f, 1.f, 0.5f ), params );
   std::vector< LevelView< 2 > > ordered = { coarse, fine };
   std::vector< LevelView< 2 > > reversed = { fine, coarse };

   const auto orderedBoxes = TNL::SPH::makeMeasurementOwnerBoxArray< Config2 >( ordered );
   const auto reversedBoxes = TNL::SPH::makeMeasurementOwnerBoxArray< Config2 >( reversed );
   expect( failures,
           TNL::SPH::findOwningMeasurementLevel( Point{ 1.5f, 0.f }, orderedBoxes.getConstView() ) == 1,
           "nested query selects finest level in coarse-to-fine storage" );
   expect( failures,
           TNL::SPH::findOwningMeasurementLevel( Point{ 1.5f, 0.f }, reversedBoxes.getConstView() ) == 0,
           "nested query selects finest level independent of vector order" );
   expect( failures,
           TNL::SPH::findOwningMeasurementLevel( Point{ 3.f, 0.f }, orderedBoxes.getConstView() ) == 0,
           "nested coarse level owns its uncovered region" );
   expect( failures,
           TNL::SPH::findOwningMeasurementLevel( Point{ 9.f, 0.f }, orderedBoxes.getConstView() ) == -1,
           "uncovered nested query has signed minus-one owner" );
   expect( failures, TNL::SPH::validateMeasurementLevels( ordered ), "different-rf nested boxes are valid" );

   LevelView< 2 > disjoint = makeLevel< 2 >(
      fluidB, boundaryB, noPatches, noInterfaces, 2, makeOwnerBox< 2 >( 5.f, -1.f, 6.f, 1.f, 1.f ), params );
   std::vector< LevelView< 2 > > disjointLevels = { coarse, disjoint };
   const auto disjointBoxes = TNL::SPH::makeMeasurementOwnerBoxArray< Config2 >( disjointLevels );
   expect( failures,
           TNL::SPH::findOwningMeasurementLevel( Point{ 5.5f, 0.f }, disjointBoxes.getConstView() ) == 1,
           "disjoint boxes retain their independent owners" );
   expect( failures, TNL::SPH::validateMeasurementLevels( disjointLevels ), "disjoint same-rf boxes are valid" );
}

void
checkOverlapValidation( int& failures )
{
   using Point = Simulation< 2 >::SPHTraitsType::VectorType;
   using FluidPointer = Simulation< 2 >::FluidPointer;
   using BoundaryPointer = Simulation< 2 >::BoundaryPointer;
   using PatchPointer = Simulation< 2 >::MultiresolutionBoundaryPointer;
   auto fluid = makeSet< FluidPointer >( { Point{ 0.5f, 0.f } }, { 2.f }, 1.f );
   auto boundary = makeSet< BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 1.f );
   const std::vector< PatchPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   ModelParams< 2 > params;

   int nextLevelId = 0;
   auto makeDummyLevel = [ & ]( const OwnerBox< 2 >& box )
   {
      return makeLevel< 2 >( fluid, boundary, noPatches, noInterfaces, nextLevelId++, box, params );
   };
   const auto left = makeOwnerBox< 2 >( 0.f, -1.f, 1.f, 1.f, 1.f );
   const auto touching = makeOwnerBox< 2 >( 1.f, -1.f, 2.f, 1.f, 1.f );
   const auto overlapping = makeOwnerBox< 2 >( 0.5f, -1.f, 2.f, 1.f, 1.f );
   expect( failures,
           TNL::SPH::validateMeasurementLevels( std::vector< LevelView< 2 > >{ makeDummyLevel( left ),
                                                                                makeDummyLevel( touching ) } ),
           "same-rf boxes sharing only a face are valid" );
   expect( failures,
            ! TNL::SPH::validateMeasurementLevels( std::vector< LevelView< 2 > >{ makeDummyLevel( left ),
                                                                                   makeDummyLevel( overlapping ) } ),
            "positive-volume same-rf overlap is rejected" );
   expect( failures,
            ! TNL::SPH::validateMeasurementLevels( std::vector< LevelView< 2 > >{
               makeLevel< 2 >( fluid, boundary, noPatches, noInterfaces, 7, left, params ),
               makeLevel< 2 >( fluid, boundary, noPatches, noInterfaces, 7, touching, params ) } ),
            "duplicate level ids are rejected even for disjoint boxes" );
}

void
checkViewScalingAndPatchSelection( int& failures )
{
   using Point2 = Simulation< 2 >::SPHTraitsType::VectorType;
   using Simulation2 = Simulation< 2 >;
   auto fluid2 = makeSet< Simulation2::FluidPointer >( { Point2{ 1.f, 0.f } }, { 2.f }, 2.f );
   auto boundary2 = makeSet< Simulation2::BoundaryPointer >( { Point2{ 12.f, 12.f } }, { 2.f }, 2.f );
   auto patch0 = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point2{ 1.f, 0.f } }, { 4.f }, 2.f );
   auto patch1 = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point2{ 1.f, 0.f } }, { 8.f }, 2.f );
   auto reverse = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point2{ 1.f, 0.f } }, { 16.f }, 1.f );
   std::vector< Simulation2::MultiresolutionBoundaryPointer > patches = { patch0, patch1, reverse };
   const std::vector< InterfaceDescriptor > interfaces = { { 0, 1 }, { 0, 1 }, { 1, 0 } };
   ModelParams< 2 > params2;
   auto level = makeLevel< 2 >(
      fluid2, boundary2, patches, interfaces, 0, makeOwnerBox< 2 >( 0.f, -1.f, 4.f, 1.f, 2.f ), params2 );

   // rf=2/(2*0.5)=2 in 2D: h=2*0.5=1 and mass=2^2*4=16.
   expectNear( failures, level.h, 1.f, 1.e-6f, "2D view derives h from rf" );
   expectNear( failures, level.mass, 16.f, 1.e-6f, "2D view derives mass from rf^2" );
   expect( failures, level.ownPatches.size() == 2, "view includes every directed ownIdx patch" );
   expect( failures,
           level.ownPatches.size() == 2 && level.ownPatches[ 0 ] == &patches[ 0 ] && level.ownPatches[ 1 ] == &patches[ 1 ],
           "view preserves ownIdx patch order and excludes reverse direction" );
   expect( failures, level.fluid == &fluid2, "view non-owning fluid reference addresses the supplied set" );
   expect( failures, level.boundary == &boundary2, "view non-owning boundary reference addresses the supplied set" );

   using Point3 = Simulation< 3 >::SPHTraitsType::VectorType;
   using Simulation3 = Simulation< 3 >;
   auto fluid3 = makeSet< Simulation3::FluidPointer >( { Point3{ 0.f, 0.f, 0.f } }, { 2.f }, 0.25f );
   auto boundary3 = makeSet< Simulation3::BoundaryPointer >( { Point3{ 12.f, 12.f, 12.f } }, { 2.f }, 0.25f );
   const std::vector< Simulation3::MultiresolutionBoundaryPointer > noPatches3;
   const std::vector< InterfaceDescriptor > noInterfaces3;
   ModelParams< 3 > params3;
   params3.mass = 8.f;
   OwnerBox< 3 > box3;
   box3.minimum = Point3{ -1.f, -1.f, -1.f };
   box3.maximum = Point3{ 1.f, 1.f, 1.f };
   box3.refinementFactor = 0.25f;
   auto level3 = makeLevel< 3 >( fluid3, boundary3, noPatches3, noInterfaces3, 0, box3, params3 );

   // rf=0.25/(2*0.5)=0.25 in 3D: h=0.125 and mass=0.25^3*8=0.125.
   expectNear( failures, level3.h, 0.125f, 1.e-6f, "3D view derives h from rf" );
   expectNear( failures, level3.mass, 0.125f, 1.e-7f, "3D view derives mass from rf^3" );
}

struct GridViewFixture
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;

   ModelParams< 2 > params;
   Simulation2::FluidPointer coarseFluid =
      makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 2.f }, { Point{ 3.f, 1.f } }, 2.f );
   Simulation2::FluidPointer fineFluid =
      makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 100.f }, { Point{ 100.f, 100.f } }, 1.f );
   Simulation2::BoundaryPointer coarseBoundary =
      makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 200.f }, { Point{ 200.f, 200.f } }, 2.f );
   Simulation2::BoundaryPointer fineBoundary =
      makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 300.f }, { Point{ 300.f, 300.f } }, 1.f );
   Simulation2::MultiresolutionBoundaryPointer patchA = makeSet< Simulation2::MultiresolutionBoundaryPointer >(
      { Point{ 1.f, 0.f } }, { 4.f }, { Point{ 5.f, 2.f } }, 2.f );
   Simulation2::MultiresolutionBoundaryPointer patchB = makeSet< Simulation2::MultiresolutionBoundaryPointer >(
      { Point{ 1.f, 0.f } }, { 8.f }, { Point{ 11.f, 4.f } }, 2.f );
   Simulation2::MultiresolutionBoundaryPointer reverse = makeSet< Simulation2::MultiresolutionBoundaryPointer >(
      { Point{ 1.f, 0.f } }, { 1000.f }, { Point{ 1000.f, 1000.f } }, 1.f );
   std::vector< Simulation2::MultiresolutionBoundaryPointer > patches = { patchA, patchB, reverse };
   std::vector< InterfaceDescriptor > interfaces = { { 0, 1 }, { 0, 1 }, { 1, 0 } };
   OwnerBox< 2 > coarseBox = makeOwnerBox< 2 >( 0.f, -1.f, 4.f, 1.f, 2.f );
   OwnerBox< 2 > fineBox = makeOwnerBox< 2 >( 1.f, -1.f, 2.f, 1.f, 1.f );

   LevelView< 2 >
   coarseLevel()
   {
      return makeLevel< 2 >( coarseFluid, coarseBoundary, patches, interfaces, 0, coarseBox, params );
   }

   LevelView< 2 >
   fineLevel()
   {
      return makeLevel< 2 >( fineFluid, fineBoundary, patches, interfaces, 1, fineBox, params );
   }
};

void
checkGridNormalizationAndSources( int& failures )
{
   GridViewFixture fixture;
   std::vector< LevelView< 2 > > levels = { fixture.coarseLevel(), fixture.fineLevel() };
   TestGrid< 2 > grid;
   // The asserted cell has linear index 3, avoiding the legacy same-set i==j
   // index collision without changing the independent source expectation.
   initializeCellGrid< 2 >( grid, GridViewFixture::Point{ -2.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 4, 1 } );
   grid.interpolate( levels, fixture.params );

   // Query (1,0) is on fine.min x and therefore belongs to the coarse level.
   // Coarse h=1,m=16,W=2. Own fluid contributes gamma=16, own patches 8 and 4.
   // Thus gamma=28, density numerator=3*32=96, velocity x numerator=132 and
   // velocity y numerator=48. Fine/reverse/boundary sources must not contribute.
   expectNear( failures,
               grid.outputGamma().getElement( 3 ),
               28.f,
               1.e-5f,
               "grid accumulates own fluid and all ownIdx patches" );
   expectNear( failures,
               grid.outputVariables().rho.getArray().getElement( 3 ),
               24.f / 7.f,
               1.e-5f,
               "grid hand-normalized density excludes other-level sources" );
   expectNear( failures,
               grid.outputVariables().v.getArray().getElement( 3 )[ 0 ],
               33.f / 7.f,
               1.e-5f,
               "grid hand-normalized velocity x excludes boundary" );
   expectNear( failures,
               grid.outputVariables().v.getArray().getElement( 3 )[ 1 ],
               12.f / 7.f,
               1.e-5f,
               "grid hand-normalized velocity y" );
}

struct DryFineFixture
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;

   ModelParams< 2 > params;
   Simulation2::FluidPointer coarseFluid =
      makeSet< Simulation2::FluidPointer >( { Point{ 0.5f, 0.f } }, { 2.f }, { Point{ 3.f, 4.f } }, 2.f );
   Simulation2::FluidPointer fineFluid = makeSet< Simulation2::FluidPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 1.f );
   Simulation2::BoundaryPointer coarseBoundary =
      makeSet< Simulation2::BoundaryPointer >( { Point{ 0.5f, 0.f } }, { 2.f }, 2.f );
   Simulation2::BoundaryPointer fineBoundary =
      makeSet< Simulation2::BoundaryPointer >( { Point{ 0.5f, 0.f } }, { 2.f }, 1.f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > patches;
   const std::vector< InterfaceDescriptor > interfaces;
   OwnerBox< 2 > coarseBox = makeOwnerBox< 2 >( 0.f, -1.f, 2.f, 1.f, 2.f );
   OwnerBox< 2 > fineBox = makeOwnerBox< 2 >( 0.f, -1.f, 1.f, 1.f, 1.f );

   std::vector< LevelView< 2 > >
   levels()
   {
      return { makeLevel< 2 >( coarseFluid, coarseBoundary, patches, interfaces, 0, coarseBox, params ),
               makeLevel< 2 >( fineFluid, fineBoundary, patches, interfaces, 1, fineBox, params ) };
   }
};

void
checkDryFineAndUncoveredGrid( int& failures )
{
   DryFineFixture fixture;
   auto levels = fixture.levels();
   TestGrid< 2 > dryGrid;
   initializeCellGrid< 2 >( dryGrid, DryFineFixture::Point{ -0.5f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 2, 1 } );
   dryGrid.interpolate( levels, fixture.params );
   expectNear( failures, dryGrid.outputGamma().getElement( 1 ), 0.f, 1.e-6f, "dry fine owner has no coarse fallback gamma" );
   expectNear( failures,
               dryGrid.outputVariables().rho.getArray().getElement( 1 ),
               0.f,
               1.e-6f,
               "dry fine owner has no coarse fallback density" );
   expectNear( failures,
               dryGrid.outputVariables().v.getArray().getElement( 1 )[ 0 ],
               0.f,
               1.e-6f,
               "dry fine owner has no coarse fallback velocity" );

   TestGrid< 2 > uncoveredGrid;
   // Query (-1,0) is outside both boxes despite the co-located coarse source in support.
   initializeCellGrid< 2 >(
      uncoveredGrid, DryFineFixture::Point{ -2.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 2, 1 } );
   uncoveredGrid.interpolate( levels, fixture.params );
   expectNear( failures, uncoveredGrid.outputGamma().getElement( 1 ), 0.f, 1.e-6f, "uncovered grid gamma is zero" );
   expectNear( failures,
               uncoveredGrid.outputVariables().rho.getArray().getElement( 1 ),
               0.f,
               1.e-6f,
               "uncovered grid density is zero" );
   expectNear( failures,
               uncoveredGrid.outputVariables().v.getArray().getElement( 1 )[ 1 ],
               0.f,
               1.e-6f,
               "uncovered grid velocity is zero" );
}

void
checkGammaThresholdWithoutFallback( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   struct Case
   {
      float fineMass;
      float expectedGamma;
      float expectedDensity;
   };
   const Case cases[] = {
      { 0.99f, 0.495f, 0.f }, // below 0.5
      { 1.f, 0.5f, 0.f },     // exactly 0.5 remains dry
      { 1.01f, 0.505f, 3.f }, // immediately above 0.5 normalizes locally
   };
   for( const auto& testCase : cases ) {
      ModelParams< 2 > params;
      params.mass = testCase.fineMass;
      auto coarseFluid = makeSet< Simulation2::FluidPointer >( { Point{ 0.5f, 0.f } }, { 100.f }, 2.f );
      auto fineFluid =
         makeSet< Simulation2::FluidPointer >( { Point{ 0.5f, 0.f } }, { 3.f }, { Point{ 7.f, 0.f } }, 1.f );
      auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 12.f, 12.f } }, { 2.f }, 1.f );
      const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
      const std::vector< InterfaceDescriptor > noInterfaces;
      std::vector< LevelView< 2 > > levels = {
         makeLevel< 2 >(
            coarseFluid, boundary, noPatches, noInterfaces, 0, makeOwnerBox< 2 >( 0.f, -1.f, 2.f, 1.f, 2.f ), params ),
         makeLevel< 2 >(
            fineFluid,
            boundary,
            noPatches,
            noInterfaces,
            1,
            makeOwnerBox< 2 >( 0.f, -1.f, 1.f, 1.f, 1.f ),
            params )
      };
      TestGrid< 2 > grid;
      initializeCellGrid< 2 >( grid, Point{ -0.5f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 2, 1 } );
      grid.interpolate( levels, params );

      // Fine W=1.5 and mass/rho=fineMass/3, so gamma=0.5*fineMass.
      expectNear( failures, grid.outputGamma().getElement( 1 ), testCase.expectedGamma, 1.e-6f, "gamma threshold value" );
      expectNear( failures,
                  grid.outputVariables().rho.getArray().getElement( 1 ),
                  testCase.expectedDensity,
                  1.e-5f,
                  "gamma threshold never falls back to coarse density" );
      expectNear( failures,
                  grid.outputVariables().v.getArray().getElement( 1 )[ 0 ],
                  testCase.expectedDensity == 0.f ? 0.f : 7.f,
                  1.e-5f,
                  "gamma threshold never falls back to coarse velocity" );
   }
}

void
checkPressureBoundaryAndHistory( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto coarseFluid = makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 4.f }, 2.f );
   auto fineFluid = makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 100.f }, 1.f );
   auto coarseBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 12.f }, 2.f );
   auto fineBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 300.f }, 1.f );
   auto ownPatch = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point{ 1.f, 0.f } }, { 8.f }, 2.f );
   auto reverse = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point{ 1.f, 0.f } }, { 1000.f }, 1.f );
   std::vector< Simulation2::MultiresolutionBoundaryPointer > patches = { ownPatch, reverse };
   const std::vector< InterfaceDescriptor > interfaces = { { 0, 1 }, { 1, 0 } };
   auto coarse = makeLevel< 2 >(
      coarseFluid, coarseBoundary, patches, interfaces, 0, makeOwnerBox< 2 >( 0.f, -1.f, 4.f, 1.f, 2.f ), params );
   auto fine = makeLevel< 2 >(
      fineFluid, fineBoundary, patches, interfaces, 1, makeOwnerBox< 2 >( 1.f, -1.f, 2.f, 1.f, 1.f ), params );
   std::vector< Point > positions = { Point{ 1.f, 0.f }, Point{ 20.f, 20.f } };

   TestPressureSensor< 2 > sensor;
   sensor.init( positions, 2, 2, false );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse }, params );

   // With linear EOS p=rho, coarse own fluid and own patch produce gamma=8+4=12
   // and pressure numerator=32+32=64, independently giving p=16/3.
   expectNear( failures, sensor.data().getElement( 0, 0 ), 16.f / 3.f, 1.e-5f, "single-view pressure normalization" );
   expectNear( failures, sensor.data().getElement( 0, 1 ), 0.f, 1.e-6f, "single-view uncovered pressure is zero" );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse, fine }, params );
   expect( failures, sensor.historySize() == 2, "single- and multi-view events append exactly one pressure row each" );
   expectNear( failures, sensor.data().getElement( 1, 0 ), 16.f / 3.f, 1.e-5f, "other-level pressure sources are excluded" );
   expectNear( failures, sensor.data().getElement( 1, 1 ), 0.f, 1.e-6f, "multi-view uncovered pressure is zero" );

   TestPressureSensor< 2 > withBoundary;
   withBoundary.init( positions, 2, 1, true );
   withBoundary.interpolate( std::vector< LevelView< 2 > >{ coarse, fine }, params );
   // Same-level boundary adds gamma=8/3 and numerator=32: p=96/(44/3)=72/11.
   expectNear( failures,
               withBoundary.data().getElement( 0, 0 ),
               72.f / 11.f,
               1.e-5f,
               "pressure optionally adds only own boundary" );
   expect( failures, withBoundary.historySize() == 1, "boundary-enabled pressure event appends one row" );
}

void
checkWaterOwnershipSourcesAndHistory( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   params.mass = 16.f;
   auto coarseFluid = makeSet< Simulation2::FluidPointer >(
      { Point{ 0.5f, -1.f }, Point{ 0.5f, 0.f } }, { 2.f, 2.f }, 0.5f );
   auto fineFluid = makeSet< Simulation2::FluidPointer >(
      { Point{ 0.5f, 1.f }, Point{ 0.5f, 2.f } }, { 1.f, 1.f }, 0.25f );
   auto coarseBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 0.5f, 1.f } }, { 2.f }, 0.5f );
   auto fineBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 0.5f, 2.f } }, { 2.f }, 0.25f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   auto coarse = makeLevel< 2 >(
      coarseFluid, coarseBoundary, noPatches, noInterfaces, 0, makeOwnerBox< 2 >( 0.f, -2.f, 1.f, 1.f, 0.5f ), params );
   auto fine = makeLevel< 2 >(
      fineFluid, fineBoundary, noPatches, noInterfaces, 1, makeOwnerBox< 2 >( 0.f, 1.f, 1.f, 3.f, 0.25f ), params );

   std::vector< Point > positions = { Point{ 0.5f, -1.f } };
   TestWaterSensor< 2 > sensor;
   sensor.init( positions, 2, 1.f, Point{ 0.f, 1.f }, -1.f, 3.f );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse }, params );
   // Coarse-only wet samples are y=-1 and y=0: two filled levels give water=1.
   expectNear( failures, sensor.getSensorData().getElement( 0, 0 ), 1.f, 1.e-6f, "single-view water reference" );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse, fine }, params );

   // Per-sample ownership: coarse fills y=-1,0; y=1 is coarse-owned but dry;
   // fine fills y=2. Fine's particle at y=1 and coarse boundary at y=1 are not
   // sources. Three filled levels plus sensorProjection=-1 give water=2.
   expectNear( failures,
               sensor.getSensorData().getElement( 1, 0 ),
               2.f,
               1.e-6f,
               "water uses per-vertical-sample finest ownership" );
   expect( failures, sensor.historySize() == 2, "single- and multi-view water events append exactly one row each" );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse, fine }, params );
   expectNear( failures, sensor.getSensorData().getElement( 2, 0 ), 2.f, 1.e-6f, "second multi-view water value" );
   expect( failures, sensor.historySize() == 3, "each multi-view water event appends one contiguous history row" );
}

void
checkWaterDirectedPatchSources( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   params.mass = 16.f;
   auto coarseFluid = makeSet< Simulation2::FluidPointer >( { Point{ -15.5f, -8.5f } }, { 2.f }, 0.5f );
   auto fineFluid = makeSet< Simulation2::FluidPointer >( { Point{ -15.5f, -12.f } }, { 1.f }, 0.25f );
   auto coarseBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ -15.5f, -8.5f } }, { 2.f }, 0.5f );
   auto fineBoundary = makeSet< Simulation2::BoundaryPointer >( { Point{ -15.5f, -8.5f } }, { 2.f }, 0.25f );
   auto ownPatch = makeSet< Simulation2::MultiresolutionBoundaryPointer >( { Point{ -15.5f, -15.f } }, { 2.f }, 0.5f );
   auto reversePatch = makeSet< Simulation2::MultiresolutionBoundaryPointer >(
      { Point{ -15.5f, -15.f }, Point{ -15.5f, -14.f }, Point{ -15.5f, -13.f } }, { 1.f, 1.f, 1.f }, 0.25f );
   std::vector< Simulation2::MultiresolutionBoundaryPointer > patches = { ownPatch, reversePatch };
   const std::vector< InterfaceDescriptor > interfaces = { { 0, 1 }, { 1, 0 } };
   auto coarse = makeLevel< 2 >( coarseFluid,
                                 coarseBoundary,
                                 patches,
                                 interfaces,
                                 0,
                                 makeOwnerBox< 2 >( -16.f, -15.5f, -14.f, -13.f, 0.5f ),
                                 params );
   auto fine = makeLevel< 2 >( fineFluid,
                               fineBoundary,
                               patches,
                               interfaces,
                               1,
                               makeOwnerBox< 2 >( -16.f, -13.f, -14.f, -11.f, 0.25f ),
                               params );

   std::vector< Point > positions = { Point{ -15.5f, -15.f } };
   TestWaterSensor< 2 > sensor;
   sensor.init( positions, 1, 1.f, Point{ 0.f, 1.f }, 0.f, 4.f );
   sensor.interpolate( std::vector< LevelView< 2 > >{ coarse, fine }, params );

   // Coarse own fluid is dry. Its directed patch supplies gamma=2.5 only at
   // y=-15; fine fluid supplies gamma=1.125 at y=-12. Exactly two samples are
   // wet, giving water=2-15=-13. The reverse patch is collocated at y=-15 and
   // also occupies y=-14,-13, which must remain dry. Omitting the own patch
   // gives water=-14; including reverse/other-level patches gives water=-11.
   expectNear( failures,
               sensor.getSensorData().getElement( 0, 0 ),
               -13.f,
               1.e-6f,
               "water includes ownIdx patches and excludes reverse/other-level patches" );
   expect( failures, sensor.historySize() == 1, "directed-patch water event appends one history row" );
}

void
printGridValues( const char* label, TestGrid< 2 >& grid )
{
   std::cout << "CHAR " << label << std::setprecision( 9 );
   for( int i = 0; i < grid.outputGamma().getSize(); ++i ) {
      const auto velocity = grid.outputVariables().v.getArray().getElement( i );
      std::cout << " " << i << ":(g=" << grid.outputGamma().getElement( i )
                << ",r=" << grid.outputVariables().rho.getArray().getElement( i ) << ",vx=" << velocity[ 0 ]
                << ",vy=" << velocity[ 1 ] << ")";
   }
   std::cout << '\n';
}

void
printPressureHistory( const char* label, const TestPressureSensor< 2 >& sensor )
{
   std::cout << "CHAR " << label << std::setprecision( 9 );
   for( int row = 0; row < sensor.historySize(); ++row )
      std::cout << " " << sensor.data().getElement( row, 0 );
   std::cout << '\n';
}

void
printWaterHistory( const char* label, const TestWaterSensor< 2 >& sensor )
{
   std::cout << "CHAR " << label << std::setprecision( 9 );
   for( int row = 0; row < sensor.historySize(); ++row )
      std::cout << " " << sensor.getSensorData().getElement( row, 0 );
   std::cout << '\n';
}

void
checkDirectGridOneLevelParity( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 2.f, 0.f } }, { 2.f }, { Point{ 3.f, -4.f } }, 1.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 20.f, 20.f } }, { 100.f }, 1.f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   auto level = makeLevel< 2 >(
      fluid, boundary, noPatches, noInterfaces, 17, makeOwnerBox< 2 >( -1.f, -1.f, 3.f, 1.f, 1.f ), params );
   expect( failures, level.levelId == 17, "parity fixture deliberately uses nonidentity domain level id" );

   TestGrid< 2 > direct;
   TestGrid< 2 > viewed;
   initializeCellGrid< 2 >( direct, Point{ 0.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 3, 1 } );
   initializeCellGrid< 2 >( viewed, Point{ 0.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 3, 1 } );
   direct.interpolateUsingGrid( fluid, boundary, params );
   viewed.interpolate( std::vector< LevelView< 2 > >{ level }, params );

   // W=1+h=1.5 and rf=searchRadius/(2h)=1. At query indices 1 and 2,
   // gamma=W*m/rho=3, rho numerator=W*m=6, and v numerator=v*gamma.
   // Query zero is outside support. Its same-set index collision is harmless;
   // both wet query indices differ from the sole source index zero.
   auto checkExpected = [ & ]( TestGrid< 2 >& grid, const std::string& path )
   {
      for( int i = 0; i < 3; ++i ) {
         const bool wet = i != 0;
         expectNear( failures, grid.outputGamma().getElement( i ), wet ? 3.f : 0.f, 1.e-6f, path + " raw gamma" );
         expectNear( failures,
                     grid.outputVariables().rho.getArray().getElement( i ),
                     wet ? 2.f : 0.f,
                     1.e-6f,
                     path + " density" );
         const auto velocity = grid.outputVariables().v.getArray().getElement( i );
         expectNear( failures, velocity[ 0 ], wet ? 3.f : 0.f, 1.e-6f, path + " velocity x" );
         expectNear( failures, velocity[ 1 ], wet ? -4.f : 0.f, 1.e-6f, path + " velocity y" );
      }
   };
   checkExpected( direct, "direct one-level grid" );
   checkExpected( viewed, "view one-level grid" );
   for( int i = 0; i < 3; ++i ) {
      expect( failures,
              direct.outputGamma().getElement( i ) == viewed.outputGamma().getElement( i ),
              "direct/view grid raw gamma exact parity " + std::to_string( i ) );
      expect( failures,
              direct.outputVariables().rho.getArray().getElement( i )
                 == viewed.outputVariables().rho.getArray().getElement( i ),
              "direct/view grid density exact parity " + std::to_string( i ) );
      expect( failures,
              direct.outputVariables().v.getArray().getElement( i )
                 == viewed.outputVariables().v.getArray().getElement( i ),
              "direct/view grid velocity exact parity " + std::to_string( i ) );
   }
   printGridValues( "grid parity direct", direct );
   printGridValues( "grid parity one-level(levelId=17,slot=0)", viewed );
}

void
checkDirectPressureUsesBaseParameters( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   params.mass = 0.2f;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 0.f, 0.f } }, { 2.f }, 2.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 20.f, 20.f } }, { 2.f }, 2.f );
   std::vector< Point > positions = { Point{ 0.f, 0.f } };
   TestPressureSensor< 2 > sensor;
   sensor.init( positions, 1, 1, false );
   sensor.interpolate( fluid, boundary, params );

   // Base h=0.5 and mass=0.2 give gamma=1.5*0.2/2=0.15, which is dry.
   // Radius-derived h=1 and mass=0.8 instead give gamma=0.8 and pressure=2.
   expectNear( failures,
               sensor.data().getElement( 0, 0 ),
               0.f,
               1.e-6f,
               "direct pressure uses base h and mass without refinement" );
   expect( failures, sensor.historySize() == 1, "base-pressure interpolation appends one history row" );
}

void
checkDirectPressureOneLevelParity( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 1.f, 0.f } }, { 4.f }, 1.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 1.f, 0.f } }, { 8.f }, 1.f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   auto level = makeLevel< 2 >(
      fluid, boundary, noPatches, noInterfaces, 23, makeOwnerBox< 2 >( 0.f, -1.f, 2.f, 1.f, 1.f ), params );
   std::vector< Point > positions = { Point{ 1.f, 0.f } };

   for( const bool includeBoundary : { false, true } ) {
      TestPressureSensor< 2 > direct;
      TestPressureSensor< 2 > viewed;
      direct.init( positions, 1, 1, includeBoundary );
      viewed.init( positions, 1, 1, includeBoundary );
      for( int call = 0; call < 2; ++call ) {
         direct.interpolate( fluid, boundary, params );
         viewed.interpolate( std::vector< LevelView< 2 > >{ level }, params );
      }

      // Direct path uses only base h=0.5, mass=4, W=1.5 with no rf-scaling:
      // fluid gamma=1.5 and pressure numerator=6 give p=4 without boundary.
      // Adding same-set boundary raises gamma by 0.75 to 2.25 and numerator by 6
      // to 12, preserving p=12/2.25=16/3.
      const float expected = includeBoundary ? 16.f / 3.f : 4.f;
      expect( failures, direct.historySize() == 2, "direct pressure appends exactly one row per call" );
      expect( failures, viewed.historySize() == 2, "one-level pressure appends exactly one row per call" );
      for( int row = 0; row < 2; ++row ) {
         const float directValue = direct.data().getElement( row, 0 );
         const float viewValue = viewed.data().getElement( row, 0 );
         expectNear( failures, directValue, expected, 1.e-6f, "direct pressure independent value" );
         expectNear( failures, viewValue, expected, 1.e-6f, "one-level pressure independent value" );
         expect( failures, directValue == viewValue, "direct/view pressure exact parity" );
      }
      const char* mode = includeBoundary ? "with-boundary" : "fluid-only";
      printPressureHistory( ( std::string( "pressure parity direct " ) + mode ).c_str(), direct );
      printPressureHistory( ( std::string( "pressure parity one-level " ) + mode ).c_str(), viewed );
   }
}

void
checkDirectWaterOneLevelParity( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto fluid = makeSet< Simulation2::FluidPointer >(
      { Point{ 0.5f, -1.f }, Point{ 0.5f, 0.f } }, { 2.f, 2.f }, 1.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 0.5f, 2.f } }, { 2.f }, 1.f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   auto level = makeLevel< 2 >(
      fluid, boundary, noPatches, noInterfaces, 29, makeOwnerBox< 2 >( 0.f, -2.f, 1.f, 3.f, 1.f ), params );
   std::vector< Point > positions = { Point{ 0.5f, -1.f } };
   TestWaterSensor< 2 > direct;
   TestWaterSensor< 2 > viewed;
   direct.init( positions, 1, 1.f, Point{ 0.f, 1.f }, -1.f, 3.f );
   viewed.init( positions, 1, 1.f, Point{ 0.f, 1.f }, -1.f, 3.f );
   for( int call = 0; call < 2; ++call ) {
      direct.interpolate( fluid, boundary, params );
      viewed.interpolate( std::vector< LevelView< 2 > >{ level }, params );
   }

   // The fluid fills samples y=-1,0,1 and leaves y=2 dry; the coincident
   // boundary at y=2 is excluded in both paths. The comma-overloaded dot
   // product is (sensorLocation,direction)=-1, hence water=3*1-1=2.
   expect( failures, direct.historySize() == 2, "direct water appends exactly one row per call" );
   expect( failures, viewed.historySize() == 2, "one-level water appends exactly one row per call" );
   for( int row = 0; row < 2; ++row ) {
      const float directValue = direct.getSensorData().getElement( row, 0 );
      const float viewValue = viewed.getSensorData().getElement( row, 0 );
      expectNear( failures, directValue, 2.f, 1.e-6f, "direct water independent projected result" );
      expectNear( failures, viewValue, 2.f, 1.e-6f, "one-level water independent projected result" );
      expect( failures, directValue == viewValue, "direct/view water exact parity" );
   }
   printWaterHistory( "water parity direct", direct );
   printWaterHistory( "water parity one-level", viewed );
}

void
checkDirectViewGammaThresholds( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   struct Case
   {
      const char* name;
      float mass;
      float gamma;
      float rho;
      float velocityX;
      float velocityY;
      float pressure;
      float water;
   };
   const Case cases[] = {
      { "below", 0.99f, 0.495f, 0.f, 0.f, 0.f, 0.f, 0.f },
      { "equal", 1.f, 0.5f, 0.f, 0.f, 0.f, 0.f, 0.f },
      { "above", 1.01f, 0.505f, 3.f, 4.f, -2.f, 3.f, 2.f },
   };

   for( const auto& testCase : cases ) {
      ModelParams< 2 > params;
      params.mass = testCase.mass;
      auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 2.f, 0.f } }, { 3.f }, { Point{ 4.f, -2.f } }, 1.f );
      auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 20.f, 20.f } }, { 30.f }, 1.f );
      const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
      const std::vector< InterfaceDescriptor > noInterfaces;
      auto level = makeLevel< 2 >(
         fluid, boundary, noPatches, noInterfaces, 31, makeOwnerBox< 2 >( -1.f, -1.f, 3.f, 1.f, 1.f ), params );

      TestGrid< 2 > directGrid;
      TestGrid< 2 > viewedGrid;
      initializeCellGrid< 2 >( directGrid, Point{ 0.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 3, 1 } );
      initializeCellGrid< 2 >( viewedGrid, Point{ 0.f, 0.f }, typename TestGrid< 2 >::CoordinatesType{ 3, 1 } );
      directGrid.interpolateUsingGrid( fluid, boundary, params );
      viewedGrid.interpolate( std::vector< LevelView< 2 > >{ level }, params );

      // With rho=3 and W=1.5, gamma=mass/2 exactly by construction. The raw
      // gamma survives even when normalized fields use the strict >0.5 test.
      for( TestGrid< 2 >* grid : { &directGrid, &viewedGrid } ) {
         expectNear( failures, grid->outputGamma().getElement( 2 ), testCase.gamma, 1.e-6f, "grid threshold raw gamma" );
         expectNear( failures,
                     grid->outputVariables().rho.getArray().getElement( 2 ),
                     testCase.rho,
                     1.e-6f,
                     "grid threshold density" );
         const auto velocity = grid->outputVariables().v.getArray().getElement( 2 );
         expectNear( failures, velocity[ 0 ], testCase.velocityX, 1.e-6f, "grid threshold velocity x" );
         expectNear( failures, velocity[ 1 ], testCase.velocityY, 1.e-6f, "grid threshold velocity y" );
      }
      expect( failures,
              directGrid.outputGamma().getElement( 2 ) == viewedGrid.outputGamma().getElement( 2 ),
              "threshold direct/view raw grid gamma exact parity" );
      expect( failures,
              directGrid.outputVariables().rho.getArray().getElement( 2 )
                 == viewedGrid.outputVariables().rho.getArray().getElement( 2 ),
              "threshold direct/view grid density exact parity" );
      expect( failures,
              directGrid.outputVariables().v.getArray().getElement( 2 )
                 == viewedGrid.outputVariables().v.getArray().getElement( 2 ),
              "threshold direct/view grid velocity exact parity" );

      std::vector< Point > positions = { Point{ 2.f, 0.f } };
      TestPressureSensor< 2 > directPressure;
      TestPressureSensor< 2 > viewedPressure;
      directPressure.init( positions, 1, 1, false );
      viewedPressure.init( positions, 1, 1, false );
      directPressure.interpolate( fluid, boundary, params );
      viewedPressure.interpolate( std::vector< LevelView< 2 > >{ level }, params );
      const float directP = directPressure.data().getElement( 0, 0 );
      const float viewedP = viewedPressure.data().getElement( 0, 0 );
      expectNear( failures, directP, testCase.pressure, 1.e-6f, "direct pressure strict gamma threshold" );
      expectNear( failures, viewedP, testCase.pressure, 1.e-6f, "one-level pressure strict gamma threshold" );
      expect( failures, directP == viewedP, "threshold direct/view pressure exact parity" );

      // Two vertical samples at y=0,1 have the same gamma. Their projection is
      // zero, so dry gives zero and above-threshold gives two filled levels.
      TestWaterSensor< 2 > directWater;
      TestWaterSensor< 2 > viewedWater;
      directWater.init( positions, 1, 1.f, Point{ 0.f, 1.f }, 0.f, 2.f );
      viewedWater.init( positions, 1, 1.f, Point{ 0.f, 1.f }, 0.f, 2.f );
      directWater.interpolate( fluid, boundary, params );
      viewedWater.interpolate( std::vector< LevelView< 2 > >{ level }, params );
      const float directWl = directWater.getSensorData().getElement( 0, 0 );
      const float viewedWl = viewedWater.getSensorData().getElement( 0, 0 );
      expectNear( failures, directWl, testCase.water, 1.e-6f, "direct water strict gamma threshold" );
      expectNear( failures, viewedWl, testCase.water, 1.e-6f, "one-level water strict gamma threshold" );
      expect( failures, directWl == viewedWl, "threshold direct/view water exact parity" );

      std::cout << "CHAR gamma " << testCase.name << std::setprecision( 9 ) << " grid direct=(g="
                << directGrid.outputGamma().getElement( 2 ) << ",r="
                << directGrid.outputVariables().rho.getArray().getElement( 2 ) << ") one-level=(g="
                << viewedGrid.outputGamma().getElement( 2 ) << ",r="
                << viewedGrid.outputVariables().rho.getArray().getElement( 2 ) << ") pressure=(" << directP << ","
                << viewedP << ") water=(" << directWl << "," << viewedWl << ")\n";
   }
}

void
checkMismatchedRadiusPressureSemantics( int& failures )
{
   using Simulation2 = Simulation< 2 >;
   using Point = Simulation2::SPHTraitsType::VectorType;
   ModelParams< 2 > params;
   auto fluid = makeSet< Simulation2::FluidPointer >( { Point{ 0.f, 0.f } }, { 4.f }, 1.f );
   auto boundary = makeSet< Simulation2::BoundaryPointer >( { Point{ 1.5f, 0.f } }, { 8.f }, 2.f );
   const std::vector< Simulation2::MultiresolutionBoundaryPointer > noPatches;
   const std::vector< InterfaceDescriptor > noInterfaces;
   auto level = makeLevel< 2 >(
      fluid, boundary, noPatches, noInterfaces, 37, makeOwnerBox< 2 >( -1.f, -1.f, 1.f, 1.f, 1.f ), params );
   std::vector< Point > positions = { Point{ 0.f, 0.f } };
   TestPressureSensor< 2 > direct;
   TestPressureSensor< 2 > viewed;
   direct.init( positions, 1, 1, true );
   viewed.init( positions, 1, 1, true );
   direct.interpolate( fluid, boundary, params );
   viewed.interpolate( std::vector< LevelView< 2 > >{ level }, params );

   // Legacy direct pressure applies the fluid radius=1 to both sources and
   // excludes the boundary at distance 1.5: p=4. MR applies the boundary's own
   // radius=2 and includes it: p=(6+6)/(1.5+0.75)=16/3. This intentional
   // semantic difference must not be unified accidentally during deduplication.
   const float directValue = direct.data().getElement( 0, 0 );
   const float viewValue = viewed.data().getElement( 0, 0 );
   expectNear( failures, directValue, 4.f, 1.e-6f, "direct pressure retains fluid cutoff for boundary source" );
   expectNear( failures, viewValue, 16.f / 3.f, 1.e-6f, "one-level pressure uses boundary source cutoff" );
   expect( failures, directValue != viewValue, "mismatched-radius pressure results are intentionally distinct" );
   expect( failures, direct.historySize() == 1 && viewed.historySize() == 1, "mismatched-radius calls append one row" );
   std::cout << "CHAR pressure mismatched-radius direct=" << std::setprecision( 9 ) << directValue
             << " one-level=" << viewValue << '\n';
}

void
runWave1MeasurementCharacterizationTests( int& failures )
{
   checkDirectGridOneLevelParity( failures );
   checkDirectPressureOneLevelParity( failures );
   checkDirectPressureUsesBaseParameters( failures );
   checkDirectWaterOneLevelParity( failures );
   checkDirectViewGammaThresholds( failures );
   checkMismatchedRadiusPressureSemantics( failures );
}

void
runApprovedMeasurementTests( int& failures )
{
   checkHalfOpenOwnership( failures );
   checkNestedDisjointAndReversedOwnership( failures );
   checkOverlapValidation( failures );
   checkViewScalingAndPatchSelection( failures );
   checkGridNormalizationAndSources( failures );
   checkDryFineAndUncoveredGrid( failures );
   checkGammaThresholdWithoutFallback( failures );
   checkPressureBoundaryAndHistory( failures );
   checkWaterOwnershipSourcesAndHistory( failures );
   checkWaterDirectedPatchSources( failures );
   checkSolverMeasurementBuilder( failures );
   checkSolverMeasurementDispatch( failures );
}

}  // namespace

int
main()
{
   int failures = 0;
   MultiresolutionMeasurementTest::runLegacyMeasurementTests( failures );
   runWave1MeasurementCharacterizationTests( failures );
   runApprovedMeasurementTests( failures );
   if( failures != 0 ) {
      std::cerr << failures << " multiresolution measurement checks failed\n";
      return 1;
   }
   std::cout << "PASS: legacy and approved multiresolution measurement behavior\n";
   return 0;
}
