#pragma once

#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

#include "MultiresolutionBufferCapacityTestSupport.h"
#include <TNL/Particles/Readers/VTKReader.h>
#include <SPH/solvers/SolverMultiSetBlockMultiresolutionLocalTimeStepping.h>

namespace MultiresolutionSearchSchedulingTest {

using namespace ParticleSetCapacityTest;

struct EmptyScheme {};
struct EmptyOpenBoundaryModel {};

struct SchedulingParams
{
   using SPHConfig = CaseConfig;
   using TimeStepping = TNL::SPH::ConstantTimeStep< CaseConfig >;
   using KernelFunction = TNL::SPH::KernelFunctions::WendlandKernel< CaseConfig >;
   float h = 0.5f, dp = 0.5f, mass = 250.f, rho0 = 1000.f;
   float mdbcExtrapolationDetTreshold = 0.001f;
};

struct SchedulingModel
{
   using DeviceType = Device;
   using ModelConfigType = SchedulingParams;
   using ModelParams = SchedulingParams;
   using ParticlesType = Particles;
   using IntegrationSchemeType = EmptyScheme;
   using IntegrationSchemeVariables = IntegratorVariables;
   using SPHConfig = CaseConfig;
   using FluidVariables = ParticleSetCapacityTest::FluidVariables;
   using BoundaryVariables = ParticleSetCapacityTest::BoundaryVariables;
   using OpenBoundaryVariables = ParticleSetCapacityTest::OpenBoundaryVariables;
   using OpenBoundaryConfig = TNL::SPH::OpenBoundaryConfig< CaseConfig >;
   using OpenBoundaryModel = EmptyOpenBoundaryModel;

   template< typename PhysicalObjectPointer >
   void
   computePressureFromDensity( PhysicalObjectPointer& physicalObject, SchedulingParams& modelParams )
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

using UniformSolver = TNL::SPH::SolverMultiSetBlockMultiresolution< SchedulingModel >;
using LtsSolver = TNL::SPH::SolverMultiSetBlockMultiresolutionLocalTimeStepping< SchedulingModel >;

template< typename Solver >
struct Access;

template<>
struct Access< UniformSolver >
{
   static auto& patches( UniformSolver& solver ) { return solver.multiresolutionBoundaryPatches; }
   static void updateFirst( UniformSolver& solver )
   {
      solver.multiresolutionBoundaryPatches[ 0 ]->updateInterfaceBuffer(
         solver.fluidSets[ 0 ], solver.fluidSets[ 1 ], solver.modelParams, 0.01f, 0 );
   }
   static void updateSecond( UniformSolver& solver )
   {
      solver.multiresolutionBoundaryPatches[ 1 ]->updateInterfaceBuffer(
         solver.fluidSets[ 1 ], solver.fluidSets[ 0 ], solver.modelParams, 0.01f, 1 );
   }
   static void update( UniformSolver& solver ) { solver.multiresolutionUpdate(); }
};

template<>
struct Access< LtsSolver >
{
   static auto& patches( LtsSolver& solver ) { return solver.multiresolutionBoundaryPatches; }
   static void updateFirst( LtsSolver& solver )
   {
      solver.multiresolutionBoundaryPatches[ 0 ]->syncUpdateInterfaceBuffer(
         solver.fluidSets[ 0 ], solver.fluidSets[ 1 ], solver.modelParams, 0.01f );
   }
   static void updateSecond( LtsSolver& solver )
   {
      solver.multiresolutionBoundaryPatches[ 1 ]->syncUpdateInterfaceBuffer(
         solver.fluidSets[ 1 ], solver.fluidSets[ 0 ], solver.modelParams, 0.01f );
   }
   static void update( LtsSolver& solver ) { solver.syncMultiresolutionUpdate(); }
};

template< typename SetPointer >
SetPointer
makeDomainSet( const Coordinates& dimensions,
               const Coordinates& globalOrigin,
               const std::vector< Point >& points,
               const int capacity = 16 )
{
   auto set = MultiresolutionBufferCapacityTest::makeSet< SetPointer >( points.size(), capacity );
   auto particles = set->getParticles();
   particles->setDimensions( dimensions );
   particles->setOrigin( Point( globalOrigin ) );
   particles->setReferentialOrigin( Point{ 0.f, 0.f } );
   particles->setGlobalOriginCoordinates( globalOrigin );
   particles->setOverlapWidth( 1 );
   for( int i = 0; i < int( points.size() ); ++i )
      particles->getPoints().setElement( i, points[ i ] );
   set->searchForNeighbors();
   return set;
}

inline int
findStableId( const FluidPointer& fluid, const int stableId )
{
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i )
      if( fluid->getVariables()->referentialIdx.getElement( i ) == stableId )
         return i;
   return -1;
}

template< typename Solver >
void
setupTwoWay( Solver& solver )
{
   solver.numberOfSubsets = 2;
   solver.fluidSets.resize( 2 );
   solver.fluidSets[ 0 ] = makeDomainSet< FluidPointer >(
      { 10, 10 }, { 0, 0 }, { { 4.5f, 2.5f }, { 0.5f, 0.5f } } );
   solver.fluidSets[ 1 ] = makeDomainSet< FluidPointer >(
      { 6, 6 }, { 2, 2 }, { { 7.5f, 7.5f }, { 1.5f, 1.5f } } );
   const int moved0 = findStableId( solver.fluidSets[ 0 ], 1037 );
   const int moved1 = findStableId( solver.fluidSets[ 1 ], 1000 );
   solver.fluidSets[ 0 ]->getVariables()->rho.getArray().setElement( moved0, 2222.f );
   solver.fluidSets[ 0 ]->getVariables()->v.getArray().setElement( moved0, Point{ 84.f, 0.f } );
   solver.fluidSets[ 1 ]->getVariables()->rho.getArray().setElement( moved1, 1111.f );
   solver.fluidSets[ 1 ]->getVariables()->v.getArray().setElement( moved1, Point{ 42.f, 0.f } );

   auto& patches = Access< Solver >::patches( solver );
   patches.resize( 2 );
   patches[ 0 ] = makeDomainSet< std::decay_t< decltype( patches[ 0 ] ) > >(
      { 10, 10 }, { 0, 0 }, { { 2.5f, 3.5f } } );
   patches[ 1 ] = makeDomainSet< std::decay_t< decltype( patches[ 1 ] ) > >(
      { 6, 6 }, { 2, 2 }, { { 8.5f, 8.5f } } );
   patches[ 0 ]->initZones( solver.fluidSets[ 0 ]->getParticles(), solver.fluidSets[ 1 ]->getParticles(), 1.f );
   patches[ 1 ]->initZones( solver.fluidSets[ 1 ]->getParticles(), solver.fluidSets[ 0 ]->getParticles(), 1.f );
   for( int p = 0; p < 2; ++p ) {
      patches[ p ]->getVariables()->referentialIdx.setElement( 0, 9000 + p );
      patches[ p ]->getVariables()->rho.getArray().setElement( 0, -7.f );
      patches[ p ]->getVariables()->v.getArray().setElement( 0, Point{ -0.01f, 0.f } );
   }
   solver.multiresolutionBoundaryPatchInterfaces = { { 0, 1 }, { 1, 0 } };
   solver.timeMeasurement.addTimer( "multiresolution-update" );
   solver.timeStepping.setTimeStep( 0.01f );

   // Real cell motion after all initial searches. Patch zero captures fluid 0's
   // other particle before patch one needs this moved source.
   solver.fluidSets[ 1 ]->getPoints().setElement( moved1, Point{ 2.5f, 3.5f } );
   solver.fluidSets[ 0 ]->getPoints().setElement( moved0, Point{ 8.5f, 8.5f } );
}

inline int
findStableParticle( const FluidPointer& fluid, const Point& point )
{
   for( int i = 0; i < fluid->getNumberOfParticles(); ++i )
      if( fluid->getPoints().getElement( i ) == point )
         return i;
   return -1;
}

template< typename BufferPointer >
float
interpolatedDensity( const BufferPointer& buffer, const Point& query )
{
   int nearest = -1;
   float distance = FLT_MAX;
   for( int i = 0; i < buffer->getNumberOfParticles(); ++i ) {
      const Point difference = buffer->getPoints().getElement( i ) - query;
      const float candidate = std::sqrt( ( difference, difference ) );
      if( candidate < distance ) {
         nearest = i;
         distance = candidate;
      }
   }
   return nearest >= 0 && distance < 0.01f ? buffer->getVariables()->rho.getArray().getElement( nearest ) : -999.f;
}

template< typename Solver >
void
checkSingleDirection( int& failures, const std::string& label )
{
   Solver solver;
   setupTwoWay( solver );
   Access< Solver >::updateFirst( solver );
   const float density = interpolatedDensity( Access< Solver >::patches( solver )[ 0 ], Point{ 2.5f, 3.5f } );
   expect( failures, std::abs( density - 1111.f ) < 0.01f, label + ": entry source rebuilt before interpolation" );
}

template< typename Solver >
void
checkSecondDirectionCompactsSource( int& failures, const std::string& label )
{
   Solver solver;
   setupTwoWay( solver );
   Access< Solver >::updateFirst( solver );
   expect( failures, solver.fluidSets[ 0 ]->getNumberOfParticles() == 2,
           label + ": first patch defers its captured source placeholder" );
   Access< Solver >::updateSecond( solver );
   expect( failures, solver.fluidSets[ 0 ]->getNumberOfParticles() == 1,
           label + ": second directed entry compacts the mutated source" );
   expect( failures, findStableParticle( solver.fluidSets[ 0 ], Point{ 8.5f, 8.5f } ) >= 0,
           label + ": second directed entry retains the surviving source" );
   expect( failures, solver.fluidSets[ 0 ]->getParticles()->getNumberOfParticlesToRemove() == 0,
           label + ": second directed entry clears the source removal count" );
}

template< typename Solver >
void
checkTwoWay( int& failures, const std::string& label )
{
   Solver solver;
   setupTwoWay( solver );
   Access< Solver >::update( solver );
   auto& patches = Access< Solver >::patches( solver );
   const float density0 = interpolatedDensity( patches[ 0 ], Point{ 2.5f, 3.5f } );
   const float density1 = interpolatedDensity( patches[ 1 ], Point{ 8.5f, 8.5f } );
   expect( failures, std::abs( density0 - 1111.f ) < 0.01f, label + ": first moved source interpolated" );
   expect( failures, std::abs( density1 - 2222.f ) < 0.01f, label + ": captured-source successor rebuilt" );

   const Point remaining[] = { { 8.5f, 8.5f }, { 2.5f, 3.5f } };
   for( int level = 0; level < 2; ++level ) {
      auto& fluid = solver.fluidSets[ level ];
      expect( failures, fluid->getNumberOfParticles() == 1, label + ": final active count excludes captures" );
      expect( failures, findStableParticle( fluid, remaining[ level ] ) >= 0, label + ": surviving particle retained" );
      expect( failures, fluid->getParticles()->getNumberOfParticlesToRemove() == 0, label + ": removal count reset" );
      const int index = findStableParticle( fluid, remaining[ level ] );
      if( index >= 0 ) {
         auto token = fluid->getParticles()->getSearchToken( fluid->getParticles() );
         int visits = 0;
         auto count = [] __cuda_callable__( int, int, Point&, int* n ) mutable { ++( *n ); };
         Particles::NeighborsLoopAnotherSet::exec( index, remaining[ level ], token, count, &visits );
         expect( failures, visits == 1, label + ": final search token visits surviving particle" );
      }
   }
}

}  // namespace MultiresolutionSearchSchedulingTest
