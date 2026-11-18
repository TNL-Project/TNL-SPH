#include <filesystem>
#include <iostream>
#include <string>
#include <type_traits>

#include "MultiresolutionSearchSchedulingTestSupport.h"

namespace {

using namespace MultiresolutionSearchSchedulingTest;

template< typename Solver >
void
checkSavePressureCurrentDensity( int& failures, const std::string& label )
{
   Solver solver;
   setupTwoWay( solver );
   solver.numberOfSubsets = 0;
   solver.outputDirectory = "results";
   auto& patch = Access< Solver >::patches( solver )[ 0 ];
   expect( failures, patch->getNumberOfParticles() == 1, label + ": fixture has one active patch particle" );
   for( const float density : { 1111.5f, 1222.5f } ) {
      patch->getVariables()->rho.getArray().setElement( 0, density );
      patch->getVariables()->p.getArray().setElement( 0, -777.f );
      solver.save( false );
      expect( failures, patch->getVariables()->p.getArray().getElement( 0 ) == density - solver.modelParams.rho0, label );
   }
}

void
checkLtsAccumulation( int& failures )
{
   LtsSolver solver;
   solver.numberOfSubsets = 2;
   solver.fluidSets.resize( 2 );
   solver.fluidSets[ 0 ] = makeDomainSet< FluidPointer >( { 10, 10 }, { 0, 0 }, {} );
   solver.fluidSets[ 1 ] = makeDomainSet< FluidPointer >( { 6, 6 }, { 2, 2 }, { { 7.5f, 7.5f } } );
   solver.multiresolutionBoundaryPatches.resize( 1 );
   solver.multiresolutionBoundaryPatchInterfaces = { { 0, 1 } };
   auto& patch = solver.multiresolutionBoundaryPatches[ 0 ];
   patch = makeDomainSet< LtsSolver::MultiresolutionBoundaryPointer >( { 10, 10 }, { 0, 0 }, {} );
   patch->initZones( solver.fluidSets[ 0 ]->getParticles(), solver.fluidSets[ 1 ]->getParticles(), 1.f );
   patch->initMassNodes( solver.modelParams, 0, 1, 1.f );
   expect( failures, patch->getMassNodes().numberOfMassNodes > 0, "LTS accumulation fixture has mass nodes" );
   if( patch->getMassNodes().numberOfMassNodes == 0 )
      return;

   const Point moved = patch->getMassNodes().points.getElement( 0 );
   const Point normal = patch->getMassNodes().normal.getElement( 0 );
   solver.fluidSets[ 1 ]->getPoints().setElement( 0, moved );
   solver.fluidSets[ 1 ]->getVariables()->v.getArray().setElement( 0, normal );
   solver.accumulateInterfaceMassesFromNeighbor( 1, 0.01f );

   float flux = 0.f;
   for( int i = 0; i < patch->getMassNodes().numberOfMassNodes; ++i )
      flux += std::abs( patch->getMassNodes().massFlux.getElement( i ) );
   expect( failures, flux > 0.f, "LTS accumulation rebuilds the integrated source before all matching patches" );
}

void
checkEmptyAccumulation( int& failures )
{
   LtsSolver solver;
   solver.fluidSets.resize( 1 );
   solver.fluidSets[ 0 ] = makeDomainSet< FluidPointer >( { 4, 4 }, { 0, 0 }, { { 0.5f, 0.5f } } );
   solver.fluidSets[ 0 ]->getPoints().setElement( 0, Point{ FLT_MAX, FLT_MAX } );
   solver.accumulateInterfaceMassesFromNeighbor( 0, 0.01f );
   expect( failures,
           solver.fluidSets[ 0 ]->getParticles()->getNumberOfParticlesToRemove() == 0,
           "LTS accumulation skips all source work when no patch uses the level" );
}

void
checkLtsCompatibilityAliases( int& failures )
{
   using BufferAlias = typename LtsSolver::MultiresolutionBoundaryLTS;
   using PointerAlias = typename LtsSolver::MultiresolutionBoundaryLTSPointer;
   static_assert( std::is_same_v< BufferAlias, typename LtsSolver::MultiresolutionBoundary > );
   static_assert( std::is_same_v< PointerAlias, typename LtsSolver::MultiresolutionBoundaryPointer > );

   LtsSolver solver;
   solver.multiresolutionBoundaryPatchesLTS.resize( 1 );
   expect( failures, &solver.multiresolutionBoundaryPatchesLTS == &solver.multiresolutionBoundaryPatches,
           "LTS compatibility vector aliases canonical storage" );
   expect( failures, solver.multiresolutionBoundaryPatches.size() == 1,
           "LTS compatibility vector updates canonical storage" );
}

}  // namespace

int
main()
{
   std::filesystem::create_directories( "results" );
   int failures = 0;
   checkSingleDirection< UniformSolver >( failures, "uniform" );
   checkSecondDirectionCompactsSource< UniformSolver >( failures, "uniform" );
   checkTwoWay< UniformSolver >( failures, "uniform" );
   checkSavePressureCurrentDensity< UniformSolver >( failures, "uniform patch save pressure uses current density" );
   checkSingleDirection< LtsSolver >( failures, "LTS" );
   checkSecondDirectionCompactsSource< LtsSolver >( failures, "LTS" );
   checkTwoWay< LtsSolver >( failures, "LTS" );
   checkSavePressureCurrentDensity< LtsSolver >( failures, "LTS patch save pressure uses current density" );
   checkLtsAccumulation( failures );
   checkEmptyAccumulation( failures );
   checkLtsCompatibilityAliases( failures );
   if( failures != 0 ) {
      std::cerr << failures << " multiresolution search-scheduling checks failed\n";
      return 1;
   }
   std::cout << "PASS: moved sources are rebuilt per direction and final fluid searches exclude captures\n";
   return 0;
}
