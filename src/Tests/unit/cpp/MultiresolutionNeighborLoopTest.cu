#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include <TNL/Algorithms/Segments/Ellpack.h>
#include <TNL/Devices/Host.h>
#include <TNL/Particles/CellIndexer.h>
#include <TNL/Particles/ParticlesCellList.h>
#include <TNL/Particles/Writers/VTKWriter.h>
#include <TNL/Pointers/SharedPointer.h>

#include <SPH/shared/thrustExecPolicySelector.h>
#include <SPH/Fluid.h>
#include <SPH/Kernels.h>
#include <SPH/MultiresolutionRectangleBuffer.h>
#include <SPH/MultiresolutionRectangleBufferLocalTimeStepping.h>
#include <SPH/OpenBoundaryConfig.h>

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
using Point = Traits::VectorType;
using Coordinates = Traits::CoordinatesType;

class Variables
{
public:
   using ScalarArrayType = Traits::ScalarArrayType;
   using VectorArrayType = Traits::VectorArrayType;

   void
   setSize( const int size )
   {
      rho.setSize( size );
      v.setSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& )
   {}

   ScalarArrayType rho;
   VectorArrayType v;
};

class IntegratorVariables
{
public:
   void
   setSize( const int size )
   {
      rho_old.setSize( size );
      v_old.setSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& )
   {}

   Traits::ScalarArrayType rho_old;
   Traits::VectorArrayType v_old;
};

struct ModelParams
{
   using KernelFunction = TNL::SPH::KernelFunctions::WendlandKernel< CaseConfig >;

   float h = 0.5f;
   float dp = 0.5f;
   // rho0 * dp^2: m/rho = dp^2 keeps the real interpolation support test above its M00 threshold.
   float mass = 250.f;
   float mdbcExtrapolationDetTreshold = 0.001f;
};

using Particles = TNL::Particles::ParticlesCellList< ParticlesConfig, Device >;
using Fluid = TNL::SPH::Fluid< Particles, CaseConfig, Variables, IntegratorVariables >;
using OpenBoundaryConfig = TNL::SPH::OpenBoundaryConfig< CaseConfig >;
using Buffer = TNL::SPH::MultiresolutionBoundary< Particles, CaseConfig, Variables, IntegratorVariables, OpenBoundaryConfig, ModelParams >;
using BufferLTS = TNL::SPH::MultiresolutionBoundaryLTS< Particles, CaseConfig, Variables, IntegratorVariables, OpenBoundaryConfig, ModelParams >;
using FluidPointer = TNL::Pointers::SharedPointer< Fluid, Device >;

template< typename SetPointer >
void
initializeSet( SetPointer& set )
{
   auto& particles = set->getParticles();
   particles->setSize( 1 );
   particles->setSearchRadius( 1.f );
   particles->setDimensions( Coordinates{ 2, 2 } );
   particles->setOrigin( Point{ 0.f, 0.f } );
   particles->setGlobalOriginCoordinates( Coordinates{ 0, 0 } );
   particles->setNumberOfParticles( 1 );
   set->getVariables()->setSize( 1 );
   set->getIntegratorVariables()->setSize( 1 );
   particles->getPoints().setElement( 0, Point{ 0.5f, 0.5f } );
   set->searchForNeighbors();
}

template< typename BufferPointer >
struct Fixture
{
   using BufferPointerType = BufferPointer;

   Fixture()
   {
      initializeSet( fluid );
      initializeSet( buffer );
      fluid->getVariables()->rho.setElement( 0, 1000.f );
      fluid->getVariables()->v.setElement( 0, Point{ 1.f, 0.f } );
      buffer->getVariables()->rho.setElement( 0, -7.f );
      buffer->getVariables()->v.setElement( 0, Point{ -7.f, -7.f } );
      buffer->getIntegratorVariables()->rho_old.setElement( 0, 0.f );
      buffer->getIntegratorVariables()->v_old.setElement( 0, Point{ 0.f, 0.f } );
   }

   FluidPointer fluid;
   BufferPointer buffer;
   ModelParams modelParams;
};

template< typename BufferType >
void
checkLoopVisits( int& failures )
{
   using BufferPointer = TNL::Pointers::SharedPointer< BufferType, Device >;
   Fixture< BufferPointer > fixture;
   auto searchInFluid = fixture.buffer->getParticles()->getSearchToken( fixture.fluid->getParticles() );
   const Point query = fixture.buffer->getPoints().getElement( 0 );
   int sameSetVisits = 0;
   int crossSetVisits = 0;
   auto countVisit = [] __cuda_callable__( int, int, Point&, int* visits ) mutable
   {
      ++( *visits );
   };
   Particles::NeighborsLoop::exec( 0, query, searchInFluid, countVisit, &sameSetVisits );
   Particles::NeighborsLoopAnotherSet::exec( 0, query, searchInFluid, countVisit, &crossSetVisits );
   if( sameSetVisits != 0 || crossSetVisits != 1 ) {
      ++failures;
      std::cerr << "FAIL: real TNL visits at i=j=0: NeighborsLoop=" << sameSetVisits
                << " (expected 0), NeighborsLoopAnotherSet=" << crossSetVisits << " (expected 1)\n";
   }
   else
      std::cout << "real TNL visits at i=j=0: NeighborsLoop=0, NeighborsLoopAnotherSet=1\n";
}

void
expectNear( int& failures, const std::string& method, const char* quantity, float actual, float expected )
{
   if( ! std::isfinite( actual ) || std::abs( actual - expected ) > 1.e-3f ) {
      ++failures;
      std::cerr << "FAIL: " << method << " skipped the valid cross-set neighbor at query i=0, source j=0: " << quantity
                << "=" << actual << ", expected " << expected << '\n';
   }
}

template< typename BufferType >
void
checkInterpolation( int& failures, const std::string& method )
{
   using BufferPointer = TNL::Pointers::SharedPointer< BufferType, Device >;
   Fixture< BufferPointer > fixture;
   fixture.buffer->interpolateVariables( fixture.fluid, fixture.modelParams );
   expectNear( failures, method, "rho", fixture.buffer->getVariables()->rho.getElement( 0 ), 1000.f );
   expectNear( failures, method, "v[0]", fixture.buffer->getVariables()->v.getElement( 0 )[ 0 ], 1.f );
}

void
checkLTSRefresh( int& failures )
{
   using BufferPointer = TNL::Pointers::SharedPointer< BufferLTS, Device >;
   Fixture< BufferPointer > fixture;
   fixture.buffer->refreshVariablesFromNeighbor( fixture.fluid, fixture.modelParams, 1.f );
   expectNear( failures, "LTS refreshVariablesFromNeighbor", "rho", fixture.buffer->getVariables()->rho.getElement( 0 ), 1000.f );
   expectNear( failures, "LTS refreshVariablesFromNeighbor", "v[0]", fixture.buffer->getVariables()->v.getElement( 0 )[ 0 ], 1.f );
}

}  // namespace

int
main()
{
   int failures = 0;
   checkLoopVisits< Buffer >( failures );
   checkInterpolation< Buffer >( failures, "buffer interpolateVariables" );
   checkInterpolation< BufferLTS >( failures, "LTS interpolateVariables" );
    checkLTSRefresh( failures );

    int nonfiniteFailures = 0;
    expectNear( nonfiniteFailures, "nonfinite guard", "rho", std::numeric_limits< float >::quiet_NaN(), 1000.f );
    if( nonfiniteFailures != 1 )
       ++failures;

    if( failures != 0 ) {
      std::cerr << failures << " cross-set neighbor regression checks failed\n";
      return 1;
   }
   std::cout << "PASS: both buffer variants visited the valid source particle at j == i == 0\n";
   return 0;
}
