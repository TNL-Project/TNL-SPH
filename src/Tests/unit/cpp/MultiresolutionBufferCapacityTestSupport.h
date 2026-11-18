#pragma once

#include <cfloat>
#include <string>

#include "ParticleSetCapacityTestSupport.h"

#include <TNL/Particles/Writers/VTKWriter.h>

#include <SPH/shared/thrustExecPolicySelector.h>
#include <SPH/Kernels.h>
#include <SPH/MultiresolutionRectangleBuffer.h>
#include <SPH/MultiresolutionRectangleBufferLocalTimeStepping.h>
#include <SPH/OpenBoundaryConfig.h>

namespace MultiresolutionBufferCapacityTest {

using namespace ParticleSetCapacityTest;
using BoundarySettings = TNL::SPH::OpenBoundaryConfig< CaseConfig >;

struct ModelParams
{
   using KernelFunction = TNL::SPH::KernelFunctions::WendlandKernel< CaseConfig >;
   float h = 0.5f;
   float dp = 0.5f;
   float mass = 250.f;
   float rho0 = 1000.f;
};

template< typename Base >
struct BufferProbe : public Base
{
   using Base::frameDimensions;
   using Base::frameOrientation;
   using Base::frameOriginGlobalCoordinates;
   using Base::inner_overlap;
   using Base::massNodes;
   using Base::numberOfParticlesToBuffer;
   using Base::numberOfParticlesToCreate;
   using Base::numberOfParticlesToRemove;
   using Base::numberOfParticlesToRetype;
   using Base::outer_overlap;
   using Base::particlesToBuffer;
   using Base::particlesToFluid;
   using Base::particlesToRemove;
   using Base::retypeMarker;
   using Base::zone;
};

using UniformBase = TNL::SPH::
   MultiresolutionBoundary< Particles, CaseConfig, FluidVariables, IntegratorVariables, BoundarySettings, ModelParams >;
using LTSBase = TNL::SPH::
   MultiresolutionBoundaryLTS< Particles, CaseConfig, FluidVariables, IntegratorVariables, BoundarySettings, ModelParams >;
using Buffer = BufferProbe< UniformBase >;
using BufferLTS = BufferProbe< LTSBase >;
using BufferPointer = TNL::Pointers::SharedPointer< Buffer, Device >;
using BufferLTSPointer = TNL::Pointers::SharedPointer< BufferLTS, Device >;

template< typename SetPointer >
SetPointer
makeSet( const int activeCount, const int capacity )
{
   SetPointer set;
   auto particles = set->getParticles();
   particles->setSize( capacity );
   particles->setSearchRadius( 1.f );
   particles->setDimensions( Coordinates{ 4, 4 } );
   particles->setOrigin( Point{ 0.f, 0.f } );
   particles->setGlobalOriginCoordinates( Coordinates{ 0, 0 } );
   particles->setNumberOfParticles( activeCount );
   set->getVariables()->setSize( capacity );
   set->getIntegratorVariables()->setSize( capacity );

   auto& variables = *set->getVariables();
   auto& history = *set->getIntegratorVariables();
   for( int i = 0; i < activeCount; ++i ) {
      const float value = static_cast< float >( i + 1 );
      const Point point{ 0.5f + static_cast< float >( i % 2 ), 0.5f + static_cast< float >( ( i / 2 ) % 2 ) };
      particles->getPoints().setElement( i, point );
      variables.referentialIdx.setElement( i, 1000 + 37 * i );
      variables.rho.getArray().setElement( i, 1000.f + value );
      variables.drho.getArray().setElement( i, 10.f + value );
      variables.p.getArray().setElement( i, 100.f + value );
      variables.v.getArray().setElement( i, Point{ value, -value } );
      variables.a.getArray().setElement( i, Point{ 2.f * value, 3.f * value } );
      variables.gamma.getArray().setElement( i, 0.5f + value );
      variables.marker.getArray().setElement( i, 100 + i );
      history.rho_old.setElement( i, 900.f + value );
      history.v_old.setElement( i, Point{ -2.f * value, 0.25f * value } );
   }
   variables.highestReferentialIdx = 4321;
   set->searchForNeighbors();
   return set;
}

template< typename BufferPointerType >
BufferPointerType
makeBuffer( const int activeCount, const int capacity )
{
   BufferPointerType buffer = makeSet< BufferPointerType >( activeCount, capacity );
   buffer->particlesToFluid.setSize( capacity );
   buffer->particlesToRemove.setSize( capacity );
   buffer->retypeMarker.setSize( capacity );
   buffer->particlesToBuffer.setSize( 1 );
   buffer->massNodes.setSize( 1 );
   if( capacity > 0 ) {
      buffer->particlesToFluid.setElement( 0, 71 );
      buffer->particlesToRemove.setElement( 0, 72 );
      buffer->retypeMarker.setElement( 0, 73 );
   }
   buffer->inner_overlap = false;
   buffer->outer_overlap = true;
   buffer->frameOrientation = 1;
   buffer->frameOriginGlobalCoordinates = Coordinates{ 0, 0 };
   buffer->frameDimensions = Coordinates{ 4, 4 };
   return buffer;
}

template< typename SetPointer >
SnapshotMap
snapshotOf( const SetPointer& set )
{
   SnapshotMap result;
   const auto& particles = *set->getParticles();
   const auto& variables = *set->getVariables();
   const auto& history = *set->getIntegratorVariables();
   for( int i = 0; i < set->getNumberOfParticles(); ++i ) {
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

template< typename SetPointer >
void
checkPreserved( int& failures,
                const SetPointer& set,
                const SnapshotMap& expected,
                const std::string& label,
                const int expectedHighestReferentialIdx = 4321 )
{
   const auto& variables = *set->getVariables();
   const auto& history = *set->getIntegratorVariables();
   for( const auto& [ id, old ] : expected ) {
      int found = -1;
      for( int i = 0; i < set->getNumberOfParticles() && found < 0; ++i )
         if( variables.referentialIdx.getElement( i ) == id )
            found = i;
      expect( failures, found >= 0, label + ": stable ID " + std::to_string( id ) + " preserved" );
      if( found < 0 )
         continue;
      const std::string at = label + ": stable ID " + std::to_string( id ) + " ";
      expect( failures, set->getPoints().getElement( found ) == old.point, at + "position" );
      expect( failures, variables.rho.getArray().getElement( found ) == old.density, at + "rho" );
      expect( failures, variables.drho.getArray().getElement( found ) == old.densityRate, at + "drho" );
      expect( failures, variables.p.getArray().getElement( found ) == old.pressure, at + "pressure" );
      expect( failures, variables.v.getArray().getElement( found ) == old.velocity, at + "velocity" );
      expect( failures, variables.a.getArray().getElement( found ) == old.acceleration, at + "acceleration" );
      expect( failures, variables.gamma.getArray().getElement( found ) == old.gamma, at + "gamma" );
      expect( failures, variables.marker.getArray().getElement( found ) == old.marker, at + "marker" );
      expect( failures, history.rho_old.getElement( found ) == old.oldDensity, at + "rho_old" );
      expect( failures, history.v_old.getElement( found ) == old.oldVelocity, at + "v_old" );
   }
   expect( failures, variables.highestReferentialIdx == expectedHighestReferentialIdx,
           label + ": highest referential ID has the expected next-ID value" );
}

template< typename SetPointer >
void
checkAbsent( int& failures, const SetPointer& set, const int removedId, const std::string& label )
{
   bool found = false;
   for( int i = 0; i < set->getNumberOfParticles(); ++i )
      found |= set->getVariables()->referentialIdx.getElement( i ) == removedId;
   expect( failures, ! found, label + ": marked removal compacted" );
}

template< typename Pointer >
void
checkScratchGrowth( int& failures, const Pointer& buffer, const int capacity, const int zoneSize )
{
   expect( failures, buffer->particlesToFluid.getSize() == capacity, "particlesToFluid follows capacity" );
   expect( failures, buffer->particlesToRemove.getSize() == capacity, "particlesToRemove follows capacity" );
   expect( failures, buffer->retypeMarker.getSize() == capacity, "retypeMarker follows capacity" );
   expect( failures, buffer->particlesToFluid.getElement( 0 ) == 71, "particlesToFluid resize preserves contents" );
   expect( failures, buffer->particlesToRemove.getElement( 0 ) == 72, "particlesToRemove resize preserves contents" );
   expect( failures, buffer->retypeMarker.getElement( 0 ) == 73, "retypeMarker resize preserves contents" );
   expect( failures, buffer->particlesToBuffer.getSize() == zoneSize, "capture scratch remains zone-sized" );
   expect( failures, buffer->massNodes.points.getSize() == buffer->massNodes.numberOfMassNodes,
           "mass nodes remain independent" );
}

}  // namespace MultiresolutionBufferCapacityTest
