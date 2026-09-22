#pragma once

#include "SolverMultiSetBlockMultiresolution.h"
#include "../MultiresolutionRectangleBufferLocalTimeStepping.h"

namespace TNL {
namespace SPH {

/**
 * \brief Local-timestepping variant of the block multiresolution solver.
 *
 * Levels advance with per-level time steps derived from the refinement-factor ratios
 * (finest level uses the configured time step). The fine side of an interface gets a
 * two-layer buffer band (SolverMultiSetBlockMultiresolutionLocalTimestepping uses
 * MultiresolutionBoundaryLTS patches), the coarse side keeps one layer; the global
 * clock advances with the coarse step and fine substeps are driven from the case
 * (see Documentation/multiresolution-local-timestepping.md for the schedule).
 */
template< typename Model >
class SolverMultiSetBlockMultiresolutionLocalTimestepping : public SolverMultiSetBlockMultiresolution< Model >
{
public:

   using BaseType = SolverMultiSetBlockMultiresolution< Model >;

   using typename BaseType::DeviceType;
   using typename BaseType::ModelType;
   using typename BaseType::ModelParams;
   using typename BaseType::ParticlesType;
   using typename BaseType::SPHConfig;
   using typename BaseType::GlobalIndexType;
   using typename BaseType::RealType;
   using typename BaseType::CoordinatesType;
   using typename BaseType::VectorType;
   using typename BaseType::TimeStepping;

   using typename BaseType::FluidVariables;
   using typename BaseType::FluidPointer;
   using typename BaseType::BoundaryPointer;
   using typename BaseType::OpenBoundaryConfigType;

   using MultiresolutionBoundaryLTS = MultiresolutionBoundaryLTS<
      ParticlesType, SPHConfig, FluidVariables, typename BaseType::IntegrationSchemeVariablesType, OpenBoundaryConfigType, ModelParams >;
   using MultiresolutionBoundaryLTSPointer = Pointers::SharedPointer< MultiresolutionBoundaryLTS, DeviceType >;

   SolverMultiSetBlockMultiresolutionLocalTimestepping( std::ostream& out = std::cout ) : BaseType( out ) {};

   void
   init( int argc, char* argv[] );

   void
   initializeBlockBasedMultiResolutionSimulation();

   void
   initParticleSets();

   void
   initMultiResolutionBoundaryPatches();

   float
   getLevelRefinementFactor( int idx ) const;

   // two layers on the fine side of an interface, one on the coarse side
   int
   getLevelOverlapWidth( int idx ) const;

    /* Per-level step primitives driven from the case exec loop. coarse/fine indices
      follow the subdomain order; levelSubsteps[i] gives the number of steps level i
      takes per global (coarsest) step, levelTimeStepping[i] its time-stepping state. */
    void
    interactCoarseLevel( int coarseIdx );

    void
    interactFineSubstep( int fineIdx, bool computeActiveLayer );

    void
    integrateLevel( int idx );

    /* Buffer operations restricted to the interface whose finer side is fineLevelIdx
      (a chain of N zones creates N-1 interfaces; the band of iface(i,i+1) follows the
      clock of level i+1). */
    void
    advectSubstepBuffers( RealType dt, bool skipLayer1, int fineLevelIdx );

    void
    integrateBuffersLayer1( RealType dt, bool useEulerRestart, int fineLevelIdx );

   void
   accumulateInterfaceMassesFromNeighbor( int neighborIdx, RealType dt );

   void
   refreshSubstepSearches( int fineIdx );

   void
   syncMultiresolutionUpdate();

   void
   save( bool writeParticleCellIndex = false );

   void
   writeProlog( bool writeSystemInformation = true ) noexcept;

   void
   writeInfo() noexcept;

   std::vector< MultiresolutionBoundaryLTSPointer > multiresolutionBoundaryPatchesLTS;

   std::vector< TimeStepping > levelTimeStepping;
   std::vector< int > levelSubsteps;

};

} // namespace SPH
} // namespace TNL

#include "SolverMultiSetBlockMultiresolutionLocalTimestepping.hpp"
