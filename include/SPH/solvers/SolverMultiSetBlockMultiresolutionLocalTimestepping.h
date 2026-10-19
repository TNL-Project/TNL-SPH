#pragma once

#include "SolverMultiSetBlockMultiresolution.h"
#include "../MultiresolutionRectangleBufferLocalTimeStepping.h"
#include "../TimeStep.h"

namespace TNL {
namespace SPH {

/**
 * \brief Local-timestepping variant of the block multiresolution solver.
 *
 * Levels advance with per-level time steps derived from the refinement-factor ratios
 * (finest level uses the configured time step). The fine side of an interface gets a
 * two-layer buffer band (MultiresolutionBoundaryLTS patches), the coarse side keeps
 * one layer. The global clock advances with the coarse step and fine substeps are
 * driven from the case exec loop as a Berger-Colella recursion over the levels with
 * the coarsest level included (see the damBreak2D multiresolution local-timestepping
 * example for the schedule). Buffer bands drift ballistically during the substeps of
 * the level that owns them and their state is re-grounded from the coarser neighbor
 * inside the window (refreshBuffersFromNeighbor): exactly at level-clock coincidences
 * and as a theta-blend at level-clock midpoints. The full buffer lifecycle runs only
 * at sync points, when all level clocks coincide again.
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

    /* Per-level step primitives driven from the case exec loop. Level indices follow
      the subdomain order; getNumberOfSubsteps(i) gives the number of substeps level i
      takes per global (coarsest) step, levelTimeStepping[i] its time-stepping state. */
    void
    interact( int subdomainIdx );

    void
    integrateLevel( int idx );

    int
    getNumberOfSubsteps( int level ) const;

    /* Re-derive the window time step and the per-level steps from the current CFL
      conditions, one governing step per global window. Every level proposes its own
      CFL-limited step (with that level's smoothing length), each is converted to the
      implied finest-level step, the minimum becomes the window's finest step, and all
      level steps are re-snapped from it with the refinement-factor ratios so that
      levelSubsteps(i) * dt_i == dt_window holds exactly - the level-clock coincidence
      the R1/R2 band refreshes and the sync-point buffer lifecycle rely on. The
      substep counts themselves stay fixed. No-op for ConstantTimeStep (the configured
      initial steps are kept). Call once per synchronized global step, before the
      Berger-Colella recursion. */
    void
    updateLevelTimeSteps();

    /* Ballistic advection of the buffer bands restricted to the interface whose finer
      side is fineLevelIdx (a chain of N zones creates N-1 interfaces; the band of
      iface(i,i+1) follows the clock of level i+1). Both layers drift. */
    void
    advectSubstepBuffers( RealType dt, int fineLevelIdx );

    /* Band state refresh of the interface clocked by fineLevelIdx: interpolates rho
      and v from the (already advanced) coarser neighbor and blends them with the
      sync-time anchors kept in rho_old/v_old,
         q <- theta * q_interp( now ) + ( 1 - theta ) * q_anchor( sync ).
      theta = 1 is the exact coincidence refresh used when the two level clocks meet
      inside the window. Rebuilds the neighbor's search structure first - its particles
      moved since the last sync. */
    void
    refreshBuffersFromNeighbor( RealType theta, int fineLevelIdx );

    /* Re-stamp the band anchors (rho_old/v_old) to the current band state. Called
      after an exact coincidence refresh so the following midpoint blends stay
      time-bracketed (the anchor becomes the band state at the coincidence time).
      Restricted to the fine-side band of the interface clocked by fineLevelIdx. */
    void
    stampBuffersAnchor( int fineLevelIdx );

    void
    accumulateInterfaceMassesFromNeighbor( int neighborIdx, RealType dt );

    void
    refreshSubstepSearches( int fineIdx );

    /* Berger-Colella recursion over the refinement levels with the coarsest level
      included: drives the whole local-timestepping schedule of one global window and
      invokes subdomainSubstep( level, substep, globalSubstep, dt ) for every level
      substep (the case supplies the per-substep physics). Level L takes
      getNumberOfSubsteps(L) substeps of dt_L per window and recurses into L+1 after
      each of them, so the finer level always advances within the just-integrated
      window of its parent. Substeps are indexed from zero; the window-global index
      g = base + substep is threaded through the recursion. Searches are refreshed
      before every substep except the genuine first of the window (g == 0, covered by
      the sync-time performNeighborSearch). The state of the band of iface(level,
      level+1) is re-grounded inside the window: exactly (theta = 1, anchors re-stamped)
      after every non-final parent substep where the two clocks coincide, and blended
      with theta = modelParams.ltsTheta at each level's own midpoint substep. */
    template< typename SubstepUpdate >
    void
    recursiveSubdomainsLoop( SubstepUpdate&& subdomainSubstep );

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

