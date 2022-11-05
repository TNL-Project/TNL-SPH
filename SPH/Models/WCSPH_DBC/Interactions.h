#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "../../SPHTraits.h"
#include "Variables.h"

/**
 * Modules used as default.
 **/
#include "../EquationOfState.h"
#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
  using DeviceType = typename Particles::Device; //?

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType; //Particles::?
  using RealType = typename SPHFluidTraitsType::RealType; //Particles::?

  using PointType = typename Particles::PointType;
  using PointArrayType = typename Particles::PointArrayType;

  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using Integrator = VerletIntegrator< Particles, SPHFluidConfig, Variables >;

  /**
   * Constructor.
   **/
  WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref )
  : vars( size ), points( points_ref ), particles( particles_ref ), integrator( size, vars, points ) {};

  /**
   * Process one general parcile with index i.
   */
  __cuda_callable__
  void
  ProcessOneParticle( GlobalIndexType index_i );

  /**
   * Process one fluid particle with index i.
   */
  __cuda_callable__
  void
  ProcessOneFluidParticle( GlobalIndexType index_i  );

  /**
   * Process one boundary particle with index i.
   */
  __cuda_callable__
  void
  ProcessOneBoundaryParticle( GlobalIndexType index_i );

  /**
   * Evaluate interaction between fluid particle i and fluid particle j.
   */
  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j );

  /**
   * Evaluate interaction between fluid particle i and boundary particle j.
   */
  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j );

  /**
   * Evaluate interaction between boundary particle i and fluid particle j.
   */
  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionBF( GlobalIndexType i, GlobalIndexType j );

  void sortParticlesAndVariables();

  template< typename EquationOfState = TaitWeaklyCompressibleEOS >
  void
  ComputePressureFromDensity();

//protected:

  Variables vars;
  PointArrayType& points;

  //PointArrayTypeView pointss; //mby like this?
  //PointArrayType_ptr points; //mby like this?

  Particles& particles;

  Integrator integrator; //temp

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"
