#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "../../SPHTraits.h"
#include "Variables.h"

/**
 * Modules used as default.
 **/
#include "../EquationOfState.h"
#include "../DiffusiveTerms.h"
#include "../VisousTerms.h"
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

  using Integrator = VerletIntegrator< Particles, SPHFluidConfig, Variables >; //-> template
  using DiffusiveTerm = MolteniDiffusiveTerm< SPHFluidConfig >; //-> template
  using ViscousTerm = ArtificialViscosity< SPHFluidConfig >; //-> template

	using ParticlePointer = typename Pointers::SharedPointer< Particles, DeviceType >;

	/* VARIABLES FIELDS */
  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;

  /**
   * Constructor.
   **/
  //IDLT: WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref )
  //IDLT: : vars( size ), points( points_ref ), particles( particles_ref ), integrator( size, vars, points ) {};
  WCSPH_DBC( GlobalIndexType size, ParticlePointer& particles )
  //: type( size ), rho( size ), drho( size ), p( size ), v( size ), a( size ) , particles( particles ), integrator( size, vars, particles ) {}; //add integrator
  : type( size ), rho( size ), drho( size ), p( size ), v( size ), a( size ) , particles( particles ) {}; //add integrator

	/* NEW */
  /**
   * Get fileds with variables.
   */
  const ParticleTypeArrayType&
  getParticleType() const;

  ParticleTypeArrayType&
  getParticleType();

  const ScalarArrayType&
  getRho() const;

  ScalarArrayType&
  getRho();

  const ScalarArrayType&
  getDrho() const;

  ScalarArrayType&
  getDrho();

  const ScalarArrayType&
  getPress() const;

  ScalarArrayType&
  getPress();

  const VectorArrayType&
  getVel() const;

  VectorArrayType&
  getVel();

  const VectorArrayType&
  getAcc() const;

  VectorArrayType&
  getAcc();

	/* OLD */

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
  template< typename SPHKernelFunction, typename DiffusiveTerm, typename VisousTerm >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j );

  /**
   * Evaluate interaction between fluid particle i and boundary particle j.
   */
  template< typename SPHKernelFunction, typename DiffusiveTerm, typename VisousTerm >
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

  template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
  void
  ComputePressureFromDensity();

//protected:

  //Variables vars;
  //IDLT: PointArrayType& points;

  /* Variables - Fields */
  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;

	/* Constants */
  RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

  //PointArrayTypeView pointss; //mby like this?
  //PointArrayType_ptr points; //mby like this?

  Particles particles;

  //Integrator integrator; //temp

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"

