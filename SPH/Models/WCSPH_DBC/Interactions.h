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
//#include "Integrator.h"

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

  using ScalarType = typename SPHFluidTraitsType::ScalarType;
  using VectorType = typename SPHFluidTraitsType::VectorType;
  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using DiffusiveTerm = MolteniDiffusiveTerm< SPHFluidConfig >; //-> template
  using ViscousTerm = ArtificialViscosity< SPHFluidConfig >; //-> template

	using ParticlePointer = typename Pointers::SharedPointer< Particles, DeviceType >;

	/* VARIABLES FIELDS */
  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;

  /**
   * Constructor.
   */
  WCSPH_DBC( GlobalIndexType size, ParticlePointer& particles )
  : type( size ), rho( size ), drho( size ), p( size ), v( size ), a( size ) , rhoO( size ), vO( size ), particles( particles )
	{
		vO = 0.; //remove
		rhoO = 1000.; //remove
	}

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

  /**
   * Sort particles and all variables based on particle cell index.
	 * TODO: Move this on the side of nbsearch/particles.
   */
  void sortParticlesAndVariables();

  /**
   * Compute pressure from density.
	 * TODO: Move out.
   */
  template< typename EquationOfState = TaitWeaklyCompressibleEOS< SPHFluidConfig > >
  void
  ComputePressureFromDensity();

	/* TEMP INTEGRATORS, MOVE OUT */
	void
	IntegrateVerlet( RealType dt );

	void
	IntegrateEuler( RealType dt );

//protected:

  /* Variables - Fields */
  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;

	/* TEMP INTEGRATORS, MOVE OUT */
  ScalarArrayType rhoO;
  VectorArrayType vO;

	/* Constants */
  RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

  ParticlePointer particles;
  //Integrator integrator; //temp

};

} // SPH
} // ParticleSystem
} // TNL

#include "Interactions_impl.h"

