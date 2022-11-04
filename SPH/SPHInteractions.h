#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "SPHFluidTraits.h"
#include "SPHFluidVariables.h"

#include "SPHFluidVariables.h"
#include "../Particles/ParticlesTraits.h"

#include "SPHequationOfState.h"

#include "SPHKernels.h"

#include "SPHInteractionsIntegrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename Variables::SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;
  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using DeviceType = typename Particles::Device;
  using PointType = typename Particles::PointType;
  using PointArrayType = typename Particles::PointArrayType;

  using Integrator = VerletIntegrator< Particles, SPHFluidConfig, Variables >;

  //using PointArrayTypeView = Containers::ArrayView< PointType >; /*temp-test*/

  WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref )
  : vars( size ), points( points_ref ), particles( particles_ref ), integrator( size, vars, points ) {};
  //WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref,  PointArrayTypeView points_view ) : vars( size ), points( points_ref ), particles( particles_ref ),
  //points_view( points_view ) {};

  __cuda_callable__
  void
  ProcessOneParticle( GlobalIndexType index_i );

  __cuda_callable__
  void
  ProcessOneFluidParticle( GlobalIndexType index_i  );

  __cuda_callable__
  void
  ProcessOneBoundaryParticle( GlobalIndexType index_i );

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j );

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j );

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType
  PerformParticleInteractionBF( GlobalIndexType i, GlobalIndexType j );

  void sortParticlesAndVariables();

  template< typename EquationOfState = TaitWeaklyCompressibleEOS >
  void
  ComputePressureFromDensity();

  Variables vars;
  //PointArrayTypeView pointss;
  PointArrayType& points;
  //PointArrayType_ptr points;

  Particles& particles;

  //DRAFT Integrator
  Integrator integrator;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPHInteractions_impl.h"
