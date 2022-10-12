#pragma once

#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
class SPHFluidVariables
{

  public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;
  using InteractionResultTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;

  SPHFluidVariables(GlobalIndexType size)
  : rho(size), p(size), v(size), type(size), DrhoDv(size) {}

  /* Variables - Fields */

  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType p;
  VectorArrayType v;

  InteractionResultTypeArrayType DrhoDv;

  /* Variables - constans */
  RealType h;


};

template< typename SPHFluidConfig >
class SPHFluidConstantVariables
{

};

} // SPH
} // ParticleSystem
} // TNL
