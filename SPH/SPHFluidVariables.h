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
  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;

  SPHFluidVariables(GlobalIndexType size)
  : rho(size), p(size), a(size), v(size) {}

  ScalarArrayType rho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;

};

} // SPH
} // ParticleSystem
} // TNL
