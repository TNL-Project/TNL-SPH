#pragma once

#include "../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
class SPHFluidVariables
{
  public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;

  SPHFluidVariables( GlobalIndexType size )
  : rho( size ), p( size ), v( size ), type( size ) {}

  /* Variables - Fields */
  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType p;
  VectorArrayType v;

  /* Variables - constans */
  RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

};

} // SPH
} // ParticleSystem
} // TNL

