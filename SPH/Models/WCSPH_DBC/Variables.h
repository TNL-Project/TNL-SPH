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
  using InteractionResultTypeArrayType = typename SPHFluidTraitsType::InteractionResultTypeArray;


  SPHFluidVariables( GlobalIndexType size )
  : rho( size ), p( size ), v( size ), type( size ), DrhoDv( size ) {}

  /* Variables - Fields */
  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType p;
  VectorArrayType v;

  InteractionResultTypeArrayType DrhoDv;

  /* Variables - constans */
  RealType h, m, speedOfSound, coefB, rho0;

};

} // SPH
} // ParticleSystem
} // TNL
