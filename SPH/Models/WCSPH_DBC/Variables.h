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

  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  using ParticleTypeArrayType = typename SPHFluidTraitsType::ParticleTypeArrayType;
  using InteractionResultTypeArrayType = typename SPHFluidTraitsType::InteractionResultTypeArray;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;

  SPHFluidVariables( GlobalIndexType size )
  : rho( size ), p( size ), v( size ), type( size ), DrhoDv( size ) {}

  /* Variables - Fields */

  ParticleTypeArrayType type;

  ScalarArrayType rho;
  ScalarArrayType p;
  VectorArrayType v;

  InteractionResultTypeArrayType DrhoDv;

  /* Integrator */

  /* Variables - constans */
  RealType h, m, speedOfSound, coefB, rho0;


};

template< typename SPHFluidConfig >
class SPHFluidConstantVariables
{

};

} // SPH
} // ParticleSystem
} // TNL
