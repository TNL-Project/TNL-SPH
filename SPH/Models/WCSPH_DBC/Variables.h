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

  SPHFluidVariables( GlobalIndexType size )
  : rho( size ), drho ( size ), p( size ), v( size ), a( size ) {}

  /* Variables - Fields */
  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;

  /* Variables - constans */
  //RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

};

template< typename SPHFluidConfig >
class SPHBoundaryVariables
{
  public:
  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  SPHBoundaryVariables( GlobalIndexType size )
  : rho( size ), drho ( size ), p( size ), v( size ), a( size ) {}

  /* Variables - Fields */
  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;

  /* Variables - constans */
  //RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

};

} // SPH
} // ParticleSystem
} // TNL

