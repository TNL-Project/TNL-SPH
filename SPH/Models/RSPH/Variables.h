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
  : rho( size ), drho ( size ), p( size ), v( size ), a( size ), n( size ) {}

  /* Variables - Fields */
  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;
  VectorArrayType n;

  /* Variables - constans */
  //RealType h, m, speedOfSound, coefB, rho0, delta, alpha;

};

//temp
template< typename SPHFluidConfig, typename PointArrayType >
class SWAPFluidVariables
{
  public:
  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType; //REDO

  SWAPFluidVariables( GlobalIndexType size )
  : indicesMap( size ), rho_swap( size ), v_swap( size ), points_swap( size ) {}

  /* Variables - Fields */
  //TEMP
  IndexArrayType indicesMap;

  /* Copy of arrays needed for resort */
  ScalarArrayType rho_swap;
  VectorArrayType v_swap;
  PointArrayType points_swap;

};

template< typename SPHFluidConfig, typename PointArrayType >
class SWAPBoundaryVariables
{
  public:
  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType; //REDO

  SWAPBoundaryVariables( GlobalIndexType size )
  : indicesMap( size ) {}

  /* Variables - Fields */
  IndexArrayType indicesMap;


};

} // SPH
} // ParticleSystem
} // TNL

