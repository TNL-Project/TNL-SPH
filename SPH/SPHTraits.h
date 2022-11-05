#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Array.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
class SPHFluidTraits
{
  public:
  static constexpr int spaceDimension = SPHFluidConfig::spaceDimension;

  using DeviceType = typename SPHFluidConfig::DeviceType;
  using GlobalIndexType = typename SPHFluidConfig::GlobalIndexType;
  using LocalIndexType = typename SPHFluidConfig::LocalIndexType;
  using CellIndexType = typename SPHFluidConfig::CellIndexType;
  using RealType = typename SPHFluidConfig::RealType;

  /* particle related */
  using ParticleType = unsigned short int; //kuk
  using ParticleTypeArrayType = Containers::Array< ParticleType, DeviceType, GlobalIndexType >;

  using ScalarType = RealType;
  using ScalarArrayType = Containers::Array< ScalarType, DeviceType, GlobalIndexType >;

  using VectorType = Containers::StaticVector< spaceDimension, RealType >;
  using VectorArrayType = Containers::Array< VectorType, DeviceType, GlobalIndexType >;

  /* ? */
  using InteractionResultType = Containers::StaticVector< spaceDimension + 1, RealType >;
  using InteractionResultTypeArray = Containers::Array< InteractionResultType, DeviceType, GlobalIndexType >;
};

} // SPH
} // ParticleSystem
} // TNL
