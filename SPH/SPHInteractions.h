#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "SPHFluidTraits.h"
#include "SPHFluidVariables.h"

#include "SPHFluidVariables.h"
#include "../Particles/ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using RealType = typename SPHFluidTraitsType::RealType;

  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;
  using GlobalIndexType = typename Variables::SPHFluidTraitsType::GlobalIndexType;

  using PointArrayType = typename Particles::PointArrayType;
  //using PointArrayTypeView = typename Containers::ArrayView< PointArrayType >;

  WCSPH_DBC(GlobalIndexType size, PointArrayType& points_ref) : vars(size), points(points_ref) {};

  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionFF(GlobalIndexType i, GlobalIndexType j)
  {

  }

  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionFB(GlobalIndexType i, GlobalIndexType j)
  {

  }

  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionBF(GlobalIndexType i, GlobalIndexType j)
  {

  }


  Variables vars;
  //PointArrayTypeView points;
  PointArrayType& points;

};

} // SPH
} // ParticleSystem
} // TNL
