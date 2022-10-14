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

  using GlobalIndexType = typename Variables::SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;
  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using PointType = typename Particles::PointArrayType;
  using PointArrayType = typename Particles::PointArrayType;

  WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref ) : vars( size ), points( points_ref ), particles( particles_ref ) {};

  /*
  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j )
  {
    const PointType dr = points[ i ] - points[ j ];
    const PointType dv = vars.v[ i ] - vars.v[ j ];

    const RealType drs = l2Norm( dr );
    const RealType F = SPHKernel::F( drs, vars.h );
    const PointType gradW = dr*F;

    const RealType drho = ( dv, gradW )*vars.m;

    const RealType p_term = ( vars.p[ i ] + vars.p[ j ] ) / (vars.rho [ i ] * vars.rho[ j ]);
    const PointType a = p_term * gradW * vars.m;

    return { drs.x, a[ 0 ], a[ 1 ] };
  }
  */

  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j )
  {

  }

  template< typename SPHKernel >
  InteractionResultType PerformParticleInteractionBF( GlobalIndexType i, GlobalIndexType j )
  {

  }


  Variables vars;
  //PointArrayTypeView points;
  PointArrayType& points;

  Particles& particles;


};

} // SPH
} // ParticleSystem
} // TNL
