#include "SPHInteractions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariables()
{
   auto view_particleCellIndices = particles.getParticleCellIndices().getView();
   auto view_points = points.getView();

   auto view_type = vars.type.getView();
   auto view_rho = vars.rho.getView();
   auto view_p = vars.p.getView();
   auto view_v = vars.v.getView();
   auto view_DrhoDv = vars.DrhoDv.getView();

   Algorithms::sort< DeviceType, GlobalIndexType >(
       0, particles.getNumberOfParticles(),
       [=] __cuda_callable__ ( int i, int j ) -> bool {
         return view_particleCellIndices[ i ] < view_particleCellIndices[ j ]; },
       [=] __cuda_callable__ ( int i, int j ) mutable {
         swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
         swap( view_points[ i ], view_points[ j ] );
         swap( view_type[ i ], view_type[ j ] );
         swap( view_rho[ i ], view_rho[ j ] );
         swap( view_p[ i ], view_p[ j ] );
         swap( view_v[ i ], view_v[ j ] );
         swap( view_DrhoDv[ i ], view_DrhoDv[ j ] );
         } );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename SPHKernelFunction >
__cuda_callable__
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::InteractionResultType
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j )
{
  const PointType dr = this->points[ i ] - this->points[ j ];
  const PointType dv = vars.v[ i ] - vars.v[ j ];

  const RealType drs = l2Norm( dr );
  const RealType F = SPHKernelFunction::F( drs, vars.h );
  const PointType gradW = dr*F;

  const RealType drho = ( dv, gradW )*vars.m;

  const RealType p_term = ( vars.p[ i ] + vars.p[ j ] ) / (vars.rho [ i ] * vars.rho[ j ]);
  const PointType a = p_term * gradW * vars.m;

  return { drho, a[ 0 ], a[ 1 ] };
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename SPHKernelFunction >
__cuda_callable__
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::InteractionResultType
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j )
{
  const PointType dr = this->points[ i ] - this->points[ j ];
  const PointType dv = vars.v[ i ] - vars.v[ j ];

  const RealType drs = l2Norm( dr );
  const RealType F = SPHKernelFunction::F( drs, vars.h );
  const PointType gradW = dr*F;

  const RealType drho = ( dv, gradW )*vars.m;

  const RealType p_term = ( vars.p[ i ] + vars.p[ j ] ) / (vars.rho [ i ] * vars.rho[ j ]);
  const PointType a = p_term * gradW * vars.m;

  return { drho, a[ 0 ], a[ 1 ] };
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename SPHKernelFunction >
__cuda_callable__
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::InteractionResultType
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::PerformParticleInteractionBF( GlobalIndexType i, GlobalIndexType j )
{
  const PointType dr = this->points[ i ] - this->points[ j ];
  const PointType dv = vars.v[ i ] - vars.v[ j ];

  const RealType drs = l2Norm( dr );
  const RealType F = SPHKernelFunction::F( drs, vars.h );
  const PointType gradW = dr*F;

  const RealType drho = ( dv, gradW )*vars.m;
  const PointType a = { 0., 0. };

  return { drho, a[ 0 ], a[ 1 ] };
}


} // SPH
} // ParticleSystem
} // TNL
