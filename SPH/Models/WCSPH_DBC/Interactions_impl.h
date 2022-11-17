#include "Interactions.h"

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
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ProcessOneParticle( GlobalIndexType index_i )
{
  if( vars.type[ index_i ] == 0 )
    ProcessOneFluidParticle( index_i );
  else
    ProcessOneBoundaryParticle( index_i );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ProcessOneFluidParticle( GlobalIndexType index_i )
{
  auto fetch = [=] __cuda_callable__ ( int i ) -> InteractionResultType
  {
    GlobalIndexType index_j = particles.getNeighbor( index_i, i );

    if( vars.type[ index_j ] == 0 )
      return PerformParticleInteractionFF< WendlandKernel, DiffusiveTerm, ViscousTerm >( index_i , index_j );
    else
      return PerformParticleInteractionFB< WendlandKernel, DiffusiveTerm, ViscousTerm >( index_i , index_j );
  };

  auto reduction = [] __cuda_callable__ ( const InteractionResultType& a, const InteractionResultType& b )
  {
     return a + b;
  };

  vars.DrhoDv[ index_i ] = Algorithms::reduce< DeviceType, int, InteractionResultType >( 0, particles.getNeighborsCount( index_i ), fetch, reduction, { 0., 0., -9.81 } );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ProcessOneBoundaryParticle( GlobalIndexType index_i )
{
  auto fetch = [=] __cuda_callable__ ( int i ) -> InteractionResultType
  {
    GlobalIndexType index_j = particles.getNeighbor( index_i, i );

    if( vars.type[ index_j ] == 0 )
      return PerformParticleInteractionBF< WendlandKernel >( index_i , index_j );
    else
      return {0., 0., 0.};
  };

  auto reduction = [] __cuda_callable__ ( const InteractionResultType& a, const InteractionResultType& b )
  {
     return a + b;
  };

  vars.DrhoDv[ index_i ] = Algorithms::reduce< DeviceType, int, InteractionResultType >( 0, particles.getNeighborsCount( index_i ), fetch, reduction, { 0., 0., 0. } );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
__cuda_callable__
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::InteractionResultType
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j )
{
  const PointType dr = this->points[ i ] - this->points[ j ];
  const PointType dv = vars.v[ i ] - vars.v[ j ];

  const RealType drs = l2Norm( dr );
  const RealType F = SPHKernelFunction::F( drs, vars.h );
  const PointType gradW = dr * F;

	const RealType psi = DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );
	const RealType diffTerm =  psi * ( dr, gradW ) * vars.m / vars.rho[ j ];
  const RealType drho = ( dv, gradW ) * vars.m - diffTerm;
  //const RealType drho = ( dv, gradW ) * vars.m + DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );

  const RealType p_term = ( vars.p[ i ] + vars.p[ j ] ) / ( vars.rho [ i ] * vars.rho[ j ] );
  const RealType visco =  ViscousTerm::Pi( vars.rho[ i ], vars.rho[ j ], drs, ( dr, dv ) );
  const PointType a = ( -1 ) * ( p_term + visco )* gradW * vars.m;
  //const PointType a = ( p_term + visco )* gradW * vars.m;

  return { drho, a[ 0 ], a[ 1 ] };
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
__cuda_callable__
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::InteractionResultType
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j )
{
  const PointType dr = this->points[ i ] - this->points[ j ];
  const PointType dv = vars.v[ i ] - vars.v[ j ];

  const RealType drs = l2Norm( dr );
  const RealType F = SPHKernelFunction::F( drs, vars.h );
  const PointType gradW = dr*F;

	const RealType psi = DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );
	const RealType diffTerm =  psi * ( dr, gradW ) * vars.m / vars.rho[ j ];
  const RealType drho = ( dv, gradW ) * vars.m - diffTerm;
  //const RealType drho = ( dv, gradW ) * vars.m + DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );

  const RealType p_term = ( vars.p[ i ] + vars.p[ j ] ) / ( vars.rho [ i ] * vars.rho[ j ] );
  const RealType visco =  ViscousTerm::Pi( vars.rho[ i ], vars.rho[ j ], drs, ( dr, dv ) );
  const PointType a = ( -1 ) * ( p_term + visco ) * gradW * vars.m;

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

	const RealType psi = DiffusiveTerm::Psi( vars.rho[ i ], vars.rho[ j ], drs );
	const RealType diffTerm =  psi * ( dr, gradW ) * vars.m / vars.rho[ j ];
  const RealType drho = ( dv, gradW ) * vars.m - diffTerm;
  //const RealType drho = ( dv, gradW )*vars.m;
  const PointType a = { 0., 0. };

  return { drho, a[ 0 ], a[ 1 ] };
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename EquationOfState >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ComputePressureFromDensity()
{
  auto view_rho = vars.rho.getView();
  auto view_p = vars.p.getView();

  auto init = [=] __cuda_callable__ ( int i ) mutable
  {
     view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ] );
  };
  Algorithms::ParallelFor< DeviceType >::exec( 0, particles.getNumberOfParticles(), init );
}


} // SPH
} // ParticleSystem
} // TNL

