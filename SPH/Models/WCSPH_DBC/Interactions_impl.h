#include "Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

	/* NEW */

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ParticleTypeArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getParticleType() const
{
	return this->type;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ParticleTypeArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getParticleType()
{
	return this->type;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getRho() const
{
	return this->rho;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getRho()
{
	return this->rho;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getDrho() const
{
	return this->drho;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getDrho()
{
	return this->drho;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getPress() const
{
	return this->p;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ScalarArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getPress()
{
	return this->p;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::VectorArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getVel() const
{
	return this->v;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::VectorArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getVel()
{
	return this->v;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::VectorArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getAcc() const
{
	return this->a;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::VectorArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getAcc()
{
	return this->a;
}

	/* OLD */

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariables()
{
   auto view_particleCellIndices = particles->getParticleCellIndices().getView();
   auto view_points = particles->getPoints().getView();

   auto view_type = type.getView();
   auto view_rho = rho.getView();
   auto view_drho = drho.getView();
   auto view_p = p.getView();
   auto view_v = v.getView();
   auto view_a = a.getView();

	 //temp, this needs to be done in better way - sort integrator arrays
	 auto view_rhoO = rhoO.getView();
	 auto view_vO = vO.getView();

   Algorithms::sort< DeviceType, GlobalIndexType >(
       0, particles->getNumberOfParticles(),
       [=] __cuda_callable__ ( int i, int j ) -> bool {
         return view_particleCellIndices[ i ] <= view_particleCellIndices[ j ]; },
       [=] __cuda_callable__ ( int i, int j ) mutable {
         swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
         swap( view_points[ i ], view_points[ j ] );
         swap( view_type[ i ], view_type[ j ] );
         swap( view_rho[ i ], view_rho[ j ] );
         swap( view_drho[ i ], view_drho[ j ] );
         swap( view_p[ i ], view_p[ j ] );
         swap( view_v[ i ], view_v[ j ] );
         swap( view_a[ i ], view_a[ j ] );
	 			//temp, this needs to be done in better way - sort integrator arrays
         swap( view_rhoO[ i ], view_rhoO[ j ] );
         swap( view_vO[ i ], view_vO[ j ] );
         } );
}

/*
template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ProcessOneParticle( GlobalIndexType index_i )
{
  if( type[ index_i ] == 0 )
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

    if( type[ index_j ] == 0 )
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
*/

/*
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
  const PointType a = ( -1.0 ) * ( p_term + visco )* gradW * vars.m;
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
  const PointType a = ( -1.0 ) * ( p_term + visco ) * gradW * vars.m;

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
*/

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename EquationOfState >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ComputePressureFromDensity()
{
  auto view_rho = this->getRho().getView();
  auto view_p = this->getPress().getView();

  auto init = [=] __cuda_callable__ ( int i ) mutable
  {
     view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ] );
  };
  Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), init );
}

/* TEMP INTEGRATION, REMOVE LATER */
template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IntegrateVerlet( RealType dt )
{
    auto rhoO_view = this->rhoO.getView();
    auto vO_view = this->vO.getView();

    auto rho_view = this->getRho().getView();
    auto v_view = this->getVel().getView();
    auto r_view = this->particles->getPoints().getView();

    auto drho_view = this->getDrho().getView();
    auto a_view = this->getAcc().getView();

    RealType dtdt05 = 0.5 * dt * dt;
    RealType dt2 = 2 * dt;

    auto init = [=] __cuda_callable__ ( int i ) mutable
    {

			/*
       r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
       rho_view[ i ] = rhoOO_view[ i ] + drho_view[ i ] * dt2;
       v_view[ i ] = vOO_view[ i ] + a_view[ i ] * dt2;

       vOO_view[ i ] = vO_view[ i ];
       vO_view[ i ] = v_view[ i ];

       rhoOO_view[ i ] = rhoO_view[ i ];
       rhoO_view[ i ] = rho_view[ i ];
			 */

       r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
       vO_view[ i ] += a_view[ i ] * dt2;

       rhoO_view[ i ] += drho_view[ i ] * dt2;

    };
    Algorithms::ParallelFor< DeviceType >::exec( 0, this->particles->getNumberOfParticles(), init );

		//auto velocity_auxData = getVel().getData();
		//v_view = vO_view;
		//vO_view = velocity_auxView;

		//auto rho_auxView = getRho().getView();
		//rho_view = rhoO_view;
		//rhoO_view = rho_auxView;
		std::swap( this->v, this->vO );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IntegrateEuler( RealType dt )
{

  auto rhoO_view = this->rhoO.getView();
  auto vO_view = this->vO.getView();

  auto rho_view = this->getRho().getView();
  auto v_view = this->getVel().getView();
  auto r_view = this->particles->getPoints().getView();

  auto drho_view = this->getDrho().getView();
  auto a_view = this->getAcc().getView();

  RealType dtdt05 = 0.5 * dt * dt;

  auto init = [=] __cuda_callable__ ( int i ) mutable
  {

		/*
			v_view[ i ] += a_view[ i ] * dt;
     	r_view[ i ] += vO_view[ i ] * dt + a_view[ i ] * dtdt05;
			rho_view[ i ] += drho_view[ i ] * dt;

		 	vOO_view[ i ] = vO_view[ i ];
     	vO_view[ i ] = v_view[ i ];

     	rhoOO_view[ i ] = rhoO_view[ i ];
     	rhoO_view[ i ] = rho_view[ i ];
			*/

      r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
			vO_view[ i ] = v_view[ i ];
			v_view[ i ] += a_view[ i ] * dt;
			rhoO_view[ i ] = rho_view[ i ];
			rho_view[ i ] += drho_view[ i ] * dt;

  };

  Algorithms::ParallelFor< DeviceType >::exec( 0, this->particles->getNumberOfParticles(), init );
}


} // SPH
} // ParticleSystem
} // TNL

