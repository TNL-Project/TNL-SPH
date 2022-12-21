#include "Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

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

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariables()
{
   auto view_particleCellIndices = particles->getParticleCellIndices().getView();
   auto view_points = particles->getPoints().getView();
   auto view_type = type.getView();
   auto view_rho = rho.getView();
   auto view_v = v.getView();
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
         swap( view_v[ i ], view_v[ j ] );
         swap( view_rhoO[ i ], view_rhoO[ j ] );
         swap( view_vO[ i ], view_vO[ j ] );
         } );
}

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

/* TEMP INTEGRATION, MOVE OUT */
template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IntegrateVerlet( RealType dt )
{
   auto rho_view = this->getRho().getView();
   auto v_view = this->getVel().getView();
   auto r_view = this->particles->getPoints().getView();

   auto rhoO_view = this->rhoO.getView();
   auto vO_view = this->vO.getView();

   const auto drho_view = this->getDrho().getView();
   const auto a_view = this->getAcc().getView();

   RealType dtdt05 = 0.5 * dt * dt;
   RealType dt2 = 2 * dt;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
      vO_view[ i ] += a_view[ i ] * dt2;
      rhoO_view[ i ] += drho_view[ i ] * dt2;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, this->particles->getNumberOfParticles(), init );

   std::swap( this->v, this->vO ); //swap pointers instead
   std::swap( this->rho, this->rhoO ); //swap pointers instead
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IntegrateEuler( RealType dt )
{
   auto rho_view = this->getRho().getView();
   auto v_view = this->getVel().getView();
   auto r_view = this->particles->getPoints().getView();

   auto rhoO_view = this->rhoO.getView();
   auto vO_view = this->vO.getView();

   auto drho_view = this->getDrho().getView();
   auto a_view = this->getAcc().getView();

   RealType dtdt05 = 0.5 * dt * dt;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
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

