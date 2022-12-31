#include "Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
const Variables&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getFluidVariables() const
{
   return FluidVariables;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
Variables&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getFluidVariables()
{
   return this->FluidVariables;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariables()
{
   auto view_particleCellIndices = particles->getParticleCellIndices().getView();
   auto view_points = particles->getPoints().getView();
   auto view_rho = this->FluidVariables.rho.getView();
   auto view_v = this->FluidVariables.v.getView();
   auto view_rhoO = this->rhoO.getView();
   auto view_vO = this->vO.getView();

   Algorithms::sort< DeviceType, GlobalIndexType >(
       0, particles->getNumberOfParticles(),
       [=] __cuda_callable__ ( int i, int j ) -> bool {
         return view_particleCellIndices[ i ] <= view_particleCellIndices[ j ]; },
       [=] __cuda_callable__ ( int i, int j ) mutable {
         swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
         swap( view_points[ i ], view_points[ j ] );
         swap( view_rho[ i ], view_rho[ j ] );
         swap( view_v[ i ], view_v[ j ] );
         swap( view_rhoO[ i ], view_rhoO[ j ] );
         swap( view_vO[ i ], view_vO[ j ] );
         } );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariablesThrust()
{
   GlobalIndexType numberOfParticle = particles->getNumberOfParticles();
   auto view_particleCellIndices = particles->getParticleCellIndices().getView();
   auto view_points = particles->getPoints().getView();

   auto view_rho = this->FluidVariables.rho.getView();
   auto view_v = this->FluidVariables.v.getView();
   auto view_rhoO = rhoO.getView();
   auto view_vO = vO.getView();

#ifdef PREFER_SPEED_OVER_MEMORY
   //Reset indices:
   indicesMap.forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );

   auto view_indicesMap = indicesMap.getView();
   auto view_points_swap = points_swap.getView();
   auto view_rho_swap = rho_swap.getView();
   auto view_v_swap = v_swap.getView();
   auto view_rhoO_swap = rhoO_swap.getView();
   auto view_vO_swap = vO_swap.getView();

   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, view_indicesMap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_points.getArrayData(), view_points_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_rho.getArrayData(), view_rho_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_v.getArrayData(), view_v_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_rhoO.getArrayData(), view_rhoO_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_vO.getArrayData(), view_vO_swap.getArrayData() );

   particles->getPoints().swap( points_swap );
   FluidVariables.rho.swap( rho_swap );
   FluidVariables.v.swap( v_swap );
   rhoO.swap( rhoO_swap );
   vO.swap( vO_swap );
#else
   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, thrust::make_zip_iterator( thrust::make_tuple( view_points.getArrayData(), view_rho.getArrayData(), view_v.getArrayData(), view_rhoO.getArrayData(), view_vO.getArrayData() ) ) );
#endif
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename EquationOfState >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::ComputePressureFromDensity()
{
   auto view_rho = this->FluidVariables.rho.getView();
   auto view_p = this->FluidVariables.p.getView();

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
   auto rho_view = this->FluidVariables.rho.getView();
   auto v_view = this->FluidVariables.v.getView();
   auto r_view = this->particles->getPoints().getView();

   auto rhoO_view = this->rhoO.getView();
   auto vO_view = this->vO.getView();

   const auto drho_view = this->FluidVariables.drho.getView();
   const auto a_view = this->FluidVariables.a.getView();

   RealType dtdt05 = 0.5 * dt * dt;
   RealType dt2 = 2 * dt;

   auto init = [=] __cuda_callable__ ( int i ) mutable
   {
      r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
      vO_view[ i ] += a_view[ i ] * dt2;
      rhoO_view[ i ] += drho_view[ i ] * dt2;
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, this->particles->getNumberOfParticles(), init );

   FluidVariables.v.swap( vO );
   FluidVariables.rho.swap( rhoO );
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IntegrateEuler( RealType dt )
{
   auto rho_view = this->FluidVariables.rho.getView();
   auto v_view = this->FluidVariables.v.getView();
   auto r_view = this->particles->getPoints().getView();

   auto rhoO_view = this->rhoO.getView();
   auto vO_view = this->vO.getView();

   const auto drho_view = this->FluidVariables.drho.getView();
   const auto a_view = this->FluidVariables.a.getView();

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

