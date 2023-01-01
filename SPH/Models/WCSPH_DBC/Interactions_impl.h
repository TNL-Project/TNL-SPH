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
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IndexArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getIndicesForReoder() const
{
   return indicesMap;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IndexArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getIndicesForReoder()
{
   return indicesMap;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariables()
{
   auto view_particleCellIndices = particles->getParticleCellIndices().getView();
   auto view_points = particles->getPoints().getView();
   auto view_rho = this->FluidVariables.rho.getView();
   auto view_v = this->FluidVariables.v.getView();

   Algorithms::sort< DeviceType, GlobalIndexType >(
       0, particles->getNumberOfParticles(),
       [=] __cuda_callable__ ( int i, int j ) -> bool {
         return view_particleCellIndices[ i ] <= view_particleCellIndices[ j ]; },
       [=] __cuda_callable__ ( int i, int j ) mutable {
         swap( view_particleCellIndices[ i ], view_particleCellIndices[ j ] );
         swap( view_points[ i ], view_points[ j ] );
         swap( view_rho[ i ], view_rho[ j ] );
         swap( view_v[ i ], view_v[ j ] );
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

#ifdef PREFER_SPEED_OVER_MEMORY
   //Reset indices:
   indicesMap.forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );

   auto view_indicesMap = indicesMap.getView();
   auto view_points_swap = points_swap.getView();
   auto view_rho_swap = rho_swap.getView();
   auto view_v_swap = v_swap.getView();

   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, view_indicesMap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_points.getArrayData(), view_points_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_rho.getArrayData(), view_rho_swap.getArrayData() );
   thrust::gather( thrust::device, indicesMap.getArrayData(), indicesMap.getArrayData() + numberOfParticle, view_v.getArrayData(), view_v_swap.getArrayData() );

   particles->getPoints().swap( points_swap );
   FluidVariables.rho.swap( rho_swap );
   FluidVariables.v.swap( v_swap );
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

} // SPH
} // ParticleSystem
} // TNL

