#include "Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//const Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getFluidVariables() const
//{
//   return this->fluidVariables;
//}
//
//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getFluidVariables()
//{
//   return this->fluidVariables;
//}
//
//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//const Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getBoundaryVariables() const
//{
//   return this->boundaryVariables;
//}
//
//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getBoundaryVariables()
//{
//   return this->boundaryVariables;
//}
//
//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//const Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getInletVariables() const
//{
//   return this->inletVariables;
//}
//
//template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//Variables&
//WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::getInletVariables()
//{
//   return this->inletVariables;
//}

template< typename Particles, typename SPHFluidConfig, typename Variables >
const typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IndexArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getIndicesForReoder() const
{
   return swapFluid.indicesMap;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
typename WCSPH_DBC< Particles, SPHFluidConfig, Variables >::IndexArrayType&
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::getIndicesForReoder()
{
   return swapFluid.indicesMap;
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortParticlesAndVariablesThrust( ParticlePointer& particleSys, VariablesPointer& variables, SwapVariables& variables_swap )
{
   GlobalIndexType numberOfParticle = particleSys->getNumberOfParticles();
   auto view_particleCellIndices = particleSys->getParticleCellIndices().getView();
   auto view_points = particleSys->getPoints().getView();

   //auto view_rho = this->FluidVariables.rho.getView();
   //auto view_v = this->FluidVariables.v.getView();
   auto view_rho = variables->rho.getView();
   auto view_v = variables->v.getView();

#ifdef PREFER_SPEED_OVER_MEMORY
   //Reset indices:
   variables_swap.indicesMap.forAllElements( [] __cuda_callable__ ( int i, int& value ) { value = i; } );

   auto view_indicesMap = variables_swap.indicesMap.getView();
   auto view_points_swap = variables_swap.points_swap.getView();
   auto view_rho_swap = variables_swap.rho_swap.getView();
   auto view_v_swap = variables_swap.v_swap.getView();

   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, view_indicesMap.getArrayData() );
   thrust::gather( thrust::device, variables_swap.indicesMap.getArrayData(), variables_swap.indicesMap.getArrayData() + numberOfParticle, view_points.getArrayData(), view_points_swap.getArrayData() );
   thrust::gather( thrust::device, variables_swap.indicesMap.getArrayData(), variables_swap.indicesMap.getArrayData() + numberOfParticle, view_rho.getArrayData(), view_rho_swap.getArrayData() );
   thrust::gather( thrust::device, variables_swap.indicesMap.getArrayData(), variables_swap.indicesMap.getArrayData() + numberOfParticle, view_v.getArrayData(), view_v_swap.getArrayData() );

   particleSys->getPoints().swap( variables_swap.points_swap );
   variables->rho.swap( variables_swap.rho_swap );
   variables->v.swap( variables_swap.v_swap );
#else
   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, thrust::make_zip_iterator( thrust::make_tuple( view_points.getArrayData(), view_rho.getArrayData(), view_v.getArrayData(), view_rhoO.getArrayData(), view_vO.getArrayData() ) ) );
#endif
}

template< typename Particles, typename SPHFluidConfig, typename Variables >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::sortBoundaryParticlesAndVariablesThrust( ParticlePointer& particleSys, VariablesPointer& variables, SwapVariables& variables_swap )
{
   GlobalIndexType numberOfParticle = particleSys->getNumberOfParticles();
   auto view_particleCellIndices = particleSys->getParticleCellIndices().getView();
   auto view_points = particleSys->getPoints().getView();

   auto view_indicesMap = variables_swap.indicesMap.getView();
   auto view_rho = variables->rho.getView();
   auto view_v = variables->v.getView();

   thrust::sort_by_key( thrust::device, view_particleCellIndices.getArrayData(), view_particleCellIndices.getArrayData() + numberOfParticle, thrust::make_zip_iterator( thrust::make_tuple( view_points.getArrayData(), view_rho.getArrayData(), view_v.getArrayData(), view_indicesMap.getArrayData() ) ) );
}



//USETHIS:template< typename Particles, typename OpenBoundary, typename SPHFluidConfig, typename Variables >
//USETHIS:template< typename EquationOfState >
//USETHIS:void
//USETHIS:WCSPH_DBC< Particles, OpenBoundary, SPHFluidConfig, Variables >::ComputePressureFromDensity()
//USETHIS:{
//USETHIS:   auto view_rho = this->getFluidVariables().rho.getView();
//USETHIS:   auto view_p = this->getFluidVariables().p.getView();
//USETHIS:
//USETHIS:   auto init = [=] __cuda_callable__ ( int i ) mutable
//USETHIS:   {
//USETHIS:      view_p[ i ] = EquationOfState::DensityToPressure( view_rho[ i ] );
//USETHIS:   };
//USETHIS:   Algorithms::ParallelFor< DeviceType >::exec( 0, particles->getNumberOfParticles(), init );
//USETHIS:}

} // SPH
} // ParticleSystem
} // TNL

