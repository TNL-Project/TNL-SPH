#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "Variables.h"
#include "../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ModelPointer, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig > >
class VerletIntegrator
{
public:

   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename SPHFluidConfig::DeviceType;

   using RealType = typename Variables::RealType;
   using GlobalIndexType = typename Variables::GlobalIndexType;

   using ScalarType = typename SPHFluidTraitsType::ScalarType;
   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorType = typename SPHFluidTraitsType::VectorType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

   VerletIntegrator( ModelPointer model, GlobalIndexType size, GlobalIndexType size_boundary )
   : rho_old( size ), v_old( size ), rho_old_swap( size ), v_old_swap( size ),
   rhoBoundary_old( size_boundary ), vBoundary_old( size_boundary ), rhoBoundary_old_swap( size_boundary ), vBoundary_old_swap( size_boundary ),
   model( model )
   {
      //rho_old = model->getFluidVariables().rho;
      //v_old = model->getFluidVariables().v;
      rho_old = 1000.f; rho_old_swap = 1000.f;
      v_old = 0.f; v_old_swap = 0.f;

      rhoBoundary_old = 1000.f; rhoBoundary_old_swap = 1000.f;
      vBoundary_old = 0.f; vBoundary_old_swap = 0.f;
   }

   void
   IntegrateVerlet( RealType dt )
   {
      auto rho_view = model->getFluidVariables().rho.getView();
      auto v_view = model->getFluidVariables().v.getView();
      auto r_view = model->particles->getPoints().getView();

      auto rho_old_view = this->rho_old.getView();
      auto v_old_view = this->v_old.getView();

      const auto drho_view = model->getFluidVariables().drho.getView();
      const auto a_view = model->getFluidVariables().a.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] += a_view[ i ] * dt2;
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, model->particles->getNumberOfParticles(), init );

      model->getFluidVariables().v.swap( v_old );
      model->getFluidVariables().rho.swap( rho_old );
   }

   void
   IntegrateVerletBoundary( RealType dt )
   {
      auto rho_view = model->getBoundaryVariables().rho.getView();
      auto rho_old_view = this->rhoBoundary_old.getView();

      const auto drho_view = model->getBoundaryVariables().drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, model->boundaryParticles->getNumberOfParticles(), init );

      model->getBoundaryVariables().rho.swap( rhoBoundary_old );
   }

   void
   IntegrateEuler( RealType dt )
   {
      auto rho_view = model->getFluidVariables().rho.getView();
      auto v_view = model->getFluidVariables().v.getView();
      auto r_view = model->particles->getPoints().getView();

      auto rho_old_view = this->rho_old.getView();
      auto v_old_view = this->v_old.getView();

      const auto drho_view = model->getFluidVariables().drho.getView();
      const auto a_view = model->getFluidVariables().a.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] = v_view[ i ];
         v_view[ i ] += a_view[ i ] * dt;
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, model->particles->getNumberOfParticles(), init );
   }

   void
   IntegrateEulerBoundary( RealType dt )
   {
      auto rho_view = model->getBoundaryVariables().rho.getView();
      auto rho_old_view = this->rhoBoundary_old.getView();

      const auto drho_view = model->getBoundaryVariables().drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, model->boundaryParticles->getNumberOfParticles(), init );
   }

   void
   sortIntegratorArrays()
   {
      GlobalIndexType numberOfParticle = model->particles->getNumberOfParticles();
      auto view_indicesMap = model->swapFluid.indicesMap.getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();

      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      thrust::gather( thrust::device, view_indicesMap.getArrayData(), view_indicesMap.getArrayData() + numberOfParticle, view_rho_old.getArrayData(), view_rho_old_swap.getArrayData() );
      thrust::gather( thrust::device, view_indicesMap.getArrayData(), view_indicesMap.getArrayData() + numberOfParticle, view_v_old.getArrayData(), view_v_old_swap.getArrayData() );

      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   void
   sortIntegratorBoundaryArrays()
   {
      GlobalIndexType numberOfParticle = model->boundaryParticles->getNumberOfParticles();
      auto view_indicesMap = model->swapBoundary.indicesMap.getView();

      auto view_rho_old = rhoBoundary_old.getView();
      auto view_v_old = vBoundary_old.getView();

      auto view_rho_old_swap = rhoBoundary_old_swap.getView();
      auto view_v_old_swap = vBoundary_old_swap.getView();

      thrust::gather( thrust::device, view_indicesMap.getArrayData(), view_indicesMap.getArrayData() + numberOfParticle, view_rho_old.getArrayData(), view_rho_old_swap.getArrayData() );
      thrust::gather( thrust::device, view_indicesMap.getArrayData(), view_indicesMap.getArrayData() + numberOfParticle, view_v_old.getArrayData(), view_v_old_swap.getArrayData() );

      rhoBoundary_old.swap( rhoBoundary_old_swap );
      vBoundary_old.swap( vBoundary_old_swap );
   }

protected:

   ModelPointer model;

   ScalarArrayType rho_old;
   VectorArrayType v_old;

   ScalarArrayType rhoBoundary_old;
   VectorArrayType vBoundary_old;

#ifdef PREFER_SPEED_OVER_MEMORY
   ScalarArrayType rho_old_swap;
   VectorArrayType v_old_swap;

   ScalarArrayType rhoBoundary_old_swap;
   VectorArrayType vBoundary_old_swap;
#endif
};

} // SPH
} // ParticleSystem
} // TNL

