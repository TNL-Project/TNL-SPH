#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#include "Variables.h"
#include "../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
class IntegratorVariables
{
   public:
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;

   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

   using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHFluidConfig::DeviceType >;

   IntegratorVariables( GlobalIndexType size )
   : rho_old( size ), v_old( size ), rho_old_swap( size ), v_old_swap( size )
   {
      rho_old = 1000.; //TODO: Fix this.
      v_old = 0.;
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();

      auto view_rho_old_swap = rho_old_swap.getView();
      auto view_v_old_swap = v_old_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho_old.getArrayData(), view_rho_old_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v_old.getArrayData(), view_v_old_swap.getArrayData() );

      rho_old.swap( rho_old_swap );
      v_old.swap( v_old_swap );
   }

   ScalarArrayType rho_old;
   VectorArrayType v_old;

   ScalarArrayType rho_old_swap;
   VectorArrayType v_old_swap;
};

template< typename ModelPointer, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig > >
class VerletIntegrator
{
public:

   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
   using DeviceType = typename SPHFluidConfig::DeviceType;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;
   using VectorType = typename SPHFluidTraitsType::VectorType;

   using IntegratorVariablesType = IntegratorVariables< SPHFluidConfig >;
   using IntegratorVariablesPointer = typename Pointers::SharedPointer< IntegratorVariablesType, DeviceType >;

   VerletIntegrator() = default;

   template< typename FluidPointer >
   void
   IntegrateVerlet( RealType dt, FluidPointer& fluid )
   {
      auto rho_view = fluid->variables->rho.getView();
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();

      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();

      const auto drho_view = fluid->variables->drho.getView();
      const auto a_view = fluid->variables->a.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] += a_view[ i ] * dt2;
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, fluid->particles->getNumberOfParticles(), init );

      fluid->variables->v.swap( fluid->integratorVariables->v_old );
      fluid->variables->rho.swap( fluid->integratorVariables->rho_old );
   }

   template< typename BoundaryPointer >
   void
   IntegrateVerletBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] += drho_view[ i ] * dt2;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, boundary->particles->getNumberOfParticles(), init );

      boundary->variables->rho.swap( boundary->integratorVariables->rho_old );
   }

   template< typename FluidPointer >
   void
   IntegrateEuler( RealType dt, FluidPointer& fluid )
   {
      auto rho_view = fluid->variables->rho.getView();
      auto v_view = fluid->variables->v.getView();
      auto r_view = fluid->particles->getPoints().getView();

      auto rho_old_view = fluid->integratorVariables->rho_old.getView();
      auto v_old_view = fluid->integratorVariables->v_old.getView();

      const auto drho_view = fluid->variables->drho.getView();
      const auto a_view = fluid->variables->a.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         r_view[ i ] += v_view[ i ] * dt + a_view[ i ] * dtdt05;
         v_old_view[ i ] = v_view[ i ];
         v_view[ i ] += a_view[ i ] * dt;
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, fluid->particles->getNumberOfParticles(), init );
   }

   template< typename BoundaryPointer >
   void
   IntegrateEulerBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         rho_old_view[ i ] = rho_view[ i ];
         rho_view[ i ] += drho_view[ i ] * dt;
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, boundary->particles->getNumberOfParticles(), init );
   }

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   updateBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
   {
      const GlobalIndexType numberOfParticle = fluid->particles->getNumberOfParticles();
      const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

      //Buffer
      auto view_r_buffer = openBoundary->particles->getPoints().getView();
      auto view_v_buffer = openBoundary->variables->v.getView();
      auto view_rho_buffer = openBoundary->variables->rho.getView();

      auto view_inletMark = openBoundary->particleMark->getView();
      view_inletMark = 1; //TODO: this can be avoided

      const VectorType inletOrientation = openBoundary->parameters.orientation;
      const VectorType inletConstVelocity = openBoundary->parameters.velocity;
      const RealType inletConstDensity = openBoundary->parameters.density;
      const RealType bufferEdge = openBoundary->parameters.bufferEdge;
      const VectorType bufferWidth = openBoundary->parameters.bufferWidth;

      //Fluid ( variable arrays and integrator variable arrays )
      auto view_r_fluid = fluid->particles->getPoints().getView();
      auto view_v_fluid = fluid->variables->v.getView();
      auto view_rho_fluid = fluid->variables->rho.getView();

      auto view_rho_old = fluid->integratorVariables->rho_old.getView();
      auto view_v_old = fluid->integratorVariables->v_old.getView();

      auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
         //TODO: Ugly, ugly stuff.. Keep the scope inside, remove everything else.
         if( inletOrientation[ 0 ] > 0 ) //This isn't.
         {
            if( ( view_r_buffer[ i ], inletOrientation ) > bufferEdge ) //It is possible to do for each direction. This if statemens is good.
            {
               view_inletMark[ i ] = 0;
            }
         }
         else if( inletOrientation[ 0 ] < 0 ) //This isn't.
         {
            if( ( view_r_buffer[ i ], inletOrientation ) > ( -1.f * bufferEdge ) ) //It is possible to do for each direction. This if statemens is good.
            {
               view_inletMark[ i ] = 0;
            }
         }
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfBufferParticles, moveBufferParticles );

      auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_inletMark[ i ]; };
      auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      const GlobalIndexType numberOfRetyped = numberOfBufferParticles - Algorithms::reduce< DeviceType >( 0, view_inletMark.getSize(), fetch, reduction, 0.0 ); //I like zeros baceause sort.
      std::cout << "... InletBuffer - number of retyped particles: " << numberOfRetyped << std::endl;

      if( numberOfRetyped == 0 )
         return;

      //Sort particles by mark //TODO: can be this avoided?
      thrust::sort_by_key( thrust::device, view_inletMark.getArrayData(), view_inletMark.getArrayData() + numberOfBufferParticles,
            thrust::make_zip_iterator( thrust::make_tuple( view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );
      std::cout << "... InletBuffer - particles sorted." << std::endl;

      std::cout << ".................. numberOfParticles: " << fluid->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfAllocatedParticles: " << fluid->particles->getNumberOfAllocatedParticles() << std::endl;
      std::cout << ".................. numberOfParticles: " << fluid->particles->getNumberOfParticles() << std::endl;

      auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         if( view_inletMark[ i ] == 0 )
         {
            //Retype buffer particle to fluid
            view_r_fluid[ numberOfParticle + i ] = view_r_buffer[ i ];
            view_rho_fluid[ numberOfParticle + i ] = view_rho_buffer[ i ];
            view_v_fluid[ numberOfParticle + i ] = view_v_buffer[ i ];
            view_rho_old[ numberOfParticle + i ] = view_rho_buffer[ i ];
            view_v_old[ numberOfParticle + i ] = view_v_buffer[ i ];

            //Generate new bufffer partice
            //const VectorType newBufferParticle = view_r_buffer[ i ] - bufferWidth ; //This works in 1d
            const VectorType newBufferParticle = view_r_buffer[ i ] - bufferWidth * inletOrientation[ 0 ] ; //This works in 1d
            view_r_buffer[ i ] = newBufferParticle;
            view_v_buffer[ i ] = inletConstVelocity;
            view_rho_buffer[ i ] = inletConstDensity;
         }
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfRetyped, createNewFluidParticles );
      fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
      std::cout << "... InletBuffer - system updated." << std::endl;
   }
};

} // SPH
} // ParticleSystem
} // TNL

