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

   typename SPHFluidTraitsType::IndexArrayType markOutlet;

   template< typename FluidPointer >
   void
   integrateVerlet( RealType dt, FluidPointer& fluid )
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
      Algorithms::parallelFor< DeviceType >( 0, fluid->particles->getNumberOfParticles(), init );

      fluid->variables->v.swap( fluid->integratorVariables->v_old );
      fluid->variables->rho.swap( fluid->integratorVariables->rho_old );
   }

   template< typename BoundaryPointer >
   void
   integrateVerletBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;
      RealType dt2 = 2 * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         //rho_old_view[ i ] += drho_view[ i ] * dt2;

         const RealType rho = rho_old_view[ i ];
         const RealType rho_new = rho + drho_view[ i ] * dt2;
         rho_old_view[ i ] = ( rho_new > 1000.f ? rho_new : 1000.f ); //TODO: Use rho ref instead of 1000.f
      };
      Algorithms::parallelFor< DeviceType >( 0, boundary->particles->getNumberOfParticles(), init );

      boundary->variables->rho.swap( boundary->integratorVariables->rho_old );
   }

   template< typename FluidPointer >
   void
   integrateEuler( RealType dt, FluidPointer& fluid )
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
      Algorithms::parallelFor< DeviceType >( 0, fluid->particles->getNumberOfParticles(), init );
   }

   template< typename BoundaryPointer >
   void
   integrateEulerBoundary( RealType dt, BoundaryPointer& boundary )
   {
      auto rho_view = boundary->variables->rho.getView();
      auto rho_old_view = boundary->integratorVariables->rho_old.getView();

      const auto drho_view = boundary->variables->drho.getView();

      RealType dtdt05 = 0.5 * dt * dt;

      auto init = [=] __cuda_callable__ ( int i ) mutable
      {
         //rho_old_view[ i ] = rho_view[ i ];
         //rho_view[ i ] += drho_view[ i ] * dt;

         const RealType rho = rho_view[ i ];
         rho_old_view[ i ] = rho;
         const RealType rho_new = rho + drho_view[ i ] * dt;
         rho_view[ i ] = ( rho_new > 1000.f ? rho_new : 1000.f ); //TODO: Use rho ref instead of 1000.f
      };
      Algorithms::parallelFor< DeviceType >( 0, boundary->particles->getNumberOfParticles(), init );
   }

   template< typename FluidPointer, typename BoundaryPointer, typename TimeStepping >
   void
   integratStepVerlet( FluidPointer& fluid, BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      if( timeStepping.getStep() % 20 == 0 ) {
         integrateEuler( timeStepping.getTimeStep(), fluid ); //TODO: Timer!
         integrateEulerBoundary( timeStepping.getTimeStep(), boundary );
      }
      else {
         integrateVerlet( timeStepping.getTimeStep(), fluid );
         integrateVerletBoundary( timeStepping.getTimeStep(), boundary );
      }
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

      auto view_inletMark = openBoundary->variables->particleMark.getView();
      view_inletMark = 1; //TODO: this can be avoided

      const VectorType inletOrientation = openBoundary->parameters.orientation;
      const VectorType inletConstVelocity = openBoundary->parameters.velocity;
      const RealType inletConstDensity = openBoundary->parameters.density;
      const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
      const VectorType bufferPosition = openBoundary->parameters.position;

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
         const VectorType r = view_r_buffer[ i ];
         //const VectorType r_relative = r - bufferPosition;
         const VectorType r_relative = bufferPosition - r;

         //if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] ) //It is possible to do for each direction. This if statemens is good.
         //if( ( r_relative, inletOrientation ) > 0 ) //It is possible to do for each direction. This if statemens is good.
         if( ( r_relative, inletOrientation ) <= 0.f ) //It is possible to do for each direction. This if statemens is good.
         {
            view_inletMark[ i ] = 0;
            printf("sc: %f ", ( r_relative, inletOrientation ));
         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles );

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
            //:const VectorType r_relative = view_r_buffer[ i ] - bufferPosition;
            //:const VectorType newBufferParticle = view_r_buffer[ i ] - ( r_relative, inletOrientation ) * inletOrientation; //This works in 1d
            const VectorType r_relative = bufferPosition - view_r_buffer[ i ];
            const VectorType newBufferParticle = view_r_buffer[ i ] - ( r_relative, inletOrientation ) * inletOrientation - bufferWidth[ 0 ] * inletOrientation; //This works in 1d

            view_r_buffer[ i ] = newBufferParticle;
            view_v_buffer[ i ] = inletConstVelocity;
            view_rho_buffer[ i ] = inletConstDensity;
         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfRetyped, createNewFluidParticles );
      fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );

      //added
      view_v_buffer = inletConstVelocity;
      std::cout << "... InletBuffer - system updated." << std::endl;
   }

   template< typename FluidPointer, typename OpenBoundaryPointer >
   void
   updateOutletBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
   {
      const GlobalIndexType numberOfParticle = fluid->particles->getNumberOfParticles();
      GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

      //Buffer
      auto view_r_buffer = openBoundary->particles->getPoints().getView();
      auto view_v_buffer = openBoundary->variables->v.getView();
      auto view_rho_buffer = openBoundary->variables->rho.getView();

      auto view_inletMark = openBoundary->variables->particleMark.getView();
      view_inletMark = 0; //TODO: this can be avoided

      const VectorType inletOrientation = openBoundary->parameters.orientation;
      const VectorType inletConstVelocity = openBoundary->parameters.velocity;
      const RealType inletConstDensity = openBoundary->parameters.density;
      const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
      const VectorType bufferPosition = openBoundary->parameters.position;

      //Fluid ( variable arrays and integrator variable arrays )
      auto view_r_fluid = fluid->particles->getPoints().getView();
      auto view_v_fluid = fluid->variables->v.getView();
      auto view_rho_fluid = fluid->variables->rho.getView();

      auto view_rho_old = fluid->integratorVariables->rho_old.getView();
      auto view_v_old = fluid->integratorVariables->v_old.getView();

      std::cout << "************************************OUTLET-START******************************************** " << std::endl;

      std::cout << "... OutletBuffer - move buffer particles, " << std::endl;
      auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
         const VectorType r = view_r_buffer[ i ];
         const VectorType r_relative = bufferPosition - r;

         if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] ) //It is possible to do for each direction. This if statemens is good.
         //if( ( r_relative, inletOrientation ) > 0 ) //It is possible to do for each direction. This if statemens is good.
         //if( ( r_relative, inletOrientation ) <= 0.f ) //It is possible to do for each direction. This if statemens is good.
         {
            view_inletMark[ i ] = 1;
            printf("--> moveBufferParticles: sc: %f r: [ %f, %f ]", ( r_relative, inletOrientation ), view_r_buffer[ i ][ 0 ], view_r_buffer[ i ][ 1 ] );
         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles );
      std::cout << std::endl << ".................. moved." << std::endl;

      auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_inletMark[ i ]; };
      auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      //const GlobalIndexType numberOfRetyped = numberOfBufferParticles - Algorithms::reduce< DeviceType >( 0, view_inletMark.getSize(), fetch, reduction, 0.0 ); //I like zeros baceause sort.
      const GlobalIndexType numberOfRetyped = Algorithms::reduce< DeviceType >( 0, view_inletMark.getSize(), fetch, reduction, 0.0 ); //I like zeros baceause sort.
      std::cout << "... OutletBuffer - number of retyped particles: " << numberOfRetyped << std::endl;

      //if( numberOfRetyped == 0 )
      //   return;

      //Sort particles by mark //TODO: can be this avoided?
      thrust::sort_by_key( thrust::device, view_inletMark.getArrayData(), view_inletMark.getArrayData() + numberOfBufferParticles,
            thrust::make_zip_iterator( thrust::make_tuple( view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );
      std::cout << "... OutletBuffer - particles sorted." << std::endl;

      std::cout << ".................. numberOfParticles: " << fluid->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfAllocatedParticles: " << fluid->particles->getNumberOfAllocatedParticles() << std::endl;
      std::cout << ".................. numberOfParticles: " << fluid->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfBufferParticles: " << openBoundary->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfBufferAlocatedParticles: " << openBoundary->particles->getNumberOfAllocatedParticles() << std::endl;

      auto discardBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         if( view_inletMark[ i ] == 1 )
         {

            //Discard old buffer particles
            view_r_buffer[ i ] = FLT_MAX;
            view_rho_buffer[ i ] = 0.f;
            view_v_buffer[ i ] = 0.f;

         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, discardBufferParticles ); //this can loop over retyped range

      openBoundary->particles->setNumberOfParticles( numberOfBufferParticles - numberOfRetyped );
      numberOfBufferParticles = numberOfBufferParticles - numberOfRetyped;
      std::cout << "... OutletBuffer - move and discard buffer particles - system updated." << std::endl;
      std::cout << ".................. numberOfBufferParticles: " << openBoundary->particles->getNumberOfParticles() << std::endl;

      //At this moment, oulet buffer is updated, we need to just resize the fluid.
      auto view_fluidOutMark = fluid->variables->fluidOutMark.getView();
      auto view_fluidOutMarkIndex = fluid->variables->fluidOutMarkIndex.getView();
      view_fluidOutMark = 0;
      view_fluidOutMarkIndex = 0;

      auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         const VectorType r = view_r_fluid[ i ];
         const VectorType r_relative = bufferPosition - r;

         if( ( r_relative, inletOrientation ) > 0 )
         {
            view_fluidOutMark[ i ] = 1;
            view_fluidOutMarkIndex[ i ] = i;
            printf("-> checkFluidParticles: %f for r: [ %f, %f ]\n", ( r_relative, inletOrientation ), r[ 0 ], r[ 1 ] );
         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfParticle, checkFluidParticles );

      auto fetchFluid = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_fluidOutMark[ i ]; };
      auto reductionFluid = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      //const GlobalIndexType numberOfFluidToBuffer = Algorithms::reduce< DeviceType >( 0, view_fluidOutMark.getSize(), fetchFluid, reductionFluid, 0.0 ); //I like zeros baceause sort.
      const GlobalIndexType numberOfFluidToBuffer = Algorithms::reduce< DeviceType >( 0, numberOfParticle, fetchFluid, reductionFluid, 0.0 ); //I like zeros baceause sort.

      //Sort particles by mark //TODO: can be this avoided?
      thrust::sort_by_key( thrust::device, view_fluidOutMark.getArrayData(), view_fluidOutMark.getArrayData() + numberOfParticle,
            view_fluidOutMarkIndex.getArrayData() );


      const GlobalIndexType numberOfParticleNew = numberOfParticle - numberOfFluidToBuffer;

      std::cout << "... OutletBuffer - fluid to buffer - particles sorted." << std::endl;
      std::cout << ".................. numberOfParticles: " << fluid->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfAllocatedParticles: " << fluid->particles->getNumberOfAllocatedParticles() << std::endl;
      std::cout << ".................. newNumberOfParticles (n - retyped): " << numberOfParticleNew << std::endl;
      std::cout << ".................. numberOfFluidToBuffer: " << numberOfFluidToBuffer << std::endl;

      auto retypeFluidToOutlet = [=] __cuda_callable__ ( int i ) mutable
      {
         const GlobalIndexType p = view_fluidOutMarkIndex[ numberOfParticle - i - 1];
         if(  view_fluidOutMark[ numberOfParticle - i - 1 ] == 1 )
         {
         printf( "-> retypeFluidToOutlet: p: %d \n", p );
         view_r_buffer[ numberOfBufferParticles + i ] = view_r_fluid[ p ];
         view_rho_buffer[ numberOfBufferParticles + i ] = view_rho_fluid[ p ];
         view_v_buffer[ numberOfBufferParticles + i ] = view_v_fluid[ p ];
         printf( "-> retypeFluidToOutlet: r_bufferNew: [ %f, %f ]\n", view_r_buffer[ numberOfBufferParticles + i ][ 0 ], view_r_buffer[ numberOfBufferParticles + i ][ 1 ]);


         //Deactivate the old fluid particles
         //swap( view_r_fluid[ p ], view_r_fluid[ numberOfParticles - i ]);

         view_r_fluid[ p ] = FLT_MAX;
         swap( view_r_fluid[ p ], view_r_fluid[ numberOfParticle - i - 1] );
         swap( view_rho_fluid[ p ], view_rho_fluid[ numberOfParticle - i -1] );
         swap( view_v_fluid[ p ], view_v_fluid[ numberOfParticle - i -1] );
         swap( view_rho_old[ p ], view_rho_old[ numberOfParticle - i -1] );
         swap( view_v_old[ p ], view_v_old[ numberOfParticle - i -1] );
         }


      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfFluidToBuffer, retypeFluidToOutlet );

      openBoundary->particles->setNumberOfParticles( numberOfBufferParticles + numberOfFluidToBuffer );
      fluid->particles->setNumberOfParticles( numberOfParticleNew );
      std::cout << "... OutletBuffer - upradetd." << std::endl;
      std::cout << ".................. numberOfBufferParticles: " << openBoundary->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. numberOfBufferAlocatedParticles: " << openBoundary->particles->getNumberOfAllocatedParticles() << std::endl;
      //std::cout << ".................. bufferParticles: " << openBoundary->particles->getPoints() << std::endl;
      std::cout << "... OutletBuffer - finished." << std::endl;
      std::cout << "************************************OUTLET-FINISH********************************************* " << std::endl;

   }



};

} // SPH
} // ParticleSystem
} // TNL

