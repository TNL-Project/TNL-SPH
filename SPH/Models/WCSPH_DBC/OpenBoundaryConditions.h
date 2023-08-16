#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< SPHFluidConfig >::updateBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfParticle = fluid->particles->getNumberOfParticles();
   const GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   //Load buffer fields
   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_rho_buffer = openBoundary->variables->rho.getView();
   auto view_inletMark = openBoundary->variables->particleMark.getView();
   view_inletMark = 1;

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType inletConstVelocity = openBoundary->parameters.velocity;
   const RealType inletConstDensity = openBoundary->parameters.density;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   const VectorType bufferPosition = openBoundary->parameters.position;

   //Load fluid fileds (variables + integrator variables)
   auto view_r_fluid = fluid->particles->getPoints().getView();
   auto view_v_fluid = fluid->variables->v.getView();
   auto view_rho_fluid = fluid->variables->rho.getView();

   auto view_rho_old = fluid->integratorVariables->rho_old.getView();
   auto view_v_old = fluid->integratorVariables->v_old.getView();

   //Move buffer particles
   auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      //Check if particle leaves the buffer
      if( ( r_relative, inletOrientation ) <= 0.f )
         view_inletMark[ i ] = 0;

      return view_inletMark[ i ];
   };
   const GlobalIndexType numberOfNotRetyped = Algorithms::reduce< DeviceType >(
         0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );
   const GlobalIndexType numberOfRetyped = numberOfBufferParticles - numberOfNotRetyped;

   if( numberOfRetyped == 0 )
      return;

   //Sort particles by mark
   thrust::sort_by_key( thrust::device, view_inletMark.getArrayData(),
                        view_inletMark.getArrayData() + numberOfBufferParticles,
                        thrust::make_zip_iterator( thrust::make_tuple( view_r_buffer.getArrayData(),
                                                                       view_v_buffer.getArrayData(),
                                                                       view_rho_buffer.getArrayData() ) ) );

   //Retype fluid to buffer
   auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      if( view_inletMark[ i ] == 0 )
      {
         view_r_fluid[ numberOfParticle + i ] = view_r_buffer[ i ];
         view_rho_fluid[ numberOfParticle + i ] = view_rho_buffer[ i ];
         view_v_fluid[ numberOfParticle + i ] = view_v_buffer[ i ];
         view_rho_old[ numberOfParticle + i ] = view_rho_buffer[ i ];
         view_v_old[ numberOfParticle + i ] = view_v_buffer[ i ];

         //const VectorType r_relative = bufferPosition - view_r_buffer[ i ];
         //const VectorType newBufferParticle = view_r_buffer[ i ] - ( r_relative, inletOrientation ) * inletOrientation - bufferWidth[ 0 ] * inletOrientation;

         const VectorType newBufferParticle = view_r_buffer[ i ] - bufferWidth[ 0 ] * inletOrientation;

         view_r_buffer[ i ] = newBufferParticle;
         view_v_buffer[ i ] = inletConstVelocity;
         //view_rho_buffer[ i ] = inletConstDensity; //TODO: Keep the density same.
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfRetyped, createNewFluidParticles );

   //Update number of particles
   fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() + numberOfRetyped );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() + numberOfRetyped );
}

template< typename SPHFluidConfig >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< SPHFluidConfig >::updateOutletBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
{
   const GlobalIndexType numberOfParticle = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfBufferParticles = openBoundary->particles->getNumberOfParticles();

   //Load buffer fields
   auto view_r_buffer = openBoundary->particles->getPoints().getView();
   auto view_v_buffer = openBoundary->variables->v.getView();
   auto view_rho_buffer = openBoundary->variables->rho.getView();
   auto view_inletMark = openBoundary->variables->particleMark.getView();
   view_inletMark = 0;

   const VectorType inletOrientation = openBoundary->parameters.orientation;
   const VectorType inletConstVelocity = openBoundary->parameters.velocity;
   const RealType inletConstDensity = openBoundary->parameters.density;
   const VectorType bufferWidth = openBoundary->parameters.bufferWidth;
   const VectorType bufferPosition = openBoundary->parameters.position;

   //check for particles leaving the buffer
   auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      //view_r_buffer[ i ] += ( view_v_buffer[ i ], orientation ) * orientation * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      //buffer particle leaves the entire computation space
      if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] )
         view_inletMark[ i ] = 1;

      return view_inletMark[ i ];
   };
   const GlobalIndexType removeFromBufferCount = Algorithms::reduce< DeviceType >(
         0, numberOfBufferParticles, moveBufferParticles );

   //TODO: This should work with view_r_buffer[ i ] = FLT_MAX; instead of using view_inletMark[ i ] = 1.
   //      Maybe it should be done in splited functions. After sortParticles() I need to reload the vector views.
   //openBoundary->particles->resetListWithIndices();
   //openBoundary->particles->computeParticleCellIndices();
   //openBoundary->sortParticles();

   //TODO: Sort particle based on mark. This should be replaced with something.
   thrust::sort_by_key( thrust::device,
                        view_inletMark.getArrayData(),
                        view_inletMark.getArrayData() + numberOfBufferParticles,
                        thrust::make_zip_iterator( thrust::make_tuple(
                              view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );

   numberOfBufferParticles = numberOfBufferParticles - removeFromBufferCount;
   openBoundary->particles->setNumberOfParticles( numberOfBufferParticles - removeFromBufferCount );
   openBoundary->particles->setLastActiveParticle( openBoundary->getLastActiveParticle() - removeFromBufferCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() - removeFromBufferCount );

   //----- FLUID TO BUFFER -------------------------------------------------
   //Load all fluid fields
   auto view_r_fluid = fluid->particles->getPoints().getView();
   auto view_v_fluid = fluid->variables->v.getView();
   auto view_rho_fluid = fluid->variables->rho.getView();
   auto view_rho_old = fluid->integratorVariables->rho_old.getView();
   auto view_v_old = fluid->integratorVariables->v_old.getView();

   //TODO: Ugly ugly ugly temp workaround.
   const typename SPHFluidTraitsType::IndexVectorType gridIndex = TNL::floor(
         ( bufferPosition[ 0 ] - openBoundary->particles->getGridOrigin() ) / openBoundary->particles->getSearchRadius() );
   const GlobalIndexType gridColumnAuxTrick = gridIndex[ 0 ];

   //const PairIndexType particleRangeToCheck = fluid->particles->getFirstLastParticleInColumnOfCells( gridColumnAuxTrick );
   PairIndexType particleRangeToCheck;
   if constexpr( SPHFluidConfig::spaceDimension == 2 )
      particleRangeToCheck = fluid->particles->getFirstLastParticleInColumnOfCells( gridColumnAuxTrick );
   else if constexpr( SPHFluidConfig::spaceDimension == 3 )
      particleRangeToCheck = fluid->particles->getFirstLastParticleInBlockOfCells( gridColumnAuxTrick );

   auto receivingParticleMark_view = openBoundary->variables->receivingParticleMark.getView();
   receivingParticleMark_view = INT_MAX;

   //obtain fluid particles that should be retyped to buffer particles
   auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const VectorType r = view_r_fluid[ i ];
      const VectorType r_relative = bufferPosition - r;

      //fluid particle should be retyped to the buffer particle
      if( ( r_relative, inletOrientation ) > 0 ){
         printf(" place to: %d, place i: %d",  i - particleRangeToCheck[ 0 ], i );
         receivingParticleMark_view[ i - particleRangeToCheck[ 0 ] ] = i;
         return 1;
      }
      return 0;
   };
   const GlobalIndexType fluidToBufferCount = Algorithms::reduce< DeviceType >(
         particleRangeToCheck[ 0 ], particleRangeToCheck[ 1 ] + 1, checkFluidParticles, TNL::Plus() );

   thrust::sort( thrust::device,
                 receivingParticleMark_view.getArrayData(),
                 receivingParticleMark_view.getArrayData() + numberOfBufferParticles );

   //retype fluid particles to buffer particles
   auto retypeFluidToOutlet = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = receivingParticleMark_view[ i ];

      view_r_buffer[ numberOfBufferParticles + i ] = view_r_fluid[ p ];
      view_rho_buffer[ numberOfBufferParticles + i ] = view_rho_fluid[ p ];
      view_v_buffer[ numberOfBufferParticles + i ] = view_v_fluid[ p ];

      view_r_fluid[ p ] = FLT_MAX;
   };
   Algorithms::parallelFor< DeviceType >( 0, fluidToBufferCount, retypeFluidToOutlet );

   //update particle ranges for buffer and fluid
   openBoundary->particles->setNumberOfParticles( numberOfBufferParticles + fluidToBufferCount );
   openBoundary->particles->setLastActiveParticle( openBoundary->particles->getLastActiveParticle() + fluidToBufferCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() + fluidToBufferCount );

   openBoundary->numberOfFluidParticlesToRemove = fluidToBufferCount;

}



} // SPH
} // ParticleSystem
} // TNL

