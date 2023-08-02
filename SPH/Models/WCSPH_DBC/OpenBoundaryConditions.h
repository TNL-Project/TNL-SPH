#include "Integrator.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ModelPointer, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< ModelPointer, SPHFluidConfig, Variables >::updateBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
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
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      if( ( r_relative, inletOrientation ) <= 0.f )
         view_inletMark[ i ] = 0;

      return view_inletMark[ i ];
   };
   const GlobalIndexType numberOfNotRetyped = Algorithms::reduce< DeviceType >( 0, numberOfBufferParticles, moveBufferParticles, TNL::Plus() );
   const GlobalIndexType numberOfRetyped = numberOfBufferParticles - numberOfNotRetyped;

   //auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_inletMark[ i ]; };
   //auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
   //const GlobalIndexType numberOfRetyped = numberOfBufferParticles - Algorithms::reduce< DeviceType >( 0, view_inletMark.getSize(), fetch, reduction, 0.0 ); //I like zeros baceause sort.
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
         //view_rho_buffer[ i ] = inletConstDensity; //TODO: Keep the density same.
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfRetyped, createNewFluidParticles );


   //fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
   fluid->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() + numberOfRetyped );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() + numberOfRetyped );

   //------EXERIMENT-//update//profile-------------------------------------------
   view_v_buffer = inletConstVelocity;

   auto bufferVelocityProfile = [=] __cuda_callable__ ( int i ) mutable
   {
      const VectorType r = view_r_buffer[ i ];

      const RealType eps = 0.0001;
      if( ( r[ 1 ] > ( 0.002  - eps ) ) && ( r[ 1 ] < ( 0.002  + eps ) ) ){
         view_v_buffer[ i ][ 0 ] = 0.3523f;
         //view_rho_buffer[ i ] = 1001.2;
      }
      else if( ( r[ 1 ] > ( 0.004  - eps ) ) && ( r[ 1 ] < ( 0.004  + eps ) ) ){
         view_v_buffer[ i ][ 0 ] = 0.6325f;
         //view_rho_buffer[ i ] = 1001.0;
      }
      else if( ( r[ 1 ] > ( 0.006  - eps ) ) && ( r[ 1 ] < ( 0.006  + eps ) ) ){
         view_v_buffer[ i ][ 0 ] = 0.8f;
         //view_rho_buffer[ i ] = 1000.68;
      }
      else if( ( r[ 1 ] > ( 0.008  - eps ) ) && ( r[ 1 ] < ( 0.008  + eps ) ) )
         view_v_buffer[ i ][ 0 ] = 0.9f;
         //view_rho_buffer[ i ] = 1000.68;

   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfBufferParticles, bufferVelocityProfile );



   //added
   std::cout << "... InletBuffer - system updated." << std::endl;
}

template< typename ModelPointer, typename SPHFluidConfig, typename Variables >
template< typename FluidPointer, typename OpenBoundaryPointer >
void
VerletIntegrator< ModelPointer, SPHFluidConfig, Variables >::updateOutletBuffer( RealType dt, FluidPointer& fluid, OpenBoundaryPointer& openBoundary )
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

   //check for particles leaving the buffer
   auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
      const VectorType r = view_r_buffer[ i ];
      const VectorType r_relative = bufferPosition - r;

      //buffer particle leaves the entire computation space
      if( ( r_relative, inletOrientation ) > bufferWidth[ 0 ] )
         view_inletMark[ i ] = 1;

      return view_inletMark[ i ];
   };
   const GlobalIndexType removeFromBufferCount = Algorithms::reduce< DeviceType >(
         0, numberOfBufferParticles, moveBufferParticles );

   //sort particles by mark
   //TODO: replace this with something
   thrust::sort_by_key( thrust::device,
                        view_inletMark.getArrayData(),
                        view_inletMark.getArrayData() + numberOfBufferParticles,
                        thrust::make_zip_iterator( thrust::make_tuple(
                              view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );

  // Due to sort, this can be simply ignored
  // auto discardBufferParticles = [=] __cuda_callable__ ( int i ) mutable
  // {
  //       view_r_buffer[ i ] = FLT_MAX;
  //       view_rho_buffer[ i ] = 0.f;
  //       view_v_buffer[ i ] = 0.f;
  // };
  // Algorithms::parallelFor< DeviceType >(
  //       numberOfBufferParticles - removeFromBufferCount , numberOfBufferParticles, discardBufferParticles );

   numberOfBufferParticles = numberOfBufferParticles - removeFromBufferCount;
   openBoundary->particles->setNumberOfParticles( numberOfBufferParticles - removeFromBufferCount );
   openBoundary->particles->setNumberOfParticles( openBoundary->particles->getLastActiveParticle() - removeFromBufferCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() - removeFromBufferCount );

   //----- FLUID TO BUFFER -------------------------------------------------
   //Load all fluid fields
   auto view_r_fluid = fluid->particles->getPoints().getView();
   auto view_v_fluid = fluid->variables->v.getView();
   auto view_rho_fluid = fluid->variables->rho.getView();
   auto view_rho_old = fluid->integratorVariables->rho_old.getView();
   auto view_v_old = fluid->integratorVariables->v_old.getView();

   //TODO: Ugly ugly ugly temp workaround.
   const typename SPHFluidTraitsType::IndexVectorType gridIndex = TNL::floor( ( view_r_buffer.getElement( 0 ) - openBoundary->particles->getGridOrigin() ) / openBoundary->particles->getSearchRadius() );
   const GlobalIndexType gridColumnAuxTrick = gridIndex[ 0 ];

   const PairIndexType particleRangeToCheck = fluid->particles->getFirstLastParticleInColumnOfCells( gridColumnAuxTrick );
   auto receivingParticleMark_view = openBoundary->variables->receivingParticleMark.getView();

   //obtain fluid particles that should be retyped to buffer particles
   auto checkFluidParticles = [=] __cuda_callable__ ( int i ) mutable
   {
      const VectorType r = view_r_fluid[ i ];
      const VectorType r_relative = bufferPosition - r;

      //fluid particle should be retyped to the buffer particle
      if( ( r_relative, inletOrientation ) > 0 ){
         receivingParticleMark_view[ i - particleRangeToCheck[ 0 ] ] = i;
         return 1;
      }
      return 0;
   };
   const GlobalIndexType fluidToBufferCount = Algorithms::reduce< DeviceType >(
         0, numberOfParticle, checkFluidParticles, TNL::Plus() );

   thrust::sort( thrust::device,
                 receivingParticleMark_view.getArrayData(),
                 receivingParticleMark_view.getArrayData() + numberOfParticle );

   //retype fluid particles to buffer particles
   auto retypeFluidToOutlet = [=] __cuda_callable__ ( int i ) mutable
   {
      const GlobalIndexType p = receivingParticleMark_view[ i ];

      view_r_buffer[ numberOfBufferParticles + i ] = view_r_fluid[ p ];
      view_rho_buffer[ numberOfBufferParticles + i ] = view_rho_fluid[ p ];
      view_v_buffer[ numberOfBufferParticles + i ] = view_v_fluid[ p ];

      //shift the particles to keep compact set TODO: This can be avoided by just particle sort if FLT_MAX.
      view_r_fluid[ p ] = FLT_MAX;
      swap( view_r_fluid[ p ], view_r_fluid[ numberOfParticle - i - 1 ] );
      swap( view_rho_fluid[ p ], view_rho_fluid[ numberOfParticle - i -1 ] );
      swap( view_v_fluid[ p ], view_v_fluid[ numberOfParticle - i -1 ] );
      swap( view_rho_old[ p ], view_rho_old[ numberOfParticle - i -1 ] );
      swap( view_v_old[ p ], view_v_old[ numberOfParticle - i -1 ] );
   };
   Algorithms::parallelFor< DeviceType >( 0, fluidToBufferCount, retypeFluidToOutlet );

   //update particle ranges for buffer and fluid
   openBoundary->particles->setNumberOfParticles( numberOfBufferParticles + fluidToBufferCount );
   openBoundary->particles->setNumberOfParticles( openBoundary->particles->getLastActiveParticle() + fluidToBufferCount );
   openBoundary->setLastActiveParticle( openBoundary->getLastActiveParticle() + fluidToBufferCount );

   fluid->particles->setNumberOfParticles( numberOfParticle - fluidToBufferCount );
   fluid->particles->setLastActiveParticle( fluid->particles->getLastActiveParticle() - fluidToBufferCount );
   fluid->setLastActiveParticle( fluid->getLastActiveParticle() - fluidToBufferCount );
}



} // SPH
} // ParticleSystem
} // TNL

