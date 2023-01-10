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

   using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;

   VerletIntegrator( ModelPointer model, GlobalIndexType size, GlobalIndexType size_boundary, GlobalIndexType size_inlet )
   : rho_old( size ), v_old( size ), rho_old_swap( size ), v_old_swap( size ),
   rhoBoundary_old( size_boundary ), vBoundary_old( size_boundary ), rhoBoundary_old_swap( size_boundary ), vBoundary_old_swap( size_boundary ),
   rhoInlet_old( size_inlet ), vInlet_old( size_inlet ), rhoInlet_old_swap( size_inlet ), vInlet_old_swap( size_inlet ), inletMark( size_inlet ),
   model( model )
   {
      //rho_old = model->getFluidVariables().rho;
      //v_old = model->getFluidVariables().v;
      rho_old = 1000.f; rho_old_swap = 1000.f;
      v_old = 0.f; v_old_swap = 0.f;

      rhoBoundary_old = 1000.f; rhoBoundary_old_swap = 1000.f;
      vBoundary_old = 0.f; vBoundary_old_swap = 0.f;

      inletMark = 1;
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

   void
   updateBuffer( RealType dt, RealType bufferEdge )
   {
      const GlobalIndexType numberOfParticle = model->particles->getNumberOfParticles();
      const GlobalIndexType numberOfBufferParticles = model->inletParticles->getNumberOfParticles();

      //Buffer
      auto view_r_buffer = model->inletParticles->getPoints().getView();
      auto view_v_buffer = model->getInletVariables().v.getView();
      auto view_rho_buffer = this->model->getInletVariables().rho.getView();

      auto view_inletMark = inletMark.getView();
      view_inletMark = 1;

      //Fluid
      auto view_r_fluid = model->particles->getPoints().getView();
      auto view_v_fluid = model->getFluidVariables().v.getView();
      auto view_rho_fluid = model->getFluidVariables().rho.getView();

      auto view_rho_old = rho_old.getView();
      auto view_v_old = v_old.getView();


      //std::cout << "Buffer points:" << view_r_buffer << std::endl;

      //Move buffer particles
      auto moveBufferParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         view_r_buffer[ i ] += view_v_buffer[ i ] * dt;
         if( view_r_buffer[ i ][ 0 ] > bufferEdge )
         {
            view_inletMark[ i ] = 0;
         }
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, model->inletParticles->getNumberOfParticles(), moveBufferParticles );

      std::cout << "---> updateBuffer: Particle marked." << std::endl;


      //Get number of retyped particles
      auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_inletMark[ i ]; };
      auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      const GlobalIndexType numberOfRetyped = numberOfBufferParticles - Algorithms::reduce< DeviceType >( 0, view_inletMark.getSize(), fetch, reduction, 0.0 ); //I like zeros baceause sort.

      std::cout << ".................. numberOfRetyped: " << numberOfRetyped << std::endl;

      //Sort particles by mark
      thrust::sort_by_key( thrust::device, view_inletMark.getArrayData(), view_inletMark.getArrayData() + numberOfBufferParticles,
            thrust::make_zip_iterator( thrust::make_tuple( view_r_buffer.getArrayData(), view_v_buffer.getArrayData(), view_rho_buffer.getArrayData() ) ) );

      std::cout << "---> updateBuffer: Particle sorted." << std::endl;
      //std::cout << "Buffer points:" << view_r_buffer << std::endl;


      std::cout << ".................. numberOfParticles: " << model->particles->getNumberOfParticles() << std::endl;
      std::cout << "---> updateBuffer: Update number of particles in particle system." << std::endl;
      std::cout << ".................. numberOfAllocatedParticles: " << model->particles->getNumberOfAllocatedParticles() << std::endl;
      std::cout << ".................. numberOfParticles: " << model->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. fieldSize: " << view_r_fluid.getSize()  << " | " << view_v_fluid.getSize() << " | " << view_rho_fluid.getSize() << std::endl;

      //std::cout<< "density" << view_rho_buffer << std::endl;
      std::cout<< "velocity:" << view_v_buffer << std::endl;

      //Assign new fluid particles
      auto createNewFluidParticles = [=] __cuda_callable__ ( int i ) mutable
      {
         if( view_inletMark[ i ] == 0 )
         {
            //Create new fluid
            view_r_fluid[ numberOfParticle + i ] = view_r_buffer[ i ];
            //view_rho_fluid[ numberOfParticle + i ] = view_rho_buffer[ i ];
            //view_v_fluid[ numberOfParticle + i ] = view_v_buffer[ i ];

            view_rho_fluid[ numberOfParticle + i ] = 1000.f;
            view_rho_old[ numberOfParticle + i ] = 1000.f;
            view_v_fluid[ numberOfParticle + i ] = { 0.5f, 0.f };
            view_v_old[ numberOfParticle + i ] = { 0.5f, 0.f };

            printf( "___RETYPE___ i: %d, numberOfParticle + i: %d.\n", i, numberOfParticle + i );

            //Generate new bufffer particle
            // TODO: Generalize!
            const RealType bufferLength = 0.007f;
            const VectorType newBufferParticle = { view_r_buffer[ i ][ 0 ] - bufferLength, view_r_buffer[ i ][ 1 ] };
            view_r_buffer[ i ] = newBufferParticle;
            view_v_buffer[ i ] = { 0.5f, 0.f };
            view_rho_buffer[ i ] = 1000.f;

         }
      };
      Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfRetyped, createNewFluidParticles );

      //Change number of fluid particles
      model->particles->setNumberOfParticles( numberOfParticle + numberOfRetyped );
      std::cout << ".................. numberOfParticles: " << model->particles->getNumberOfParticles() << std::endl;
      std::cout << "---> updateBuffer: Update number of particles in particle system." << std::endl;
      std::cout << ".................. numberOfAllocatedParticles: " << model->particles->getNumberOfAllocatedParticles() << std::endl;
      std::cout << ".................. numberOfParticles: " << model->particles->getNumberOfParticles() << std::endl;
      std::cout << ".................. fieldSize: " << view_r_fluid.getSize()  << " | " << view_v_fluid.getSize() << " | " << view_rho_fluid.getSize() << std::endl;
      std::cout << ".................. fieldSize integrator: " << view_r_fluid.getSize()  << " | " << view_v_fluid.getSize() << " | " << view_rho_fluid.getSize() << std::endl;


      std::cout << "---> updateBuffer: Particle transfer done." << std::endl;

   }

protected:

   ModelPointer model;

   ScalarArrayType rho_old;
   VectorArrayType v_old;

   ScalarArrayType rhoBoundary_old;
   VectorArrayType vBoundary_old;

   ScalarArrayType rhoInlet_old;
   VectorArrayType vInlet_old;

   IndexArrayType inletMark;

#ifdef PREFER_SPEED_OVER_MEMORY
   ScalarArrayType rho_old_swap;
   VectorArrayType v_old_swap;

   ScalarArrayType rhoBoundary_old_swap;
   VectorArrayType vBoundary_old_swap;

   ScalarArrayType rhoInlet_old_swap;
   VectorArrayType vInlet_old_swap;
#endif
};

} // SPH
} // ParticleSystem
} // TNL

