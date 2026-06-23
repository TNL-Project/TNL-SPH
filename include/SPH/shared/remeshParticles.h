#pragma once

#include <TNL/Algorithms/reduce.h>
#include <TNL/Functional.h>
#include "../SPHTraits.h"
#include <cfloat>
#include <cmath>

namespace TNL {
namespace SPH {
namespace remeshing {

template< typename RealType >
__cuda_callable__
RealType
Lambda21( RealType s )
{
   const RealType abs_s = ( s < 0 ) ? -s : s;
   if( abs_s < 1.f )
      return 1.f - 2.5f * abs_s * abs_s + 1.5f * abs_s * abs_s * abs_s;
   else if( abs_s < 2.f )
      return 2.f - 4.f * abs_s + 2.5f * abs_s * abs_s - 0.5f * abs_s * abs_s * abs_s;
   else
      return 0.f;
}

template< typename SPHConfig, typename FluidPointer >
void
getParticleCoordsLimits( FluidPointer& fluid,
                         typename SPHFluidTraits< SPHConfig >::VectorType& lowerCorner,
                         typename SPHFluidTraits< SPHConfig >::VectorType& upperCorner )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHTraitsType::DeviceType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using VectorType = typename SPHTraitsType::VectorType;

   const GlobalIndexType n = fluid->getNumberOfParticles();
   auto view_points = fluid->getPoints().getView();

   for( int d = 0; d < VectorType::getSize(); d++ ) {
      auto fetch = [=] __cuda_callable__( GlobalIndexType i )
      {
         return view_points[ i ][ d ];
      };
      lowerCorner[ d ] = Algorithms::reduce< DeviceType >( 0, n, fetch, TNL::Min() );
      upperCorner[ d ] = Algorithms::reduce< DeviceType >( 0, n, fetch, TNL::Max() );
   }
}

template< typename SPHConfig, typename FluidPointer, typename ModelParams >
void
createParticlesOnGrid( FluidPointer& inputSet, FluidPointer& outputSet, ModelParams& modelParams )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHTraitsType::DeviceType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;

   VectorType lowerCorner = 0;
   VectorType upperCorner = 0;
   getParticleCoordsLimits< SPHConfig >( inputSet, lowerCorner, upperCorner );

   const RealType dp = modelParams.dp;
   const RealType rho0 = modelParams.rho0;

   IndexVectorType gridCounts = 0;
   GlobalIndexType totalParticles = 1;
   for( int d = 0; d < IndexVectorType::getSize(); d++ ) {
      gridCounts[ d ] = static_cast< GlobalIndexType >( std::floor( ( upperCorner[ d ] - lowerCorner[ d ] ) / dp ) ) + 1;
      totalParticles *= gridCounts[ d ];
   }

   const GlobalIndexType allocatedParticles = outputSet->getNumberOfAllocatedParticles();
   if( totalParticles > allocatedParticles ) {
      std::cerr << "Remeshing error: number of grid particles (" << totalParticles
                << ") exceeds allocated size (" << allocatedParticles << ")." << std::endl;
      return;
   }

   outputSet->getParticles()->setNumberOfParticles( totalParticles );

   auto view_r_out = outputSet->getParticles()->getPoints().getView();
   auto view_rho_out = outputSet->getVariables()->rho.getView();
   auto view_v_out = outputSet->getVariables()->v.getView();
   auto view_marker_out = outputSet->getVariables()->marker.getView();

   if constexpr( SPHConfig::spaceDimension == 2 ) {
      auto create = [ = ] __cuda_callable__( const IndexVectorType& idx ) mutable
      {
         const GlobalIndexType i = idx[ 1 ] * gridCounts[ 0 ] + idx[ 0 ];
         view_r_out[ i ] = lowerCorner + VectorType( idx[ 0 ] * dp, idx[ 1 ] * dp );
         view_rho_out[ i ] = rho0;
         view_v_out[ i ] = 0.f;
         view_marker_out[ i ] = 0;
      };
      IndexVectorType begin = 0;
      Algorithms::parallelFor< DeviceType >( begin, gridCounts, create );
   }
   else if constexpr( SPHConfig::spaceDimension == 3 ) {
      auto create = [ = ] __cuda_callable__( const IndexVectorType& idx ) mutable
      {
         const GlobalIndexType i = ( idx[ 2 ] * gridCounts[ 1 ] + idx[ 1 ] ) * gridCounts[ 0 ] + idx[ 0 ];
         view_r_out[ i ] = lowerCorner + VectorType( idx[ 0 ] * dp, idx[ 1 ] * dp, idx[ 2 ] * dp );
         view_rho_out[ i ] = rho0;
         view_v_out[ i ] = 0.f;
         view_marker_out[ i ] = 0;
      };
      IndexVectorType begin = 0;
      Algorithms::parallelFor< DeviceType >( begin, gridCounts, create );
   }
}


template< typename ParticlesType, typename FluidPointer, typename ModelParams >
void
interpolateParticlesToGrid( FluidPointer& inputSet, FluidPointer& outputSet, ModelParams& modelParams )
{
   using SPHConfig = typename ModelParams::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHTraitsType::DeviceType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using VectorType = typename SPHTraitsType::VectorType;
   using RealType = typename SPHTraitsType::RealType;

   inputSet->searchForNeighbors();

   const RealType searchRadius = inputSet->getParticles()->getSearchRadius();
   typename ParticlesType::NeighborsLoopParams searchInputSet( inputSet->getParticles() );
   const GlobalIndexType n_output = outputSet->getNumberOfParticles();

   const RealType dp = modelParams.dp;
   const RealType m = modelParams.mass;

   auto view_r_in = inputSet->getParticles()->getPoints().getView();
   auto view_rho_in = inputSet->getVariables()->rho.getView();
   auto view_v_in = inputSet->getVariables()->v.getView();

   auto view_r_out = outputSet->getParticles()->getPoints().getView();
   auto view_rho_out = outputSet->getVariables()->rho.getView();
   auto view_v_out = outputSet->getVariables()->v.getView();
   auto view_rho_old_out = outputSet->getIntegratorVariables()->rho_old.getView();
   auto view_v_old_out = outputSet->getIntegratorVariables()->v_old.getView();

   auto pairInterpolation = [ = ] __cuda_callable__(
         LocalIndexType i,
         LocalIndexType j,
         VectorType& r_i,
         RealType* rho_acc,
         VectorType* rho_v_acc,
         RealType* gamma_acc ) mutable
   {
      const VectorType r_j = view_r_in[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         RealType W = 1.f;
         for( int d = 0; d < VectorType::getSize(); d++ ) {
            const RealType diff = r_i[ d ] - r_j[ d ];
            const RealType s = ( diff < 0.f ) ? -diff : diff;
            W *= Lambda21( s / dp );
         }
         if( W > 0.f ) {
            const RealType rho_j = view_rho_in[ j ];
            const VectorType v_j = view_v_in[ j ];
            *rho_acc += W * rho_j * m;
            *rho_v_acc += W * rho_j * v_j * m;
            *gamma_acc += W * rho_j * m;
         }
      }
   };

   auto interpolate = [ = ] __cuda_callable__( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_r_out[ i ];

      RealType rho_acc = 0.f;
      VectorType rho_v_acc = 0.f;
      RealType gamma_acc = 0.f;

      ParticlesType::NeighborsLoopAnotherSet::exec(
            i, r_i, searchInputSet, pairInterpolation, &rho_acc, &rho_v_acc, &gamma_acc );

      if( gamma_acc > 0.01f ) {
         view_rho_out[ i ] = rho_acc / gamma_acc;
         view_v_out[ i ] = rho_v_acc / rho_acc;
         view_rho_old_out[ i ] = view_rho_out[ i ];
         view_v_old_out[ i ] = view_v_out[ i ];
         return 0;
      }
      else {
         view_r_out[ i ] = FLT_MAX;
         return 1;
      }
   };

   const GlobalIndexType n_discarded = Algorithms::reduce< DeviceType >(
         0, n_output, interpolate, TNL::Plus() );

   if( n_discarded > 0 )
      outputSet->getParticles()->setNumberOfParticlesToRemove(
            outputSet->getParticles()->getNumberOfParticlesToRemove() + n_discarded );
}

template< typename ParticlesType, typename SPHConfig, typename FluidPointer, typename ModelParams >
void
remeshParticles( FluidPointer& inputSet, FluidPointer& outputSet, ModelParams& modelParams )
{
   createParticlesOnGrid< SPHConfig >( inputSet, outputSet, modelParams );

   interpolateParticlesToGrid< ParticlesType >( inputSet, outputSet, modelParams );

   outputSet->searchForNeighbors();

   std::swap( inputSet, outputSet );
}

}  //namespace remeshing
}  //namespace SPH
}  //namespace TNL
