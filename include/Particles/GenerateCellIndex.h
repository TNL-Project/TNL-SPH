#pragma once

#include "ParticlesTraits.h"
#include <bits/utility.h>

using namespace TNL;
using namespace TNL::ParticleSystem;

//note:
//row-major - std::index_sequence< 0, 1 >
//column-major - std::index_sequence< 1, 0 >

template< int Dimension >
struct DefaultPermutation;

template <>
struct DefaultPermutation< 2 >
{
   using type = std::index_sequence< 0, 1 >;
};

template <>
struct DefaultPermutation< 3 >
{
   using type = std::index_sequence< 0, 1, 2 >;
};

template< int Dimension, typename ParticleConfig, typename Permutation = typename DefaultPermutation< Dimension >::value >
class SimpleCellIndex
{};

template< typename ParticleConfig, typename Permutation  >
class SimpleCellIndex< 2, ParticleConfig, Permutation >
{
public:
   using DeviceType = typename ParticleConfig::DeviceType;
   using ParticleTraitsType = ParticlesTraits< ParticleConfig, DeviceType >;
   using CellIndexView = TNL::Containers::ArrayView< typename ParticleTraitsType::CellIndexType, DeviceType >;
   using PointsView = TNL::Containers::ArrayView< typename ParticleTraitsType::PointType, DeviceType >;

   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using RealType = typename ParticleTraitsType::RealType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using PointType = typename ParticleTraitsType::PointType;

   static void ComputeCellIndex( CellIndexView cells, PointsView points, IndexVectorType gridSize )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1 > > )
      {
         auto f = [=] __cuda_callable__ ( const IndexVectorType& i ) mutable
         {
            cells[ i[ 1 ] * gridSize[ 0 ] + i[ 0 ] ] = i[ 1 ] * gridSize[ 0 ] + i[ 0 ];
         };

         IndexVectorType begin{ 0, 0 };
         Algorithms::parallelFor< DeviceType >( begin, gridSize,  f );
      }

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 1, 0 > > )
      {
         auto f = [=] __cuda_callable__ ( const IndexVectorType& i ) mutable
         {
            cells[ i[ 0 ] * gridSize[ 1 ] + i[ 1 ] ] = i[ 0 ] * gridSize[ 1 ] + i[ 1 ];
         };

         IndexVectorType begin{ 0, 0 };
         Algorithms::parallelFor< DeviceType >( begin, gridSize,  f );
      }
   }

   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices,
                                         const PointsView view_points,
                                         const GlobalIndexType firstActiveParticle,
                                         const GlobalIndexType lastActiveParticle,
                                         const IndexVectorType gridSize,
                                         const PointType gridOrigin,
                                         const RealType searchRadius  )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1 > > )
      {
         auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) + \
                                           TNL::floor( ( view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridSize[ 0 ];
         };
         Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle + 1, f );
      }

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 1, 0 > > )
      {
         auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            /**
             * Due to rounding errors, there is difference between this form of expression and the actually used.
             * view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) + \
             *                              TNL::floor( ( view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) * gridSize[ 1 ];
             */
            view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 1 ] / searchRadius - gridOrigin[ 1 ] / searchRadius ) ) + \
                                           TNL::floor( ( view_points[ i ][ 0 ] / searchRadius - gridOrigin[ 0 ] / searchRadius ) ) * gridSize[ 1 ];
         };
         Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle + 1, f );
      }
   }

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const GlobalIndexType& i, const GlobalIndexType& j, const IndexVectorType& gridSize )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1 > > )
         return j * gridSize[ 0 ]  + i;

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 1, 0 > > )
         return i * gridSize[ 1 ]  + j;
   }

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const IndexVectorType& i, const IndexVectorType& gridSize )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1 > > )
         return i[ 1 ] * gridSize[ 0 ]  + i[ 0 ];

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 1, 0 > > )
         return i[ 0 ] * gridSize[ 1 ]  + i[ 1 ];
   }

};

template< typename ParticleConfig, typename Permutation >
class SimpleCellIndex< 3, ParticleConfig, Permutation >
{
public:
   using DeviceType = typename ParticleConfig::DeviceType;
   using ParticleTraitsType = ParticlesTraits< ParticleConfig, DeviceType >;
   using CellIndexView = TNL::Containers::ArrayView< typename ParticleTraitsType::CellIndexType, DeviceType >;
   using PointsView = TNL::Containers::ArrayView< typename ParticleTraitsType::PointType, DeviceType >;

   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using RealType = typename ParticleTraitsType::RealType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using PointType = typename ParticleTraitsType::PointType;


   static void ComputeCellIndex( CellIndexView cells, PointsView points, IndexVectorType gridSize  )
   {
      auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, LocalIndexType k ) mutable
      {
         cells[ k * ( gridSize[ 0 ] * gridSize[ 1 ] ) + j * gridSize[ 0 ] + i ] = k * ( gridSize[ 0 ] * gridSize[ 1 ] ) + j * gridSize[ 0 ] + i;
      };
      IndexVectorType begin{ 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, gridSize, f );
   }

   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices,
                                         const PointsView view_points,
                                         const GlobalIndexType firstActiveParticle,
                                         const GlobalIndexType lastActiveParticle,
                                         const IndexVectorType gridSize,
                                         const PointType gridOrigin,
                                         const RealType searchRadius  )
   {

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1, 2 > > )
      {
         auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) +
                                           TNL::floor( ( view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridSize[ 0 ] +
                                           TNL::floor( ( view_points[ i ][ 2 ] - gridOrigin[ 2 ] ) / searchRadius ) * gridSize[ 0 ] * gridSize[ 1 ];
         };
         Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle + 1, f );
      }

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 2, 1, 0 > > )
      {
         auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 2 ] - gridOrigin[ 2 ] ) / searchRadius ) +
                                           TNL::floor( ( view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridSize[ 1 ] +
                                           TNL::floor( ( view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) * gridSize[ 1 ] * gridSize[ 2 ];
         };
         Algorithms::parallelFor< DeviceType >( firstActiveParticle, lastActiveParticle + 1, f );
      }
   }

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const GlobalIndexType& i, const GlobalIndexType& j, const GlobalIndexType& k, const IndexVectorType& gridSize )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1, 2 > > )
         return k * gridSize[ 0 ] * gridSize[ 1 ] + j * gridSize[ 0 ] + i;

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 2, 1, 0 > > )
         return i * gridSize[ 1 ] * gridSize[ 2 ] + j * gridSize[ 1 ] + k;
   }

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const IndexVectorType& i, const IndexVectorType& gridSize )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1, 2 > > )
         return i[ 2 ] * gridSize[ 0 ] * gridSize[ 1 ] + i[ 1 ] * gridSize[ 0 ] + i[ 0 ];

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 2, 1, 0 > > )
         return i[ 0 ] * gridSize[ 1 ] * gridSize[ 2 ] + i[ 1 ] * gridSize[ 1 ] + i[ 2 ];
   }
};

