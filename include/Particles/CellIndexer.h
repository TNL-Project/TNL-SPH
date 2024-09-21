#pragma once

#include "ParticlesTraits.h"
#include <cfloat>
#include <cstdint>

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
   using value = std::index_sequence< 0, 1 >;
};

template <>
struct DefaultPermutation< 3 >
{
   using value = std::index_sequence< 0, 1, 2 >;
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

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const PointType& r,
                      const PointType& gridOrigin,
                      const IndexVectorType& gridDimension,
                      const RealType& searchRadius )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1 > > )
         return TNL::floor( ( r[ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) + \
                TNL::floor( ( r[ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridDimension[ 0 ];

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 1, 0 > > )
         return TNL::floor( ( r[ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) + \
                TNL::floor( ( r[ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) * gridDimension[ 1 ];
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

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const PointType& r,
                      const PointType& gridOrigin,
                      const IndexVectorType& gridDimension,
                      const RealType& searchRadius )
   {
      if constexpr( std::is_same_v< Permutation, std::index_sequence< 0, 1, 2 > > )
         return TNL::floor( ( r[ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) +
                TNL::floor( ( r[ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridDimension[ 0 ] +
                TNL::floor( ( r[ 2 ] - gridOrigin[ 2 ] ) / searchRadius ) * gridDimension[ 0 ] * gridDimension[ 1 ];

      if constexpr( std::is_same_v< Permutation, std::index_sequence< 2, 1, 0 > > )
         return TNL::floor( ( r[ 2 ] - gridOrigin[ 2 ] ) / searchRadius ) +
                TNL::floor( ( r[ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridDimension[ 1 ] +
                TNL::floor( ( r[ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) * gridDimension[ 1 ] * gridDimension[ 2 ];
   }
};

