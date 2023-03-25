#pragma once

#include "ParticlesTraits.h"

using namespace TNL;
using namespace TNL::ParticleSystem;

template< int Dimension, typename ParticleConfig >
class SimpleCellIndex
{};

template< typename ParticleConfig >
class SimpleCellIndex< 2, ParticleConfig >
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

   //static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize;

   static void ComputeCellIndex( CellIndexView cells, PointsView points, IndexVectorType gridSize )
   {
      auto f = [=] __cuda_callable__ ( const IndexVectorType& i ) mutable
      {
         cells[ i[ 1 ] * gridSize[ 0 ] + i[ 0 ] ] = i[ 1 ] * gridSize[ 0 ] + i[ 0 ];
      };

      IndexVectorType begin{ 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, gridSize,  f );
      //Algorithms::ParallelFor2D< DeviceType >::exec(
      //   ( LocalIndexType ) 0,
      //   ( LocalIndexType ) 0,
      //   ( LocalIndexType ) gridSize[ 0 ],
      //   ( LocalIndexType ) gridSize[ 1 ],
      //   f );
   }

   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices, const PointsView view_points, const GlobalIndexType _numberOfParticles, const IndexVectorType gridSize, const PointType gridOrigin, const RealType searchRadius  )
   {

      auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         //is necessary to norm particle coordinates to cell index, moreover 2D,3D,...
         view_particeCellIndices[ i ] = TNL::floor((view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius) + TNL::floor((view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius)* gridSize[ 0 ];
      };
      Algorithms::parallelFor< DeviceType >( 0, _numberOfParticles, f );
      }

   __cuda_callable__
   static uint32_t
   //EvaluateCellIndex( GlobalIndexType& i, GlobalIndexType& j )
   EvaluateCellIndex( const GlobalIndexType& i, const GlobalIndexType& j, const IndexVectorType& gridSize )
   {
      //return j*_numberOfCells + i;
      return j * gridSize[ 0 ]  + i;
   }

};

template< typename ParticleConfig >
class SimpleCellIndex< 3, ParticleConfig >
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

   //static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize;
   //static constexpr GlobalIndexType _numberOfCellsY = ParticleConfig::gridYsize;

   static void ComputeCellIndex( CellIndexView cells, PointsView points, IndexVectorType gridSize  )
   {
      auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, LocalIndexType k ) mutable
      {
         cells[ k * ( gridSize[ 0 ] * gridSize[ 1 ] ) + j * gridSize[ 0 ] + i ] = k * ( gridSize[ 0 ] * gridSize[ 1 ] ) + j * gridSize[ 0 ] + i;
      };
      IndexVectorType begin{ 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, gridSize, f );
   }

   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices, const PointsView view_points, const GlobalIndexType _numberOfParticles, const IndexVectorType gridSize, const PointType gridOrigin, const RealType searchRadius  )
   {

      auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         view_particeCellIndices[ i ] = TNL::floor( ( view_points[ i ][ 0 ] - gridOrigin[ 0 ] ) / searchRadius ) +
                                        TNL::floor( ( view_points[ i ][ 1 ] - gridOrigin[ 1 ] ) / searchRadius ) * gridSize[ 0 ] +
                                        TNL::floor( ( view_points[ i ][ 2 ] - gridOrigin[ 2 ] ) / searchRadius ) * gridSize[ 0 ] * gridSize[ 1 ];
      };
      Algorithms::parallelFor< DeviceType >( 0, _numberOfParticles, f );
   }

   __cuda_callable__
   static uint32_t
   EvaluateCellIndex( const GlobalIndexType& i, const GlobalIndexType& j, const GlobalIndexType& k, const IndexVectorType& gridSize )
   {
      //return k * _numberOfCells * _numberOfCellsY + j *_numberOfCells + i;
      return k * gridSize[ 0 ] * gridSize[ 1 ] + j * gridSize[ 0 ] + i;
   }

};

//template< typename ParticleConfig, typename DeviceType >
//class ZOrderCurve
//{
//public:
//   using ParticleTraitsType = ParticlesTraits< ParticleConfig, DeviceType >;
//   using CellIndexView = TNL::Containers::ArrayView< typename ParticleTraitsType::CellIndexType, DeviceType >;
//   using PointsView = TNL::Containers::ArrayView< typename ParticleTraitsType::PointType, DeviceType >;
//
//   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
//   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
//
//   static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize;
//
//   static void ComputeCellIndex( CellIndexView cells, PointsView points )
//   {
//      auto f = [=] __cuda_callable__ (  const IndexVectorType& i  ) mutable
//      {
//         //https://graphics.stanford.edu/~seander/bithacks.html
//         static const uint32_t MASKS[] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF };
//         static const uint32_t SHIFTS[] = { 1, 2, 4, 8 };
//
//         uint32_t x = i;
//         uint32_t y = j;
//
//         x = ( x | ( x << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//         x = ( x | ( x << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//         x = ( x | ( x << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//         x = ( x | ( x << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//         y = ( y | ( y << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//         y = ( y | ( y << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//         y = ( y | ( y << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//         y = ( y | ( y << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//         const uint32_t result = x | ( y << 1 );
//         cells[ i[ 1 ] * _numberOfCells + i[ 0 ] ] = result;
//      };
//      IndexVectorType begin{ 0, 0 };
//      Algorithms::parallelFor< DeviceType >( begin, gridSize,  f );
//      //Algorithms::ParallelFor2D< DeviceType >::exec(
//      //   ( LocalIndexType ) 0,
//      //   ( LocalIndexType ) 0,
//      //   ( LocalIndexType ) ParticleConfig::gridXsize,
//      //   ( LocalIndexType ) ParticleConfig::gridYsize,
//      //   f );
//   }
//
//   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices, PointsView view_points, GlobalIndexType _numberOfParticles )
//   {
//      auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
//      {
//         //https://graphics.stanford.edu/~seander/bithacks.html
//         static const uint32_t MASKS[] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF };
//         static const uint32_t SHIFTS[] = { 1, 2, 4, 8 };
//
//         //is necessary to norm particle coordinates to cell index, moreover 2D,3D,...
//         uint32_t x = TNL::floor( ( view_points[ i ][ 0 ] - ParticleConfig::gridXbegin )/ ParticleConfig::searchRadius );
//         uint32_t y = TNL::floor( ( view_points[ i ][ 1 ] - ParticleConfig::gridYbegin )/ ParticleConfig::searchRadius );
//
//         x = ( x | ( x << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//         x = ( x | ( x << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//         x = ( x | ( x << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//         x = ( x | ( x << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//         y = ( y | ( y << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//         y = ( y | ( y << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//         y = ( y | ( y << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//         y = ( y | ( y << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//         const uint32_t result = x | ( y << 1 );
//
//         view_particeCellIndices[ i ] = result;
//      };
//      Algorithms::ParallelFor< DeviceType >::exec(
//         ( LocalIndexType ) 0,
//          _numberOfParticles,
//         f );
//      }
//
//   __cuda_callable__
//   static uint32_t
//   EvaluateCellIndex( unsigned int i, unsigned int j )
//   {
//      ////https://graphics.stanford.edu/~seander/bithacks.html
//      static const uint32_t MASKS[] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF };
//      static const uint32_t SHIFTS[] = { 1, 2, 4, 8 };
//
//      uint32_t x = i;
//      uint32_t y = j;
//
//      x = ( x | ( x << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//      x = ( x | ( x << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//      x = ( x | ( x << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//      x = ( x | ( x << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//      y = ( y | ( y << SHIFTS[ 3 ] ) ) & MASKS[ 3 ];
//      y = ( y | ( y << SHIFTS[ 2 ] ) ) & MASKS[ 2 ];
//      y = ( y | ( y << SHIFTS[ 1 ] ) ) & MASKS[ 1 ];
//      y = ( y | ( y << SHIFTS[ 0 ] ) ) & MASKS[ 0 ];
//
//      const uint32_t result = x | ( y << 1 );
//      return result;
//   }
//
//};

