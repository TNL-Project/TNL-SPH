#pragma once

#include "ParticlesTraits.h"

using namespace TNL;
using namespace TNL::ParticleSystem;

template < typename ParticleConfig, typename DeviceType >
class SimpleCellIndex
{
public:
   using ParticleTraitsType = ParticlesTraits< ParticleConfig, DeviceType >;
   using CellIndexView = TNL::Containers::ArrayView< typename ParticleTraitsType::CellIndexType, DeviceType >;
   using PointsView = TNL::Containers::ArrayView< typename ParticleTraitsType::PointType, DeviceType >;

   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;

   static constexpr GlobalIndexType _numberOfCells = ParticleConfig::gridXsize;

   static void ComputeCellIndex( CellIndexView cells, PointsView points )
   {
      auto f = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j ) mutable
      {
         cells[ j*_numberOfCells + i ] = j*_numberOfCells + i;
      };
      Algorithms::ParallelFor2D< DeviceType >::exec(
         ( LocalIndexType ) 0,
         ( LocalIndexType ) 0,
         ( LocalIndexType ) ParticleConfig::gridXsize,
         ( LocalIndexType ) ParticleConfig::gridYsize,
         f );
   }

   static void ComputeParticleCellIndex( CellIndexView view_particeCellIndices, PointsView view_points, GlobalIndexType _numberOfParticles )
   {

      auto f = [=] __cuda_callable__ ( LocalIndexType i ) mutable
      {
         //is necessary to norm particle coordinates to cell index, moreover 2D,3D,...
         view_particeCellIndices[ i ] = TNL::floor((view_points[ i ][ 0 ] - ParticleConfig::gridXbegin)/ ParticleConfig::searchRadius) + TNL::floor((view_points[ i ][ 1 ] - ParticleConfig::gridYbegin)/ ParticleConfig::searchRadius)*_numberOfCells;
      };
      Algorithms::ParallelFor< DeviceType >::exec(
         ( LocalIndexType ) 0,
          _numberOfParticles,
         f );
      }

};

