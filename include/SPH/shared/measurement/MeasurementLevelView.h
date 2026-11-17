// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Array.h>

#include <SPH/DecompositionTopology.h>
#include <SPH/SPHTraits.h>

namespace TNL::SPH {

template< typename SPHConfig >
struct MeasurementOwnerBox
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   VectorType minimum{};
   VectorType maximum{};
   RealType refinementFactor = 0;
};

template< typename SPHConfig, typename SPHSimulation, typename AuxiliaryFluidPointer >
struct MeasurementLevelView
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using OwnerBox = MeasurementOwnerBox< SPHConfig >;

   const FluidPointer* fluid = nullptr;
   const BoundaryPointer* boundary = nullptr;
   std::vector< const AuxiliaryFluidPointer* > ownPatches;
   bool hasMatchingPatchInterfaces = false;
   OwnerBox ownerBox;
   RealType h = 0;
   RealType mass = 0;
   int levelId = -1;

   MeasurementLevelView() = default;

   MeasurementLevelView(
      const FluidPointer& fluid,
      const BoundaryPointer& boundary,
      const std::vector< AuxiliaryFluidPointer >& allPatches,
      const std::vector< InterfaceDescriptor >& interfaces,
      const int ownIdx,
      const OwnerBox& ownerBox,
      const auto& modelParams )
   : fluid( &fluid ),
     boundary( &boundary ),
     hasMatchingPatchInterfaces( allPatches.size() == interfaces.size() ),
     ownerBox( ownerBox ),
     levelId( ownIdx )
   {
      const std::size_t patchCount = std::min( allPatches.size(), interfaces.size() );
      for( std::size_t patch = 0; patch < patchCount; ++patch )
         if( interfaces[ patch ].ownIdx == ownIdx )
            ownPatches.push_back( &allPatches[ patch ] );

      const RealType searchRadius = fluid->getParticles()->getSearchRadius();
      const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
      h = refinementFactor * modelParams.h;
      mass = std::pow( refinementFactor, SPHConfig::spaceDimension ) * modelParams.mass;
   }
};

template< typename VectorType, typename OwnerBoxArrayView >
__cuda_callable__
int
findOwningMeasurementLevel( const VectorType& point, const OwnerBoxArrayView& ownerBoxes )
{
   int owner = -1;
   for( int level = 0; level < ownerBoxes.getSize(); ++level ) {
      const auto& box = ownerBoxes[ level ];
      if( !( box.refinementFactor > 0 ) )
         continue;
      bool contains = true;
      for( int d = 0; d < VectorType::getSize(); ++d )
         contains &= point[ d ] > box.minimum[ d ] && point[ d ] <= box.maximum[ d ];
      if( contains && ( owner < 0 || box.refinementFactor < ownerBoxes[ owner ].refinementFactor ) )
         owner = level;
   }
   return owner;
}

template< typename SPHConfig, typename LevelView >
auto
makeMeasurementOwnerBoxArray( const std::vector< LevelView >& levels )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHTraitsType::DeviceType;
   using OwnerBox = MeasurementOwnerBox< SPHConfig >;
   Containers::Array< OwnerBox, DeviceType > ownerBoxes;
   ownerBoxes.setSize( levels.size() );
   for( std::size_t level = 0; level < levels.size(); ++level )
      ownerBoxes.setElement( level, levels[ level ].ownerBox );
   return ownerBoxes;
}

template< typename PointArray, typename OwnerBoxArray, typename OwnerIndexArray >
void
computeMeasurementOwners( const PointArray& points, const OwnerBoxArray& ownerBoxes, OwnerIndexArray& owners )
{
   using DeviceType = typename OwnerBoxArray::DeviceType;
   using IndexType = typename OwnerIndexArray::IndexType;
   if( owners.getSize() != points.getSize() )
      owners.setSize( points.getSize() );

   const auto pointsView = points.getConstView();
   const auto ownerBoxesView = ownerBoxes.getConstView();
   auto ownersView = owners.getView();
   auto assignOwner = [ = ] __cuda_callable__( IndexType i ) mutable
   {
      ownersView[ i ] = findOwningMeasurementLevel( pointsView[ i ], ownerBoxesView );
   };
   Algorithms::parallelFor< DeviceType >( 0, points.getSize(), assignOwner );
}

namespace detail {

template< typename RealType >
bool
isPositiveFinite( const RealType value )
{
   return std::isfinite( value ) && value > 0;
}

template< typename LevelView >
bool
isValidMeasurementLevel( const LevelView& level )
{
   if( level.levelId < 0 || level.fluid == nullptr || level.boundary == nullptr || ! level.hasMatchingPatchInterfaces )
      return false;
   if( ! isPositiveFinite( level.h ) || ! isPositiveFinite( level.mass )
       || ! isPositiveFinite( level.ownerBox.refinementFactor ) )
      return false;
   for( int d = 0; d < LevelView::SPHTraitsType::spaceDimension; ++d ) {
      if( ! std::isfinite( level.ownerBox.minimum[ d ] ) || ! std::isfinite( level.ownerBox.maximum[ d ] )
          || !( level.ownerBox.minimum[ d ] < level.ownerBox.maximum[ d ] ) )
         return false;
   }
   for( const auto* patch : level.ownPatches )
      if( patch == nullptr )
         return false;
   return true;
}

template< typename OwnerBox >
bool
hasPositiveVolumeOverlap( const OwnerBox& left, const OwnerBox& right )
{
   for( int d = 0; d < OwnerBox::SPHTraitsType::spaceDimension; ++d ) {
      const auto overlapMinimum = std::max( left.minimum[ d ], right.minimum[ d ] );
      const auto overlapMaximum = std::min( left.maximum[ d ], right.maximum[ d ] );
      if( !( overlapMinimum < overlapMaximum ) )
         return false;
   }
   return true;
}

}  // namespace detail

template< typename LevelView >
bool
validateMeasurementLevels( const std::vector< LevelView >& levels )
{
   if( levels.empty() )
      return false;
   for( const auto& level : levels )
      if( ! detail::isValidMeasurementLevel( level ) )
         return false;

   for( std::size_t left = 0; left < levels.size(); ++left ) {
      for( std::size_t right = left + 1; right < levels.size(); ++right ) {
         if( levels[ left ].levelId == levels[ right ].levelId )
            return false;
         if( levels[ left ].ownerBox.refinementFactor == levels[ right ].ownerBox.refinementFactor
             && detail::hasPositiveVolumeOverlap( levels[ left ].ownerBox, levels[ right ].ownerBox ) )
            return false;
      }
   }
   return true;
}

}  // namespace TNL::SPH
