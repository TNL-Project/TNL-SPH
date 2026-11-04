// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Devices/Host.h>
#include <cmath>
#include <algorithm>
#include <utility>

namespace TNL::SPH {

namespace detail {

template< typename VectorType >
inline bool
isInsideBox( const VectorType& point, const VectorType& boxOrigin, const VectorType& boxSize )
{
   bool isInside = true;
   for( int d = 0; d < VectorType::getSize(); d++ )
      if( ( point[ d ] < boxOrigin[ d ] ) || ( point[ d ] > ( boxOrigin[ d ] + boxSize[ d ] ) ) )
         isInside = false;
   return isInside;
}

}
template< int spaceDimension, typename VectorType, typename RealType >
TNL::Containers::Array< VectorType, TNL::Devices::Host >
computeExcludedDomainFaceNormals( const VectorType& faceMin,
                                  const VectorType& faceMax,
                                  const VectorType& domMin,
                                  const VectorType& domMax,
                                  const RealType tol,
                                  const bool flipNormals )
{
   VectorType excludedTmp[ 2 * spaceDimension ];
   int nExcluded = 0;
   for( int d = 0; d < spaceDimension; d++ )
      for( int side = 0; side < 2; side++ ) {
         if( std::fabs( ( side == 0 ? faceMin[ d ] : faceMax[ d ] ) - ( side == 0 ? domMin[ d ] : domMax[ d ] ) ) > tol )
            continue;
         VectorType normal = 0.f;
         normal[ d ] = ( side == 0 ) ? 1.f : -1.f;
         if( flipNormals )
            normal[ d ] = -normal[ d ];
         excludedTmp[ nExcluded++ ] = normal;
      }
   return TNL::Containers::Array< VectorType, TNL::Devices::Host >( excludedTmp, nExcluded );
}

/* Domain extent is the union of the subdomain boxes (origins and dimensions are
   in world coordinates); the topology's own global grid cannot be trusted here. */
template< int spaceDimension, typename Topology >
decltype( auto )
computeSubdomainUnionBox( const Topology& topology )
{
   using VectorType = std::decay_t< decltype( topology.getSubdomain( 0 ).grid.getOrigin() ) >;
   using RealType = typename Topology::RealType;
   const RealType globalCell = topology.getGlobalGrid().getSpaceSteps()[ 0 ];
   const int nSubdomains = topology.getNumberOfSubdomains();
   VectorType domMin = 0.f, domMax = 0.f;
   for( int sd = 0; sd < nSubdomains; sd++ ) {
      for( int d = 0; d < spaceDimension; d++ ) {
         const RealType lo = topology.getSubdomain( sd ).grid.getOrigin()[ d ];
         const RealType hi =
            lo + topology.getSubdomain( sd ).grid.getDimensions()[ d ] * topology.getSubdomain( sd ).refinementFactor * globalCell;
         if( sd == 0 || lo < domMin[ d ] ) domMin[ d ] = lo;
         if( sd == 0 || hi > domMax[ d ] ) domMax[ d ] = hi;
      }
   }
   return std::make_pair( domMin, domMax );
}

template< typename VectorType, typename RealType >
struct MassNodeFaceContext {
   VectorType faceMin, faceMax;
   RealType tol;
   bool flipNormals;
};

template< typename SubdomainDescriptor, typename RealType >
MassNodeFaceContext< std::decay_t< decltype( SubdomainDescriptor::grid.getOrigin() ) >, RealType >
computeMassNodeFaceContext( const SubdomainDescriptor& ownSd,
                            const SubdomainDescriptor& nbSd,
                            const RealType globalCell )
{
   using VectorType = std::decay_t< decltype( ownSd.grid.getOrigin() ) >;
   const RealType ownCell = ownSd.refinementFactor * globalCell;
   const RealType nbCell = nbSd.refinementFactor * globalCell;
   const RealType tol = (RealType) 1.5f * std::max( ownCell, nbCell );

   const VectorType ownMin = ownSd.grid.getOrigin();
   const VectorType ownSize = ownSd.grid.getDimensions() * ownCell;
   const VectorType nbMin = nbSd.grid.getOrigin();

   const auto& fineSd = ownCell <= nbCell ? ownSd : nbSd;
   const VectorType faceMin = fineSd.grid.getOrigin();
   const VectorType faceSize = fineSd.grid.getDimensions() * std::min( ownCell, nbCell );

   return { faceMin, faceMin + faceSize, tol, detail::isInsideBox( nbMin, ownMin, ownSize ) };
}

}  // namespace TNL::SPH
