#include <limits>

namespace TNL {
namespace ParticleSystem {
namespace detail {

//This was the original idea
template< int spaceDimension >
struct checkParticle
{};

template <>
struct checkParticle< 2 >
{
   template< typename PointType >
   static bool
   __cuda_callable__
   isValid( const PointType& point )
   {
      if( point[ 0 ] == FLT_MAX || point[ 1 ] == FLT_MAX )
         return false;
      else
         return true;
   }

   template< typename PointType >
   __cuda_callable__
   static bool
   isInvalid( const PointType& point )
   {
      if( point[ 0 ] == FLT_MAX || point[ 1 ] == FLT_MAX )
         return true;
      else
         return false;
   }
};

template <>
struct checkParticle< 3 >
{
   template< typename PointType >
   __cuda_callable__
   static bool
   isValid( const PointType& point )
   {
      if( point[ 0 ] == FLT_MAX || point[ 1 ] == FLT_MAX ||  point[ 2 ] == FLT_MAX )
         return false;
      else
         return true;
   }

   template< typename PointType >
   __cuda_callable__
   static bool
   isInvalid( const PointType& point )
   {
      if( point[ 0 ] == FLT_MAX || point[ 1 ] == FLT_MAX ||  point[ 1 ] == FLT_MAX )
         return true;
      else
         return false;
   }
};

// This is maybe much better solution
// FIXME: I would like to check the point directly using if constexpr for the dimension,
//        but since it is udes in lambda function, nvcc disaegree

template< typename PointType >
__cuda_callable__
static bool
isValid( PointType point )
{
   static constexpr int dim = PointType::getSize();
   return checkParticle< dim >::isValid( point );
}

template< typename PointType >
__cuda_callable__
static bool
isInvalid( PointType point )
{
   static constexpr int dim = PointType::getSize();
   return checkParticle< dim >::isInvalid( point );
}

} //detail
} //namespace Particles
} //namespace TNL

