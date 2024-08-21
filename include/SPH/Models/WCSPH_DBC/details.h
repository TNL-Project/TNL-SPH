#pragma once

namespace TNL {
namespace SPH {

template< typename Matrix, typename RealType, typename VectorType >
__cuda_callable__
static Matrix
matrixCorrection2D( RealType W, VectorType gradW, VectorType r_ij, RealType V )
{
   Matrix A = 0.f;

   A( 0, 0 ) = W * V;
   A( 0, 1 ) = r_ij[ 0 ] * W * V;
   A( 0, 2 ) = r_ij[ 1 ] * W * V;

   A( 1, 0 ) = gradW[ 0 ] * V;
   A( 1, 1 ) = r_ij[ 0 ] * gradW[ 0 ] * V;
   A( 1, 2 ) = r_ij[ 1 ] * gradW[ 0 ] * V;

   A( 2, 0 ) = gradW[ 1 ] * V;
   A( 2, 1 ) = r_ij[ 0 ] * gradW[ 1 ] * V;
   A( 2, 2 ) = r_ij[ 1 ] * gradW[ 1 ] * V;

   return A;
}

template< typename Matrix, typename RealType, typename VectorType >
__cuda_callable__
static Matrix
matrixCorrection3D( RealType W, VectorType gradW, VectorType r_ij, RealType V )
{
   Matrix A = 0.f;

   A( 0, 0 ) = W * V;
   A( 0, 1 ) = r_ij[ 0 ] * W * V;
   A( 0, 2 ) = r_ij[ 1 ] * W * V;
   A( 0, 3 ) = r_ij[ 2 ] * W * V;

   A( 1, 0 ) = gradW[ 0 ] * V;
   A( 1, 1 ) = r_ij[ 0 ] * gradW[ 0 ] * V;
   A( 1, 2 ) = r_ij[ 1 ] * gradW[ 0 ] * V;
   A( 1, 3 ) = r_ij[ 2 ] * gradW[ 0 ] * V;

   A( 2, 0 ) = gradW[ 1 ] * V;
   A( 2, 1 ) = r_ij[ 0 ] * gradW[ 1 ] * V;
   A( 2, 2 ) = r_ij[ 1 ] * gradW[ 1 ] * V;
   A( 2, 3 ) = r_ij[ 2 ] * gradW[ 1 ] * V;

   A( 3, 0 ) = gradW[ 2 ] * V;
   A( 3, 1 ) = r_ij[ 0 ] * gradW[ 2 ] * V;
   A( 3, 2 ) = r_ij[ 1 ] * gradW[ 2 ] * V;
   A( 3, 3 ) = r_ij[ 2 ] * gradW[ 2 ] * V;

   return A;
}

template< typename VectorExtendedType, typename RealType, typename VectorType, typename VariableType >
__cuda_callable__
static VectorExtendedType
getVariableValueAndGradient2D( RealType W, VectorType gradW, VariableType f, RealType V )
{
   VectorExtendedType f_gradf = 0.f;

   f_gradf[ 0 ] = f * W * V;
   f_gradf[ 1 ] = gradW[ 0 ] * f * V;
   f_gradf[ 2 ] = gradW[ 1 ] * f * V;

   return f_gradf;
}

template< typename VectorExtendedType, typename RealType, typename VectorType, typename VariableType >
__cuda_callable__
static VectorExtendedType
getVariableValueAndGradient3D( RealType W, VectorType gradW, VariableType f, RealType V )
{
   VectorExtendedType f_gradf = 0.f;

   f_gradf[ 0 ] = f * W * V;
   f_gradf[ 1 ] = gradW[ 0 ] * f * V;
   f_gradf[ 2 ] = gradW[ 1 ] * f * V;
   f_gradf[ 3 ] = gradW[ 2 ] * f * V;

   return f_gradf;
}

//template< typename VectorExtendedTyep, typename VectorType >
//__cuda_callable__
//static typename VectorType::RealType
//interpolateGhostNode2D( VectorExtendedType& variableValAndGrad, VectorType& dist )
//{
//   const RealType valInterpolated = variableValAndGrad[ 0 ] + \
//                                    variableValAndGrad[ 1 ] * dist[ 0 ] + \
//                                    variableValAndGrad[ 2 ] * dist[ 1 ];
//   return valInterpolated;
//}

namespace ghostNodeDetail {

template< typename VectorExtendedType, typename VectorType >
__cuda_callable__
static typename VectorType::RealType
interpolateGhostNode( VectorExtendedType& variableValAndGrad, VectorType& dist )
{
   using RealType = typename VectorType::RealType;

   if constexpr( VectorType::getSize() == 2 ){
      const RealType valInterpolated = variableValAndGrad[ 0 ] + \
                                       variableValAndGrad[ 1 ] * dist[ 0 ] + \
                                       variableValAndGrad[ 2 ] * dist[ 1 ];
      return valInterpolated;
   }

   if constexpr( VectorType::getSize() == 3 ){
      const RealType valInterpolated = variableValAndGrad[ 0 ] + \
                                       variableValAndGrad[ 1 ] * dist[ 0 ] + \
                                       variableValAndGrad[ 2 ] * dist[ 1 ] + \
                                       variableValAndGrad[ 3 ] * dist[ 2 ];
      return valInterpolated;
   }
}

} // ghostNodeDetails

} // SPH
} // TNL

