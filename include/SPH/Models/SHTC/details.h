#pragma once

namespace TNL {
namespace SPH {
namespace details {

template< typename MatrixType  >
__cuda_callable__
static MatrixType
deviator2D( const MatrixType& A )
{
   const MatrixType eye( { 1, 0, 0, 1 } );
   return A - 1.f / 3.f * ( A( 0, 0 ) + A( 1, 1 ) ) * eye;
}

template< typename MatrixType  >
__cuda_callable__
static MatrixType
deviator3D( MatrixType& A )
{
   const MatrixType eye( { 1, 0, 0, 0, 1, 0, 0, 0, 1 } );
   return A - 1.f / 3.f * ( A( 0, 0 ) + A( 1, 1 ) + A( 2, 2 ) ) * eye;
}

template< typename MatrixType  >
__cuda_callable__
static MatrixType
deviator( const MatrixType& A )
{
   if constexpr( MatrixType::getRows() == 2 )
      return deviator2D( A );
   if constexpr( MatrixType::getRows() == 3 )
      return deviator3D( A );
}


template< typename MatrixType  >
__cuda_callable__
static MatrixType
unitMatrix()
{
   //TODO: Allow only for square matrices
   if constexpr( MatrixType::getRows() == 2 )
      return MatrixType( { 1, 0, 0, 1 } );
   if constexpr( MatrixType::getRows() == 3 )
      return MatrixType( { 1, 0, 0, 0, 1, 0, 0, 0, 1 } );
}

} // details
} // SPH
} // TNL

