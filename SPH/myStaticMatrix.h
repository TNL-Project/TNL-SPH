// Copyright (c) 2004-2023 Tomáš Oberhuber et al.
//
// This file is part of TNL - Template Numerical Library (https://tnl-project.org/)
//
// SPDX-License-Identifier: MIT

// Implemented by: Jakub Klinkovsky

#pragma once

#include <iomanip>
#include <TNL/Containers/NDArray.h>
#include <TNL/Containers/StaticVector.h>

//test
#include <TNL/Containers/detail/VectorAssignment.h>

namespace TNL {
namespace Matrices {

template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation = std::index_sequence< 0, 1 > >  // identity by default
class StaticMatrix : public Containers::StaticNDArray< Value,
                                                       // note that using std::size_t in SizesHolder does not make sense, since
                                                       // the StaticNDArray is based on StaticArray, which uses int as IndexType
                                                       Containers::SizesHolder< int, Rows, Columns >,
                                                       Permutation >
{
   using Base = Containers::StaticNDArray< Value, Containers::SizesHolder< int, Rows, Columns >, Permutation >;

public:
   //temp gargabe
   using IndexType = int;


   ////> inherit all constructors
   using Base::Base;

   // inherit all assignment operators
   using Base::operator=;

	////> add value based constructor
   /**
    * \brief Default constructor.
    */
   __cuda_callable__
   constexpr StaticMatrix() = default;

	template< typename T >
	__cuda_callable__
	//constexpr StaticMatrix( const T& v );
	StaticMatrix( const T& v );

	//... or just like this
	//__cuda_callable__
	//constexpr StaticMatrix( const Value& v )

   static constexpr std::size_t
   getRows()
   {
      return Rows;
   }

   __cuda_callable__
   static constexpr std::size_t
   getColumns()
   {
      return Columns;
   }

	////> add overloaded assign operator
	template< typename T >
	__cuda_callable__
	constexpr StaticMatrix&
	operator=( const T& v );

   __cuda_callable__
   Containers::StaticVector< Rows, Value >
   operator*( const Containers::StaticVector< Columns, Value >& vector ) const
   {
      Containers::StaticVector< Rows, Value > result;
      for( std::size_t i = 0; i < Rows; i++ ) {
         Value v = 0;
         for( std::size_t j = 0; j < Columns; j++ )
            v += ( *this )( i, j ) * vector[ j ];
         result[ i ] = v;
      }
      return result;
   }

   ////> add overload some other usefull operators
   __cuda_callable__
   StaticMatrix< Value, Rows, Columns >
   operator*( const Value& value ) const
   {
      StaticMatrix< Value, Rows, Columns > result;
      for( std::size_t i = 0; i < Rows; i++ ) {
         for( std::size_t j = 0; j < Columns; j++ )
            result( i, j ) = ( *this )( i, j ) * value;
      }
      return result;
   }

    __cuda_callable__
   constexpr StaticMatrix< Value, Rows, Columns >&
   operator*=( const Value& value )
   {
      for( std::size_t i = 0; i < Rows; i++ )
         for( std::size_t j = 0; j < Columns; j++ )
            ( *this )( i, j ) *= value;

      return *this;
   }

   //__cuda_callable__
   //StaticMatrix< Value, Rows, Columns >
   //operator+( const StaticMatrix< Value, Rows, Columns >& matrix ) const
   //{
   //   StaticMatrix< Value, Rows, Columns > result;
   //   for( std::size_t i = 0; i < Rows; i++ )
   //      for( std::size_t j = 0; j < Columns; j++ )
   //        result( i, j ) += matrix( i, j );

   //   return result;
   //}

   __cuda_callable__
   constexpr StaticMatrix< Value, Rows, Columns >&
   operator+=( const StaticMatrix< Value, Rows, Columns >& matrix )
   {
      for( std::size_t i = 0; i < Rows; i++ )
         for( std::size_t j = 0; j < Columns; j++ )
           ( *this )( i, j ) += matrix( i, j );

      return *this;
   }

   __cuda_callable__
   constexpr StaticMatrix< Value, Rows, Columns >&
   operator-=( const StaticMatrix< Value, Rows, Columns >& matrix )
   {
      for( std::size_t i = 0; i < Rows; i++ )
         for( std::size_t j = 0; j < Columns; j++ )
           ( *this )( i, j ) -= matrix( i, j );

      return *this;
   }

   void
   print( std::ostream& str ) const;

};

template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >
template< typename T >
__cuda_callable__
//constexpr StaticMatrix< Value, Rows, Columns, Permutation >::StaticMatrix( const T& v )
StaticMatrix< Value, Rows, Columns, Permutation >::StaticMatrix( const T& v )
{
	//This works but it complaints about calliign __host__ function from __host__ __device__ function.
	//this->setValue( v );
	this->array = v;
}

template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >  // identity by default
void
StaticMatrix< Value, Rows, Columns, Permutation >::print( std::ostream& str ) const
{
   for( IndexType row = 0; row < this->getRows(); row++ ) {
      str << "Row: " << row << " -> ";
      for( IndexType column = 0; column < this->getColumns(); column++ ) {
         std::stringstream str_;
         str_ << std::setw( 4 ) << std::right << column << ":" << std::setw( 4 ) << std::left
              << ( *this )( row, column );
         str << std::setw( 10 ) << str_.str();
      }
      if( row < this->getRows() - 1 )
         str << std::endl;
   }
}

template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >
std::ostream&
operator<<( std::ostream& str, const StaticMatrix< Value, Rows, Columns, Permutation >& matrix )
{
   matrix.print( str );
   return str;
}


template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >
__cuda_callable__
StaticMatrix< Value, Rows, Columns, Permutation >
operator+( StaticMatrix< Value, Rows, Columns, Permutation > a, const StaticMatrix< Value, Rows, Columns >& b )
{
   a +=b;
   return a;
}

template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >
__cuda_callable__
StaticMatrix< Value, Rows, Columns, Permutation >
operator-( StaticMatrix< Value, Rows, Columns, Permutation > a, const StaticMatrix< Value, Rows, Columns >& b )
{
   a -=b;
   return a;
}


template< typename Value,
          std::size_t Rows,
          std::size_t Columns,
          typename Permutation >
template< typename T >
__cuda_callable__
constexpr StaticMatrix< Value, Rows, Columns, Permutation >&
StaticMatrix< Value, Rows, Columns, Permutation >::operator=( const T& v )
{
	//This works but it complaints about calliign __host__ function from __host__ __device__ function.
	//this->setValue( v );
	this->array = v;
	return *this;
}

//Special operations for special cases:

template< typename Real >
__cuda_callable__
Real
Determinant( const StaticMatrix< Real, 3, 3 > &A )
{
	Real det;
   det = A( 0,  0 ) * A( 1, 1 ) * A( 2, 2 ) + \
			A( 0,  1 ) * A( 1, 2 ) * A( 2, 0 ) + \
    		A( 0,  2 ) * A( 1, 0 ) * A( 2, 1 ) - \
    		A( 2,  0 ) * A( 1, 1 ) * A( 0, 2 ) - \
    		A( 2,  1 ) * A( 1, 2 ) * A( 0, 0 ) - \
    		A( 2,  2 ) * A( 1, 0 ) * A( 0, 1 ) ;
	return det;
}

template< typename Real >
__cuda_callable__
StaticMatrix< Real, 3, 3 >
Inverse( const StaticMatrix< Real, 3, 3 > &A )
{
   Real det = Determinant( A );
   StaticMatrix< Real, 3, 3 > invA;

   invA( 0, 0 ) =    ( A( 1, 1 ) * A( 2, 2 ) - A( 1, 2 ) * A( 2, 1 ) ) / det,
   invA( 0, 1 ) =  - ( A( 0, 1 ) * A( 2, 2 ) - A( 0, 2 ) * A( 2, 1 ) ) / det,
   invA( 0, 2 ) =    ( A( 0, 1 ) * A( 1, 2 ) - A( 0, 2 ) * A( 1, 1 ) ) / det,
   invA( 1, 0 ) =  - ( A( 1, 0 ) * A( 2, 2 ) - A( 1, 2 ) * A( 2, 0 ) ) / det,
   invA( 1, 1 ) =    ( A( 0, 0 ) * A( 2, 2 ) - A( 0, 2 ) * A( 2, 0 ) ) / det,
   invA( 1, 2 ) =  - ( A( 0, 0 ) * A( 1, 2 ) - A( 0, 2 ) * A( 1, 0 ) ) / det,
   invA( 2, 0 ) =    ( A( 1, 0 ) * A( 2, 1 ) - A( 1, 1 ) * A( 2, 0 ) ) / det,
   invA( 2, 1 ) =  - ( A( 0, 0 ) * A( 2, 1 ) - A( 0, 1 ) * A( 2, 0 ) ) / det,
   invA( 2, 2 ) =    ( A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 ) ) / det;

   return invA;
}

template< typename Real >
__cuda_callable__
Containers::StaticVector< 3, Real >
Solve( const StaticMatrix< Real, 3, 3 > &A, const Containers::StaticVector< 3, Real > &b )
{
   return Inverse( A ) * b;
}

















}  // namespace Matrices
}  // namespace TNL

