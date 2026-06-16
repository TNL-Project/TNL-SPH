#pragma once

namespace TNL {
namespace SPH {
namespace Interpolation {

constexpr int factorial( int n ) {
    if ( n <= 1 ) return 1.0;
    double result = 1.0;
    for ( int i = 2; i <= n; ++i ) {
        result *= i;
    }
    return result;
}

// p = sum from k=0 to m of (k + dim - 1)! / (k! * (dim - 1)!) - 1
template < int Dim, int Order >
constexpr int get_p() {
    static_assert( Dim >= 1 && Dim <= 3, "Dim must be 2 or 3." );
    static_assert( Order >= 0, "Order (m) must be non-negative," );

    int result = 0;
    for ( int k = 0; k <= Order; ++k ) {
        result += factorial( k + Dim - 1 ) / ( factorial( k ) * factorial( Dim - 1 ) );
    }
    return result; // add -1 to subtract the zeroth order term
}

}  //namespace Interpolation
}  //namespace SPH
}  //namespace TNL

