#include <TNL/Containers/StaticVector.h>

namespace TNL {
namespace SPH {
namespace Interpolation {

template< int Dim, int Order, typename RealType >
class TaylorMonomials
{
public:
    static constexpr int p = get_p< Dim, Order >();
    using VectorType = Containers::StaticVector< Dim, RealType >;
    using BaseVectorType = Containers::StaticVector< p, RealType >;

    __cuda_callable__
    static BaseVectorType
    eval( const VectorType& r_ij, RealType h )
    {
       BaseVectorType result;
       if constexpr ( Dim == 2 && Order == 1 ) {
           result = { 1.f,
                      r_ij[ 0 ],
                      r_ij[ 1 ] };
       }

       else if constexpr ( Dim == 3 && Order == 1 ) {
           result = { 1.f,
                      r_ij[ 0 ],
                      r_ij[ 1 ],
                      r_ij[ 2 ] };
       }

       else if constexpr ( Dim == 2 && Order == 2 ) {
           result = { 1.f,
                      r_ij[ 0 ],
                      r_ij[ 1 ],
                      0.5f * r_ij[ 0 ] * r_ij[ 0 ],
                      r_ij[ 0 ] * r_ij[ 1 ],
                      0.5f * r_ij[ 1 ] * r_ij[ 1 ] };
       }
       else if constexpr ( Dim == 3 && Order == 2 ) {
           result = { 1.f,
                      r_ij[ 0 ],
                      r_ij[ 1 ],
                      r_ij[ 2 ],
                      0.5f * r_ij[ 0 ] * r_ij[ 0 ],
                      r_ij[ 0 ] * r_ij[ 1 ],
                      0.5f * r_ij[ 1 ] * r_ij[ 1 ],
                      r_ij[ 0 ] * r_ij[ 2 ],
                      0.5f * r_ij[ 2 ] * r_ij[ 2 ],
                      r_ij[ 1 ] * r_ij[ 2 ] };
       }
       else {
           static_assert( Dim >= 2 && Order >= 0, "Unsupported dimension or order" );
       }

       return result;
    }
};

}  //namespace Interpolation
}  //namespace SPH
}  //namespace TNL

