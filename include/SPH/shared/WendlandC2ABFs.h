#include <TNL/Containers/StaticVector.h>
#include "./../Kernels.h"

namespace TNL {
namespace SPH {
namespace Interpolation {

template< int Dim, int Order, typename SolverConfig >
class WendlandC2ABFs
{
public:
    static constexpr int p = get_p< Dim, Order >();
    //FIXME:
    using RealType = typename SolverConfig::RealType;
    using WendlandKernel = KernelFunctions::WendlandKernel< SolverConfig >; //FIXME
    using VectorType = Containers::StaticVector< Dim, RealType >;
    using BaseVectorType = Containers::StaticVector< p, RealType >;


    __cuda_callable__
    static BaseVectorType
    eval( const VectorType& r_ij, RealType h )
    {
       const RealType drs = l2Norm( r_ij );
       const RealType drs2 = drs * drs;
       BaseVectorType result;

       if constexpr ( Dim == 2 && Order == 1 ) {
          const RealType W = WendlandKernel::W( drs, h );
          const RealType F = WendlandKernel::F( drs, h );

          result = { W,
                     r_ij[ 0 ] * F ,
                     r_ij[ 1 ] * F  };
       }

       else if constexpr ( Dim == 3 && Order == 1 ) {
          const RealType W = WendlandKernel::W( drs, h );
          const RealType F = WendlandKernel::F( drs, h );

          result = { W,
                     r_ij[ 0 ] * F,
                     r_ij[ 1 ] * F,
                     r_ij[ 2 ] * F };
       }

       else if constexpr ( Dim == 2 && Order == 2 ) {
          const RealType W = WendlandKernel::W( drs, h );
          const RealType F = WendlandKernel::F( drs, h );
          const RealType F2 = WendlandKernel::F2( drs, h );

          result = { W,
                     r_ij[ 0 ] * F,
                     r_ij[ 1 ] * F,
                     r_ij[ 0 ] * r_ij[ 0 ] * F2 / drs + ( r_ij[ 1 ] * r_ij[ 1 ] ) * F / drs2,
                     r_ij[ 0 ] * r_ij[ 1 ] * F2 / drs - ( r_ij[ 0 ] * r_ij[ 1 ] ) * F / drs2,
                     r_ij[ 1 ] * r_ij[ 1 ] * F2 / drs + ( r_ij[ 1 ] * r_ij[ 1 ] ) * F / drs2 };
       }
       else if constexpr ( Dim == 3 && Order == 2 ) {
          const RealType W = WendlandKernel::W( drs, h );
          const RealType F = WendlandKernel::F( drs, h );
          const RealType F2 = WendlandKernel::F2( drs, h );

          result = { W,
                     r_ij[ 0 ] * F,
                     r_ij[ 1 ] * F,
                     r_ij[ 2 ] * F,
                     r_ij[ 0 ] * r_ij[ 0 ] * F2 / drs + ( r_ij[ 1 ] * r_ij[ 1 ] + r_ij[ 2 ] * r_ij[ 2 ] ) * F / drs2,
                     r_ij[ 0 ] * r_ij[ 1 ] * F2 / drs - ( r_ij[ 0 ] * r_ij[ 1 ] ) * F / drs2,
                     r_ij[ 1 ] * r_ij[ 1 ] * F2 / drs + ( r_ij[ 2 ] * r_ij[ 2 ] + r_ij[ 0 ] * r_ij[ 0 ] ) * F / drs2,
                     r_ij[ 0 ] * r_ij[ 2 ] * F2 / drs - ( r_ij[ 0 ] * r_ij[ 2 ] ) * F / drs2,
                     r_ij[ 2 ] * r_ij[ 2 ] * F2 / drs + ( r_ij[ 0 ] * r_ij[ 0 ] + r_ij[ 1 ] * r_ij[ 1 ] ) * F / drs2,
                     r_ij[ 1 ] * r_ij[ 2 ] * F2 / drs - ( r_ij[ 1 ] * r_ij[ 2 ] ) * F / drs2 };
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

