#include <TNL/Containers/StaticVector.h>

#include "mfdConstants.h"
#include "TaylorMonomials.h"

namespace TNL {
namespace SPH {
namespace Interpolation {

// TODO: Pass this real type to ABFs directisimo!
template< int Dim, int Order, typename RealType, template< int, int, typename > typename ABFsType, typename KernelFunction, typename SolverConfig >
class MFD
{
   public:
   constexpr static int p = get_p< Dim, Order >();
   using VectorType = Containers::StaticVector< Dim, RealType >; //FIXME: Should not be there
   using BaseVectorType = Containers::StaticVector< p, RealType >;
   using BaseMatrixType = Matrices::StaticMatrix< RealType, p, p >;
   using ABFs = ABFsType< Dim, Order, SolverConfig >;
   using TaylorMonomials = TaylorMonomials< Dim, Order, RealType >;

   __cuda_callable__
   static const BaseMatrixType
   getPairCorrectionMatrix( VectorType r_ij, RealType h )
   {
      const BaseVectorType W_ji = ABFs::eval( ( -1 ) * r_ij, h ); //FIXME!
      const BaseVectorType X_ji = TaylorMonomials::eval(( -1 ) * r_ij, h ); //FIXME
      return Matrices::cross( X_ji, W_ji );
   }

   __cuda_callable__
   static const BaseVectorType
   getPairVariableAndDerivatives( VectorType r_ij, RealType h )
   {
      return ABFs::eval( r_ij, h );
   }
};

} //namespace Interpolation
} //namespace SPH
} //namespace TNL

