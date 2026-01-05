#pragma once

#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace PressureGradients {

template< typename SPHConfig >
class Symmetric
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   __cuda_callable__
   static VectorType
   grad_p( const RealType& p_i, const RealType& p_j, const VectorType& gradWV_j )
   {
      return ( p_i + p_j ) * gradWV_j;
   }

   __cuda_callable__
   static VectorType
   BI_grad_p( const RealType& p_i, const RealType& p_j, const VectorType& gradWS_k )
   {
      return ( p_i + p_j ) * gradWS_k;
   }
};

template< typename SPHConfig >
class TIC
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   __cuda_callable__
   static VectorType
   grad_p( const RealType& p_i, const RealType& p_j, const VectorType& gradWV_j )
   {
      const RealType grad_p_core = ( p_i >= 0 ) ? ( p_j + p_i ) : ( p_j - p_i );
      return grad_p_core * gradWV_j;
   }

   __cuda_callable__
   static VectorType
   BI_grad_p( const RealType& p_i, const RealType& p_j, const VectorType& gradWS_k )
   {
      const RealType grad_p_core = ( p_i >= 0 ) ? ( p_j + p_i ) : ( p_j - p_i );
      return grad_p_core * gradWS_k;
   }
};


} // Pressure gradient
} // SPH
} // TNL

