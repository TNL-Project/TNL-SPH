#pragma once

#include "SPHTraits.h"

namespace TNL {
namespace SPH {
namespace KernelFunctions {

template< typename SPHConfig, int Dimension = SPHConfig::spaceDimension >
class WendlandKernel
{};

template< typename SPHConfig >
class WendlandKernel< SPHConfig, 2 >
{
public:

   __cuda_callable__
   static float F( float r, float h )
   {
      const float wConst = -0.3482f / ( h * h * h * h ); // 7/(4*PI*16*h^4)
      const float q = r / h;
      return ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * wConst;
   }

   __cuda_callable__
   static float W( float r, float h )
   {
      const float wConst = 0.03482f / ( h * h ); // 7/(4*PI*16*h^2)*(5/8)
      const float q = r / h;
      return ( 1.f + 2.f * q) * ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * wConst;
   }

   __cuda_callable__
   static float F2( float r, float h )
   {
      const float wConst = -2.7852f /( h * h * h * h );
      const float q = r / h;
      return ( 1.f - 0.5f * q ) * ( 1.f - 0.5f * q ) * ( 2.f * q - 1.f ) * wConst;
   }
};

template< typename SPHConfig >
class WendlandKernel< SPHConfig, 3 >
{
public:

   __cuda_callable__
   static float F( float r, float h )
   {
      const float wConst =  -0.2611136f / ( h * h * h * h * h) ; // 21/(16*PI*h^5)*(5/8)
      const float q = r / h;
      return ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * wConst;
   }

   __cuda_callable__
   static float W( float r, float h )
   {
      const float wConst = 0.02611136f / ( h * h * h ); // 21/(16*PI*h^3)
      const float q = r / h;
      return ( 1.f + 2.f * q) * ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * ( 2.f - q ) * wConst;
   }
};

} // KernelFunctions
} // SPH
} // TNL

