#pragma once

#include "SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHCaseConfig >
class WendlandKernel2D //2D
{
public:

   __cuda_callable__
   static float F( float r, float h )
   {
      const float bwen = float( -2.7852 / ( h * h * h ) ); //2D move to const
      const float qq = r / h; //2D
      const float wqq1 = 1.f - 0.5f * qq;
      return bwen * qq * wqq1 * wqq1 * wqq1 / r;
      /*
      if ( qq <= 2.0 )
         return bwen * qq * wqq1 * wqq1 * wqq1 / r;
      else
         return 0.0;
      */
   }

   __cuda_callable__
   static float W( float r, float h )
   {
      float W = 0;

      const float awen = float( 0.557 / ( h * h ) ); //2D
      const float bwen = float( -2.7852 / ( h * h * h) ); //2D
      const float qq = r / h;
      const float wqq1 = 1.f - 0.5f * qq;
      const float wqq2 = wqq1 * wqq1;
      const float wqq = qq + qq + 1.f;
      if ( qq <= 2.0 ){
         W = awen * wqq * wqq2 * wqq2;
      }
      else{
         W = 0.0;
      }
      return W;
   }
};

class WendlandKernel3D //2D
{
public:

   __cuda_callable__
   static float F( float r, float h )
   {
      const float bwen = float( -2.08891 / ( h * h * h * h ) ); //2D move to const
      const float qq = r / h; //2D
      const float wqq1 = 1.f - 0.5f * qq;
      return bwen * qq * wqq1 * wqq1 * wqq1 / r;
      /*
      if ( qq <= 2.0 )
         return bwen * qq * wqq1 * wqq1 * wqq1 / r;
      else
         return 0.0;
      */
   }

   __cuda_callable__
   static float W( float r, float h )
   {
      float W = 0;

      const float awen = float( 0.41778 / ( h * h *  h ) ); //can be constexpr
      const float bwen = float( -2.08891 / ( h * h * h * h ) ); //2D
      const float qq = r / h;
      const float wqq1 = 1.f - 0.5f * qq;
      const float wqq2 = wqq1 * wqq1;
      const float wqq = qq + qq + 1.f;
      if ( qq <= 2.0 ){
         W = awen * wqq * wqq2 * wqq2;
      }
      else{
         W = 0.0;
      }
      return W;
   }
};

} // SPH
} // ParticleSystem
} // TNL

