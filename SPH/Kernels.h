#pragma once

#include "SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

//template< typename SPHCaseConfig >
class WendlandKernel //2D
{
public:
   __cuda_callable__
   static double F( double r, double h )
   {
      double F = 0;
      const double awen = float( 0.557 / ( h * h ) ); //2D
      const double bwen = float( -2.7852 / ( h * h * h ) ); //2D
      const double qq = r / h; //2D
      const float wqq1 = 1.f - 0.5f * qq;
      const float wqq2 = wqq1 * wqq1;
      const float wqq = qq + qq + 1.f;

      if ( qq <= 2.0 ){
         F = bwen * qq * wqq2 * wqq1 / r;
      }
      else{
         F = 0.0;
      }
      return F;
   }

   __cuda_callable__
   static double W( double r, double h )
   {
      double W = 0;

      const double awen = float( 0.557 / ( h * h ) ); //2D
      const double bwen = float( -2.7852 / ( h * h * h) ); //2D
      const double qq = r / h;
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

