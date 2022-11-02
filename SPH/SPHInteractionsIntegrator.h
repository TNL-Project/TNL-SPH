#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "SPHFluidVariables.h"
#include "SPHFluidTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig > >
class VerletIntegrator
{
  public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;
  using DeviceType = typename SPHFluidConfig::DeviceType;

  using RealType = typename Variables::RealType;
  using GlobalIndexType = typename Variables::GlobalIndexType;

  using ScalarArrayType = typename Variables::ScalarArrayType;
  using VectorArrayType = typename Variables::VectorArrayType;

  using ScalarType = typename SPHFluidTraitsType::ScalarType;
  using VectorType = typename SPHFluidTraitsType::VectorType;
  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using ScalarArrayView = typename Containers::ArrayView< ScalarType >;
  using VectorArrayView = Containers::ArrayView< VectorType >;
  using InteractionResultView = Containers::ArrayView< InteractionResultType >;

  VerletIntegrator( GlobalIndexType size, Variables& variables_ref ) : rhoO( size ), vO( size ), rhoOO( size ), vOO( size ), variables( variables_ref )
  {

  }

  void IntegrateVerlet( GlobalIndexType numberOfParticles, RealType dt )
  {
    ScalarArrayView rhoO_view = rhoO.getView();
    VectorArrayView vO_view = vO.getView();

    ScalarArrayView rhoOO_view = rhoOO.getView();
    VectorArrayView vOO_view = vOO.getView();

    ScalarArrayView rho_view = variables.rho.getView();
    VectorArrayView v_view = variables.v.getView();
    VectorArrayView r_view = variables.r.getView();

    InteractionResultView DrhoDv_view = variables.DrhoDv.getView();

    RealType dtdt05 = 0.5 * dt * dt;
    RealType dt2 = 2 * dt;

    auto init = [=] __cuda_callable__ ( int i ) mutable
    {
       const RealType drho = { DrhoDv_view[ i ][ 0 ] };
       const VectorType a = { DrhoDv_view[ i ][ 1 ], DrhoDv_view[ i ][ 2 ] };

       r_view[ i ] += v_view[ i ] * dt + a * dtdt05;
       rho_view[ i ] = rhoOO_view + drho * dt2;
       vOO_view[ i ] = vOO_view + drho * dt2;
    };

    Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, init );

  }

  ScalarArrayType rhoO;
  VectorArrayType vO;

  ScalarArrayType rhoOO;
  VectorArrayType vOO;

  Variables& variables;

};

} // SPH
} // ParticleSystem
} // TNL
