#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "SPHFluidTraits.h"
#include "SPHFluidVariables.h"

#include "SPHFluidVariables.h"
#include "../Particles/ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables = SPHFluidVariables< SPHFluidConfig> >
class WCSPH_DBC
{
public:

  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename Variables::SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;
  using InteractionResultType = typename SPHFluidTraitsType::InteractionResultType;

  using DeviceType = typename Particles::Device;
  using PointType = typename Particles::PointType;
  using PointArrayType = typename Particles::PointArrayType;

  //using PointArrayTypeView = Containers::ArrayView< PointType >; /*temp-test*/

  WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref ) : vars( size ), points( points_ref ), particles( particles_ref ) {};
  //WCSPH_DBC( GlobalIndexType size, PointArrayType& points_ref, Particles& particles_ref,  PointArrayTypeView points_view ) : vars( size ), points( points_ref ), particles( particles_ref ),
  //points_view( points_view ) {};

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType PerformParticleInteractionFF( GlobalIndexType i, GlobalIndexType j );

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType PerformParticleInteractionFB( GlobalIndexType i, GlobalIndexType j );

  template< typename SPHKernelFunction >
  __cuda_callable__
  InteractionResultType PerformParticleInteractionBF( GlobalIndexType i, GlobalIndexType j );

  void sortParticlesAndVariables();

  /* TEMP, TEST */
  void FillParticleTypeWithInts()
  {
     auto view_type = vars.type.getView();

     for( int i = 0; i < particles.getNumberOfParticles(); i++ )
        view_type[ i ] = i;
  }

  Variables vars;
  //PointArrayTypeView pointss;
  PointArrayType& points;
  //PointArrayType_ptr points;

  Particles& particles;

};

} // SPH
} // ParticleSystem
} // TNL

#include "SPHInteractions_impl.h"
