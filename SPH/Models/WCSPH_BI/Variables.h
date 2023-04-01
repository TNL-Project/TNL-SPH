#pragma once

#include "../../SPHTraits.h"
#include <thrust/gather.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

//TODO: How should I use this?
template< typename SPHFluidConfig >
class SPHFluidConstants
{
   public:
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

   using RealType = typename SPHFluidTraitsType::RealType;
   using VectorType = typename SPHFluidTraitsType::VectorType;

};

template< typename SPHFluidConfig >
class SPHFluidVariables
{
   public:
   using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;

   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

   using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHFluidConfig::DeviceType >;

   SPHFluidVariables( GlobalIndexType size )
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ),
     rho_swap( size ), v_swap( size ) {}

   /* Variables - Fields */
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData(), view_rho_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData(), view_v_swap.getArrayData() );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

#ifdef PREFER_SPEED_OVER_MEMORY
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
#endif
};

template< typename SPHFluidConfig >
class SPHBoundaryVariables
{
  public:
  using SPHFluidTraitsType = SPHFluidTraits< SPHFluidConfig >;

  using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
  using RealType = typename SPHFluidTraitsType::RealType;

  using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
  using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

  using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;
  using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHFluidConfig::DeviceType >;

  SPHBoundaryVariables( GlobalIndexType size )
  : rho( size ), drho ( size ), p( size ), v( size ), a( size ), n( size ),
    rho_swap( size ), v_swap( size ), n_swap( size ) {}

  /* Variables - Fields */
  ScalarArrayType rho;
  ScalarArrayType drho;
  ScalarArrayType p;
  VectorArrayType v;
  VectorArrayType a;
  VectorArrayType n;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();
      auto view_n = n.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();
      auto view_n_swap = n_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData(), view_rho_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData(), view_v_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_n.getArrayData(), view_n_swap.getArrayData() );

      rho.swap( rho_swap );
      v.swap( v_swap );
      n.swap( n_swap );
   }

#ifdef PREFER_SPEED_OVER_MEMORY
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   VectorArrayType n_swap;
#endif

};

} // SPH
} // ParticleSystem
} // TNL

