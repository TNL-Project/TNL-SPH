#pragma once

#include "../../SPHTraits.h"
#include <thrust/gather.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

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
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ), rho_swap( size ), v_swap( size ) {}

   /* Variables - Fields */
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;

   /* Additional variable fields to avoid inmpace sort. */
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData() + firstActiveParticle, view_rho_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData() + firstActiveParticle, view_v_swap.getArrayData() + firstActiveParticle );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      reader.template readParticleVariable< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }
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
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ), ghostNodes( size ),
     rho_swap( size ), v_swap( size ), ghostNodes_swap( size ) {}

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;
   VectorArrayType ghostNodes;

   //Additional variable fields to avoid inmpace sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   VectorArrayType ghostNodes_swap;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();
      auto view_ghostNodes = ghostNodes.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();
      auto view_ghostNodes_swap = ghostNodes_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData() + firstActiveParticle, view_rho_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData() + firstActiveParticle, view_v_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_ghostNodes.getArrayData() + firstActiveParticle, view_ghostNodes_swap.getArrayData() + firstActiveParticle );

      rho.swap( rho_swap );
      v.swap( v_swap );
      ghostNodes.swap( ghostNodes_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      reader.template readParticleVariable< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
      reader.template readParticleVariable2D< VectorArrayType, typename ScalarArrayType::ValueType >( ghostNodes, "GhostNodes" ); //FIXME!
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }
};

} // SPH
} // ParticleSystem
} // TNL

