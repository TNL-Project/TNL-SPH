#pragma once

#include "../../SPHTraits.h"
#include <thrust/gather.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig >
class SPHFluidVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   SPHFluidVariables( GlobalIndexType size )
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ), gamma( size ), rho_swap( size ), v_swap( size ) {}

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;
   ScalarArrayType gamma;

   //Additional variable fields to avoid inmpace sort
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
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType& firstActiveParticle )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 ); //TODO: Obvious.
   }

};

template< typename SPHState >
class SPHOpenBoundaryVariables : public SPHFluidVariables< SPHState >
{
   public:
   using BaseType = SPHFluidVariables< SPHState >;
   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   SPHOpenBoundaryVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHState >( size ), particleMark( size ), receivingParticleMark( size ) {};

   IndexArrayType particleMark;
   IndexArrayType receivingParticleMark;
};

template< typename SPHConfig >
class SPHBoundaryVariables
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   SPHBoundaryVariables( GlobalIndexType size )
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ), n( size ),
     rho_swap( size ), v_swap( size ), n_swap( size ) {}

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;
   VectorArrayType n;

   //Additional variable fields to avoid inmpace sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   VectorArrayType n_swap;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();
      auto view_n = n.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();
      auto view_n_swap = n_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData() + firstActiveParticle, view_rho_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData() + firstActiveParticle, view_v_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_n.getArrayData() + firstActiveParticle, view_n_swap.getArrayData() + firstActiveParticle );

      rho.swap( rho_swap );
      v.swap( v_swap );
      n.swap( n_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      reader.template readParticleVariable< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
      reader.template readParticleVariable2D< VectorArrayType, typename ScalarArrayType::ValueType >( n, "Normals" ); //FIXME!
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", 3, numberOfParticles ); //TODO: Obvious.
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType& firstActiveParticle )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }

};

} // SPH
} // ParticleSystem
} // TNL

