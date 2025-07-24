#pragma once

#include "../../SPHTraits.h"
#include <TNL/Particles/details/thrustExecPolicySelector.h>
#include <thrust/gather.h>

#ifdef HAVE_MPI
#include "../../shared/utils.h"
#endif

namespace TNL {
namespace SPH {

template< typename SPHState >
class FluidVariables
{
   public:
   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;

   //Additional variable fields to avoid inmpace sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;

   void
   setSize( const GlobalIndexType& size )
   {
      rho.setSize( size );
      drho.setSize( size );
      p.setSize( size );
      v.setSize( size );
      a.setSize( size );
      rho_swap.setSize( size );
      v_swap.setSize( size );
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      particles->reorderArray( rho, rho_swap );
      particles->reorderArray( v, v_swap );
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

template< typename SPHState >
class BoundaryVariables : public FluidVariables< SPHState >
{
public:
   using Base = FluidVariables< SPHState >;
   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayTypePointer = typename Base::IndexArrayTypePointer;

   void
   setSize( const GlobalIndexType& size )
   {
      Base::setSize( size );
      n.setSize( size );
      n_swap.setSize( size );
   }

   VectorArrayType n;
   VectorArrayType n_swap;

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      Base::sortVariables( particles );
      particles->reorderArray( n, n_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      Base::readVariables( reader );
      //FIXME
      if constexpr( SPHConfig::spaceDimension == 2 )
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" );
      if constexpr( SPHConfig::spaceDimension == 3 )
         reader.template readParticleVariable< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" );
   }
};


} // SPH
} // TNL

