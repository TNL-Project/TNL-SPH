// Variables.h  –  AUTO-GENERATED
#pragma once

#include "../../SPHTraits.h"
#include <TNL/Particles/details/thrustExecPolicySelector.h>
#include <thrust/gather.h>
#include "../WCSPH_DBC/BoundaryConditionsTypes.h"

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

   ScalarArrayType rho;
   ScalarArrayType rho_swap;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType v_swap;
   VectorArrayType a;

   void setSize( const GlobalIndexType& size )
   {
      rho.setSize( size );
      rho_swap.setSize( size );
      drho.setSize( size );
      p.setSize( size );
      v.setSize( size );
      v_swap.setSize( size );
      a.setSize( size );
   }

   template< typename ParticlesPointer >
   void sortVariables( ParticlesPointer& particles )
   {
      particles->reorderArray( rho, rho_swap );
      particles->reorderArray( v, v_swap );
   }

   template< typename ReaderType >
   void readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      if constexpr( SPHConfig::spaceDimension == 2 )
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( v, "Velocity" );
      if constexpr( SPHConfig::spaceDimension == 3 )
         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( v, "Velocity" );
   }

   template< typename WriterType >
   void writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }

#ifdef HAVE_MPI
   template< typename Synchronizer, typename DistributedParticlesPointer >
   void synchronizeVariables( Synchronizer& synchronizer,
                             DistributedParticlesPointer& distributedParticles )
   {
      synchronizer.synchronize( rho, distributedParticles );
      synchronizer.synchronize( v, distributedParticles );
   }
#endif

};

template< typename SPHState >
class BoundaryVariables : public FluidVariables< SPHState >
{};

template< typename SPHState >
class OpenBoundaryVariables : public FluidVariables< SPHState >
{
public:
   using BaseType = FluidVariables< SPHState >;
   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   void setSize( const GlobalIndexType& size )
   {
      BaseType::setSize( size );
      particleMark.setSize( size );
      receivingParticleMark.setSize( size );
   }

   typename SPHTraitsType::IndexArrayType particleMark;
   typename SPHTraitsType::IndexArrayType receivingParticleMark;
};

} // SPH
} // TNL

// end Variables.h