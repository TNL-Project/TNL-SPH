#pragma once

#include "../../SPHTraits.h"
#include <TNL/Particles/details/thrustExecPolicySelector.h>
#include <thrust/gather.h>
#include "BoundaryConditionsTypes.h"

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
   using MarkerArrayType = typename SPHTraitsType::MarkerArrayType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;
   ScalarArrayType gamma;
   MarkerArrayType marker;

   GlobalIndexType highestReferentialIdx;
   IndexArrayType  referentialIdx;

   //Additional variable fields to avoid in-place sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   MarkerArrayType marker_swap;
   IndexArrayType  referentialIdx_swap;

   void
   setSize( const GlobalIndexType& size )
   {
      rho.setSize( size );
      drho.setSize( size );
      p.setSize( size );
      v.setSize( size );
      a.setSize( size );
      gamma.setSize( size );
      marker.setSize( size );
      rho_swap.setSize( size );
      v_swap.setSize( size );
      marker_swap.setSize( size );
      referentialIdx.setSize( size );
      referentialIdx_swap.setSize( size );

      referentialIdx.forAllElements( [] __cuda_callable__( GlobalIndexType i, GlobalIndexType& value ) { value = i; } );
      highestReferentialIdx = size;
   }

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      particles->reorderArray( rho, rho_swap );
      particles->reorderArray( v, v_swap );
      particles->reorderArray( marker, marker_swap );
      particles->reorderArray( referentialIdx, referentialIdx_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      reader.template readParticleVariable< MarkerArrayType, typename MarkerArrayType::ValueType >( marker, "Ptype" );
      //FIXME
      if constexpr( SPHConfig::spaceDimension == 2 )
         reader.template readParticleVariable2D< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
      if constexpr( SPHConfig::spaceDimension == 3 )
         reader.template readParticleVariable3D< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >(
         v, "Velocity", numberOfParticles, firstActiveParticle, 3 );  //TODO: Obvious.
      writer.template writePointData< IndexArrayType >(
            referentialIdx, "ReferentialIndex", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( gamma, "Gamma", numberOfParticles, firstActiveParticle, 1 );
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
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   void
   setSize( const GlobalIndexType& size )
   {
      Base::setSize( size );
      n.setSize( size );
      n_swap.setSize( size );
      elementSize.setSize( size );
      elementSize_swap.setSize( size );
   }

   VectorArrayType n;
   VectorArrayType n_swap;
   ScalarArrayType elementSize;
   ScalarArrayType elementSize_swap;

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      Base::sortVariables( particles );
      particles->reorderArray( n, n_swap );
      particles->reorderArray( elementSize, elementSize_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      Base::readVariables( reader );
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( elementSize, "ElementSize" );
      //FIXME
      if constexpr( SPHConfig::spaceDimension == 2 )
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" );
      if constexpr( SPHConfig::spaceDimension == 3 )
         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" ); //FIXME!
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      Base::writeVariables( writer, numberOfParticles );
      writer.template writeVector< VectorArrayType, typename Base::RealType >(
            n, "Normals", numberOfParticles, firstActiveParticle, 3 );  //TODO: Obvious.
   }
};

template< typename SPHState >
class OpenBoundaryVariables : public BoundaryVariables< SPHState >
{
public:
   using BaseType = BoundaryVariables< SPHState >;
   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   //SPHOpenBoundaryVariables( GlobalIndexType size )
   //: SPHFluidVariables< SPHState >( size ), particleMark( size ), receivingParticleMark( size ) {};
   void
   setSize( const GlobalIndexType& size )
   {
      BaseType::setSize( size );
      particleMark.setSize( size );
      receivingParticleMark.setSize( size );
   }

   IndexArrayType particleMark;
   IndexArrayType receivingParticleMark;
};

}  //namespace SPH
}  //namespace TNL

