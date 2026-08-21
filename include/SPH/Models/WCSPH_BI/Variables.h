#pragma once

#include "../../SPHTraits.h"
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
   IndexArrayType referentialIdx;

   //Additional variable fields to avoid in-place sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   MarkerArrayType marker_swap;
   IndexArrayType referentialIdx_swap;

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

      referentialIdx.forAllElements(
         [] __cuda_callable__( GlobalIndexType i, GlobalIndexType & value )
         {
            value = i;
         } );
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
      reader.template readParticleVariable< ScalarArrayType >( rho, "Density" );
      reader.template readParticleVariable< MarkerArrayType >( marker, "Ptype" );
      reader.template readParticleVariable< VectorArrayType >( v, "Velocity" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure" );
      writer.template writePointData< ScalarArrayType >( rho, "Density" );
      writer.template writePointData< VectorArrayType >( v, "Velocity" );
      writer.template writePointData< IndexArrayType >( referentialIdx, "ReferentialIndex" );
      writer.template writePointData< ScalarArrayType >( gamma, "Gamma" );
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
      reader.template readParticleVariable< ScalarArrayType >( elementSize, "ElementSize" );
      reader.template readParticleVariable< VectorArrayType >( n, "Normals" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer )
   {
      Base::writeVariables( writer );
      writer.template writePointData< VectorArrayType >( n, "Normals" );
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

