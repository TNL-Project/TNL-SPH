#pragma once

#include "../../SPHTraits.h"
#include <TNL/Particles/details/thrustExecPolicySelector.h>
#include <thrust/gather.h>
#include "BoundaryConditionsTypes.h"

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

   //SPHFluidVariables() = default;

   //SPHFluidVariables( GlobalIndexType size )
   //: rho( size ), drho ( size ), p( size ), v( size ), a( size ), rho_swap( size ), v_swap( size ) {}

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
      //FIXME
      if constexpr( SPHConfig::spaceDimension == 2 )
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( v, "Velocity" );
      if constexpr( SPHConfig::spaceDimension == 3 )
         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( v, "Velocity" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 ); //FIXME This part of code is slow and wrong.
   }

#ifdef HAVE_MPI
   template< typename Synchronizer, typename DistributedParticlesPointer >
   void
   synchronizeVariables( Synchronizer& synchronizer,
                         DistributedParticlesPointer& distributedParticles )
   {
      synchronizer.synchronize( rho, distributedParticles );
      synchronizer.synchronize( v, distributedParticles );
   }
#endif

protected:

};

template< typename SPHState >
class OpenBoundaryVariables : public FluidVariables< SPHState >
{
   public:
   using BaseType = FluidVariables< SPHState >;
   using SPHTraitsType = typename BaseType::SPHTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

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

template< typename SPHState, typename Enable = void >
class BoundaryVariables : public FluidVariables< SPHState >
{};

template< typename SPHState >
class BoundaryVariables< SPHState,
                         typename std::enable_if_t< std::is_same_v< typename SPHState::BCType, WCSPH_BCTypes::DBC > > >
: public FluidVariables< SPHState >
{
public:
   using BaseType = FluidVariables< SPHState >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;

   //SPHBoundaryVariables( GlobalIndexType size )
   //: SPHFluidVariables< SPHState >( size ) {};
};

template< typename SPHState >
class BoundaryVariables< SPHState,
                         typename std::enable_if_t< std::is_same_v< typename SPHState::BCType, WCSPH_BCTypes::MDBC > > >
: public FluidVariables< SPHState >
{
public:
   using Base = FluidVariables< SPHState >;
   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayTypePointer = typename Base::IndexArrayTypePointer;
   using VectorExtendedArrayType = typename SPHTraitsType::VectorExtendedArrayType;
   using MatrixExtendedArrayType = typename SPHTraitsType::MatrixExtendedArrayType;

   void
   setSize( const GlobalIndexType& size )
   {
      Base::setSize( size );
      ghostNodes.setSize( size );
      ghostNodes_swap.setSize( size );
      n.setSize( size );
      n_swap.setSize( size );

      rhoGradRho_gn.setSize( size );
      cMatrix_gn.setSize( size );
   }

   VectorArrayType ghostNodes;
   VectorArrayType ghostNodes_swap;
   VectorArrayType n;
   VectorArrayType n_swap;

   VectorExtendedArrayType rhoGradRho_gn;
   MatrixExtendedArrayType cMatrix_gn;

   template< typename ParticlesPointer >
   void
   sortVariables( ParticlesPointer& particles )
   {
      Base::sortVariables( particles );
      particles->reorderArray( ghostNodes, ghostNodes_swap );
      particles->reorderArray( n, n_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      Base::readVariables( reader );
      //FIXME
      if constexpr( SPHConfig::spaceDimension == 2 ){
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( ghostNodes, "GhostNodes" ); //FIXME!
         reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" ); //FIXME!
      }
      if constexpr( SPHConfig::spaceDimension == 3 ){
         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( ghostNodes, "GhostNodes" ); //FIXME!
         reader.template readParticleVariable3D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( n, "Normals" ); //FIXME!
      }
   }
};

} // SPH
} // TNL

