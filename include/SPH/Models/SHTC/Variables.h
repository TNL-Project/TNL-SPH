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
class Empty
{};

template< typename SPHState >
class SHTCVariables
{
   public:
   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using MatrixArrayType = typename SPHTraitsType::MatrixArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   //Variables - Fields
   ScalarArrayType rho;
   ScalarArrayType drhodt;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType dvdt;
   MatrixArrayType A;
   MatrixArrayType dAdt;
   MatrixArrayType stress;

   //Additional variable fields to avoid inmpace sort
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;
   MatrixArrayType A_swap;

   void
   setSize( const GlobalIndexType& size )
   {
      rho.setSize( size );
      drhodt.setSize( size );
      p.setSize( size );
      v.setSize( size );
      dvdt.setSize( size );
      A.setSize( size );
      dAdt.setSize( size );
      stress.setSize( size );
      rho_swap.setSize( size );
      v_swap.setSize( size );
      A_swap.setSize( size );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();
      auto view_A = A.getView();
      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();
      auto view_A_swap = A_swap.getView();


      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData(), view_rho_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData(), view_v_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_A.getArrayData(), view_A_swap.getArrayData() );

      rho.swap( rho_swap );
      v.swap( v_swap );
      A.swap( A_swap );

      //particles->reorderArray( rho, rho_swap );
      //particles->reorderArray( v, v_swap );
      //particles->reorderArray( A, A_swap );
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
class SHTCBoundaryVariables : public SHTCVariables< SPHState >
{
public:
   using Base = SHTCVariables< SPHState >;
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

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      Base::sortVariables( map, numberOfParticles );

      auto view_map = map->getView();
      auto view_n = n.getView();
      auto view_n_swap = n_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_n.getArrayData(), view_n_swap.getArrayData() );

      n.swap( n_swap );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      Base::sortVariables( map, numberOfParticles, firstActiveParticle );

      auto view_map = map->getView();
      auto view_n = n.getView();
      auto view_n_swap = n_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_n.getArrayData() + firstActiveParticle, view_n_swap.getArrayData() + firstActiveParticle );

      n.swap( n_swap );
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

