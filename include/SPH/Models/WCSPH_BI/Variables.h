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
   ScalarArrayType gamma;

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
      gamma.setSize( size );
      rho_swap.setSize( size );
      v_swap.setSize( size );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();
      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
                      view_rho.getArrayData(), view_rho_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
                      view_v.getArrayData(), view_v_swap.getArrayData() );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrust::device,
                      view_map.getArrayData(),
                      view_map.getArrayData() + numberOfParticles,
                      view_rho.getArrayData() + firstActiveParticle,
                      view_rho_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device,
                      view_map.getArrayData(),
                      view_map.getArrayData() + numberOfParticles,
                      view_v.getArrayData() + firstActiveParticle,
                      view_v_swap.getArrayData() + firstActiveParticle );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
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
      //writer.template writePointData< ScalarArrayType >( gamma, "Gamma", numberOfParticles, firstActiveParticle, 1 );
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
   using IndexArrayTypePointer = typename Base::IndexArrayTypePointer;

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

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      Base::sortVariables( map, numberOfParticles );

      auto view_map = map->getView();
      auto view_n = n.getView();
      auto view_n_swap = n_swap.getView();
      auto view_elementSize = elementSize.getView();
      auto view_elementSize_swap = elementSize_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
                      view_n.getArrayData(), view_n_swap.getArrayData() );
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
                      view_elementSize.getArrayData(), view_elementSize_swap.getArrayData() );

      n.swap( n_swap );
      elementSize.swap( elementSize_swap );
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

