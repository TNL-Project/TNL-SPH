#include "TNL/Containers/Array.h"
#include <string>

namespace TNL {
namespace ParticleSystem {

template< typename ValueType, typename DeviceType, typename GlobalIndexType >
class ParticleField
{
public:

   using ArrayType = TNL::Containers::Array< ValueType, DeviceType, GlobalIndexType  >;

   ParticleField() = default;

   ParticleField( const GlobalIndexType& size, const std::string& name )
   : field( size ), field_swap( size ), fieldName( name ) {}

   void
   setName( const std::string& name )
   {
      this->fieldName = name;
   }

   std::string
   getName()
   {
      return this->name;
   }

   void
   setSize( const GlobalIndexType& size )
   {
      field.setSize( size );
      field_swap.setSize( size );
   }

   const GlobalIndexType
   getSize() const
   {
      return field.getSize();
   }

   void
   reorder( IndexArrayTypePointer& map, const GlobalIndexType& numberOfParticles )
   {
      auto view_map = map->getView();
      auto view_field = field.getView();
      auto view_field_swap = field_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_field.getArrayData(), view_field_swap.getArrayData() );
      field.swap( field_swap );
   }

   void
   reorder( IndexArrayTypePointer& map,
            const GlobalIndexType& numberOfParticles,
            const GlobalIndexType& firstActiveParticle )
   {
      auto view_map = map->getView();
      auto view_field = field.getView();
      auto view_field_swap = field_swap.getView();

      using ThrustDeviceType = TNL::Thrust::ThrustExecutionPolicy< typename SPHConfig::DeviceType >;
      ThrustDeviceType thrustDevice;
      thrust::gather( thrustDevice, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_field.getArrayData() + firstActiveParticle, view_field_swap.getArrayData() + firstActiveParticle );
      field.swap( field_swap );
   }

   template< typename ReaderType >
   void
   read( ReaderType& reader )
   {
      reader.template readParticleVariable< ArrayType, ValueType >( field, fieldName );
   }

   template< typename WriterType >
   void
   save( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      //writer.template writePointData< ArrayType >( field, fieldName, numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }

#ifdef HAVE_MPI
   template< typename Synchronizer, typename SimulationSubdomainInfo >
   void
   synchronizeVariables( Synchronizer& synchronizer, SimulationSubdomainInfo& subdomainInfo )
   {
      synchronizer.template synchronizeArray< ValueType >( field, field_swap, subdomainInfo, 1 );
   }

   void
   centerInMemory( const GlobalIndexType& firstActiveParticle,
                   const GlobalIndexType& shiftInMemory,
                   const GlobalIndexType& numberOfParticles )
   {
      utils::shiftArray( field, field_swap, firstActiveParticle, shiftInMemory, numberOfParticles );
   }
#endif

protected:
   ParticlesPointer particlesPointer;

   std::string fieldName;
   ArrayType field;
   ArrayType field_swap;

};

} //namespace ParticleSystem
} //namespace TNL

