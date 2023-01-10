#pragma once

#include <limits>

#include "VTKWriter.h"

namespace TNL {
namespace ParticleSystem {
namespace Writers {

template< typename ParticleSystem >
VTKWriter< ParticleSystem >::VTKWriter( std::ostream& str, VTK::FileFormat format ) : str( str.rdbuf() ), format( format )
{
   if( format != VTK::FileFormat::ascii && format != VTK::FileFormat::binary )
      throw std::domain_error( "The Legacy VTK file formats support only ASCII and BINARY formats." );
}

template< typename ParticleSystem >
void
VTKWriter< ParticleSystem >::writeMetadata( int cycle, double time )
{
   if( ! headerWritten )
      writeHeader();

   int n_metadata = 0;
   if( cycle >= 0 )
      ++n_metadata;
   if( time >= 0 )
      ++n_metadata;
   if( n_metadata > 0 )
      str << "FIELD FieldData " << n_metadata << "\n";
   if( cycle >= 0 ) {
      str << "CYCLE 1 1 int\n";
      writeValue( format, str, cycle );
      str << "\n";
   }
   if( time >= 0 ) {
      str << "TIME 1 1 double\n";
      writeValue( format, str, time );
      str << "\n";
   }
}

template< typename ParticleSystem >
template< int EntityDimension >
void
VTKWriter< ParticleSystem >::writeParticles( const ParticleSystem& particles )
{
   if( ! headerWritten )
      writeHeader();
   //writePoints( particles );
   writePointsTemp( particles );

   const std::uint64_t verticesCount = particles.getNumberOfParticles();

   // legacy VTK files always have fixed integer width, even in the BINARY format
   // - DataFormat version 2.0: 32-bit
   // - DataFormat version 5.1: 64-bit (vtktypeint64)
   str << std::endl << "VERTICES " << verticesCount + 1 << " " << verticesCount << std::endl;
   str << "OFFSETS vtktypeint64" << std::endl;
   const int particlesCount = particles.getNumberOfParticles();
      for( int i = 0; i < particlesCount + 1; i++ ) {
         writeValue< signed long >( format, str, ( signed long ) i );
         if( format == VTK::FileFormat::ascii )
             str << "\n";
      }
   str << "CONNECTIVITY vtktypeint64" << std::endl;
   //detail::VTKMeshEntitiesWriter< Mesh, EntityDimension >::template writeConnectivity< std::int64_t >( mesh, str, format );
      for( int i = 0; i < particlesCount; i++ ) {
         writeValue< signed long >( format, str, ( signed long ) i );
         if( format == VTK::FileFormat::ascii )
             str << "\n";
      }

}

template< typename ParticleSystem >
template< typename Array >
void
VTKWriter< ParticleSystem >::writePointData( const Array& array, const std::string& name, const int numberOfComponents )
{
   if( array.getSize() / numberOfComponents != typename Array::IndexType( pointsCount ) )
      throw std::length_error( "Mismatched array size for POINT_DATA section: " + std::to_string( array.getSize() )
                               + " (there are " + std::to_string( pointsCount ) + " points in the file)" );

   // check that we won't start the section second time
   if( currentSection != VTK::DataType::PointData && cellDataArrays * pointDataArrays != 0 )
      throw std::logic_error( "The requested data section is not the current section and it has already been written." );
   currentSection = VTK::DataType::PointData;

   // start the appropriate section if necessary
   if( pointDataArrays == 0 )
      str << std::endl << "POINT_DATA " << pointsCount << std::endl;
   ++pointDataArrays;

   writeDataArray( array, name, numberOfComponents );
}


template< typename ParticleSystem >
template< typename Array >
void
VTKWriter< ParticleSystem >::writeDataArray( const Array& array, const std::string& name, const int numberOfComponents )
{
   // use a host buffer if direct access to the array elements is not possible
   if( std::is_same< typename Array::DeviceType, Devices::Cuda >::value ) {
      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;
      HostArray hostBuffer;
      hostBuffer = array;
      writeDataArray( hostBuffer, name, numberOfComponents );
      return;
   }

   if( numberOfComponents != 1 && numberOfComponents != 3 )
      throw std::logic_error( "Unsupported numberOfComponents parameter: " + std::to_string( numberOfComponents ) );

   // write DataArray header
   if( numberOfComponents == 1 ) {
      str << "SCALARS " << name << " " << getType< typename Array::ValueType >() << std::endl;
      str << "LOOKUP_TABLE default" << std::endl;
   }
   else {
      str << "VECTORS " << name << " " << getType< typename Array::ValueType >() << std::endl;
   }

   for( typename Array::IndexType i = 0; i < array.getSize(); i++ ) {
      writeValue( format, str, array[ i ] );
      if( format == VTK::FileFormat::ascii )
         str << "\n";
   }
}

template< typename ParticleSystem >
void
VTKWriter< ParticleSystem >::writePoints( const ParticleSystem& particles )
{
   pointsCount = particles.getNumberOfParticles();
   str << "POINTS " << pointsCount << " " << getType< typename ParticleSystem::RealType >() << std::endl;
   for( std::uint64_t i = 0; i < pointsCount; i++ ) {
      const auto& point = particles.getPoint( i );
      for( int j = 0; j < point.getSize(); j++ )
         writeValue( format, str, point[ j ] );
      // VTK needs zeros for unused dimensions
      for( int j = point.getSize(); j < 3; j++ )
         writeValue( format, str, (typename ParticleSystem::PointType::RealType) 0 );
      if( format == VTK::FileFormat::ascii )
         str << "\n";
   }
}

template< typename ParticleSystem >
void
VTKWriter< ParticleSystem >::writePointsTemp( const ParticleSystem& particles )
{
   //I need to use host buffer also for points,
   /* write data array */
   // use a host buffer if direct access to the array elements is not possible
   if( std::is_same< typename ParticleSystem::Device, Devices::Cuda >::value ) {
      using Array = typename ParticleSystem::PointArrayType;
      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;
     HostArray hostBuffer;
     hostBuffer = particles.getPoints();

      pointsCount = particles.getNumberOfParticles();
      str << "POINTS " << pointsCount << " " << getType< typename ParticleSystem::RealType >() << std::endl;
      for( std::uint64_t i = 0; i < pointsCount; i++ ) {
         const auto& point = hostBuffer.getElement( i );
         for( int j = 0; j < point.getSize(); j++ )
            writeValue( format, str, point[ j ] );
         // VTK needs zeros for unused dimensions
         for( int j = point.getSize(); j < 3; j++ )
            writeValue( format, str, (typename ParticleSystem::PointType::RealType) 0 );
         if( format == VTK::FileFormat::ascii )
            str << "\n";
      }
   }
   else
      writePoints( particles );
}

template< typename ParticleSystem >
void
VTKWriter< ParticleSystem >::writeHeader()
{
   str << "# vtk DataFile Version 5.1\n"
       << "TNL DATA\n"
       << ( ( format == VTK::FileFormat::ascii ) ? "ASCII\n" : "BINARY\n" ) << "DATASET POLYDATA\n";
   headerWritten = true;
}

template< typename ParticleSystem >
template< typename T >
void
VTKWriter< ParticleSystem >::writeValue( VTK::FileFormat format, std::ostream& str, T value )
{
   if( format == VTK::FileFormat::binary ) {
      value = forceBigEndian( value );
      str.write( reinterpret_cast< const char* >( &value ), sizeof( T ) );
   }
   else {
      // precision affects only floating-point types, not integers
      str.precision( std::numeric_limits< T >::digits10 );
      str << value << " ";
   }
}

/* TEMPORARY, AWFUL, AWFUL WAY HOW TO WRITE VECTORS */
template< typename ParticleSystem >
template< typename Array, typename Type >
void
VTKWriter< ParticleSystem >::writeVector( const Array& array, const std::string& name, const int numberOfComponents, const int numberOfElements )
{
   /* write point data */
   //: if( array.getSize() / numberOfComponents != typename Array::IndexType( pointsCount ) )
   //:    throw std::length_error( "Mismatched array size for POINT_DATA section: " + std::to_string( array.getSize() )
   //:                             + " (there are " + std::to_string( pointsCount ) + " points in the file)" );

   // check that we won't start the section second time
   if( currentSection != VTK::DataType::PointData && cellDataArrays * pointDataArrays != 0 )
      throw std::logic_error( "The requested data section is not the current section and it has already been written." );
   currentSection = VTK::DataType::PointData;

   // start the appropriate section if necessary
   if( pointDataArrays == 0 )
      str << std::endl << "POINT_DATA " << pointsCount << std::endl;
   ++pointDataArrays;

   /* write data array */
   // use a host buffer if direct access to the array elements is not possible
   if( std::is_same< typename Array::DeviceType, Devices::Cuda >::value ) {
      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;
      HostArray hostBuffer;
      hostBuffer = array;
      hostBuffer.resize( numberOfElements );

   /* write points */
      pointsCount = hostBuffer.getSize();
      str << "VECTORS " << name << " " << getType< Type >() << std::endl;
      for( std::uint64_t i = 0; i < pointsCount; i++ ) {
         const auto& point = array.getElement( i );
         for( int j = 0; j < point.getSize(); j++ )
            writeValue( format, str, point[ j ] );
         // VTK needs zeros for unused dimensions
         for( int j = point.getSize(); j < 3; j++ )
            writeValue( format, str, ( Type ) 0 );
         if( format == VTK::FileFormat::ascii )
            str << "\n";
      }
   }
}

/* TEMP FOR OPEN SYSTEM */
template< typename ParticleSystem >
template< typename Array >
void
VTKWriter< ParticleSystem >::writePointData( const Array& array, const std::string& name, const int numberOfElements, const int numberOfComponents )
{
   //if( array.getSize() / numberOfComponents != typename Array::IndexType( pointsCount ) )
   //   throw std::length_error( "Mismatched array size for POINT_DATA section: " + std::to_string( array.getSize() )
   //                            + " (there are " + std::to_string( pointsCount ) + " points in the file)" );

   // check that we won't start the section second time
   if( currentSection != VTK::DataType::PointData && cellDataArrays * pointDataArrays != 0 )
      throw std::logic_error( "The requested data section is not the current section and it has already been written." );
   currentSection = VTK::DataType::PointData;

   // start the appropriate section if necessary
   if( pointDataArrays == 0 )
      str << std::endl << "POINT_DATA " << pointsCount << std::endl;
   ++pointDataArrays;

   writeDataArray( array, name, numberOfElements, numberOfComponents );
}


template< typename ParticleSystem >
template< typename Array >
void
VTKWriter< ParticleSystem >::writeDataArray( const Array& array, const std::string& name, const int numberOfElements, const int numberOfComponents )
{
   // use a host buffer if direct access to the array elements is not possible
   if( std::is_same< typename Array::DeviceType, Devices::Cuda >::value ) {
      using HostArray = typename Array::
         template Self< std::remove_const_t< typename Array::ValueType >, Devices::Host, typename Array::IndexType >;
      HostArray hostBuffer;
      hostBuffer = array;
      hostBuffer.resize( numberOfElements );
      writeDataArray( hostBuffer, name, numberOfElements, numberOfComponents );
      return;
   }

   if( numberOfComponents != 1 && numberOfComponents != 3 )
      throw std::logic_error( "Unsupported numberOfComponents parameter: " + std::to_string( numberOfComponents ) );

   // write DataArray header
   if( numberOfComponents == 1 ) {
      str << "SCALARS " << name << " " << getType< typename Array::ValueType >() << std::endl;
      str << "LOOKUP_TABLE default" << std::endl;
   }
   else {
      str << "VECTORS " << name << " " << getType< typename Array::ValueType >() << std::endl;
   }

   for( typename Array::IndexType i = 0; i < array.getSize(); i++ ) {
      writeValue( format, str, array[ i ] );
      if( format == VTK::FileFormat::ascii )
         str << "\n";
   }
}

}  // namespace Writers
}  // namespace Meshes
}  // namespace TNL
