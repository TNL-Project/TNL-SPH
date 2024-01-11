#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <TNL/Endianness.h>
#include "particleReader.h"

namespace TNL {
namespace ParticleSystem {
namespace Readers {

class VTKReader : public ParticleReader
{
public:

   VTKReader() = default;

   VTKReader( const std::string& fileName ) : ParticleReader( fileName ) {}

   void
   detectParticleSystem() override
   {
      //reset();

      std::ifstream inputFile( fileName );
      if( ! inputFile )
        throw ParticleReaderError( "VTKReader", "failed to open the file '" + fileName +"'" );

      parseHeader( inputFile );

      if( dataset != "POLYDATA" )
         throw ParticleReaderError( "VTKReader", "the dataset '" + dataset + "' is not supported" );

      // parse the file, find the starting positions of all relevant sections
      findSections( inputFile );

      std::string line;
      std::string aux;
      std::istringstream iss;

      // parse points section
      if( sectionPositions.count( "POINTS" ) == 0 )
         throw ParticleReaderError( "VTKReader", "unable to find the POINTS section, the file may be invalid or corrupted" );
      inputFile.seekg( sectionPositions[ "POINTS" ] );
      if(getline( inputFile, line )) {}//{std::cout << "Reading done." << std::endl;}
      iss.clear();
      iss.str( line );
      iss >> aux;
      iss >> NumberOfPoints;
      iss >> pointsType;
      if( pointsType != "float" && pointsType != "double" )
         throw ParticleReaderError( "VTKReader", "unsupported data type for POINTS: " + pointsType );

      // only std::uint8_t makes sense for entity types
      typesType = "std::uint8_t";

      // arrays holding the data from the VTK file
      std::vector< double > pointsArray;
      //std::vector< std::int64_t > cellConnectivityArray;
      //std::vector< std::int64_t > cellOffsetsArray;
      std::vector< std::uint8_t > typesArray;

      // read points
      spaceDimension = 0;
      for( std::size_t pointIndex = 0; pointIndex < NumberOfPoints; pointIndex++ ) {
         if( ! inputFile )
            throw ParticleReaderError( "VTKReader", "unable to read enough vertices, the file may be invalid or corrupted" );

         // read the coordinates and compute the space dimension
         for( int i = 0; i < 3; i++ ) {
            double aux = 0;
            if( pointsType == "float" )
               aux = readValue< float >( dataFormat, inputFile );
            else
               aux = readValue< double >( dataFormat, inputFile );
            if( ! inputFile )
               throw ParticleReaderError( "VTKReader",
                                          "unable to read " + std::to_string( i ) + "th component of the vertex number "
                                          + std::to_string( pointIndex ) );
            if( aux != 0.0 )
               spaceDimension = std::max( spaceDimension, i + 1 );
            pointsArray.push_back( aux );
         }
      }

      // set the arrays to the base class
      this->pointsArray = std::move( pointsArray );

      // indicate success by setting the mesh type
      //meshType = "Meshes::Mesh";

   }

   VariantVector
   readPointData( const std::string& arrayName ) override
   {
      return readPointOrCellData( "POINT_DATA", arrayName );
   }

protected:
   // output of parseHeader
   std::string formatVersion;
   VTK::FileFormat dataFormat = VTK::FileFormat::ascii;
   std::string dataset;

   // output of findSection
   std::map< std::string, std::ios::pos_type > sectionPositions;

   // particleSystem properties - needed for reading POINT_DATA and CELL_DATA
   std::int32_t points_count = 0;
   std::int32_t cells_count = 0; //temp, remove

   void
   parseHeader( std::istream& str )
   {
     std::string line;
     std::istringstream iss;

     // check header
     getline( str, line );
     static const std::string prefix = "# vtk DataFile Version ";
     if( line.size() < prefix.size() )
        throw ParticleReaderError( "VTKReader", "failed to parse the VTK file header: unsupported VTK header '" + line + "'" );
     formatVersion = line.substr( prefix.length() );
     if( line.substr( 0, prefix.length() ) != prefix )
        throw ParticleReaderError( "VTKReader", "failed to parse the VTK file header: unsupported VTK header '" + line + "'" );
     if( formatVersion != "2.0" && formatVersion != "5.1" )
         throw ParticleReaderError( "VTKReader", "unsupported VTK DataFile Version: '" + formatVersion + "'" );

     // skip title
      if( ! str )
         throw ParticleReaderError( "VTKReader", "failed to parse the VTK file header" );
      getline( str, line );

     // parse data type
      if( ! str )
         throw ParticleReaderError( "VTKReader", "failed to parse the VTK file header" );
      std::string format;
      getline( str, format );
      if( format == "ASCII" )
         dataFormat = VTK::FileFormat::ascii;
      else if( format == "BINARY" )
         dataFormat = VTK::FileFormat::binary;
      else
         throw ParticleReaderError( "VTKReader", "unknown data format: '" + format + "'" );

      // parse dataset
      if( ! str )
         throw ParticleReaderError( "VTKReader", "failed to parse the VTK file header" );
      getline( str, line );
      iss.clear();
      iss.str( line );
      std::string tmp;
      iss >> tmp;
      if( tmp != "DATASET" )
         throw ParticleReaderError( "VTKReader", "wrong dataset specification: '" + line + "'" );
      iss >> dataset;

   }

   static void
   skip_meta( std::istream& str )
   {
      // skip possible metadata
      // https://vtk.org/doc/nightly/html/IOLegacyInformationFormat.html
      std::string line;
      while( true ) {
         getline( str, line );
         if( ! str )
            throw ParticleReaderError( "VTKReader", "failed to parse a METADATA section: is it terminated by a blank line?" );
         // strip whitespace
         line.erase( std::remove_if( line.begin(), line.end(), isspace ), line.end() );
         if( line.empty() )
            break;
      }
   }


   void
   findSections( std::istream& str )
   {
      while( str ) {
         // drop all whitespace (empty lines etc) before saving a position and reading a line

         str >> std::ws;
         if( str.eof() ){
            str.clear();                   //debug
            //str.seekg(0, std::ios::beg); //debug
            break;
         }

         // read a line which should contain the following section header
         const std::ios::pos_type currentPosition = str.tellg();
         std::string line;
         getline( str, line );
         if( ! str )
            throw ParticleReaderError( "VTKReader", "failed to parse sections of the VTK file" );

         // parse the section name
         std::istringstream iss( line );
         std::string name;
         iss >> name;

         if( name == "FIELD" ) {
            sectionPositions.insert( { "FIELD", currentPosition } );
            // parse the rest of the line: FIELD FieldData <count>
            std::string aux;
            int count = 0;
            iss >> aux >> count;
            // skip the FieldData arrays
            for( int i = 0; i < count; i++ ) {
               getline( str, line );
               iss.clear();
               iss.str( line );
               // <name> <components> <tuples> <datatype>
               std::int32_t components = 0;
               std::int32_t tuples = 0;
               std::string datatype;
               iss >> aux;
               if( aux == "METADATA" ) {
                  // skip metadata and read again
                  skip_meta( str );
                  i--;
                  continue;
               }
               iss >> components >> tuples >> datatype;
               if( ! iss )
                  throw ParticleReaderError( "VTKReader", "failed to extract FieldData information from line '" + line + "'" );
               // skip the points coordinates
               for( std::int32_t j = 0; j < components * tuples; j++ )
                  skipValue( dataFormat, str, datatype );
               // skip end of line (or any whitespace)
               str >> std::ws;
            }
         }
         else if( name == "POINTS" ) {
            sectionPositions.insert( { "POINTS", currentPosition } );
            // parse the rest of the line: POINTS <points_count> <datatype>
            std::string datatype;
            iss >> points_count >> datatype;
            // skip the values
            for( std::int32_t j = 0; j < 3 * points_count; j++ )
               skipValue( dataFormat, str, datatype );
            // skip end of line (or any whitespace)
            str >> std::ws;
         }
         // METADATA is a thing since version 5.1 of the file format (or something else newer than 2.0)
         else if( name == "METADATA" ) {
            sectionPositions.insert( { "METADATA", currentPosition } );
            skip_meta( str );
         }
         else if( name == "VERTICES" ) {
            sectionPositions.insert( { "VERTICES", currentPosition } );

            //OFFSETS follow VERTICES, skip them
            getline( str, line );
            iss.clear();
            iss.str( line );
            std::string aux;
            std::string datatype;
            iss >> aux >> datatype;

            if( aux != "OFFSETS" )
              throw ParticleReaderError( "VTKReader", "expected OFFSETS section, found '" + aux + "'" );

            //next, skip VERTICES, we don't need them for particles
            for( std::int32_t j = 0; j < 1 * points_count + 1; j++ )
               skipValue( dataFormat, str, "vtktypeint64" );

            str >> std::ws;

            //follows connectivity, we don't need that either
            str >> std::ws;
            getline( str, line );
            iss.clear();
            iss.str( line );
            iss >> aux >> datatype;

            if( aux != "CONNECTIVITY" )
               throw ParticleReaderError( "VTKReader", "expected CONNECTIVITY section, found '" + aux + "'" );

            for( std::int32_t j = 0; j < points_count; j++ )
               skipValue( dataFormat, str, datatype );

         }
         else if( name == "POINT_DATA" ) { //we can remove cell data
            if( points_count == 0 )
               throw ParticleReaderError( "VTKReader",
                                          "encountered a " + name
                                           + " section, but the particle system was not parsed yet "
                                             "(cells count = "
                                           + std::to_string( cells_count ) + ", points count = " + std::to_string( points_count )
                                           + ")" );

            while( str ) {
               // drop all whitespace (empty lines etc) before saving a position and reading a line
               //str >> std::ws; //debug, this should be on
               if( str.eof() )
                  break;

               // read a line which should contain the following array metadata
               const std::ios::pos_type currentPosition = str.tellg();
               std::string line;
               getline( str, line );
               if( ! str )
                  throw ParticleReaderError( "VTKReader", "failed to parse sections of the VTK file" );

               // parse the array type
               std::istringstream iss( line );
               std::string type;
               iss >> type;

               // handle switching between CELL_DATA and POINT_DATA
               if( ( name == "CELL_DATA" && type == "POINT_DATA" ) || ( name == "POINT_DATA" && type == "CELL_DATA" ) ) {
                  name = type;
                  continue;
               }

               const std::int32_t elements = ( name == "CELL_DATA" ) ? cells_count : points_count;

               // scalars: 1 value per cell/point
               // vectors: 3 values per cell/point
               // fields: arbitrary number of values per cell/point
               int values_per_element = 1;

               // additional metadata
               std::string array_name;
               std::string datatype;

               if( type == "SCALARS" ) {
                  // parse the rest of the line: SCALARS <array_name> <datatype>
                  iss >> array_name >> datatype;
                  sectionPositions.insert( { name + "::" + array_name, currentPosition } );
                  // skip the LOOKUP_TABLE line
                  getline( str, line );
               }
               else if( type == "VECTORS" ) {
                  values_per_element = 3;
                  // parse the rest of the line: VECTORS <array_name> <datatype>
                  iss >> array_name >> datatype;
                  sectionPositions.insert( { name + "::" + array_name, currentPosition } );
               }
               else if( type == "TENSORS" ) {
                  values_per_element = 9;
                  // parse the rest of the line: TENSORS <array_name> <datatype>
                  iss >> array_name >> datatype;
                  sectionPositions.insert( { name + "::" + array_name, currentPosition } );
               }
               else if( type == "FIELD" ) {
                  // parse the rest of the line: FIELD FieldData <count>
                  std::string aux;
                  int count = 0;
                  iss >> aux >> count;
                  // skip the FieldData arrays
                  for( int i = 0; i < count; i++ ) {
                     // drop all whitespace (empty lines etc) before saving a position and reading a line
                     str >> std::ws;
                     const std::ios::pos_type currentPosition = str.tellg();
                     getline( str, line );
                     iss.clear();
                     iss.str( line );
                     // <array_name> <components> <tuples> <datatype>
                     std::int32_t components = 0;
                     std::int32_t tuples = 0;
                     std::string datatype;
                     iss >> array_name;
                     if( array_name == "METADATA" ) {
                        // skip metadata and read again
                        skip_meta( str );
                        i--;
                        continue;
                     }
                     iss >> components >> tuples >> datatype;
                     if( ! iss )
                        throw ParticleReaderError( "VTKReader",
                                               "failed to extract FieldData information from line '" + line + "'" );
                     sectionPositions.insert( { name + "::" + array_name, currentPosition } );
                     // skip the points coordinates
                     for( std::int32_t j = 0; j < components * tuples; j++ )
                        skipValue( dataFormat, str, datatype );
                     // skip end of line (or any whitespace)
                     str >> std::ws;
                  }
                  continue;
               } // FIELD
               else if( type == "METADATA" ) {
                  skip_meta( str );
                  continue;
               }
               else {
                  std::cerr << "VTKReader: encountered an unsupported CELL_DATA array type: " << type
                            << ". Ignoring the rest of the file." << std::endl;
                  return;
               }

               // skip the values
               for( std::int32_t j = 0; j < elements * values_per_element; j++ )
                  skipValue( dataFormat, str, datatype );
               // skip end of line (or any whitespace)
               str >> std::ws;
            }
         } //pointdata
         else
            throw ParticleReaderError( "VTKReader",
                                   "parsing error: unexpected section start at byte " + std::to_string( currentPosition )
                                      + " (section name is '" + name + "')" );

      } //while loop
   } //find Sections

   VariantVector
   readPointOrCellData( std::string sectionName, const std::string& arrayName )
   {
      // NOTE: we must open the file in binary mode to prevent CR/CRLF conversions on Windows
      std::ifstream inputFile( fileName, std::ios_base::binary );
      if( ! inputFile )
         throw ParticleReaderError( "VTKReader", "failed to open the file '" + fileName + "'" );

      std::int32_t elements = ( sectionName == "CELL_DATA" ) ? cells_count : points_count;
      int values_per_element = 1;

      sectionName += "::" + arrayName;
      if( sectionPositions.count( sectionName ) == 0 )
         throw ParticleReaderError( "VTKReader", "array " + arrayName + " was not found in the CELL_DATA section" );
      inputFile.seekg( sectionPositions[ sectionName ] );

      // type: SCALARS, VECTORS, etc.
      // datatype: int, float, double
      std::string type;
      std::string datatype;

      // parse the metadata line
      std::string line;
      getline( inputFile, line );
      line = rstrip_cr( line );
      std::istringstream iss( line );
      iss >> type;

      // if the line starts with the array name, it must be a FIELD
      if( type == arrayName ) {
         // parse <array_name> <components> <tuples> <datatype>
         iss >> values_per_element >> elements >> datatype;
      }
      else {
         // parse the rest of the line: <type> <array_name> <datatype>
         std::string array_name;
         iss >> array_name >> datatype;
         if( type == "SCALARS" ) {
            values_per_element = 1;
            // skip the LOOKUP_TABLE line
            getline( inputFile, line );
         }
         else if( type == "VECTORS" )
            values_per_element = 3;
         else if( type == "TENSORS" )
            values_per_element = 9;
         else
            throw ParticleReaderError( "VTKReader", "requested array type " + type + " is not implemented in the reader" );
      }

      if( datatype == "int" )
         return readDataArray< std::int32_t >( inputFile, elements * values_per_element );
      else if( datatype == "float" )
         return readDataArray< float >( inputFile, elements * values_per_element );
      else if( datatype == "double" )
         return readDataArray< double >( inputFile, elements * values_per_element );
      else
         throw ParticleReaderError( "VTKReader", "found data type which is not implemented in the reader: " + datatype );
   }

   template< typename T >
   std::vector< T >
   readDataArray( std::istream& str, std::int32_t values )
   {
      std::vector< T > vector( values );
      for( std::int32_t i = 0; i < values; i++ )
         vector[ i ] = readValue< T >( dataFormat, str );
      return vector;
   }

   static inline std::string
   rstrip_cr( const std::string& string )
   {
      const auto end = string.find_last_not_of( '\r' );
      return string.substr( 0, end + 1 );
   }


   static void
   skipValue( VTK::FileFormat format, std::istream& str, const std::string& datatype )
   {
      if( datatype == "int" )  // implicit in vtk DataFile Version 2.0
         readValue< std::int32_t >( format, str );
      else if( datatype == "vtktypeint32" )  // vtk DataFile Version 5.1
         readValue< std::int32_t >( format, str );
      else if( datatype == "vtktypeint64" )  // vtk DataFile Version 5.1
         readValue< std::int64_t >( format, str );
      else if( datatype == "unsigned_int" )
         readValue< unsigned int >( format, str );
      else if( datatype == "float" )
         readValue< float >( format, str );
      else if( datatype == "double" )
         readValue< double >( format, str );
      else if( datatype == "unsigned_char" )
         readValue< unsigned int >( format, str );
      else
         throw ParticleReaderError( "VTKReader", "found data type which is not implemented in the reader: " + datatype );
   }

   template< typename T >
   static T
   readValue( VTK::FileFormat format, std::istream& str )
   {
      T value;
      if( format == VTK::FileFormat::binary ) {
         str.read( reinterpret_cast< char* >( &value ), sizeof( T ) );
         // forceBigEndian = swapIfLittleEndian, i.e. here it forces a big-endian
         // value to the correct system endianness
         value = forceBigEndian( value );
      }
      else {
         str >> value;
      }
      return value;
   }


};

} // Readers
} // ParticleSystem
} // TNL

