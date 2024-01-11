#pragma once

#include <ostream>

#include <stdexcept>
#include <string>
#include <cstdint>

#include "../Readers/VTKTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace Writers {

template< typename ParticleSystem >
class VTKWriter
{
   static_assert( ParticleSystem::getParticleDimension() <= 3, "The VTK format supports only 1D, 2D and 3D meshes." );
   // TODO: check also space dimension when grids allow it
   //   static_assert( Mesh::getSpaceDimension() <= 3, "The VTK format supports only 1D, 2D and 3D meshes." );

public:

   VTKWriter() = delete;

   VTKWriter( std::ostream& str, VTK::FileFormat format = VTK::FileFormat::binary );

   void
   writeMetadata( std::int32_t cycle = -1, double time = -1 );

   template< int EntityDimension = ParticleSystem::getParticleDimension() >
   void
   writeParticles( const ParticleSystem& particles );

   template< typename Array >
   void
   writePointData( const Array& array, const std::string& name, int numberOfComponents = 1 );

   template< typename Array >
   void
   writeDataArray( const Array& array, const std::string& name, int numberOfComponents = 1 );

   /* TEMPORARY, AWFUL, AWFUL WAY HOW TO WRITE VECTORS */
   template< typename Array, typename Type >
   void
   writeVector( const Array& array, const std::string& name,  const int numberOfElementsm, const int writeFromElement, const int numberOfComponents );

   /* TEMP for opensystem */
   template< typename Array >
   void
   writePointData( const Array& array, const std::string& name,  const int numberOfElements, const int writeFromElement, int numberOfComponents = 1 );

   template< typename Array >
   void
   writeDataArray( const Array& array, const std::string& name,  const int numberOfElements, const int writeFromElement, int numberOfComponents = 1 );

protected:
   void
   writePoints( const ParticleSystem& particles );

   void
   writePointsTemp( const ParticleSystem& particles );

   void
   writeHeader();

   std::ostream str;

   VTK::FileFormat format;

   // number of cells (in the VTK sense) written to the file
   std::uint64_t cellsCount = 0;

   // number of points written to the file
   std::uint64_t pointsCount = 0;

   // indicator if the header has been written
   bool headerWritten = false;

   // number of data arrays written in each section
   int cellDataArrays = 0;
   int pointDataArrays = 0;

   // indicator of the current section
   VTK::DataType currentSection = VTK::DataType::CellData;

   template< typename T >
   void
   writeValue( VTK::FileFormat format, std::ostream& str, T value );

};

}  // namespace Writers
}  // namespace Meshes
}  // namespace TNL

#include "VTKWriter.hpp"

