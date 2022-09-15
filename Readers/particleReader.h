#pragma once

#include <string>
#include <utility>
#include <vector>
//#include <TNL/3rdparty/mpark/variant.hpp>  // backport of std::variant from C++17
//#include "variant.hpp"
#include <mpark/variant.hpp>

namespace TNL {
namespace ParticleSystem {
namespace Readers {

class ParticleReader
{
public:

   using VariantVector = mpark::variant< std::vector< std::int8_t >,
                                         std::vector< std::uint8_t >,
                                         std::vector< std::int16_t >,
                                         std::vector< std::uint16_t >,
                                         std::vector< std::int32_t >,
                                         std::vector< std::uint32_t >,
                                         std::vector< std::int64_t >,
                                         std::vector< std::uint64_t >,
                                         std::vector< float >,
                                         std::vector< double > >;

   ParticleReader() = default;

   ParticleReader( std::string fileName ) : fileName( std::move( fileName ) ) {}

   virtual ~ParticleReader() = default;

   void
   setFileName( const std::string& fileName )
   {
      reset();
      this->fileName = fileName;
   }

   /**
    * \brief This method resets the reader to an empty state.
    *
    * In particular, implementations should call the \ref resetBase method to
    * reset the arrays holding the intermediate particle representation.
    */
   virtual void
   reset()
   {
      resetBase();
   }

   virtual void
   detectParticleSystem() = 0;

   template< typename ParticleType >
   void // ... -> std::enable_if_t< isGrid< MeshType >::value >
   loadParticle( ParticleType& particles )
   {

      // skip empty particleSystem (the cell shape is indeterminate)
      if( NumberOfPoints == 0 ) {
         particles = ParticleType{};
         return;
      }

      using PointType = typename ParticleType::PointType;

      // assign points
      using mpark::visit;
      visit(
         [ &particles ]( auto&& array )
         {
            PointType p;
            std::size_t i = 0;
            for( auto c : array ) {
               int dim = i++ % 3;
               if( dim >= PointType::getSize() )
                  continue;
               p[ dim ] = c;
               if( dim == PointType::getSize() - 1 )
                  particles.setPoint( ( i - 1 ) / 3, p );
            }
         },
         pointsArray );



   }

   int
   getSpaceDimension() const
   {
      return spaceDimension;
   }


protected:
   // input file name
   std::string fileName;

   // attributes of the particle system
   std::size_t NumberOfPoints;
   int spaceDimension;

   // string representation of ptcs types (forced means specified by the user, otherwise
   // the type detected by detectMesh takes precedence)
   std::string forcedRealType;
   std::string forcedGlobalIndexType;
   std::string forcedLocalIndexType = "short int";  // not stored in any file format

   // points
   std::vector< std::int32_t > pointsArray; // ... -> VariantVector pointsArray;

   // string representation of each array's value type
   std::string pointsType, connectivityType, offsetsType, typesType;

   void
   resetBase()
   {
      NumberOfPoints = 0;
      spaceDimension = 0;

      pointsArray = {};
      pointsType = connectivityType = offsetsType = typesType = "";
   }






};

} // Readers
} // ParticleSystem
} // TNL
