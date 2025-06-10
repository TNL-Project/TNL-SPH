#pragma once

#include <TNL/Containers/StaticVector.h>
#include <TNL/Matrices/StaticMatrix.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/Vector.h>

namespace TNL {
namespace SPH {

template< typename SPHFluidConfig >
class SPHFluidTraits
{
   public:
   static constexpr int spaceDimension = SPHFluidConfig::spaceDimension;

   using DeviceType = typename SPHFluidConfig::DeviceType;
   using GlobalIndexType = typename SPHFluidConfig::GlobalIndexType;
   using LocalIndexType = typename SPHFluidConfig::LocalIndexType;
   using CellIndexType = typename SPHFluidConfig::CellIndexType;
   using RealType = typename SPHFluidConfig::RealType;

   /* particle related */
   using MarkerType = int; //FIXME
   using MarkerArrayType = Containers::Vector< MarkerType, DeviceType, GlobalIndexType >;
   using ScalarType = RealType;
   using ScalarArrayType = Containers::Vector< ScalarType, DeviceType, GlobalIndexType >;
   using VectorType = Containers::StaticVector< spaceDimension, RealType >;
   using VectorArrayType = Containers::Vector< VectorType, DeviceType, GlobalIndexType >;
   using IndexArrayType = Containers::Array< GlobalIndexType, DeviceType >;
   using IndexVectorType = Containers::StaticVector< spaceDimension, GlobalIndexType >;

   //types for correction matrices related to MDBC
   using MatrixType = Matrices::StaticMatrix< RealType, SPHFluidConfig::spaceDimension, SPHFluidConfig::spaceDimension >;
   using MatrixArrayType = Containers::Vector< MatrixType, DeviceType, GlobalIndexType >;

   using VectorExtendedType = Containers::StaticVector< SPHFluidConfig::spaceDimension + 1, RealType >;
   using VectorExtendedArrayType = Containers::Array< VectorExtendedType, DeviceType, GlobalIndexType >;
   using MatrixExtendedType = Matrices::StaticMatrix< RealType, SPHFluidConfig::spaceDimension + 1, SPHFluidConfig::spaceDimension + 1 >;
   using MatrixExtendedArrayType = Containers::Array< MatrixExtendedType, DeviceType, GlobalIndexType >;

};

} // SPH
} // TNL

