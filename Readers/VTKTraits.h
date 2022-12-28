#pragma once

#include <stdexcept>
#include <string>
#include <cstdint>

namespace TNL {
namespace ParticleSystem {
namespace VTK {

// VTK file formats
enum class FileFormat : std::uint8_t
{
   ascii,
   binary,
   zlib_compressed
};

// VTK data types
enum class DataType : std::uint8_t
{
   CellData,
   PointData
};


}  // namespace VTK
}  // namespace ParticleSystem
}  // namespace TNL

