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


}  // namespace VTK
}  // namespace ParticleSystem
}  // namespace TNL


