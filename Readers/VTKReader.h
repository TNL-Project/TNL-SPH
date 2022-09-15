#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <map>

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

   }


};

} // Readers
} // ParticleSystem
} // TNL
