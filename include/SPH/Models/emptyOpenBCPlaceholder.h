#include "../OpenBoundaryConfig.h"

namespace TNL {
namespace SPH {

template< typename SPHConfig >
void
configSetupOpenBoundaryModelPatch( TNL::Config::ConfigDescription& config, std::string prefix )
{}

template< typename SPHConfig >
class NoOpenBC : public OpenBoundaryConfig< SPHConfig >
{
   public:
   using Base = OpenBoundaryConfig< SPHConfig >;

   NoOpenBC() = default;

   void
   init( TNL::Config::ParameterContainer& parameters,
         TNL::Config::ParameterContainer& parametersOpenBoundary,
         std::string prefix )
   {}

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {}
};

} //namespace SPH
} //namespace TNL

