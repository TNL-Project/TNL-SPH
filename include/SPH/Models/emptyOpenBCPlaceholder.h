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

template< typename SPHState >
class EmptyVariables
{
public:

   using SPHConfig = typename SPHState::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;

   void
   setSize( const GlobalIndexType& size ) {};

   template< typename IndexArrayTypePointer >
   void
   sortVariables( IndexArrayTypePointer& map, const GlobalIndexType& numberOfParticles ) {};

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader ) {};

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType& startingIdx = 0 ) {};

};

} //namespace SPH
} //namespace TNL

