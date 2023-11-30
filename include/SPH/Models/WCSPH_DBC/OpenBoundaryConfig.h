#pragma once
#include "../../OpenBoundaryConfig.h"
#include "BoundaryConditionsTypes.h"
#include <SPH/SPHTraits.h>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig >
bool SPHModelInit( TNL::Config::ConfigDescription config )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;


   config.addEntry< std::string >( "identifier", "Identifier of the open boundary patch.", "empty" );
   config.addEntry< VectorType >( "orientation", "Orientation of the open boundary buffer.", 0 );
   config.addEntry< VectorType >( "position", "Position of the open boundary buffer.", 0 );
   config.addEntry< VectorType >( "bufferWidth", "Width of the open boundary buffer.", 0 );
   config.addEntry< VectorType >( "bufferHeight", "Height of the open boundary buffer.", 0 );
   config.addEntry< int >( "numberOfParticlesPerCell", "Maximum allowed number of particles per grid cell.", 0 );

   config.addEntry< std::string >( "rho_bc", "Define type of open boundary condition for density.", "undefine" );
   config.addEntryEnum( "fiexd" );
   config.addEntryEnum( "extrapolated" );
   config.addEntry< RealType >( "density", "Open boundary value for density.", 0 );
   config.addEntry< std::string >( "v_bc", "Define type of open boundary condition for velocity.", "undefine" );
   config.addEntry< RealType >( "velocity", "Open boundary value for density.", 0 );
   config.addEntry< RealType >( "extrapolationDetTreshold", "Velocity treshold required in case of extrapolation.", 0 );
   config.addEntryEnum( "fiexd" );
   config.addEntryEnum( "extrapolated" );
}

template< typename SPHConfig >
class DBCOpenBoundaryConfig : public OpenBoundaryConfig< SPHConfig >
{
   public:
   using Base = OpenBoundaryConfig< SPHConfig >;
   using RealType = typename Base::RealType;
   using VectorType = typename Base::VectorType;

   DBCOpenBoundaryConfig() = default;
   DBCOpenBoundaryConfig( DBCOpenBoundaryConfig&& config );

   WCSPH_BCTypes::OpenBoundaryConditionsType type;

   /**
    * Define value handling on the buffer. Three options are possible:
    *
    * - "fixed" sets constat value for given variable. Together with
    *   this, values corresponding with given model needs to be specified.
    *
    * - "profile" sets profile for given variable. Together with this
    *   lambda function specifing the profil along the buffer needs to be specified.
    *
    * - "extrapolated" extrapolates given variable from fluid to boundary
    *   buffer. Corresponding values are set based on initial condition.
    *   Together with this, treshold for the interpolation matrix needs to be specify:
    *   use 1e-3 first order, 1e3 zero order.
    */
   std::string rho_bc = "undefine";
   RealType density = 0.f;
   RealType densityProfile = 0.f;

   std::string v_bc = "undefine";
   VectorType velocity = 0.f;
   VectorType velocityProfile = 0.f;

   RealType extrapolationDetTreshold = 0.f;
};

} // SPH
} // ParticleSystem
} // TNL

