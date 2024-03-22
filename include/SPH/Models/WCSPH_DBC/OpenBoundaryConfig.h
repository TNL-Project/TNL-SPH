#pragma once
#include "../../OpenBoundaryConfig.h"
#include "BoundaryConditionsTypes.h"
#include "TNL/Containers/Expressions/ExpressionTemplates.h"
#include "TNL/Containers/Vector.h"
#include "TNL/Logger.h"
#include <SPH/SPHTraits.h>
#include <iterator>

namespace TNL {
namespace SPH {

template< typename SPHConfig >
void
configSetupOpenBoundaryModelPatch( TNL::Config::ConfigDescription& config, std::string prefix )
{
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   // setup base open boundary confit
   conigSetupOpenBoundaryPatch< SPHConfig >( config, prefix );

   // parameters for model-dependent open boundary conditions
   config.addEntry< std::string >( prefix + "type", "Type of open boundary patch.", "none" );
   config.addEntry< std::string >( prefix + "rho_bc", "Define type of open boundary condition for density.", "undefine" );
      config.addEntryEnum( "fixed" );
      config.addEntryEnum( "extrapolated" );
      config.addEntryEnum( "do-nothing" );
   config.addEntry< RealType >( prefix + "density", "Open boundary value for density.", 0 );
   config.addEntry< std::string >( prefix + "v_bc", "Define type of open boundary condition for velocity.", "undefine" );
      config.addEntryEnum( "fixed" );
      config.addEntryEnum( "extrapolated" );
      config.addEntryEnum( "do-nothing" );
   config.addEntry< RealType >( prefix + "velocity-x", "Open boundary value for density.", 0 );
   config.addEntry< RealType >( prefix + "velocity-y", "Open boundary value for density.", 0 );
   config.addEntry< RealType >( prefix + "velocity-z", "Open boundary value for density.", 0 );
   config.addEntry< RealType >( prefix + "extrapolationDetTreshold", "Velocity treshold required in case of extrapolation.", 0 );
}

template< typename SPHConfig >
class DBCOpenBoundaryConfig : public OpenBoundaryConfig< SPHConfig >
{
   public:
   using Base = OpenBoundaryConfig< SPHConfig >;
   using RealType = typename Base::RealType;
   using VectorType = typename Base::VectorType;

   DBCOpenBoundaryConfig() = default;

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

   void
   init( TNL::Config::ParameterContainer& parameters, std::string prefix )
   {
      Base::init( parameters, prefix );

      //identify buffer type
      const std::string openBoundaryPatchType = parameters.getParameter< std::string >( prefix + "type" );
      if( openBoundaryPatchType == "periodic" )
         type = WCSPH_BCTypes::OpenBoundaryConditionsType::PeriodicBoundary;
      else if( openBoundaryPatchType == "inlet" )
         type = WCSPH_BCTypes::OpenBoundaryConditionsType::Inlet;
      else if( openBoundaryPatchType == "outlet" )
         type = WCSPH_BCTypes::OpenBoundaryConditionsType::Outlet;
      else
         std::cerr << "Open boundary patch " << prefix << " has invalid type " << openBoundaryPatchType << "." << std::endl;

      //: //identify density boundary condition
      //: const std::string densityBoundaryCondition = parameters.getParameter< std::string >( prefix + "rho_bc" );
      //: if( densityBoundaryCondition == "fixed" )
      //:    rho_bc = densityBoundaryCondition;
      //: else if( densityBoundaryCondition == "extrapolated" ){
      //:    rho_bc = densityBoundaryCondition;
      //: }
      //: else
      //:    std::cerr << "Open boundary patch " << prefix << " has invalid type " << openBoundaryPatchType << "." << std::endl;

      this->rho_bc = parameters.getParameter< std::string >( prefix + "rho_bc" );
      density = parameters.getParameter< RealType >( prefix + "density" );
      this->v_bc = parameters.getParameter< std::string >( prefix + "v_bc" );
      velocity = parameters.getXyz< VectorType >( prefix + "velocity" );
      this->extrapolationDetTreshold = parameters.getParameter< RealType >( prefix + "extrapolationDetTreshold" );

   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      Base::writeProlog( logger );
      logger.writeSeparator();
      std::string openBoundaryPatchType;
      if( type == WCSPH_BCTypes::OpenBoundaryConditionsType::PeriodicBoundary ){
         openBoundaryPatchType = "periodic boundary";
         logger.writeParameter( "Open boundary patch type: ", openBoundaryPatchType );
         logger.writeParameter( "Paired periodic boundary patch index: ", this->pairedPeriodicBuffer );
         logger.writeParameter( "Periodic boundary shift vector: ", this->shift );
         return;
      }
      if( type == WCSPH_BCTypes::OpenBoundaryConditionsType::Inlet )
         openBoundaryPatchType = "inlet";
      if( type == WCSPH_BCTypes::OpenBoundaryConditionsType::Outlet )
         openBoundaryPatchType = "outlet";
      logger.writeParameter( "Open boundary patch type: ", openBoundaryPatchType );
      logger.writeParameter( "Density boundary condition:", rho_bc );
      logger.writeParameter( "Prescribed boundary density:", density );
      logger.writeParameter( "Velocity boundary conditon:", v_bc );
      logger.writeParameter( "Prescribed boundary velocity:", velocity );

   }
};

} // SPH
} // TNL

