#include "Fluid.h"
#include "Boundary.h"

namespace TNL {
namespace SPH {

template < typename Model >
class SPHMultiset_CFDDomain
{
   using DeviceType = typename Model::DeviceType;
   using ModelType = Model;
   using ModelParams = typename ModelType::ModelParams;
   using ParticlesType = typename ModelType::ParticlesType;;
   using IntegrationSchemeType = typename ModelType::IntegrationSchemeType;
   using IntegratorPointer = typename Pointers::SharedPointer< IntegrationSchemeType, DeviceType >;
   using IntegrationSchemeVariablesType = typename Model::IntegrationSchemeVariables;
   using TimeStepping = typename Model::ModelParams::TimeStepping;

   using FluidVariables = typename Model::FluidVariables;
   using Fluid = Fluid< ParticlesType, SPHConfig, FluidVariables, IntegrationSchemeVariablesType >;
   using FluidPointer = Pointers::SharedPointer< Fluid, DeviceType >;
   using BoundaryVariables = typename Model::BoundaryVariables;
   using Boundary = Boundary< ParticlesType, SPHConfig, BoundaryVariables, IntegrationSchemeVariablesType >;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;

   SPHMultiset_CFDDomain() = default;

   /**
    * \brief Set parameters of the simulation domain.
    *
    * \param parameters global parameters config with data about simulation.
    * \param logger global logger to write log from initialization.
    */
   void
   init( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   /**
    * \brief Perform neighbors search for all particle objects in the domain.
    */
   void
   performNeighborSearch( TNL::Logger& log );

   /**
    * \brief Perform interaction between all particles and all particle objects
    * in the simulation.
    */
   void
   interact();

protected:

   FluidPointer fluid;
   FluidPointer fluidOverlap;

   BoundaryPointer boundary;
   BoundaryPointer boundaryOverlap;

   std::vector< OpenBoundaryPointer > openBoundaryPatches;
};

} // SPH
} // TNL

