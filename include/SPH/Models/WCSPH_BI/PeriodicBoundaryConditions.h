#include "../../shared/PeriodicBoundaryConditions.h"

#include "../../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Model >
class PeriodicBoundaryConditions : public PeriodicBoundaryConditionsShared< Model >
{
public:

   template< typename FluidPointer, typename BoudaryPointer, typename ParticleParams >
   void applyPeriodicBoundaryConditionForFields( FluidPointer& fluid, BoudaryPointer& boundary, ParticleParams& particleParams )
   {
      this->initializePeriodicBoundaryTransfer( fluid, particleParams );
      this->periodicUpdateOfParticleField( fluid->getFluidVariables()->rho );
      this->periodicUpdateOfParticleField( fluid->getFluidVariables()->v );
      this->finalizePeriodicBoundaryTransfer();

      this->initializePeriodicBoundaryTransfer( boundary, particleParams );
      this->periodicUpdateOfParticleField( boundary->getBoundaryVariables()->rho );
      this->periodicUpdateOfParticleField( boundary->getBoundaryVariables()->v );
      this->periodicUpdateOfParticleField( boundary->getBoundaryVariables()->n );
      this->finalizePeriodicBoundaryTransfer();
   }

   //TODO: Remove this.
   template< typename FluidPointer, typename BoudaryPointer, typename ParticleParams >
   void initializePeriodicBoundaryConditionForField( FluidPointer& fluid, BoudaryPointer& boundary, ParticleParams& particleParams )
   {
      this->initialize( fluid, particleParams );
      this->initializeParticleField( fluid, fluid->getFluidVariables()->rho, fluid->getFluidVariables()->rho_swap );
      this->initializeParticleField( fluid, fluid->getFluidVariables()->v, fluid->getFluidVariables()->v_swap  );

      this->initialize( boundary, particleParams );
      this->initializeParticleField( boundary, boundary->getBoundaryVariables()->rho, boundary->getBoundaryVariables()->rho_swap );
      this->initializeParticleField( boundary, boundary->getBoundaryVariables()->v, boundary->getBoundaryVariables()->v_swap );
      this->initializeParticleField( boundary, boundary->getBoundaryVariables()->n, boundary->getBoundaryVariables()->n_swap );
   }
};

} // SPH
} // ParticleSystem
} // TNL

