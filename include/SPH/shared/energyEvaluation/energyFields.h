# pragma once

#include <fstream>
#include "../../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename SPHDefs >
class WCSPHEnergyFields
{
public:
   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;

   // kinetic energy
   ScalarArrayType dekindt;
   RealType dEkindt = 0.f;
   RealType Ekin = 0.f;

   // potential energy
   ScalarArrayType depotdt;
   RealType dEpotdt = 0.f;
   RealType Epot = 0.f;

   // potential energy
   ScalarArrayType decompdt;
   RealType dEcompdt = 0.f;
   RealType Ecomp = 0.f;

   template< typename FluidPointer >
   void
   init( FluidPointer& fluid )
   {
      const IndexType size = fluid->getParticles()->getNumberOfAllocatedParticles();

      dekindt.setSize( size );
      depotdt.setSize( size );
      decompdt.setSize( size );

      dekindt = 0.;
      depotdt = 0.;
      decompdt = 0.;
   }

   void
   integrate( const RealType& dt )
   {
      Ekin += dt * dEkindt;
      Epot += dt * dEpotdt;
      Ecomp += dt * dEcompdt;
   }

   void
   output( const std::string& outputPath, const int step, const RealType time )
   {
      std::ofstream outfile;
      outfile.open(outputPath, std::ios_base::app );
      outfile << step << " " << time << " " << this->Ekin << " " << Epot << " " << Ecomp << std::endl;
   }

   template< typename FluidPointer, typename ModelParams  >
   void
   computeEnergyDerivatives( FluidPointer& fluid, ModelParams& modelParams )
   {
      using EOS = typename ModelParams::EOS;

      const RealType m = modelParams.mass;
      const VectorType gravity = modelParams.gravity;
      typename EOS::ParamsType eosParams( modelParams );

      auto dekindt_view = dekindt.getView();
      auto depotdt_view = depotdt.getView();
      auto decompdt_view = decompdt.getView();

      const auto v_view = fluid->variables->v.getConstView();
      const auto dvdt_view = fluid->variables->a.getConstView();
      const auto rho_view = fluid->variables->rho.getConstView();
      const auto drhodt_view = fluid->variables->drho.getConstView();

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const RealType rho_i = rho_view[ i ];
         const RealType drhodt_i = drhodt_view[ i ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
         const VectorType v_i = v_view[ i ];
         const VectorType dvdt_i = dvdt_view[ i ];

         dekindt_view[ i ] += ( -1.f ) * m * ( gravity, v_i );
         depotdt_view[ i ] += m * ( v_i, dvdt_i );
         decompdt_view[ i ] += m * p_i / ( rho_i * rho_i ) * drhodt_i;
      };
      fluid->particles->forAll( particleLoop );

      dEkindt = TNL::sum( dekindt );
      dEpotdt = TNL::sum( depotdt );
      dEcompdt = TNL::sum( decompdt );
   }
};

} // SPH
} // TNL

