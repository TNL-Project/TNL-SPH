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

   // kinetic energy - derivatives
   ScalarArrayType dekindt;
   RealType dEkindt = 0.f;
   RealType Ekin = 0.f;

   // potential energy - derivatives
   ScalarArrayType depotdt;
   RealType dEpotdt = 0.f;
   RealType Epot = 0.f;

   // compresibillity energy - derivatives
   ScalarArrayType decompdt;
   RealType dEcompdt = 0.f;
   RealType Ecomp = 0.f;

   // kinetic energy - sum
   ScalarArrayType ekin;
   RealType EkinSum = 0.f;

   // potential energy - sum
   ScalarArrayType epot;
   RealType EpotSum = 0.f;

   // compresibillity energy - sum
   ScalarArrayType ecomp;
   RealType EcompSum = 0.f;

   bool energySnapshots = false;

   /**
   bool includeOpenBoundary = false;

   // kinet energy related to open boudaries
   ScalarArrayType dekindt_inlet;
   ScalarArrayType dekindt_outlet;
   RealType dEkindt_inlet = 0.f;
   RealType dEkindt_outlet = 0.f;

   // potential energy related to open boudaries
   ScalarArrayType depotdt_inlet;
   ScalarArrayType depotdt_outlet;
   RealType dEpotdt_inlet = 0.f;
   RealType dEpotdt_outlet = 0.f;

   // compresibillity energy related to open boudaries
   ScalarArrayType decompdt_inlet;
   ScalarArrayType decompdt_outlet;
   RealType dEcompdt_inlet = 0.f;
   RealType dEcompdt_outlet = 0.f;
   */

   WCSPHEnergyFields() = default;

   template< typename FluidPointer >
   WCSPHEnergyFields( FluidPointer& fluid, bool initEnergySnapshots = false )
   {
      init( fluid, initEnergySnapshots );
   }

   template< typename FluidPointer >
   void
   init( FluidPointer& fluid, bool initEnergySnapshots = false )
   {
      const IndexType size = fluid->getParticles()->getNumberOfAllocatedParticles();

      dekindt.setSize( size );
      depotdt.setSize( size );
      decompdt.setSize( size );
      dekindt = 0.;
      depotdt = 0.;
      decompdt = 0.;

      energySnapshots = initEnergySnapshots;
      if( energySnapshots )
      {
         ekin.setSize( size );
         epot.setSize( size );
         ecomp.setSize( size );
         ekin = 0.;
         epot = 0.;
         ecomp = 0.;
      }
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
      if( energySnapshots )
         outfile << step << " " << time << " " << Ekin << " " << Epot << " " << Ecomp << " " << EkinSum << " " << EpotSum << " " << EcompSum << std::endl;
      else
         outfile << step << " " << time << " " << Ekin << " " << Epot << " " << Ecomp << std::endl;
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

      const auto v_view = fluid->getVariables()->v.getConstView();
      const auto dvdt_view = fluid->getVariables()->a.getConstView();
      const auto rho_view = fluid->getVariables()->rho.getConstView();
      const auto drhodt_view = fluid->getVariables()->drho.getConstView();

      // reset energy derivatives
      dekindt_view = 0.f;
      depotdt_view = 0.f;
      decompdt_view = 0.f;

      // compute current rates
      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const RealType rho_i = rho_view[ i ];
         const RealType drhodt_i = drhodt_view[ i ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
         const VectorType v_i = v_view[ i ];
         const VectorType dvdt_i = dvdt_view[ i ];

         dekindt_view[ i ] = m * ( v_i, dvdt_i );
         depotdt_view[ i ] = ( -1.f ) * m * ( gravity, v_i );
         decompdt_view[ i ] = m * p_i / ( rho_i * rho_i ) * drhodt_i;
      };
      fluid->getParticles()->forAll( particleLoop );

      dEkindt = TNL::sum( dekindt );
      dEpotdt = TNL::sum( depotdt );
      dEcompdt = TNL::sum( decompdt );
   }

   template< typename FluidPointer, typename ModelParams  >
   void
   computeEnergyLevels( FluidPointer& fluid, ModelParams& modelParams )
   {
      using EOS = typename ModelParams::EOS;

      const RealType m = modelParams.mass;
      const VectorType gravity = modelParams.gravity;
      typename EOS::ParamsType eosParams( modelParams );

      auto ekin_view = ekin.getView();
      auto epot_view = epot.getView();
      auto ecomp_view = ecomp.getView();

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto v_view = fluid->getVariables()->v.getConstView();
      const auto rho_view = fluid->getVariables()->rho.getConstView();

      // reset energy derivatives
      ekin_view = 0.f;
      epot_view = 0.f;
      ecomp_view = 0.f;

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const VectorType r_i = r_view[ i ];
         const RealType rho_i = rho_view[ i ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
         const VectorType v_i = v_view[ i ];

         ekin_view[ i ] = 0.5 * m * ( v_i, v_i );
         epot_view[ i ] = ( -1.f ) * m * ( gravity, r_i );
         ecomp_view[ i ] = m * m * p_i / rho_i;
      };
      fluid->getParticles()->forAll( particleLoop );

      EkinSum = TNL::sum( ekin );
      EpotSum = TNL::sum( epot );
      EcompSum = TNL::sum( ecomp );
   }

   /*
   template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams >
   void
   computeEnergyDerivativesOpenBoundary( FluidPointer& fluid, OpenBoundaryPointer& openBoundary, ModelParams& modelParams )
   {
      auto dekindt_view = dekindt.getView();
      auto depotdt_view = depotdt.getView();
      auto decompdt_view = decompdt.getView();

      auto dekindt_inlet_view = dekindt_inlet.getView();
      auto depotdt_inlet_view = depotdt_inlet.getView();
      auto decompdt_inlet_view = decompdt_inlet.getView();

      const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();
      const VectorType inletOrientation = openBoundary->parameters.orientation;
      const VectorType bufferPosition = openBoundary->parameters.position;

      const auto r_view = fluid->getParticles().getPoints().getConstView();

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r = r_view[ p ];
         const VectorType r_relative = r - bufferPosition;
         if( (r_relative, inletOrientation ) > 0 )

         dekindt_inlet_view[ p ] = 0.f;
         depotdt_inlet_view[ p ] = -dekindt_view[ p ];
         decompdt_inlet_view[ p ] = decompdt_view[ p ];
         dekindt_view[ p ] = -0.f;
         depotdt_view[ p ] = 0.f;
         decomdt[ p ] = 0.f;
      };
      fluid->getParticles()->forAll( particleLoop );
   }
   */

};

} // SPH
} // TNL

