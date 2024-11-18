# pragma once

#include <fstream>
#include "../../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename SPHDefs >
class WCSPHOpenBoundaryEnergyFields
{
public:
   using SPHConfig = typename SPHDefs::SPHConfig;
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;

   // kinetic energy - sum
   ScalarArrayType ekin_inlet;
   ScalarArrayType ekin_outlet;
   RealType EkinSum_inlet = 0.f;
   RealType EkinSum_outlet = 0.f;

   // potential energy - sum
   ScalarArrayType epot_inlet;
   ScalarArrayType epot_outlet;
   RealType EpotSum_inlet = 0.f;
   RealType EpotSum_outlet = 0.f;

   // compresibillity energy - sum
   ScalarArrayType ecomp_inlet;
   ScalarArrayType ecomp_outlet;
   RealType EcompSum_inlet = 0.f;
   RealType EcompSum_outlet = 0.f;

   template< typename OpenBoundaryPointer >
   void
   addInlet( OpenBoundaryPointer& openBoundary )
   {
      const IndexType patchSize = openBoundary->getParticles()->getNumberOfAllocatedParticles();
      const IndexType currentSize = ekin_inlet.getSize();
      ekin_inlet.resize( currentSize + patchSize );
      epot_inlet.setSize( currentSize + patchSize );
      ecomp_inlet.setSize( currentSize + patchSize );
      ekin_inlet = 0.;
      epot_inlet = 0.;
      ecomp_inlet = 0.;
   }

   template< typename OpenBoundaryPointer >
   void
   addOutlet( OpenBoundaryPointer& openBoundary )
   {
      // This should be allocated based on the number of particles that can be contained in the
      // buffer zone.
      const IndexType patchNumberOfCells = openBoundary->zone.getNumberOfCells();
      const IndexType particlesPerCell = openBoundary->zone.getNumberOfParticlesPerCell();
      const IndexType patchSize = patchNumberOfCells * particlesPerCell;
      const IndexType currentSize = ekin_inlet.getSize();
      ekin_outlet.setSize( currentSize + patchSize );
      epot_outlet.setSize( currentSize + patchSize );
      ecomp_outlet.setSize( currentSize +  patchSize );
      ekin_outlet = 0.;
      epot_outlet = 0.;
      ecomp_outlet = 0.;
   }

   void
   output( const std::string& outputPath, const int step, const RealType time )
   {
      EkinSum_inlet = TNL::sum( ekin_inlet );
      EpotSum_inlet = TNL::sum( epot_inlet );
      EcompSum_inlet = TNL::sum( ecomp_inlet );
      EkinSum_outlet = TNL::sum( ekin_outlet );
      EpotSum_outlet = TNL::sum( epot_outlet );
      EcompSum_outlet = TNL::sum( ecomp_outlet );

      std::ofstream outfile;
      outfile.open(outputPath, std::ios_base::app );
      outfile << step << " " << time << " " << EkinSum_inlet << " " << EpotSum_inlet << " " << EcompSum_inlet << " " << EkinSum_outlet << " " << EpotSum_outlet << " " << EcompSum_outlet << std::endl;
   }

   template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams >
   void
   computeInflowEnergyLevels( FluidPointer& fluid,
                              OpenBoundaryPointer& openBoundary,
                              ModelParams& modelParams,
                              const RealType& dt )
   {
      using EOS = typename ModelParams::EOS;
      const RealType m = modelParams.mass;
      const VectorType gravity = modelParams.gravity;
      typename EOS::ParamsType eosParams( modelParams );

      auto ekin_inlet_view = ekin_inlet.getView();
      auto epot_inlet_view = epot_inlet.getView();
      auto ecomp_inlet_view = ecomp_inlet.getView();

      const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();
      const VectorType inletOrientation = openBoundary->parameters.orientation;
      const VectorType bufferPosition = openBoundary->parameters.position;

      const auto r_view = openBoundary->getParticles()->getPoints().getConstView();
      const auto v_view = openBoundary->variables->v.getConstView();
      const auto rho_view = openBoundary->variables->rho.getConstView();

      // reset energy levels
      ekin_inlet_view = 0.f;
      epot_inlet_view = 0.f;
      ecomp_inlet_view = 0.f;

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const VectorType r_i = r_view[ i ];
         const VectorType v_i = v_view[ i ];
         const RealType rho_i = rho_view[ i ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

         const VectorType r_new = r_i + dt * v_i;
         const VectorType r_relative = r_new - bufferPosition;

         // if particle moves out of the buffer, compute its energy contribution
         if( ( r_relative, inletOrientation ) > 0 ){
            ekin_inlet_view[ i ] = 0.5 * m * ( v_i, v_i );
            epot_inlet_view[ i ] = ( -1.f ) * m * ( gravity, r_new );
            ecomp_inlet_view[ i ] = m * m * p_i / rho_i;
         }
      };
      openBoundary->particles->forAll( particleLoop );
   }

   template< typename FluidPointer, typename OpenBoundaryPointer, typename ModelParams >
   void
   computeOutletEnergyLevels( FluidPointer& fluid,
                              OpenBoundaryPointer& openBoundary,
                              ModelParams& modelParams,
                              const RealType& dt )
   {
      using EOS = typename ModelParams::EOS;
      const RealType m = modelParams.mass;
      const VectorType gravity = modelParams.gravity;
      typename EOS::ParamsType eosParams( modelParams );

      auto ekin_outlet_view = ekin_outlet.getView();
      auto epot_outlet_view = epot_outlet.getView();
      auto ecomp_outlet_view = ecomp_outlet.getView();

      const auto zoneParticleIndices_view = openBoundary->zone.getParticlesInZone().getConstView();
      const IndexType numberOfZoneParticles = openBoundary->zone.getNumberOfParticles();
      const VectorType outletOrientation = openBoundary->parameters.orientation;
      const VectorType bufferPosition = openBoundary->parameters.position;

      const auto r_view = fluid->getParticles()->getPoints().getConstView();
      const auto v_view = fluid->variables->v.getConstView();
      const auto rho_view = fluid->variables->rho.getConstView();

      // reset energy levels
      ekin_outlet_view = 0.f;
      epot_outlet_view = 0.f;
      ecomp_outlet_view = 0.f;

      auto particleLoop = [ = ] __cuda_callable__( IndexType i ) mutable
      {
         const IndexType p = zoneParticleIndices_view[ i ];
         const VectorType r_i = r_view[ p ];
         const VectorType v_i = v_view[ p ];
         const RealType rho_i = rho_view[ p ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

         const VectorType r_relative = bufferPosition - r_i;

         // if particle entered the outlet buffer, compute its energy
         if( ( r_relative, outletOrientation ) > 0 ){
            ekin_outlet_view[ i ] = 0.5 * m * ( v_i, v_i );
            epot_outlet_view[ i ] = ( -1.f ) * m * ( gravity, r_i );
            ecomp_outlet_view[ i ] = m * m * p_i / rho_i;
         }
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, particleLoop );
   }
};

} // SPH
} // TNL

