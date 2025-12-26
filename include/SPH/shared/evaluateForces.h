# pragma once

#include <fstream>
#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename SPHDefs >
class EvaluateForces
{

public:

   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   using KernelFunction = typename SPHDefs::KernelFunction;
   using EOS = typename SPHDefs::EOS;
   using ViscousTerm = typename SPHDefs::ViscousTerm;

   EvaluateForces() = default;

   template< typename BoudaryPointer >
   EvaluateForces( BoudaryPointer& boundary )
   {
      init( boundary );
   }

   template< typename FluidPointer >
   void
   init( FluidPointer& fluid )
   {
      const IndexType size = fluid->getParticles()->getNumberOfAllocatedParticles();
      F_pressure.setSize( size );
      F_visco.setSize( size );
      F_pressure = 0.f;
      F_visco = 0.f;
   }

   template< typename FluidPointer, typename BoundaryPointer, typename ModelParams >
   void
   computeForces( FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams, const int markerValueToProcess = 0 )
   {
      auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->getParticles() );
      const RealType searchRadius = boundary->getParticles()->getSearchRadius();
      const RealType h = modelParams.h;
      const RealType m = modelParams.mass;
      typename ModelParams::ViscousTerm::ParamsType viscousTermTermsParams( modelParams );
      typename ModelParams::EOS::ParamsType eosParams( modelParams );

      const auto view_points = fluid->getParticles()->getPoints().getConstView();
      const auto view_rho = fluid->getVariables()->rho.getConstView();
      const auto view_v = fluid->getVariables()->v.getConstView();
      const auto view_points_bound = boundary->getParticles()->getPoints().getConstView();
      const auto view_rho_bound = boundary->getVariables()->rho.getConstView();
      const auto view_v_bound = boundary->getVariables()->v.getConstView();
      const auto view_marker_bound = boundary->getVariables()->marker.getConstView();

      auto view_F_pressure = F_pressure.getView();
      auto view_F_visco = F_visco.getView();

      auto comptueForce = [=] __cuda_callable__ (
            IndexType i,
            IndexType j,
            VectorType& r_i,
            VectorType& v_i,
            RealType& rho_i,
            RealType& p_i,
            VectorType* F_pressure_i,
            VectorType* F_visco_i ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
            const VectorType v_j = view_v[ j ];
            const VectorType v_ij = v_i - v_j;


            const RealType WV_j = KernelFunction::W( drs, h ) * m / rho_j;
            const VectorType gradWV_j = r_ij * KernelFunction::F( drs, h ) * m / rho_j;

            const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );

            *F_pressure_i += ( p_i + p_j ) * gradWV_j;
            *F_visco_i += ViscousTerm::Pi( rho_i, rho_j, drs, ( r_ij, v_ij ), viscousTermTermsParams ) * ( gradWV_j * rho_j );
         }
      };

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_i = view_points_bound[ i ];
         const VectorType v_i = view_v_bound[ i ];
         const RealType rho_i = view_rho_bound[ i ];
         const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );
         const int maker_i  = view_marker_bound[ i ];
         VectorType F_pressure_i = 0.f;
         VectorType F_visco_i = 0.f;

         ParticlesType::NeighborsLoop::exec( i, r_i, searchInFluid, comptueForce, v_i, rho_i, p_i, &F_pressure_i, &F_visco_i );

         F_pressure_i *= m / rho_i;
         F_visco_i *= m / rho_i;

         if( maker_i == markerValueToProcess ){
            view_F_pressure[ i ] = F_pressure_i;
            view_F_visco[ i ] = F_visco_i;
         }
         else{
            view_F_pressure[ i ] = 0.f;
            view_F_visco[ i ] = 0.f;
         }
      };
      boundary->getParticles()->forAll( particleLoop );

      F_pressure_sum = TNL::l2Norm( TNL::sum( F_pressure ) );
      F_visco_sum = TNL::l2Norm( TNL::sum( F_visco ) );
   }

   void
   output( const std::string& outputPath, const int step, const RealType time )
   {
      std::ofstream outfile;
      outfile.open(outputPath, std::ios_base::app );
      outfile << step << " " << time << " " << F_pressure_sum  << " " << F_visco_sum << std::endl;
   }

   VectorArrayType F_pressure;
   VectorArrayType F_visco;

   RealType F_pressure_sum;
   RealType F_visco_sum;

};

} // SPH
} // TNL

