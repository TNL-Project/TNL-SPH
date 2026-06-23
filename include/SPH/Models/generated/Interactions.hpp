// Interactions.hpp  –  AUTO-GENERATED
#include "Interactions.h"

namespace TNL {
namespace SPH {

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
void
WCSPH_DBC< Particles, ModelConfig >::interaction(
   FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
{
   auto searchInFluid = fluid->getParticles()->getSearchToken( fluid->getParticles() );
   auto searchInBound = fluid->getParticles()->getSearchToken( boundary->getParticles() );

   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType dp = modelParams.dp;
   const RealType mass = modelParams.mass;
   const RealType delta = modelParams.delta;
   const RealType alpha = modelParams.alpha;
   const RealType dynamicViscosity = modelParams.dynamicViscosity;
   const RealType speedOfSound = modelParams.speedOfSound;
   const RealType rho0 = modelParams.rho0;
   const RealType dtInit = modelParams.dtInit;
   const RealType cfl = modelParams.cfl;
   const RealType dtMin = modelParams.dtMin;
   const RealType eps = modelParams.eps;
   const VectorType gravity = modelParams.gravity;
   typename EOS::ParamsType eosParams( modelParams );

   const auto view_points_bound = boundary->getParticles()->getPoints().getConstView();
   const auto view_rho_bound = boundary->getVariables()->rho.getConstView();
   const auto view_v_bound = boundary->getVariables()->v.getConstView();
   const auto view_points = fluid->getParticles()->getPoints().getConstView();
   const auto view_rho = fluid->getVariables()->rho.getConstView();
   auto view_Drho = fluid->getVariables()->drho.getView();
   const auto view_v = fluid->getVariables()->v.getConstView();
   auto view_A = fluid->getVariables()->a.getView();

   auto fluid_fluid = [=] __cuda_callable__ (
         LocalIndexType i, LocalIndexType j,
         VectorType& r_i,
         VectorType& v_i,
         RealType& rho_i,
         RealType& p_i,
         RealType* drho_i,
         VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];
         const VectorType v_j = view_v[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const auto eps_h2 = h * h * eps;
         const auto v_ij = v_i - v_j;
         const auto drdv = ( r_ij, v_ij );
         const auto mu = h * drdv / ( ( drs * drs + eps_h2 ) );
         const auto visco = ( drdv < 0.f ) ? ( -2.f * alpha * speedOfSound * mu / ( ( rho_i + rho_j ) ) ) : ( 0.f );
         const auto psi = 2.f * h * delta * speedOfSound * ( rho_j - rho_i ) / ( ( drs * drs + eps_h2 ) );

         *drho_i += ( v_ij, gradW ) * mass - psi * ( r_ij, gradW ) * mass / rho_j;
         *a_i += -( ( p_i + p_j ) / ( rho_i * rho_j ) + visco ) * gradW * mass;
      }
   };

   auto fluid_boundary = [=] __cuda_callable__ (
         LocalIndexType i, LocalIndexType j,
         VectorType& r_i,
         VectorType& v_i,
         RealType& rho_i,
         RealType& p_i,
         RealType* drho_i,
         VectorType* a_i ) mutable
   {
      const VectorType r_j = view_points_bound[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho_bound[ j ];
         const VectorType v_j = view_v_bound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const auto eps_h2 = h * h * eps;
         const auto v_ij = v_i - v_j;
         const auto drdv = ( r_ij, v_ij );
         const auto mu = h * drdv / ( ( drs * drs + eps_h2 ) );
         const auto visco = ( drdv < 0.f ) ? ( -2.f * alpha * speedOfSound * mu / ( ( rho_i + rho_j ) ) ) : ( 0.f );
         const auto psi = 2.f * h * delta * speedOfSound * ( rho_j - rho_i ) / ( ( drs * drs + eps_h2 ) );

         *drho_i += ( v_ij, gradW ) * mass - psi * ( r_ij, gradW ) * mass / rho_j;
         *a_i += -( ( p_i + p_j ) / ( rho_i * rho_j ) + visco ) * gradW * mass;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points[ i ];
      const VectorType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      RealType drho_i = 0.f;
      VectorType a_i = VectorType( 0.f );

      Particles::NeighborsLoop::exec( i, r_i, searchInFluid, fluid_fluid,
                        v_i, rho_i, p_i, &drho_i, &a_i );
      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInBound, fluid_boundary,
                        v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      a_i += gravity;
      view_A[ i ] = a_i;
   };
   fluid->getParticles()->forAll( particleLoop );

   if( fluid->periodicPatches.size() > 0 ) {
      for( long unsigned int p = 0; p < std::size( fluid->periodicPatches ); p++ ) {

         const auto zoneParticleIndices_view = fluid->periodicPatches[ p ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = fluid->periodicPatches[ p ]->particleZone.getNumberOfParticles();
         const VectorType shift = fluid->periodicPatches[ p ]->config.shift;

         auto periodicParticleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            const GlobalIndexType pi = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points[ pi ] + shift;
            const VectorType v_i = view_v[ pi ];
            const RealType rho_i = view_rho[ pi ];
            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

            RealType drho_i = 0.f;
            VectorType a_i = VectorType( 0.f );

            Particles::NeighborsLoop::exec( pi, r_i, searchInFluid, fluid_fluid,
                              v_i, rho_i, p_i, &drho_i, &a_i );
            Particles::NeighborsLoopAnotherSet::exec( pi, r_i, searchInBound, fluid_boundary,
                              v_i, rho_i, p_i, &drho_i, &a_i );

            view_Drho[ pi ] += drho_i;
            view_A[ pi ] += a_i;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename FluidPointer, typename BoundaryPointer >
void
WCSPH_DBC< Particles, ModelConfig >::updateSolidBoundary(
   FluidPointer& fluid, BoundaryPointer& boundary, ModelParams& modelParams )
{
   auto searchInFluid = boundary->getParticles()->getSearchToken( fluid->getParticles() );

   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType dp = modelParams.dp;
   const RealType mass = modelParams.mass;
   const RealType delta = modelParams.delta;
   const RealType alpha = modelParams.alpha;
   const RealType dynamicViscosity = modelParams.dynamicViscosity;
   const RealType speedOfSound = modelParams.speedOfSound;
   const RealType rho0 = modelParams.rho0;
   const RealType dtInit = modelParams.dtInit;
   const RealType cfl = modelParams.cfl;
   const RealType dtMin = modelParams.dtMin;
   const RealType eps = modelParams.eps;
   const VectorType gravity = modelParams.gravity;
   typename EOS::ParamsType eosParams( modelParams );

   const auto view_points_bound = boundary->getParticles()->getPoints().getConstView();
   const auto view_rho_bound = boundary->getVariables()->rho.getConstView();
   auto view_Drho_bound = boundary->getVariables()->drho.getView();
   const auto view_v_bound = boundary->getVariables()->v.getConstView();
   const auto view_points = fluid->getParticles()->getPoints().getConstView();
   const auto view_rho = fluid->getVariables()->rho.getConstView();
   const auto view_v = fluid->getVariables()->v.getConstView();

   auto boundary_fluid = [=] __cuda_callable__ (
         LocalIndexType i, LocalIndexType j,
         VectorType& r_i,
         VectorType& v_i,
         RealType& rho_i,
         RealType& p_i,
         RealType* drho_i ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];
         const VectorType v_j = view_v[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType F = KernelFunction::F( drs, h );
         const VectorType gradW = r_ij * F;

         const auto v_ij = v_i - v_j;

         *drho_i += ( v_ij, gradW ) * mass * rho_i / rho_j;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
   {
      const VectorType r_i = view_points_bound[ i ];
      const VectorType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

      RealType drho_i = 0.f;

      Particles::NeighborsLoopAnotherSet::exec( i, r_i, searchInFluid, boundary_fluid,
                        v_i, rho_i, p_i, &drho_i );

      view_Drho_bound[ i ] = drho_i;
   };
   boundary->getParticles()->forAll( particleLoop );

   if( boundary->periodicPatches.size() > 0 ) {
      for( long unsigned int p = 0; p < std::size( boundary->periodicPatches ); p++ ) {

         const auto zoneParticleIndices_view = boundary->periodicPatches[ p ]->particleZone.getParticlesInZone().getConstView();
         const GlobalIndexType numberOfZoneParticles = boundary->periodicPatches[ p ]->particleZone.getNumberOfParticles();
         const VectorType shift = boundary->periodicPatches[ p ]->config.shift;

         auto periodicParticleLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
         {
            const GlobalIndexType pi = zoneParticleIndices_view[ i ];
            const VectorType r_i = view_points_bound[ pi ] + shift;
            const VectorType v_i = view_v_bound[ pi ];
            const RealType rho_i = view_rho_bound[ pi ];
            const RealType p_i = EOS::DensityToPressure( rho_i, eosParams );

            RealType drho_i = 0.f;

            Particles::NeighborsLoopAnotherSet::exec( pi, r_i, searchInFluid, boundary_fluid,
                              v_i, rho_i, p_i, &drho_i );

            view_Drho_bound[ pi ] += drho_i;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfZoneParticles, periodicParticleLoop );
      }
   }
}

template< typename Particles, typename ModelConfig >
template< typename EOS_, typename PhysicalObjectPointer >
void
WCSPH_DBC< Particles, ModelConfig >::computePressureFromDensity(
   PhysicalObjectPointer& physicalObject, ModelParams& modelParams )
{
   auto view_rho = physicalObject->getVariables()->rho.getView();
   auto view_p = physicalObject->getVariables()->p.getView();
   typename EOS_::ParamsType eosParams( modelParams );

   auto evalPressure = [=] __cuda_callable__ ( int i ) mutable
   {
      view_p[ i ] = EOS_::DensityToPressure( view_rho[ i ], eosParams );
   };
   physicalObject->getParticles()->forAll( evalPressure );
}

} // SPH
} // TNL

// end Interactions.hpp