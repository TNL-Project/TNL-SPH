#include "Interactions.h"
#include "../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename NeighborSearchPointer, typename ModelPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::Interaction( NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound, ModelPointer& boundary )
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = this->particles->getSearchRadius();

   static constexpr RealType gridXbegin = Particles::Config::gridXbegin; //FIXIT
   static constexpr RealType gridYbegin = Particles::Config::gridYbegin; //FIXIT

   const auto view_firstLastCellParticle = neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   const auto view_points = particles->getPoints().getView();

   const auto view_firstLastCellParticle_bound = neighborSearch_bound->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex_bound = boundary->particles->getParticleCellIndices().getView();
   const auto view_points_bound = boundary->particles->getPoints().getView();

   /* CONSTANT VARIABLES */
   const RealType h = this->h;
   const RealType m = this->m;
   const RealType speedOfSound = this->speedOfSound;
   const RealType coefB = this->coefB;
   const RealType rho0 = this->rho0;
   const RealType delta = this->delta;
   const RealType alpha = this->alpha;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_rho = this->FluidVariables.rho.getView();
   auto view_Drho = this->FluidVariables.drho.getView();
   const auto view_v = this->FluidVariables.v.getView();
   auto view_a = this->FluidVariables.a.getView();

   const auto view_rho_bound = boundary->getFluidVariables().rho.getView();
   auto view_Drho_bound = boundary->getFluidVariables().drho.getView();
   const auto view_v_bound = boundary->getFluidVariables().v.getView();

   auto FluidFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, PointType& r_i, PointType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, PointType* a_i ) mutable
   {
      /* This should be some interaction structure - properties of particle B:*/
      const PointType r_j = view_points[ j ];
      const PointType dr = r_i - r_j;
      const RealType drs = l2Norm( dr );
      if (drs <= searchRadius )
      {
         const PointType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Interaction: */
         const PointType dv = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const PointType gradW = dr * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
         *drho_i += ( dv, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( dr, dv ) );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;
      }
   };

   auto FluidBound = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, PointType& r_i, PointType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, PointType* a_i ) mutable
   {
      /* This should be some interaction structure - properties of particle B:*/
      const PointType r_j = view_points_bound[ j ];
      const PointType dr = r_i - r_j;
      const RealType drs = l2Norm( dr );
      if (drs <= searchRadius )
      {
         const PointType v_j = view_v_bound[ j ];
         const RealType rho_j = view_rho_bound[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );

         /* Interaction: */
         const PointType dv = v_i - v_j;

         const RealType F = SPHKernelFunction::F( drs, h );
         const PointType gradW = dr * F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
         *drho_i += ( dv, gradW ) * m - diffTerm;

         const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
         const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( dr, dv ) );
         *a_i += ( -1.0f ) * ( p_term + visco )* gradW * m;
      }
   };

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, PointType& r_i, PointType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      //if( view_particleType[ j ] == 0 ){
         /* This should be some interaction structure, mby. - properties of particle B: */
         const PointType r_j = view_points[ j ];
         const PointType dr = r_i - r_j;
         const RealType drs = l2Norm( dr );
         if( drs <= searchRadius )
         {
            const PointType v_j = view_v[ j ];
            const RealType rho_j = view_rho[ j ];
            const RealType p_j = EOS::DensityToPressure( rho_j );

            /* Interaction */
            const PointType dv = v_i - v_j;

            const RealType F = SPHKernelFunction::F( drs, h );
            const PointType gradW = dr*F;

            const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
            const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
            *drho_i += ( dv, gradW ) * m - diffTerm;
         }
      //}
   };



   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const PointType r_i = view_points[ i ];
      const PointType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i );

      const int gridIndexI = TNL::floor( ( r_i[ 0 ] - gridXbegin ) / searchRadius );
      const int gridIndexJ = TNL::floor( ( r_i[ 1 ] - gridYbegin ) / searchRadius );

      PointType a_i = {0.f, 0.f};
      RealType drho_i = 0.f;

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndexI, gridIndexJ, view_firstLastCellParticle, view_particleCellIndex, FluidFluid, r_i, v_i, rho_i, p_i, &drho_i, &a_i );
      neighborSearch_bound->loopOverNeighbors( i, numberOfParticles_bound, gridIndexI, gridIndexJ, view_firstLastCellParticle_bound, view_particleCellIndex, FluidBound, r_i, v_i, rho_i, p_i, &drho_i, &a_i );

      view_Drho[ i ] = drho_i;
      a_i[ 1 ] -= 9.81f ;
      view_a[ i ] = a_i;
   };
   SPHParallelFor::exec( 0, numberOfParticles, particleLoop, neighborSearch, neighborSearch_bound );

   auto particleLoopBoundary = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const PointType r_i = view_points_bound[ i ];
      const PointType v_i = view_v_bound[ i ];
      const RealType rho_i = view_rho_bound[ i ];
      const RealType p_i = EOS::DensityToPressure( rho_i );

      const int gridIndexI = TNL::floor( ( r_i[ 0 ] - gridXbegin ) / searchRadius );
      const int gridIndexJ = TNL::floor( ( r_i[ 1 ] - gridYbegin ) / searchRadius );

      PointType a_i = {0.f, 0.f};
      RealType drho_i = 0.;

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndexI, gridIndexJ, view_firstLastCellParticle, view_particleCellIndex_bound, BoundFluid, r_i, v_i, rho_i, p_i, &drho_i );

      view_Drho_bound[ i ] = drho_i;

   };
   SPHParallelFor::exec( 0, numberOfParticles_bound, particleLoopBoundary, neighborSearch );
}

} // SPH
} // ParticleSystem
} // TNL

