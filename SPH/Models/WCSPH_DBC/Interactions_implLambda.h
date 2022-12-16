#include "Interactions.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename Particles, typename SPHFluidConfig, typename Variables >
template< typename NeighborSearchPointer, typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm >
void
WCSPH_DBC< Particles, SPHFluidConfig, Variables >::Interaction( NeighborSearchPointer& neighborSearch )
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   GlobalIndexType numberOfParticles = particles->getNumberOfParticles();
   static constexpr GlobalIndexType _numberOfCells = Particles::Config::gridXsize; //FIXIT
   const auto view_firstCellParticle = neighborSearch->getCellFirstParticleList().getView();
   const auto view_lastCellParticle = neighborSearch->getCellLastParticleList().getView(); // DEBUG
   const auto view_particleCellIndex = particles->getParticleCellIndices().getView();
   const auto view_points = particles->getPoints().getView();
   const RealType searchRadius = this->particles->getSearchRadius();

   /* CONSTANT VARIABLES */
   const RealType h = this->h;
   const RealType m = this->m;
   const RealType speedOfSound = this->speedOfSound;
   const RealType coefB = this->coefB;
   const RealType rho0 = this->rho0;
   const RealType delta = this->delta;
   const RealType alpha = this->alpha;

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_particleType = this->getParticleType().getView();
   const auto view_rho = this->getRho().getView();
   auto view_Drho = this->getDrho().getView();
   const auto view_p = this->getPress().getView();
   const auto view_v = this->getVel().getView();
   auto view_a = this->getAcc().getView();

   auto FluidFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, PointType& r_i, PointType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i, PointType* a_i ) mutable
   {
      /* This should be some interaction structure - properties of particle B:*/
      const PointType r_j = view_points[ j ];
      const PointType v_j = view_v[ j ];
      const RealType rho_j = view_rho[ j ];
      const RealType p_j = view_p[ j ];

      /* Interaction: */
      const PointType dr = r_i - r_j;
      const PointType dv = v_i - v_j;

      const RealType drs = l2Norm( dr );
      const RealType F = SPHKernelFunction::F( drs, h );
      const PointType gradW = dr * F;

      const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
      const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
      *drho_i += ( dv, gradW ) * m - diffTerm;

      const RealType p_term = ( p_i + p_j ) / ( rho_i * rho_j );
      const RealType visco =  ViscousTerm::Pi( rho_i, rho_j, drs, ( dr, dv ) );
      *a_i += ( -1.0 ) * ( p_term + visco )* gradW * m;
   };

   auto BoundFluid = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, PointType& r_i, PointType& v_i, RealType& rho_i, RealType& p_i, RealType* drho_i ) mutable
   {
      if( view_particleType[ j ] == 0 ){
         /* This should be some interaction structure, mby. - properties of particle B: */
         const PointType r_j = view_points[ j ];
         const PointType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = view_p[ j ];

         /* Interaction */
         const PointType dr = r_i - r_j;
         const PointType dv = v_i - v_j;

         const RealType drs = l2Norm( dr );
         const RealType F = SPHKernelFunction::F( drs, h );
         const PointType gradW = dr*F;

         const RealType psi = DiffusiveTerm::Psi( rho_i, rho_j, drs );
         const RealType diffTerm =  psi * ( dr, gradW ) * m / rho_j;
         *drho_i += ( dv, gradW ) * m - diffTerm;
      }
   };

   auto particleLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch ) mutable
   {
      const unsigned int activeCell = view_particleCellIndex[ i ];

      /*TODO: This should be some interaction structure  - properties of particle A:*/
      const PointType r_i = view_points[ i ];
      const PointType v_i = view_v[ i ];
      const RealType rho_i = view_rho[ i ];
      const RealType p_i = view_p[ i ];

      PointType a_i = {0., 0.};
      RealType drho_i = 0.;

      // Process fluid particle
      if( view_particleType[ i ] == 0 )
      {
         neighborSearch->loopOverNeighbors( i, numberOfParticles, view_firstCellParticle, view_particleCellIndex, FluidFluid, r_i, v_i, rho_i, p_i, &drho_i, &a_i );

         view_Drho[ i ] = drho_i;
         a_i[ 1 ] -= 9.81 ;
         view_a[ i ] = a_i;
      }
      // Process boundary particle
      else if( view_particleType[ i ] == 1 )
      {
         neighborSearch->loopOverNeighbors( i, numberOfParticles, view_firstCellParticle, view_particleCellIndex, BoundFluid, r_i, v_i, rho_i, p_i, &drho_i );

         view_Drho[ i ] = drho_i;
         a_i = { 0., 0. };
         view_a[ i ] = a_i;
      }
      else
      {
         printf( "INVALID PARTICLE TYPE!" );
      }
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfParticles, particleLoop, neighborSearch );
}


} // SPH
} // ParticleSystem
} // TNL

