#include "../../SPHTraits.h"
#include ""

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename Variables >
class Interpolation
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   template< typename VariableType >
   VariableType InterpolatePoint();

   template< typename VariableType, typename SPHKernelFunction >
   void InterpolateGrid();

   using GridType = typename ParticleTraitsType::GridType;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   //protected:
   /*
      - grid
      - variables (grid size)
    */
   GridType interpolationGrid;
   VariablesPointer interpolation;

};

} // SPH
} // ParticleSystem
} // TNL

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig >
template< typename VariableType, typename SPHKernelFunction >
void
Interpolation< SPHConfig >::InterpolateGrid()
{
   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );

         const RealType V = m / rho_j;

         v += v_j * W * V;
         rho += W * m;

         gamma += W * V;
      }
   };

   auto gridLoop = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j ) mutable
   {
   };
   Algorithms::ParallelFor2D< DeviceType >::exec(
      ( LocalIndexType ) 0,
      ( LocalIndexType ) 0,
      ( LocalIndexType ) gridDimX,
      ( LocalIndexType ) gridDimY,
      f );
}

template< typename SPHConfig >
template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
void
Interpolation< SPHConfig >::InterpolatePoint( const VectorType r, NeighborSearchPointer neighborSearch )
{

   /* CONSTANT VARIABLES */ //TODO: Do this like a human.
   const RealType h = SPHConfig::h;
   const RealType m = SPHConfig::m;

   /* INTERPOLATED VARIABLES */
   VectorType v = 0.f;
   RealType rho = 0.f;

   RealType gamma = 0.f;

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );

         const RealType V = m / rho_j;

         v += v_j * W * V;
         rho += W * m;

         gamma += W * V;
      }
   };
   neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, interpolate, &rho, &v, &gamma );


}

} // SPH
} // ParticleSystem
} // TNL

