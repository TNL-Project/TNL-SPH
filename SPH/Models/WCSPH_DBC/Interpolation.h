//#include <TNL/Meshes/Writers/VTUWriter.h>
#include <TNL/Meshes/Writers/VTKWriter.h>
//#include <TNL/Meshes/Writers/VTIWriter.h>

#include "../../SPHTraits.h"
#include "../../customParallelFor.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename Variables >
class Interpolation
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;

   //using GridType = typename ParticleTraitsType::GridType;
   using GridType = Meshes::Grid< 2, RealType, typename SPHConfig::DeviceType, GlobalIndexType >;

   using CoordinatesType = typename GridType::CoordinatesType;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, typename SPHConfig::DeviceType >;

   Interpolation(VectorType gridOrigin, CoordinatesType gridDimension, VectorType spaceSteps,
      GlobalIndexType numberOfSavedSteps, GlobalIndexType numberOfDensitySensors )
   : variables( gridDimension[ 0 ] * gridDimension[ 1 ] ), gridDimension( gridDimension )
   {
      interpolationGrid.setOrigin( gridOrigin );
      interpolationGrid.setDimensions( gridDimension );
      interpolationGrid.setSpaceSteps( spaceSteps );

      densitySensors.setSize( numberOfDensitySensors );
      densitySensorsPositions.setSize( numberOfDensitySensors );
   }

   template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
   void InterpolateSensors( FluidPointer& fluid, BoudaryPointer& boundary );

   template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
   void InterpolateGrid( FluidPointer& fluid, BoudaryPointer& boundary );

   void saveInterpolation( std::string outputFileName );

   //protected:
   GridType interpolationGrid;
   CoordinatesType gridDimension;

   VariablesPointer variables;

   GlobalIndexType numberOfSensors;
   //std::vector< VariablesPointer > sensorsHistory;

   using VariablesArrayType = Containers::Array< typename Variables::ScalarArrayType, DeviceType >;
   typename Variables::VectorArrayType densitySensorsPositions;
   VariablesArrayType densitySensors;
   GlobalIndexType numberOfDensitySensors;
   GlobalIndexType densitySensorsTimeIndexer = 0;

};

} // SPH
} // ParticleSystem
} // TNL

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename Variables >
template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
void
Interpolation< SPHConfig, Variables >::InterpolateGrid( FluidPointer& fluid, BoudaryPointer& boundary )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */ //TODO: Do this like a human.
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   const VectorType gridOrigin = fluid->particles->getGridOrigin();
   const IndexVectorType gridSize = fluid->particles->getGridSize();

   const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();

   /* CONSTANT VARIABLES */ //TODO: Do this like a human.
   const RealType h = SPHConfig::h;
   const RealType m = SPHConfig::mass;

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

         *v += v_j * W * V;
         *rho += W * m;

         *gamma += W * V;
      }
   };

   auto gridLoop = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound ) mutable
   {
      VectorType v = 0.f;
      RealType rho = 0.f;
      RealType gamma = 0.f;

      VectorType r = { ( i + 1 ) * searchRadius , ( j + 1) * searchRadius };
      const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );
      const GlobalIndexType idx =  j * 45 + i;
      printf("%d ", idx);

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, interpolate, r, &rho, &v, &gamma );


     if( gamma > 0.5f ){
        view_v_interpolation[ idx ] = v / gamma;
        view_rho_interpolation[ idx ] = rho /gamma;
     }
     else{
        view_v_interpolation[ idx ] = 0.f;
        view_rho_interpolation[ idx ] = 0.f;
     }
   };
   Algorithms::ParallelFor2D< DeviceType >::exec(
      ( LocalIndexType ) 0,
      ( LocalIndexType ) 0,
      ( LocalIndexType ) gridDimension[ 0 ],
      ( LocalIndexType ) gridDimension[ 1 ],
      gridLoop, fluid->neighborSearch, boundary->neighborSearch );
}

template< typename SPHConfig, typename Variables >
void
Interpolation< SPHConfig, Variables >::saveInterpolation( const std::string outputFileName )
{

   std::cout << "Saving interpolation ........... ";
   std::cout << "Interpolation grid: " << std::endl << interpolationGrid << std::endl;
   std::ofstream file( outputFileName );
   using Writer = TNL::Meshes::Writers::VTKWriter< GridType >;
   //using Writer = TNL::Meshes::Writers::VTIWriter< GridType >;
   Writer writer( file );
   writer.writeEntities( interpolationGrid );
   //writer.writeImageData( interpolationGrid );
   //writer.template writeDataArray< typename Variables::ScalarArrayType >( variables->rho.getView(), "Density" );
   std::string densityName = "Density";
   //writer.template writeDataArray< typename Variables::ScalarArrayType >( variables->rho, densityName, 1 );
   writer.template writeCellData< typename Variables::ScalarArrayType >( variables->rho, densityName, 1 );
   //writer.template writeCellData< typename Variables::VectorArrayType >( variables->v, "Velocity", 2 );
   //std::cout << variables->rho << std::endl;
   std::cout << variables->rho.getSize() << std::endl;
   std::cout << "MY size: "<< gridDimension[ 0 ] * gridDimension[ 1 ] << std::endl;
   std::cout << " SAVED." << std::endl;

}


template< typename SPHConfig, typename Variables >
template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
void
Interpolation< SPHConfig, Variables >::InterpolateSensors( FluidPointer& fluid, BoudaryPointer& boundary )
{

   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */ //TODO: Do this like a human.
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   const VectorType gridOrigin = fluid->particles->getGridOrigin();
   const IndexVectorType gridSize = fluid->particles->getGridSize();

   const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();

   /* CONSTANT VARIABLES */ //TODO: Do this like a human.
   const RealType h = SPHConfig::h;
   const RealType m = SPHConfig::m;


   auto view_densitySensorsPositions = densitySensorsPositions.getView();
   auto view_densitySensors = densitySensors.getView();

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );

         const RealType V = m / rho_j;

         rho += W * m;

         gamma += W * V;
      }
   };

   auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound ) mutable
   {
      RealType rho = 0.f;
      VectorType v = 0.f;
      RealType gamma = 0.f;

      VectorType r = view_densitySensorsPositions[ i ];
      const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, interpolate, r, &rho, &v, &gamma );


     auto view_densitySensor = view_densitySensors[ i ].getView();

     if( gamma > 0.5f ){
        view_densitySensor[ densitySensorsTimeIndexer ] = rho /gamma;
     }
     else{
        view_densitySensor[ densitySensorsTimeIndexer ] = 0.f;
     }
   };
   Algorithms::ParallelFor< DeviceType >::exec( 0, numberOfDensitySensors, sensorsLoop, fluid->neighborSearch, boundary->neighborSearch );

   densitySensorsTimeIndexer++;
}

} // SPH
} // ParticleSystem
} // TNL

