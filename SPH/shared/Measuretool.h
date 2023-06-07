#include <TNL/Meshes/Writers/VTKWriter.h>
#include <TNL/Containers/NDArray.h>

#include "../SPHTraits.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename SPHSimulation >
class InterpolateToGrid
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using NeighborSearch = typename SPHSimulation::NeighborSearchType;
   using NeighborSearchPointer = typename SPHSimulation::NeighborSearchPointer;
   using Variables = typename SPHSimulation::FluidVariables;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, DeviceType >;

   using GridType = Meshes::Grid< 2, RealType, DeviceType, GlobalIndexType >;
   using CoordinatesType = typename GridType::CoordinatesType;

   template< typename InterpolationConfigType >
   InterpolateToGrid( InterpolationConfigType configInit )
   : variables( configInit.gridSize[ 0 ] * configInit.gridSize[ 1 ] ),
     gridDimension( configInit.gridSize )
   {
      interpolationGrid.setOrigin( configInit.gridOrigin );
      interpolationGrid.setDimensions( configInit.gridSize );
      interpolationGrid.setSpaceSteps( configInit.gridStep );
   }

   template< typename SPHKernelFunction, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:

   GridType interpolationGrid;
   CoordinatesType gridDimension;
   VariablesPointer variables;


};

template< typename SPHConfig, typename SPHSimulation >
class SensorInterpolation
{
   public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using NeighborSearchPointer = typename SPHSimulation::NeighborSearchPointer;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 DeviceType >;         // Device

   SensorInterpolation( GlobalIndexType numberOfSavedSteps, std::vector< VectorType >& sensorsPoints )
   : sensorPositions( sensorsPoints ),
     numberOfSensors( sensorsPoints.size() ),
     numberOfSavedSteps( numberOfSavedSteps )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );
   }

   template<typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState );

   void
   save( const std::string outputFileName );

protected:

   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;

   GlobalIndexType numberOfSavedSteps;
   GlobalIndexType sensorIndexer = 0;
};

template< typename SPHConfig, typename SPHSimulation >
class SensorWaterLevel
{
   public:
   using DeviceType = typename SPHConfig::DeviceType;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;

   using LocalIndexType = typename SPHTraitsType::LocalIndexType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   using FluidPointer = typename SPHSimulation::FluidPointer;
   using BoundaryPointer = typename SPHSimulation::BoundaryPointer;
   using NeighborSearchPointer = typename SPHSimulation::NeighborSearchPointer;

   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 Devices::Host >;         // Device - store data on HOST

   SensorWaterLevel( GlobalIndexType numberOfSavedSteps,
                     std::vector< VectorType >& sensorsPoints,
                     RealType levelIncrement,
                     VectorType direction,
                     RealType startLevel,
                     RealType endLevel )
   : sensorPositions( sensorsPoints ), numberOfSensors( sensorsPoints.size() ), numberOfSavedSteps( numberOfSavedSteps ),
     levelIncrement( levelIncrement ), direction( direction )
   {
      sensors.setSizes( numberOfSavedSteps, numberOfSensors );

      //New
      numberOfLevels = TNL::ceil( ( endLevel - startLevel ) / levelIncrement );
      levels.setSize( numberOfLevels );
   }

   //protected:
   GlobalIndexType numberOfSensors;
   SensorsDataArray sensors;
   VectorArrayType sensorPositions;

   GlobalIndexType numberOfSavedSteps;
   GlobalIndexType sensorIndexer = 0;

   //New
   IndexArrayType levels;
   GlobalIndexType numberOfLevels; //numberOfLevels = TNL::ceil( TNL::norm( endPoints - startPoint ) / levelIncrement );
   RealType levelIncrement;
   VectorType direction;

   RealType startLevel;
   RealType endLevel;

   template< typename SPHKernelFunction, typename EOS, typename SPHState >
   void
   interpolateSensors( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
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
      const RealType h = sphState.h;
      const RealType m = sphState.mass;

      auto view_sensorsPositions = sensorPositions.getView();
      auto view_sensors = sensors.getView();
      auto view_levels = levels.getView();

      auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* gamma ) mutable
      {
         const VectorType r_j = view_points[ j ];
         const VectorType r_ij = r_i - r_j;
         const RealType drs = l2Norm( r_ij );
         if( drs <= searchRadius )
         {
            const RealType rho_j = view_rho[ j ];
            const RealType W = SPHKernelFunction::W( drs, h );
            const RealType V = m / rho_j;
            *gamma += W * V;
         }
      };

      std::cout << "-----------------------> Number of levels: " << numberOfLevels << std::endl;
      for( int s = 0; s < numberOfSensors; s++ )
      {
         //const VectorType startPoint = view_sensorsPositions[ 0 ];
         const VectorType startPoint = 0.f;
         view_levels = 0.f;

         auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound, GlobalIndexType sensorIndexer ) mutable
         {
            RealType gamma = 0.f;
            const VectorType zax = { 0.f, 1.f };
            //VectorType r = view_sensorsPositions[ 0 ] + i * levelIncrement * zax;
            VectorType r = view_sensorsPositions[ s ] + i * h * zax;
            const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );

            neighborSearch->loopOverNeighbors(
                  i,
                  numberOfParticles,
                  gridIndex,
                  gridSize,
                  view_firstLastCellParticle,
                  interpolate, r, &gamma );

            if( gamma > 0.5f )
               view_levels[ i ] = 1;
            else
               view_levels[ i ] = 0;
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, sensorsLoop, fluid->neighborSearch, boundary->neighborSearch, this->sensorIndexer );

         auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_levels[ i ]; };
         auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
         const GlobalIndexType numberOfFilledLevels = Algorithms::reduce< DeviceType >( 0, view_levels.getSize(), fetch, reduction, 0.0 );
         const RealType waterLevel = h * numberOfFilledLevels + ( startPoint, direction );

         view_sensors( sensorIndexer, s ) = waterLevel;
      }

      sensorIndexer++;
   }

   void
   saveSensors( const std::string outputFileName )
   {
      using HostSensorsDataArray = Containers::NDArray< RealType,
                                                        Containers::SizesHolder< int, 0, 0 >,
                                                        std::index_sequence< 0, 1 >,
                                                        Devices::Host >;
      HostSensorsDataArray sensorsDataHost;
      sensorsDataHost = sensors;
      const auto view_sensorsDataHost = sensorsDataHost.getView();

      std::cout << std::endl << "MEASURETOOL - SAVING DATA :" << std::endl;
      std::cout << "Number of sensors............................. " << numberOfSensors << " ." << std::endl;
      std::cout << "Number of saved time steps.................... " << sensorIndexer << " ." << std::endl;
      std::cout << "Ouput filename................................ " << outputFileName << " ." << std::endl;

      std::ofstream sensorsFile;
      sensorsFile.open( outputFileName );
      sensorsFile << "#Data from measuretool:\n";
      for( int i = 0; i < sensorIndexer; i++ ){
         sensorsFile << i;
         for( int j = 0; j < numberOfSensors; j++ ){
            sensorsFile << " " <<  view_sensorsDataHost( i, j );
         }
         sensorsFile << std::endl;
      }
      sensorsFile.close();
   }

};

} // SPH
} // ParticleSystem
} // TNL

#include "Measuretool.hpp"

