//#include <TNL/Meshes/Writers/VTUWriter.h>
#include <TNL/Meshes/Writers/VTKWriter.h>
//#include <TNL/Meshes/Writers/VTIWriter.h>
#include <TNL/Containers/NDArray.h>

#include "../../SPHTraits.h"

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

   using GridType = Meshes::Grid< 2, RealType, typename SPHConfig::DeviceType, GlobalIndexType >;

   using CoordinatesType = typename GridType::CoordinatesType;
   using VariablesPointer = typename Pointers::SharedPointer< Variables, typename SPHConfig::DeviceType >;

   Interpolation( VectorType gridOrigin, CoordinatesType gridDimension, VectorType spaceSteps, //grid
                  GlobalIndexType numberOfSavedSteps, GlobalIndexType numberOfDensitySensors, typename SPHTraitsType::VectorArrayType& sensorsPoints ) //sensors
   : variables( gridDimension[ 0 ] * gridDimension[ 1 ] ), gridDimension( gridDimension ), sensorPositions( numberOfDensitySensors )
   {
      interpolationGrid.setOrigin( gridOrigin );
      interpolationGrid.setDimensions( gridDimension );
      interpolationGrid.setSpaceSteps( spaceSteps );

      //test
      sensorsDensity.setSizes( numberOfSavedSteps, numberOfDensitySensors );
      numberOfSensors = numberOfDensitySensors;
      numberOfSavedSteps = numberOfSavedSteps;
      sensorPositions = sensorsPoints;

   }

   template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer, typename EOS >
   void InterpolateSensors( FluidPointer& fluid, BoudaryPointer& boundary );

   template< typename FluidPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
   void SensorWaterLevel( FluidPointer& fluid );

   template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
   void InterpolateGrid( FluidPointer& fluid, BoudaryPointer& boundary );

   void saveInterpolation( std::string outputFileName );
   void saveSensors( std::string outputFileName );

   //protected:
   GridType interpolationGrid;
   CoordinatesType gridDimension;

   VariablesPointer variables;

   //*********************SENSORS************************

   typename Variables::ScalarArrayType sensorDensity;
   using SensorsDataArray = Containers::NDArray< RealType,  // Value
                                                 Containers::SizesHolder< int, 0, 0 >,     // SizesHolder
                                                 std::index_sequence< 0, 1 >,  // Permutation
                                                 DeviceType >;         // Device

   SensorsDataArray sensorsDensity;
   //typename Variables::VectorArrayType sensorPositions = { 0.01f, 0.01f, 0.1f, 0.1f };
   typename Variables::VectorArrayType sensorPositions;
   GlobalIndexType sensorIndexer = 0;
   GlobalIndexType numberOfSensors;
   GlobalIndexType numberOfSavedSteps;

   ////*********************WATER-LEVEL-SENSORS************************
   ///*
   // * Init with  - Number of sensors
   // *            - Start + end position of the sensor.
   // *
   // *
   // */

   //SensorsDataArray sensorsWaterLevel;
   //typename Variables::VectorArrayType sensorWaterLevelPositions;
   //typename Variables::VectorArrayType sensorWaterLevelPositions_end;
   //GlobalIndexType numberOfWaterLevelSensors;

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

   auto gridLoop = [=] __cuda_callable__ ( const IndexVectorType& i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound ) mutable
   {
      VectorType v = 0.f;
      RealType rho = 0.f;
      RealType gamma = 0.f;

      VectorType r = { ( i + 1 ) * searchRadius , ( j + 1) * searchRadius };
      const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );
      const GlobalIndexType idx =  j * 45 + i;

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
   IndexVectorType begin{ 0, 0 };
   Algorithms::parallelFor< DeviceType >( begin, gridDimension, gridLoop, fluid->neighborSearch, boundary->neighborSearch );
}

template< typename SPHConfig, typename Variables >
void
Interpolation< SPHConfig, Variables >::saveInterpolation( const std::string outputFileName )
{
   using Writer = TNL::Meshes::Writers::VTKWriter< GridType >;

   std::cout << std::endl << "MEASURETOOL - SAVING INTERPOLATION :" << std::endl;
   std::cout << "Ouput filename................................ " << outputFileName << " ." << std::endl;
   std::cout << "Grid origin................................... " << outputFileName << " ." << std::endl;

   std::ofstream file( outputFileName );
   Writer writer( file );
   writer.writeEntities( interpolationGrid );
   //writer.writeImageData( interpolationGrid );


   //writer.template writeDataArray< typename Variables::ScalarArrayType >( variables->rho, densityName, 1 );
   writer.template writeCellData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );
   //writer.template writeCellData< typename Variables::VectorArrayType >( variables->v, "Velocity", 2 );

   //std::cout << variables->rho.getSize() << std::endl;
   //std::cout << "MY size: "<< gridDimension[ 0 ] * gridDimension[ 1 ] << std::endl;

}

template< typename SPHConfig, typename Variables >
void
Interpolation< SPHConfig, Variables >::saveSensors( const std::string outputFileName )
{
   using HostSensorsDataArray = Containers::NDArray< RealType,
                                                     Containers::SizesHolder< int, 0, 0 >,
                                                     std::index_sequence< 0, 1 >,
                                                     Devices::Host >;
   HostSensorsDataArray sensorsDataHost;
   sensorsDataHost = sensorsDensity;
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


template< typename SPHConfig, typename Variables >
template< typename FluidPointer, typename BoudaryPointer, typename SPHKernelFunction, typename NeighborSearchPointer, typename EOS >
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
   const RealType m = SPHConfig::mass;

   auto view_sensorsPositions = sensorPositions.getView();
   auto view_densitySensors = sensorsDensity.getView();

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j );
         const RealType W = SPHKernelFunction::W( drs, h );

         const RealType V = m / rho_j;

         //*rho += W * m;
         *rho += p_j * W * V;
         *gamma += W * V;

      }
   };

   auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, NeighborSearchPointer& neighborSearch_bound, GlobalIndexType sensorIndexer ) mutable
   {
      RealType rho = 0.f;
      VectorType v = 0.f;
      RealType gamma = 0.f;

      VectorType r = view_sensorsPositions[ i ];
      const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );

      neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, interpolate, r, &rho, &v, &gamma );

      if( gamma > 0.5f ){
         view_densitySensors( sensorIndexer, i ) = rho / gamma;
      }
      else{
         view_densitySensors( sensorIndexer, i ) = 0.f;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, 2, sensorsLoop, fluid->neighborSearch, boundary->neighborSearch, this->sensorIndexer );

   sensorIndexer++;
}

//: template< typename SPHConfig, typename Variables >
//: template< typename FluidPointer, typename SPHKernelFunction, typename NeighborSearchPointer >
//: void
//: Interpolation< SPHConfig, Variables >::SensorWaterLevel( FluidPointer& fluid )
//: {
//:
//:    /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */ //TODO: Do this like a human.
//:    GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
//:    const RealType searchRadius = fluid->particles->getSearchRadius();
//:
//:    const VectorType gridOrigin = fluid->particles->getGridOrigin();
//:    const IndexVectorType gridSize = fluid->particles->getGridSize();
//:
//:    const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
//:    const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();
//:
//:    /* VARIABLES AND FIELD ARRAYS */
//:    const auto view_points = fluid->particles->getPoints().getView();
//:    const auto view_rho = fluid->variables->rho.getView();
//:
//:    /* CONSTANT VARIABLES */ //TODO: Do this like a human.
//:    const RealType h = SPHConfig::h;
//:    const RealType m = SPHConfig::mass;
//:
//:    auto view_sensorsPositions = sensorPositions.getView();
//:    auto view_densitySensors = sensorsDensity.getView();
//:
//:    auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* gamma ) mutable
//:    {
//:       const VectorType r_j = view_points[ j ];
//:       const VectorType r_ij = r_i - r_j;
//:       const RealType drs = l2Norm( r_ij );
//:       if( drs <= searchRadius )
//:       {
//:          const RealType rho_j = view_rho[ j ];
//:          const RealType p_j = EOS::DensityToPressure( rho_j );
//:          const RealType W = SPHKernelFunction::W( drs, h );
//:
//:          const RealType V = m / rho_j;
//:
//:          *gamma += W * V;
//:
//:       }
//:    };
//:
//:    for( int s = 0; s < numberOfWaterLevelSensors; s++ )
//:    {
//:       auto pointsInSensorLoop = [=] __cuda_callable__ ( LocalIndexType i, NeighborSearchPointer& neighborSearch, GlobalIndexType sensorIndexer ) mutable
//:       {
//:          RealType rho = 0.f;
//:          VectorType v = 0.f;
//:          RealType gamma = 0.f;
//:
//:          VectorType r = view_sensorsPositions[ i ];
//:          const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );
//:
//:          neighborSearch->loopOverNeighbors( i, numberOfParticles, gridIndex, gridSize, view_firstLastCellParticle, view_particleCellIndex, interpolate, r, &gamma );
//:
//:          if( gamma > 0.5f ){
//:             view_densitySensors( sensorIndexer, i ) = rho / gamma;
//:          }
//:          else{
//:             view_densitySensors( sensorIndexer, i ) = 0.f;
//:          }
//:       };
//:       Algorithms::ParallelFor< DeviceType >::exec( 0, 2, sensorsLoop, fluid->neighborSearch, this->sensorIndexer );
//:    }
//:
//:    sensorIndexer++;
//: }

} // SPH
} // ParticleSystem
} // TNL

