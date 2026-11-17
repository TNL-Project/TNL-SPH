// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

namespace TNL::SPH {

template< typename SPHConfig, typename SPHSimulation >
template< typename NeighborLoop, typename ModelParams, typename SourcePointer >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::accumulateGrid(
   const SourcePointer& source,
   const ModelParams& modelParams,
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   using SPHKernelFunction = typename ModelParams::KernelFunction;
   const RealType h = modelParams.h;
   const RealType mass = modelParams.mass;
   const RealType searchRadius = source->getParticles()->getSearchRadius();
   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   // specify the offset to distinguish cell-based and vertex-based interpolation
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
   const VectorType gridOrigin = this->interpolationGrid.getOrigin() + gridOriginOffsetFactor * gridSpaceSteps;

   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const auto view_v = source->getVariables()->v.getConstView();
   auto view_rho_interpolation = variables->rho.getView();
   auto view_v_interpolation = variables->v.getView();
   auto view_gamma_interpolation = gamma.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto interpolate =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = mass / rho_j;
         *rho += W * mass;
         *v += v_j * W * V;
         *gamma += W * V;
      }
   };

   auto processPoint = [ = ] __cuda_callable__( const VectorType r, const GlobalIndexType idx ) mutable
   {
      if( ownersCount != 0 && ownersView[ idx ] != levelSlot )
         return;
      VectorType v = view_v_interpolation[ idx ];
      RealType rho = view_rho_interpolation[ idx ];
      RealType gamma = view_gamma_interpolation[ idx ];
      NeighborLoop::exec( idx, r, searchInSource, interpolate, &rho, &v, &gamma );
      view_v_interpolation[ idx ] = v;
      view_rho_interpolation[ idx ] = rho;
      view_gamma_interpolation[ idx ] = gamma;
   };

   if( interpolatedGridEntity == 0 )
      interpolationGrid.template forAllEntities< 0 >(
         [ = ] __cuda_callable__( const GridVertex& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            processPoint( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );
   else
      interpolationGrid.template forAllEntities< SPHConfig::spaceDimension >(
         [ = ] __cuda_callable__( const GridCell& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            processPoint( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );
}

template< typename SPHConfig, typename SPHSimulation >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::normalizeGrid(
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   const GlobalIndexType entitiesCount = this->interpolationGrid.getEntitiesCount( this->interpolatedGridEntity );
   auto view_rho_interpolation = variables->rho.getView();
   auto view_v_interpolation = variables->v.getView();
   auto view_gamma_interpolation = gamma.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto normalize = [ = ] __cuda_callable__( GlobalIndexType idx ) mutable
   {
      if( ownersCount != 0 && ownersView[ idx ] != levelSlot )
         return;
      const RealType gamma = view_gamma_interpolation[ idx ];
      if( gamma > 0.5f ) {
         view_v_interpolation[ idx ] /= gamma;
         view_rho_interpolation[ idx ] /= gamma;
      }
      else {
         view_v_interpolation[ idx ] = 0.f;
         view_rho_interpolation[ idx ] = 0.f;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, entitiesCount, normalize );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingGrid(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   const ModelParams& modelParams )
{
   (void) boundary;
   this->variables->rho = 0.f;
   this->variables->v = VectorType( 0.f );
   this->gamma = 0.f;

   accumulateGrid< typename ParticlesType::NeighborsLoop >( fluid, modelParams );
   normalizeGrid();
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingParallelFor(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   const ModelParams& modelParams )
{
   using SPHKernelFunction = typename ModelParams::KernelFunction;

   // neighbor search objects
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->getParticles() );

   // loaded arrays
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   const auto view_v = fluid->getVariables()->v.getView();

   auto view_rho_interpolation = this->variables()->rho.getView();
   auto view_v_interpolation = this->variables()->v.getView();

   // loaded constants
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();

   // specify the offset to distinguish cell-based and vertex-based interpolation
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
   const CoordinatesType gridDimensionsOffset = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 1.0f );

   // interpolation grid constatns
   const CoordinatesType dimensions = this->interpolationGrid.getDimensions() + gridDimensionsOffset;
   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   const VectorType gridOrigin = this->interpolationGrid.getOrigin() + gridOriginOffsetFactor * gridSpaceSteps;

   auto interpolate =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType & r_i, RealType * rho, VectorType * v, RealType * gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const VectorType v_j = view_v[ j ];
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );

         const RealType V = m / rho_j;

         *v += v_j * W * V;
         *rho += W * m;
         *gamma += W * V;
      }
   };

   auto loopOverGrid = [ = ] __cuda_callable__( const CoordinatesType& i ) mutable
   {
      const VectorType r = i * gridSpaceSteps + gridOrigin;
      const GlobalIndexType idx = CellIndexer::evaluateCellIndex( i, dimensions );

      VectorType v = 0.f;
      RealType rho = 0.f;
      RealType gamma = 0.f;

      ParticlesType::NeighborsLoop::exec( idx, r, searchInFluid, interpolate, &rho, &v, &gamma );

      if( gamma > 0.5f ) {
         view_v_interpolation[ idx ] = v / gamma;
         view_rho_interpolation[ idx ] = rho / gamma;
      }
      else {
         view_v_interpolation[ idx ] = 0.f;
         view_rho_interpolation[ idx ] = 0.f;
      }
   };
   const CoordinatesType beg = 0;
   Algorithms::parallelFor< DeviceType >( beg, dimensions, interpolate );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   const ModelParams& modelParams )
{
   interpolateUsingGrid( fluid, boundary, modelParams );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams, typename LevelView >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   const ModelParams& modelParams )
{
   if( ! validateMeasurementLevels( levelViews ) )
      throw std::invalid_argument( "Invalid multiresolution measurement levels" );

   const GlobalIndexType entitiesCount = this->interpolationGrid.getEntitiesCount( this->interpolatedGridEntity );
   this->variables->rho = 0.f;
   this->variables->v = VectorType( 0.f );
   this->gamma = 0.f;

   auto ownerBoxes = makeMeasurementOwnerBoxArray< SPHConfig >( levelViews );
   IndexArrayType owners;
   owners.setSize( entitiesCount );
   const auto ownerBoxesView = ownerBoxes.getConstView();
   auto ownersView = owners.getView();

   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
   const VectorType gridOrigin = this->interpolationGrid.getOrigin() + gridOriginOffsetFactor * gridSpaceSteps;

   auto assignOwner = [ = ] __cuda_callable__( const VectorType& r, const GlobalIndexType idx ) mutable
   {
      ownersView[ idx ] = findOwningMeasurementLevel( r, ownerBoxesView );
   };
   if( this->interpolatedGridEntity == 0 )
      this->interpolationGrid.template forAllEntities< 0 >(
         [ = ] __cuda_callable__( const GridVertex& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            assignOwner( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );
   else
      this->interpolationGrid.template forAllEntities< SPHConfig::spaceDimension >(
         [ = ] __cuda_callable__( const GridCell& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            assignOwner( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );

   for( std::size_t level = 0; level < levelViews.size(); ++level ) {
      const int levelSlot = static_cast< int >( level );
      auto levelParams = modelParams;
      levelParams.h = levelViews[ level ].h;
      levelParams.mass = levelViews[ level ].mass;

      accumulateGrid< typename ParticlesType::NeighborsLoopAnotherSet >(
         *levelViews[ level ].fluid, levelParams, owners, levelSlot );
      for( const auto* patch : levelViews[ level ].ownPatches )
         accumulateGrid< typename ParticlesType::NeighborsLoopAnotherSet >( *patch, levelParams, owners, levelSlot );

      normalizeGrid( owners, levelSlot );
   }
}

template< typename SPHConfig, typename SPHSimulation >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   // write interpolation grid
   using Writer = TNL::Meshes::Writers::VTKWriter< GridType >;
   std::ofstream file( outputFileName );
   Writer writer( file );
   writer.writeEntities( interpolationGrid );

   //FIXME: Workaround for vectorArray;
   const int entitiesCount = this->interpolationGrid.getEntitiesCount( interpolatedGridEntity );
   using BufferType = Containers::Array< RealType, Devices::Host, GlobalIndexType >;
   BufferType buffer( 3 * entitiesCount );
   const auto velocityView = variables->v.getView();
   GlobalIndexType k = 0;
   for( GlobalIndexType i = 0; i < entitiesCount; i++ ) {
      const VectorType vector = velocityView.getElement( i );
      for( int j = 0; j < 3; j++ )
         buffer[ k++ ] = ( j < vector.getSize() ? vector[ j ] : 0 );
   }

   // write data
   if( this->interpolatedGridEntity > 0 ) {
      writer.template writeCellData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );
      writer.template writeCellData< typename Variables::ScalarArrayType >( gamma, "Gamma", 1 );
      writer.template writeCellData( buffer, "Velocity", 3 );
   }
   else {
      writer.template writePointData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );
      writer.template writePointData< typename Variables::ScalarArrayType >( gamma, "Gamma", 1 );
      writer.template writePointData( buffer, "Velocity", 3 );
   }
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams, typename SourcePointer >
void
SensorInterpolation< SPHConfig, SPHSimulation >::accumulatePressure(
   const SourcePointer& source,
   const ModelParams& modelParams,
   const RealType searchRadius,
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   using SPHKernelFunction = typename ModelParams::KernelFunction;
   using EOS = typename ModelParams::EOS;
   const RealType h = modelParams.h;
   const RealType mass = modelParams.mass;
   const typename EOS::ParamsType eosParams( modelParams );

   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const auto sensorPositionsView = sensorPositions.getConstView();
   auto pressureView = pressure.getView();
   auto gammaView = gamma.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto interpolate =
      [ = ] __cuda_callable__( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* p, RealType* g ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = mass / rho_j;
         *p += p_j * W * V;
         *g += W * V;
      }
   };

   auto sensorsLoop = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( ownersCount != 0 && ownersView[ i ] != levelSlot )
         return;
      VectorType r = sensorPositionsView[ i ];
      RealType p = pressureView[ i ];
      RealType g = gammaView[ i ];
      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInSource, interpolate, &p, &g );
      pressureView[ i ] = p;
      gammaView[ i ] = g;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, sensorsLoop );
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorInterpolation< SPHConfig, SPHSimulation >::resetPressure( const GlobalIndexType outputRow )
{
   if( pressure.getSize() != numberOfSensors )
      pressure.setSize( numberOfSensors );
   if( gamma.getSize() != numberOfSensors )
      gamma.setSize( numberOfSensors );
   pressure = 0.f;
   gamma = 0.f;

   auto sensorsView = sensors.getView();
   auto resetRow = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      sensorsView( outputRow, i ) = 0.f;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, resetRow );
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorInterpolation< SPHConfig, SPHSimulation >::normalizePressure(
   const GlobalIndexType outputRow,
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   const auto pressureView = pressure.getConstView();
   const auto gammaView = gamma.getConstView();
   auto sensorsView = sensors.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto normalize = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( ownersCount != 0 && ownersView[ i ] != levelSlot )
         return;
      const RealType g = gammaView[ i ];
      if( g > 0.5f )
         sensorsView( outputRow, i ) = pressureView[ i ] / g;
      else
         sensorsView( outputRow, i ) = 0.f;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, normalize );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   const ModelParams& modelParams )
{
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const bool includeBoundary = this->includeBoundary;
   const GlobalIndexType outputRow = this->sensorIndexer;

   resetPressure( outputRow );
   accumulatePressure( fluid, modelParams, searchRadius );
   if( includeBoundary )
      accumulatePressure( boundary, modelParams, searchRadius );

   normalizePressure( outputRow );
   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams, typename LevelView >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   const ModelParams& modelParams )
{
   if( ! validateMeasurementLevels( levelViews ) )
      throw std::invalid_argument( "Invalid multiresolution measurement levels" );

   IndexArrayType owners;
   auto ownerBoxes = makeMeasurementOwnerBoxArray< SPHConfig >( levelViews );
   computeMeasurementOwners( sensorPositions, ownerBoxes, owners );

   const GlobalIndexType outputRow = this->sensorIndexer;
   resetPressure( outputRow );

   const bool includeBoundary = this->includeBoundary;
   for( std::size_t level = 0; level < levelViews.size(); ++level ) {
      const int levelSlot = static_cast< int >( level );
      auto levelParams = modelParams;
      levelParams.h = levelViews[ level ].h;
      levelParams.mass = levelViews[ level ].mass;

      accumulatePressure(
         *levelViews[ level ].fluid,
         levelParams,
         ( *levelViews[ level ].fluid )->getParticles()->getSearchRadius(),
         owners,
         levelSlot );
      for( const auto* patch : levelViews[ level ].ownPatches )
         accumulatePressure(
            *patch, levelParams, ( *patch )->getParticles()->getSearchRadius(), owners, levelSlot );
      if( includeBoundary )
         accumulatePressure(
            *levelViews[ level ].boundary,
            levelParams,
            ( *levelViews[ level ].boundary )->getParticles()->getSearchRadius(),
            owners,
            levelSlot );

      normalizePressure( outputRow, owners, levelSlot );
   }

   ++this->sensorIndexer;
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorInterpolation< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   using HostSensorsDataArray =
      Containers::NDArray< RealType, Containers::SizesHolder< int, 0, 0 >, std::index_sequence< 0, 1 >, Devices::Host >;
   HostSensorsDataArray sensorsDataHost;
   sensorsDataHost = sensors;
   const auto view_sensorsDataHost = sensorsDataHost.getView();

   std::ofstream sensorsFile;
   sensorsFile.open( outputFileName );
   sensorsFile << "#Data from measuretool:\n";
   for( int i = 0; i < sensorIndexer; i++ ) {
      sensorsFile << i;
      for( int j = 0; j < numberOfSensors; j++ ) {
         sensorsFile << " " << view_sensorsDataHost( i, j );
      }
      sensorsFile << std::endl;
   }
   sensorsFile.close();
}

//template< typename SPHConfig, typename SPHSimulation >
//template<typename ModelParams >
//void
//SensorGeneralInterpolation< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid,
//                                                                     BoundaryPointer& boundary,
//                                                                     ModelParams& modelParams,
//                                                                     bool includeBoundary )
//{
//   using SPHKernelFunction = typename ModelParams::KernelFunction;
//   using EOS = typename ModelParams::EOS;
//
//   GlobalIndexType numberOfParticles = fluid->getParticles()->getNumberOfParticles();
//   GlobalIndexType numberOfParticles_bound = boundary->getParticles()->getNumberOfParticles();
//   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
//
//   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->getParticles() );
//   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->getParticles() );
//
//   /* VARIABLES AND FIELD ARRAYS */
//   const auto view_points = fluid->getParticles()->getPoints().getView();
//   const auto view_rho = fluid->getVariables()->rho.getView();
//
//   const auto view_points_boundary = boundary->getParticles()->getPoints().getView();
//   const auto view_rho_boundary = boundary->getVariables()->rho.getView();
//
//   const auto view_field = field.getView();
//
//   /* CONSTANT VARIABLES */
//   const RealType h = modelParams.h;
//   const RealType m = modelParams.mass;
//
//   typename EOS::ParamsType eosParams( modelParams );
//
//   auto view_sensorsPositions = sensorPositions.getView();
//   auto view_pressureSensors = sensors.getView();
//
//   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
//         VectorType& r_i, RealType* p, VectorType* v, RealType* gamma ) mutable
//   {
//      const VectorType r_j = view_points[ j ];
//      const VectorType r_ij = r_i - r_j;
//      const RealType drs = l2Norm( r_ij );
//      if( drs <= searchRadius )
//      {
//         const RealType rho_j = view_rho[ j ];
//         const RealType value_j = view_field[ j ];
//         const RealType W = SPHKernelFunction::W( drs, h );
//         const RealType V = m / rho_j;
//
//         *value += value_j * W * V;
//         *gamma += W * V;
//      }
//   };
//
//   auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, GlobalIndexType sensorIndexer ) mutable
//   {
//      VectorType r = view_sensorsPositions[ i ];
//      ValueType value = 0.f;
//      RealType gamma = 0.f;
//
//      NeighborsLoop::exec( i, r, searchInFluid, interpolate, &value, &gamma );
//      //if( includeBoundary ){
//      //   NeighborsLoop::exec( i, r, searchInBound, interpolateBoundary, &value, &gamma );
//      //}
//
//      if( gamma > 0.5f ){
//         view_pressureSensors( sensorIndexer, i ) = value / gamma;
//      }
//      else{
//         view_pressureSensors( sensorIndexer, i ) = 0.f;
//      }
//   };
//   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, sensorsLoop, this->sensorIndexer );
//
//   sensorIndexer++;
//}
//
//template< typename SPHConfig, typename SPHSimulation >
//void
//SensorGeneralInterpolation< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
//{
//   using HostSensorsDataArray = Containers::NDArray< RealType,
//                                                     Containers::SizesHolder< int, 0, 0 >,
//                                                     std::index_sequence< 0, 1 >,
//                                                     Devices::Host >;
//   HostSensorsDataArray sensorsDataHost;
//   sensorsDataHost = sensors;
//   const auto view_sensorsDataHost = sensorsDataHost.getView();
//
//   std::cout << std::endl << "MEASURETOOL - SAVING DATA :" << std::endl;
//   std::cout << "Number of sensors............................. " << numberOfSensors << " ." << std::endl;
//   std::cout << "Number of saved time steps.................... " << sensorIndexer << " ." << std::endl;
//   std::cout << "Ouput filename................................ " << outputFileName << " ." << std::endl;
//
//   std::ofstream sensorsFile;
//   sensorsFile.open( outputFileName );
//   sensorsFile << "#Data from measuretool:\n";
//   for( int i = 0; i < sensorIndexer; i++ ){
//      sensorsFile << i;
//      for( int j = 0; j < numberOfSensors; j++ ){
//         sensorsFile << " " <<  view_sensorsDataHost( i, j );
//      }
//      sensorsFile << std::endl;
//   }
//   sensorsFile.close();
//}

template< typename SPHConfig, typename SPHSimulation >
const typename SensorWaterLevel< SPHConfig, SPHSimulation >::SensorsDataArray&
SensorWaterLevel< SPHConfig, SPHSimulation >::getSensorData() const
{
   return sensors;
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams, typename SourcePointer >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::accumulateWater(
   const SourcePointer& source,
   const ModelParams& modelParams,
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   using SPHKernelFunction = typename ModelParams::KernelFunction;
   const RealType h = modelParams.h;
   const RealType mass = modelParams.mass;
   const RealType searchRadius = source->getParticles()->getSearchRadius();

   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const auto samplePositionsView = samplePositions.getConstView();
   auto sampleGammaView = sampleGamma.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto interpolateWaterLevel =
      [ = ] __cuda_callable__( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = mass / rho_j;
         *gamma += W * V;
      }
   };

   auto samplesLoop = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( ownersCount != 0 && ownersView[ i ] != levelSlot )
         return;
      VectorType r = samplePositionsView[ i ];
      RealType gamma = sampleGammaView[ i ];
      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInSource, interpolateWaterLevel, &gamma );
      sampleGammaView[ i ] = gamma;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, samplesLoop );
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::resetWaterSamples( const VectorType& sensorLocation )
{
   if( samplePositions.getSize() != numberOfLevels )
      samplePositions.setSize( numberOfLevels );
   if( sampleGamma.getSize() != numberOfLevels )
      sampleGamma.setSize( numberOfLevels );
   if( levels.getSize() != numberOfLevels )
      levels.setSize( numberOfLevels );
   sampleGamma = 0.f;
   levels = 0;

   const RealType sampleLevelIncrement = levelIncrement;
   const VectorType sampleDirection = direction;
   auto samplePositionsView = samplePositions.getView();
   auto setSamplePosition = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      samplePositionsView[ i ] = sensorLocation + i * sampleLevelIncrement * sampleDirection;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, setSamplePosition );
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::classifyWaterSamples(
   const IndexArrayType& owners,
   const GlobalIndexType levelSlot )
{
   const auto sampleGammaView = sampleGamma.getConstView();
   auto levelsView = levels.getView();
   const auto ownersView = owners.getConstView();
   const GlobalIndexType ownersCount = owners.getSize();

   auto classifySample = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( ownersCount == 0 || ownersView[ i ] == levelSlot )
         levelsView[ i ] = sampleGammaView[ i ] > 0.5f ? 1 : 0;
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, classifySample );
}

template< typename SPHConfig, typename SPHSimulation >
typename SensorWaterLevel< SPHConfig, SPHSimulation >::RealType
SensorWaterLevel< SPHConfig, SPHSimulation >::reduceWaterLevel( const VectorType& sensorLocation )
{
   const auto levelsView = levels.getConstView();
   auto fetch = [ = ] __cuda_callable__( GlobalIndexType i ) -> GlobalIndexType
   {
      return levelsView[ i ];
   };
   auto reduction = [] __cuda_callable__( const GlobalIndexType& a, const GlobalIndexType& b )
   {
      return a + b;
   };
   const GlobalIndexType numberOfFilledLevels =
      Algorithms::reduce< DeviceType >( 0, levelsView.getSize(), fetch, reduction, 0 );
   // TNL overloads the comma operator as a dot product.
   return levelIncrement * numberOfFilledLevels + ( sensorLocation, direction );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   const ModelParams& modelParams )
{
   (void) boundary;
   const GlobalIndexType outputRow = this->sensorIndexer;

   for( GlobalIndexType sensor = 0; sensor < this->numberOfSensors; ++sensor ) {
      const VectorType sensorLocation = sensorPositions.getElement( sensor );
      resetWaterSamples( sensorLocation );
      accumulateWater( fluid, modelParams );
      classifyWaterSamples();
      sensors( outputRow, sensor ) = reduceWaterLevel( sensorLocation );
   }

   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
template< typename ModelParams, typename LevelView >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   const ModelParams& modelParams )
{
   if( ! validateMeasurementLevels( levelViews ) )
      throw std::invalid_argument( "Invalid multiresolution measurement levels" );

   auto ownerBoxes = makeMeasurementOwnerBoxArray< SPHConfig >( levelViews );
   IndexArrayType sampleOwners;
   const GlobalIndexType outputRow = this->sensorIndexer;

   for( GlobalIndexType sensor = 0; sensor < this->numberOfSensors; ++sensor ) {
      const VectorType sensorLocation = sensorPositions.getElement( sensor );
      resetWaterSamples( sensorLocation );
      computeMeasurementOwners( samplePositions, ownerBoxes, sampleOwners );

      for( std::size_t level = 0; level < levelViews.size(); ++level ) {
         const int levelSlot = static_cast< int >( level );
         auto levelParams = modelParams;
         levelParams.h = levelViews[ level ].h;
         levelParams.mass = levelViews[ level ].mass;

         accumulateWater( *levelViews[ level ].fluid, levelParams, sampleOwners, levelSlot );
         for( const auto* patch : levelViews[ level ].ownPatches )
            accumulateWater( *patch, levelParams, sampleOwners, levelSlot );

         classifyWaterSamples( sampleOwners, levelSlot );
      }

      sensors( outputRow, sensor ) = reduceWaterLevel( sensorLocation );
   }

   ++this->sensorIndexer;
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   using HostSensorsDataArray =
      Containers::NDArray< RealType, Containers::SizesHolder< int, 0, 0 >, std::index_sequence< 0, 1 >, Devices::Host >;
   HostSensorsDataArray sensorsDataHost;
   sensorsDataHost = sensors;
   const auto view_sensorsDataHost = sensorsDataHost.getView();

   std::ofstream sensorsFile;
   sensorsFile.open( outputFileName );
   sensorsFile << "#Data from measuretool:\n";
   for( int i = 0; i < sensorIndexer; i++ ) {
      sensorsFile << i;
      for( int j = 0; j < numberOfSensors; j++ ) {
         sensorsFile << " " << view_sensorsDataHost( i, j );
      }
      sensorsFile << std::endl;
   }
   sensorsFile.close();
}

}  // namespace TNL::SPH
