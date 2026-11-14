// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

namespace TNL::SPH {

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingGrid(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& modelParams )
{
   // neighbor search objects
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->getParticles() );

   // loaded arrays
   const auto view_points = fluid->getParticles()->getPoints().getView();
   const auto view_rho = fluid->getVariables()->rho.getView();
   const auto view_v = fluid->getVariables()->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();
   auto view_gamma_interpolation = this->gamma.getView();

   // loaded constants
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();

   // interpolation grid constatns
   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   // specify the offset to distinguish cell-based and vertex-based interpolation
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
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

   auto processPoint = [ = ] __cuda_callable__( const VectorType r, const GlobalIndexType idx ) mutable
   {
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
      view_gamma_interpolation[ idx ] = gamma;
   };

   if( this->interpolatedGridEntity == 0 )
      this->interpolationGrid.template forAllEntities< 0 >(
         [ = ] __cuda_callable__( const GridVertex& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            const VectorType r = coords * gridSpaceSteps + gridOrigin;
            const GlobalIndexType idx = gridEntity.getIndex();
            processPoint( r, idx );
         } );
   else
      this->interpolationGrid.template forAllEntities< SPHConfig::spaceDimension >(
         [ = ] __cuda_callable__( const GridCell& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            const VectorType r = coords * gridSpaceSteps + gridOrigin;
            const GlobalIndexType idx = gridEntity.getIndex();
            processPoint( r, idx );
         } );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingParallelFor(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& modelParams )
{
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
template< typename SPHKernelFunction, typename ModelParams >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& modelParams )
{
   interpolateUsingGrid< SPHKernelFunction, ModelParams >( fluid, boundary, modelParams );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SourcePointer >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::accumulateGridSource(
   const SourcePointer& source,
   const IndexArrayType& owners,
   const GlobalIndexType levelId,
   const RealType h,
   const RealType m,
   const VectorType& gridOrigin,
   const VectorType& gridSpaceSteps )
{
   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const auto view_v = source->getVariables()->v.getConstView();
   const RealType searchRadius = source->getParticles()->getSearchRadius();
   const auto ownersView = owners.getConstView();
   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();
   auto view_gamma_interpolation = this->gamma.getView();

   auto interpolate =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;
         *v += view_v[ j ] * W * V;
         *rho += W * m;
         *gamma += W * V;
      }
   };

   auto processPoint = [ = ] __cuda_callable__( const VectorType r, const GlobalIndexType idx ) mutable
   {
      if( ownersView[ idx ] != levelId )
         return;
      VectorType v = view_v_interpolation[ idx ];
      RealType rho = view_rho_interpolation[ idx ];
      RealType gamma = view_gamma_interpolation[ idx ];
      ParticlesType::NeighborsLoopAnotherSet::exec( idx, r, searchInSource, interpolate, &rho, &v, &gamma );
      view_v_interpolation[ idx ] = v;
      view_rho_interpolation[ idx ] = rho;
      view_gamma_interpolation[ idx ] = gamma;
   };

   if( this->interpolatedGridEntity == 0 )
      this->interpolationGrid.template forAllEntities< 0 >(
         [ = ] __cuda_callable__( const GridVertex& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            processPoint( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );
   else
      this->interpolationGrid.template forAllEntities< SPHConfig::spaceDimension >(
         [ = ] __cuda_callable__( const GridCell& gridEntity ) mutable
         {
            const CoordinatesType coords = gridEntity.getCoordinates();
            processPoint( coords * gridSpaceSteps + gridOrigin, gridEntity.getIndex() );
         } );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename ModelParams, typename LevelView >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   ModelParams& modelParams )
{
   (void) modelParams;
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

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();
   auto view_gamma_interpolation = this->gamma.getView();

   for( std::size_t level = 0; level < levelViews.size(); ++level ) {
      const int levelId = static_cast< int >( level );
      const RealType h = levelViews[ level ].h;
      const RealType m = levelViews[ level ].mass;

      accumulateGridSource< SPHKernelFunction >(
         *levelViews[ level ].fluid, owners, levelId, h, m, gridOrigin, gridSpaceSteps );
      for( const auto* patch : levelViews[ level ].ownPatches )
         accumulateGridSource< SPHKernelFunction >( *patch, owners, levelId, h, m, gridOrigin, gridSpaceSteps );

      auto normalize = [ = ] __cuda_callable__( GlobalIndexType idx ) mutable
      {
         if( ownersView[ idx ] != levelId )
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
template< typename SPHKernelFunction, typename EOS, typename ModelParams >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->getParticles() );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getConstView();
   const auto view_rho = fluid->getVariables()->rho.getConstView();
   const auto view_points_boundary = boundary->getParticles()->getPoints().getConstView();
   const auto view_rho_boundary = boundary->getVariables()->rho.getConstView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   //const RealType h = modelParams.h;
   //const RealType m = modelParams.mass;
   // ----
   const RealType refinementFactor = searchRadius / ( 2.f * modelParams.h );
   const unsigned int dim = SPHConfig::spaceDimension;
   const RealType h = refinementFactor * modelParams.h;
   const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
   // ----
   const bool includeBoundary = this->includeBoundary;

   typename EOS::ParamsType eosParams( modelParams );

   auto view_sensorsPositions = sensorPositions.getView();
   auto view_pressureSensors = sensors.getView();

   auto interpolate =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType & r_i, RealType * p, VectorType * v, RealType * gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;

         //*rho += W * m;
         *p += p_j * W * V;
         *gamma += W * V;
      }
   };

   auto interpolateBoundary =
      [ = ] __cuda_callable__(
         LocalIndexType i, LocalIndexType j, VectorType & r_i, RealType * p, VectorType * v, RealType * gamma ) mutable
   {
      const VectorType r_j = view_points_boundary[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho_boundary[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;

         //*rho += W * m;
         *p += p_j * W * V;
         *gamma += W * V;
      }
   };

   auto sensorsLoop = [ = ] __cuda_callable__( LocalIndexType i, GlobalIndexType sensorIndexer ) mutable
   {
      const VectorType r = view_sensorsPositions[ i ];
      RealType p = 0.f;
      VectorType v = 0.f;
      RealType gamma = 0.f;

      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInFluid, interpolate, &p, &v, &gamma );
      if( includeBoundary ) {
         ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInBound, interpolateBoundary, &p, &v, &gamma );
      }

      if( gamma > 0.5f ) {
         view_pressureSensors( sensorIndexer, i ) = p / gamma;
      }
      else {
         view_pressureSensors( sensorIndexer, i ) = 0.f;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, sensorsLoop, this->sensorIndexer );

   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename EOS, typename SourcePointer >
void
SensorInterpolation< SPHConfig, SPHSimulation >::accumulatePressureSource(
   const SourcePointer& source,
   const IndexArrayType& owners,
   ScalarArrayType& pressure,
   ScalarArrayType& gamma,
   const GlobalIndexType levelId,
   const RealType h,
   const RealType m,
   const typename EOS::ParamsType& eosParams )
{
   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const RealType searchRadius = source->getParticles()->getSearchRadius();
   const auto ownersView = owners.getConstView();
   const auto sensorPositionsView = sensorPositions.getConstView();
   auto pressureView = pressure.getView();
   auto gammaView = gamma.getView();

   auto interpolate =
      [ = ] __cuda_callable__( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* p, RealType* g ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;
         *p += p_j * W * V;
         *g += W * V;
      }
   };

   auto sensorsLoop = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( ownersView[ i ] != levelId )
         return;
      VectorType r = sensorPositionsView[ i ];
      RealType p = pressureView[ i ];
      RealType g = gammaView[ i ];
      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInSource, interpolate, &p, &g );
      pressureView[ i ] = p;
      gammaView[ i ] = g;
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfSensors, sensorsLoop );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename EOS, typename ModelParams, typename LevelView >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   ModelParams& modelParams )
{
   if( ! validateMeasurementLevels( levelViews ) )
      throw std::invalid_argument( "Invalid multiresolution measurement levels" );

   IndexArrayType owners;
   auto ownerBoxes = makeMeasurementOwnerBoxArray< SPHConfig >( levelViews );
   computeMeasurementOwners( sensorPositions, ownerBoxes, owners );
   auto ownersView = owners.getView();

   ScalarArrayType pressure;
   ScalarArrayType gamma;
   pressure.setSize( this->numberOfSensors );
   gamma.setSize( this->numberOfSensors );
   pressure = 0.f;
   gamma = 0.f;
   auto pressureView = pressure.getView();
   auto gammaView = gamma.getView();
   auto sensorsView = sensors.getView();
   const GlobalIndexType outputRow = this->sensorIndexer;

   auto resetRow = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      sensorsView( outputRow, i ) = 0.f;
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfSensors, resetRow );

   typename EOS::ParamsType eosParams( modelParams );
   const bool includeBoundary = this->includeBoundary;
   for( std::size_t level = 0; level < levelViews.size(); ++level ) {
      const int levelId = static_cast< int >( level );
      const RealType h = levelViews[ level ].h;
      const RealType m = levelViews[ level ].mass;

      accumulatePressureSource< SPHKernelFunction, EOS >(
         *levelViews[ level ].fluid, owners, pressure, gamma, levelId, h, m, eosParams );
      for( const auto* patch : levelViews[ level ].ownPatches )
         accumulatePressureSource< SPHKernelFunction, EOS >( *patch, owners, pressure, gamma, levelId, h, m, eosParams );
      if( includeBoundary )
         accumulatePressureSource< SPHKernelFunction, EOS >(
            *levelViews[ level ].boundary, owners, pressure, gamma, levelId, h, m, eosParams );

      auto normalize = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
      {
         if( ownersView[ i ] != levelId )
            return;
         if( gammaView[ i ] > 0.5f )
            sensorsView( outputRow, i ) = pressureView[ i ] / gammaView[ i ];
         else
            sensorsView( outputRow, i ) = 0.f;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfSensors, normalize );
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
//template<typename SPHKernelFunction, typename EOS, typename ModelParams >
//void
//SensorGeneralInterpolation< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid,
//                                                                     BoundaryPointer& boundary,
//                                                                     ModelParams& modelParams,
//                                                                     bool includeBoundary )
//{
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
template< typename SPHKernelFunction, typename EOS, typename ModelParams >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::interpolate(
   FluidPointer& fluid,
   BoundaryPointer& boundary,
   ModelParams& modelParams )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->getParticles() );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->getParticles() );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->getParticles()->getPoints().getConstView();
   const auto view_rho = fluid->getVariables()->rho.getConstView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->getParticles()->getSearchRadius();
   const RealType h = modelParams.h;
   const RealType m = modelParams.mass;
   const RealType levelIncrement = this->levelIncrement;

   auto view_sensors = sensors.getView();
   auto view_levels = levels.getView();

   auto interpolateWaterLevel =
      [ = ] __cuda_callable__( LocalIndexType i, LocalIndexType j, VectorType & r_i, RealType * gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;
         *gamma += W * V;
      }
   };

   for( int s = 0; s < numberOfSensors; s++ ) {
      const VectorType direction = this->direction;
      const VectorType sensorLocation = sensorPositions.getElement( s );
      view_levels = 0;

      auto sensorsLoop = [ = ] __cuda_callable__( LocalIndexType i ) mutable
      {
         RealType gamma = 0.f;
         const VectorType r = sensorLocation + i * levelIncrement * direction;
         ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInFluid, interpolateWaterLevel, &gamma );

         if( gamma > 0.5f )
            view_levels[ i ] = 1;
         else
            view_levels[ i ] = 0;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, sensorsLoop );

      auto fetch = [ = ] __cuda_callable__( GlobalIndexType i ) -> GlobalIndexType
      {
         return view_levels[ i ];
      };
      auto reduction = [] __cuda_callable__( const GlobalIndexType& a, const GlobalIndexType& b )
      {
         return a + b;
      };
      const GlobalIndexType numberOfFilledLevels =
         Algorithms::reduce< DeviceType >( 0, view_levels.getSize(), fetch, reduction, 0 );
      //const RealType waterLevel = h * numberOfFilledLevels + this->startLevel; //start level can be written as ( startPoint,
      //direction )
      const RealType waterLevel =
         levelIncrement * numberOfFilledLevels
         + ( sensorLocation,
             direction );  //TODO: in case start level is included, the result is nonsense //TODO: added start water level
      view_sensors( sensorIndexer, s ) = waterLevel;
   }

   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SourcePointer >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::accumulateWaterSource(
   const SourcePointer& source,
   const VectorArrayType& samplePositions,
   const IndexArrayType& sampleOwners,
   ScalarArrayType& sampleGamma,
   const GlobalIndexType levelId,
   const RealType h,
   const RealType m )
{
   typename ParticlesType::NeighborsLoopParams searchInSource( source->getParticles() );
   const auto view_points = source->getParticles()->getPoints().getConstView();
   const auto view_rho = source->getVariables()->rho.getConstView();
   const RealType searchRadius = source->getParticles()->getSearchRadius();
   const auto samplePositionsView = samplePositions.getConstView();
   const auto sampleOwnersView = sampleOwners.getConstView();
   auto sampleGammaView = sampleGamma.getView();

   auto interpolateWaterLevel =
      [ = ] __cuda_callable__( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const RealType drs = l2Norm( r_i - r_j );
      if( drs <= searchRadius ) {
         const RealType rho_j = view_rho[ j ];
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;
         *gamma += W * V;
      }
   };

   auto samplesLoop = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
   {
      if( sampleOwnersView[ i ] != levelId )
         return;
      VectorType r = samplePositionsView[ i ];
      RealType gamma = sampleGammaView[ i ];
      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInSource, interpolateWaterLevel, &gamma );
      sampleGammaView[ i ] = gamma;
   };
   Algorithms::parallelFor< DeviceType >( 0, this->numberOfLevels, samplesLoop );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename EOS, typename ModelParams, typename LevelView >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::interpolate(
   const std::vector< LevelView >& levelViews,
   ModelParams& modelParams )
{
   (void) modelParams;
   if( ! validateMeasurementLevels( levelViews ) )
      throw std::invalid_argument( "Invalid multiresolution measurement levels" );

   auto ownerBoxes = makeMeasurementOwnerBoxArray< SPHConfig >( levelViews );
   VectorArrayType samplePositions;
   IndexArrayType sampleOwners;
   ScalarArrayType sampleGamma;
   samplePositions.setSize( this->numberOfLevels );
   sampleOwners.setSize( this->numberOfLevels );
   sampleGamma.setSize( this->numberOfLevels );
   if( this->levels.getSize() != this->numberOfLevels )
      this->levels.setSize( this->numberOfLevels );

   auto samplePositionsView = samplePositions.getView();
   auto sampleOwnersView = sampleOwners.getView();
   auto sampleGammaView = sampleGamma.getView();
   auto wetLevelsView = this->levels.getView();
   auto sensorsView = sensors.getView();
   const RealType levelIncrement = this->levelIncrement;
   const VectorType direction = this->direction;
   const GlobalIndexType outputRow = this->sensorIndexer;

   for( GlobalIndexType sensor = 0; sensor < this->numberOfSensors; ++sensor ) {
      const VectorType sensorLocation = sensorPositions.getElement( sensor );
      auto setSamplePosition = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
      {
         samplePositionsView[ i ] = sensorLocation + i * levelIncrement * direction;
      };
      Algorithms::parallelFor< DeviceType >( 0, this->numberOfLevels, setSamplePosition );
      computeMeasurementOwners( samplePositions, ownerBoxes, sampleOwners );
      sampleGamma = 0.f;
      this->levels = 0;

      for( std::size_t level = 0; level < levelViews.size(); ++level ) {
         const int levelId = static_cast< int >( level );
         const RealType h = levelViews[ level ].h;
         const RealType m = levelViews[ level ].mass;

         accumulateWaterSource< SPHKernelFunction >(
            *levelViews[ level ].fluid, samplePositions, sampleOwners, sampleGamma, levelId, h, m );
         for( const auto* patch : levelViews[ level ].ownPatches )
            accumulateWaterSource< SPHKernelFunction >(
               *patch, samplePositions, sampleOwners, sampleGamma, levelId, h, m );

         auto classifySample = [ = ] __cuda_callable__( GlobalIndexType i ) mutable
         {
            if( sampleOwnersView[ i ] == levelId )
               wetLevelsView[ i ] = sampleGammaView[ i ] > 0.5f ? 1 : 0;
         };
         Algorithms::parallelFor< DeviceType >( 0, this->numberOfLevels, classifySample );
      }

      auto fetch = [ = ] __cuda_callable__( GlobalIndexType i ) -> GlobalIndexType
      {
         return wetLevelsView[ i ];
      };
      auto reduction = [] __cuda_callable__( const GlobalIndexType& a, const GlobalIndexType& b )
      {
         return a + b;
      };
      const GlobalIndexType numberOfFilledLevels =
         Algorithms::reduce< DeviceType >( 0, wetLevelsView.getSize(), fetch, reduction, 0 );
      const RealType waterLevel = levelIncrement * numberOfFilledLevels + ( sensorLocation, direction );
      sensorsView( outputRow, sensor ) = waterLevel;
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
