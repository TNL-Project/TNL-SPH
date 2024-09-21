namespace TNL {
namespace SPH {


template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SPHState >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingGrid( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   // neighbor search objects
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->particles );

   // loaded arrays
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();

   // loaded constants
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType searchRadius = fluid->particles->getSearchRadius();

   // interpolation grid constatns
   const IndexVectorType gridDimensions = this->interpolationGrid.getDimensions();
   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   // specify the offset to distinguish cell-based and vertex-based interpolation
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
   const VectorType gridOrigin = this->interpolationGrid.getOrigin() + gridOriginOffsetFactor * gridSpaceSteps;

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
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

   auto processPoint = [=] __cuda_callable__ ( const VectorType r, const GlobalIndexType idx ) mutable
   {
      VectorType v = 0.f;
      RealType rho = 0.f;
      RealType gamma = 0.f;

      ParticlesType::NeighborsLoop::exec( idx, r, searchInFluid, interpolate, &rho, &v, &gamma );

      if( gamma > 0.5f ){
         view_v_interpolation[ idx ] = v / gamma;
         view_rho_interpolation[ idx ] = rho /gamma;
      }
      else{
         view_v_interpolation[ idx ] = 0.f;
         view_rho_interpolation[ idx ] = 0.f;
      }
   };

   if( this->interpolatedGridEntity == 0 )
      this->interpolationGrid.template forAllEntities< 0 >(
         [=] __cuda_callable__ ( const GridVertex& gridEntity ) mutable
         {
            const IndexVectorType coords = gridEntity.getCoordinates();
            const VectorType r = coords * gridSpaceSteps + gridOrigin;
            const GlobalIndexType idx = gridEntity.getIndex();
            processPoint( r, idx );
         } );
   else
      this->interpolationGrid.template forAllEntities< SPHConfig::spaceDimension >(
         [=] __cuda_callable__ ( const GridCell& gridEntity ) mutable
         {
            const IndexVectorType coords = gridEntity.getCoordinates();
            const VectorType r = coords * gridSpaceSteps + gridOrigin;
            const GlobalIndexType idx = gridEntity.getIndex();
            processPoint( r, idx );
         } );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SPHState >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolateUsingParallelFor( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   // neighbor search objects
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->particles );

   // loaded arrays
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();

   // loaded constants
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType searchRadius = fluid->particles->getSearchRadius();

   // specify the offset to distinguish cell-based and vertex-based interpolation
   const VectorType gridOriginOffsetFactor = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 0.5f );
   const IndexVectorType gridDimensionsOffset = ( this->interpolatedGridEntity ) ? ( 0.f ) : ( 1.0f );

   // interpolation grid constatns
   const IndexVectorType gridDimensions = this->interpolationGrid.getDimensions() + gridDimensionsOffset;
   const VectorType gridSpaceSteps = this->interpolationGrid.getSpaceSteps();
   const VectorType gridOrigin = this->interpolationGrid.getOrigin() + gridOriginOffsetFactor * gridSpaceSteps;

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* rho, VectorType* v, RealType* gamma ) mutable
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

   auto loopOverGrid = [=] __cuda_callable__ ( const IndexVectorType& i  ) mutable
   {
      const VectorType r = i * gridSpaceSteps + gridOrigin;
      const GlobalIndexType idx = CellIndexer::EvaluateCellIndex( i, gridDimensions );

      VectorType v = 0.f;
      RealType rho = 0.f;
      RealType gamma = 0.f;

      ParticlesType::NeighborsLoop::exec( idx, r, searchInFluid, interpolate, &rho, &v, &gamma );

      if( gamma > 0.5f ){
         view_v_interpolation[ idx ] = v / gamma;
         view_rho_interpolation[ idx ] = rho /gamma;
      }
      else{
         view_v_interpolation[ idx ] = 0.f;
         view_rho_interpolation[ idx ] = 0.f;
      }
   };
   const IndexVectorType beg = 0;
   Algorithms::parallelFor< DeviceType >( beg, gridDimensions, interpolate );
}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SPHState >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   interpolateUsingGrid< SPHKernelFunction, SPHState >( fluid, boundary, sphState );
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
   if( this->interpolatedGridEntity >  0 ){
      writer.template writeCellData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );
      writer.template writeCellData( buffer, "Velocity", 3 );
   }
   else{
      writer.template writePointData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );
      writer.template writePointData( buffer, "Velocity", 3 );
   }

}

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename EOS, typename SPHState >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid,
                                                              BoundaryPointer& boundary,
                                                              SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->particles );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getConstView();
   const auto view_rho = fluid->variables->rho.getConstView();
   const auto view_points_boundary = boundary->particles->getPoints().getConstView();
   const auto view_rho_boundary = boundary->variables->rho.getConstView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const bool includeBoundary = this->includeBoundary;

   typename EOS::ParamsType eosParams( sphState );

   auto view_sensorsPositions = sensorPositions.getView();
   auto view_pressureSensors = sensors.getView();

   auto interpolate = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* p, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;

         //*rho += W * m;
         *p += p_j * W * V;
         *gamma += W * V;
      }
   };

   auto interpolateBoundary = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j,
         VectorType& r_i, RealType* p, VectorType* v, RealType* gamma ) mutable
   {
      const VectorType r_j = view_points_boundary[ j ];
      const VectorType r_ij = r_i - r_j;
      const RealType drs = l2Norm( r_ij );
      if( drs <= searchRadius )
      {
         const RealType rho_j = view_rho_boundary[ j ];
         const RealType p_j = EOS::DensityToPressure( rho_j, eosParams );
         const RealType W = SPHKernelFunction::W( drs, h );
         const RealType V = m / rho_j;

         //*rho += W * m;
         *p += p_j * W * V;
         *gamma += W * V;
      }
   };

   auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, GlobalIndexType sensorIndexer ) mutable
   {
      const VectorType r = view_sensorsPositions[ i ];
      RealType p = 0.f;
      VectorType v = 0.f;
      RealType gamma = 0.f;

      ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInFluid, interpolate, &p, &v, &gamma );
      if( includeBoundary ){
         ParticlesType::NeighborsLoopAnotherSet::exec( i, r, searchInBound, interpolateBoundary, &p, &v, &gamma );
      }

      if( gamma > 0.5f ){
         view_pressureSensors( sensorIndexer, i ) = p / gamma;
      }
      else{
         view_pressureSensors( sensorIndexer, i ) = 0.f;
      }
   };
   Algorithms::parallelFor< DeviceType >( 0, numberOfSensors, sensorsLoop, this->sensorIndexer );

   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorInterpolation< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   using HostSensorsDataArray = Containers::NDArray< RealType,
                                                     Containers::SizesHolder< int, 0, 0 >,
                                                     std::index_sequence< 0, 1 >,
                                                     Devices::Host >;
   HostSensorsDataArray sensorsDataHost;
   sensorsDataHost = sensors;
   const auto view_sensorsDataHost = sensorsDataHost.getView();

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

//template< typename SPHConfig, typename SPHSimulation >
//template<typename SPHKernelFunction, typename EOS, typename SPHState >
//void
//SensorGeneralInterpolation< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid,
//                                                                     BoundaryPointer& boundary,
//                                                                     SPHState& sphState,
//                                                                     bool includeBoundary )
//{
//   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
//   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
//   const RealType searchRadius = fluid->particles->getSearchRadius();
//
//   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
//   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );
//
//   /* VARIABLES AND FIELD ARRAYS */
//   const auto view_points = fluid->particles->getPoints().getView();
//   const auto view_rho = fluid->variables->rho.getView();
//
//   const auto view_points_boundary = boundary->particles->getPoints().getView();
//   const auto view_rho_boundary = boundary->variables->rho.getView();
//
//   const auto view_field = field.getView();
//
//   /* CONSTANT VARIABLES */
//   const RealType h = sphState.h;
//   const RealType m = sphState.mass;
//
//   typename EOS::ParamsType eosParams( sphState );
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
template< typename SPHKernelFunction, typename EOS, typename SPHState >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticlesType::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticlesType::NeighborsLoopParams searchInBound( boundary->particles );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getConstView();
   const auto view_rho = fluid->variables->rho.getConstView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType levelIncrement = this->levelIncrement;

   auto view_sensors = sensors.getView();
   auto view_levels = levels.getView();

   auto interpolateWaterLevel = [=] __cuda_callable__ ( LocalIndexType i, LocalIndexType j, VectorType& r_i, RealType* gamma ) mutable
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

   for( int s = 0; s < numberOfSensors; s++ )
   {
      const VectorType direction = this->direction;
      const VectorType sensorLocation = sensorPositions.getElement( s );
      view_levels = 0;

      auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i ) mutable
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

      auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_levels[ i ]; };
      auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      const GlobalIndexType numberOfFilledLevels = Algorithms::reduce< DeviceType >( 0, view_levels.getSize(), fetch, reduction, 0 );
      //const RealType waterLevel = h * numberOfFilledLevels + this->startLevel; //start level can be written as ( startPoint, direction )
      const RealType waterLevel = levelIncrement * numberOfFilledLevels; //TODO: in case start level is included, the result is nonsense
      view_sensors( sensorIndexer, s ) = waterLevel;
   }

   sensorIndexer++;
}

template< typename SPHConfig, typename SPHSimulation >
void
SensorWaterLevel< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   using HostSensorsDataArray = Containers::NDArray< RealType,
                                                     Containers::SizesHolder< int, 0, 0 >,
                                                     std::index_sequence< 0, 1 >,
                                                     Devices::Host >;
   HostSensorsDataArray sensorsDataHost;
   sensorsDataHost = sensors;
   const auto view_sensorsDataHost = sensorsDataHost.getView();

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


} // SPH
} // TNL
