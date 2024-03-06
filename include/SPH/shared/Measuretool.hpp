namespace TNL {
namespace SPH {

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SPHState >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();

   /* CONSTANT VARIABLES */
   const RealType h = sphState.h;
   const RealType m = sphState.mass;
   const RealType searchRadius = fluid->particles->getSearchRadius();

   const IndexVectorType _gridDimension = this->gridDimension;
   const VectorType _gridStep = this->gridStep;
   const VectorType _gridOrigin = this->interpolationGrid.getOrigin();

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

   if constexpr( SPHState::SPHConfig::spaceDimension == 2 ){
      auto gridLoop = [=] __cuda_callable__ ( const IndexVectorType& i ) mutable
      {
         VectorType v = 0.f;
         RealType rho = 0.f;
         RealType gamma = 0.f;

         //TODO: Replace following lines with some reasonable function
         //NOTE: I used to increase the indices i[ 0 ], i[ 1 ] by + 1 but I dont know why.
         const VectorType  r = { ( i[ 0 ] ) * _gridStep[ 0 ] + _gridOrigin[ 0 ],
                                 ( i[ 1 ] ) * _gridStep[ 1 ] + _gridOrigin[ 1 ]  };
         const GlobalIndexType  idx =  i[ 1 ] * _gridDimension[ 0 ] + i[ 0 ];

         TNL::ParticleSystem::NeighborsLoop::exec( i[ 0 ], r, searchInFluid, interpolate, &rho, &v, &gamma );

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
      Algorithms::parallelFor< DeviceType >( begin, gridDimension, gridLoop );
   }
   if constexpr( SPHState::SPHConfig::spaceDimension == 3 ){

      auto gridLoop = [=] __cuda_callable__ ( const IndexVectorType& i ) mutable
      {
         VectorType v = 0.f;
         RealType rho = 0.f;
         RealType gamma = 0.f;

         //TODO: Replace following lines with some reasonable function
         //NOTE: I used to increase the indices i[ 0 ], i[ 1 ], i[ 2 ]  by + 1 but I dont know why.
         const VectorType r = { ( i[ 0 ] + 1 ) * _gridStep[ 0 ] + _gridOrigin[ 0 ],
                                ( i[ 1 ] + 1 ) * _gridStep[ 1 ] + _gridOrigin[ 1 ],
                                ( i[ 2 ] + 1 ) * _gridStep[ 2 ] + _gridOrigin[ 2 ]};
         const GlobalIndexType  idx =  i[ 2 ] * _gridDimension[ 0 ] * _gridDimension[ 1 ] + \
                                       i[ 1 ] * _gridDimension[ 0 ] + \
                                       i[ 0 ];

         TNL::ParticleSystem::NeighborsLoop::exec( i[ 0 ], r, searchInFluid, interpolate, &rho, &v, &gamma );

         if( gamma > 0.5f ){
            view_v_interpolation[ idx ] = v / gamma;
            view_rho_interpolation[ idx ] = rho /gamma;
         }
         else{
            view_v_interpolation[ idx ] = 0.f;
            view_rho_interpolation[ idx ] = 0.f;
         }
      };
      IndexVectorType begin{ 0, 0, 0 };
      Algorithms::parallelFor< DeviceType >( begin, gridDimension, gridLoop );
   }

}

template< typename SPHConfig, typename SPHSimulation >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::save( const std::string outputFileName )
{
   using Writer = TNL::Meshes::Writers::VTKWriter< GridType >;

   std::ofstream file( outputFileName );
   Writer writer( file );
   writer.writeEntities( interpolationGrid );

   writer.template writeCellData< typename Variables::ScalarArrayType >( variables->rho, "Density", 1 );

   //FIXME: Workaround for vectorArray;
   using BufferType = Containers::Array< RealType, Devices::Host, GlobalIndexType >;
   BufferType buffer( 3 * gridDimension[ 0 ] * gridDimension[ 1 ] );

   const auto velocityView = variables->v.getView();
   GlobalIndexType k = 0;
   GlobalIndexType entitiesCount = gridDimension[ 0 ] * gridDimension[ 1 ];
   for( GlobalIndexType i = 0; i < entitiesCount; i++ ) {
      const VectorType vector = velocityView.getElement( i );
      for( int j = 0; j < 3; j++ )
         buffer[ k++ ] = ( j < vector.getSize() ? vector[ j ] : 0 );
   }

   writer.template writeCellData( buffer, "Velocity", 3 );
}

template< typename SPHConfig, typename SPHSimulation >
template<typename SPHKernelFunction, typename EOS, typename SPHState >
void
SensorInterpolation< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid,
                                                              BoundaryPointer& boundary,
                                                              SPHState& sphState,
                                                              bool includeBoundary )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_points_boundary = boundary->particles->getPoints().getView();
   const auto view_rho_boundary = boundary->variables->rho.getView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
   const RealType h = sphState.h;
   const RealType m = sphState.mass;

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

      TNL::ParticleSystem::NeighborsLoop::exec( i, r, searchInFluid, interpolate, &p, &v, &gamma );
      if( includeBoundary ){
         TNL::ParticleSystem::NeighborsLoop::exec( i, r, searchInBound, interpolateBoundary, &p, &v, &gamma );
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
   typename ParticleSystem::NeighborsLoopParams searchInFluid( fluid->particles );
   typename ParticleSystem::NeighborsLoopParams searchInBound( boundary->particles );

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();

   /* CONSTANT VARIABLES */
   const RealType searchRadius = fluid->particles->getSearchRadius();
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

   for( int s = 0; s < numberOfSensors; s++ )
   {
      const VectorType direction = this->direction;
      view_levels = 0.f;

      auto sensorsLoop = [=] __cuda_callable__ ( LocalIndexType i, GlobalIndexType sensorIndexer ) mutable
      {
         RealType gamma = 0.f;
         const VectorType r = view_sensorsPositions[ s ] + i * h * direction;
         TNL::ParticleSystem::NeighborsLoop::exec( i, r, searchInFluid, interpolate, &gamma );

         if( gamma > 0.5f )
            view_levels[ i ] = 1;
         else
            view_levels[ i ] = 0;
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfLevels, sensorsLoop, this->sensorIndexer );

      auto fetch = [=] __cuda_callable__ ( GlobalIndexType i ) -> GlobalIndexType { return view_levels[ i ]; };
      auto reduction = [] __cuda_callable__ ( const GlobalIndexType& a, const GlobalIndexType& b ) { return a + b; };
      const GlobalIndexType numberOfFilledLevels = Algorithms::reduce< DeviceType >( 0, view_levels.getSize(), fetch, reduction, 0.0 );
      //const RealType waterLevel = h * numberOfFilledLevels + this->startLevel; //start level can be written as ( startPoint, direction )
      const RealType waterLevel = h * numberOfFilledLevels; //TODO: in case start level is included, the result is nonsense
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
