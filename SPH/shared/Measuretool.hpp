namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig, typename SPHSimulation >
template< typename SPHKernelFunction, typename SPHState >
void
InterpolateToGrid< SPHConfig, SPHSimulation >::interpolate( FluidPointer& fluid, BoundaryPointer& boundary, SPHState& sphState )
{
   /* PARTICLES AND NEIGHBOR SEARCH ARRAYS */ //TODO: Do this like a human.
   GlobalIndexType numberOfParticles = fluid->particles->getNumberOfParticles();
   GlobalIndexType numberOfParticles_bound = boundary->particles->getNumberOfParticles();
   const RealType searchRadius = fluid->particles->getSearchRadius();

   const VectorType gridOrigin = fluid->particles->getGridOrigin();
   const IndexVectorType gridSize = fluid->particles->getGridSize();

   //typename NeighborSearch::NeighborsLoopParams fluidLoopParams( fluid->neighborSearch );
   //typename NeighborSearch::NeighborsLoopParams boundaryLoopParams( boundary->neighborSearch );

   const auto view_firstLastCellParticle = fluid->neighborSearch->getCellFirstLastParticleList().getView();
   const auto view_particleCellIndex = fluid->particles->getParticleCellIndices().getView();

   /* VARIABLES AND FIELD ARRAYS */
   const auto view_points = fluid->particles->getPoints().getView();
   const auto view_rho = fluid->variables->rho.getView();
   const auto view_v = fluid->variables->v.getView();

   auto view_rho_interpolation = this->variables->rho.getView();
   auto view_v_interpolation = this->variables->v.getView();

   /* CONSTANT VARIABLES */ //TODO: Do this like a human.
   const RealType h = sphState.h;
   const RealType m = sphState.mass;

   const IndexVectorType gridDimension = this->gridDimension;

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

      VectorType r = { ( i[ 0 ] + 1 ) * searchRadius , ( i[ 1 ] + 1 ) * searchRadius };
      const IndexVectorType gridIndex = TNL::floor( ( r - gridOrigin ) / searchRadius );
      const GlobalIndexType idx =  i[ 1 ] * gridDimension[ 0 ] + i[ 0 ];

      //fluidLoopParams.i = i[ 0 ];
      //fluidLoopParams.gridIndex = gridIndex;

      neighborSearch->loopOverNeighbors(
            //fluidLoopParams,
            i[ 0 ],
            numberOfParticles, //TODO: Remove.
            gridIndex,
            gridSize, //TODO: Remove
            view_firstLastCellParticle,
            interpolate, r, &rho, &v, &gamma );

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

template< typename SPHConfig, typename SPHSimulation >
void
InterpolateToGrid< SPHConfig, SPHSimulation >save( const std::string outputFileName )
{
   using Writer = TNL::Meshes::Writers::VTKWriter< GridType >;

   std::cout << std::endl << "MEASURETOOL - SAVING INTERPOLATION :" << std::endl;
   std::cout << "Ouput filename................................ " << outputFileName << " ." << std::endl;
   std::cout << "Grid origin................................... " << interpolationGrid.getOrigin() << " ." << std::endl;
   std::cout << "Grid size..................................... " << gridDimension[ 0 ] * gridDimension[ 1 ] << " ." << std::endl;

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


} // SPH
} // ParticleSystem
} // TNL
