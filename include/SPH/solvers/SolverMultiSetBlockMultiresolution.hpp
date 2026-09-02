#include "SolverMultiSetBlockMultiresolution.h"
#include "../tempFunctionsToConfigMultiresolution.h"

namespace TNL {
namespace SPH {

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::init( int argc, char* argv[] )
{
   auto& params = this->parameters;
   auto& log = this->logger;

   try {
      initialize< SimulationType >( argc, argv, this->cliParams, this->cliConfig, params, this->config );
   }
   catch ( ... ) {
      std::cerr << std::endl;
   }

   log.writeHeader( "SPH simulation initialization." );
   this->caseName = params.template getParameter< std::string >( "case-name" );
   this->verbose = params.template getParameter< std::string >( "verbose-intensity" );
   this->outputDirectory = params.template getParameter< std::string >( "output-directory" );
   this->particlesFormat = params.template getParameter< std::string >( "particles-format" );

#ifdef HAVE_MPI
   initializeDistributedSimulation();
#else
   initializeBlockBasedMultiResolutionSimulation();
#endif
   log.writeHeader( "SPH simulation successfully initialized." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initializeBlockBasedMultiResolutionSimulation()
{
   auto& params = this->parameters;
   auto& log = this->logger;

   this->numberOfSubsets = params.template getParameter< int >( "numberOfSubdomains" );
   const std::string configSubdomainsPath = params.template getParameter< std::string >( "subdomains-config" );
   const std::string subdomainsInline = params.template getParameter< std::string >( "subdomains" );
    for( int subset = 0; subset < this->numberOfSubsets; subset++ )
       TNL::SPH::configSubdomain( subset, this->configSubdomains );
    for( int bufferIdx = 0; bufferIdx < 2 * this->numberOfSubsets; bufferIdx++ ) {
       TNL::SPH::configMultiresolutionBuffer( bufferIdx, this->configSubdomains );
       TNL::SPH::configBoundaryGhostBuffer( bufferIdx, this->configSubdomains );
    }
   parseDistributedConfig( configSubdomainsPath, parametersSubdomains, configSubdomains, log, subdomainsInline );

   topology.loadFromConfig( params, parametersSubdomains );
   topology.finalizeLinear();
   initParticleSets();
   initMultiResolutionBoundaryPatches();
   this->timeMeasurement.addTimer( "multiresolution-update" );

   const bool hasOpenBcFile = params.template getParameter< std::string >( "open-boundary-config" ) != "";
   const bool hasOpenBcInline = params.template getParameter< std::string >( "open-boundary" ) != "";
   if( hasOpenBcFile || hasOpenBcInline ){
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

    const bool hasPeriodicBcFile = params.template getParameter< std::string >( "periodic-boundary-config" ) != "";
    const bool hasPeriodicBcInline = params.template getParameter< std::string >( "periodic-boundary" ) != "";
    if( hasPeriodicBcFile || hasPeriodicBcInline ){
       this->initPeriodicBoundaryPatches( params, log );

       this->timeMeasurement.addTimer( "enforce-periodic-bc" );
       this->timeMeasurement.addTimer( "transfer-periodic-bc" );
       this->timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
       this->timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
    }

    this->modelParams.init( params );
    for( long unsigned int p = 0; p < multiresolutionBoundaryPatches.size(); p++ ){
      const auto& iface = multiresolutionBoundaryPatchInterfaces[ p ];
      std::string subdomainKey = "subdomain-" + std::to_string( iface.ownIdx ) + "-";
      const float refinementFactor = parametersSubdomains.getParameter< float >( subdomainKey + "refinement-factor" );
      multiresolutionBoundaryPatches[ p ]->initMassNodes( this->modelParams, iface.ownIdx, iface.neighborIdx, refinementFactor );
   }

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initial-time-step" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

    readParticlesFiles();
    initBoundaryGhosts();

    log.writeSeparator();
    const bool hasMeasuretoolFile = params.template getParameter< std::string >( "measuretool-config" ) != "";
    const bool hasMeasuretoolInline = params.template getParameter< std::string >( "measuretool" ) != "";
    if( hasMeasuretoolFile || hasMeasuretoolInline ) {
       log.writeParameter( "Simulation monitor initialization.", "" );
       this->simulationMonitor.init( params, this->timeStepping, log );
       this->simulationMonitor.setupVolumetricFlowRateZones( this->fluidSets[ 0 ] );
       log.writeParameter( "Simulation monitor initialization.", "Done." );
    }
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initializeDistributedSimulation()
{
   auto& params = this->parameters;
   auto& log = this->logger;

   this->numberOfSubsets = 1;
   this->fluidSets.resize( 1 );
   this->boundarySets.resize( 1 );

   log.writeParameter( "Configuration of distributed simulation:", "" );
   Containers::StaticVector< 2, int > numberOfSubdomains = params.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   log.writeParameter( "Number of subdomains:", numberOfSubdomains );
   for( int x = 0; x < numberOfSubdomains[ 0 ]; x++ )
      for( int y = 0; y < numberOfSubdomains[ 1 ]; y++ )
         TNL::SPH::configSetupDistributedSubdomain( x, y, this->configDistributed );
    std::string configDistributedPath = params.template getParameter< std::string >( "distributed-config" );
    const std::string distributedDomainInline = params.template getParameter< std::string >( "distributed-domain" );
    parseDistributedConfig( configDistributedPath, this->parametersDistributed, this->configDistributed, log, distributedDomainInline );

   initDistributedParticleSets( params, this->parametersDistributed, log );

   this->loadBalancingMeasure = params.template getParameter< std::string >( "load-balancing-measure" );
   this->loadBalancingStepInterval = params.template getParameter< int >( "load-balancing-step-inteval" );
   this->fluid()->getDistributedParticles()->setParticlesCountResizeThreshold(
         params.template getParameter< float >( "number-of-particles-balancing-coef" ) );
   this->fluid()->getDistributedParticles()->setCompTimeResizePercentageThreshold(
         params.template getParameter< float >( "computational-time-balancing-coef" ) );

   this->timeMeasurement.addTimer( "synchronize", false );
   this->timeMeasurement.addTimer( "rebalance", false );

   const bool hasOpenBcFile = params.template getParameter< std::string >( "open-boundary-config" ) != "";
   const bool hasOpenBcInline = params.template getParameter< std::string >( "open-boundary" ) != "";
   if( hasOpenBcFile || hasOpenBcInline ){
      this->initOpenBoundaryPatches( params, log );

      this->timeMeasurement.addTimer( "extrapolate-openbc" );
      this->timeMeasurement.addTimer( "apply-openbc" );
   }

    const bool hasPeriodicBcFile = params.template getParameter< std::string >( "periodic-boundary-config" ) != "";
    const bool hasPeriodicBcInline = params.template getParameter< std::string >( "periodic-boundary" ) != "";
    if( hasPeriodicBcFile || hasPeriodicBcInline ){
       this->initPeriodicBoundaryPatches( params, log );

       this->timeMeasurement.addTimer( "enforce-periodic-bc" );
       this->timeMeasurement.addTimer( "transfer-periodic-bc" );
       this->timeMeasurement.addTimer( "periodicity-fluid-updateZone", false );
       this->timeMeasurement.addTimer( "periodicity-boundary-updateZone", false );
    }

   this->modelParams.init( params );

   this->timeStepping.setTimeStep( params.template getParameter< RealType >( "initial-time-step" ) );
   this->timeStepping.setEndTime( params.template getParameter< RealType >( "final-time" ) );
   this->timeStepping.addOutputTimer( "save_results", params.template getParameter< RealType >( "snapshot-period" ) );

   readParticleFilesDistributed( params, this->parametersDistributed, log );

   log.writeSeparator();
    if( params.template getParameter< std::string >( "measuretool-config" ) != "" ) {
       log.writeParameter( "Simulation monitor initialization.", "" );
       this->simulationMonitor.init( params, this->timeStepping, log );
       this->simulationMonitor.setupVolumetricFlowRateZones( this->fluidSets[ 0 ] );
       log.writeParameter( "Simulation monitor initialization.", "Done." );
    }
}
#endif

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initParticleSets()
{
   const int nSubsets = topology.getNumberOfSubdomains();
   this->numberOfSubsets = nSubsets;
   this->fluidSets.resize( nSubsets );
   this->boundarySets.resize( nSubsets );
   multiresolutionBoundaryPatches.resize( nSubsets );

   for( int i = 0; i < nSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      this->fluidSets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->fluidSets[ i ]->getParticles()->setParticleSetLabel( 0 );

      this->boundarySets[ i ]->initializeAsDistributed(
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n" ),
            parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );
      if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
         this->boundarySets[ i ]->getParticles()->setParticleSetLabel( 1 );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initMultiResolutionBoundaryPatches()
{
   const int mrbCount = topology.getNumberOfSubdomainInterfaces();
   multiresolutionBoundaryPatches.resize( mrbCount );

   int mrbIdx = 0;
   for( int i = 0; i < topology.getNumberOfSubdomains(); i++ ) {
      const std::string key = "subdomain-" + std::to_string( i ) + "-";
      const float rf = parametersSubdomains.getParameter< float >( key + "refinement-factor" );

      for( const auto& iface : topology.getInterfacesOfSubdomain( i ) ) {

         const std::string bufferKey = "multiresolution-buffer-" + std::to_string( mrbIdx ) + "-";
         const int mrbNumOfPtcs = parametersSubdomains.getParameter< int >( bufferKey + "n" );
         // TODO: use some reasonable estimate
         const RealType resolutionScale = std::pow( rf, -int( SPHConfig::spaceDimension ) );
         const int mrbNumOfAllocPtcs =
            std::max( int( 0.1 * this->getTotalFluidParticlesCount() * resolutionScale ), 2 * mrbNumOfPtcs );

         multiresolutionBoundaryPatches[ mrbIdx ]->initializeAsDistributed(
            mrbNumOfPtcs,
            mrbNumOfAllocPtcs,
            topology.getLocalGrid( i ),
            topology.getLocalOriginCoordinates( i ),
            topology.getGlobalGrid(),
            topology.getNumberOfOverlapLayers() );

         multiresolutionBoundaryPatches[ mrbIdx ]->initZones(
            this->fluidSets[ iface.ownIdx ]->getParticles(),
            this->fluidSets[ iface.neighborIdx ]->getParticles(),
            rf );

         multiresolutionBoundaryPatchInterfaces.push_back( iface );
         mrbIdx++;
      }
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::readParticlesFiles()
{
   auto& params = this->parameters;
   auto& log = this->logger;
   auto& paramsOB = this->parametersOpenBoundary;

   const int nSubsets = params.template getParameter< int >( "numberOfSubdomains" );

   for( int i = 0; i < nSubsets; i++ ){
      std::string subdomainKey = "subdomain-" + std::to_string( i ) + "-";
      if( parametersSubdomains.getParameter< int >( subdomainKey + "fluid_n" ) != 0 ){
         const std::string fluidFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "fluid-particles" );
         log.writeParameter( "Reading fluid particles:", fluidFileName );
         this->fluidSets[ i ]->template readParticlesAndVariables< typename BaseType::Reader >( fluidFileName );
      }
      if( parametersSubdomains.getParameter< int >( subdomainKey + "boundary_n" ) != 0 ){
         const std::string boundaryFileName = parametersSubdomains.getParameter< std::string >( subdomainKey + "boundary-particles" );
         log.writeParameter( "Reading boundary particles:", boundaryFileName );
         this->boundarySets[ i ]->template readParticlesAndVariables< typename BaseType::Reader >( boundaryFileName );
      }
   }

   const int numberOfBoundaryPatches = params.template getParameter< int >( "openBoundaryPatches" );
   if( numberOfBoundaryPatches > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         log.writeParameter( "Reading open boundary particles:", paramsOB.template getParameter< std::string >( prefix + "particles" ) );
         this->openBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::Reader >(
            paramsOB.template getParameter< std::string >( prefix + "particles" ) );
      }
   }

   const int numberOfMultiresolutionBuffers = multiresolutionBoundaryPatches.size();
   for( int i = 0; i < numberOfMultiresolutionBuffers; i++ ) {
      std::string bufferKey = "multiresolution-buffer-" + std::to_string( i ) + "-";
      if( parametersSubdomains.getParameter< int >( bufferKey + "n" ) != 0 ) {
         const std::string mrbFileName = parametersSubdomains.getParameter< std::string >( bufferKey + "particles" );
         log.writeParameter( "Reading multiresolution buffer particles:", mrbFileName );
         multiresolutionBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::Reader >( mrbFileName );
      }
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initBoundaryGhosts()
{
   auto& log = this->logger;

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchInterfaces.size(); p++ ) {
      const auto& iface = multiresolutionBoundaryPatchInterfaces[ p ];
      const auto& patch = multiresolutionBoundaryPatches[ p ];
      BoundaryPointer& ownBoundary = this->boundarySets[ iface.ownIdx ];
      BoundaryPointer& srcBoundary = this->boundarySets[ iface.neighborIdx ];

      // own-resolution ghost lattice pre-generated by the init script
      const std::string ghostKey = "boundary-ghost-buffer-" + std::to_string( p ) + "-";
      const int ghostsInFile = parametersSubdomains.getParameter< int >( ghostKey + "n" );
      if( ghostsInFile > 0 ) {
         BoundaryGhostParticles rec;
         rec.ownIdx = iface.ownIdx;
         rec.neighborIdx = iface.neighborIdx;
         rec.ghostBegin = ownBoundary->getNumberOfParticles();
         rec.useInterpolation = true;

         BoundaryPointer ghostGeometry;
         ghostGeometry->initializeAsDistributed( ghostsInFile,
                                                 ghostsInFile,
                                                 topology.getLocalGrid( iface.ownIdx ),
                                                 topology.getLocalOriginCoordinates( iface.ownIdx ),
                                                 topology.getGlobalGrid(),
                                                 topology.getNumberOfOverlapLayers() );
         if constexpr( ParticlesType::specifySearchedSetExplicitly() == true )
            ghostGeometry->getParticles()->setParticleSetLabel( 1 );

         const std::string ghostFileName = parametersSubdomains.getParameter< std::string >( ghostKey + "particles" );
         log.writeParameter( "Reading boundary ghost particles:", ghostFileName );
         ghostGeometry->template readParticlesAndVariables< typename BaseType::Reader >( ghostFileName );

         const GlobalIndexType numberOfGhosts = ghostGeometry->getNumberOfParticles();
         rec.numberOfGhostParticles = numberOfGhosts;

         const GlobalIndexType allocated = ownBoundary->getNumberOfAllocatedParticles();
         if( rec.ghostBegin + numberOfGhosts > allocated )
            throw std::runtime_error( "SolverMultiSetBlockMultiresolution: boundary set " + std::to_string( iface.ownIdx )
                                      + " has no room for " + std::to_string( numberOfGhosts )
                                      + " ghost particles (allocated " + std::to_string( allocated )
                                      + ", needed " + std::to_string( rec.ghostBegin + numberOfGhosts )
                                      + ") - increase the boundary allocation factor in the init script." );

         const GlobalIndexType ghostBegin = rec.ghostBegin;

         auto dst_r = ownBoundary->getPoints().getView();
         const auto src_r = ghostGeometry->getPoints().getConstView();
         auto dst_rho = ownBoundary->getVariables()->rho.getView();
         const auto src_rho = ghostGeometry->getVariables()->rho.getConstView();
         auto dst_drho = ownBoundary->getVariables()->drho.getView();
         auto dst_p = ownBoundary->getVariables()->p.getView();
         auto dst_v = ownBoundary->getVariables()->v.getView();
         const auto src_v = ghostGeometry->getVariables()->v.getConstView();
         auto dst_a = ownBoundary->getVariables()->a.getView();
         auto dst_gamma = ownBoundary->getVariables()->gamma.getView();
         auto dst_marker = ownBoundary->getVariables()->marker.getView();
         const auto src_marker = ghostGeometry->getVariables()->marker.getConstView();
         auto dst_n = ownBoundary->getVariables()->n.getView();
         const auto src_n = ghostGeometry->getVariables()->n.getConstView();
         auto dst_elementSize = ownBoundary->getVariables()->elementSize.getView();
         const auto src_elementSize = ghostGeometry->getVariables()->elementSize.getConstView();

         // the VTK input carries no values for drho, p, a and gamma - start them at zero
         auto copyGhost = [ = ] __cuda_callable__( GlobalIndexType k ) mutable
         {
            const GlobalIndexType d = ghostBegin + k;
            dst_r[ d ] = src_r[ k ];
            dst_rho[ d ] = src_rho[ k ];
            dst_drho[ d ] = 0.f;
            dst_p[ d ] = 0.f;
            dst_v[ d ] = src_v[ k ];
            dst_a[ d ] = 0.f;
            dst_gamma[ d ] = 0.f;
            dst_marker[ d ] = src_marker[ k ];
            dst_n[ d ] = src_n[ k ];
            dst_elementSize[ d ] = src_elementSize[ k ];
         };
         Algorithms::parallelFor< DeviceType >( 0, numberOfGhosts, copyGhost );

         ownBoundary->getParticles()->setNumberOfParticles( ghostBegin + numberOfGhosts );

         boundaryGhostInterfaces.push_back( rec );
         log.writeParameter( "Boundary ghosts for set " + std::to_string( iface.ownIdx )
                                + " from set " + std::to_string( iface.neighborIdx ) + ":",
                             numberOfGhosts );
         continue;
      }

      const CoordinatesType frameBackOrigin = patch->getFrameBackOriginGlobalCoordinates();
      const CoordinatesType frameBackDims = patch->getFrameBackDimensions();
      const CoordinatesType frameFrontOrigin = patch->getFrameFrontOriginGlobalCoordinates();
      const CoordinatesType frameFrontDims = patch->getFrameFrontDimensions();
      const VectorType refOrig = ownBoundary->getParticles()->getReferentialOrigin();
      const RealType inv_sr = 1.f / ownBoundary->getParticles()->getSearchRadius();

      // host-side band selection (init-time only, keeps it deterministic)
      const GlobalIndexType srcCount = srcBoundary->getNumberOfParticles();
      Containers::Array< VectorType, Devices::Host > srcPointsHost;
      srcPointsHost = srcBoundary->getPoints();
      std::vector< GlobalIndexType > srcListHost;
      for( GlobalIndexType j = 0; j < srcCount; j++ ) {
         const CoordinatesType gc = TNL::floor( ( srcPointsHost[ j ] - refOrig ) * inv_sr );
         // band = symmetric difference frameBack ⊖ frameFront, works for both inner and outer frames
         const bool inBack = MultiresolutionBoundary::isInsideBox( gc - frameBackOrigin, frameBackDims );
         const bool inFront = MultiresolutionBoundary::isInsideBox( gc - frameFrontOrigin, frameFrontDims );
         if( inBack != inFront )
            srcListHost.push_back( j );
      }

      const GlobalIndexType numberOfGhosts = srcListHost.size();
      if( numberOfGhosts == 0 )
         continue;

      BoundaryGhostParticles rec;
      rec.ownIdx = iface.ownIdx;
      rec.neighborIdx = iface.neighborIdx;
      rec.ghostBegin = ownBoundary->getNumberOfParticles();
      rec.numberOfGhostParticles = numberOfGhosts;

      const GlobalIndexType allocated = ownBoundary->getNumberOfAllocatedParticles();
      if( rec.ghostBegin + numberOfGhosts > allocated )
         throw std::runtime_error( "SolverMultiSetBlockMultiresolution: boundary set " + std::to_string( iface.ownIdx )
                                   + " has no room for " + std::to_string( numberOfGhosts )
                                   + " ghost particles (allocated " + std::to_string( allocated )
                                   + ", needed " + std::to_string( rec.ghostBegin + numberOfGhosts )
                                   + ") - increase the boundary allocation factor in the init script." );

      Containers::Array< GlobalIndexType, Devices::Host > srcListHostArray( numberOfGhosts );
      for( GlobalIndexType k = 0; k < numberOfGhosts; k++ )
         srcListHostArray[ k ] = srcListHost[ k ];
      rec.srcIndexList = srcListHostArray;

      const GlobalIndexType ghostBegin = rec.ghostBegin;
      const auto srcIdx = rec.srcIndexList.getConstView();

      auto dst_r = ownBoundary->getPoints().getView();
      const auto src_r = srcBoundary->getPoints().getConstView();
      auto dst_rho = ownBoundary->getVariables()->rho.getView();
      const auto src_rho = srcBoundary->getVariables()->rho.getConstView();
      auto dst_drho = ownBoundary->getVariables()->drho.getView();
      const auto src_drho = srcBoundary->getVariables()->drho.getConstView();
      auto dst_p = ownBoundary->getVariables()->p.getView();
      const auto src_p = srcBoundary->getVariables()->p.getConstView();
      auto dst_v = ownBoundary->getVariables()->v.getView();
      const auto src_v = srcBoundary->getVariables()->v.getConstView();
      auto dst_a = ownBoundary->getVariables()->a.getView();
      const auto src_a = srcBoundary->getVariables()->a.getConstView();
      auto dst_gamma = ownBoundary->getVariables()->gamma.getView();
      const auto src_gamma = srcBoundary->getVariables()->gamma.getConstView();
      auto dst_marker = ownBoundary->getVariables()->marker.getView();
      const auto src_marker = srcBoundary->getVariables()->marker.getConstView();
      auto dst_n = ownBoundary->getVariables()->n.getView();
      const auto src_n = srcBoundary->getVariables()->n.getConstView();
      auto dst_elementSize = ownBoundary->getVariables()->elementSize.getView();
      const auto src_elementSize = srcBoundary->getVariables()->elementSize.getConstView();

      auto copyGhost = [ = ] __cuda_callable__( GlobalIndexType k ) mutable
      {
         const GlobalIndexType s = srcIdx[ k ];
         const GlobalIndexType d = ghostBegin + k;
         dst_r[ d ] = src_r[ s ];
         dst_rho[ d ] = src_rho[ s ];
         dst_drho[ d ] = src_drho[ s ];
         dst_p[ d ] = src_p[ s ];
         dst_v[ d ] = src_v[ s ];
         dst_a[ d ] = src_a[ s ];
         dst_gamma[ d ] = src_gamma[ s ];
         dst_marker[ d ] = src_marker[ s ];
         dst_n[ d ] = src_n[ s ];
         dst_elementSize[ d ] = src_elementSize[ s ];
      };
      Algorithms::parallelFor< DeviceType >( 0, numberOfGhosts, copyGhost );

      ownBoundary->getParticles()->setNumberOfParticles( ghostBegin + numberOfGhosts );

      boundaryGhostInterfaces.push_back( rec );
      log.writeParameter( "Boundary ghosts for set " + std::to_string( iface.ownIdx )
                             + " from set " + std::to_string( iface.neighborIdx ) + ":",
                          numberOfGhosts );
   }

   // scratch for the referential-index inverse used by the per-step refresh
   GlobalIndexType scratchSize = 0;
   for( const auto& rec : boundaryGhostInterfaces ) {
      scratchSize = std::max( scratchSize, this->boundarySets[ rec.ownIdx ]->getNumberOfParticles() );
      scratchSize = std::max( scratchSize, this->boundarySets[ rec.neighborIdx ]->getNumberOfParticles() );
   }
   ghostIndexInverseOwn.setSize( scratchSize );
   ghostIndexInverseSrc.setSize( scratchSize );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::refreshBoundaryGhostValues()
{
   for( const auto& rec : boundaryGhostInterfaces ) {
      BoundaryPointer& ownBoundary = this->boundarySets[ rec.ownIdx ];
      BoundaryPointer& srcBoundary = this->boundarySets[ rec.neighborIdx ];

      const GlobalIndexType ownCount = ownBoundary->getNumberOfParticles();

      const auto ownRefIdx = ownBoundary->getVariables()->referentialIdx.getConstView();
      auto invOwn = ghostIndexInverseOwn.getView();
      auto invSrc = ghostIndexInverseSrc.getView();

      Algorithms::parallelFor< DeviceType >( 0, ownCount, [ = ] __cuda_callable__( GlobalIndexType j ) mutable
      {
         invOwn[ ownRefIdx[ j ] ] = j;
      } );

      const GlobalIndexType ghostBegin = rec.ghostBegin;
      auto dst_rho = ownBoundary->getVariables()->rho.getView();
      auto dst_gamma = ownBoundary->getVariables()->gamma.getView();
      const auto src_rho = srcBoundary->getVariables()->rho.getConstView();
      const auto src_gamma = srcBoundary->getVariables()->gamma.getConstView();

      if( rec.useInterpolation ) {
         using KernelFunction = typename ModelParams::KernelFunction;
         using MFD = Interpolation::
            MFD< ParticlesType::getParticlesDimension(), 1, RealType, Interpolation::WendlandC2ABFs, KernelFunction, SPHConfig >;
         using MfdMatrixType = typename MFD::BaseMatrixType;
         using MfdVectorType = typename MFD::BaseVectorType;

         const RealType searchRadius = srcBoundary->getParticles()->getSearchRadius();
         const auto srcParams = this->modelParams.refined( searchRadius );
         const RealType h = srcParams.h;
         const RealType m = srcParams.mass;
         typename ParticlesType::NeighborsLoopParams searchInSrcBoundary( srcBoundary->getParticles() );

         const auto src_points = srcBoundary->getPoints().getConstView();
         const auto own_points = ownBoundary->getPoints().getConstView();

         auto interpolateFromSrc = [ = ] __cuda_callable__( GlobalIndexType i,
                                                            GlobalIndexType j,
                                                            VectorType& r_x,
                                                            MfdMatrixType* M_x,
                                                            MfdVectorType* brho_x,
                                                            MfdVectorType* bgamma_x ) mutable
         {
            const RealType rho_j = src_rho[ j ];
            const VectorType r_j = src_points[ j ];
            const VectorType r_xj = r_x - r_j;
            const RealType drs = l2Norm( r_xj );
            // source particles without fluid support carry rho = 0 and would make V_j infinite
            if( drs <= searchRadius && rho_j > 0.f ) {
               const RealType gamma_j = src_gamma[ j ];
               const RealType V_j = m / rho_j;

               *M_x += MFD::getPairCorrectionMatrix( r_xj, h ) * V_j;
               const MfdVectorType b_x = MFD::getPairVariableAndDerivatives( r_xj, h ) * V_j;
               *brho_x += rho_j * b_x;
               *bgamma_x += gamma_j * b_x;
            }
         };

         auto interpolateGhost = [ = ] __cuda_callable__( GlobalIndexType k ) mutable
         {
            const GlobalIndexType slot = invOwn[ ghostBegin + k ];
            const VectorType r_x = own_points[ slot ];
            MfdMatrixType M_x = 0.f;
            MfdVectorType brho_x = 0.f;
            MfdVectorType bgamma_x = 0.f;

            ParticlesType::NeighborsLoopAnotherSet::exec(
               slot, r_x, searchInSrcBoundary, interpolateFromSrc, &M_x, &brho_x, &bgamma_x );

            // keep the last values where the interpolation has no reliable support
            const RealType detM = Matrices::determinant( M_x );
            if( M_x( 0, 0 ) > 0.05 && detM > 0.001f ) {
               dst_rho[ slot ] = Matrices::solve( M_x, brho_x )[ 0 ];
               dst_gamma[ slot ] = Matrices::solve( M_x, bgamma_x )[ 0 ];
            }
            else if( M_x( 0, 0 ) > 0.05f ) {
               dst_rho[ slot ] = brho_x[ 0 ] / M_x( 0, 0 );
               dst_gamma[ slot ] = bgamma_x[ 0 ] / M_x( 0, 0 );
            }
         };
         Algorithms::parallelFor< DeviceType >( 0, rec.numberOfGhostParticles, interpolateGhost );
         continue;
      }

      const GlobalIndexType srcCount = srcBoundary->getNumberOfParticles();
      const auto srcRefIdx = srcBoundary->getVariables()->referentialIdx.getConstView();

      Algorithms::parallelFor< DeviceType >( 0, srcCount, [ = ] __cuda_callable__( GlobalIndexType j ) mutable
      {
         invSrc[ srcRefIdx[ j ] ] = j;
      } );

      const auto srcIdx = rec.srcIndexList.getConstView();

      Algorithms::parallelFor< DeviceType >( 0, rec.numberOfGhostParticles, [ = ] __cuda_callable__( GlobalIndexType k ) mutable
      {
         const GlobalIndexType ghostSlot = invOwn[ ghostBegin + k ];
         const GlobalIndexType srcSlot = invSrc[ srcIdx[ k ] ];
         dst_rho[ ghostSlot ] = src_rho[ srcSlot ];
         dst_gamma[ ghostSlot ] = src_gamma[ srcSlot ];
      } );
   }
}

#ifdef HAVE_MPI
template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                                                          TNL::Config::ParameterContainer& parametersDistributed,
                                                                          TNL::Logger& logger )
{
   logger.writeHeader( "SPH simulation initialization." );
   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   Containers::StaticVector< 2, int > subdomainCoordinates = distributed::restoreSubdomainCoordinatesFromRank( rank, numberOfSubdomains );
   const std::string subdomainKey = distributed::getSubdomainKey( rank, numberOfSubdomains );
   logger.writeParameter( "Initializing rank: ", rank );
   logger.writeParameter( "Initializing subdomain: ", subdomainCoordinates );
   logger.writeParameter( "Subdomain key: ", subdomainKey );
   logger.writeSeparator();

   const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
   const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
   const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
   const CoordinatesType domainGridDimension = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );
   const int numberOfOverlapLayers = parameters.getParameter< int >( "overlapWidth" );

   const VectorType subdomainOrigin = parametersDistributed.getXyz< VectorType >( subdomainKey + "origin" );
   const CoordinatesType subdomainGridDimension = parametersDistributed.getXyz< CoordinatesType >( subdomainKey + "grid-dimensions" );
   const CoordinatesType subdomainGridOriginGlobalCoords = parametersDistributed.getXyz< CoordinatesType >( subdomainKey + "origin-global-coords" );

   logger.writeParameter( "initDistributed:", "fluid->initialize" );
   this->fluid()->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "fluid_n" ),
                                            parametersDistributed.getParameter< int >( subdomainKey + "fluid_n_allocated" ),
                                            searchRadius,
                                            domainGridDimension,
                                            domainOrigin,
                                            subdomainGridDimension,
                                            subdomainGridOriginGlobalCoords,
                                            numberOfOverlapLayers,
                                            numberOfSubdomains,
                                            subdomainOrigin,
                                            logger );
   this->fluid()->getDistributedParticles()->writeProlog( logger );
   this->fluid()->setCommunicator( this->communicator );
   this->fluid()->getDistributedParticlesSynchronizer().initialize( this->fluid()->getDistributedParticles() );

   logger.writeParameter( "initDistributed:", "boundary->initialize" );
   this->boundary()->initializeAsDistributed( parametersDistributed.getParameter< int >( subdomainKey + "boundary_n" ),
                                               parametersDistributed.getParameter< int >( subdomainKey + "boundary_n_allocated" ),
                                               searchRadius,
                                               domainGridDimension,
                                               domainOrigin,
                                               subdomainGridDimension,
                                               subdomainGridOriginGlobalCoords,
                                               numberOfOverlapLayers,
                                               numberOfSubdomains,
                                               subdomainOrigin,
                                               logger );
   this->boundary()->getDistributedParticles()->writeProlog( logger );
   this->boundary()->setCommunicator( this->communicator );
   this->boundary()->getDistributedParticlesSynchronizer().initialize( this->boundary()->getDistributedParticles() );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
                                                                           TNL::Config::ParameterContainer& parametersDistributed,
                                                                           TNL::Logger& logger )
{
   int rank = TNL::MPI::GetRank();
   Containers::StaticVector< 2, int > numberOfSubdomains = parameters.getXyz< Containers::StaticVector< 2, int > >( "subdomains" );
   const std::string subdomainKey = distributed::getSubdomainKey( rank, numberOfSubdomains );

   logger.writeParameter( "Reading fluid particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   this->fluid()->template readParticlesAndVariables< typename BaseType::Reader >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "fluid-particles" ) );
   logger.writeParameter( "Reading boundary particles:", parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );
   this->boundary()->template readParticlesAndVariables< typename BaseType::Reader >(
      parametersDistributed.getParameter< std::string >( subdomainKey + "boundary-particles" ) );

   const int numberOfBoundaryPatches = parameters.getParameter< int >( "openBoundaryPatches" );
   if( this->openBoundaryPatches.size() > 0 ) {
      for( int i = 0; i < numberOfBoundaryPatches; i++ ) {
         std::string prefix = subdomainKey + "buffer-" + std::to_string( i + 1 ) + "-";
         this->openBoundaryPatches[ i ]->template readParticlesAndVariables< typename BaseType::Reader >(
            parametersDistributed.getParameter< std::string >( prefix + "particles" ) );
      }
   }
}
#endif

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::multiresolutionUpdate()
{
   this->timeMeasurement.start( "multiresolution-update" );

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatches.size(); p++ ) {
      const auto& iface = multiresolutionBoundaryPatchInterfaces[ p ];
      multiresolutionBoundaryPatches[ p ]->updateInterfaceBuffer(
            this->fluidSets[ iface.ownIdx ],
            this->fluidSets[ iface.neighborIdx ],
            this->modelParams,
            this->timeStepping.getTimeStep(),
            iface.ownIdx );
   }

   this->timeMeasurement.stop( "multiresolution-update" );
   this->writeLog( "Update multiresolution BC...", "Done.");
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::interact()
{
    this->timeMeasurement.start( "interact" );

    for( int i = 0; i < this->numberOfSubsets; i++ ){
       this->model.updateSolidBoundary( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
       this->model.finalizeBoundaryInteraction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
    }

    refreshBoundaryGhostValues();

    for( int i = 0; i < this->numberOfSubsets; i++ ){
       this->model.interaction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
       for( long unsigned int p = 0; p < multiresolutionBoundaryPatchInterfaces.size(); p++ )
          if( multiresolutionBoundaryPatchInterfaces[ p ].ownIdx == i )
             this->model.interactionWithOpenBoundary(
                   this->fluidSets[ i ], multiresolutionBoundaryPatches[ p ], this->modelParams );
       this->model.finalizeInteraction( this->fluidSets[ i ], this->boundarySets[ i ], this->modelParams );
    }

   this->timeMeasurement.stop( "interact" );
   this->writeLog( "Interact...", "Done." );
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::save( bool writeParticleCellIndex )
{
   if( ( this->verbose == "with-snapshot" ) || ( this->verbose == "full" ) )
      this->writeInfo();

   const RealType time = this->timeStepping.getTime();

   for( int i = 0; i < this->numberOfSubsets; i++ ){
#ifdef HAVE_MPI
      std::string outputFileNameFluid = this->outputDirectory + "/fluid_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameFluid = this->outputDirectory + "/fluid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->fluidSets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >( outputFileNameFluid, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameFluid );

#ifdef HAVE_MPI
      std::string outputFileNameBound = this->outputDirectory + "/boundary_rank" + std::to_string( TNL::MPI::GetRank() ) + "_" + std::to_string( time ) + "_particles.vtk";
#else
      std::string outputFileNameBound = this->outputDirectory + "/boundary_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + "_particles.vtk";
#endif
      this->boundarySets[ i ]->template writeParticlesAndVariables< typename BaseType::Writer >( outputFileNameBound, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameBound );

#ifdef HAVE_MPI
      std::string outputFileNameGrid = this->outputDirectory + "/grid_rank" + std::to_string( TNL::MPI::GetRank() + 1 ) + "_" + std::to_string( time ) + ".vtk";
#else
      std::string outputFileNameGrid = this->outputDirectory + "/grid_subdomain" + std::to_string( i ) + "_" + std::to_string( time ) + ".vtk";
#endif
      TNL::Particles::Writers::writeBackgroundGrid( outputFileNameGrid, this->fluidSets[ i ]->getParticles()->getDimensions(), this->fluidSets[ i ]->getParticles()->getOrigin(), this->fluidSets[ i ]->getParticles()->getSearchRadius() );
      this->logger.writeParameter( "Saved:", outputFileNameGrid );

      this->simulationMonitor.save( this->logger );
   }

   for( long unsigned int p = 0; p < multiresolutionBoundaryPatchInterfaces.size(); p++ ) {
      const auto& iface = multiresolutionBoundaryPatchInterfaces[ p ];
      std::string outputFileNameMultiresolutionBound =
         this->outputDirectory + "/multiresolutionBoundaryPatch_subdomain" + std::to_string( iface.ownIdx ) + "_to_"
            + std::to_string( iface.neighborIdx ) + "_" + std::to_string( time ) + "_particles.vtk";
      multiresolutionBoundaryPatches[ p ]->template writeParticlesAndVariables< typename BaseType::Writer >(
            outputFileNameMultiresolutionBound, writeParticleCellIndex );
      this->logger.writeParameter( "Saved:", outputFileNameMultiresolutionBound );
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::writeProlog( bool writeSystemInformation ) noexcept
{
   auto& log = this->logger;

   log.writeHeader( "SPH simulation configuration." );
   log.writeParameter( "Case name:", this->caseName );
   log.writeSeparator();

   const bool printGPUs = std::is_same< DeviceType, TNL::Devices::Cuda >::value;

   if( TNL::MPI::isInitialized() )
      log.writeParameter( "MPI processes:", TNL::MPI::GetSize() );
   log.writeParameter( "Device type:", TNL::getType< DeviceType >() );
   if( ! printGPUs ) {
      if( TNL::Devices::Host::isOMPEnabled() ) {
         log.writeParameter( "OMP enabled:", "yes", 1 );
         log.writeParameter( "OMP threads:", TNL::Devices::Host::getMaxThreadsCount(), 1 );
      }
      else
         log.writeParameter( "OMP enabled:", "no", 1 );
   }

   log.writeParameter( "Particle system type:", Model::ParticlesType::writeModelType() );
   log.writeParameter( "SPH model:", Model::writeModelType() );
   log.writeParameter( "Verbose:", this->verbose );
   log.writeParameter( "Output directory:", this->outputDirectory );
   log.writeParameter( "Particles format", this->particlesFormat );
   writePrologModel( log, this->modelParams );

   for( int i = 0; i < this->numberOfSubsets; i++ ){
      log.writeHeader( "Fluid " + std::to_string( i ) + " object information." );
      this->fluidSets[ i ]->writeProlog( log );
      log.writeHeader( "Boundary " + std::to_string( i ) + " object information:" );
      this->boundarySets[ i ]->writeProlog( log );
   }

   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ ) {
         log.writeHeader( "Open boundary buffer" + std::to_string( i + 1 ) + "." );
         this->openBoundaryPatches[ i ]->writeProlog( log );
         log.writeSeparator();
         this->openBoundaryPatches[ i ]->config.writeProlog( log );
      }
   }

   if( multiresolutionBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ ) {
         log.writeHeader( "Multiresolution boundary buffer" + std::to_string( i ) + "." );
         multiresolutionBoundaryPatches[ i ]->writeProlog( log, i );
      }
   }

   log.writeHeader( "System information." );
   if( writeSystemInformation ) {
      log.writeSystemInformation( printGPUs );
      log.writeSeparator();
      log.writeCurrentTime( "Started at:" );
      log.writeSeparator();
   }
}

template< typename Model >
void
SolverMultiSetBlockMultiresolution< Model >::writeInfo() noexcept
{
   auto& log = this->logger;

   log.writeSeparator();
   log.writeParameter( "Simulation time: " + std::to_string( this->timeStepping.getTime() )
                                  + " s, simulation step: " + std::to_string( this->timeStepping.getStep() ),
                               "" );
   log.writeCurrentTime( "Current time:" );
   for( int i = 0; i < this->numberOfSubsets; i++ ){
      log.writeParameter( "Number of fluid particles:", this->fluidSets[ i ]->getNumberOfParticles() );
      log.writeParameter( "Number of allocated fluid particles:", this->fluidSets[ i ]->getNumberOfAllocatedParticles() );
      log.writeParameter( "Number of boundary particles:", this->boundarySets[ i ]->getNumberOfParticles() );
      log.writeParameter( "Number of allocated boundary particles:", this->boundarySets[ i ]->getNumberOfAllocatedParticles() );
   }
   if( this->openBoundaryPatches.size() > 0 ) {
      for( long unsigned int i = 0; i < this->openBoundaryPatches.size(); i++ )
         log.writeParameter( "Number of buffer" + std::to_string( i + 1 ) + " particles:",
                                this->openBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   if( multiresolutionBoundaryPatches.size() > 0 ) {
       for( long unsigned int i = 0; i < multiresolutionBoundaryPatches.size(); i++ )
         log.writeParameter( "Number of mr-buffer " + std::to_string( i ) + " particles:",
               multiresolutionBoundaryPatches[ i ]->getNumberOfParticles() );
   }
   log.writeSeparator();
}

} // namespace SPH
} // namespace TNL
