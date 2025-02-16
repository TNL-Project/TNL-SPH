#pragma once

#include "SPHTraits.h"
#include "TNL/Logger.h"
#include <memory>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#include "OpenBoundaryConfig.h"
#include "PeriodicBoundaryBuffers.h"

#if HAVE_MPI
#include "DistributedSPHSynchronizer.h"
#include "shared/utils.h"
#endif

#include <Particles/DistributedParticles.h>
#include <Particles/DistributedParticlesSynchronizer.h>

namespace TNL {
namespace SPH {

class ParticleSetMetada
{

};

template< typename ParticleSystem, typename SPHCaseConfig, typename Variables, typename IntegratorVariables >
class ParticleSet
{
   public:
   using DeviceType = typename ParticleSystem::Device;
   using ParticlePointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using VariablesPointerType = typename Pointers::SharedPointer< Variables, DeviceType >;
   using IntegratorVariablesPointerType = typename Pointers::SharedPointer< IntegratorVariables, DeviceType >;

   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using IndexVectorType = typename SPHTraitsType::IndexVectorType;
   using VectorType = typename SPHTraitsType::VectorType;

   using OpenBoundaryConfig = OpenBoundaryConfig< SPHCaseConfig >;
   using PeriodicBoundary = PeriodicBoundary< ParticleSystem, OpenBoundaryConfig >;
   using PeriodicBoundaryPointer = typename Pointers::SharedPointer< PeriodicBoundary, DeviceType >;

#ifdef  HAVE_MPI
   using DistributedParticlesType = TNL::ParticleSystem::DistributedParticleSystem< ParticleSystem >;
   using DistributedParticlesPointerType = typename Pointers ::SharedPointer< DistributedParticlesType, DeviceType >;
   using Synchronizer = TNL::ParticleSystem::DistributedParticlesSynchronizer< DistributedParticlesType >;
#endif

   ParticleSet() : particles(), variables(), integratorVariables() {}

   ParticleSet( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h, numberOfCells ),
     variables( sizeAllocated ),
     integratorVariables( sizeAllocated ) {};

   void
   initialize( const unsigned int numberOfParticles,
               const unsigned int numberOfAllocatedParticles,
               const RealType& searchRadius,
               const IndexVectorType& gridDimension,
               const VectorType& gridOrigin )
   {
      this->getParticles()->setSize( numberOfAllocatedParticles );
      this->getParticles()->setSearchRadius( searchRadius );
      this->getParticles()->setGridDimensions( gridDimension );
      this->getParticles()->setGridOrigin( gridOrigin );
      this->getParticles()->setNumberOfParticles( numberOfParticles );
      this->variables->setSize( numberOfAllocatedParticles );
      this->integratorVariables->setSize( numberOfAllocatedParticles );
      const VectorType zeroVector = 0;
      this->getParticles()->setGridOriginGlobalCoords( zeroVector );
   }

#ifdef HAVE_MPI
   void
   initializeAsDistributed( const unsigned int numberOfParticles,
                            const unsigned int numberOfAllocatedParticles,
                            const RealType& searchRadius,
                            const IndexVectorType& domainGridDimension,
                            const VectorType& domainOrigin,
                            const IndexVectorType& subdomainGridDimension,
                            const IndexVectorType& subdomainGridOriginGlobalCoords,
                            const int& numberOfOverlapLayers,
                            const Containers::StaticVector< 2, int >& numberOfSubdomains,
                            TNL::Logger& logger )
   {
      const VectorType shiftOriginDueToOverlaps =  searchRadius * numberOfOverlapLayers;

      this->getParticles()->setSize( numberOfAllocatedParticles );
      this->getParticles()->setSearchRadius( searchRadius );
      this->getParticles()->setGridDimensions( domainGridDimension );
      this->getParticles()->setGridOrigin( domainOrigin );
      this->getParticles()->setOverlapWidth( numberOfOverlapLayers );
      this->getParticles()->setNumberOfParticles( numberOfParticles );
      this->getParticles()->setGridReferentialOrigin( domainOrigin - shiftOriginDueToOverlaps );
      this->getParticles()->setGridOriginGlobalCoords( subdomainGridOriginGlobalCoords );
      this->variables->setSize( numberOfAllocatedParticles );
      this->integratorVariables->setSize( numberOfAllocatedParticles );

      this->distributedParticles->setDistributedGridParameters( searchRadius,
                                                                domainGridDimension,
                                                                domainOrigin,
                                                                subdomainGridDimension,
                                                                subdomainGridOriginGlobalCoords,
                                                                numberOfOverlapLayers,
                                                                numberOfSubdomains );
   }
#endif

   void
   initializePeriodicity( TNL::Config::ParameterContainer& parameters )
   {
      //TODO: I don't like the compute domain properties here, this class should not take parameters as arg.
      const VectorType domainOrigin = parameters.getXyz< VectorType >( "domainOrigin" );
      const VectorType domainSize = parameters.getXyz< VectorType >( "domainSize" );
      const RealType searchRadius = parameters.getParameter< RealType >( "searchRadius" );
      const IndexVectorType gridSize = TNL::ceil( ( domainSize - domainOrigin ) / searchRadius );

      const int numberOfPeriodicPatches = parameters.getParameter< int >( "periodicBoundaryPatches" );
      std::cout << "Number of periodic patches: " << numberOfPeriodicPatches << std::endl;
      periodicPatches.resize( numberOfPeriodicPatches );
      for( int i = 0; i < numberOfPeriodicPatches; i++ ) {
         std::string prefix = "buffer-" + std::to_string( i + 1 ) + "-";
         periodicPatches[ i ]->initialize( parameters,
                                           prefix,
                                           searchRadius,
                                           gridSize,
                                           domainOrigin );
                                           //parameters.getParameter< int >( prefix + "numberOfParticlesPerCell" ) );
      }
   }

   const GlobalIndexType
   getNumberOfParticles() const
   {
      return this->particles->getNumberOfParticles();
   }

   const GlobalIndexType
   getNumberOfAllocatedParticles() const
   {
      return this->particles->getNumberOfAllocatedParticles();
   }

   ParticlePointerType&
   getParticles()
   {
      #ifdef HAVE_MPI
         return this->distributedParticles->getLocalParticles();
      #else
         return this->particles;
      #endif
   }

   const ParticlePointerType&
   getParticles() const
   {
      return this->particles;
   }

   typename ParticleSystem::PointArrayType&
   getPoints()
   {
      return this->particles->getPoints();
   }

   const typename ParticleSystem::PointArrayType&
   getPoints() const
   {
      return this->particles->getPoints();
   }

   virtual VariablesPointerType&
   getVariables()
   {
      return this->variables;
   }

   virtual const VariablesPointerType&
   getVariables() const
   {
      return this->variables;
   }

   void
   sortParticles()
   {
      particles->sortParticles();
      variables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles());
      integratorVariables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles() );
   }

   void
   sortVariables( const GlobalIndexType numberOfParticlesToRemove = 0 )
   {
      variables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles() + numberOfParticlesToRemove );
      integratorVariables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles() + numberOfParticlesToRemove );
   }

   void
   searchForNeighbors()
   {
      const GlobalIndexType numberOfParticlesToRemove = particles->getNumberOfParticlesToRemove();
      this->particles->searchForNeighbors();
      this->sortVariables( numberOfParticlesToRemove );

   }

   void
   makeSetSearchable()
   {
       if constexpr( ParticleSystem::specifySearchedSetExplicitly() == true ){
         const GlobalIndexType numberOfParticlesToRemove = particles->getNumberOfParticlesToRemove();
         this->particles->makeSetSearchable();
         this->sortVariables( numberOfParticlesToRemove );
       }
       else if constexpr( ParticleSystem::specifySearchedSetExplicitly() == false ){
         const GlobalIndexType numberOfParticlesToRemove = particles->getNumberOfParticlesToRemove();
         this->particles->searchForNeighbors();
         this->sortVariables( numberOfParticlesToRemove );
      }
   }

   void
   enforcePeriodicPatches()
   {
      for( long unsigned int i = 0; i < std::size( periodicPatches ); i++ ){
         periodicPatches[ i ]->particleZone.updateParticlesInZone( particles );
      }
   }

   template< typename ReaderType >
   void
   readParticlesAndVariables( const std::string& inputFileName )
   {
      ReaderType reader( inputFileName, particles->getNumberOfParticles(), particles->getNumberOfAllocatedParticles() );
      reader.template readParticles< typename ParticleSystem::PointArrayType >( particles->getPoints() ) ;
      variables->readVariables( reader );
   }

   template< typename WriterType >
   void
   writeParticlesAndVariables( const std::string& outputFileName, bool writeParticleCellIndex = false )
   {
      std::ofstream outputFileFluid ( outputFileName, std::ofstream::out );
      WriterType writer( outputFileFluid );
      writer.writeParticles( *particles );
      variables->writeVariables( writer, particles->getNumberOfParticles() );

      if( writeParticleCellIndex == true )
         writer.template writePointData< typename ParticleSystem::CellIndexArrayType >(
               particles->getParticleCellIndices(),
               "GridIndex",
               particles->getNumberOfParticles(),
               1 );
   }

   void
   writeProlog( TNL::Logger& logger )
   {
      logger.writeParameter( "Number of particles:", this->particles->getNumberOfParticles() );
      logger.writeParameter( "Number of allocated particles:", this->particles->getNumberOfAllocatedParticles() );
      logger.writeParameter( "Search radius:", this->particles->getSearchRadius() );
      logger.writeParameter( "Grid size:", this->particles->getSearchRadius() );
   }

#ifdef HAVE_MPI
   template< typename OverlapSetPointer >
   void
   synchronizeObject( OverlapSetPointer& overlapSet, TNL::Logger& logger )
   {
      this->distributedParticles->collectParticlesInInnerOverlaps( particles ); //TODO: Merge ptcs and distPtcs
      this->synchronizer.synchronizeOverlapSizes( distributedParticles, particles );
      // check numberOfParitlces, numberOfAllocatedParticles and numberOfRecvParticles

      // sychronize
      this->synchronizer.synchronize( this->getPoints(), overlapSet->getPoints(), distributedParticles );
      this->variables->synchronizeVariables( synchronizer, overlapSet->getVariables(), distributedParticles );
      this->integratorVariables->synchronizeVariables( synchronizer, overlapSet->integratorVariables, distributedParticles );

      // update the number of particles inside subdomain
      const GlobalIndexType numberOfRecvParticles = this->synchronizer.getNumberOfRecvParticles();
      particles->setNumberOfParticles( particles->getNumberOfParticles() + numberOfRecvParticles );
   }

   void
   synchronizeBalancingMeasures()
   {
      this->synchronizer.synchronizeBalancingMeasures( distributedParticles );
   }
#endif

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      logger.writeParameter( "Particle system parameters:", "" );
      particles->writeProlog( logger );
   }

   //Properties of physical object
   ParticlePointerType particles;
#ifdef HAVE_MPI
   DistributedParticlesPointerType distributedParticles;
   Synchronizer synchronizer;
#endif
   VariablesPointerType variables;
   IntegratorVariablesPointerType integratorVariables;

   std::vector< PeriodicBoundaryPointer > periodicPatches;


};

}
}

