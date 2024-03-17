#pragma once

#include "SPHTraits.h"
#include "TNL/Logger.h"
#include <memory>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#if HAVE_MPI
#include "DistributedSPHSynchronizer.h"
#include "shared/utils.h"
#endif

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

   ParticleSet() : particles(), variables(), integratorVariables() {}

   ParticleSet( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h, numberOfCells ),
     variables( sizeAllocated ),
     integratorVariables( sizeAllocated ),
     firstActiveParticle( 0 ),
     lastActiveParticle( size - 1 ) {};

   void
   initialize( unsigned int numberOfParticles,
               unsigned int numberOfAllocatedParticles,
               RealType searchRadius,
               IndexVectorType gridSize,
               VectorType gridOrigin )
   {
      this->particles->setSize( numberOfAllocatedParticles );
      this->particles->setSearchRadius( searchRadius );
      this->particles->setGridSize( gridSize );
      this->particles->setGridOrigin( gridOrigin );
      this->particles->setNumberOfParticles( numberOfParticles );
      this->particles->setFirstActiveParticle( 0 );
      this->particles->setLastActiveParticle( numberOfParticles - 1 );
      this->firstActiveParticle = 0;
      this->lastActiveParticle = numberOfParticles - 1;
      this->variables->setSize( numberOfAllocatedParticles );
      this->integratorVariables->setSize( numberOfAllocatedParticles );
   }

   const GlobalIndexType
   getFirstActiveParticle() const
   {
      return this->firstActiveParticle;
   }

   void
   setFirstActiveParticle( GlobalIndexType firstActiveParticle )
   {
      this->firstActiveParticle = firstActiveParticle;
   }

   const GlobalIndexType
   getLastActiveParticle() const
   {
      return this->lastActiveParticle;
   }

   void
   setLastActiveParticle( GlobalIndexType lastActiveParticle )
   {
      this->lastActiveParticle = lastActiveParticle;
   }

   const GlobalIndexType
   getNumberOfActiveParticles() const
   {
      return ( this->lastActiveParticle - this->firstActiveParticle + 1 );
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
      return this->particles;
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

   void sortParticles()
   {
      particles->sortParticles();
      variables->sortVariables(
            particles->getSortPermutations(), particles->getNumberOfParticles(), particles->getFirstActiveParticle() );
      integratorVariables->sortVariables(
            particles->getSortPermutations(), particles->getNumberOfParticles(), particles->getFirstActiveParticle() );
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
      variables->writeVariables( writer, particles->getNumberOfParticles(), particles->getFirstActiveParticle() );

      if( writeParticleCellIndex == true )
         writer.template writePointData< typename ParticleSystem::CellIndexArrayType >(
               particles->getParticleCellIndices(),
               "GridIndex",
               particles->getNumberOfParticles(),
               particles->getFirstActiveParticle(),
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
   template< typename Synchronzier >
   void
   synchronizeObject( Synchronzier& synchronizer )
   {
      variables->synchronizeVariables( synchronizer, subdomainInfo );
      integratorVariables->synchronizeVariables( synchronizer, subdomainInfo );
      synchronizer.template synchronizeArray< typename ParticleSystem::PointArrayType >(
            particles->getPoints(), particles->getPointsSwap(), subdomainInfo, 1 );
   }
#endif

   void
   writeProlog( TNL::Logger& logger ) const noexcept
   {
      logger.writeParameter( "First active particle index:", this->firstActiveParticle );
      logger.writeParameter( "Last active particle index:", this->lastActiveParticle );
      logger.writeParameter( "Particle system parameters:", "" );
      particles->writeProlog( logger );
   }

   //Some additional informations
   GlobalIndexType firstActiveParticle = 0;
   GlobalIndexType lastActiveParticle = 0;

   //Properties of physical object
   ParticlePointerType particles;
   VariablesPointerType variables;
   IntegratorVariablesPointerType integratorVariables;

#ifdef HAVE_MPI
   using SimulationSubdomainInfo = DistributedParticleSetInfo< typename ParticleSystem::Config >;
   SimulationSubdomainInfo subdomainInfo;
#endif

};

}
}

