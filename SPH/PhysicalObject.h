#pragma once

#include "SPHTraits.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

#if HAVE_MPI
#include "DistributedSPHSynchronizer.h"
#include "shared/utils.h"
#endif

namespace TNL {
namespace ParticleSystem {
namespace SPH {

class PhysicalObjectMetada
{

};

template< typename ParticleSystem, typename SPHCaseConfig, typename Variables, typename IntegratorVariables >
class PhysicalObject
{
   public:
   using DeviceType = typename ParticleSystem::Device;
   using ParticlePointerType = typename Pointers::SharedPointer< ParticleSystem, DeviceType >;
   using VariablesPointerType = typename Pointers::SharedPointer< Variables, DeviceType >;
   using IntegratorVariablesPointerType = typename Pointers::SharedPointer< IntegratorVariables, DeviceType >;

   using SPHTraitsType = SPHFluidTraits< SPHCaseConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   PhysicalObject( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType h, GlobalIndexType numberOfCells )
   : particles( size, sizeAllocated, h, numberOfCells ),
     variables( sizeAllocated ),
     integratorVariables( sizeAllocated ),
     firstActiveParticle( 0 ),
     lastActiveParticle( size - 1 ) {};

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

   void
   centerObjectArraysInMemory()
   {
      //: const GlobalIndexType numberOfParticles = this->getNumberOfActiveParticles();
      //: const GlobalIndexType numberOfAllocatedParticles = this->getNumberOfAllocatedParticles();
      //: const GlobalIndexType shiftInMemory = static_cast< int >( ( numberOfAllocatedParticles - numberOfParticles ) / 2 );

      //: variables->centerVariablesInMemory( this->firstActiveParticle, shiftInMemory, numberOfParticles );
      //: integratorVariables->centerVariablesInMemory( this->firstActiveParticle, shiftInMemory, numberOfParticles );

      //: utils::shiftArray(
      //:       particles->getPoints(), particles->getPointsSwap(), this->firstActiveParticle, shiftInMemory, numberOfParticles );

      //: this->firstActiveParticle = shiftInMemory;
      //: this->lastActiveParticle = shiftInMemory + numberOfParticles - 1 ;
      //: this->particles->setFirstActiveParticle( shiftInMemory ); //FIXME: is this OK?
      //: this->particles->setLastActiveParticle( shiftInMemory + numberOfParticles - 1 ); //FIXME: is this OK?

      //edit
      const GlobalIndexType particlesStart = this->particles->getFirstActiveParticle();
      const GlobalIndexType numberOfParticlesToCopy = this->particles->getLastActiveParticle() -
                                                      this->particles->getFirstActiveParticle() + 1;
      //----- debug ------------------------------------------------------
      std::cout << "| particles->getNumberOfParticles(): " << particles->getNumberOfParticles() << " numberOfParticlesToCopy: " << numberOfParticlesToCopy << std::endl;
      if( particles->getNumberOfParticles() != numberOfParticlesToCopy )
         exit(1);
      //----- end-debug --------------------------------------------------

      const GlobalIndexType numberOfAllocatedParticles = this->getNumberOfAllocatedParticles();
      const GlobalIndexType shiftInMemory = static_cast< int >( ( numberOfAllocatedParticles - numberOfParticlesToCopy ) / 2 );

      variables->centerVariablesInMemory( particlesStart, shiftInMemory, numberOfParticlesToCopy );
      integratorVariables->centerVariablesInMemory( particlesStart, shiftInMemory, numberOfParticlesToCopy );

      utils::shiftArray(
            particles->getPoints(), particles->getPointsSwap(), particlesStart, shiftInMemory, numberOfParticlesToCopy );

      //experiment
      //this->firstActiveParticle = firstActiveParticle + shiftInMemory;
      //this->lastActiveParticle = lastActiveParticle + shiftInMemory;
      this->firstActiveParticle =  shiftInMemory + subdomainInfo.receivedBegin;
      this->lastActiveParticle = shiftInMemory + numberOfParticlesToCopy + subdomainInfo.receivedEnd - 1;

      this->particles->setFirstActiveParticle( shiftInMemory ); //FIXME: is this OK?
      this->particles->setLastActiveParticle( shiftInMemory + numberOfParticlesToCopy - 1 ); //FIXME: is this OK?

      //----- debug ------------------------------------------------------
      //TNL::MPI::Barrier();
      if( TNL::MPI::GetRank() == 0 ){
         std::cout << "rank 0:" << std::endl;
         std::cout << "| shift in memory: " << shiftInMemory << std::endl;
         std::cout << "| first active particle: " << firstActiveParticle << std::endl;
         std::cout << "| last active particle: " << lastActiveParticle << std::endl;
         std::cout << "| particles - first active particle: " << particles->getFirstActiveParticle() << std::endl;
         std::cout << "| particles - last active particle: " << particles->getLastActiveParticle() << std::endl;
      }
      //TNL::MPI::Barrier();
      if( TNL::MPI::GetRank() == 1 ){
         std::cout << "rank 1:" << std::endl;
         std::cout << "| shift in memory: " << shiftInMemory << std::endl;
         std::cout << "| first active particle: " << firstActiveParticle << std::endl;
         std::cout << "| last active particle: " << lastActiveParticle << std::endl;
         std::cout << "| particles - first active particle: " << particles->getFirstActiveParticle() << std::endl;
         std::cout << "| particles - last active particle: " << particles->getLastActiveParticle() << std::endl;
      }
      //TNL::MPI::Barrier();
      if( TNL::MPI::GetRank() == 2 ){
         std::cout << "rank 2:" << std::endl;
         std::cout << "| shift in memory: " << shiftInMemory << std::endl;
         std::cout << "| first active particle: " << firstActiveParticle << std::endl;
         std::cout << "| last active particle: " << lastActiveParticle << std::endl;
         std::cout << "| particles - first active particle: " << particles->getFirstActiveParticle() << std::endl;
         std::cout << "| particles - last active particle: " << particles->getLastActiveParticle() << std::endl;
      }
      //TNL::MPI::Barrier( distributedSPHSimulation.communicator );
      //----- end-debug --------------------------------------------------
   }

   void
   completeSynchronization()
   {
      //Update the particle ranges
      //GlobalIndexType updatedNumberOfParticles = synchronizer.getNewParticleCount(
      //      subdomainInfo, particles->getNumberOfParticles() );

      const GlobalIndexType numberOfParticlesToSet = subdomainInfo.lastParticleInLastGridColumn -
                                                     subdomainInfo.firstParticleInFirstGridColumn +
                                                     subdomainInfo.receivedEnd +
                                                     subdomainInfo.receivedBegin + 1;

      particles->setNumberOfParticles( numberOfParticlesToSet );
      particles->setFirstActiveParticle( subdomainInfo.firstParticleInFirstGridColumn - subdomainInfo.receivedBegin );
      particles->setLastActiveParticle( subdomainInfo.lastParticleInLastGridColumn + subdomainInfo.receivedEnd );
   }
#endif

   //Some additional informations
   GlobalIndexType firstActiveParticle = 0;
   GlobalIndexType lastActiveParticle = 0;

   //Properties of physical object
   ParticlePointerType particles;
   VariablesPointerType variables;
   IntegratorVariablesPointerType integratorVariables;

#ifdef HAVE_MPI
   using SimulationSubdomainInfo = DistributedPhysicalObjectInfo< typename ParticleSystem::Config >;
   SimulationSubdomainInfo subdomainInfo;
#endif

};

}
}
}

