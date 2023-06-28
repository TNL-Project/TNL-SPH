#pragma once

#include "SPHTraits.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>

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
   getNumberOfParticles() const
   {
      return this->particles->getNumberOfParticles();
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
      variables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles() );
      integratorVariables->sortVariables( particles->getSortPermutations(), particles->getNumberOfParticles() );
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
   writeParticlesAndVariables( const std::string& outputFileName )
   {
      std::ofstream outputFileFluid ( outputFileName, std::ofstream::out );
      WriterType writer( outputFileFluid );
      writer.writeParticles( *particles );
      variables->writeVariables( writer, particles->getNumberOfParticles() );
   }

   //Some additional informations
   GlobalIndexType firstActiveParticle = 0;
   GlobalIndexType lastActiveParticle = 0;

   //Properties of physical object
   ParticlePointerType particles;
   VariablesPointerType variables;
   IntegratorVariablesPointerType integratorVariables;
};

}
}
}

