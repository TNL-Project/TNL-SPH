#pragma once

#include <TNL/Logger.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Algorithms/reduce.h>

#include <limits> //UINT_MAX

#include <TNL/MPI/ScopedInitializer.h>

#include "../Particles/Particles.h"
#include "../Particles/ParticlesTraits.h"

#include "SPH.h"
#include "SPHTraits.h"

#include "DistributedSPHSynchronizer.h"

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename ParticleConfig >
class SimulationSubdomainInfo
{
   public:
   using ParticleTraitsType = ParticlesTraits< ParticleConfig, typename ParticleConfig::DeviceType >;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using RealType = typename ParticleTraitsType::RealType;
   using PointType = typename ParticleTraitsType::PointType;

   /**
    * Store the index of grid column (in 2D) or slince (in 3D), which
    * represents the origin and end of the current subdomain with respect
    * to the global grid.
    */
   GlobalIndexType gridIdxOverlapBegin;
   GlobalIndexType gridIdxBegin;
   GlobalIndexType gridIdxOverlapEnd;
   GlobalIndexType gridIdxEnd;

   /**
    * Store index of first and last particle of first and last
    * grid column (in 2D) or slice (in 3D) of the subdomain.
    */
   GlobalIndexType firstParticleInFirstGridColumn;
   GlobalIndexType lastParticleInFirstGridColumn;
   GlobalIndexType firstParticleInLastGridColumn;
   GlobalIndexType lastParticleInLastGridColumn;

   /**
    * Number of particles to send from current rank to its
    * left subdomain (numberOfParticlesToSet) and right
    * subdomain (numberOfParticlesToSendEnd).
    */
   GlobalIndexType numberOfParticlesToSendBegin = 0;
   GlobalIndexType numberOfParticlesToSendEnd = 0;

   /**
    * Number of particles received on current rank from its
    * left subdomain (receivedBegin) and right subdomain (receivedEnd).
    * This is necessary to incorporte the recieved data to
    * arrays that exits in the subdomain.
    */
   GlobalIndexType receivedBegin = 0;
   GlobalIndexType receivedEnd = 0;

   /**
    * Load balancing metadata - domainds are resized based on the
    * computation time on current, left and right subdomain and number of
    * particles in current, left and right subdomain.
    */
   GlobalIndexType numberOfParticlesInPreviousSubdomain;
   GlobalIndexType numberOfParticlesInNextSubdomain;
   GlobalIndexType numberOfParticlesInThisSubdomain;

   RealType solutionTimeInPreviousSubdomain;
   RealType solutionTimeInNextSubdomain;
   RealType solutionTimeInThisSubdomain;

   template< typename SubdomainInfoParams >
   void loadParameters( SubdomainInfoParams params )
   {
      this->gridIdxOverlapBegin = params.gridIdxOverlapStar;
      this->gridIdxBegin = params.gridIdxStart;

      this->gridIdxOverlapEnd = params.gridIdxOverlapEnd;
      this->gridIdxEnd = params.gridIdxEnd;

      //Local simulation info:
      this->firstParticleInFirstGridColumn = 0;
      this->lastParticleInFirstGridColumn = 0;
      this->firstParticleInLastGridColumn = 0;
      this->lastParticleInLastGridColumn = 0;
   }

   void
   writeProlog( TNL::Logger& logger ) const noexcept;
};


template< typename SPHSimulation >
class DistributedSPHSimpleFluid
{
public:

   using SPHSimulationType = SPHSimulation;
   using DeviceType = typename SPHSimulationType::DeviceType;
   using ParticleSystem = typename SPHSimulationType::ParticleSystemType;
   using ParticleTraitsType = typename SPHSimulation::ParticleSystemType::ParticleTraitsType;

   using SPHConfig = typename SPHSimulationType::ModelType::SPHConfig;
   using ParticleConfig = typename SPHSimulationType::ParticleSystemType::Config;

   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using RealType = typename ParticleTraitsType::RealType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using PairIndexType = typename ParticleTraitsType::PairIndexType;

   //using SimulationSubdomainInfo = SimulationSubdomainInfo< ParticleConfig >;
   using SimulationSubdomainInfo = DistributedPhysicalObjectInfo< ParticleConfig >;
   using Synchronizer = DistributedSPHSimulationSynchronizer< ParticleConfig, SPHConfig >;

   //Synchronizer types
   using SPHTriatsType = SPHFluidTraits< SPHConfig >;
   using IndexArrayType = Containers::Vector< GlobalIndexType, DeviceType, GlobalIndexType >;
   using ByteArrayView = Containers::ArrayView< std::uint8_t, DeviceType, GlobalIndexType >;
   using RequestsVector = std::vector< MPI_Request >;

   DistributedSPHSimpleFluid() = default;

   DistributedSPHSimpleFluid( SPHSimulation&& localSimulation )
   : localSimulation( std::move( localSimulation ) ) {}

   void
   setParticlesDecomposition( const IndexVectorType& domainDecomposition );

   const IndexVectorType&
   getParticleDimension() const;

   void
   setCommunicator( const MPI::Comm& communicator );

   //Steps we need to perform before synchronization:
   /**
    * Returns first and last particle in given column of cells.
    */
   template< typename SPHObjectPointer >
   PairIndexType
   getFirstLastParticleInColumnOfCells( const GlobalIndexType& gridColumn, const SPHObjectPointer& sphObject );

   /**
    * Obtains particle indices, ranges and limits necessary for
    * simulation synchronization.
    */
   void
   performLoadBalancing();

   void
   updateLocalSubdomain();

   void
   synchronize();

   template< typename SPHObjectPointer >
   void
   updateLocalSimulationInfo( SimulationSubdomainInfo& subdomainInfo, SPHObjectPointer& sphObject );

   template< typename Array >
   void
   synchronizeArray( Array& arraySend, Array& arrayReceive, SimulationSubdomainInfo& subdomainInfo, int valuesPerElement = 1 );

   void
   synchronizeByteArray( ByteArrayView arraySend,
                         ByteArrayView arrayReceive,
                         SimulationSubdomainInfo& subdomainInfo,
                         int bytesPerValue )
   {
      auto requests = synchronizeByteArrayAsyncWorker( arraySend, arrayReceive, subdomainInfo, bytesPerValue );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeByteArrayAsyncWorker( ByteArrayView arraySend,
                                    ByteArrayView arrayReceive,
                                    SimulationSubdomainInfo& subdomainInfo,
                                    int bytesPerValue );

   template< typename Array, typename SPHObjectPointer >
   void
   arrangeRecievedAndLocalData( Array& arraySend,
                                Array& arrayReceive,
                                SPHObjectPointer& sphObject,
                                SimulationSubdomainInfo& subdomainInfo,
                                bool tempSetNumberOfPtcs );

   //Load balancing
   void
   synchronizeSubdomainMetaData( SimulationSubdomainInfo& subdomainInfo )
   {
      auto requests = synchronizeSubdomainMetaDataArrayAsyncWorker( subdomainInfo );
      MPI::Waitall( requests.data(), requests.size() );
   }

   RequestsVector
   synchronizeSubdomainMetaDataArrayAsyncWorker( SimulationSubdomainInfo& subdomainInfo );

   void
   updateSubdomainSize( SimulationSubdomainInfo& subdomainInfo, SimulationSubdomainInfo& subdomainInfo_boundary );

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

   //implement standard interaction functions aswell for distributed
   template< typename SPHKernelFunction, typename DiffusiveTerm, typename ViscousTerm, typename EOS, typename SPHState >
   void
   interact( SPHState& sphState );

   template< typename Writer >
   void
   save( const std::string& outputFilename, const int step );


//protected:

   SPHSimulation localSimulation;
   SimulationSubdomainInfo localSimulationInfo;
   SimulationSubdomainInfo localSimulationInfo_boundary;

   MPI::Comm communicator = MPI_COMM_WORLD;
   Synchronizer synchronizer;


   //Control
   IndexVectorType domainDecomposition = 0;
   bool distributed = false;
   bool isSet = false;

};

} // SPH
} // ParticleSystem
} // TNL

template< typename SPHSimulation >
std::ostream&
operator<<( std::ostream& str, const SPH::DistributedSPHSimpleFluid< SPHSimulation >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

//DEBUG
template< typename ParticleConfig >
std::ostream&
operator<<( std::ostream& str, const SPH::SimulationSubdomainInfo< ParticleConfig >& subdomainInfo )
{
   TNL::Logger logger( 100, str );

   subdomainInfo.writeProlog( logger );

   return str;
}


#include "DistributedSPH_impl.h"

