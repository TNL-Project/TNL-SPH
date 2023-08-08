#pragma once

#include "../../SPHTraits.h"
#include <thrust/gather.h>
#include "BoundaryConditionsTypes.h"

#ifdef HAVE_MPI
#include "../../shared/utils.h"
#endif

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHState >
class SPHFluidVariables
{
   public:
   using SPHConfig = typename SPHState::SPHConfig;
   using SPHFluidTraitsType = SPHFluidTraits< SPHConfig >;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using RealType = typename SPHFluidTraitsType::RealType;
   using ScalarArrayType = typename SPHFluidTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;
   using IndexArrayType = typename SPHFluidTraitsType::IndexArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, typename SPHConfig::DeviceType >;

   SPHFluidVariables( GlobalIndexType size )
   : rho( size ), drho ( size ), p( size ), v( size ), a( size ),
     rho_swap( size ), v_swap( size ) {}

   /* Variables - Fields */
   ScalarArrayType rho;
   ScalarArrayType drho;
   ScalarArrayType p;
   VectorArrayType v;
   VectorArrayType a;

   /* Additional variable fields to avoid inmpace sort. */
   ScalarArrayType rho_swap;
   VectorArrayType v_swap;

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData(), view_rho_swap.getArrayData() );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData(), view_v_swap.getArrayData() );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      auto view_map = map->getView();

      auto view_rho = rho.getView();
      auto view_v = v.getView();

      auto view_rho_swap = rho_swap.getView();
      auto view_v_swap = v_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_rho.getArrayData() + firstActiveParticle, view_rho_swap.getArrayData() + firstActiveParticle );
      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_v.getArrayData() + firstActiveParticle, view_v_swap.getArrayData() + firstActiveParticle );

      rho.swap( rho_swap );
      v.swap( v_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      reader.template readParticleVariable< ScalarArrayType, typename ScalarArrayType::ValueType >( rho, "Density" );
      reader.template readParticleVariable< VectorArrayType, typename ScalarArrayType::ValueType >( v, "Velocity" );
   }

   template< typename WriterType >
   void
   writeVariables( WriterType& writer, const GlobalIndexType& numberOfParticles, const GlobalIndexType firstActiveParticle = 0 )
   {
      writer.template writePointData< ScalarArrayType >( p, "Pressure", numberOfParticles, firstActiveParticle, 1 );
      writer.template writePointData< ScalarArrayType >( rho, "Density", numberOfParticles, firstActiveParticle, 1 );
      writer.template writeVector< VectorArrayType, RealType >( v, "Velocity", numberOfParticles, firstActiveParticle, 3 );
   }

#ifdef HAVE_MPI
   template< typename Synchronizer, typename SimulationSubdomainInfo >
   void
   synchronizeVariables( Synchronizer& synchronizer, SimulationSubdomainInfo& subdomainInfo )
   {
      synchronizer.template synchronizeArray< ScalarArrayType >( rho, rho_swap, subdomainInfo, 1 );
      synchronizer.template synchronizeArray< VectorArrayType >( v, v_swap, subdomainInfo, 1 );
   }

   void
   centerVariablesInMemory( const GlobalIndexType firstActiveParticle,
                            const GlobalIndexType shiftInMemory,
                            const GlobalIndexType numberOfParticles )
   {
      utils::shiftArray( rho, rho_swap, firstActiveParticle, shiftInMemory, numberOfParticles );
      utils::shiftArray( v, v_swap, firstActiveParticle, shiftInMemory, numberOfParticles );
   }
#endif

};

template< typename SPHState >
class SPHOpenBoundaryVariables : public SPHFluidVariables< SPHState >
{
   public:
   using BaseType = SPHFluidVariables< SPHState >;
   using SPHTraitsType = typename BaseType::SPHFluidTraitsType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using IndexArrayType = typename SPHTraitsType::IndexArrayType;

   SPHOpenBoundaryVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHState >( size ), particleMark( size ), receivingParticleMark( size ) {};

   IndexArrayType particleMark;
   IndexArrayType receivingParticleMark;
};

template< typename SPHState, typename Enable = void >
class SPHBoundaryVariables : public SPHFluidVariables< SPHState >
{};

template< typename SPHState >
class SPHBoundaryVariables< SPHState,
                            typename std::enable_if_t< std::is_same_v< typename SPHState::BCType, WCSPH_BCTypes::DBC > > >
: public SPHFluidVariables< SPHState >
{
public:
   using BaseType = SPHFluidVariables< SPHState >;
   using GlobalIndexType = typename BaseType::GlobalIndexType;

   SPHBoundaryVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHState >( size ) {};
};

template< typename SPHState >
class SPHBoundaryVariables< SPHState,
                            typename std::enable_if_t< std::is_same_v< typename SPHState::BCType, WCSPH_BCTypes::MDBC > > >
: public SPHFluidVariables< SPHState >
{
public:
   using Base = SPHFluidTraits< SPHState >;
   using SPHFluidTraitsType = SPHFluidTraits< typename SPHState::SPHConfig >;

   using GlobalIndexType = typename SPHFluidTraitsType::GlobalIndexType;
   using VectorArrayType = typename SPHFluidTraitsType::VectorArrayType;

   SPHBoundaryVariables( GlobalIndexType size )
   : SPHFluidVariables< SPHState >( size ), ghostNodes( size ), ghostNodes_swap( size ) {}

   VectorArrayType ghostNodes;
   VectorArrayType ghostNodes_swap;

   template< typename IndexArrayTypePointer >
   void
   sortVariables( IndexArrayTypePointer& map, GlobalIndexType numberOfParticles, GlobalIndexType firstActiveParticle )
   {
      Base::sortVariables( map, numberOfParticles, firstActiveParticle );

      auto view_map = map->getView();
      auto view_ghostNodes = ghostNodes.getView();
      auto view_ghostNodes_swap = ghostNodes_swap.getView();

      thrust::gather( thrust::device, view_map.getArrayData(), view_map.getArrayData() + numberOfParticles,
            view_ghostNodes.getArrayData() + firstActiveParticle, view_ghostNodes_swap.getArrayData() + firstActiveParticle );

      ghostNodes.swap( ghostNodes_swap );
   }

   template< typename ReaderType >
   void
   readVariables( ReaderType& reader )
   {
      Base::readVariables( reader );
      reader.template readParticleVariable2D< VectorArrayType, typename VectorArrayType::ValueType::ValueType >( ghostNodes, "GhostNodes" ); //FIXME!
   }

};

} // SPH
} // ParticleSystem
} // TNL

