// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <SPH/solvers/SolverMultiSetBase.h>
#include <SPH/shared/remeshParticles.h>

namespace TNL::SPH {

template< typename Model >
class SolverMultiSetRemeshed : public SolverMultiSetBase< Model >
{
public:

   using BaseType = SolverMultiSetBase< Model >;
   using SimulationType = SolverMultiSetRemeshed< Model >;

   using typename BaseType::DeviceType;
   using typename BaseType::ModelType;
   using typename BaseType::ModelParams;
   using typename BaseType::ParticlesType;
   using typename BaseType::SPHConfig;
   using typename BaseType::GlobalIndexType;
   using typename BaseType::RealType;
   using typename BaseType::CoordinatesType;
   using typename BaseType::VectorType;

   using typename BaseType::FluidPointer;
   using typename BaseType::BoundaryPointer;
   using typename BaseType::SimulationMonitor;

   SolverMultiSetRemeshed( std::ostream& out = std::cout ) : BaseType( out ) {};

   void
   initRemeshedSimulation( int argc, char* argv[] );

   void
   interact();

   void
   remeshParticles();

};

}  // namespace TNL::SPH

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const TNL::SPH::SolverMultiSetRemeshed< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include <SPH/solvers/SolverMultiSetRemeshed.hpp>

