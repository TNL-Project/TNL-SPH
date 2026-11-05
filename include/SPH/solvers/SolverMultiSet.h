// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once

#include <SPH/solvers/SolverMultiSetBase.h>

namespace TNL::SPH {

template< typename Model >
class SolverMultiSet : public SolverMultiSetBase< Model >
{
public:

   using BaseType = SolverMultiSetBase< Model >;
   using SimulationType = SolverMultiSet< Model >;

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
   using typename BaseType::OpenBoundaryPointer;

   using typename BaseType::SimulationMonitor;

   SolverMultiSet( std::ostream& out = std::cout ) : BaseType( out ) {};

   void
   init( int argc, char* argv[] );

   void
   initParticleSets( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

   void
   readParticlesFiles( TNL::Config::ParameterContainer& parameters, TNL::Logger& logger );

#ifdef HAVE_MPI
   void
   initDistributedParticleSets( TNL::Config::ParameterContainer& parameters,
                                TNL::Config::ParameterContainer& parametersDistributed,
                                TNL::Logger& logger );

   void
   readParticleFilesDistributed( TNL::Config::ParameterContainer& parameters,
                                 TNL::Config::ParameterContainer& parametersDistributed,
                                 TNL::Logger& logger );
#endif

   void
   interact();

};

}  // namespace TNL::SPH

template< typename Model >
std::ostream&
operator<<( std::ostream& str, const TNL::SPH::SolverMultiSet< Model >& sphSimulation )
{
   TNL::Logger logger( 100, str );

   sphSimulation.writeProlog( logger );

   return str;
}

#include <SPH/solvers/SolverMultiSet.hpp>
