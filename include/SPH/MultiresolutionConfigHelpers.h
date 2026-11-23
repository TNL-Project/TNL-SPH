// SPDX-FileComment: This file is part of TNL-SPH (https://gitlab.com/tnl-project/tnl-sph)
// SPDX-License-Identifier: MIT

#pragma once
#include <TNL/Devices/Cuda.h>
#include <TNL/Devices/Host.h>
#include <TNL/MPI.h>
#include <TNL/Config/ConfigDescription.h>

namespace TNL::SPH {

//TODO: Merge configSetupDistributedSubdomain and configSubdomain, since these are basically the same
void
configSubdomain( int subdomain, TNL::Config::ConfigDescription& config )
{
   std::string subdomainKey = "subdomain" + std::to_string( subdomain );
   config.addRequiredEntry< std::string >( subdomainKey + "FluidParticles", "Input fluid particles file path." );
   config.addRequiredEntry< std::string >( subdomainKey + "BoundaryParticles", "Input boundary particles file path." );
   config.addEntry< int >( subdomainKey + "Fluid_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Fluid_n_allocated", "The allocated number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Boundary_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "Boundary_n_allocated", "The allocated number of fluid particles.", 0 );

   config.addEntry< double >( subdomainKey + "RefinementFactor", "Refinement factor of current subdomain.", 0. );  //different

   config.addEntry< double >( subdomainKey + "Origin-x", "The origin of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "Origin-y", "The origin of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "Origin-z", "The origin of domain in z direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoords-x", "The origin of domain in global cell coords. in x direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoords-y", "The origin of domain in global cell coords. in y direction.", 0. );
   config.addEntry< int >(
      subdomainKey + "OriginGlobalCoords-z", "The origin of domain in global cell coords. in z direction.", 0. );

   config.addEntry< double >( subdomainKey + "Size-x", "The size of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "Size-y", "The size of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "Size-z", "The size of domain in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensions-x", "The size of domain in cells in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensions-y", "The size of domain in cells in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "GridDimensions-z", "The size of domain in cells in z direction.", 0. );
}

}  // namespace TNL::SPH
