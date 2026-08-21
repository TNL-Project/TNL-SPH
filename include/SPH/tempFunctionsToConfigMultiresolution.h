#pragma once
#include <TNL/Devices/Cuda.h>
#include <TNL/Devices/Host.h>
#include <TNL/MPI.h>
#include <TNL/Config/ConfigDescription.h>

namespace TNL {
namespace SPH {

//TODO: Merge configSetupDistributedSubdomain and configSubdomain, since these are basically the same
void
configSubdomain( int subdomain, TNL::Config::ConfigDescription& config )
{
   std::string subdomainKey = "subdomain-" + std::to_string( subdomain ) + "-";
   config.addRequiredEntry< std::string >( subdomainKey + "fluid-particles", "Input fluid particles file path." );
   config.addRequiredEntry< std::string >( subdomainKey + "boundary-particles", "Input boundary particles file path." );
   config.addEntry< int >( subdomainKey + "fluid_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "fluid_n_allocated", "The allocated number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "boundary_n", "The initial number of fluid particles.", 0 );
   config.addEntry< int >( subdomainKey + "boundary_n_allocated", "The allocated number of fluid particles.", 0 );

   config.addEntry< double >( subdomainKey + "refinement-factor", "Refinement factor of current subdomain.", 0. ); //different

   config.addEntry< double >( subdomainKey + "origin-x", "The origin of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "origin-y", "The origin of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "origin-z", "The origin of domain in z direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-x", "The origin of domain in global cell coords. in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-y", "The origin of domain in global cell coords. in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "origin-global-coords-z", "The origin of domain in global cell coords. in z direction.", 0. );

   config.addEntry< double >( subdomainKey + "size-x", "The size of domain in x direction.", 0. );
   config.addEntry< double >( subdomainKey + "size-y", "The size of domain in y direction.", 0. );
   config.addEntry< double >( subdomainKey + "size-z", "The size of domain in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-x", "The size of domain in cells in x direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-y", "The size of domain in cells in y direction.", 0. );
   config.addEntry< int >( subdomainKey + "grid-dimensions-z", "The size of domain in cells in z direction.", 0. );
}

} // SPH
} // TNL

