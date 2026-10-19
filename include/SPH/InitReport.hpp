#include "InitReport.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <TNL/Containers/StaticVector.h>
#include <TNL/MPI.h>

namespace TNL {
namespace SPH {

inline nlohmann::json
makeInitReportDocument()
{
   nlohmann::json doc;
   doc[ "schema_version" ] = 1;
   doc[ "config" ] = nlohmann::json::object();
   doc[ "fluids" ] = nlohmann::json::array();
   doc[ "boundaries" ] = nlohmann::json::array();
   doc[ "open_boundary_patches" ] = nlohmann::json::array();
   doc[ "periodic_patches" ] = nlohmann::json::array();
   doc[ "multiresolution_buffers" ] = nlohmann::json::array();
   doc[ "decomposition" ] = nlohmann::json::object();
   return doc;
}

template< typename Params >
void
addConfigSection( nlohmann::json& doc, const Params& parameters )
{
   nlohmann::json config = nlohmann::json::object();

   // Domain geometry - parameter access mirrors DecompositionTopology.h.
   // The vectors always have three components (z is 0 for 2D configurations).
   using PointType = Containers::StaticVector< 3, double >;
   const PointType domainOrigin = parameters.template getXyz< PointType >( "domainOrigin" );
   const PointType domainSize = parameters.template getXyz< PointType >( "domainSize" );
   config[ "domainOrigin" ] = { domainOrigin[ 0 ], domainOrigin[ 1 ], domainOrigin[ 2 ] };
   config[ "domainSize" ] = { domainSize[ 0 ], domainSize[ 1 ], domainSize[ 2 ] };
   config[ "searchRadius" ] = parameters.template getParameter< double >( "searchRadius" );

   // Model resolution parameters - spelled exactly as in the JSONC config
   if( parameters.checkParameter( "dp" ) )
      config[ "dp" ] = parameters.template getParameter< double >( "dp" );
   if( parameters.checkParameter( "h" ) )
      config[ "h" ] = parameters.template getParameter< double >( "h" );

   // Decomposition: an absent entry (or the "unset" default 0 of configSetup)
   // means a single subdomain
   int numberOfSubdomains = 1;
   if( parameters.checkParameter( "numberOfSubdomains" ) )
      numberOfSubdomains = std::max( 1, parameters.template getParameter< int >( "numberOfSubdomains" ) );
   config[ "numberOfSubdomains" ] = numberOfSubdomains;

   doc[ "config" ] = std::move( config );
}

inline void
writeInitReportFile( const nlohmann::json& doc, const std::string& outputDirectory )
{
   std::string fileName = "init_report.json";
#ifdef HAVE_MPI
   if( TNL::MPI::isInitialized() && TNL::MPI::GetSize() > 1 )
      fileName = "init_report_subdomain" + std::to_string( TNL::MPI::GetRank() ) + ".json";
#endif
   const std::string reportPath = outputDirectory + "/" + fileName;

   std::ofstream reportFile( reportPath );
   if( ! reportFile.is_open() )
      throw std::runtime_error( "Failed to open the init report file: " + reportPath );
   reportFile << doc.dump( 2 ) << '\n';
}

}  //namespace SPH
}  //namespace TNL
