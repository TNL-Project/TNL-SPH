#pragma once

#include <vector>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

/*
   All grid-level parameters for one subdomain, extracted from config once
   and reused by both initParticleSets and initZones.
*/
template< typename GridType >
struct SubdomainDescriptor
{
   using CoordinatesType = typename GridType::CoordinatesType;

   int index;
   float refinementFactor;
   GridType grid;
   CoordinatesType originCoordinates; //TODO: We need coordinates of the dmoain in terms of global grid. (Do we?)
};


struct InterfaceDescriptor
{
   int ownIdx;
   int neighborIdx;
};

/*
   Base topology.
*/
class DecompositionTopology
{
public:
   using Interface = InterfaceDescriptor;

   virtual ~DecompositionTopology() = default;
   virtual int getNumberOfSubdomains() const = 0;
   virtual std::vector< Interface > getInterfacesOfSubdomain( int idx ) const = 0;
};

/*
   Flat topology - linear and nested both supported via explicit addInterface()
*/
template< typename GridType>
class DecompositionTopologyFlat : public DecompositionTopology
{
public:

   using RealType = typename GridType::RealType;
   using PointType = typename GridType::PointType;
   using CoordinatesType = typename GridType::CoordinatesType;
   using Subdomain = SubdomainDescriptor< GridType >;
   using Interface = InterfaceDescriptor;

   // Construction from config - call once, reuse everywhere
   void loadFromConfig( TNL::Config::ParameterContainer& parameters,
                        TNL::Config::ParameterContainer& parametersSubdomains )
   {
      const RealType radius = parameters.getParameter< RealType >( "searchRadius" );
      const PointType origin = parameters.getXyz< PointType >( "domainOrigin" );
      const PointType size = parameters.getXyz< PointType >( "domainSize" );

      // Global domain as a TNL grid
      GridType global;
      global.setDimensions( TNL::ceil( ( size - origin ) / radius ) );
      global.setOrigin( origin );
      global.setSpaceSteps( PointType( radius ) );
      setGlobalGrid( global );
      numberOfOverlapLayers = 1; // TODO: from config, but for pure particles, it is 1 by default

      // Build local domains a s TNL grids
      const int n = parameters.getParameter< int >( "numberOfSubdomains" );
      for( int i = 0; i < n; i++ ) {
         const std::string key = "subdomain-" + std::to_string( i ) + "-";
         const float rf = parametersSubdomains.getParameter< float >( key + "refinement-factor" );
         const RealType local_radius = rf * radius;

         GridType local;
         local.setDimensions( parametersSubdomains.getXyz< CoordinatesType >( key + "grid-dimensions" ) );
         local.setOrigin( parametersSubdomains.getXyz< PointType >( key + "origin" ) ); //TODO: Remove?
         local.setSpaceSteps( PointType( local_radius ) );

         SubdomainDescriptor< GridType > sd;
         sd.index = i;
         sd.refinementFactor = rf;
         sd.grid = local;
         sd.originCoordinates = parametersSubdomains.getXyz< CoordinatesType >( key + "origin-global-coords" );
         subdomains.push_back( sd );
      }
   }

   void addInterface( int subdomain_idx_a, int subdomain_idx_b )
   {
      interfaces.push_back( { subdomain_idx_a, subdomain_idx_b } );
      interfaces.push_back( { subdomain_idx_b, subdomain_idx_a } );
   }

   void finalizeLinear()
   {
      for( int i = 0; i + 1 < ( int )subdomains.size(); ++i )
         addInterface( subdomains[ i ].index, subdomains[ i + 1 ].index );
   }

   // DecompositionTopology interface
   int getNumberOfSubdomains() const override
   {
      return subdomains.size();
   }

   std::vector< Interface > getInterfacesOfSubdomain( int idx ) const override
   {
      std::vector< Interface > result;
      for( const auto& iface : interfaces )
         if( iface.ownIdx == idx ) result.push_back( iface );
      return result;
   }

   const Subdomain& getSubdomain( int idx ) const
   {
      for( const auto& sd : subdomains )
         if( sd.index == idx ) return sd;
      throw std::out_of_range( "DecompositionTopologyFlat: unknown subdomain index" );
   }

   int
   getNumberOfSubdomainInterfaces() const
   {
      int totalPatches = 0;
      for( int i = 0; i < getNumberOfSubdomains(); i++ )
         totalPatches += getInterfacesOfSubdomain( i ).size();
      return totalPatches;
   }

   // Global domain accessors
   void
   setGlobalGrid( const GridType& grid )
   {
      globalGrid = grid;
   }

   const GridType&
   getGlobalGrid() const
   {
      return globalGrid;
   }

   RealType
   getGlobalSearchRadius() const
   {
      return globalGrid.getSpaceSteps()[ 0 ];
   }

   PointType
   getDomainOrigin() const
   {
      return globalGrid.getOrigin();
   }

   PointType
   getDomainSize() const
   {
      return globalGrid.getProportions();
   }

   CoordinatesType
   getDomainGridDimension() const
   {
      return globalGrid.getDimensions();
   }

   int
   getNumberOfOverlapLayers() const
   {
      return numberOfOverlapLayers;
   }


   // Local domain accessors
   const GridType&
   getLocalGrid( int subdomain_idx ) const
   {
      return getSubdomain( subdomain_idx ).grid;
   }

   CoordinatesType
   getLocalOriginCoordinates( int subdomain_idx )
   {
      return getSubdomain( subdomain_idx ).originCoordinates;
   }

   /**
    * \brief Collects the decomposition geometry into the initialization report.
    *
    * Emits the global (coarse) grid built by loadFromConfig and one entry per
    * subdomain of the flat topology - a plain read-back of the values as loaded
    * (this header has no writeProlog sibling - the report is a new standalone
    * method, not a mirror of an existing log section).
    *
    * The method exists only on the flat topology class: the %DecompositionTopology
    * base interface exposes neither the global grid nor the per-subdomain
    * descriptors, so there is no meaningful base-class override to give.
    *
    * \b region_dimensions is the footprint of each subdomain expressed in GLOBAL
    * coarse-grid coordinates (cells): `grid-dimensions * refinement-factor` per
    * axis. For a 1D x-split chain this reproduces the consecutive-origin
    * differences (the last region extends to the end of the global grid); for
    * nested refinement chains it is the patch footprint in coarse cells. A
    * consecutive-origin shortcut is intentionally NOT used - it misreports
    * nested patches, which are the actual consumers of this topology.
    * Vector entries are three-component JSON arrays (z = 0 in 2D), mirroring
    * the config-section serialization of InitReport.hpp.
    *
    * \returns A JSON object with the keys \e global (\e origin, \e dimensions)
    *          and \e subdomains (entries with \e index, \e origin_global_coordinates,
    *          \e region_dimensions).
    */
   [[nodiscard]] nlohmann::json
   writeInitReport() const
   {
      nlohmann::json report = nlohmann::json::object();
      nlohmann::json global = nlohmann::json::object();
      global[ "origin" ] = vectorToInitReportArray( getDomainOrigin() );
      global[ "dimensions" ] = vectorToInitReportArray( getDomainGridDimension() );
      report[ "global" ] = std::move( global );

      report[ "subdomains" ] = nlohmann::json::array();
      for( const Subdomain& sd : subdomains ) {
         // footprint in global coarse-grid cells, materialized via the scalar-left promotion used by MrRectangleBuffer's frameDimensions
         const PointType regionDimensions = sd.refinementFactor * sd.grid.getDimensions();
         nlohmann::json subdomain = nlohmann::json::object();
         subdomain[ "index" ] = sd.index;
         subdomain[ "origin_global_coordinates" ] = vectorToInitReportArray( sd.originCoordinates );
         subdomain[ "region_dimensions" ] = vectorToInitReportArray( regionDimensions );
         report[ "subdomains" ].push_back( std::move( subdomain ) );
      }
      return report;
   }

protected:

   /**
    * \brief Serializes a space-dimension-sized static vector as a
    *        three-component JSON array.
    *
    * The z component is emitted as 0 in 2D, so the reported arrays always
    * have three components - the same convention as the \e domainOrigin /
    * \e domainSize serialization in InitReport.hpp. When the input vector
    * has been scaled to doubles (e.g. the \e region_dimensions footprint),
    * the z component is double 0 as well.
    */
   template< typename Vector >
   static nlohmann::json
   vectorToInitReportArray( const Vector& v )
   {
      if constexpr( Vector::getSize() > 2 )
         return nlohmann::json{ v[ 0 ], v[ 1 ], v[ 2 ] };
      else
         return nlohmann::json{ v[ 0 ], v[ 1 ], typename Vector::RealType( 0 ) };
   }

private:

   std::vector< Subdomain > subdomains;
   std::vector< Interface > interfaces;
   GridType globalGrid;
   int numberOfOverlapLayers = 1;
};

