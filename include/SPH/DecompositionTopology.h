#pragma once

#include <vector>
#include <stdexcept>
#include <string>

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


private:

   std::vector< Subdomain > subdomains;
   std::vector< Interface > interfaces;
   GridType globalGrid;
   int numberOfOverlapLayers = 1;
};

