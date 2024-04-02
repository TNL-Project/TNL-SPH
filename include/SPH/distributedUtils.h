namespace TNL {
namespace SPH {
namespace distributed {

Containers::StaticVector< 2, int >
restoreSubdomainCoordinatesFromRank( int rank, Containers::StaticVector< 2, int > numberOfSubdomains )
{
   Containers::StaticVector< 2, int > subdomainCoordinates = 0.;
   subdomainCoordinates[ 1 ] = std::floor( rank / numberOfSubdomains[ 0 ] );
   subdomainCoordinates[ 0 ] = rank - std::floor( rank / numberOfSubdomains[ 0 ] );
   return subdomainCoordinates;
}

std::string
getSubdomainKey( int rank, Containers::StaticVector< 2, int > numberOfSubdomains )
{
   Containers::StaticVector< 2, int > subdomainCoordinates = distributed::restoreSubdomainCoordinatesFromRank(
         rank, numberOfSubdomains );
   //NOTE: Only 2D decomposition is allowed
   std::string subdomainKey = "subdomain-x-" + std::to_string( subdomainCoordinates[ 0 ] ) +
                              "-y-" + std::to_string( subdomainCoordinates[ 1 ] ) + "-";
   return subdomainKey;
}

}  //namespace distributed
}  //namespace SPH
}  //namespace TNL

