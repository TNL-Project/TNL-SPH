#include <vector>
#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

template< typename SPHDefs >
class MassConservationMonitor
{
public:
   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   template< typename FluidSetsVector, typename ModelParams >
   void
   init( const FluidSetsVector& fluidSets, const ModelParams& modelParams )
   {
      numberOfSubsets = fluidSets.size();
      particleMassPerSet.resize( numberOfSubsets );
      massPerSet.resize( numberOfSubsets );

      for( int i = 0; i < numberOfSubsets; i++ ){
            const RealType rcut = fluidSets[ i ]->getParticles()->getSearchRadius();
            const RealType refinementFactor = rcut / ( 2.f * modelParams.h );
            const unsigned int dim = SPHConfig::spaceDimension;
            const RealType m = std::pow( refinementFactor, dim ) * modelParams.mass;
            particleMassPerSet[ i ] = m;
      }
   }

   template< typename FluidSetsVector >
   void
   sumTotalMass( const FluidSetsVector& fluidSets )
   {
      totalMass = 0.;

      for( int i = 0; i < numberOfSubsets; i++ ){
         const RealType m = particleMassPerSet[ i ];
         const IndexType np = fluidSets[ i ]->getNumberOfParticles();
         massPerSet[ i ] = m * np;
         totalMass += m * np;
      }
   }

   void
   output( const std::string& outputPath, const int step, const RealType time )
   {
      std::ofstream outfile;
      outfile.open(outputPath, std::ios_base::app );
      outfile << step << " " << time << " " << totalMass;
      for( int i = 0; i < numberOfSubsets; i++ )
         outfile <<  " " << massPerSet[ i ];
      outfile << std::endl;
   }

protected:

   int numberOfSubsets;
   std::vector< RealType > particleMassPerSet;
   std::vector< RealType > massPerSet;
   RealType totalMass;
};

} // SPH
} // TNL

