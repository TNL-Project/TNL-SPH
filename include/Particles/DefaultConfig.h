#include <type_traits>
#include <Particles/CellIndexer.h>
#include <Particles/ParticlesTraits.h>

namespace TNL {
namespace ParticleSystem {

template< int SpaceDimension, typename Device, typename Real = float, typename GlobalIndex = int, typename LocalIndex = int >
class DefaultConfig
{
   public:
   using DeviceType = Device;

   using GlobalIndexType = int;
   using LocalIndexType = int;
   using CellIndexType = int;
   using RealType = float;

   static constexpr int spaceDimension = SpaceDimension;

   using UseWithDomainDecomposition = std::false_type;
   using CoordinatesType = Containers::StaticVector< spaceDimension, int >;
   using CellIndexerSequence = std::conditional_t< SpaceDimension == 2,
                                                   std::index_sequence< 0, 1 >,
                                                   std::index_sequence< 0, 1, 2 > >;
   using CellIndexerType = SimpleCellIndex< spaceDimension, DefaultConfig, CellIndexerSequence >;
   using NeighborListType = typename Algorithms::Segments::Ellpack< DeviceType, int >;
};

} //namespace Particles
} //namespace TNL

