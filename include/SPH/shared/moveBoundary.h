#include <cstdio>
#include <cstdlib>
#include "../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace features {

template< typename Boundary >
class Motion
{
public:

   using DeviceType = typename Boundary::DeviceType;
   using BoundaryPointer = Pointers::SharedPointer< Boundary, DeviceType >;

   using TraitsType = TNL::SPH::SPHFluidTraits< typename SPHDefs::SPHConfig >;
   using RealType = typename TraitsType::RealType;
   using VectorType = typename TraitsType::VectorType;
   using MatrixType = typename TraitsType::MatrixType;

   Motion() = default;
   Motion( BoundaryPointer& boundary )
   {
      boundaryRefConfiguration = boundary;
   }

   void
   init( BoundaryPointer& boundary )
   {
      boundaryRefConfiguration = boundary;
   }

   void
   transform( BoundaryPointer& boundary, const VectorType rotVect, const VectorType translVect )
   {

      const RealType sx = sin( rotVect.x() );
      const RealType cx = cos( rotVect.x() );
      const RealType sy = sin( rotVect.y() );
      const RealType cy = cos( rotVect.y() );
      const RealType sz = sin( rotVect.z() );
      const RealType cz = cos( rotVect.z() );
      const MatrixType rotMatrix = { ( cz * cy ), ( cz * sy * sx - sz * cx ), ( cz * sy * cx + sz * sx ),
                                     ( sz * cy ), ( sz * sy * sx + cz * cx ), ( sz * sy * cx - cz * sx ),
                                     ( -sy     ), ( cy * sx                ), ( cy * cx                ) };

      auto r_view = boundary->getPoints().getView();
      auto v_view = boundary->getVariables()->v.getView();
      auto n_view = boundary->getVariables()->n.getView();
      const auto marker_view = boundary->getVariables()->marker.getConstView();

      auto rotate  = [ = ] __cuda_callable__( int i ) mutable
      {
         if( marker_view[ i ] == movingBoundary )
         {
            const VectorType r = r_view[ i ];
            const VectorType n = n_view[ i ];

            const VectorType r_new = rotMatrix * r + translVect;
            r_view[ i ] = r_new;
            const VectorType n_new = rotMatrix * n;
            n_view[ i ] = n_new;
            v_view[ i ] = VectorProduct( rotVect, r );
         }
      };
      boundary->getParticles()->forAll( rotate );
   }

   void
   move( BoundaryPointer& boundary, const VectorType rotVect, const VectorType translVect )
   {
      // synchronize positions, normals and element sizes (the rest can be ignored)
      boundary->getParticles()->reorderArray(
            boundaryRefConfiguration->getParticles()->getPoints(), boundaryRefConfiguration->getParticles()->getPointsSwap() );
      boundaryRefConfiguration->getVariables()->sortVariables( boundary->getParticles() );

      // copy geometridal data from the referential boundary and transfer them to actual boundary
      boundary->getPoints() = boundaryRefConfiguration->getPoints();
      boundary->getVariables()->n = boundaryRefConfiguration->getVariables()->n;
      boundary->getVariables()->elementSize = boundaryRefConfiguration->getVariables()->elementSize;

      // update the live object
      transform( boundary, rotVect, translVect );
   }

protected:

   BoundaryPointer boundaryRefConfiguration;

};

} // features
} // SPH
} // TNL

