#pragma once

#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Pointers/SharedPointer.h>
#include <TNL/Algorithms/sort.h>

#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>
#include <type_traits>
#include "../SPH/shared/thrustExecPolicySelector.h"

#include "ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {

template < typename ParticleConfig, typename DeviceType >
class Particles
{
public:
   using Device = DeviceType;
   using Config = ParticleConfig;
   using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;

   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using LocalIndexType = typename ParticleTraitsType::LocalIndexType;
   using RealType = typename ParticleTraitsType::RealType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using PointType = typename ParticleTraitsType::PointType;

   using IndexArrayType = typename ParticleTraitsType::CellIndexArrayType;
   using PointArrayType = typename ParticleTraitsType::PointArrayType;
   using IndexArrayTypePointer = typename Pointers::SharedPointer< IndexArrayType, Device >; //TODO: Do I need this?

   static constexpr int spaceDimension = Config::spaceDimension;

   static std::string
   writeModelType()
   {
      return "TNL::ParticleSystem::Particles";
   }

   /**
    * Constructors.
    */
   Particles() : points( 0 ), points_swap( 0 ), sortPermutations( 0 ) {}

   Particles( GlobalIndexType size, GlobalIndexType sizeAllocated )
   : numberOfParticles( size ), numberOfAllocatedParticles( sizeAllocated ), points( sizeAllocated ) { }

   Particles( GlobalIndexType size, GlobalIndexType sizeAllocated, RealType radius )
   : numberOfParticles( size ),
     numberOfAllocatedParticles( sizeAllocated ),
     points( sizeAllocated ),
     points_swap( sizeAllocated ),
     sortPermutations( sizeAllocated ),
     mark( sizeAllocated ),
     radius( radius ) { }

   /**
    * Get dimension of particle system.
    */
   [[nodiscard]]
   static constexpr int
   getParticlesDimension();

   void
   setSize( const GlobalIndexType& size );

   /**
    * Get number of particles in particle system.
    */
   [[nodiscard]] __cuda_callable__
   const GlobalIndexType
   getNumberOfParticles() const;

   [[nodiscard]] __cuda_callable__
   const GlobalIndexType
   getNumberOfAllocatedParticles() const;

   /**
    * \brief Set number of active particles in particle system.
    *
    * \param new number of active particles in the system
    */
   void
   setNumberOfParticles( const GlobalIndexType& newNumberOfParticles );

   /**
    * \brief Get value of the search (cutoff) radius.
    *
    * This method can be called from device kernels.
    */
   [[nodiscard]] __cuda_callable__
   const RealType
   getSearchRadius() const;

   /**
    * \brief Set value of the search (cutoff) radius.
    *
    * \param new value of the search radius
    */
   void
   setSearchRadius( const RealType& searchRadius );

   /**
    * Grid referential origin.
    */
   [[nodiscard]] const PointType
   getGridReferentialOrigin() const;

   void
   setGridReferentialOrigin( const PointType& origin );

   /**
    *
    */
   [[nodiscard]] const IndexVectorType
   getGridOriginGlobalCoords() const;

   void
   setGridOriginGlobalCoords( const IndexVectorType& origin );

   [[nodiscard]] const IndexVectorType
   getGridOriginGlobalCoordsWithOverlap() const;

   /**
    * \brief Returns origin of the implicit linked list grid.
    *
    * \return the origin of the grid.
    */
   [[nodiscard]] const PointType
   getGridOrigin() const;

   /**
    * \brief Set origin of the implicit linked list grid.
    *
    * \param the origin of the grid.
    */
   void
   setGridOrigin( const PointType& origin );

   /**
    * \brief Returns origin of the implicit linked list grid including overlap.
    *
    * \return the origin of the grid.
    */
   [[nodiscard]] const PointType
   getGridOriginWithOverlap() const;

   /**
    * \brief Returns dimensions of the implicit linked list grid.
    *
    * \return Coordinate vector with number of edges along each axis.
    */
   [[nodiscard]] const IndexVectorType
   getGridDimensions() const;

   /**
    * \brief Set dimensions of the implicit linked list grid.
    *
    * \param gridSize grid dimensions given in a form of coordinate vector.
    */
   virtual void
   setGridDimensions( const IndexVectorType& dimensions );

   /**
    * \brief Returns dimensions of the implicit linked list grid including overlap.
    *
    * \return Coordinate vector with number of edges along each axis.
    */
   [[nodiscard]] const IndexVectorType
   getGridDimensionsWithOverlap() const;

   /**
    * \brief Set width of overlap expressed in number of cells.
    *
    * \param Integer expressing the overlap width in number of cells
    */
   virtual void
   setOverlapWidth( const GlobalIndexType width );

   /**
    * \brief Set width of overlap expressed in number of cells.
    *
    * \return Integer expressing the overlap width in number of cells
    */
   [[nodiscard]] const GlobalIndexType
   getOverlapWidth() const;

   /**
    * \brief Return const-qualified array of spatial coordinates of particles.
    */
   [[nodiscard]] const PointArrayType&
   getPoints() const;

   /**
    * \brief Return array of spatial coordinates of particles.
    */
   PointArrayType&
   getPoints();

   /**
    * \brief Return const-qualified array of spatial coordinates of particles used as swap field.
    * Points are defined in two duplicit arrays to avoid some inplace operations.
    */
   [[nodiscard]] const PointArrayType&
   getPointsSwap() const;

   /**
    * \brief Return array of spatial coordinates of particles used as swap field.
    * Points are defined in two duplicit arrays to avoid some inplace operations.
    */
   PointArrayType&
   getPointsSwap();

   /**
    * \brief Get list with new particle ordering. This is used for reodering of fields with variables
    * defined on particles.
    *
    * \return pointer to array with permutations.
    */
   [[nodiscard]] const IndexArrayTypePointer&
   getSortPermutations() const;

   /**
    * \brief Get list with new particle ordering. This is used for reodering of fields with variables
    * defined on particles.
    *
    * \return pointer to array with permutations.
    */
   IndexArrayTypePointer&
   getSortPermutations();

   /**
    * \brief Get number of particles that are set to be removed during next search procedure.
    *
    * \return number of particles
    */
   [[nodiscard]] const GlobalIndexType
   getNumberOfParticlesToRemove() const;

   /**
    * \brief Set number of particles to remove. Since we can remove particles from any position in the set,
    * in order to remove particles, their position should be invalidated (set as FLT_MAX) and corresponding
    * number of particles has to be saved. Particles are erased during next search procedure.
    *
    * \param number of invalidated particles that should be erased during next search procedure.
    */
   void
   setNumberOfParticlesToRemove( const GlobalIndexType& removeCount );

   /**
    * \brief Execute function \e f in parallel for all particles
    *
    * The function \e f is executed as `f(i)`, where `GlobalIndexType i` is the global index of the
    * particle to be processed. The particle set itself is not passed to the function `f`, it is the user's
    * responsibility to ensure proper access to the mesh if needed, e.g. by the means of lambda capture
    * and/or using a \ref TNL::Pointers::SharedPointer "SharedPointer".
    */
   template< typename Device2 = DeviceType,
             typename Func,
             typename UseWithDomainDecomposition = typename Config::UseWithDomainDecomposition,
             std::enable_if_t< !UseWithDomainDecomposition::value, bool > Enabled = true >
   void
   forAll( Func f ) const;

   /**
    * \brief Execute function \e f in parallel for all particles
    *
    * The function \e f is executed as `f(i)`, where `GlobalIndexType i` is the global index of the
    * particle to be processed. The particle set itself is not passed to the function `f`, it is the user's
    * responsibility to ensure proper access to the mesh if needed, e.g. by the means of lambda capture
    * and/or using a \ref TNL::Pointers::SharedPointer "SharedPointer".
    */
   template< typename Device2 = DeviceType,
             typename Func,
             typename UseWithDomainDecomposition = typename Config::UseWithDomainDecomposition,
             std::enable_if_t< UseWithDomainDecomposition::value, bool > Enabled = true >
   void
   forAll( Func f ) const;

   /**
    * \brief Test if is particle with given index located inside domain
    *
    * \param index of particle to test
    *
    * \return True or False based on the presence of particles in given domain.
    */
   __cuda_callable__
   bool
   isInsideDomain( const GlobalIndexType& index );

   __cuda_callable__
   bool
   isInsideDomain( const PointType& point, const PointType& domainOrigin, const PointType& domainSize ) const;

   __cuda_callable__
   bool
   isInsideDomain( const IndexVectorType& particleCell, const IndexVectorType& gridDimensionWithOverlap ) const;

   /**
    * \brief Reorder particles based on permutation vector.
    */
   void
   reorderParticles();

   void
   writeProlog( TNL::Logger& logger ) const noexcept;

protected:

   /**
    * Number of active paticles inside the domain.
    */
   GlobalIndexType numberOfParticles;

   /**
    * Total number of allocated particles inside domain. In case we assume
    * that number of particles is not constant, we pre-allocate the fields with
    * bigger size.
    */
   GlobalIndexType numberOfAllocatedParticles;

   /**
    * Number of particles to remove. Remove procedure is not trivial and it is performed
    * during the search procedure.
    */
   GlobalIndexType numberOfParticlesToRemove = 0;

   /**
    * Search radius (cut off raidus) for nearest neighbors search.
    */
   RealType radius;

   /**
    * Origin of the global domain where particles live.
    */
   bool gridReferentialOriginSet = false;
   PointType gridReferentialOrigin;

   /**
    * In case that referential origin is used, store global coordinates of current origin.
    */
   IndexVectorType gridOriginGlobalCoords;

   /**
    * Origin of the domain where particles live. The domain is assumed in form
    * of uniform cartesian rectangular grid.
    */
   PointType gridOrigin;

   /**
    * Size of the domain where particles live expressed as a grid dimension
    * with cell size corresponding to search radius. The domain s assumed in form of
    * uniform cartezian grid.
    */
   IndexVectorType gridDimension;

   /**
    * Overlap width expressed by the number of cells.
    */
   LocalIndexType overlapWidth = 0;

   /**
    * Array contaning the points (particle positions). Size is numberOfAllocatedParticles.
    */
   PointArrayType points;

   /**
    * Duplicit array with points (particle positions). This one is used to avoid
    * in-place sort. During the sort, this array is filled with the particles using new
    * ordering and then swaped with the original points array.
    */
   PointArrayType points_swap;

   /**
    * Array to store permutation of new particle ordering during particle sort.
    */
   IndexArrayTypePointer sortPermutations;

   /**
    * Array marking particles (if allocated).
    */
   IndexArrayType mark;
};

} //namespace Particles
} //namespace TNL

#include "Particles.hpp"

