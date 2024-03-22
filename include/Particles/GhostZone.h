#pragma once

#include <TNL/Algorithms/scan.h>

#include "ParticlesTraits.h"

namespace TNL {
namespace ParticleSystem {


template< typename ParticleConfig >
class ParticleZone
{
public:

   using DeviceType = typename ParticleConfig::DeviceType;
   using Config = ParticleConfig;
   using ParticleTraitsType = ParticlesTraits< Config, DeviceType >;
   using CellIndexer = typename ParticleConfig::CellIndexerType;

   using GlobalIndexType = typename ParticleTraitsType::GlobalIndexType;
   using PairIndexType = typename ParticleTraitsType::PairIndexType;
   using IndexVectorType = typename ParticleTraitsType::IndexVectorType;
   using IndexArrayType = typename ParticleTraitsType::CellIndexArrayType; //FIXME
   using RealType = typename ParticleTraitsType::RealType;
   using PointType = typename ParticleTraitsType::PointType;


   /**
    * Constructor.
    */
   ParticleZone() = default;

   /**
    * Constructor.
    */
   ParticleZone( GlobalIndexType numberOfParticlesPerCell ) : numberOfParticlesPerCell( numberOfParticlesPerCell ) {}

   /**
    * Constructor - allocate only the field.
    */
   ParticleZone( GlobalIndexType numberOfCells, GlobalIndexType numerOfParticlesPerCells )
   : numberOfCellsInZone( numberOfCells ),
     cellsInZone( numberOfCells ),
     numberOfParticlesInCell( numberOfCells ),
     particlesInZone( numberOfCells * numerOfParticlesPerCells ) {}


   /**
    * Assign cells from point and direction for grid-base orthogonal zones
    */
   void
   assignCells( IndexVectorType startingPoint, IndexVectorType size, IndexVectorType gridSize );

   /**
    * Assign cells from point and direction for grid-base orthogonal zones
    */
   void
   assignCells( const PointType firstPoint,
                const PointType secondPoint,
                IndexVectorType gridSize,
                PointType gridOrigin,
                RealType searchRadius );

   /**
    * Assign cells from another array.
    */
   template< typename Array >
   void
   assignCells( Array& inputCells );

   /**
    * Get indices of cells contained in the zone.
    */
   const IndexArrayType&
   getCellsInZone() const;

   IndexArrayType&
   getCellsInZone();

   /**
    * Get indices of particles contained in the zone..
    */
   const IndexArrayType&
   getParticlesInZone() const;

   IndexArrayType&
   getParticlesInZone();

   /**
    * Get number of particles in zone.
    */
   const GlobalIndexType
   getNumberOfParticles() const;

   /**
    * Get number of cells in zone.
    */
   const GlobalIndexType
   getNumberOfCells() const;

   ///**
   // * Get number of cells in zone.
   // */
   //const GlobalIndexType
   //getMaxNumberOfParticlesPerCell() const;

   void
   setNumberOfParticlesPerCell( const GlobalIndexType numberOfParticlesPerCell );

   void
   resetParticles();

   void
   resetZoneCells();

   /**
    * Collect particles contained in the zone
    */
   template< typename ParticlesPointer >
   void
   collectNumbersOfParticlesInCells( const ParticlesPointer& particles );

   /**
    * Collect particles contained in the zone
    */
   template< typename ParticlesPointer >
   void
   buildParticleList( const ParticlesPointer& particles );

   /**
    * Collect particles contained in the zone
    */
   template< typename ParticlesPointer >
   void updateParticlesInZone( const ParticlesPointer& particles );

   /**
    * Collect particles contained in the zone
    */
   template< typename ParticlesPointer, typename TimeMeasurement >
   void updateParticlesInZone( const ParticlesPointer& particles, TimeMeasurement& timeMeasurement );

   void
   writeProlog( TNL::Logger& logger ) const noexcept;


protected:

   GlobalIndexType numberOfParticlesPerCell;

   GlobalIndexType numberOfCellsInZone;

   IndexArrayType cellsInZone;

   IndexArrayType numberOfParticlesInCell; //FIXME: Rename.

   GlobalIndexType numberOfParticlesInZone = 0;
   IndexArrayType particlesInZone;


};

} // Particles
} // TNL

#include "GhostZone.hpp"

