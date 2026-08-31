#pragma once

#include <cstdio>
#include <fstream>
#include "../SPHTraits.h"

namespace TNL {
namespace SPH {

/**
 * \brief Evaluates torque (moment of forces) acting on a boundary with respect
 * to a given center and rotation axis.
 *
 * The evaluation reuses the per-particle forces stored by a forces monitor
 * (EvaluateForces or EvaluateForces_BoundaryIntegrals), so the forces monitor
 * must be updated by its computeForces() method before calling computeTorque().
 * The marker-based selection of boundary parts is inherited from the forces
 * monitor - particles outside the selection have zero force and hence
 * contribute zero torque.
 *
 * The torque is computed as
 *   T = sum_i ( r_i - center ) x F_i
 * and the moment about the given axis is the projection M = ( axis_n, T ),
 * where axis_n is the normalized axis vector (axis must be non-zero).
 */
template< typename SPHDefs >
class EvaluateTorque
{

public:

   using SPHConfig = typename SPHDefs::SPHConfig;
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using IndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using VectorType = typename SPHTraitsType::VectorType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   EvaluateTorque() = default;

   template< typename BoundaryPointer >
   EvaluateTorque( BoundaryPointer& boundary )
   {
      init( boundary );
   }

   template< typename BoundaryPointer >
   void
   init( BoundaryPointer& boundary )
   {
      const IndexType size = boundary->getParticles()->getNumberOfAllocatedParticles();
      T_pressure.setSize( size );
      T_visco.setSize( size );
      T_pressure = 0.f;
      T_visco = 0.f;
      T_pressure_sum = 0.f;
      T_visco_sum = 0.f;
      M_pressure = 0.f;
      M_visco = 0.f;
   }

   template< typename BoundaryPointer, typename ForcesMonitor >
   void
   computeTorque( BoundaryPointer& boundary,
                  ForcesMonitor& forcesMonitor,
                  const VectorType& center,
                  const VectorType& axis )
   {
      const auto view_points_bound = boundary->getParticles()->getPoints().getConstView();
      const auto view_F_pressure = forcesMonitor.F_pressure.getConstView();
      const auto view_F_visco = forcesMonitor.F_visco.getConstView();

      auto view_T_pressure = T_pressure.getView();
      auto view_T_visco = T_visco.getView();

      auto particleLoop = [=] __cuda_callable__ ( IndexType i ) mutable
      {
         const VectorType r_i = view_points_bound[ i ] - center;
         view_T_pressure[ i ] = VectorProduct( r_i, view_F_pressure[ i ] );
         view_T_visco[ i ] = VectorProduct( r_i, view_F_visco[ i ] );
      };
      boundary->getParticles()->forAll( particleLoop );

      T_pressure_sum = TNL::sum( T_pressure );
      T_visco_sum = TNL::sum( T_visco );

      // moments with respect to the rotation axis (axis must be non-zero)
      const VectorType axis_n = axis / TNL::l2Norm( axis );
      M_pressure = ( axis_n, T_pressure_sum );
      M_visco = ( axis_n, T_visco_sum );
   }

   void
   output( const std::string& outputPath, const int step, const RealType time )
   {
      const VectorType T_total = T_pressure_sum + T_visco_sum;
      const RealType M_total = M_pressure + M_visco;
      std::ofstream outfile;
      outfile.open( outputPath, std::ios_base::app );
      outfile << step << " " << time << " "
              << T_pressure_sum.x() << " " << T_pressure_sum.y() << " " << T_pressure_sum.z() << " "
              << T_visco_sum.x() << " " << T_visco_sum.y() << " " << T_visco_sum.z() << " "
              << T_total.x() << " " << T_total.y() << " " << T_total.z() << " "
              << M_pressure << " " << M_visco << " " << M_total << std::endl;
   }

   VectorArrayType T_pressure;
   VectorArrayType T_visco;

   VectorType T_pressure_sum;
   VectorType T_visco_sum;

   RealType M_pressure;
   RealType M_visco;

};

} // SPH
} // TNL
