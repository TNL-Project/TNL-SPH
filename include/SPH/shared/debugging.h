#pragma once

#include <cuda_runtime.h>
#include <TNL/Algorithms/reduce.h>
#include <exception>
#include <iostream>
#include <string>

// Debug macros — centralized here so all translation units share one definition.
// Use #ifndef guards so files that defined their own copy still compile.

#ifndef SPH_CUDA_CHECK
   #define SPH_CUDA_CHECK()                                                                  \
      do {                                                                                   \
         cudaDeviceSynchronize();                                                            \
         const cudaError_t _sph_err = cudaGetLastError();                                    \
         if( _sph_err != cudaSuccess ) {                                                     \
            std::cerr << "CUDA error in " << __func__ << " (" << __FILE__ << ":" << __LINE__ \
                      << "): " << cudaGetErrorString( _sph_err ) << std::endl;               \
            std::abort();                                                                    \
         }                                                                                   \
      } while( 0 )
#endif

// Wrap a call with try/catch + CUDA sync check. Aborts on failure.
#ifndef SPH_DEBUG_CALL
   #define SPH_DEBUG_CALL( call, name )                                                                          \
      do {                                                                                                       \
         try {                                                                                                   \
            call;                                                                                                \
         }                                                                                                       \
         catch( const std::exception& _sph_e ) {                                                                 \
            std::cerr << "FAILED in " name " (" __FILE__ ":" << __LINE__ << "): " << _sph_e.what() << std::endl; \
            std::abort();                                                                                        \
         }                                                                                                       \
         SPH_CUDA_CHECK();                                                                                       \
      } while( 0 )
#endif

// Same as SPH_DEBUG_CALL, but runs a debug action (e.g. inspect particles)
// in the catch block before aborting. Use this to diagnose the failure state.
//
// Example:
//   SPH_DEBUG_CALL_D(
//      fluid_own->searchForNeighbors(),
//      "fluid_own->searchForNeighbors",
//      Debugging::checkParticlesOutDomain( fluid_own, "fluid_own" ) );
//
#ifndef SPH_DEBUG_CALL_D
   #define SPH_DEBUG_CALL_D( call, name, debug_action )                                                          \
      do {                                                                                                       \
         try {                                                                                                   \
            call;                                                                                                \
         }                                                                                                       \
         catch( const std::exception& _sph_e ) {                                                                 \
            std::cerr << "FAILED in " name " (" __FILE__ ":" << __LINE__ << "): " << _sph_e.what() << std::endl; \
            debug_action;                                                                                        \
            std::abort();                                                                                        \
         }                                                                                                       \
         SPH_CUDA_CHECK();                                                                                       \
      } while( 0 )
#endif

// Debug functions
namespace TNL {
namespace SPH {
namespace Debugging {

// Inspect a particle set for out-of-domain or NaN/Inf particles.
// Prints per-particle details (index, position, cell coords) and a summary.
//
// Works on any type that exposes getPoints(), getNumberOfParticles(),
// getNumberOfAllocatedParticles(), and getParticles() — i.e. ParticleSet
// and its subclasses (FluidPointer, MultiresolutionBoundaryPointer, etc.).
//
// Pass the dereferenced object: checkParticlesOutDomain( *fluid_ptr, "fluid" )
// or from inside a method: checkParticlesOutDomain( *this, "buffer" )
//
// isBuffer: labels the output as "(buffer)" — purely cosmetic.
template< typename ParticleSetType >
void
checkParticlesOutDomain( const ParticleSetType& particles, const std::string& name, bool isBuffer = false )
{
   using DeviceType = typename ParticleSetType::DeviceType;
   using IndexVectorType = typename ParticleSetType::IndexVectorType;
   using VectorType = typename ParticleSetType::VectorType;
   using RealType = typename ParticleSetType::RealType;

   const auto r_view = particles.getPoints().getConstView();
   const int n = particles.getNumberOfParticles();
   const int n_alloc = particles.getNumberOfAllocatedParticles();
   const int n_remove = particles.getParticles()->getNumberOfParticlesToRemove();

   const VectorType gridRefOrigin = particles.getParticles()->getGridReferentialOrigin();
   const RealType searchRadius = particles.getParticles()->getSearchRadius();
   const RealType invSearchRadius = 1.0f / searchRadius;
   const IndexVectorType gridOriginGlob = particles.getParticles()->getGridOriginGlobalCoordsWithOverlap();
   const IndexVectorType gridDims = particles.getParticles()->getGridDimensionsWithOverlap();

   const int dim = IndexVectorType::getSize();

   // --- Header ---
   std::cout << "\n=== checkParticlesOutDomain ===" << std::endl;
   std::cout << "  set:         " << name << ( isBuffer ? " (buffer)" : "" ) << std::endl;
   std::cout << "  particles:   " << n << " / " << n_alloc << " allocated"
             << " (toRemove: " << n_remove << ")" << std::endl;
   std::cout << "  searchRadius: " << searchRadius << std::endl;
   std::cout << "  gridRefOrigin: (";
   for( int d = 0; d < dim; d++ )
      std::cout << gridRefOrigin[ d ] << ( d < dim - 1 ? ", " : "" );
   std::cout << ")" << std::endl;
   std::cout << "  gridOriginGlob: (";
   for( int d = 0; d < dim; d++ )
      std::cout << gridOriginGlob[ d ] << ( d < dim - 1 ? ", " : "" );
   std::cout << ")" << std::endl;
   std::cout << "  gridDims:    (";
   for( int d = 0; d < dim; d++ )
      std::cout << gridDims[ d ] << ( d < dim - 1 ? ", " : "" );
   std::cout << ")" << std::endl;

   // --- Device kernel: check each particle ---
   auto check = [ = ] __cuda_callable__( int i ) -> int
   {
      const VectorType r = r_view[ i ];

      // Skip already-removed particles (marked with FLT_MAX)
      if( r[ 0 ] == FLT_MAX )
         return 0;

      // Check for NaN (NaN != NaN is the only reliable test on GPU)
      for( int d = 0; d < dim; d++ ) {
         if( r[ d ] != r[ d ] ) {
            printf( "  [NaN] i=%d pos=(%f,%f,%f)\n", i, r[ 0 ], r[ 1 ], r[ 2 ] );
            return 1;
         }
      }

      // Check for Inf (anything beyond a sane threshold)
      for( int d = 0; d < dim; d++ ) {
         if( r[ d ] > 1e30f || r[ d ] < -1e30f ) {
            printf( "  [Inf] i=%d pos=(%f,%f,%f)\n", i, r[ 0 ], r[ 1 ], r[ 2 ] );
            return 1;
         }
      }

      // Check domain bounds (extended domain including overlap)
      const IndexVectorType cg = TNL::floor( ( r - gridRefOrigin ) * invSearchRadius );
      const IndexVectorType cc = cg - gridOriginGlob;
      for( int d = 0; d < dim; d++ ) {
         if( cc[ d ] < 0 || cc[ d ] >= gridDims[ d ] ) {
            printf( "  [OUT] i=%d pos=(%.3f,%.3f,%.3f) cell=(%d,%d,%d) local=(%d,%d,%d) dims=(%d,%d,%d)\n",
                    i,
                    r[ 0 ],
                    r[ 1 ],
                    r[ 2 ],
                    cg[ 0 ],
                    cg[ 1 ],
                    cg[ 2 ],
                    cc[ 0 ],
                    cc[ 1 ],
                    cc[ 2 ],
                    gridDims[ 0 ],
                    gridDims[ 1 ],
                    gridDims[ 2 ] );
            return 1;
         }
      }
      return 0;
   };

   const int bad = TNL::Algorithms::reduce< DeviceType >( 0, n, check, TNL::Plus() );
   cudaDeviceSynchronize();

   // --- Summary ---
   if( bad == 0 )
      std::cout << "  result:      OK — all " << n << " particles inside domain" << std::endl;
   else
      std::cout << "  result:      " << bad << " bad particle(s) found" << std::endl;
   std::cout << "=== end ===" << std::endl;
}

}  // namespace Debugging
}  // namespace SPH
}  // namespace TNL

