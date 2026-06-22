#pragma once

#include <array>
#include <cmath>
#include <algorithm>

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>

#include "../../../SPHTraits.h"

namespace TNL {
namespace SPH {
namespace IntegrationSchemes {

namespace AndersonDetail {

// Tiny dense solver for the (regularized) normal equations of the AA
// least-squares problem. n is always small (history depth, typically 2-6),
// so plain Gaussian elimination with partial pivoting is more than fast
// enough - this costs nothing next to a single neighbor search + interact().
template< typename Real, int MaxHistory >
void
solveSmallSystem( int n,
                   std::array< std::array< Real, MaxHistory >, MaxHistory >& A,
                   std::array< Real, MaxHistory >& b,
                   std::array< Real, MaxHistory >& x )
{
   for( int col = 0; col < n; col++ ) {
      int piv = col;
      Real maxAbs = std::abs( A[ col ][ col ] );
      for( int row = col + 1; row < n; row++ )
         if( std::abs( A[ row ][ col ] ) > maxAbs ) {
            maxAbs = std::abs( A[ row ][ col ] );
            piv = row;
         }
      if( piv != col ) {
         std::swap( A[ piv ], A[ col ] );
         std::swap( b[ piv ], b[ col ] );
      }
      const Real diag = A[ col ][ col ];
      if( std::abs( diag ) < Real( 1e-30 ) )
         continue;  // leave x[col] = 0 via back-substitution below
      for( int row = col + 1; row < n; row++ ) {
         const Real factor = A[ row ][ col ] / diag;
         for( int k = col; k < n; k++ )
            A[ row ][ k ] -= factor * A[ col ][ k ];
         b[ row ] -= factor * b[ col ];
      }
   }
   for( int row = n - 1; row >= 0; row-- ) {
      Real s = b[ row ];
      for( int k = row + 1; k < n; k++ )
         s -= A[ row ][ k ] * x[ k ];
      x[ row ] = ( std::abs( A[ row ][ row ] ) < Real( 1e-30 ) ) ? Real( 0 ) : s / A[ row ][ row ];
   }
}

}  // namespace AndersonDetail


/**
 * Anderson acceleration (Type-II / Walker & Ni 2011, Fang & Saad 2009) of the
 * midpoint sub-iteration fixed point
 *
 *    x_{k+1} = g( x_k ),     x = ( a, drho )
 *
 * It is a drop-in replacement for the convex-combination relax() step:
 * it consumes the same x_in = (dvdt_in, drhodt_in) and g(x_in) = (a, drho)
 * (the field just produced by interact()) that relax() already uses, but
 * extrapolates from a short history of past iterates/residuals instead of
 * mixing only the last two. No Jacobian is needed - interact() is treated
 * as a black box, exactly like in the Picard scheme it replaces.
 *
 * Reference: x_{k+1} = x_k + beta*f_k - (DX_k + beta*DF_k) * gamma, where
 * f_k = g(x_k) - x_k, DX_k/DF_k hold the last mk increments of x/f, and
 * gamma minimizes || f_k - DF_k gamma ||_2 (solved via the normal equations,
 * Tikhonov-regularized for robustness).
 */
template< typename SPHConfig, int MaxHistory = 6 >
class AndersonAccelerator
{
public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
   using ScalarArrayType = typename SPHTraitsType::ScalarArrayType;
   using VectorArrayType = typename SPHTraitsType::VectorArrayType;

   void
   setSize( const GlobalIndexType& n )
   {
      size = n;
      for( int i = 0; i <= MaxHistory; i++ ) {
         X_a[ i ].setSize( n );
         G_a[ i ].setSize( n );
         X_rho[ i ].setSize( n );
         G_rho[ i ].setSize( n );
      }
      for( int i = 0; i < MaxHistory; i++ ) {
         dXa[ i ].setSize( n );
         dFa[ i ].setSize( n );
         dXrho[ i ].setSize( n );
         dFrho[ i ].setSize( n );
      }
      Fk_a.setSize( n );
      Fk_rho.setSize( n );
      dotScratch.setSize( n );
   }

   // m            - history length actually used (clamped to MaxHistory; 3-5 is a good start)
   // beta         - damping applied to the residual itself (1.0 = none)
   // reg          - relative Tikhonov regularization of the small Gram matrix
   // gammaLimit   - safeguard: if sum|gamma_i| exceeds this, fall back to a
   //                plain damped Picard step and restart the history
   void
   configure( int m, RealType beta_, RealType reg_ = RealType( 1e-8 ), RealType gammaLimit_ = RealType( 1e2 ) )
   {
      depth = std::min( m, MaxHistory );
      beta = beta_;
      regularization = reg_;
      gammaSafeguard = gammaLimit_;
   }

   // Call once per outer time step, before entering the sub-iteration loop
   // (e.g. at the end of predictor()). Clears the history and freezes the
   // reference scales used to combine the (a, drho) fields into a single
   // inner product for the duration of this step's sub-iterations.
   //
   // Lazily (re)allocates the history buffers if the particle count changed
   // (including the very first call) - this means no separate setSize()
   // call has to be wired into the simulation's initialization code, and
   // particle-count changes from open/periodic boundaries are handled
   // automatically.
   template< typename ParticlesPointer >
   void
   reset( ParticlesPointer& particles, const VectorArrayType& a0, const ScalarArrayType& drho0 )
   {
      if( size != a0.getSize() )
         setSize( a0.getSize() );
      count = 0;
      scaleA2 = meanSquareVec( particles, a0 ) + RealType( 1e-12 );
      scaleRho2 = meanSquareScalar( drho0 ) + RealType( 1e-12 );
   }

   // x_in       = (a_in, rho_in): the iterate that produced the interaction below
   // a_io/rho_io: on input = g(x_in) (freshly computed by interact()),
   //              on output = the accelerated next iterate (written in place,
   //              same convention as MidpointScheme::relax())
   template< typename ParticlesPointer >
   void
   mix( ParticlesPointer& particles,
        const VectorArrayType& a_in, const ScalarArrayType& rho_in,
        VectorArrayType& a_io, ScalarArrayType& rho_io )
   {
      const int slot = count % ( MaxHistory + 1 );
      X_a[ slot ] = a_in;
      X_rho[ slot ] = rho_in;
      G_a[ slot ] = a_io;
      G_rho[ slot ] = rho_io;

      fillDifference( particles, G_a[ slot ], X_a[ slot ], G_rho[ slot ], X_rho[ slot ], Fk_a, Fk_rho );
      std::cout << "mixing - count: " << count << std::endl;

      const int mk = std::min( count, depth );
      if( mk == 0 ) {
         a_io = X_a[ slot ] + beta * Fk_a;
         rho_io = X_rho[ slot ] + beta * Fk_rho;
         count++;
         return;
      }

      // increments dX_j = X_{i+1}-X_i, dF_j = F_{i+1}-F_i over the last mk pairs
      for( int j = 0; j < mk; j++ ) {
         const int i0 = ( count - mk + j ) % ( MaxHistory + 1 );
         const int i1 = ( count - mk + j + 1 ) % ( MaxHistory + 1 );
         fillIncrement( particles,
                         X_a[ i1 ], X_a[ i0 ], G_a[ i1 ], G_a[ i0 ],
                         X_rho[ i1 ], X_rho[ i0 ], G_rho[ i1 ], G_rho[ i0 ],
                         dXa[ j ], dFa[ j ], dXrho[ j ], dFrho[ j ] );
      }

      // Gram matrix / rhs of the small least-squares problem:  A gamma = b
      std::array< std::array< RealType, MaxHistory >, MaxHistory > A{};
      std::array< RealType, MaxHistory > b{};
      std::array< RealType, MaxHistory > gamma{};

      RealType traceA = 0;
      for( int j = 0; j < mk; j++ ) {
         for( int l = j; l < mk; l++ ) {
            const RealType val = dotCombined( particles, dFa[ j ], dFa[ l ], dFrho[ j ], dFrho[ l ] );
            A[ j ][ l ] = val;
            A[ l ][ j ] = val;
         }
         traceA += A[ j ][ j ];
         b[ j ] = dotCombined( particles, dFa[ j ], Fk_a, dFrho[ j ], Fk_rho );
      }
      const RealType reg = regularization * ( traceA / mk + RealType( 1e-30 ) );
      for( int j = 0; j < mk; j++ )
         A[ j ][ j ] += reg;

      AndersonDetail::solveSmallSystem< RealType, MaxHistory >( mk, A, b, gamma );

      RealType gammaNorm = 0;
      for( int j = 0; j < mk; j++ )
         gammaNorm += std::abs( gamma[ j ] );

      if( ! std::isfinite( gammaNorm ) || gammaNorm > gammaSafeguard ) {
         // history looks unreliable (near-collinear increments, blow-up, ...):
         // fall back to a plain damped Picard step and restart the history
         a_io = X_a[ slot ] + beta * Fk_a;
         rho_io = X_rho[ slot ] + beta * Fk_rho;
         count = 0;
         return;
      }

      a_io = X_a[ slot ] + beta * Fk_a;
      rho_io = X_rho[ slot ] + beta * Fk_rho;
      for( int j = 0; j < mk; j++ ) {
         a_io = a_io - gamma[ j ] * ( dXa[ j ] + beta * dFa[ j ] );
         rho_io = rho_io - gamma[ j ] * ( dXrho[ j ] + beta * dFrho[ j ] );
      }

      count++;
   }

public:
   template< typename ParticlesPointer >
   RealType
   meanSquareVec( ParticlesPointer& particles, const VectorArrayType& v )
   {
      auto v_view = v.getConstView();
      auto out_view = dotScratch.getView();
      auto kernel = [=] __cuda_callable__ ( int i ) mutable
      {
         out_view[ i ] = ( v_view[ i ], v_view[ i ] );  // per-particle dot, same idiom as midpointResiduals()
      };
      particles->forAll( kernel );
      return TNL::sum( dotScratch ) / RealType( v.getSize() );
   }

   RealType
   meanSquareScalar( const ScalarArrayType& s )
   {
      return TNL::sum( s * s ) / RealType( s.getSize() );
   }

   // Fa = Ga - Xa, Frho = Grho - Xrho  (the residual f = g(x) - x of one iterate)
   template< typename ParticlesPointer >
   void
   fillDifference( ParticlesPointer& particles,
                    const VectorArrayType& Ga, const VectorArrayType& Xa,
                    const ScalarArrayType& Grho, const ScalarArrayType& Xrho,
                    VectorArrayType& Fa, ScalarArrayType& Frho )
   {
      auto Ga_v = Ga.getConstView();
      auto Xa_v = Xa.getConstView();
      auto Grho_v = Grho.getConstView();
      auto Xrho_v = Xrho.getConstView();
      auto Fa_v = Fa.getView();
      auto Frho_v = Frho.getView();
      auto kernel = [=] __cuda_callable__ ( int i ) mutable
      {
         Fa_v[ i ] = Ga_v[ i ] - Xa_v[ i ];
         Frho_v[ i ] = Grho_v[ i ] - Xrho_v[ i ];
      };
      particles->forAll( kernel );
   }

   // dXa_j = Xa1 - Xa0, dFa_j = (Ga1-Xa1) - (Ga0-Xa0)  -- and the scalar analogues
   template< typename ParticlesPointer >
   void
   fillIncrement( ParticlesPointer& particles,
                  const VectorArrayType& Xa1, const VectorArrayType& Xa0,
                  const VectorArrayType& Ga1, const VectorArrayType& Ga0,
                  const ScalarArrayType& Xrho1, const ScalarArrayType& Xrho0,
                  const ScalarArrayType& Grho1, const ScalarArrayType& Grho0,
                  VectorArrayType& dXa_j, VectorArrayType& dFa_j,
                  ScalarArrayType& dXrho_j, ScalarArrayType& dFrho_j )
   {
      auto Xa1_v = Xa1.getConstView();
      auto Xa0_v = Xa0.getConstView();
      auto Ga1_v = Ga1.getConstView();
      auto Ga0_v = Ga0.getConstView();
      auto Xrho1_v = Xrho1.getConstView();
      auto Xrho0_v = Xrho0.getConstView();
      auto Grho1_v = Grho1.getConstView();
      auto Grho0_v = Grho0.getConstView();
      auto dXa_v = dXa_j.getView();
      auto dFa_v = dFa_j.getView();
      auto dXrho_v = dXrho_j.getView();
      auto dFrho_v = dFrho_j.getView();
      auto kernel = [=] __cuda_callable__ ( int i ) mutable
      {
         dXa_v[ i ] = Xa1_v[ i ] - Xa0_v[ i ];
         dFa_v[ i ] = ( Ga1_v[ i ] - Xa1_v[ i ] ) - ( Ga0_v[ i ] - Xa0_v[ i ] );
         dXrho_v[ i ] = Xrho1_v[ i ] - Xrho0_v[ i ];
         dFrho_v[ i ] = ( Grho1_v[ i ] - Xrho1_v[ i ] ) - ( Grho0_v[ i ] - Xrho0_v[ i ] );
      };
      particles->forAll( kernel );
   }

   // Combined, non-dimensionalized inner product over the (a, drho) state,
   // using the same per-particle dot idiom as MidpointScheme::midpointResiduals().
   template< typename ParticlesPointer >
   RealType
   dotCombined( ParticlesPointer& particles,
                const VectorArrayType& va, const VectorArrayType& vb,
                const ScalarArrayType& sa, const ScalarArrayType& sb )
   {
      auto va_v = va.getConstView();
      auto vb_v = vb.getConstView();
      auto sa_v = sa.getConstView();
      auto sb_v = sb.getConstView();
      auto out_v = dotScratch.getView();
      const RealType invA = RealType( 1 ) / scaleA2;
      const RealType invRho = RealType( 1 ) / scaleRho2;
      auto kernel = [=] __cuda_callable__ ( int i ) mutable
      {
         out_v[ i ] = invA * ( va_v[ i ], vb_v[ i ] ) + invRho * sa_v[ i ] * sb_v[ i ];
      };
      particles->forAll( kernel );
      return TNL::sum( dotScratch );
   }

   GlobalIndexType size = 0;
   int depth = 3;
   int count = 0;
   RealType beta = 1.0;
   RealType regularization = 1e-8;
   RealType gammaSafeguard = 1e2;
   RealType scaleA2 = 1, scaleRho2 = 1;

   std::array< VectorArrayType, MaxHistory + 1 > X_a, G_a;
   std::array< ScalarArrayType, MaxHistory + 1 > X_rho, G_rho;
   std::array< VectorArrayType, MaxHistory > dXa, dFa;
   std::array< ScalarArrayType, MaxHistory > dXrho, dFrho;
   VectorArrayType Fk_a;
   ScalarArrayType Fk_rho;
   ScalarArrayType dotScratch;
};

}  // namespace IntegrationSchemes
}  // namespace SPH
}  // namespace TNL

