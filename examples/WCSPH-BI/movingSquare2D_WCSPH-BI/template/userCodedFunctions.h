#include <SPH/SPHTraits.h>
#include <climits>

namespace userCodedFunctions {

void
userConfigSetup( TNL::Config::ConfigDescription& config )
{
   config.addEntry< int >( "filtering-steps-interval", "The initial number of fluid particles.", INT_MAX );
}

class CustomMotion
{
   public:
   static constexpr int movingBoundary = 1;

   // load file with prescribed motion
   CustomMotion( const std::string motionFileName )
   {
      std::ifstream file( motionFileName );
      if ( !file.is_open() ) {
         std::cerr << "Cannot open file with body motion: " + motionFileName << std::endl;
         exit( 1 );
      }
      std::string line;
      // skip first two lines (header)
      std::getline( file, line );
      std::getline( file, line );

      while( std::getline( file, line ) ){
         // clean multiple spaces
         while( line.find( "  " ) != std::string::npos ) {
             line.replace( line.find( "  " ), 2, " " );
         }
         // trim leading whitespace
         line.erase( 0, line.find_first_not_of( " \t" ) );
         if ( line.empty() ) continue;

         std::istringstream iss( line );
         float t, dvdt, v, x;
         if ( !( iss >> t >> dvdt >> v >> x ) ) continue;

         time.push_back( t );
         acc.push_back( dvdt );
         vel.push_back( v );
         pos.push_back( x );
      }
      file.close();
   }

   // temp custom interpolation function
   float
   interpolate( float t, const std::vector< float >& times, const std::vector< float >& values )
   {
      if ( times.empty() ) return 0.0f;
      // clamp to range ( np.interp default behavior: extrapolate with boundary values )
      if ( t <= times.front() ) return values.front();
      if ( t >= times.back() ) return values.back();

      // Find the interval ( lower_bound returns first >= t )
      auto it = std::lower_bound( times.begin(), times.end(), t );
      if ( it == times.end() ) return values.back();
      if ( *it == t ) return values[ it - times.begin() ];

      size_t i = it - times.begin() - 1;
      float t0 = times[ i ];
      float t1 = times[ i + 1 ];
      float v0 = values[ i ];
      float v1 = values[ i + 1 ];

      return v0 + ( t - t0 ) * ( v1 - v0 ) / ( t1 - t0 );
   }

   template< typename BoundaryPointer, typename TimeStepping >
   void
   assignMotion( BoundaryPointer& boundary, TimeStepping& timeStepping )
   {
      using TraitsType = TNL::SPH::SPHFluidTraits< typename SPHDefs::SPHConfig >;
      using RealType = typename  TraitsType::RealType;
      using VectorType = typename TraitsType::VectorType;

      RealType t = timeStepping.getTime();
      RealType dt = timeStepping.getTimeStep();
      RealType vx = interpolate( t, time, vel );
      RealType accx = interpolate( t, time, acc );

      auto r_view = boundary->getPoints().getView();
      auto v_view = boundary->getVariables()->v.getView();
      auto dvdt_view = boundary->getVariables()->a.getView();
      const auto marker_view = boundary->getVariables()->marker.getConstView();

      const VectorType vel = { vx, 0 };
      const VectorType acc = { accx, 0 };

      auto moveSquare = [ = ] __cuda_callable__( int i ) mutable
      {
         if( marker_view[ i ] == movingBoundary )
         {
            r_view[ i ] += dt * vel;
            v_view[ i ] = vel;
            dvdt_view[ i ] = acc;
         }
      };
      boundary->getParticles()->forAll( moveSquare );
   }

   std::vector< float > time;
   std::vector< float > pos;
   std::vector< float > vel;
   std::vector< float > acc;

};

} // namespace userCodedFunctions

