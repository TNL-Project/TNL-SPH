#pragma once

#include "TNL/Timer.h"
#include <map>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

class TimerMeasurement
{
public:

   TimerMeasurement()
   {
      timers.insert( { "search", TNL::Timer()  } );
      timers.insert( { "search_reset", TNL::Timer()  } );
      timers.insert( { "search_cellIndices", TNL::Timer()  } );
      timers.insert( { "search_sort", TNL::Timer()  } );
      timers.insert( { "search_toCells", TNL::Timer()  } );
      timers.insert( { "interact", TNL::Timer()  } );
      timers.insert( { "integrate", TNL::Timer()  } );
      timers.insert( { "total", TNL::Timer()  } );
   }

   void
   addTimer( const std::string keyword )
   {
      timers.insert( { keyword, TNL::Timer() } );
   }

   void
   start( const std::string keyword )
   {
      this->timers[ keyword ].start();
   }

   void
   stop( const std::string keyword )
   {
      this->timers[ keyword ].stop();
   }

   void
   print( const int stepsTotal )
   {
      std::ostream &str = std::cout;
      TNL::Logger logger( 100, str );
      logger.writeHeader( "Computation time" );
      //float totalTime = this->timers[ "total" ].getRealTime();

      for ( auto const& [ key, val ] : this->timers )
      {
         logger.writeParameter( key, val.getRealTime() );
         logger.writeParameter( key + "-average", val.getRealTime() / stepsTotal );
         logger.writeParameter( key + "-percentage", val.getRealTime() );
         logger.writeSeparator();
      }
   }

protected:

   std::map< std::string, TNL::Timer > timers;

};


} // SPH
} // ParticleSystem
} // TNL

