#pragma once

#include "TNL/Timer.h"
#include <map>

namespace TNL {
namespace SPH {

class TimerMeasurement
{
public:

   TimerMeasurement()
   {
      timers.insert( { "search", std::make_pair( TNL::Timer(), true ) } );
      timers.insert( { "search_reset", std::make_pair( TNL::Timer(), false ) } );
      timers.insert( { "search_cellIndices", std::make_pair( TNL::Timer(), false ) } );
      timers.insert( { "search_sort", std::make_pair( TNL::Timer(), false ) } );
      timers.insert( { "search_toCells", std::make_pair( TNL::Timer(), false ) } );
      timers.insert( { "interact", std::make_pair( TNL::Timer(), true ) } );
      timers.insert( { "integrate", std::make_pair( TNL::Timer(), true ) } );
   }

   void
   addTimer( const std::string keyword, const bool subTimer = true )
   {
      timers.insert( { keyword, std::make_pair( TNL::Timer(), subTimer ) } );
   }

   void
   start( const std::string keyword )
   {
      this->timers[ keyword ].first.start();
   }

   void
   stop( const std::string keyword )
   {
      this->timers[ keyword ].first.stop();
   }

   void
   writeInfo( TNL::Logger& logger, const int stepsTotal ) const noexcept
   {
      logger.writeHeader( "Computation time" );
      float totalTime = 0.f;
      for ( auto const& [ key, val ] : this->timers ) {
         if ( val.second == true )
            totalTime += val.first.getRealTime();
      }

      for ( auto const& [ key, val ] : this->timers ) {
         logger.writeParameter( key, val.first.getRealTime() );
         logger.writeParameter( key + "-average", val.first.getRealTime() / stepsTotal );
         logger.writeParameter( key + "-percentage", val.first.getRealTime() / totalTime * 100 );
         logger.writeSeparator();
      }
      logger.writeParameter( "Total time:", totalTime );
      logger.writeSeparator();
   }

protected:

   std::map< std::string, std::pair< TNL::Timer, bool >> timers;

};


} // SPH
} // TNL

