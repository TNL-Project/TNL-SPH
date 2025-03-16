#pragma once

#include <TNL/Benchmarks/Benchmarks.h>
#include "TNL/Timer.h"
#include <map>
#include <string>

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

   float
   getTotalTime() const
   {
      float totalTime = 0.f;
      for ( auto const& [ key, val ] : this->timers ) {
         if ( val.second == true )
            totalTime += val.first.getRealTime();
      }
      return totalTime;
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

   void
   writeInfoToJson( const std::string resultsPath, const int stepsTotal )
   {
      std::map< std::string, std::string > timeResults;

      float totalTime = 0.f;
      for ( auto const& [ key, val ] : this->timers ) {
         if ( val.second == true )
            totalTime += val.first.getRealTime();
      }

      for ( auto const& [ key, val ] : this->timers ) {
         timeResults.insert( { key, std::to_string( val.first.getRealTime() ) } );
         timeResults.insert( { key + "-average", std::to_string( val.first.getRealTime() / stepsTotal ) } );
         timeResults.insert( { key + "-percentage", std::to_string( val.first.getRealTime() / totalTime * 100 ) } );
      }
      timeResults.insert( { "total", std::to_string( totalTime ) } );
      timeResults.insert( { "total-average", std::to_string( totalTime / stepsTotal ) } );
      Benchmarks::writeMapAsJson( timeResults, resultsPath, ".json" );
   }

protected:

   std::map< std::string, std::pair< TNL::Timer, bool >> timers;

};


} // SPH
} // TNL

