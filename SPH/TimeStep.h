#include <map>

namespace TNL {
namespace ParticleSystem {
namespace SPH {

template< typename SPHConfig >
class ConstantTimeStep
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   ConstantTimeStep( RealType initialTimeStep, RealType endTime )
   : timeStep( initialTimeStep ), endTime( endTime )
   {
      step = 0;
      time = 0;
   }

   bool
   runTheSimulation()
   {
      bool runSimulation = false;
      if( time < endTime )
         runSimulation = true;
      return runSimulation;
   }

   const RealType
   getTime() const
   {
      return time;
   }

   const GlobalIndexType
   getStep() const
   {
      return step;
   }

   const RealType
   getTimeStep() const
   {
      return timeStep;
   }

   void
   setTimeStep( RealType timeStep )
   {
      timeStep = timeStep;
   }

   void
   updateTimeStep()
   {
      time += timeStep;
      step += 1;
   }

   void
   addOutputTimer( const std::string keyword, const RealType outputTime )
   {
      outputTimers.insert({ keyword , { outputTime, this->getTime() } } );
   }

   void
   changeOutputTimer( const std::string keyword, const RealType outputTime )
   {
      this->outputTimers[ keyword ].first = outputTime;
      this->outputTimers[ keyword ].second = this->getTime();
   }

   bool
   checkOutputTimer( const std::string keyword )
   {
      if( this->getTime() > this->outputTimers[ keyword ].second )
      {
         this->outputTimers[ keyword ].second += this->outputTimers[ keyword ].first;
         return true;
      }

      return false;
   }

protected:

   GlobalIndexType step;
   RealType time;
   RealType timeStep;
   RealType endTime;

   //control timers for arbitrary outpus
   std::map< std::string, std::pair< RealType, RealType > > outputTimers;

};

template< typename SPHConfig >
class VariableTimeStep
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;
};


} // SPH
} // ParticleSystem
} // TNL

