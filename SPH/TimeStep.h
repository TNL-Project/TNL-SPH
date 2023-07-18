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

   const RealType
   getEndTime() const
   {
      return endTime;
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
class VariableTimeStep : public ConstantTimeStep< SPHConfig >
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   VariableTimeStep( RealType initialTimeStep, RealType endTime )
   : ConstantTimeStep< SPHConfig >( initialTimeStep, endTime ) {};

   template< typename FluidPointer, typename SPHState >
   void computeTimeStep( FluidPointer& fluid, SPHState& params )
   {
      const RealType delta_r_max = 0.1f * params.h;
      const RealType CFL = params.CFL;

      const auto view_v = fluid->getFluidVariables()->v.getConstView();
      const auto view_a = fluid->getFluidVariables()->a.getConstView();

      auto fetch = [=] __cuda_callable__ ( int i )
      {
         const RealType dt_vel = delta_r_max / l2Norm( view_v[ i ] );
         const RealType dt_acc = sqrt( delta_r_max / ( 0.5f * l2Norm( view_a[ i ] ) ) );

         return CFL * min( dt_vel, dt_acc );
      };
      RealType newTimeStep = Algorithms::reduce< DeviceType >( fluid->getFirstActiveParticle(),
                                                               fluid->getLastActiveParticle() + 1,
                                                               fetch,
                                                               TNL::Max() );

      std::cout << "[ VariableTimeStep.computeTimeStep ] params.dtMin: " << params.dtMin << std::endl;
      std::cout << "[ VariableTimeStep.computeTimeStep ] params.dtInit: " << params.dtInit << std::endl;
      std::cout << "[ VariableTimeStep.computeTimeStep ] newTimeStep: " << newTimeStep << std::endl;

      this->timeStep = std::max( std::min( newTimeStep, params.dtInit ), params.dtMin );
   }
};

} // SPH
} // ParticleSystem
} // TNL

