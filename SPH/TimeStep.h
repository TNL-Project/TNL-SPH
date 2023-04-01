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

   protected:
   GlobalIndexType step;
   RealType time;
   RealType timeStep;
   RealType endTime;

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

