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

   ConstantTimeStep( RealType initialTimeStep )
   : timeStep( initialTimeStep )
   {
      step = 0;
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

   protected:
   GlobalIndexType step;
   RealType time;
   RealType timeStep;

};

template< typename SPHConfig >
class VariableTimeStep
{
   public:
   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   ConstantTimeStep( RealType initialTimeStep )
   : timeStep( initialTimeStep )
   {
      step = 0;
   }

   const RealType
   getTimeStep() const
   {
      return dt;
   }

   //TODO: Compute time step form flow variables.
   void
   computeTimeStep()
   {}

   protected:
   GlobalIndexType step;
   RealType time;
   RealType timeStep;

};


} // SPH
} // ParticleSystem
} // TNL

