#pragma once

#include <iterator>
#include <map>
#include <fstream>

namespace TNL {
namespace SPH {

template< typename SPHConfig >
class ConstantTimeStep
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   ConstantTimeStep()
   {
      step = 0;
      time = 0;
   }

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
   setEndTime( const RealType& endTime )
   {
      this->endTime = endTime;
   }

   void
   setTimeStep( const RealType& timeStep )
   {
      this->timeStep = timeStep;
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
      if( this->getTime() >= this->outputTimers[ keyword ].second )
      {
         this->outputTimers[ keyword ].second += this->outputTimers[ keyword ].first;
         return true;
      }

      return false;
   }

   void
   outputTimeStep( const std::string& outputPath ) {}

   template< typename FluidPointer, typename SPHState >
   void computeTimeStep( FluidPointer& fluid, SPHState& params ) {}

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

   VariableTimeStep() : ConstantTimeStep< SPHConfig >() {};

   VariableTimeStep( RealType initialTimeStep, RealType endTime )
   : ConstantTimeStep< SPHConfig >( initialTimeStep, endTime ) {};

   void
   outputTimeStep( const std::string& outputPath )
   {
      std::ofstream outfile;
      outfile.open(outputPath, std::ios_base::app );
      outfile << this->step << " " << this->time << " " << this->timeStep << std::endl;
   }

   template< typename FluidPointer, typename SPHState >
   void computeTimeStep( FluidPointer& fluid, SPHState& params )
   {
      const RealType delta_r_max = 0.1f * params.h;
      const RealType cfl = params.cfl;

      const auto view_v = fluid->getFluidVariables()->v.getConstView();
      const auto view_a = fluid->getFluidVariables()->a.getConstView();

      auto fetch = [=] __cuda_callable__ ( int i )
      {
         const RealType dt_vel = delta_r_max / l2Norm( view_v[ i ] );
         const RealType dt_acc = sqrt( delta_r_max / ( 0.5f * l2Norm( view_a[ i ] ) ) );

         return cfl * min( dt_vel, dt_acc );
      };
      RealType newTimeStep = Algorithms::reduce< DeviceType >( 0,
                                                               fluid->getNumberOfParticles(),
                                                               fetch,
                                                               TNL::Min() );

      //std::cout << "[ VariableTimeStep.computeTimeStep ] params.dtMin: " << params.dtMin << std::endl;
      //std::cout << "[ VariableTimeStep.computeTimeStep ] params.dtInit: " << params.dtInit << std::endl;
      //std::cout << "[ VariableTimeStep.computeTimeStep ] newTimeStep: " << newTimeStep << std::endl;

      this->timeStep = std::max( std::min( newTimeStep, params.dtInit ), params.dtMin );
   }
};

template< typename SPHConfig >
class VariableTimeStepWithReduction : public ConstantTimeStep< SPHConfig >
{
public:

   using SPHTraitsType = SPHFluidTraits< SPHConfig >;
   using DeviceType = typename SPHConfig::DeviceType;
   using GlobalIndexType = typename SPHTraitsType::GlobalIndexType;
   using RealType = typename SPHTraitsType::RealType;

   VariableTimeStepWithReduction( RealType initialTimeStep, RealType endTime )
   : ConstantTimeStep< SPHConfig >( initialTimeStep, endTime ) {};

   const RealType
   getMaxViscosEffect() const
   {
      return maxViscosityEffect;
   }

   void setMaxViscosityEffect( RealType maxViscosity )
   {
      maxViscosityEffect = maxViscosity;
   }

   template< typename FluidPointer, typename SPHState >
   void computeTimeStep( FluidPointer& fluid, SPHState& params )
   {
      const RealType delta_r_max = 0.1f * params.h;
      const RealType CFL = params.CFL;
      const RealType h = params.h;

      const auto view_v = fluid->getFluidVariables()->v.getConstView();
      const auto view_a = fluid->getFluidVariables()->a.getConstView();

      auto fetch = [=] __cuda_callable__ ( int i )
      {
         return sqrt( h / ( l2Norm( view_a[ i ] ) ) );
      };
      RealType dt_acc = Algorithms::reduce< DeviceType >( fluid->getFirstActiveParticle(),
                                                          fluid->getLastActiveParticle() + 1,
                                                          fetch,
                                                          TNL::Max() );

      RealType dt_visco = params.h / ( params.speedOfSound + params.h * maxViscosityEffect );

      this->timeStep = std::max( std::min( dt_acc, dt_visco ), params.dtMin );
   }

protected:

   RealType maxViscosityEffect;

};


} // SPH
} // TNL

