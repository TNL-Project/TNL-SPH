#include <iostream>
#include <fstream> //temp, to write output

#include <TNL/Devices/Cuda.h>
#include <string>
#include <sys/types.h>

/**
 * Particle system.
 */
#include "../../../Particles/Particles.h"
#include "../../../Particles/neighbourSearch.h"

/**
 * Particle system reader.
 **/
#include "../../../Readers/VTKReader.h"
#include "../../../Writers/VTKWriter.h"
#include "../../../Readers/readSPHSimulation.h"

/**
 * Case configuration
 * One configuration for particle system, one for SPH.
 */
#include "ParticlesConfig.h"
#include "SPHCaseConfig.h"
#include "MeasuretoolConfig.h"
#include "SimulationControlConfig.h"

/**
 * SPH general toolds.
 */
#include "../../../SPH/SPHOpen.h"

/**
 * SPH model.
 */
#include "../../../SPH/Models/WCSPH_DBCopen/Variables.h"
#include "../../../SPH/Models/WCSPH_DBCopen/Interactions.h"
#include "../../../SPH/Models/EquationOfState.h"

#include "../../../SPH/Models/EquationOfState.h"
#include "../../../SPH/Models/DiffusiveTerms.h"
#include "../../../SPH/Kernels.h"

using namespace TNL;

int main( int argc, char* argv[] )
{
   /**
    * Number of particles
    */
   using Device = Devices::Cuda;

   /**
    * Load simulation configs.
    * - Particle system config:
    *   config for definition of particle system (datatypes, dimension,...)
    *   config with parameters of particle system (domain size, search radius,...)
    *
    * - Configuration of particle system.
    *   config with initial parameters of the particle system
    *
    * - SPH method config:
    *   config with parameteres and constants of the SPH method
    *
    * - Simulation control:
    *   config with path to initial condition, path to store results, end time etc.
    *
    * - Measuretool:
    *   config for pressure measurement
    *   config for water level measurement
    *   config for grid interpolation
    */
   using ParticlesConfig = ParticleSystemConfig::ParticleSystemConfig< Device >;
   using SPHSimulationConfig = TNL::ParticleSystem::SPH::SPHSimpleFluidConfig< ParticlesConfig >;

   using ParticlesInitParameters = ParticleSystem::ParticleSystemConfig::ParticleInitialSetup;

   using SPHConfig = SPH::SPHCaseConfig< Device >;

   using SimulationControl = TNL::ParticleSystem::SPH::SimulationConstrolConfiguration::SPHSimulationControl;

   using MeasuretoolInitGridInterpolation = TNL::ParticleSystem::SPH::MeasuretoolConfiguration::GridInterpolationConfig< SPHConfig >;

   /**
    * Particle and neighbor search model.
    */
   using ParticleSystem = typename ParticleSystem::Particles< ParticlesConfig, Device >;
   using NeighborSearch = typename TNL::ParticleSystem::NeighborSearch< ParticlesConfig, ParticleSystem >;

   /**
    * Define simulation SPH model and SPH formulation.
    *
    * - SPHModel: is the model of used SPH method (WCSPH_DBC, WCSPH_BI, RSPH, etc.)
    *   IMPORATANT: Constants and parameters of used model have to be defined in the SPHConfig.
    *
    * - SPHSimulation: defines the type of problem (simple fluid, problem with open or
    *   moving boundaries or multiphase flows). For the chosen type of simulation,
    *   appropriate SPH scheme is required!
    */
   using SPHModel = TNL::ParticleSystem::SPH::WCSPH_DBC< ParticleSystem, SPHConfig >;
   using SPHSimulation = TNL::ParticleSystem::SPH::SPHOpenSystem< SPHModel, ParticleSystem, NeighborSearch >;

   /**
    * Define particular schemes of the used SPH model.
    * Here in the case of WCSPH, we work with DiffusiveTerm, ViscousTerm and EquationOfState.
    *
    * IMPORATANT: Constants and parameters of used schemes have to be defined in the SPHConfig.
    */
   using DiffusiveTerm = TNL::ParticleSystem::SPH::MolteniDiffusiveTerm< SPHConfig >;
   using ViscousTerm = TNL::ParticleSystem::SPH::ArtificialViscosity< SPHConfig >;
   using EOS = TNL::ParticleSystem::SPH::TaitWeaklyCompressibleEOS< SPHConfig >;

   /**
    * Particle reader and writer.
    */
   using Reader = TNL::ParticleSystem::Readers::VTKReader;
   using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;

   /**
    * Define time step control.
    * There is const time step option and variable time step option.
    */
   using TimeStepping = TNL::ParticleSystem::SPH::ConstantTimeStep< SPHConfig >;

   /**
    * Load simulation parameters.
    */
   SPHSimulationConfig mySPHSimulationConfig;
   mySPHSimulationConfig.template loadParameters< ParticlesInitParameters >();

   /**
    * Load simulation control (file names, time steps,...)
    */
   SimulationControl mySimulationControl;

   /**
    * Create the simulation.
    */
   SPHSimulation mySPHSimulation( mySPHSimulationConfig );
   std::cout << mySPHSimulation << std::endl;

   //Resolve this

   mySPHSimulation.addOpenBoundaryPatch( ParticlesConfig_inlet::numberOfParticles, ParticlesConfig_inlet::numberOfAllocatedParticles,
         ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

   mySPHSimulation.addOpenBoundaryPatch( ParticlesConfig_inlet2::numberOfParticles, ParticlesConfig_inlet2::numberOfAllocatedParticles,
         ParticlesConfig::searchRadius, ParticlesConfig::gridXsize * ParticlesConfig::gridYsize );

   //-------------------------- up to this point
   /**
    * Read the particle file.
    */
   //:TNL::ParticleSystem::ReadParticles< ParticlesConfig, Reader > myFluidReader( inputParticleFile );
   //:myFluidReader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.fluid->particles->getPoints() ) ;

   //:myFluidReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
   //:      mySPHSimulation.fluid->getFluidVariables()->rho, "Density" );
   //:myFluidReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
   //:      mySPHSimulation.fluid->getFluidVariables()->p, "Pressure" );
   //:myFluidReader.template readParticleVariable< SPHModel::VectorArrayType, float >(
   //:      mySPHSimulation.fluid->getFluidVariables()->v, "Velocity" );

   TNL::ParticleSystem::ReadParticles< ParticlesConfig_bound, Reader > myBoundaryReader( inputParticleFile_bound );
   myBoundaryReader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.boundary->particles->getPoints() ) ;

   myBoundaryReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->rho, "Density" );
   myBoundaryReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->p, "Pressure" );
   myBoundaryReader.template readParticleVariable< SPHModel::VectorArrayType, float >(
         mySPHSimulation.boundary->getBoundaryVariables()->v, "Velocity" );

   /**
    * Setup inlet.
    */
   TNL::ParticleSystem::ReadParticles< ParticlesConfig_inlet, Reader > myInletReader( inputParticleFile_inlet );
   myInletReader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.openBoundaryPatches[ 0 ]->particles->getPoints() ) ;

   myInletReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.openBoundaryPatches[ 0 ]->getOpenBoundaryVariables()->rho, "Density" );
   myInletReader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
         mySPHSimulation.openBoundaryPatches[ 0 ]->getOpenBoundaryVariables()->p, "Pressure" );
   myInletReader.template readParticleVariable2D< SPHModel::VectorArrayType, float >(
         mySPHSimulation.openBoundaryPatches[ 0 ]->getOpenBoundaryVariables()->v, "Velocity" );

   /**
    * Setup parameters for open boundary.
    */
   mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.orientation = {
      SPHConfig::INLET::orientation_x, SPHConfig::INLET::orientation_y };
   mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.velocity = {
      SPHConfig::INLET::velocity_x, SPHConfig::INLET::velocity_y };
   mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.density = SPHConfig::INLET::inlet_density;
   mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.bufferWidth = {
      SPHConfig::INLET::bufferWidth_x, SPHConfig::INLET::bufferWidth_y };
   mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.position = {
      SPHConfig::INLET::position_x, SPHConfig::INLET::position_y };

   std::cout << "Inlet parameters: " << std::endl;
   std::cout << "Orientation ................. " << mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.orientation << std::endl;
   std::cout << "Velocity .................... " << mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.velocity << std::endl;
   std::cout << "BufferWidth ................. " << mySPHSimulation.openBoundaryPatches[ 0 ]->parameters.bufferWidth << std::endl;

  /**
   * Setup second inlet.
   */
  TNL::ParticleSystem::ReadParticles< ParticlesConfig_inlet2, Reader > myInlet2Reader( inputParticleFile_inlet2 );
  myInlet2Reader.template readParticles< ParticleSystem::PointArrayType >( mySPHSimulation.openBoundaryPatches[ 1 ]->particles->getPoints() ) ;

  myInlet2Reader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
        mySPHSimulation.openBoundaryPatches[ 1 ]->getOpenBoundaryVariables()->rho, "Density" );
  myInlet2Reader.template readParticleVariable< SPHModel::ScalarArrayType, float >(
        mySPHSimulation.openBoundaryPatches[ 1 ]->getOpenBoundaryVariables()->p, "Pressure" );
  myInlet2Reader.template readParticleVariable2D< SPHModel::VectorArrayType, float >(
        mySPHSimulation.openBoundaryPatches[ 1 ]->getOpenBoundaryVariables()->v, "Velocity" );

  /**
   * Setup parameters for open boundary.
   */
  mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.orientation = {
     SPHConfig::INLET2::orientation_x, SPHConfig::INLET2::orientation_y };
  mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.velocity = {
     SPHConfig::INLET2::velocity_x, SPHConfig::INLET2::velocity_y };
  mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.density = SPHConfig::INLET2::inlet_density;
  mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.bufferWidth = {
     SPHConfig::INLET2::bufferWidth_x, SPHConfig::INLET2::bufferWidth_y };
  mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.position = {
     SPHConfig::INLET2::position_x, SPHConfig::INLET2::position_y };

  std::cout << "Inlet2 parameters: " << std::endl;
  std::cout << "Orientation ................. " << mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.orientation << std::endl;
  std::cout << "Velocity .................... " << mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.velocity << std::endl;
  std::cout << "BufferWidth ................. " << mySPHSimulation.openBoundaryPatches[ 1 ]->parameters.bufferWidth << std::endl;

   mySPHSimulation.fluid->particles->setGridOrigin( { ParticlesConfig::gridXbegin, ParticlesConfig::gridYbegin } );
   mySPHSimulation.boundary->particles->setGridOrigin( { ParticlesConfig::gridXbegin, ParticlesConfig::gridYbegin } );
   mySPHSimulation.openBoundaryPatches[ 0 ]->particles->setGridOrigin( { ParticlesConfig::gridXbegin, ParticlesConfig::gridYbegin } );
   mySPHSimulation.openBoundaryPatches[ 1 ]->particles->setGridOrigin( { ParticlesConfig::gridXbegin, ParticlesConfig::gridYbegin } );

   mySPHSimulation.fluid->particles->setGridSize( { ParticlesConfig::gridXsize, ParticlesConfig::gridYsize } );
   mySPHSimulation.boundary->particles->setGridSize( { ParticlesConfig::gridXsize, ParticlesConfig::gridYsize } );
   mySPHSimulation.openBoundaryPatches[ 0 ]->particles->setGridSize( { ParticlesConfig::gridXsize, ParticlesConfig::gridYsize } );
   mySPHSimulation.openBoundaryPatches[ 1 ]->particles->setGridSize( { ParticlesConfig::gridXsize, ParticlesConfig::gridYsize } );

   /**
    * Define timers to measure computation time.
    */
   TNL::Timer timer_search, timer_interact, timer_integrate, timer_inlet, timer_pressure;
   TNL::Timer timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells;

   /**
    * TEMP: Determine number of interation for constant timestep.
    * Perform simulation main loop.
    */
   int steps = endTime / SPHConfig::dtInit;
   std::cout << "Number of steps: " << steps << std::endl;


   for( unsigned int iteration = 0; iteration < steps; iteration ++ )
   {
      std::cout << "STEP: " << iteration << std::endl;

      /**
       * Find neighbors within the SPH simulation.
       */
      timer_search.start();
      mySPHSimulation.PerformNeighborSearch(
            iteration, timer_search_reset, timer_search_cellIndices, timer_search_sort, timer_search_toCells );
      timer_search.stop();
      std::cout << "Search... done. " << std::endl;

      /**
       * Perform interaction with given model.
       */
      timer_interact.start();
      mySPHSimulation.template Interact< SPH::WendlandKernel2D, DiffusiveTerm, ViscousTerm, EOS >();
      timer_interact.stop();
      std::cout << "Interact... done. " << std::endl;

      //#include "outputForDebug.h"

      /**
       * Perform time integration, i.e. update particle positions.
       */
      timer_integrate.start();
      if( iteration % 20 == 0 ) {
         mySPHSimulation.integrator->IntegrateEuler< typename SPHSimulation::FluidPointer >( SPHConfig::dtInit, mySPHSimulation.fluid );
         mySPHSimulation.integrator->IntegrateEulerBoundary< typename SPHSimulation::BoundaryPointer >( SPHConfig::dtInit, mySPHSimulation.boundary );
         timer_inlet.start();
         mySPHSimulation.integrator->updateBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 0 ] );
         //mySPHSimulation.integrator->updateBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 1 ] );
         mySPHSimulation.integrator->updateOutletBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 1 ] );
         timer_inlet.stop();
      }
      else {
         mySPHSimulation.integrator->IntegrateVerlet< typename SPHSimulation::FluidPointer >( SPHConfig::dtInit, mySPHSimulation.fluid );
         mySPHSimulation.integrator->IntegrateVerletBoundary< typename SPHSimulation::BoundaryPointer >( SPHConfig::dtInit, mySPHSimulation.boundary );
         timer_inlet.start();
         mySPHSimulation.integrator->updateBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 0 ] );
         //mySPHSimulation.integrator->updateBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 1 ] );
         mySPHSimulation.integrator->updateOutletBuffer< typename SPHSimulation::FluidPointer, typename SPHSimulation::OpenBoudaryPatchPointer >( SPHConfig::dtInit, mySPHSimulation.fluid, mySPHSimulation.openBoundaryPatches[ 1 ] );
         timer_inlet.stop();
      }
      timer_integrate.stop();
      std::cout << "Integration... done. " << std::endl;

      /**
       * Output particle data
       */
      if( ( iteration % outputStep ==  0) && (iteration > 0) )
      {
         /**
          * Compute pressure from density.
          * This is not necessary since we do this localy, if pressure is needed.
          * Its useful for output anywal.
          */
         timer_pressure.start();
         mySPHSimulation.model->template ComputePressureFromDensity< EOS >( mySPHSimulation.fluid->variables, mySPHSimulation.fluid->particles->getNumberOfParticles() ); //TODO: FIX.
         mySPHSimulation.model->template ComputePressureFromDensity< EOS >( mySPHSimulation.boundary->variables, mySPHSimulation.boundary->particles->getNumberOfParticles() ); //TODO: FIX.
         timer_pressure.stop();
         std::cout << "Compute pressure... done. " << std::endl;

         std::string outputFileNameFluid = outputFileName + std::to_string( iteration ) + "_fluid.vtk";
         std::ofstream outputFileFluid ( outputFileNameFluid, std::ofstream::out );
         Writer myWriter( outputFileFluid, VTK::FileFormat::ascii );
         myWriter.writeParticles( *mySPHSimulation.fluid->particles );
         myWriter.template writePointData< SPHModel::ScalarArrayType >(
               mySPHSimulation.fluid->getFluidVariables()->rho, "Density", mySPHSimulation.fluid->particles->getNumberOfParticles(), 1 );
         myWriter.template writeVector< SPHModel::VectorArrayType, SPHConfig::RealType >(
               mySPHSimulation.fluid->getFluidVariables()->v, "Velocity", 3, mySPHSimulation.fluid->particles->getNumberOfParticles() );

         std::string outputFileNameBoundary = outputFileName + std::to_string( iteration ) + "_boundary.vtk";
         std::ofstream outputFileBoundary ( outputFileNameBoundary, std::ofstream::out );
         Writer myWriterBoundary( outputFileBoundary, VTK::FileFormat::ascii );
         myWriterBoundary.writeParticles( *mySPHSimulation.boundary->particles );
         myWriterBoundary.template writePointData< SPHModel::ScalarArrayType >(
               mySPHSimulation.boundary->getBoundaryVariables()->p, "Pressure", mySPHSimulation.boundary->particles->getNumberOfParticles(), 1 );
         //myWriterBoundary.template writeVector< SPHModel::VectorArrayType, SPHConfig::RealType >(
         //      mySPHSimulation.boundary->getBoundaryVariables()->v, "Velocity", 3, mySPHSimulation.boundary->particles->getNumberOfParticles() );
      }
   }

   /**
    * Output simulation stats.
    */
   float totalTime = ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() );

   float totalTimePerStep = totalTime / steps;

   std::cout << std::endl << "COMPUTATION TIME:" << std::endl;
   std::cout << "Search........................................ " << timer_search.getRealTime() << " sec." << std::endl;
   std::cout << "Search (average time per step)................ " << timer_search.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Search (percentage)........................... " << timer_search.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Reset ..................................... " << timer_search_reset.getRealTime() << " sec." << std::endl;
   std::cout << " - Reset (average time per step).............. " << timer_search_reset.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Reset (percentage)......................... " << timer_search_reset.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Index by cell ............................. " << timer_search_cellIndices.getRealTime() << " sec." << std::endl;
   std::cout << " - Index by cell (average time per step)...... " << timer_search_cellIndices.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Index by cell (percentage)................. " << timer_search_cellIndices.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Sort ...................................... " << timer_search_sort.getRealTime() << " sec." << std::endl;
   std::cout << " - Sort (average time per step)............... " << timer_search_sort.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Sort (percentage).......................... " << timer_search_sort.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Particle to cell .......................... " << timer_search_toCells.getRealTime() << " sec." << std::endl;
   std::cout << " - Particle to cell (average time per step)... " << timer_search_toCells.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Particle to cell (percentage).............. " << timer_search_toCells.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Interaction................................... " << timer_interact.getRealTime() << " sec." << std::endl;
   std::cout << "Interaction (average time per step)........... " << timer_interact.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Interaction (percentage)...................... " << timer_interact.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Integrate..................................... " << timer_integrate.getRealTime() << " sec." << std::endl;
   std::cout << "Integrate (average time per step)............. " << timer_integrate.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Integrate (percentage)........................ " << timer_integrate.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << " - Update inlet............................... " << timer_inlet.getRealTime() << " sec." << std::endl;
   std::cout << " - Update inlet (average time per step)....... " << timer_inlet.getRealTime() / steps << " sec." << std::endl;
   std::cout << " - Update inlet (percentage).................. " << timer_inlet.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Pressure update............................... " << timer_pressure.getRealTime() << " sec." << std::endl;
   std::cout << "Pressure update (average time per step)....... " << timer_pressure.getRealTime() / steps << " sec." << std::endl;
   std::cout << "Pressure update (percentage).................. " << timer_pressure.getRealTime() / totalTime * 100 << " %." << std::endl;
   std::cout << "Total......................................... " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) << " sec." << std::endl;
   std::cout << "Total (average time per step)................. " << ( timer_search.getRealTime() + \
   + timer_interact.getRealTime() + timer_integrate.getRealTime() + timer_pressure.getRealTime() ) / steps << " sec." << std::endl;

   std::cout << "\nDone ... " << std::endl;
}

