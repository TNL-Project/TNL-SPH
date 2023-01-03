//Create temp particle system to read data in:
using ParticleSystemToReadData = typename ParticleSystem::Particles< ParticlesConfig, Devices::Host >;
ParticleSystemToReadData particlesToRead( ParticlesConfig::numberOfParticles, ParticlesConfig::searchRadius );

using ParticleSystemToReadData_bound = typename ParticleSystem::Particles< ParticlesConfig_bound, Devices::Host >;
ParticleSystemToReadData_bound particlesToRead_bound( ParticlesConfig_bound::numberOfParticles, ParticlesConfig::searchRadius );

const std::string inputFileName = inputParticleFile;
Reader myReader( inputFileName );
myReader.detectParticleSystem();
//myReader.template loadParticle< ParticleSystem >( *mySPHSimulation.particles );
myReader.template loadParticle< ParticleSystemToReadData >( particlesToRead );

const std::string inputFileName_bound = inputParticleFile_bound;
Reader myReader_bound( inputFileName_bound );
myReader_bound.detectParticleSystem();
//myReader.template loadParticle< ParticleSystem >( *mySPHSimulation.particles );
myReader_bound.template loadParticle< ParticleSystemToReadData_bound >( particlesToRead_bound );

/**
 * Setup type for boundary particles and initial condition.
 */
//#include "asignIinitialCondition.h"
#include "asignInitialCondition_fluid.h"
std::cout << " Fluid vars loaded. " << std::endl;
#include "asignInitialCondition_boundary.h"
std::cout << " Boundary vars loaded. " << std::endl;

