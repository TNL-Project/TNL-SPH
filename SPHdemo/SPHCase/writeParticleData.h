/* WRITE PARTICLES TEMPLATE */

//using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystemToReadData >;
using Writer = TNL::ParticleSystem::Writers::VTKWriter< ParticleSystem >;
const std::string outputFileName = "output.vtk";
std::ofstream outputFile (outputFileName, std::ofstream::out);
Writer myWriter( outputFile, VTK::FileFormat::binary );
myWriter.writeParticles( *mySPHSimulation.particles );
//myWriter.template writeMetadata( 1, 0 );
myWriter.template writePointData< SPHModel::ScalarArrayType >( mySPHSimulation.model->FluidVariables.rho, "Density" );
myWriter.template writeVector< SPHModel::VectorArrayType, float >( mySPHSimulation.model->FluidVariables.v, "Velocity", 3 );

