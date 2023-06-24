namespace TNL {
namespace ParticleSystem {
namespace SPH {
namespace SimulationConstrolConfiguration {

/**
 * PARAMETERS FOR SIMULATION CONTROL (necessary)
 *
 * This class is used to store core parameters for simulation control and
 * simulation initialization. This includes variables such as paths for loading
 * and saving files or the length of the simulation and the frequency of saving outputs.
 *
 * It is necessary to enter:
 * - const std::string inputParticleFile - input file with initial particle configuration
 * - const std::string inputParticleFile_bound - input file with initial boundary configuration
 * - std::string outputFilaName - path to store results
 * - const float endTime - end time of the simulation
 * - const float outputTime - time interval for sagin the results
 */
class SPHSimulationControl
{
   public:

   const std::string inputParticleFile = "dambreak_fluid.vtk";
   const std::string inputParticleFile_bound = "dambreak_boundary.vtk";

   const std::string inputParticleFile_g1 = "dambreak_fluid_g1.vtk";
   const std::string inputParticleFile_bound_g1 = "dambreak_boundary_g1.vtk";

   const std::string inputParticleFile_g2 = "dambreak_fluid_g2.vtk";
   const std::string inputParticleFile_bound_g2 = "dambreak_boundary_g2.vtk";

   std::string outputFileName = "results/particles";

   const float endTime = 0.75;
   const float outputTime = 0.04f;
};

} //SimulationConstrolConfiguration
} //SPH
} //ParticleSystem
} //TNL

