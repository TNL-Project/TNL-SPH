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

   const std::string inputParticleFile = "openchannel_fluid.vtk";
   const std::string inputParticleFile_bound = "openchanlnel_boundary.vtk";
   const std::string inputParticleFile_inlet = "openchanlnel_inlet.vtk";
   const std::string inputParticleFile_outlet = "openchanlnel_outlet.vtk";
   std::string outputFileName = "results/particles";

   const float endTime = 0.5;
   const float outputTime = 0.04f;
};

} //SimulationConstrolConfiguration
} //SPH
} //ParticleSystem
} //TNL

