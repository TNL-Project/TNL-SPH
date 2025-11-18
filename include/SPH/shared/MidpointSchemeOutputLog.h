#include <string>
#include <fstream>

namespace TNL {
namespace SPH {
namespace Info {

template< typename MidpointTimeIntegrationSchemePointer, typename ModelParams, typename TimeStepping >
void
midpointSchemeOutputLog( const MidpointTimeIntegrationSchemePointer& midpoint,
                         const ModelParams& modelParams,
                         const TimeStepping& timeStepping,
                         const std::string& outputPath )
{
      if( midpoint->midpointIteration == modelParams.midpointMaxInterations ) {
         std::ofstream outfile;
         outfile.open( outputPath, std::ios_base::app );
         outfile << timeStepping.getStep() << " "
                 << timeStepping.getTime() << " "
                 << midpoint->residual << " "
                 << midpoint->midpointIterationBackup << " "
                 << midpoint->midpointRelaxCoef << std::endl;
      }
}

template< typename MidpointTimeIntegrationSchemePointer, typename ModelParams, typename TimeStepping >
void
mipodintSchemeOutputSubiterationsLog( const MidpointTimeIntegrationSchemePointer& midpoint,
                                      const ModelParams& modelParams,
                                      const TimeStepping& timeStepping,
                                      const std::string& outputPath )
{
         std::ofstream outfile;
         outfile.open( outputPath, std::ios_base::app );
         outfile << timeStepping.getStep() << " "
                 << timeStepping.getTime() << " "
                 << midpoint->residual << " "
                 << midpoint->midpointIteration << " "
                 << midpoint->midpointRelaxCoef << std::endl;
}

} // Info
} // SPH
} // TNL

