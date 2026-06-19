import sys
sys.path.append( '..' )
import damBreak2D_plot_comparison as dbp
from pathlib import Path

conf_dir = Path(__file__).parent
example_dir = (conf_dir / ".." /".." / ".." / "examples/WCSPH-BI/damBreak2D_WCSPH-BI").resolve()
resources_dir = (conf_dir / ".." / ".." / ".." /"examples/resources/damBreak2D/damBreak2D_experimentalDataLobovsky2014").resolve()
output_dir = conf_dir / "damBreak2D_WCSPH-BI_MR"

# LR local
dbp.plot_multiple_pressure_sensors(
    [
     "results_dp0.004",
     "results_dp0.002",
     "/home/tomas/Work/sph-projects/fresh/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution/results"
    ],
    [
      r"$\Delta x = 0.004$",
      r"$\Delta x = 0.002$",
      "multiresolution"
    ],
    'damBreak2D_WCSPH-BI_MR_dp0.004-0.002',
    example_dir,
    resources_dir,
    output_dir )

## HR hop
#dbp.plot_multiple_pressure_sensors(
#    [
#     "results_dp0.002",
#     "/home/tomas/hub/hop/playground/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI/results",
#     "/home/tomas/hub/hop/playground/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution/results"
#    ],
#    [
#      r"$\Delta x = 0.002$",
#      r"$\Delta x = 0.001$",
#      "multiresolution"
#    ],
#    'damBreak2D_WCSPH-BI_MR_dp0.002-0.001',
#    example_dir,
#    resources_dir,
#    output_dir )
