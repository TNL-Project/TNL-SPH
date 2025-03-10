#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Include globaly used postpro scripts
import sys
sys.path.append('../../../src/tools')
import plotTimeStep
import plotEnergy

#TODO: Rename the default sensors labels providet by Lobovsky et. al.
case_tag = "poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC"
example_dir = Path(__file__).parent
results_dir = ( example_dir / "results" ).resolve()

plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 24
})

if __name__ == "__main__":
    import argparse
    import os
    argparser = argparse.ArgumentParser(description="Dam break example postprocessing")
    argparser.add_argument("--with-paraview", default=False,
            help="perform postprocessing using paraview tools")
    argparser.add_argument("--config", default="sources/config.ini",
            help="path to the config file (relative to the path of this script)")

    # create folder for postprocessing results
    postproPath = r'./results/postprocessing'
    if not os.path.exists( postproPath ):
        os.makedirs( postproPath )

    #plotEnergy.plot_energy( results_dir, Epot0 = 264.87 )
    plotEnergy.plot_not_normalized_energy( results_dir )
    plotEnergy.plot_energy_snapshots( results_dir )
    plotEnergy.plot_not_normalized_open_boundary_energy_snapshots( results_dir )
