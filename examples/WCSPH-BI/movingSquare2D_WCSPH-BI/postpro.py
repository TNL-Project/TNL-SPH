#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Include globaly used postpro scripts
import sys
sys.path.append('../../../src/tools')

#TODO: Rename the default sensors labels providet by Lobovsky et. al.
case_tag = "movingSquare2D_WCSPH-BI"
example_dir = Path(__file__).parent
#resources_dir = (example_dir / ".." / ".." / "resources" / "damBreak2D" / "damBreak2D_experimentalDataLobovsky2014" ).resolve()
resources_dir = Path('/home/tomas/Downloads/tmp/SPHERIC_Benchmark')
results_dir = ( example_dir / "results" ).resolve()

# setup plot parameters
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 24
})

def plot_forces( results_dir ):


        forces = np.genfromtxt( results_dir / "force.dat", delimiter=" " )
        referential_data = np.loadtxt(resources_dir / "Force_Re100.dat", skiprows=2)
        """
        t    = referential_data[:, 0]
        Cd_p = referential_data[:, 1]
        Cd_v = referential_data[:, 2]
        """

        # superpone square acceleration manually
        motion_data = np.loadtxt(resources_dir / "Motion_Body.dat", skiprows=2)
        """
        t     = motion_data[:, 0]
        acc   = motion_data[:, 1]
        vel   = motion_data[:, 2]
        displ = motion_data[:, 3]
        """

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( forces[ :, 1 ], forces[ :, 2 ], label=r'$F_{pressure}$', linewidth=2, color='b'  )
        ax.plot( forces[ :, 1 ], forces[ :, 3 ], label=r'$F_{visco}$', linewidth=2, color='g'  )
        ax.plot( forces[ :, 1 ], forces[ :, 2 ] + forces[ :, 3 ], label=r'$F$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ F $ [N]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(results_dir) + f"/postprocessing/forces.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( forces[ :, 1 ], 2 * forces[ :, 2 ], label=r'$c_{D,\; pressure}$', linewidth=2, color='b'  )
        ax.plot( forces[ :, 1 ], 2 * forces[ :, 3 ], label=r'$c_{D,\; viscous}$', linewidth=2, color='g'  )
        ax.plot( forces[ :, 1 ], 2 * (forces[ :, 2 ] + forces[ :, 3 ]) / ( 1 * 1**2 * 1 * 1 ), label=r'$c_D$', linewidth=2, color='k'  )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 1 ] + referential_data[ :, 2 ], label=r'$c_D \; ref.$', linewidth=1, color='k', linestyle='--')
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ c_D $ [-]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(results_dir) + f"/postprocessing/cd.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( forces[ :, 1 ], 2 * forces[ :, 2 ], label=r'$c_{D,\; pressure}$', linewidth=2, color='b'  )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 1 ], label=r'$c_{D,\; pressure}\; ref.$ ', linewidth=1, color='b', linestyle='--' )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ c_{D,\; pressure} $ [-]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(results_dir) + f"/postprocessing/cd_pressure.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( forces[ :, 1 ], 2 * forces[ :, 3 ], label=r'$c_{D,\; viscous}$', linewidth=2, color='g'  )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 2 ], label=r'$c_{D,\; viscous}\; ref.$.', linewidth=1, color='g', linestyle='--'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ c_{D,\; viscous} $ [-]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(results_dir) + f"/postprocessing/cd_visous.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( motion_data[ :, 0 ], motion_data[ :, 1 ], label=r'$dv/dt$', linewidth=2, color='r'  )
        ax.plot( motion_data[ :, 0 ], motion_data[ :, 2 ], label=r'$v$', linewidth=2, color='b'  )
        ax.plot( motion_data[ :, 0 ], motion_data[ :, 3 ], label=r'$x$', linewidth=2, color='g'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$x$ [m], $v$ [m/s], $a$ [m$^2$/s]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(results_dir) + f"/postprocessing/square_motion.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":
    import argparse
    import os
    argparser = argparse.ArgumentParser(description="Dam break example postprocessing")
    argparser.add_argument("--with-paraview", default=False,
            help="perform postprocessing using paraview tools")
    argparser.add_argument("--config", default="sources/config.ini",
            help="path to the config file (relative to the path of this script)")

    # create folder for postprocessing results
    postproPath = str(results_dir) + '/postprocessing'
    if not os.path.exists( postproPath ):
        os.makedirs( postproPath )

    # Global postprocessing tools: group results
    import writeParaviewSeriesFile
    #writeParaviewSeriesFile.generate_series( results_dir )

    # Global postprocessing tools: plot time step log and energy
    import plotEnergy
    #plotTimeStep.plot_time_step( results_dir )
    plotEnergy.plot_energy( results_dir, Epot0 = 264.87 )
    plotEnergy.plot_energy_snapshots( results_dir, Epot0 = 264.87 )

    # ---
    plot_forces( results_dir )
