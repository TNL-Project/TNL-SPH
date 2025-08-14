#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Include globaly used postpro scripts
import sys
sys.path.append('../../../src/tools')
import groupResults

case_tag = "damBreak3D_WCSPH-BI"
example_dir = Path(__file__).parent
resources_dir = (example_dir / ".." / ".." / "resources" / "damBreak3D" / "damBreak3D_experimentalDataSphericIssaVioleau2006" ).resolve()
results_dir = ( example_dir / "results" ).resolve()

# setup plot parameters
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 24
})

def plot_pressure_sensors():

    experimental_data_file = resources_dir / "pressure.dat"
    experimental_data = np.genfromtxt( experimental_data_file, delimiter='\t', dtype=np.float64)

    simulation_data_file = results_dir / "sensorsPressure.dat"
    simulation_data = np.genfromtxt( simulation_data_file, delimiter=' ' )

    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = coef * sensor_snapshot_time,
    nondim_time_coef = 1. * 0.002
    # nondim_pressure_coef = coef
    nondim_pressure_coef = 1.

    for i in range( 0, 8 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        ax.plot( experimental_data[ :, 0 ],
                 experimental_data[ :, i + 1 ],
                label='Issa \& Violeau 2006', linewidth=2, color='b'  )
        ax.plot( nondim_time_coef * simulation_data[ :, 0 ],
                 nondim_pressure_coef * simulation_data[ :, i + 1 ],
                 label='WCSPH-BI', linewidth=1, color='k' )

        ax.set_ylabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_xlabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        title = f'Pressure sensor P{ i + 1 }'
        plt.title( title, fontsize=24 )
        output_plot_name = f"results/postprocessing/{case_tag}_pressure_sensor_{ i + 1 }.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_water_level_sensors():

    #Time (s),H1 (m),H2 (m),H3 (m),H4 (m)
    experimental_data_file = resources_dir / "water_level.dat"
    experimental_data = np.genfromtxt( experimental_data_file, delimiter='\t', dtype=np.float64)

    simulation_data_file = results_dir / "sensorsWaterLevel.dat"
    simulation_data = np.genfromtxt( simulation_data_file, delimiter=' ' )

    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = coef * sensor_snapshot_time,
    nondim_time_coef = 1. * 0.002
    # nondim_height_coef = coef
    nondim_height_coef = 1.

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        ax.plot( experimental_data[ :, 0 ],
                 experimental_data[ :, i + 1 ],
                 label='Issa \& Violeau 2006', linewidth=2, color='b'  )
        ax.plot( nondim_time_coef * simulation_data[ :, 0 ],
                 nondim_height_coef * simulation_data[ :, i + 1 ],
                 label='WCSPH-BI', linewidth=2, color='k' )

        ax.set_xlabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_ylabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        title = f'Water level sensor H{ i + 1 }'
        plt.title( title, fontsize=24 )
        output_plot_name = f"results/postprocessing/{case_tag}_water_level_sensor_{ i + 1 }.png"
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
    postproPath = r'./results/postprocessing'
    if not os.path.exists( postproPath ):
        os.makedirs( postproPath )

    # group results
    groupResults.make_data_series( example_dir )

    # plot results from pressure sensors
    plot_pressure_sensors()

    # plot results from pressure sensors
    plot_water_level_sensors()
