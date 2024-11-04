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
case_tag = "damBreak2D_WCSPH-BI"
example_dir = Path(__file__).parent
resources_dir = (example_dir / ".." / ".." / "resources" / "damBreak2D" / "damBreak2D_experimentalDataLobovsky2014" ).resolve()
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

    experimental_data_files = [ "Fig18_peak_event_5_sensors_5.dat",
                                "Fig18_peak_event_5_sensors_4.dat",
                                "Fig18_peak_event_5_sensors_2.dat",
                                "Fig18_peak_event_5_sensors_1.dat" ]

    simulation_data_file = results_dir / "sensorsPressure.dat"
    simulation_data = np.genfromtxt( simulation_data_file, delimiter=' ' )

    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = ( norm(g) * H )**0.5 * sensor_snapshot_time,
    nondim_time_coef = ( 9.81 / 0.3 )**0.5 *  0.002
    # nondim_pressure_coef = 1 / ( rho0 * norm(g) * H )
    nondim_pressure_coef = 1. / ( 1000. * 9.81 * 0.3 )

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "pressure" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='Lobovsky 2014', linewidth=2, color='b'  )
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
    experimental_data_files = [ "Fig16_WaterLevels_H1_2.dat",
                                "Fig16_WaterLevels_H2_3.dat",
                                "Fig16_WaterLevels_H3_3.dat",
                                "Fig16_WaterLevels_H4_4.dat" ]

    simulation_data_file = results_dir / "sensorsWaterLevel.dat"
    simulation_data = np.genfromtxt( simulation_data_file, delimiter=' ' )

    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = ( norm(g) * H )**0.5 * sensor_snapshot_time,
    nondim_time_coef = ( 9.81 / 0.3 )**0.5 *  0.002
    # nondim_pressure_coef = 1 / H
    nondim_height_coef = 1. / 0.3

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "waterLevel" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='Lobovsky 2014', linewidth=2, color='b'  )
        ax.plot( nondim_time_coef * simulation_data[ :, 0 ],
                nondim_height_coef * simulation_data[ :, i + 1 ],
                label='WCSPH-BI', linewidth=2, color='k' )

        ax.set_ylabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_xlabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
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

    # plot results from pressure sensors
    plot_pressure_sensors()

    # plot results from pressure sensors
    plot_water_level_sensors()

    # Global postprocessing tools: plot time step log and energy
    plotTimeStep.plot_time_step( results_dir )
    plotEnergy.plot_energy( results_dir, Epot0 = 264.87 )
