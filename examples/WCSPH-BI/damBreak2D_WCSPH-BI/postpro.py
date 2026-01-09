#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Include globaly used postpro scripts
import sys
sys.path.append('../../../src/tools')

#TODO: Rename the default sensors labels providet by Lobovsky et. al.
case_tag = "damBreak2D_WCSPH-BI"
example_dir = Path(__file__).parent
resources_dir = (example_dir / ".." / ".." / "resources" / "damBreak2D" / "damBreak2D_experimentalDataLobovsky2014" ).resolve()
results_dir = ( example_dir / "results--dp_0.001" ).resolve()

# setup plot parameters
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 40
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

        mask = nondim_time_coef * simulation_data[:, 0] < 8
        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='exp', linewidth=2, color='b'  )
        ax.plot( nondim_time_coef * simulation_data[ mask, 0 ],
                 nondim_pressure_coef * simulation_data[ mask, i + 1 ],
                 label='WCSPH-BI', linewidth=1, color='k' )

        ax.set_xlabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_ylabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Pressure sensor P{ i + 1 }'
        plt.title( title )
        output_plot_name = str(results_dir) + f"/postprocessing/{case_tag}_pressure_sensor_{ i + 1 }.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_multiple_pressure_sensors( path_to_data_to_compare_with, data_to_compare_labels, this_data_label ):

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

    data_to_compare_list = []
    for file in path_to_data_to_compare_with:
        print( file )
        data_to_compare_list.append(np.genfromtxt( file / "sensorsPressure.dat", delimiter=' ' ))

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "pressure" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='exp.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='_Lobovsky 2014', linewidth=0.7, color='grey', alpha=1 )

        for data, datalabel, color in zip( data_to_compare_list, data_to_compare_labels, [ 'b', 'm' ] ):
            mask = nondim_time_coef * data[:, 0] < 8
            ax.plot( nondim_time_coef * data[ mask, 0 ],
                     nondim_pressure_coef * data[ mask, i + 1 ],
                     label=datalabel, linewidth=1, color=color )

        mask = nondim_time_coef * simulation_data[:, 0] < 8
        ax.plot( nondim_time_coef * simulation_data[ mask, 0 ],
                 nondim_pressure_coef * simulation_data[ mask, i + 1 ],
                 label=this_data_label, linewidth=1.5, color='k' )

        ax.set_xlabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_ylabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Pressure sensor P{ i + 1 }'
        plt.title( title )
        output_plot_name = str(results_dir) + f"/postprocessing/{case_tag}_pressure_sensor_{ i + 1 }_comparison.png"
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
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Water level sensor H{ i + 1 }'
        plt.title( title )
        output_plot_name = str(results_dir) + f"/postprocessing/{case_tag}_water_level_sensor_{ i + 1 }.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_multiple_water_level_sensors( path_to_data_to_compare_with, data_to_compare_labels, this_data_label ):
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

    data_to_compare_list = []
    for file in path_to_data_to_compare_with:
        print( file )
        data_to_compare_list.append(np.genfromtxt( file / "sensorsWaterLevel.dat", delimiter=' ' ))

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "waterLevel" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='exp.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='_Lobovsky 2014', linewidth=0.7, color='grey', alpha=1 )

        for data, datalabel, color in zip( data_to_compare_list, data_to_compare_labels, [ 'b', 'm' ] ):
            ax.plot( nondim_time_coef * data[ :, 0 ],
                     nondim_height_coef * data[ :, i + 1 ],
                     label=datalabel, linewidth=1.5, color=color )

        ax.plot( nondim_time_coef * simulation_data[ :, 0 ],
                nondim_height_coef * simulation_data[ :, i + 1 ],
                label=datalabel, linewidth=1.5, color='k' )

        ax.set_ylabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_xlabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Water level sensor H{ i + 1 }'
        plt.title( title )
        output_plot_name = str(results_dir) + f"/postprocessing/{case_tag}_water_level_sensor_{ i + 1 }_comparison.png"
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
    #plot_multiple_pressure_sensors( [ example_dir / "results--dp_0.005", example_dir / "results--dp_0.002"], [ r'$\Delta x = 0.005$', r'$\Delta x = 0.002$' ], r'$\Delta x = 0.001$' )

    # plot results from pressure sensors
    plot_water_level_sensors()
    #plot_multiple_water_level_sensors( [ example_dir / "results--dp_0.005", example_dir / "results--dp_0.002"], [ r'$\Delta x = 0.005$', r'$\Delta x = 0.002$' ], r'$\Delta x = 0.001$' )

    # Global postprocessing tools: group results
    import writeParaviewSeriesFile
    writeParaviewSeriesFile.generate_series( results_dir )

    # Global postprocessing tools: plot midpoint iterations info
    import plotMidpointInfo
    plotMidpointInfo.plot_midpoint_info( results_dir, residual_trashold =.33166e-05, iteration_trashold = 10 )

    # Global postprocessing tools: plot time step log and energy
    import plotEnergy
    #plotTimeStep.plot_time_step( results_dir )
    plotEnergy.plot_energy( results_dir, Epot0 = 264.87 )
    plotEnergy.plot_energy_snapshots( results_dir, Epot0 = 264.87 )
