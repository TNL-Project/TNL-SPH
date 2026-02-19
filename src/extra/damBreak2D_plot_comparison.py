#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

conf_dir = Path(__file__).parent
example_dir = (conf_dir / ".." / ".." / ".." / "examples/WCSPH-DBC/damBreak2D_WCSPH-DBC").resolve()
resources_dir = (conf_dir / ".." / ".." / ".." / "examples/resources/damBreak2D/damBreak2D_experimentalDataLobovsky2014" ).resolve()

# setup plot parameters
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 40
})

def plot_multiple_pressure_sensors( samples, samples_labels, output_prefix ):

    experimental_data_files = [ "Fig18_peak_event_5_sensors_5.dat",
                                "Fig18_peak_event_5_sensors_4.dat",
                                "Fig18_peak_event_5_sensors_2.dat",
                                "Fig18_peak_event_5_sensors_1.dat" ]


    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = ( norm(g) * H )**0.5 * sensor_snapshot_time,
    nondim_time_coef = ( 9.81 / 0.3 )**0.5 *  0.002
    # nondim_pressure_coef = 1 / ( rho0 * norm(g) * H )
    nondim_pressure_coef = 1. / ( 1000. * 9.81 * 0.3 )

    #data_to_compare_list = []
    #for file in path_to_data_to_compare_with:
    #    print( file )
    #    data_to_compare_list.append(np.genfromtxt( file / "sensorsPressure.dat", delimiter=' ' ))

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "pressure" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='exp.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='_Lobovsky 2014', linewidth=0.7, color='grey', alpha=1 )

        for sample, sample_label, color in zip( samples, samples_labels, [ 'b', 'm', 'k' ] ):
            data_path = example_dir / sample / "sensorsPressure.dat"
            data = np.genfromtxt( data_path, delimiter=' ' )
            mask = nondim_time_coef * data[:, 0] < 8
            ax.plot( nondim_time_coef * data[ mask, 0 ],
                     nondim_pressure_coef * data[ mask, i + 1 ],
                     label=sample_label, linewidth=1, color=color )

        ax.set_xlabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_ylabel( r'$ p/(\rho ||\mathbf{g}|| H)^{1/2}$')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Pressure sensor P{ i + 1 }'
        plt.title( title )
        output_plot_name = str(conf_dir) + f"/{output_prefix}_pressure_sensor_{ i + 1 }_comparison.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_multiple_water_level_sensors( samples, samples_labels, output_prefix ):

    experimental_data_files = [ "Fig16_WaterLevels_H1_2.dat",
                                "Fig16_WaterLevels_H2_3.dat",
                                "Fig16_WaterLevels_H3_3.dat",
                                "Fig16_WaterLevels_H4_4.dat" ]

    #TODO: Read the parameters for normalisation from config.
    # nondim_time_coef = ( norm(g) * H )**0.5 * sensor_snapshot_time,
    nondim_time_coef = ( 9.81 / 0.3 )**0.5 *  0.002
    # nondim_pressure_coef = 1 / H
    nondim_height_coef = 1. / 0.3

    #data_to_compare_list = []
    #for file in path_to_data_to_compare_with:
    #    print( file )
    #    data_to_compare_list.append(np.genfromtxt( file / "sensorsWaterLevel.dat", delimiter=' ' ))

    for i in range( 0, 4 ):
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )

        experimental_data_file = resources_dir / "waterLevel" / experimental_data_files[ i ]
        experimental_data = np.genfromtxt( experimental_data_file )

        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='exp.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( experimental_data[ :, 0 ], experimental_data[ :, 1 ], label='_Lobovsky 2014', linewidth=0.7, color='grey', alpha=1 )

        for sample, sample_label, color in zip( samples, samples_labels, ['b', 'm', 'k'] ):
            data_path = example_dir / sample / "sensorsWaterLevel.dat"
            data = np.genfromtxt( data_path, delimiter=' ' )
            ax.plot( nondim_time_coef * data[ :, 0 ],
                     nondim_height_coef * data[ :, i + 1 ],
                     label=sample_label, linewidth=1.5, color=color )

        ax.set_xlabel( r'$ t( ||\mathbf{g}|| /H)^{1/2} $ ')
        ax.set_ylabel( r'$ h/H $')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        if i == 0:
            leg = ax.legend()
            leg.get_frame().set_edgecolor('k')
        title = f'Water level sensor H{ i + 1 }'
        plt.title( title )
        output_plot_name = str(conf_dir) + f"/{output_prefix}_water_level_sensor_{ i + 1 }_comparison.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":

    # # plot results from pressure sensors
    # plot_multiple_pressure_sensors(
    # [ "results_results_damBreak2D_WCSPH-BI_Verlet",
    #   "results_results_damBreak2D_WCSPH-BI_SymplecticVerlet",
    #   "results_results_damBreak2D_WCSPH-BI_RK4Scheme" ],
    # [ "Verlet",
    #   "Sympl. Verlet",
    #   "RK4" ],
    # 'generateExamplesToTnlSphPaperAppendix/damBreak2D_WCSPH-BI' )
    # # plot results from pressure sensors
    # plot_multiple_water_level_sensors(
    # [ "results_results_damBreak2D_WCSPH-BI_Verlet",
    #   "results_results_damBreak2D_WCSPH-BI_SymplecticVerlet",
    #   "results_results_damBreak2D_WCSPH-BI_RK4Scheme" ],
    # [ "Verlet",
    #   "Sympl. Verlet",
    #   "RK4" ],
    # 'generateExamplesToTnlSphPaperAppendix/damBreak2D_WCSPH-BI' )

    # plot results from pressure sensors
    plot_multiple_pressure_sensors(
    [ "results_damBreak2D_WCSPH-DBC_DBC",
      "results_damBreak2D_WCSPH-DBC_MDBC" ],
    [ "DBC",
      "MDBC" ],
    'generateExamplesToTnlSphPaperAppendix/damBreak2D_WCSPH-DBC' )
    # plot results from pressure sensors
    plot_multiple_water_level_sensors(
    [ "results_damBreak2D_WCSPH-DBC_DBC",
      "results_damBreak2D_WCSPH-DBC_MDBC" ],
    [ "DBC",
      "MDBC" ],
    'generateExamplesToTnlSphPaperAppendix/damBreak2D_WCSPH-DBC' )

    # # plot results from pressure sensors
    # plot_multiple_pressure_sensors(
    # [ "results_damBreak2D_RSPH" ],
    # [ "RSPH" ],
    # 'generateExamplesToTnlSphPaperAppendix/damBreak2D_RSPH' )
    # # plot results from pressure sensors
    # plot_multiple_water_level_sensors(
    # [ "results_damBreak2D_RSPH"],
    # [ "RSPH" ],
    # 'generateExamplesToTnlSphPaperAppendix/damBreak2D_RSPH' )
