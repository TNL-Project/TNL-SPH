#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

conf_dir = Path(__file__).parent
example_dir = (conf_dir / ".." / ".." / ".." / ".." / "examples/WCSPH-BI/movingSquare2D_WCSPH-BI").resolve()
resources_dir = Path('/home/tomas/Downloads/tmp/SPHERIC_Benchmark') #FIXME

# setup plot parameters
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble" : r"\usepackage{amsfonts}",
  "font.family": "Times",
  "font.serif" : "Times New Roman",
  "font.size"  : 40
})

def plot_set_of_results( samples, samples_labels, output_prefix ):

        referential_data = np.loadtxt(resources_dir / "Force_Re100.dat", skiprows=2)
        """
        t    = referential_data[:, 0]
        Cd_p = referential_data[:, 1]
        Cd_v = referential_data[:, 2]
        """
        # subsample the numerical data
        subsample = 30

        # plot pressure drag
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 1 ], label=r'ref.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 1 ], label=r'_ref.', linewidth=0.7, color='grey', alpha=1 )

        for sample, sample_label, color in zip( samples, samples_labels, ['b', 'm', 'k']):
            data_path = example_dir / sample / "force.dat"
            data = np.genfromtxt( data_path, delimiter=" " )
            ax.plot( data[ :, 1 ][::subsample], 2 * data[ :, 2 ][::subsample], label=sample_label, linewidth=1, color=color, alpha=0.9  )

        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ c_{D,\; pressure} $ [-]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        #leg = ax.legend()
        #leg.get_frame().set_edgecolor('k')
        output_plot_name = str(conf_dir) + f"/{output_prefix}_cd_pressure_comparison.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 2 ], label=r'$ref.$.', linewidth=8, color='dimgrey', alpha=0.6 )
        ax.plot( referential_data[ :, 0 ], referential_data[ :, 2 ], label=r'_ref.', linewidth=0.7, color='grey', alpha=1 )

        for sample, sample_label, color in zip( samples, samples_labels, ['b', 'm', 'k']):
            data_path = example_dir / sample / "force.dat"
            data = np.genfromtxt( data_path, delimiter=" " )
            ax.plot( data[ :, 1 ][::subsample], 2 * data[ :, 3 ][::subsample], label=sample_label, linewidth=1.5, color=color, alpha=1 )

        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$ c_{D,\; viscous} $ [-]' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = str(conf_dir) + f"/{output_prefix}_cd_visous_comparison.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":

    
    plot_set_of_results( [ "results_spheric06-dp0.04-BIConsistent-MVT-Verlet-cs50",
                           "results_spheric06-dp0.02-BIConsistent-MVT-Verlet-cs50",
                           "results_spheric06-dp0.0125-BIConsistent-MVT-Verlet-cs50" ],
                         [  r'$\Delta x = 0.04$', 
                            r'$\Delta x = 0.02$',
                            r'$\Delta x = 0.0125$' ],
                         "BIConsistent-MVT-Verlet-cs50" )

    plot_set_of_results( [ "results_spheric06-dp0.04-BIConsistent-MVT-Verlet-PST-DensityFilter-n5-cs50",
                           "results_spheric06-dp0.02-BIConsistent-MVT-Verlet-PST-DensityFilter-n5-cs50",
                           "results_spheric06-dp0.0125-BIConsistent-MVT-Verlet-PST-DensityFilter-n5-cs50" ],
                         [  r'$\Delta x = 0.04$', 
                            r'$\Delta x = 0.02$',
                            r'$\Delta x = 0.0125$' ],
                         "BIConsistent-MVT-Verlet-PST-DensityFilter-n5-cs50" )
