import re
import numpy as np
import matplotlib.pyplot as plt

def plot_vorumetlic_flow_rate(results_dir, vfr_file = "volumetricFlowRate.dat"):
    vfr = np.genfromtxt(results_dir / vfr_file, delimiter=" ")

    fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
    ax.plot( vfr[ :, 1 ], vfr[ :, 2 ], label=r'VFR', linewidth=2, color='b'  )
    ax.set_xlabel( r'$ t $ [s]' )
    ax.set_ylabel( r'$ Q $ [m$^{3}\cdot$s$^{-1}$]' )
    ax.grid( color='black', linestyle='--', linewidth=0.5 )
    leg = ax.legend()
    leg.get_frame().set_edgecolor('k')
    output_plot_name = f"results/postprocessing/{re.sub('.dat$', '.png', vfr_file)}"
    plt.savefig( output_plot_name, bbox_inches='tight' )

    fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
    ax.plot( vfr[ :, 1 ], vfr[ :, 3 ], label=r'VFR', linewidth=2, color='b'  )
    ax.set_xlabel( r'$ t $ [s]' )
    ax.set_ylabel( r'$ Q $ [m$^{3}\cdot$s$^{-1}$]' )
    ax.grid( color='black', linestyle='--', linewidth=0.5 )
    leg = ax.legend()
    leg.get_frame().set_edgecolor('k')
    output_plot_name = f"results/postprocessing/cumulative-{re.sub('.dat$', '.png', vfr_file)}"
    plt.savefig( output_plot_name, bbox_inches='tight' )
