import numpy as np
import matplotlib.pyplot as plt

def plot_vorumetlic_flow_rate(results_dir):
    vfr = np.genfromtxt(results_dir / "volumetricFlowRate.dat", delimiter=" ")

    fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
    ax.plot( vfr[ :, 1 ], vfr[ :, 2 ], label=r'VFR', linewidth=2, color='b'  )
    ax.set_xlabel( r'$ t $ [s]' )
    ax.set_ylabel( r'$ Q $ [m$^{3}\cdot$s$^{-1}$]' )
    ax.grid( color='black', linestyle='--', linewidth=0.5 )
    leg = ax.legend()
    leg.get_frame().set_edgecolor('k')
    output_plot_name = f"results/postprocessing/volumetricFlowRate.png"
    plt.savefig( output_plot_name, bbox_inches='tight' )

    fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
    ax.plot( vfr[ :, 1 ], vfr[ :, 3 ], label=r'VFR', linewidth=2, color='b'  )
    ax.set_xlabel( r'$ t $ [s]' )
    ax.set_ylabel( r'$ Q $ [m$^{3}\cdot$s$^{-1}$]' )
    ax.grid( color='black', linestyle='--', linewidth=0.5 )
    leg = ax.legend()
    leg.get_frame().set_edgecolor('k')
    output_plot_name = f"results/postprocessing/cumulativeVolumetricFlowRate.png"
    plt.savefig( output_plot_name, bbox_inches='tight' )
