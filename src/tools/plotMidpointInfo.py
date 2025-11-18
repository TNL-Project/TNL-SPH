import numpy as np
import matplotlib.pyplot as plt

def plot_midpoint_info( results_dir, residual_trashold=None, iteration_trashold=None ):

        midpoint_data = np.genfromtxt( results_dir / "midpointInfo.dat", delimiter=" " )
        print(midpoint_data)

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( midpoint_data[ :, 1 ], midpoint_data[ :, 2 ], label=r'$ R_{\Delta t}[E] $', linewidth=2, color='b'  )
        if( residual_trashold ): ax.axhline( y=residual_trashold, color='k', linestyle='-')
        ax.set_xlabel( r'$ t $ [s]')
        ax.set_ylabel( r'$ R_{\Delta t}[E] $ [J]')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/midpoint_info_residuals.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( midpoint_data[ :, 1 ], midpoint_data[ :, 3 ], label=r'$ N_{iter} $', linewidth=2, color='b'  )
        if( iteration_trashold ): ax.axhline( y=iteration_trashold, color='k', linestyle='-')
        ax.set_xlabel( r'$ t $ [s]')
        ax.set_ylabel( r'$ N_{iter} $ [-]')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/midpoint_info_iterations.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( midpoint_data[ :, 1 ], midpoint_data[ :, 4 ], label=r'$ mu_{relax} $', linewidth=2, color='b'  )
        ax.set_xlabel( r'$ t $ [s]')
        ax.set_ylabel( r'$ \mu_{relax} $ [-]')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/midpoint_info_relax_coef.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )
