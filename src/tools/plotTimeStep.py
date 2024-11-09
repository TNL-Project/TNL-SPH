import numpy as np
import matplotlib.pyplot as plt

def plot_time_step( results_dir ):

        dt_log = np.genfromtxt( results_dir / "timeStep.dat", delimiter=" " )
        print(dt_log)

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( dt_log[ :, 1 ], dt_log[ :, 2 ], label=r'$\Delta t$', linewidth=2, color='b'  )
        ax.set_xlabel( r'$ t $ [s]')
        ax.set_ylabel( r'$ \Delta t$ [s]')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/time_step.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":
    make_data_series()
