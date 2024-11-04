import numpy as np
import matplotlib.pyplot as plt

def plot_energy( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )
        print(energy)

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], energy[ :, 2 ] / Epot0, label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], energy[ :, 3 ] / Epot0, label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 ] / Epot0, label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] + energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'simulation time $ t $ [s]')
        ax.set_ylabel( r'time step $ \Delta t$ [s]')
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/time_step_log.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":
    make_data_series()
