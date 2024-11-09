import numpy as np
import matplotlib.pyplot as plt

def plot_energy( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )
        print(energy)

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], -energy[ :, 2 ] / Epot0, label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], 1 - energy[ :, 3 ] / Epot0, label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 ] / Epot0, label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], 1 + ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$\Delta E_{tot}$', linewidth=2, color='k'  )
        #ax.plot( energy[ :, 1 ], energy[ :, 4 ] / Epot0, label=r'$\Delta E_{comp}$', linewidth=2, color='r'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$\Delta E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy_change.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_not_normalized_energy( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )
        print(energy)

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], -energy[ :, 2 ], label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], energy[ :, 3 ], label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 ], label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ), label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E(t) [J]$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/energy.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ), label=r'$\Delta E_{tot}$', linewidth=2, color='k'  )
        #ax.plot( energy[ :, 1 ], energy[ :, 4 ] / Epot0, label=r'$\Delta E_{comp}$', linewidth=2, color='r'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$\Delta E(t) [J]$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/energy_change.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":
    make_data_series()
