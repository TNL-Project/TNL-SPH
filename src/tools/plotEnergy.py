import numpy as np
import matplotlib.pyplot as plt

def plot_energy( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], energy[ :, 2 ] / Epot0, label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], ( Epot0 + energy[ :, 3 ] ) / Epot0, label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 ] / Epot0, label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], ( Epot0 + energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$\Delta E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$\Delta E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy_change.png"
        print( "Plot energy: {output_plot_name}" )
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_energy_snapshots( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )
        offset = 3

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], energy[ :, 2 + offset ] / Epot0, label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], ( energy[ :, 3 + offset ] ) / Epot0, label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 + offset ] / Epot0, label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 + offset ] + energy[ :, 3 + offset ] + energy[ :, 4 + offset ] ) / Epot0, label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy_snapshots.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

        """
        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], ( energy[ :, 2 ] +  energy[ :, 3 ] + energy[ :, 4 ] ) / Epot0, label=r'$\Delta E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$\Delta E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy_change.png"
        print( "Plot energy: {output_plot_name}" )
        plt.savefig( output_plot_name, bbox_inches='tight' )
        """

def plot_not_normalized_energy( results_dir ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], energy[ :, 2 ], label=r'$E_{kin}$', linewidth=2, color='b'  )
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
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$\Delta E(t) [J]$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"{results_dir}/postprocessing/energy_change.png"
        print( "Plot energy: {output_plot_name}" )
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_not_normalized_energy_snapshots( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energy.dat", delimiter=" " )
        offset = 3

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        ax.plot( energy[ :, 1 ], energy[ :, 2 + offset ], label=r'$E_{kin}$', linewidth=2, color='b'  )
        ax.plot( energy[ :, 1 ], energy[ :, 3 + offset ], label=r'$E_{pot}$', linewidth=2, color='g'  )
        ax.plot( energy[ :, 1 ], energy[ :, 4 + offset ], label=r'$E_{comp}$', linewidth=2, color='r'  )
        ax.plot( energy[ :, 1 ], energy[ :, 2 + offset ] + energy[ :, 3 + offset ] + energy[ :, 4 + offset ], label=r'$E_{tot}$', linewidth=2, color='k'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energy_snapshots.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

def plot_not_normalized_open_boundary_energy_snapshots( results_dir, Epot0 = 1 ):

        energy = np.genfromtxt( results_dir / "energyOpenBoundary.dat", delimiter=" " )
        offset = 3

        totalEnergyIn = np.sum( energy[ :, 2 ] + energy[ :, 3 ] + energy[ :, 4 ] )
        totalMeanEnergyIn = np.mean( energy[ :, 2 ] + energy[ :, 3 ] + energy[ :, 4 ] )
        totalEnergyOut = np.sum( -( energy[ :, 5 ] + energy[ :, 6 ] + energy[ :, 7 ] ) )
        totalMeanEnergyOut = np.mean( -( energy[ :, 5 ] + energy[ :, 6 ] + energy[ :, 7 ] ) )
        print( f"Total energy in: {totalEnergyIn} [J], total energy out: {totalEnergyOut} [J]." )
        print( f"Total mean energy per step in: {totalMeanEnergyIn} [J], total mean energy per step out: {totalMeanEnergyOut} [J]." )

        fig, ax = plt.subplots( 1, 1, figsize=( 11, 8 ) )
        #ax.plot( energy[ :, 1 ], energy[ :, 2 ], label=r'$E_{kin-in}$', linewidth=2, color='b'  )
        #ax.plot( energy[ :, 1 ], energy[ :, 3 ], label=r'$E_{pot-in}$', linewidth=2, color='g'  )
        #ax.plot( energy[ :, 1 ], energy[ :, 4 ], label=r'$E_{comp-in}$', linewidth=2, color='r'  )
        ax.scatter( energy[ :, 1 ], energy[ :, 2 ] + energy[ :, 3 ] + energy[ :, 4 ], label=r'$E_{tot-in}$', linewidth=2, marker="^", color='b'  )
        #ax.plot( energy[ :, 1 ], -energy[ :, 5 ], label=r'$E_{kin-out}$', linewidth=2, color='b'  )
        #ax.plot( energy[ :, 1 ], -energy[ :, 6 ], label=r'$E_{pot-out}$', linewidth=2, color='g'  )
        #ax.plot( energy[ :, 1 ], -energy[ :, 7 ], label=r'$E_{comp-out}$', linewidth=2, color='r'  )
        ax.scatter( energy[ :, 1 ], -( energy[ :, 5 ] + energy[ :, 6 ] + energy[ :, 7 ] ), label=r'$E_{tot-out}$', linewidth=2, marker="v", color='r'  )
        ax.set_xlabel( r'$ t $ [s]' )
        ax.set_ylabel( r'$E/E_p(t=0)$' )
        ax.grid( color='black', linestyle='--', linewidth=0.5 )
        leg = ax.legend()
        leg.get_frame().set_edgecolor('k')
        output_plot_name = f"results/postprocessing/energyOpenBoundary_snapshots.png"
        plt.savefig( output_plot_name, bbox_inches='tight' )

if __name__ == "__main__":
    make_data_series()
