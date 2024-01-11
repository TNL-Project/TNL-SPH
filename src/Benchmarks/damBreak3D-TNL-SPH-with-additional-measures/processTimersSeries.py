import numpy as np
import matplotlib.pyplot as plt

# PyPlot settings
plt.rc('font', size=15)             # controls default text sizes
plt.rc('axes', titlesize=20)        # fontsize of the axes title
plt.rc('axes', labelsize=25)        # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)       # fontsize of the tick labels
plt.rc('ytick', labelsize=20)       # fontsize of the tick labels
plt.rc('legend', fontsize=15)       # legend fontsize
plt.rc('figure', titlesize=25)      # fontsize of the figure title


# print time series with time step
def printTimeSeries( timerName, title = '' ):
    inputFileName = 'results/timer_array-' + timerName + '_array.dat'
    timer_array = np.genfromtxt( inputFileName, delimiter=',' )

    fig = plt.figure( figsize=(18,5) )
    plt.plot( timer_array[:-1, 0], timer_array[:-1, 1], color='navy', linewidth=2  )
    plt.xlabel( 'step [-]', fontsize=24 )
    plt.ylabel( 'time per step [s]', fontsize=24 )
    plt.grid( which="major", color='black', linestyle='--', linewidth=0.5 )
    plt.grid( which="minor", color='grey', linestyle='--', linewidth=0.3 )
    plt.minorticks_on()
    plt.title( title, fontsize=26)
    outputFigureName = 'results/fig_timer_array-' + timerName +'_array.png'
    fig.savefig( outputFigureName, bbox_inches='tight' )

if __name__ == "__main__":

    printTimeSeries( 'time_interact', 'Comp. time per step - interaction' )
    printTimeSeries( 'time_search', 'Comp. time per step - search' )
    printTimeSeries( 'time_search_sort', 'Comp. time per step - search: sort' )
    printTimeSeries( 'time_total', 'Comp. time per step - total' )
