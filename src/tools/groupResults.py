import os
from os import listdir, makedirs
from os.path import isfile, join, exists
from pathlib import Path
from re import findall
import configparser as cfp
import numpy as np
import json

def create_series_file( identifier, results_dir, snapshot_names, snapshot_timers ):
   n_snapshots = len( snapshot_names )
   files_list = [ {} for _ in range( n_snapshots ) ]

   for i in range( 0, n_snapshots ):
       #snapshot_path = str( results_dir / "results_data" / snapshot_names[ i ])
       snapshot_path = "results_data/" + ( snapshot_names[ i ] )
       print( snapshot_path )
       files_list[ i ][ "name" ] = snapshot_path
       files_list[ i ][ "time" ] = snapshot_timers[ i ]

   series_file_dict = {
        "file-series-version" : "1.0",
        "files" : files_list
        }

   series_output_name = "results/" + identifier + ".vtk.series"
   # due to list of files, pprint looks better inside the file
   with open( series_output_name, "w") as file:
       json.dump( series_file_dict, file, indent = 4 )

# Read config files to know what we need to group
def load_case_config( path = "" ):
    config_file_path =  path + 'sources/config.ini'
    config = cfp.ConfigParser()
    # configparser is unable to read file without headers
    with open( config_file_path ) as stream:
        config.read_string( "[main]\n" + stream.read() )

    return config

def load_case_measuretool_config( path = "" ):
    measuretool_config_file_path = path + 'sources/config-measuretool.ini'
    measuretool_config = cfp.ConfigParser()
    # configparser is unable to read file without headers
    with open( measuretool_config_file_path ) as stream:
        measuretool_config.read_string( "[main]\n" + stream.read() )

    return measuretool_config

def select_data_and_build_series( files, results_dir, output_identifier ):
        # NOTE: This is numpy array only due to sorting
        plane_file_names_subroup = np.array( [ f for f in files if output_identifier in f ] )
        if np.size( plane_file_names_subroup ) == 0:
            return

        # NOTE: Version to find only decimals: r'[\d]*[.][\d]+'
        plane_timers_subgroup = np.array( [ findall( r'[-+]?(?:\d*\.*\d+)', f.replace( output_identifier, "" ) ) for f in plane_file_names_subroup ], dtype=float )

        # sort the array based on saved timers
        idx_sort = plane_timers_subgroup.ravel().argsort()
        plane_file_names_subroup = plane_file_names_subroup[ idx_sort ]
        plane_timers_subgroup = plane_timers_subgroup[ idx_sort ]

        create_series_file( output_identifier, results_dir, plane_file_names_subroup.ravel(), plane_timers_subgroup.ravel() )

        # move source data folder with raw files
        for plane_file in plane_file_names_subroup:
            os.rename( results_dir / Path( plane_file ),  results_dir / "results_data" / Path( plane_file ) )

def make_data_series_from_interpolation_planes( files, results_dir, config, measuretool_config ):
    for p in range( 1, int( config[ 'main' ][ 'interpolation-planes-count' ] ) + 1 ):
        plane_key = f'plane-{p}-identifier'
        plane_identifier= measuretool_config[ 'main' ][ plane_key ]
        select_data_and_build_series( files, results_dir, plane_identifier )

def make_data_series_from_particles( files, results_dir, config ):
    output_identifiers = [ 'fluid', 'boundary', 'inlet', 'outlet', 'grid' ]
    for output_identifier in output_identifiers:
        select_data_and_build_series( files, results_dir, output_identifier )

def make_data_series( example_dir ):
    results_dir = example_dir / "results"
    files = [ f for f in listdir( results_dir ) if isfile( join( results_dir, f ) ) ]

    results_raw_data_dir = results_dir / "results_data"
    if not os.path.exists( results_raw_data_dir ):
        os.makedirs( results_raw_data_dir )

    config = load_case_config()
    measuretool_config = load_case_measuretool_config()

    make_data_series_from_particles( files, results_dir, config )
    make_data_series_from_interpolation_planes( files, results_dir, config, measuretool_config )

if __name__ == "__main__":
    make_data_series()
