from os import listdir
from os.path import isfile, join
from pathlib import Path
from re import findall
import configparser as cfp
import numpy as np
import json

def create_series_file( identifier, snapshot_names, snapshot_timers ):
   n_snapshots = len( snapshot_names )
   files_list = [ {} for _ in range( n_snapshots ) ]

   for i in range( 0, n_snapshots ):
       files_list[ i ][ "name" ] = snapshot_names[ i ]
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

def make_data_series_from_interpolation_planes( config, measuretool_config ):
    for p in range( 1, int( config[ 'main' ][ 'interpolation-planes-count' ] ) + 1 ):
        plane_key = f'plane-{p}-identifier'
        plane_identifier= measuretool_config[ 'main' ][ plane_key ]
        print( plane_identifier )

        # NOTE: This is numpy array only due to
        plane_file_names_subroup = np.array( [ f for f in files if plane_identifier in f ])
        plane_timers_subgroup = np.array( [ findall( r'[\d]*[.][\d]+', f.replace( plane_identifier, "" ) ) for f in plane_file_names_subroup ], dtype=float )

        ## sort the arra based on saved timers
        idx_sort = plane_timers_subgroup.ravel().argsort()
        plane_file_names_subroup = plane_file_names_subroup[ idx_sort ]
        plane_timers_subgroup = plane_timers_subgroup[ idx_sort ]

        create_series_file( plane_identifier, plane_file_names_subroup.ravel(), plane_timers_subgroup.ravel() )

def make_data_series():
    path = Path(__file__).parent
    files = [ f for f in listdir( path ) if isfile( join( path, f ) ) ]

    config = load_case_config()
    measuretool_config = load_case_measuretool_config()

    make_data_series_from_interpolation_planes( config, measuretool_config )

if __name__ == "__main__":
    make_data_series()
