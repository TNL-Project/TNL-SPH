import requests
import re

tests = [ "0.005", "0.002", "0.001", "0.00025" ]

filename = "dualSPHysics_0.002_1.out"
#keyword = "CF-Forces"
keywords_search = [ "NL-Limits", "NL-PreSort", "NL-RadixSort", "NL-CellBegin", "NL-SortData", "NL-OutCheck" ]
search_time = 0.

keywords_interaction = [ "CF-Forces", "CF-PreForces", "SU-ComputeStep" ]
interaction_time = 0.

keywords_totalTime = [ "Simulation Runtime" ]
total_time = 0.

keywords_steps = [ "Steps of simulation" ]
total_steps = 0

with open( filename ) as file:
    lines = file.readlines()
    for line in lines:
        for keyword in keywords_search:
            if keyword in line:
                #print( line )
                value = float( re.findall( r"\d+\.\d+", line )[ 0 ] )
                search_time += value;
                #print( "Keyword: ", keyword, " value: ", value )
                break
        for keyword in keywords_interaction:
            if keyword in line:
                #print( line )
                value = float( re.findall( r"\d+\.\d+", line )[ 0 ] )
                interaction_time += value;
                #print( "Keyword: ", keyword, " value: ", value )
                break
        for keyword in keywords_totalTime:
            if keyword in line:
                #print( line )
                value = float( re.findall( r"\d+\.\d+", line )[ 0 ] )
                total_time += value;
                #print( "Keyword: ", keyword, " value: ", value )
                break
        for keyword in keywords_steps:
            if keyword in line:
                #print( line )
                value = int( re.findall( r"\b\d+\b", line )[ 0 ] )
                total_steps += value;
                #print( "Keyword: ", keyword, " value: ", value )
                break

print( "Search : ", search_time )
print( "Interaction time : ", interaction_time )
print( "Total time of simulation (not total runtime!): ", total_time )
print( "Total number of steps: ", total_steps )

import os
import json
import pandas as pd
from pandas.io.json import json_normalize
import matplotlib.pyplot as plt
import numpy as np
import math
from os.path import exists

filename = "tnl-sph_0.002_1.json"
parsed_lines = []
#with open( filename ) as f:
#    lines = f.readlines()
#    for line in lines:
#        parsed_line = json.loads(line)
#        parsed_lines.append( parsed_line )

frames = []
classes = [ 'dualSPHysics', 'TNL::SPH' ]
with open( filename ) as f:
    #multicolumns, df_data = get_multiindex()
    out_idx = 0
    lines = json.load( f )
    json_str = json.dumps( lines )
    resp = json.loads( json_str )

    value = float( resp['interaction'] )
    value2 = float( resp['interaction-average'] )
    data_f = { 'interaction:' :         [ interaction_time,                 float( resp[ 'interaction' ]            )  ],
               'interaction-average:' : [ interaction_time / total_steps,   float( resp[ 'interaction-average' ]    )  ],
               'search:' :              [ search_time,                      float( resp[ 'search' ]                 )  ],
               'search-average:' :      [ search_time / total_steps,        float( resp[ 'search-average' ]         )  ],
               'total:' :               [ total_time ,                      float( resp[ 'total' ]                  )  ],
               'total-average:' :       [ total_time / total_steps,         float( resp[ 'total-average' ]          )  ] }

    new_df = pd.DataFrame( data_f, classes )
    frames.append( new_df )



df = pd.DataFrame( parsed_lines )
print(df)

frames.append( df )
result = pd.concat( frames )
result.to_html( 'time_measurements.html' )
