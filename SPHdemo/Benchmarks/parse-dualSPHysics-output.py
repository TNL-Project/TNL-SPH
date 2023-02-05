import requests
import re

cases = [ "0.005_1", "0.002_1", "0.001_1", "0.0005_1", "0.00025_1" ]
results_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1><br><hr>'

for case in cases:
    filename = "./results/dualSPHysics_" + case + ".out"

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

    filename = "./results/tnl-sph_" + case + ".json"

    parsed_lines = []
    frames = []
    classes = [ 'dualSPHysics', 'TNL::SPH' ]
    with open( filename ) as f:
        #multicolumns, df_data = get_multiindex()
        out_idx = 0
        lines = json.load( f )
        json_str = json.dumps( lines )
        resp = json.loads( json_str )

        tnlSph_interactionTime = float( resp['interaction'] ) + float( resp['integrate'] )
        tnlSph_interactionTime_average = float( resp['interaction-average'] ) + float( resp['integrate-average'] )

        data_f = { 'interaction:' :         [ interaction_time,                 tnlSph_interactionTime                     ],
                   'interaction-average:' : [ interaction_time / total_steps,   tnlSph_interactionTime_average             ],
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

    results_string += '<h2> ' + 'Case ' + case + ' </h2>' + result.to_html(index=True,border=2,justify="center") + '<be><hr>'

results_string +='</center>'

with open("Result.html", 'w') as _file:
    _file.write( results_string )
