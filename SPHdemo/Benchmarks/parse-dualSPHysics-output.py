import re
import os
import json
import pandas as pd
from pandas.io.json import json_normalize
import matplotlib.pyplot as plt
import numpy as np
import math
from os.path import exists
from json2html import *

#cases = [ "0.005_1", "0.002_1", "0.001_1", "0.0005_1", "0.00025_1" ]
#folder = "results/"

cases = [ "0.002_1" ]
folder = "results_local-test/"

#results_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1><br><hr>'
results_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1>'

#device metadata
deviceFile = folder + "tnl-sph_" + cases[ 0 ] + ".device_metadata.json"
with open( deviceFile ) as f:
    lines = json.load( f )
    json_str = json.dumps( lines )
    resp = json.loads( json_str )

    deviceString = json2html.convert(json = resp)
    #deviceString = json2html.convert(json = resp)
    results_string += deviceString
results_string += "<br><hr>"


for case in cases:
    #filename = "./results/dualSPHysics_" + case + ".out"
    filename = folder + "dualSPHysics_" + case + ".out"

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


    #filename = "./results/tnl-sph_" + case + ".json"
    filename = folder + "tnl-sph_" + case + ".json"

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
    #result.to_html( 'time_measurements.html' )

    #Case details
    filename_caseMetadata = folder + "tnl-sph_" + case + ".case_metadata.json"
    with open( filename_caseMetadata ) as file_metadata:
        caseMetadata_lines = json.load( file_metadata )
        caseMetadata_json_str = json.dumps( caseMetadata_lines )
        caseMetadata_json = json.loads( caseMetadata_json_str )

        detail_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1>'
        detail_string += json2html.convert( json = caseMetadata_json ) + "<br><hr>"
        detail_string += json2html.convert( json = resp ) + "<br><hr>"
        detail_string +='</center>'

        outputFileName = folder + "case_detail.html"
        with open( outputFileName, 'w') as _file:
            _file.write( detail_string )

    detail_string = ' <a href=\"case_detail.html\"> Details </a>'
    #results_string += '<h2> ' + 'Case ' + case + ' </h2>' + result.to_html(index=True,border=2,justify="center") + '<be><hr>'
    results_string += '<h2> ' + 'Case ' + case + ' </h2>' + detail_string + result.to_html(index=True,border=2,justify="center") + '<be><hr>'


results_string +='</center>'

outputFileName = folder + "result.html"
with open( outputFileName, 'w') as _file:
    _file.write( results_string )
