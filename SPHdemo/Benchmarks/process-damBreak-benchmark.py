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

cases = [ "0.005_1", "0.002_1", "0.001_1", "0.0005_1", "0.00025_1" ]
folder = "results_A40galdor/"

#cases = [ "0.002_1" ]
#folder = "/"

def getCaseDetails( case, TNLSPHTimers, dualSPHTimers ):
    #tnl-sph details
    caseMetadataFileName = folder + "tnl-sph_" + case + ".case_metadata.json"
    with open( caseMetadataFileName ) as file_metadata:
        caseMetadata_lines = json.load( file_metadata )
        caseMetadata_json_str = json.dumps( caseMetadata_lines )
        caseMetadata_json = json.loads( caseMetadata_json_str )

        detail_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1>'
        detail_string += json2html.convert( json = caseMetadata_json ) + "<br><hr>"
        detail_string += json2html.convert( json = TNLSPHTimers ) + "<br><hr>"
        detail_string += json2html.convert( json = dualSPHTimers ) + "<br><hr>"
        detail_string +='</center>'

        caseMetadataResultName = "case_" + case + "_detail.html"
        caseMetadataResultFileName = folder + caseMetadataResultName
        with open( caseMetadataResultFileName, 'w') as _file:
            _file.write( detail_string )

        return caseMetadataResultName

def parseDualSPHysicsOutput( case ):
    filename = folder + "dualSPHysics_" + case + ".out"

    timersDictionary = {
        'VA-Init' : 0,
        'NL-Limits' : 0,
        'NL-PreSort' : 0,
        'NL-RadixSort' : 0,
        'NL-CellBegin' : 0,
        'NL-SortData' : 0,
        'NL-OutCheck' : 0,
        'CF-PreForces' : 0,
        'CF-Forces' : 0,
        'SU-Shifting' : 0,
        'SU-ComputeStep' : 0,
        'SU-Floating' : 0,
        'SU-Motion' : 0,
        'SU-Periodic' : 0,
        'SU-ResizeNp' : 0,
        'SU-DownData' : 0,
        'SU-SavePart' : 0,
        'SU-Chrono' : 0,
        'SU-BoundCorr' : 0,
        'SU-InOut' : 0,
        'Steps of simulation' : 0,
        'Steps per second' : 0,
        'Total Runtime' : 0,
        'Simulation Runtime' : 0 }

    with open( filename ) as file:
        lines = file.readlines()
        for line in lines:
            for key, value in timersDictionary.items():
                if key in line:

                    #print(key)
                    #print( re.findall( r"\d+\.\d+", line ) ) if \
                    #        re.findall( r"\d+\.\d+", line ) else print( re.findall( r"\b\d+\b", line ) )

                    parsedValue = float( re.findall( r"\d+\.\d+", line )[ 0 ] ) if \
                            re.findall( r"\d+\.\d+", line ) else  int( re.findall( r"\b\d+\b", line )[ 0 ] )
                    timersDictionary[ key ] = parsedValue

                    break

    #print( "DualSPHysics parsed timers : ", timersDictionary )
    return( timersDictionary )

def parseTNLSPHOutput( case ):
    filename = folder + "tnl-sph_" + case + ".json"
    with open( filename ) as f:
        out_idx = 0
        lines = json.load( f )
        json_str = json.dumps( lines )
        timersDictionary = json.loads( json_str )

        return timersDictionary

results_string = '<center>' +'<h1> SPH damBreak2D benchmark </h1>'

#device metadata
deviceMetadataFileName = folder + "tnl-sph_" + cases[ 0 ] + ".device_metadata.json"
with open( deviceMetadataFileName ) as f:
    deviceMetadata_lines = json.load( f )
    deviceMetadata_json_str = json.dumps( deviceMetadata_lines )
    deviceMetadata_json = json.loads( deviceMetadata_json_str )

    deviceString = json2html.convert(json = deviceMetadata_json)
    results_string += deviceString

results_string += "<br><hr>"

for case in cases:

    dualSPHTimers = parseDualSPHysicsOutput( case )
    TNLSPHTimers = parseTNLSPHOutput( case )

    frames = []
    classes = [ 'dualSPHysics', 'TNL::SPH' ]

    tnlSph_interactionTime = float( TNLSPHTimers['interaction'] ) + float( TNLSPHTimers['integrate'] )
    tnlSph_interactionTime_average = float( TNLSPHTimers['interaction-average'] ) + float( TNLSPHTimers['integrate-average'] )
    tnlSph_searchTime = float( TNLSPHTimers[ 'search' ] )
    tnlSph_searchTime_average = float( TNLSPHTimers[ 'search-average' ] )
    tnlSph_totalTime = float( TNLSPHTimers[ 'total' ] )
    tnlSph_totalTime_average = float( TNLSPHTimers[ 'total-average' ] )

    dualSph_totalSteps = int( dualSPHTimers['Steps of simulation'] )
    dualSph_interactionTime = float( dualSPHTimers['CF-PreForces'] ) + \
                              float( dualSPHTimers['CF-Forces'] ) + \
                              float( dualSPHTimers['SU-ComputeStep'] )
    dualSph_interactionTime_average = dualSph_interactionTime / dualSph_totalSteps
    dualSph_searchTime = float( dualSPHTimers[ 'NL-Limits' ] ) + \
                         float( dualSPHTimers[ 'NL-PreSort' ] ) + \
                         float( dualSPHTimers[ 'NL-RadixSort' ] ) + \
                         float( dualSPHTimers[ 'NL-CellBegin' ] ) + \
                         float( dualSPHTimers[ 'NL-SortData' ] ) + \
                         float( dualSPHTimers[ 'NL-OutCheck' ] )
    dualSph_searchTime_average = dualSph_searchTime / dualSph_totalSteps
    dualSph_totalTime = float( dualSPHTimers['Simulation Runtime'] )
    dualSph_totalTime_average = dualSph_totalTime / dualSph_totalSteps

    data_f = { 'interaction:' : [ dualSph_interactionTime, tnlSph_interactionTime ],
               'interaction-average:' : [ dualSph_interactionTime_average, tnlSph_interactionTime_average ],
               'search:' :  [ dualSph_searchTime,  tnlSph_searchTime ],
               'search-average:' : [ dualSph_searchTime_average, tnlSph_searchTime_average ],
               'total:' : [ dualSph_totalTime, tnlSph_totalTime ],
               'total-average:' : [ dualSph_totalTime_average, tnlSph_totalTime_average ] }

    new_df = pd.DataFrame( data_f, classes )
    frames.append( new_df )

    result = pd.concat( frames )
    caseDetail = getCaseDetails( case, TNLSPHTimers, dualSPHTimers )
    detail_string = ' <a href=\"'+ caseDetail + '\"> Details </a>'
    results_string += '<h2> ' + 'Case ' + case + ' </h2>' + detail_string + result.to_html(index=True,border=2,justify="center") + '<be><hr>'

results_string +='</center>'

outputFileName = folder + "result.html"
with open( outputFileName, 'w') as _file:
    _file.write( results_string )
