import requests
import re

tests = [ "0.005", "0.002", "0.001", "0.00025" ]

filename = "results.out"
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
