#! /usr/bin/env python3

import os
import subprocess
from pathlib import Path
import json
import pandas as pd
import tabulate

# initialize directories
tools_dir = Path(__file__).parent
project_dir = ( tools_dir / ".." / ".." ).resolve()
examples_dir = project_dir / "examples"

cases_list = [
        'RSPH/damBreak2D_RSPH',
        #'WCSPH-BI/damBreak2D_WCSPH-BI',
        #'WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI',
        'WCSPH-DBC/damBreak2D_WCSPH-DBC',
        'WCSPH-DBC/damBreak3D_WCSPH-DBC',
        #'WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC',
        #'WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC'
        ]

cases_list_with_mpi = [
        'WCSPH-DBC/damBreak2D_WCSPH-DBC_distributed',
        'WCSPH-DBC/damBreak3D_WCSPH-DBC_distributed',
        ]

results = []
results_fancy = []
results_returncode = []
computational_time = []

def init( case_dir ):
    args = []
    args += [ case_dir / "init.py" ]
    subprocess.run(args, check=True, cwd=case_dir)

def run( case_dir ):
    args = []
    args += [ case_dir / "run.py" ]

    # run the process and print its output as it is being executed
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          bufsize=1, cwd=case_dir, text=True) as p:
        for line in p.stdout:
            print(line, end="")
    results_returncode.append( subprocess.CalledProcessError( p.returncode, p.args ) )
    results.append( p.returncode )

def parse_tnl_sph_output( case_dir ):
    filename = case_dir / "results" / "time_measurements.json"
    with open( filename ) as f:
        lines = json.load( f )
        json_str = json.dumps( lines )
        timers_dictionary = json.loads( json_str )
        return timers_dictionary[ "total" ]

def run_cases():
    for case in cases_list:
        case_dir = examples_dir / case

        print( f"Initializing case: {case} in {case_dir}." )
        #init( case_dir )
        print( f"Initialization finished." )
        print( f"Executing case: {case} in {case_dir}." )
        run( case_dir )
        print( f"Successfully finished." )

        # get computational time
        computational_time.append( parse_tnl_sph_output( case_dir ) )

def process_results():
    # parse return codes to fancy output
    for entry in results:
        if entry == 0:
            results_fancy.append( '<span style="color:green">Success.</span>' )
        else:
            results_fancy.append( '<span style="color:red">__Failed__</span>' )

def write_results():
    for i in range( len( cases_list ) ):
        print( f"Case: { cases_list[ i ] }\n{ results[ i ] }\nComputational time: { computational_time[ i ] }\n" )

    summary = { 'Cases' : cases_list,
                'Result' : results_fancy,
                'Comp. Time' : computational_time,
                'Ref. Comp. Time' : computational_time }
    summary_df = pd.DataFrame( summary )
    with open('log.md', 'w') as f:
        f.write( summary_df.to_markdown( ) )

if __name__ == "__main__":
    run_cases()
    process_results()
    write_results()
