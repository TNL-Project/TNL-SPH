#! /usr/bin/env python3

from pathlib import Path
import subprocess
import json
import pandas as pd
from time import strftime, gmtime

# list of tests
import configurationsToTest
conf_list = configurationsToTest.test_configurations

# initialize directories
tools_dir = Path(__file__).parent
project_dir = ( tools_dir / ".." / ".." ).resolve()
examples_dir = project_dir / "examples"

# storage arrays
results = []
results_fancy = []
results_returncode = []
computational_time = []
referential_computational_time = []
computational_time_difference_formatted = []
cases_list = []
tests_passed_list = []
tests_total_list = []

def init( case_dir, conf ):
    args = []
    args += [ case_dir / "init.py" ]
    for key, value in conf.items():
        if key not in [ "case", "evaluation-function" ]:
            args += [ f"--{key}", str( value ) ]

    print( args )
    subprocess.run( args,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                    cwd=case_dir )

def run( case_dir ):
    args = []
    args += [ case_dir / "run.py" ]

    # run the process and print its output as it is being executed
    #with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                      bufsize=1, cwd=case_dir, text=True) as p:
    #    for line in p.stdout:
    #        print(line, end="")

    p =  subprocess.run( args,
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.STDOUT,
                         cwd=case_dir,
                         text=True )
    results_returncode.append( subprocess.CalledProcessError( p.returncode, p.args ) )
    results.append( p.returncode )

def evaluate_test_metrics( case_dir, conf ):
    evaluation_function = conf[ "evaluation-function" ]
    tests_passed, tests_total = evaluation_function( case_dir )
    tests_passed_list.append( tests_passed )
    tests_total_list.append( tests_total )

def parse_tnl_sph_output( case_dir ):
    filename = case_dir / "results" / "time_measurements.json"
    with open( filename ) as f:
        lines = json.load( f )
        json_str = json.dumps( lines )
        timers_dictionary = json.loads( json_str )
        return float( timers_dictionary[ "total" ] )

def run_cases():
    for conf in conf_list:
        case = conf[ "case" ]
        cases_list.append( case )
        case_dir = examples_dir / case

        print( f"Initializing case: {case} in {case_dir}." )
        init( case_dir, conf )
        print( f"Initialization finished." )
        print( f"Executing case: {case} in {case_dir}." )
        run( case_dir )
        print( f"Execution finished with return code: {results[ -1 ]}." )
        evaluate_test_metrics( case_dir, conf )

        # get computational time
        computational_time.append( parse_tnl_sph_output( case_dir ) )

def process_results( gpu_type ):
    # parse return codes to fancy output
    for entry in results:
        if entry == 0:
            results_fancy.append( '<span style="color:green">__Success__</span>' )
        else:
            results_fancy.append( '<span style="color:red">__Failed__</span>' )

    # compute timers
    referential_timers_filename = f"referentialComputationalTimes/referential_comp_times_{gpu_type}.json"
    with open( referential_timers_filename ) as f:
        lines = json.load( f )
        json_str = json.dumps( lines )
        timers_dictionary = json.loads( json_str )

        for key, value in timers_dictionary.items():
            # TODO Check with: if key not in cases_list:
            referential_computational_time.append( float( value ) )

        for i in range( len( cases_list ) ):
            t_live = computational_time[ i ]
            t_ref = referential_computational_time[ i ]
            t_dif_precentage = 100 * ( t_live - t_ref ) / t_ref

            if t_dif_precentage > 5:
                t_diff_percentage_string = f'<span style="color:red">__{t_dif_precentage:.1f} %__</span>'
            elif t_dif_precentage < -5:
                t_diff_percentage_string = f'<span style="color:green">__{t_dif_precentage:.1f} %__</span>'
            else:
                t_diff_percentage_string = f'{t_dif_precentage:.1f} %'

            computational_time_difference_formatted.append( t_diff_percentage_string )

def write_results():
    for i in range( len( cases_list ) ):
        print( f"Case: { cases_list[ i ] }\n{ results[ i ] }\nComputational time: { computational_time[ i ] }\n" )

    summary = { 'Cases' : cases_list,
                'Result' : results_fancy,
                'Comp. time' : computational_time,
                'Ref. comp. time' : referential_computational_time,
                'Comp. time dif.' : computational_time_difference_formatted }
    summary_df = pd.DataFrame( summary )
    with open(f'log_{strftime("%Y-%m-%d_%H:%M:%S")}.md', 'w') as f:
        f.write( f'Tests completed: {strftime("%Y-%m-%d %H:%M:%S", gmtime())}\n' )
        f.write( summary_df.to_markdown( ) )

if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser( description="Test all examples" )
    # TODO: Decetct the GPU automatically
    argparser.add_argument("--gpu", type=str, default="NVIDIA-A40",
            help="gpu model which runs the test")

    # parse the command line arguments
    args = argparser.parse_args()
    gpu_type = args.gpu

    run_cases()
    process_results( gpu_type )
    write_results()
