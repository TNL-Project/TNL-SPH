#! /usr/bin/env python3

from os import rename
from pathlib import Path
import subprocess
import json
import sys
import importlib.util
import pandas as pd
from time import strftime, gmtime
from rich.console import Console
from rich.table import Table


def load_configurations( conf_path ):
    spec = importlib.util.spec_from_file_location( "configurations_module", conf_path )
    if spec is None or spec.loader is None:
        raise SystemExit( f"Could not load config module from {conf_path}." )
    module = importlib.util.module_from_spec( spec )
    sys.path.insert( 0, str( Path( conf_path ).resolve().parent ) )
    spec.loader.exec_module( module )
    if not hasattr( module, "configurations" ) or not isinstance( module.configurations, list ):
        raise SystemExit( f"Config file {conf_path} must define a top-level 'configurations' list." )
    return module.configurations


# initialize directories
tools_dir = Path(__file__).parent
project_dir = ( tools_dir / ".." / ".." ).resolve()
examples_dir = project_dir / "examples"
build_dir = project_dir / "build" / "examples"

# storage arrays
results = []
results_fancy = []
results_returncode = []
computational_time = []
referential_computational_time = []
computational_time_difference_formatted = []
cases_tags_list = []
cases_list = []
tests_passed_list = []
tests_total_list = []
tests_output_formatted = []
tests_summary_formatted = []
tests_logs_list = []
errors_list = []

def init( case_dir, conf ):
    args = []
    args += [ case_dir / "init.py" ]
    for key, value in conf.items():
        if key not in [ "case", "case-tag", "evaluation-function" ]:
            args += [ f"--{key}", str( value ) ]

    print( args )
    subprocess.run( args,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                    cwd=case_dir )

def make( bin_dir ):
    args = [ "make" ]

    p = subprocess.run( args,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                        cwd=bin_dir,
                        text=True )
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

def run( case_dir ):
    args = []
    args += [ case_dir / "run.py" ]

    # run the process and print its output as it is being executed
    #with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    #                      bufsize=1, cwd=case_dir, text=True) as p:
    #    for line in p.stdout:
    #        print(line, end="")

    p = subprocess.run( args,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                        cwd=case_dir,
                        text=True )
    results_returncode.append( subprocess.CalledProcessError( p.returncode, p.args ) )
    results.append( p.returncode )

def evaluate_test_metrics( case_dir, conf ):
    evaluation_function = conf[ "evaluation-function" ]
    if evaluation_function is None:
        tests_passed_list.append( 0 )
        tests_total_list.append( 0 )
        return
    tests_passed, tests_total = evaluation_function( case_dir )
    tests_passed_list.append( tests_passed )
    tests_total_list.append( tests_total )

def parse_tnl_sph_output( case_dir ):
    filename = case_dir / "results" / "time_measurements.json"
    try:
        with open( filename ) as f:
            lines = json.load( f )
            json_str = json.dumps( lines )
            timers_dictionary = json.loads( json_str )
            return float( timers_dictionary[ "total" ] )
    except Exception as e:
        print( f"parse_tnl_sph_output: File {filename} not found." )
        return 0

def run_cases( conf_list, store_results ):
    for conf in conf_list:
        case = conf[ "case" ]
        cases_list.append( case )
        cases_tags_list.append( conf[ "case-tag" ] )
        case_dir = examples_dir / case

        results_len_before = len( results )
        tests_passed_len_before = len( tests_passed_list )
        tests_total_len_before = len( tests_total_list )
        comp_time_len_before = len( computational_time )

        try:
            print( f"Initializing case: {case} in {case_dir}." )
            init( case_dir, conf )
            print( f"Initialization finished." )
            bin_dir = build_dir / case
            print( f"Compiling case: {case} in {bin_dir}" )
            make( bin_dir )
            print( f"Executing case: {case} in {case_dir}." )
            run( case_dir )
            print( f"Execution finished with return code: {results[ -1 ]}." )
            evaluate_test_metrics( case_dir, conf )

            # get computational time
            computational_time.append( parse_tnl_sph_output( case_dir ) )

            # backup the results only if the run succeeded and store_results is set
            if store_results and results[ -1 ] == 0:
                results_dir = case_dir / "results"
                results_with_tag = "results_" + conf[ "case-tag" ]
                results_dir_renamed = case_dir / results_with_tag
                rename( results_dir, results_dir_renamed )

            errors_list.append( "" )
        except Exception as e:
            # ensure all parallel lists keep an entry for this case
            if len( results ) > results_len_before:
                results[ -1 ] = -1
            else:
                results.append( -1 )
            if len( tests_passed_list ) == tests_passed_len_before:
                tests_passed_list.append( 0 )
            if len( tests_total_list ) == tests_total_len_before:
                tests_total_list.append( 0 )
            if len( computational_time ) == comp_time_len_before:
                computational_time.append( 0.0 )
            errors_list.append( str( e ) )
            print( f"Case {case} FAILED: {e}" )
            continue

def process_results():
    # parse return codes to fancy output
    for entry in results:
        if entry == 0:
            results_fancy.append( '<span style="color:green">__Success__</span>' )
        else:
            results_fancy.append( '<span style="color:red">__Failed__</span>' )

    for i in range( len( cases_list ) ):
        # process test results
        tests_passed = tests_passed_list[ i ]
        tests_total = tests_total_list[ i ]
        if tests_passed == tests_total:
            if tests_total > 0:
                tests_output_string = f'<span style="color:green">__{tests_passed}/{tests_total}__</span>'
                tests_summary_string = f'<span style="color:green">__Passed__</span>'
            else:
                tests_output_string = '-'
                tests_summary_string = '-'
        else:
            tests_output_string = f'<span style="color:red">__{tests_passed}/{tests_total}__</span>'
            tests_summary_string = f'<span style="color:red">__Failed__</span>'

        tests_summary_formatted.append( tests_summary_string )
        tests_output_formatted.append( tests_output_string )

def write_results():
    for i in range( len( cases_list ) ):
        print( f"Case: { cases_list[ i ] }\n{ results[ i ] }\nComputational time: { computational_time[ i ] }\n" )

    summary = { 'Cases' : cases_tags_list,
                'Result' : results_fancy,
                'Comp. time' : computational_time,
                'Tests' : tests_output_formatted,
                'Tests results' : tests_summary_formatted }
    summary_df = pd.DataFrame( summary )
    with open(f'log_{strftime("%Y-%m-%d_%H:%M:%S")}.md', 'w') as f:
        f.write( f'Tests completed: {strftime("%Y-%m-%d %H:%M:%S")}\n' )
        f.write( summary_df.to_markdown( ) or "" )

def print_summary_table():
    console = Console()
    table = Table( title="Multirunner Results", show_lines=True )
    table.add_column( "Case tag" )
    table.add_column( "Result" )
    table.add_column( "Comp. time" )
    table.add_column( "Tests" )
    table.add_column( "Tests result" )
    table.add_column( "Error" )

    for i in range( len( cases_tags_list ) ):
        if results[ i ] == 0:
            result_str = "[green]Success[/green]"
        else:
            result_str = "[red]Failed[/red]"
        passed = tests_passed_list[ i ]
        total = tests_total_list[ i ]
        if total == 0:
            tests_str = "-"
            tests_result_str = "-"
        elif passed == total:
            tests_str = f"[green]{passed}/{total}[/green]"
            tests_result_str = "[green]Passed[/green]"
        else:
            tests_str = f"[red]{passed}/{total}[/red]"
            tests_result_str = "[red]Failed[/red]"
        err = errors_list[ i ]
        if err:
            err_str = err[ :60 ]
        else:
            err_str = ""
        table.add_row( cases_tags_list[ i ], result_str, str( computational_time[ i ] ), tests_str, tests_result_str, err_str )

    console.print( table )

if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser( description="Run multiple SPH examples from a config file" )
    argparser.add_argument( "--conf", type=str, required=True, help="path to Python config file defining a 'configurations' list" )
    argparser.add_argument( "--store-results", action="store_true", default=False, help="rename results/ to results_<case-tag>/ after each successful run" )
    argparser.add_argument( "--verbose", action="store_true", default=False, help="print a colored summary table to the terminal" )

    # parse the command line arguments
    args = argparser.parse_args()
    conf_list = load_configurations( args.conf )
    run_cases( conf_list, args.store_results )
    process_results()
    write_results()
    if args.verbose:
        print_summary_table()
