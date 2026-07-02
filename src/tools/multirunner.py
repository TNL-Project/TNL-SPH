#! /usr/bin/env python3

from os import rename
from pathlib import Path
import subprocess
import json
import sys
import importlib.util
import pandas as pd
from time import strftime
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

# storage: one list of result dicts, each entry represents a single case run
# keys: tag, case, returncode, comp_time, tests_passed, tests_total, error
results_list = []

def init( case_dir, conf, verbose=False ):
    args = []
    args += [ case_dir / "init.py" ]
    for key, value in conf.items():
        if key not in [ "case", "case-tag", "evaluation-function" ]:
            args += [ f"--{key}", str( value ) ]

    if verbose:
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
    return p.returncode

def evaluate_test_metrics( case_dir, conf ):
    evaluation_function = conf[ "evaluation-function" ]
    if evaluation_function is None:
        return 0, 0
    return evaluation_function( case_dir )

def parse_tnl_sph_output( case_dir ):
    filename = case_dir / "results" / "time_measurements.json"
    try:
        with open( filename ) as f:
            timers_dictionary = json.load( f )
            return float( timers_dictionary[ "total" ] )
    except Exception as e:
        print( f"parse_tnl_sph_output: Could not read {filename}: {e}" )
        return 0

def run_cases( conf_list, store_results, verbose=False ):
    for conf in conf_list:
        case = conf[ "case" ]
        tag = conf[ "case-tag" ]
        case_dir = examples_dir / case

        entry = {
            "tag": tag,
            "case": case,
            "returncode": -1,
            "comp_time": 0.0,
            "tests_passed": 0,
            "tests_total": 0,
            "error": "",
        }

        try:
            print( f"Initializing case: {case} in {case_dir}." )
            init( case_dir, conf, verbose )
            print( f"Initialization finished." )
            bin_dir = build_dir / case
            print( f"Compiling case: {case} in {bin_dir}" )
            make( bin_dir )
            print( f"Executing case: {case} in {case_dir}." )
            entry[ "returncode" ] = run( case_dir )
            print( f"Execution finished with return code: {entry[ 'returncode' ]}." )
            entry[ "tests_passed" ], entry[ "tests_total" ] = evaluate_test_metrics( case_dir, conf )
            entry[ "comp_time" ] = parse_tnl_sph_output( case_dir )

            # backup the results only if the run succeeded and store_results is set
            if store_results and entry[ "returncode" ] == 0:
                results_dir = case_dir / "results"
                results_with_tag = "results_" + tag
                results_dir_renamed = case_dir / results_with_tag
                rename( results_dir, results_dir_renamed )
        except Exception as e:
            entry[ "error" ] = str( e )
            print( f"Case {case} FAILED: {e}" )

        results_list.append( entry )

def process_results():
    # nothing to pre-compute: write_results and print_summary_table read results_list directly
    pass

def write_results():
    for entry in results_list:
        print( f"Case: {entry[ 'case' ]}\n{entry[ 'returncode' ]}\nComputational time: {entry[ 'comp_time' ]}\n" )

    summary = {
        'Cases':             [ entry[ "tag" ] for entry in results_list ],
        'Result':            [ format_result_html( entry ) for entry in results_list ],
        'Comp. time':        [ entry[ "comp_time" ] for entry in results_list ],
        'Tests':             [ format_tests_html( entry ) for entry in results_list ],
        'Tests results':     [ format_tests_summary_html( entry ) for entry in results_list ],
    }
    summary_df = pd.DataFrame( summary )
    log_path = tools_dir / f'log_{strftime("%Y-%m-%d_%H:%M:%S")}.md'
    with open( log_path, 'w' ) as f:
        f.write( f'Tests completed: {strftime("%Y-%m-%d %H:%M:%S")}\n' )
        f.write( summary_df.to_markdown( ) or "" )

def format_result_html( entry ):
    if entry[ "returncode" ] == 0:
        return '<span style="color:green">__Success__</span>'
    return '<span style="color:red">__Failed__</span>'

def format_tests_html( entry ):
    passed = entry[ "tests_passed" ]
    total = entry[ "tests_total" ]
    if passed == total:
        if total > 0:
            return f'<span style="color:green">__{passed}/{total}__</span>'
        return '-'
    return f'<span style="color:red">__{passed}/{total}__</span>'

def format_tests_summary_html( entry ):
    passed = entry[ "tests_passed" ]
    total = entry[ "tests_total" ]
    if passed == total:
        if total > 0:
            return f'<span style="color:green">__Passed__</span>'
        return '-'
    return f'<span style="color:red">__Failed__</span>'

def print_summary_table():
    console = Console()
    table = Table( title="Multirunner Results", show_lines=True )
    table.add_column( "Case tag" )
    table.add_column( "Result" )
    table.add_column( "Comp. time" )
    table.add_column( "Tests" )
    table.add_column( "Tests result" )
    table.add_column( "Error" )

    for entry in results_list:
        if entry[ "returncode" ] == 0:
            result_str = "[green]Success[/green]"
        else:
            result_str = "[red]Failed[/red]"
        passed = entry[ "tests_passed" ]
        total = entry[ "tests_total" ]
        if total == 0:
            tests_str = "-"
            tests_result_str = "-"
        elif passed == total:
            tests_str = f"[green]{passed}/{total}[/green]"
            tests_result_str = "[green]Passed[/green]"
        else:
            tests_str = f"[red]{passed}/{total}[/red]"
            tests_result_str = "[red]Failed[/red]"
        err = entry[ "error" ]
        err_str = err[ :60 ] if err else ""
        table.add_row( entry[ "tag" ], result_str, str( entry[ "comp_time" ] ), tests_str, tests_result_str, err_str )

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
    run_cases( conf_list, args.store_results, args.verbose )
    process_results()
    write_results()
    if args.verbose:
        print_summary_table()
