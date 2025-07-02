import numpy as np
from pathlib import Path

# initialize directories
tools_dir = Path(__file__).parent
project_dir = ( tools_dir / ".." / ".." / ".." ).resolve()
examples_dir = project_dir / "examples"

def damBreak2D_WCSPH_DBC( case_dir ):
    example_name = "damBreak2D_WCSPH-DBC"
    resources_dir = ( case_dir / ".." / ".." / "resources" / "damBreak2D" / "damBreak2D_experimentalDataLobovsky2014" ).resolve()

    # Referential data
    wl_avg_ref = np.array( [ 0.23116514, 0.06056554, 0.03638269, 0.02139411 ] )
    wl_std_ref = np.array( [ 0.04803559, 0.03703002, 0.02929507, 0.02061927 ] )
    ps_avg_ref = np.array( [ 1237.25545429,  914.68332571, 682.51443457, 187.89957863 ] )
    ps_std_ref = np.array( [ 1920.12293911, 1379.96908324, 984.69733283, 521.99612149 ] )

    tests_successful_count = 0
    tests_total_count = 4

    # Compare pure simulation data
    #results_dir = examples_dir / "WCSPH-DBC" / example_name  / "results"
    results_dir = case_dir / "results"

    wl_file = results_dir / "sensorsWaterLevel.dat"
    wl = np.genfromtxt( wl_file, delimiter=' ' )
    wl = wl[ :, 1: ]
    wl_avg_sim = np.mean( wl, axis=0 )
    wl_std_sim = np.std( wl, axis=0 )
    wl_avg = np.allclose( wl_avg_sim, wl_avg_ref, atol=0.1 )
    wl_std = np.allclose( wl_std_sim, wl_std_ref, atol=0.01 )

    if wl_avg:
        tests_successful_count += 1
        print( f"{example_name}: Water level average comparison - successful." )
    else:
        print( f"{example_name}: Water level average comparison - failed." )

    if wl_std:
        tests_successful_count += 1
        print( f"{example_name}: Water level std comparison - successful." )
    else:
        print( f"{example_name}: Water level std comparison - failed." )

    ps_file = results_dir / "sensorsPressure.dat"
    ps = np.genfromtxt( ps_file, delimiter=' ' )
    ps = ps[ :, 1: ]
    ps_avg_sim = np.mean( ps, axis=0 )
    ps_std_sim = np.std( ps, axis=0 )
    ps_avg = np.allclose( ps_avg_sim, ps_avg_ref, atol=100 )
    ps_std = np.allclose( ps_std_sim, ps_std_ref, atol=100 )

    if ps_avg:
        tests_successful_count += 1
        print( f"{example_name}: Pressure sensors average comparison - successful." )
    else:
        print( f"{example_name}: Pressure sensors average comparison - failed." )

    if ps_std:
        tests_successful_count += 1
        print( f"{example_name}: Pressure sensors std comparison - successful." )
    else:
        print( f"{example_name}: Pressure sensors std comparison - failed." )

    return tests_successful_count, tests_total_count

def damBreak2D_WCSPH_MDBC( case_dir ):
    return 0, 0

def damBreak3D_WCSPH_DBC( case_dir ):
    return 0, 0

def damBreak3D_WCSPH_MDBC( case_dir ):
    return 0, 0

def poiseuilleFlowWithOpenBoundary2D_WCSPH( case_dir ):
    return 0, 0

def poiseuilleFlowWithPeriodicBoundary2D_WCSPH( case_dir ):
    return 0, 0

def damBreak2D_WCSPH_BI( case_dir ):
    return 0, 0

def damBreak3D_WCSPH_BI( case_dir ):
    return 0, 0

def poiseuilleFlowWithOpenBoundary2D_WCSPH_BI( case_dir ):
    return 0, 0

def poiseuilleFlowWithPeriodicBoundary2D_WCSPH_BI( case_dir ):
    return 0, 0

def empty( case_dir ):
    return 0, 0

if __name__ == "__main__":

    print( "Evaluate metrics - main." )
    damBreak2D_WCSPH_DBC()
