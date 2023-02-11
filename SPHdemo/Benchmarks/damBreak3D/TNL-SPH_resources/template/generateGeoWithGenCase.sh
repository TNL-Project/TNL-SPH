#!/bin/bash
fail () {
    echo Execution aborted.
    exit 1
}

#jmeno uhlohy, vystupni slozky
export name=damBreak3D_WCSPH-DBC
export dirout=${name}_out
export diroutdata=${dirout}/data

#nacteni prislusnych souboru
export dirbin=/home/tomas/mount/home/tomas/Documents/DualSPHysics_v5.0.5/DualSPHysics_v5.0/bin/linux
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${dirbin}
export gencase="${dirbin}/GenCase_linux64"
export dualsphysicscpu="${dirbin}/DualSPHysics5.0CPU_linux64"
export dualsphysicsgpu="${dirbin}/DualSPHysics5.0_linux64"
export boundaryvtk="${dirbin}/BoundaryVTK_linux64"
export partvtk="${dirbin}/PartVTK_linux64"
export partvtkout="${dirbin}/PartVTKOut_linux64"
export measuretool="${dirbin}/MeasureTool_linux64"
export computeforces="${dirbin}/ComputeForces_linux64"
export isosurface="${dirbin}/IsoSurface_linux64"
export flowtool="${dirbin}/FlowTool_linux64"
export floatinginfo="${dirbin}/FloatingInfo_linux64"

# "dirout" to store results is removed if it already exists
if [ -e ${dirout} ]; then rm -r ${dirout}; fi

# CODES are executed according the selected parameters of execution in this testcase
# Executes GenCase4 to create initial files for simulation.
${gencase} ${name}_Def ${dirout}/${name} -save:vtkfluid,vtkbound
if [ $? -ne 0 ] ; then fail; fi

echo Generated.
