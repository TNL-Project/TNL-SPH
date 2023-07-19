#!/bin/bash

## benchmark setup
benchmarkName="local-test"
samples=1
resolutions="0.02"
#resolutions="0.02 0.01 0.005"
#dualSPHysicsPath=/home/tomas/mount/home/tomas/Documents/DualSPHysics_v5.0.5/DualSPHysics_v5.0/bin/linux

## results folder
resultsFolder="results_"$benchmarkName
mkdir -p $resultsFolder

for resolution in ${resolutions}
do
   for sample in $(seq 1 ${samples})
   do
      #setup and run dualSPHysics code if path to DualSPHysics is defined
      if [ -n "${dualSPHysicsPath+1}" ]
      then
         cd dualSPHysics_resources
         cp damBreak3D_WCSPH-DBC_Def_template.xml damBreak3D_WCSPH-DBC_Def.xml
         sed -i "s/resolutionPlaceholder/${resolution}/" damBreak3D_WCSPH-DBC_Def.xml
         ./run.sh ${dualSPHysicsPath}
         cd ..
         mv dualSPHysics_resources/damBreak3D_WCSPH-DBC_out/Run.out ${resultsFolder}/dualSPHysics_${resolution}_${sample}.out
      fi

      #setup and run TNL code
      cd TNL-SPH_resources
      python3 ./generateCaseWithDSPHGenCase.py -resolution=${resolution} --generateNewGeometry
      make clean
      make
      ./damBreak3D_WCSPH-DBC
      cd ..
      mv TNL-SPH_resources/results/time_measurements.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.json
      mv TNL-SPH_resources/results/device.metadata.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.device_metadata.json
      mv TNL-SPH_resources/results/case_metadata.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.case_metadata.json
   done
done
