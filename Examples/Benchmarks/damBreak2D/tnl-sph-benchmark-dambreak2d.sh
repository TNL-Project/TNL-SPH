#!/bin/bash

## benchmark setup
benchmarkName="local-test"
samples=1
resolutions="0.002"
#resolutions="0.005 0.002 0.001 0.00025"

## results folder
resultsFolder="results_"$benchmarkName
mkdir -p $resultsFolder

for resolution in ${resolutions}
do
   for sample in $(seq 1 ${samples})
   do
      #setup and run dualSPHysics code
      cd dualSPHysics_resources
      cp damBreak2D_WCSPH-DBC_Def_template.xml damBreak2D_WCSPH-DBC_Def.xml
      sed -i "s/resolutionPlaceholder/${resolution}/" damBreak2D_WCSPH-DBC_Def.xml
      ./run.sh
      cd ..
      #save results
      mv dualSPHysics_resources/damBreak2D_WCSPH-DBC_out/Run.out ${resultsFolder}/dualSPHysics_${resolution}_${sample}.out

      #setup and run TNL code
      cd TNL-SPH_resources
      python3 generateCase.py -resolution=${resolution}
      make clean
      make
      ./damBreak2D_WCSPH-DBC
      cd ..
      #save results
      mv TNL-SPH_resources/results/time_measurements.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.json
      mv TNL-SPH_resources/results/device.metadata.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.device_metadata.json
      mv TNL-SPH_resources/results/case_metadata.json ${resultsFolder}/tnl-sph_${resolution}_${sample}.case_metadata.json
   done
done
