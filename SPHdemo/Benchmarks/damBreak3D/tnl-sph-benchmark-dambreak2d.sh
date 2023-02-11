#!/bin/bash

benchmarkName="local-test"
resolutions="0.02"
#resolutions="0.005 0.002 0.001 0.00025"

## results folder
resultsFolder="results_"$benchmarkName
mkdir $resultsFolder

## setup and run dualSPHysics code
for resolution in $resolutions
do
   #for sample in {1..5}
   for sample in 1
   do
      cd dualSPHysics_resources
      cp damBreak3D_WCSPH-DBC_Def_template.xml damBreak3D_WCSPH-DBC_Def.xml
      sed -i "s/resolutionPlaceholder/${resolution}/" damBreak3D_WCSPH-DBC_Def.xml
      ./run.sh
      cd ..

      ## setup and run TNL code
      cd TNL-SPH_resources
      ./generateCase.sh $resolution
      #python3 generateCase.py -resolution=$resolution
      make clean
      make
      ./damBreak3D_WCSPH-DBC
      cd ..

      ## save results
      # cp dualSPHysics_resources/damBreak2D_WCSPH-DBC_out/Run.out results_local/dualSPHysics_${resolution}_${sample}.out
      # cp TNL-SPH_resources/time_measurements.json results_local/tnl-sph_${resolution}_${sample}.json
      cp dualSPHysics_resources/damBreak3D_WCSPH-DBC_out/Run.out $resultsFolder/dualSPHysics_${resolution}_${sample}.out
      cp TNL-SPH_resources/time_measurements.json $resultsFolder/tnl-sph_${resolution}_${sample}.json
      ## save other files
      cp TNL-SPH_resources/device.metadata.json $resultsFolder/tnl-sph_${resolution}_${sample}.device_metadata.json
      cp TNL-SPH_resources/case_metadata.json $resultsFolder/tnl-sph_${resolution}_${sample}.case_metadata.json
   done
done

## process
