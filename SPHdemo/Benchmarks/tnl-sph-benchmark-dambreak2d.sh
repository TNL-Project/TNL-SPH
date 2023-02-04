#!/bin/bash

resolutions="0.002"
#resolutions="0.005 0.002 0.001 0.00025"

## setup and run dualSPHysics code
for resolution in $resolutions
do
   #for sample in {1..5}
   for sample in 1
   do
      cd dualSPHysics_resources
      cp damBreak2D_WCSPH-DBC_Def_template.xml damBreak2D_WCSPH-DBC_Def.xml
      sed -i "s/resolutionPlaceholder/${resolution}/" damBreak2D_WCSPH-DBC_Def.xml
      ./run.sh
      cd ..

      ## setup and run TNL code
      cd TNL-SPH_resources
      python generateCase.py -resolution=$resolution
      make
      ./damBreak2D_WCSPH-DBC
      cd ..

      ## save results
      cp dualSPHysics_resources/damBreak2D_WCSPH-DBC_out/Run.out dualSPHysics_${resolution}_${sample}.out
      cp TNL-SPH_resources/time_measurements.json tnl-sph_${resolution}_${sample}.json
   done
done

## process
