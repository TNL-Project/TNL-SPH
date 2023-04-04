#!/bin/bash

#Path to data from simulation
resunltsSensorsPressure="results/particles_sensors.dat"
resunltsSensorsWaterLevel="results/particles_sensorsWaterLevel.dat"
resultsPlotPath="results/postprocessing"

mkdir ${resultsPlotPath}

caseName="damBreak2D_WCSPH-BI"

#Path to folder with experimental data and postprocessing scripts.
resourcesPath="../resources/damBreak2D"

plotPressure_s1="$resourcesPath/plotPressure_s1.gp"
plotPressure_s2="$resourcesPath/plotPressure_s2.gp"
plotPressure_s3="$resourcesPath/plotPressure_s3.gp"
plotPressure_s4="$resourcesPath/plotPressure_s4.gp"

plotWaterLevel_s1="$resourcesPath/plotWaterLevel_h1.gp"
plotWaterLevel_s2="$resourcesPath/plotWaterLevel_h2.gp"
plotWaterLevel_s3="$resourcesPath/plotWaterLevel_h3.gp"
plotWaterLevel_s4="$resourcesPath/plotWaterLevel_h4.gp"

pressureExperimentalData="$resourcesPath/damBreak2D_experimentalDataLobovsky2014/pressure/"
waterLevelExperimentalData="$resourcesPath/damBreak2D_experimentalDataLobovsky2014/waterLevel/"

#Plot pressure
gnuplot -e "num='$resunltsSensorsPressure'" \
        -e "exp='$pressureExperimentalData/Fig18_peak_event_5_sensors_5.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_pressure_sensor1.png'" \
        $plotPressure_s1

gnuplot -e "num='$resunltsSensorsPressure'" \
        -e "exp='$pressureExperimentalData/Fig18_peak_event_5_sensors_4.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_pressure_sensor2.png'" \
        $plotPressure_s2

gnuplot -e "num='$resunltsSensorsPressure'" \
        -e "exp='$pressureExperimentalData/Fig18_peak_event_5_sensors_2.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_pressure_sensor3.png'" \
        $plotPressure_s3

gnuplot -e "num='$resunltsSensorsPressure'" \
        -e "exp='$pressureExperimentalData/Fig18_peak_event_5_sensors_1.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_pressure_sensor4.png'" \
        $plotPressure_s4

#Plot water level
gnuplot -e "num='$resunltsSensorsWaterLevel'" \
        -e "exp='$waterLevelExperimentalData/Fig16_WaterLevels_H1_2.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_waterLevel_height1.png'" \
        $plotWaterLevel_s1

gnuplot -e "num='$resunltsSensorsWaterLevel'" \
        -e "exp='$waterLevelExperimentalData/Fig16_WaterLevels_H2_3.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_waterLevel_height2.png'" \
        $plotWaterLevel_s2

gnuplot -e "num='$resunltsSensorsWaterLevel'" \
        -e "exp='$waterLevelExperimentalData/Fig16_WaterLevels_H3_3.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_waterLevel_height3.png'" \
        $plotWaterLevel_s3

gnuplot -e "num='$resunltsSensorsWaterLevel'" \
        -e "exp='$waterLevelExperimentalData/Fig16_WaterLevels_H4_4.dat'" \
        -e "outputFileName='$resultsPlotPath/results_${caseName}_waterLevel_height4.png'" \
        $plotWaterLevel_s4
