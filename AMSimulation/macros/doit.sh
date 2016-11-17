#!/bin/bash

./runPatternMerging \
    ../dataFiles/patternBank_oc_tt25_sf1_nz8_pt3_400M.root \
    ../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.root >& \
    ../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.log &

./runRoadMerging \
    ../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.1.root \
    ../dataFiles/roads_TTbar_PU140.root \
    ../dataFiles/m8_roads_TTbar_PU140.root >& \
    ../dataFiles/m8_roads_TTbar_PU140.log &

./runRoadMerging \
    ../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.1.root \
    ../dataFiles/roads_TTbar_PU200.root \
    ../dataFiles/m8_roads_TTbar_PU200.root >& \
    ../dataFiles/m8_roads_TTbar_PU200.log &

