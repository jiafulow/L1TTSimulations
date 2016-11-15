#include "RoadMerging.cc"


void runRoadMerging() {
  TString patternBankFileName = "../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.0.root";
  TString roadFileName        = "../dataFiles/roads_oc_tt25_sf1_nz8_pt3_test.root";

  RoadMerging rm;
  rm.process(roadFileName, patternBankFileName);
}

#ifndef __CINT__
int main()
{
  runRoadMerging();
  return 0;
}
#endif
