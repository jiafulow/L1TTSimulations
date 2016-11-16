#include "RoadMerging.cc"


void runRoadMerging() {
  gROOT->ProcessLine(".L Loader.C+");  // fix missing dictionaries

  //TString patternBankFileName = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/620_SLHC25p3_results/jftest1/rossin_20161022/m8_patternBank_tt27_sf1L0x2L5x2_nz8_pt3_100c_250M.root";
  //TString roadFileName        = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/620_SLHC25p3_results/test1/tt27_sf1L0x2L5x2_nz8_pt3_20160308/roads_TTbar_PU140_tt27_sf1L0x2L5x2_nz8_pt3_5oo6_100c_250M.root";
  //TString outFileName         = "../dataFiles/roads_TTbar_PU140_tt27_sf1L0x2L5x2_nz8_pt3_5oo6_100c_250M_x8_FM_.root";

  TString patternBankFileName = "../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.0.root";  //FIXME
  TString roadFileName        = "../dataFiles/roads_oc_tt25_sf1_nz8_pt3.0.root";
  TString outFileName         = "../dataFiles/m8_roads_oc_tt25_sf1_nz8_pt3.root";

  RoadMerging rm;
  rm.process(roadFileName, outFileName, patternBankFileName);
}

#ifndef __CINT__
int main()
{
  runRoadMerging();
  return 0;
}
#endif
