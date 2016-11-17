#include "RoadMerging.cc"


void runRoadMerging(TString bank, TString road, TString out)
{
  gROOT->ProcessLine(".L Loader.C+");  // fix missing dictionaries

  //TString bank = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/620_SLHC25p3_results/jftest1/rossin_20161022/m8_patternBank_tt27_sf1L0x2L5x2_nz8_pt3_100c_250M.root";
  //TString road = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/620_SLHC25p3_results/test1/tt27_sf1L0x2L5x2_nz8_pt3_20160308/roads_TTbar_PU140_tt27_sf1L0x2L5x2_nz8_pt3_5oo6_100c_250M.root";
  //TString out  = "../dataFiles/roads_TTbar_PU140_tt27_sf1L0x2L5x2_nz8_pt3_5oo6_100c_250M_x8_FM_.root";

  //TString bank = "../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.0.root";
  //TString road = "../dataFiles/roads_oc_tt25_sf1_nz8_pt3.0.root";
  //TString out  = "../dataFiles/m8_roads_oc_tt25_sf1_nz8_pt3.root";

  RoadMerging rm;
  rm.process(bank, road, out);
}

#ifndef __CINT__
int main(int argc, char** argv)
{
  TString bank = "", road = "", out = "";

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " BANK ROAD OUT" << std::endl;
  } else {
    bank = argv[1];
    road = argv[2];
    out  = argv[3];
    std::cerr << "I got BANK=" << bank << " ROAD=" << road << " OUT=" << out << std::endl;
  }

  runRoadMerging(bank, road, out);
  return 0;
}
#endif
