#include "PatternMerging.cc"


void runPatternMerging(TString bank, TString out, unsigned deltaN=80000, float targetCoverage=1.00) {
  //TString bank         = "../dataFiles/patternBank_tt27_sf1L0x2_nz8_pt3_300M.root";
  //TString out          = "../dataFiles/m8_patternBank_tt27_sf1L0x2_nz8_pt3_300M.root";
  ////unsigned deltaN      = 36000;
  //unsigned deltaN      = 60000;
  //float targetCoverage = 1.00;
  ////float targetCoverage = 1.95;

  //TString bank         = "../dataFiles/patternBank_oc_tt25_sf1_nz8_pt3_400M.root";
  //TString out          = "../dataFiles/m8_patternBank_oc_tt25_sf1_nz8_pt3_400M.root";
  //unsigned deltaN      = 60000;
  //float targetCoverage = 1.00;

  PatternMerging pm;
  pm.mergePatterns(bank, out, deltaN, targetCoverage);
}

#ifndef __CINT__
int main(int argc, char** argv)
{
  TString bank = "", out = "";

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " BANK OUT" << std::endl;
  } else {
    bank = argv[1];
    out  = argv[2];
    std::cerr << "I got BANK=" << bank << " OUT=" << out << std::endl;
  }

  runPatternMerging(bank, out);
  return 0;
}
#endif
