#include "PatternMerging.cc"


void runPatternMerging() {
  TString patternBankFileName     = "../dataFiles/patternBank_tt27_sf1L0x2_nz8_pt3_300M.root";
  TString tempPatternBankFileName = "../dataFiles/m8_patternBank_tt27_sf1L0x2_nz8_pt3_300M.root";
  unsigned deltaN                 = 36000;
  float targetCoverage            = 1.00;

  PatternMerging pm;
  pm.mergePatterns(patternBankFileName, tempPatternBankFileName, deltaN, targetCoverage);
}

#ifndef __CINT__
int main()
{
  runPatternMerging();
  return 0;
}
#endif
