#include "PatternBankReader.h"

#include <iostream>

//TString filename = "/uscms_data/d3/ocerri/CMSTrigger/CMSSW_6_2_0_SLHC25_patch3/src/PatternBanks/patternBank_tt25_sf1_nz8_L5x2_pt3_200M.root";  // from Olmo
TString filename = "root://xrootd2.ihepa.ufl.edu//store/user/jiafulow/L1TrackTrigger/6_2_0_SLHC25p3_results/tt25_test/patternBank_tt25_sf1_nz8_L5x2_pt3_200M.root";  // from Olmo


void analyzePatternBank() {
    // Open the file
    std::cout << "Opening file..." << std::endl;

    PatternBankReader reader;
    reader.init(filename);

    // Get number of events
    //Long64_t nentries = reader.getEntries();
    Long64_t nentries = 10;
    std::cout << "Number of patterns: " << nentries << std::endl;

    // Get pattern bank info
    reader.getPatternInfo();

    float    coverage    = reader.pb_coverage;
    unsigned count       = reader.pb_count;
    //unsigned magicNumber = reader.pb_superstrip_nx;
    unsigned magicNumber = 0;

    std::cout << "Bank coverage: " << coverage
              << " count: "        << count
              << " magic number: " << magicNumber
              << std::endl;

    // Loop over patterns
    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        // Retrieve the pattern
        reader.getPattern(ievt);
        reader.getPatternAttr(ievt);

        // Get the variables
        // See PatternBankReader.h for more info
        std::vector<unsigned> superstripIds = *(reader.pb_superstripIds);
        unsigned              frequency     = reader.pb_frequency;
        float                 patternInvPt  = reader.pb_invPt_mean;
        float                 patternPhi    = reader.pb_phi_mean;

        std::cout << ".. ipatt: "   << ievt
                  << " frequency: " << frequency
                  << " invPt: "     << patternInvPt
                  << " phi: "       << patternPhi
                  << std::endl;
    }  // end loop over patterns

}  // end analyzePatternBank()

#ifndef __CINT__
int main()
{
  analyzePatternBank();
  return 0;
}
#endif
