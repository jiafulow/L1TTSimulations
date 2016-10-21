#include <iostream>

#include "TFile.h"
#include "TTree.h"

void printBranches() {

    TString filename = "root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/SLHC/GEN/620_SLHC25p3_results/tt25_test/tracks_TTbar_PU140_sf1_nz8_L5x2_L10x2.root";  // from Olmo
    //TString filename = "root://xrootd2.ihepa.ufl.edu//store/user/jiafulow/SLHC/GEN/620_SLHC25p3_results/tt25_test/tracks_TTbar_PU140_sf1_nz8_L5x2_L10x2.root";  // from Olmo

    // Open the file
    //std::cout << "Opening file..." << std::endl;

    TFile* tfile = TFile::Open(filename);
    TTree* ttree = (TTree*) tfile->Get("ntupler/tree");

    TBranch* branch = 0;
    //TObject* obj = 0;

    TIter nextb(ttree->GetListOfBranches());
    while ((branch = (TBranch*)nextb())) {
        if (branch)
            std::cout << branch->GetName() << std::endl;
    }

    tfile->Close();
}
