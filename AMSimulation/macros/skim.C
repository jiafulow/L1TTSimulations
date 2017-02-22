#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void skim() {

  TString oldfilename = "root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/SLHC/GEN/620_SLHC28p1_results/jftest1/stubs_oc_tt25_200M_1.root";
  TString newfilename = "stubs_oc_tt25_100K.root";
  Long64_t nentries = 100;

  TFile * oldfile = TFile::Open(oldfilename);
  TTree * oldtree = (TTree *) oldfile->Get("ntupler/tree");

  TFile * newfile = new TFile(newfilename, "recreate");
  newfile->mkdir("ntupler")->cd();
  TTree * newtree = oldtree->CloneTree(0);

  for (Long64_t i = 0; i < nentries; ++i) {
    oldtree->GetEntry(i);
    newtree->Fill();
  }

  newtree->Print();
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
