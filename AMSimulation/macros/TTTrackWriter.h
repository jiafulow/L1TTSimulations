#ifndef _TTTRACKWRITER_H_
#define _TTTRACKWRITER_H_

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

// _____________________________________________________________________________
// This is a simple wrapper around TTree

class TTTrackWriter {
public:
    TTTrackWriter(int verbose=1);
    ~TTTrackWriter();

    void init(TTree* intree, TString out, TString prefix="", TString suffix="");

    void fill() { ttree_->Fill(); }

    Long64_t writeTree();

    TTree * getTree() { return ttree_; }

protected:
    TFile* tfile_;
    TTree* ttree_;
    const int verbose_;
};
#endif  // _TTTRACKWRITER_H_

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation
#define _TTTRACKWRITER_CXX_
#ifdef _TTTRACKWRITER_CXX_

#include <cassert>
#include <stdexcept>

TTTrackWriter::TTTrackWriter(int verbose) :
    verbose_(verbose) {}

TTTrackWriter::~TTTrackWriter() {
    if (ttree_)  delete ttree_;
    if (tfile_)  delete tfile_;
}

void TTTrackWriter::init(TTree* intree, TString out, TString prefix, TString suffix) {
    if (!out.EndsWith(".root")) {
        TString msg = "Output filename must be .root";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << out << std::endl;
    tfile_ = TFile::Open(out, "RECREATE");

    if (tfile_) {
        //if (verbose_)  std::cout << "Successfully opened " << out << std::endl;
    } else {
        TString msg = "Failed to open " + out;
        throw std::invalid_argument(msg.Data());
    }

    tfile_->mkdir("ntupler")->cd();
    ttree_ = (TTree*) intree->CloneTree(0); // Do not copy the data yet
}

Long64_t TTTrackWriter::writeTree() {
    Long64_t nentries = ttree_->GetEntries();
    tfile_->Write();
    //tfile_->Close();
    return nentries;
}
#endif  // _TTTRACKWRITER_CXX_
