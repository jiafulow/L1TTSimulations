#ifndef _PATTERNBANKREADER_H_
#define _PATTERNBANKREADER_H_

#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

// _____________________________________________________________________________
// This is a simple wrapper around TTree. It sets the branch names, addresses
// etc. Its functions are essentially the same as the functions of TTree.

class PatternBankReader {
public:
    PatternBankReader(int verbose=1);
    ~PatternBankReader();

    void init(TString src);

    Int_t getPattern    (Long64_t entry)   { return ttree1_->GetEntry(entry); }
    Int_t getPatternInfo(Long64_t entry=0) { return ttree2_->GetEntry(0);     }  // only 1 entry
    Int_t getPatternAttr(Long64_t entry)   { return ttree3_->GetEntry(entry); }

    Long64_t getEntries() const { return ttree1_->GetEntries(); }

    TTree * getTree()     { return ttree1_; }
    TTree * getInfoTree() { return ttree2_; }
    TTree * getAttrTree() { return ttree3_; }

    // Typedefs
    typedef uint32_t superstrip_type;
    typedef uint16_t frequency_type;

    // Pattern attributes
    float                          pb_invPt_mean;
    float                          pb_invPt_sigma;
    float                          pb_cotTheta_mean;
    float                          pb_cotTheta_sigma;
    float                          pb_phi_mean;
    float                          pb_phi_sigma;
    float                          pb_z0_mean;
    float                          pb_z0_sigma;

    // Pattern bank statistics
    float                          pb_coverage;
    unsigned                       pb_count;
    unsigned                       pb_tower;
    std::string *                  pb_superstrip;
    //unsigned                       pb_superstrip_nx;
    //unsigned                       pb_superstrip_nz;

    // Pattern bank
    frequency_type                 pb_frequency;
    std::vector<superstrip_type> * pb_superstripIds;

protected:
    TFile* tfile_;
    TTree* ttree1_;  // for pattern bank
    TTree* ttree2_;  // for pattern bank statistics
    TTree* ttree3_;  // for pattern attributes
    const int verbose_;
};
#endif  // _PATTERNBANKREADER_H_

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation
#define _PATTERNBANKREADER_CXX_
#ifdef _PATTERNBANKREADER_CXX_

#include <cassert>
#include <stdexcept>

PatternBankReader::PatternBankReader(int verbose) :
    pb_invPt_mean     (0.),
    pb_invPt_sigma    (0.),
    pb_cotTheta_mean  (0.),
    pb_cotTheta_sigma (0.),
    pb_phi_mean       (0.),
    pb_phi_sigma      (0.),
    pb_z0_mean        (0.),
    pb_z0_sigma       (0.),
    //
    pb_coverage       (0.),
    pb_count          (0),
    pb_tower          (0),
    pb_superstrip     (0),
    //pb_superstrip_nx  (0),
    //pb_superstrip_nz  (0),
    //
    pb_frequency      (0),
    pb_superstripIds  (0),
    //
    verbose_(verbose) {}

PatternBankReader::~PatternBankReader() {
    if (ttree3_) delete ttree3_;
    if (ttree2_) delete ttree2_;
    if (ttree1_) delete ttree1_;
    if (tfile_)  delete tfile_;
}

void PatternBankReader::init(TString src) {
    if (!src.EndsWith(".root")) {
        TString msg = "Input source must be .root";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << src << std::endl;
    tfile_ = TFile::Open(src);

    if (tfile_) {
        //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
    } else {
        TString msg = "Failed to read " + src;
        throw std::invalid_argument(msg.Data());
    }

    ttree3_ = (TTree*) tfile_->Get("patternAttributes");
    assert(ttree3_ != 0);

    ttree3_->SetBranchAddress("invPt_mean"    , &(pb_invPt_mean));
    ttree3_->SetBranchAddress("invPt_sigma"   , &(pb_invPt_sigma));
    ttree3_->SetBranchAddress("cotTheta_mean" , &(pb_cotTheta_mean));
    ttree3_->SetBranchAddress("cotTheta_sigma", &(pb_cotTheta_sigma));
    ttree3_->SetBranchAddress("phi_mean"      , &(pb_phi_mean));
    ttree3_->SetBranchAddress("phi_sigma"     , &(pb_phi_sigma));
    ttree3_->SetBranchAddress("z0_mean"       , &(pb_z0_mean));
    ttree3_->SetBranchAddress("z0_sigma"      , &(pb_z0_sigma));

    ttree2_ = (TTree*) tfile_->Get("patternBankInfo");
    assert(ttree2_ != 0);

    ttree2_->SetBranchAddress("coverage"      , &(pb_coverage));
    ttree2_->SetBranchAddress("count"         , &(pb_count));
    ttree2_->SetBranchAddress("tower"         , &(pb_tower));
    ttree2_->SetBranchAddress("superstrip"    , &(pb_superstrip));
    //ttree2_->SetBranchAddress("superstrip_nx" , &(pb_superstrip_nx));
    //ttree2_->SetBranchAddress("superstrip_nz" , &(pb_superstrip_nz));

    ttree1_ = (TTree*) tfile_->Get("patternBank");
    assert(ttree1_ != 0);

    ttree1_->SetBranchAddress("frequency"     , &(pb_frequency));
    ttree1_->SetBranchAddress("superstripIds" , &(pb_superstripIds));
}
#endif  // _PATTERNBANKREADER_CXX_
