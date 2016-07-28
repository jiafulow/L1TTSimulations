#ifndef _TTSTUBREADER_H_
#define _TTSTUBREADER_H_

#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

// _____________________________________________________________________________
// This is a simple wrapper around TChain. It sets the branch names, addresses
// etc. Its functions are essentially the same as the functions of TChain.

class TTStubReader {
public:
    TTStubReader(int verbose=1);
    ~TTStubReader();

    void init(TString src);

    Long64_t loadTree(Long64_t entry) { return tchain_->LoadTree(entry); }

    Int_t getEntry(Long64_t entry) { return tchain_->GetEntry(entry); }

    Long64_t getEntries() const { return tchain_->GetEntries(); }

    TChain* getChain() { return tchain_; }

    // genParticle information
    std::vector<float> *          vp_pt;        // particle pT [GeV]
    std::vector<float> *          vp_eta;       // particle eta
    std::vector<float> *          vp_phi;       // particle phi [rad]
    std::vector<float> *          vp_vx;        // particle vertex, x [cm]
    std::vector<float> *          vp_vy;        // particle vertex, y [cm]
    std::vector<float> *          vp_vz;        // particle vertex, z [cm]
    std::vector<int> *            vp_charge;    // particle charge

    // Stub information
    std::vector<float> *          vb_z;         // stub global z [cm]
    std::vector<float> *          vb_r;         // stub global r [cm]
    std::vector<float> *          vb_eta;       // stub global eta
    std::vector<float> *          vb_phi;       // stub global phi [rad]
    std::vector<float> *          vb_coordx;    // stub local phi
    std::vector<float> *          vb_coordy;    // stub local z (barrel) or r (endcap)
    std::vector<float> *          vb_trigBend;  // stub local bend distance (a.k.a. delta-s)
    std::vector<unsigned> *       vb_modId;     // stub module ID
    std::vector<int> *            vb_tpId;      // stub tracking particle ID

protected:
    TChain* tchain_;
    int treenumber_;
    const int verbose_;
};
#endif  // _TTSTUBREADER_H_

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation
#define _TTSTUBREADER_CXX_
#ifdef _TTSTUBREADER_CXX_

#include <cassert>
#include <stdexcept>

TTStubReader::TTStubReader(int verbose) :
    vp_pt               (0),
    vp_eta              (0),
    vp_phi              (0),
    vp_vx               (0),
    vp_vy               (0),
    vp_vz               (0),
    vp_charge           (0),
    //
    vb_z                (0),
    vb_r                (0),
    vb_eta              (0),
    vb_phi              (0),
    vb_coordx           (0),
    vb_coordy           (0),
    vb_trigBend         (0),
    vb_modId            (0),
    vb_tpId             (0),
    //
    verbose_(verbose) {}

TTStubReader::~TTStubReader() {
    if (tchain_)  delete tchain_;
}

void TTStubReader::init(TString src) {
    if (!src.EndsWith(".root") && !src.EndsWith(".txt")) {
        TString msg = "Input source must be either .root or .txt";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << src << std::endl;
    tchain_ = new TChain("ntupler/tree");

    if (src.EndsWith(".root")) {
        if (tchain_->Add(src) ) {
            //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
        } else {
            TString msg = "Failed to read " + src;
            throw std::invalid_argument(msg.Data());
        }

    } else if (src.EndsWith(".txt")) {
        TFileCollection fc("fileinfolist", "", src);
        if (tchain_->AddFileInfoList((TCollection*) fc.GetList()) ) {
            //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
        } else {
            TString msg = "Failed to read " + src;
            throw std::invalid_argument(msg.Data());
        }
    }

    assert(tchain_ != 0);
    treenumber_ = tchain_->GetTreeNumber();

    tchain_->SetBranchAddress("genParts_pt"       , &(vp_pt));
    tchain_->SetBranchAddress("genParts_eta"      , &(vp_eta));
    tchain_->SetBranchAddress("genParts_phi"      , &(vp_phi));
    tchain_->SetBranchAddress("genParts_vx"       , &(vp_vx));
    tchain_->SetBranchAddress("genParts_vy"       , &(vp_vy));
    tchain_->SetBranchAddress("genParts_vz"       , &(vp_vz));
    tchain_->SetBranchAddress("genParts_charge"   , &(vp_charge));
    //
    tchain_->SetBranchAddress("TTStubs_z"         , &(vb_z));
    tchain_->SetBranchAddress("TTStubs_r"         , &(vb_r));
    tchain_->SetBranchAddress("TTStubs_eta"       , &(vb_eta));
    tchain_->SetBranchAddress("TTStubs_phi"       , &(vb_phi));
    tchain_->SetBranchAddress("TTStubs_coordx"    , &(vb_coordx));
    tchain_->SetBranchAddress("TTStubs_coordy"    , &(vb_coordy));
    tchain_->SetBranchAddress("TTStubs_trigBend"  , &(vb_trigBend));
    tchain_->SetBranchAddress("TTStubs_modId"     , &(vb_modId));
    tchain_->SetBranchAddress("TTStubs_tpId"      , &(vb_tpId));
}
#endif  // _TTSTUBREADER_CXX_
