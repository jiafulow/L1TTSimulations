#include "TTStubReader.h"

#include <iostream>

void analyzeStubs() {

    //TString filename = "root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/SLHC/GEN/Demo_Emulator/stubs_tt27_300M_emu.root";
    TString filename = "root://xrootd2.ihepa.ufl.edu//store/user/jiafulow/L1TrackTrigger/6_2_0_SLHC25p3/Demo_Emulator/stubs_tt27_300M_emu.root";

    // Open the file
    std::cout << "Opening file..." << std::endl;

    TTStubReader reader;
    reader.init(filename);

    // Get number of events
    //Long64_t nentries = reader.getEntries();
    Long64_t nentries = 100;
    std::cout << "Number of events: " << nentries << std::endl;

    // Loop over events
    // There is only one track per event
    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        // Retrieve the event
        reader.getEntry(ievt);

        // Get the variables
        // See TTStubReader.h for more info
        float pt     = reader.vp_pt    ->front();
        float eta    = reader.vp_eta   ->front();
        float phi    = reader.vp_phi   ->front();
        float vx     = reader.vp_vx    ->front();
        float vy     = reader.vp_vy    ->front();
        float vz     = reader.vp_vz    ->front();
        int   charge = reader.vp_charge->front();

        std::cout << "ievt: "     << ievt
                  << " muon pT: " << pt
                  << " eta: "     << eta
                  << " phi: "     << phi
                  << " vz: "      << vz
                  << std::endl;

        // Loop over stubs
        unsigned nstubs = reader.vb_modId->size();
        std::cout << "Num of stubs: " << nstubs << std::endl;
        for (unsigned istub = 0; istub < nstubs; ++istub) {
            // Get the variables
            // See TTStubReader.h for more info
            float    stub_z        = reader.vb_z       ->at(istub);
            float    stub_r        = reader.vb_r       ->at(istub);
            float    stub_eta      = reader.vb_eta     ->at(istub);
            float    stub_phi      = reader.vb_phi     ->at(istub);
            float    stub_coordx   = reader.vb_coordx  ->at(istub);
            float    stub_coordy   = reader.vb_coordy  ->at(istub);
            float    stub_trigBend = reader.vb_trigBend->at(istub);
            unsigned stub_modId    = reader.vb_modId   ->at(istub);
            int      stub_tpId     = reader.vb_tpId    ->at(istub);

            std::cout << ".. istub: "     << istub
                      << " module ID: "   << stub_modId
                      << " local phi: "   << stub_coordx
                      << " z: "           << stub_coordy
                      << " global phi: "  << stub_phi
                      << " z: "           << stub_z
                      << " r: "           << stub_r
                      << std::endl;

        }  // end loop over stubs

        std::cout << std::endl;

    }  // end loop over events

}  // end analyze()
