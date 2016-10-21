#include "TTTrackReader.h"

#include <iostream>

void analyzeTracks() {

    TString filename = "root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/SLHC/GEN/620_SLHC25p3_results/tt25_test/tracks_TTbar_PU140_sf1_nz8_L5x2_L10x2.root";  // from Olmo
    //TString filename = "root://xrootd2.ihepa.ufl.edu//store/user/jiafulow/SLHC/GEN/620_SLHC25p3_results/tt25_test/tracks_TTbar_PU140_sf1_nz8_L5x2_L10x2.root";  // from Olmo

    // Open the file
    std::cout << "Opening file..." << std::endl;

    TTTrackReader reader;
    reader.init(filename);

    // Get number of events
    //Long64_t nentries = reader.getEntries();
    Long64_t nentries = 10;
    std::cout << "Number of events: " << nentries << std::endl;

    // Loop over events
    // There is only one track per event
    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        // Retrieve the event
        reader.getEntry(ievt);

        std::cout << "ievt: " << ievt << std::endl;

        // Loop over tracking particles
        unsigned nparts = reader.vp2_pt->size();
        std::cout << "Num of tracking particles: " << nparts << std::endl;

        for (unsigned ipart = 0; ipart < nparts; ++ipart) {
            // Get the variables
            // See TTTrackReader.h for more info
            float pt      = reader.vp2_pt     ->at(ipart);
            float eta     = reader.vp2_eta    ->at(ipart);
            float phi     = reader.vp2_phi    ->at(ipart);
            float vx      = reader.vp2_vx     ->at(ipart);
            float vy      = reader.vp2_vy     ->at(ipart);
            float vz      = reader.vp2_vz     ->at(ipart);
            int   charge  = reader.vp2_charge ->at(ipart);
            int   pdgId   = reader.vp2_pdgId  ->at(ipart);
            bool  primary = reader.vp2_primary->at(ipart);

            if (!(primary && pt > 3.))  continue;
            if (!(M_PI/4<phi && phi<M_PI/2 && 0.<eta && eta<2.2/3))  continue;  // within trigger tower
            std::cout << ".. ipart: " << ipart
                      << " pdgId: "   << pdgId
                      << " pT: "      << pt
                      << " eta: "     << eta
                      << " phi: "     << phi
                      << " vz: "      << vz
                      << std::endl;
        }  // end loop ove tracking particles

/*
        // Loop over stubs
        unsigned nstubs = reader.vb_modId->size();
        std::cout << "Num of stubs: " << nstubs << std::endl;

        for (unsigned istub = 0; istub < nstubs; ++istub) {
            // Get the variables
            // See TTTrackReader.h for more info
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
*/

        // Loop over tracks
        unsigned ntracks = reader.vt_pt->size();
        std::cout << "Num of tracks: " << ntracks << std::endl;

        for (unsigned itrack = 0; itrack < ntracks; ++itrack) {
            // Get the variables
            // See TTTrackReader.h for more info
            float track_pt    = reader.vt_pt    ->at(itrack);
            float track_eta   = reader.vt_eta   ->at(itrack);
            float track_phi0  = reader.vt_phi0  ->at(itrack);
            float track_z0    = reader.vt_z0    ->at(itrack);
            float track_chi2  = reader.vt_chi2  ->at(itrack);
            int   track_ndof  = reader.vt_ndof  ->at(itrack);
            float track_invPt = reader.vt_invPt ->at(itrack);

            std::cout << ".. itrack: "  << itrack
                      << " pT: "        << track_pt
                      << " eta: "       << track_eta
                      << " phi0: "      << track_phi0
                      << " z0: "        << track_z0
                      << " norm chi2: " << track_chi2 / track_ndof
                      << " charge: "    << (track_invPt < 0. ? -1 : 1)
                      << std::endl;
        }  // end loop over tracks

        std::cout << std::endl;

    }  // end loop over events

}  // end analyze()
