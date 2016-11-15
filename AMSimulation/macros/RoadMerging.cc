#include "RoadMerging.h"

#include <cassert>
#include <iostream>

#include "PatternBankReader.h"
#include "TTTrackReader.h"

using namespace slhcl1tt;

RoadMerging::RoadMerging() : verbose_(1) {}

RoadMerging::~RoadMerging() {}

// _____________________________________________________________________________
void RoadMerging::process(TString src, TString bank) const {
    // Open the file
    std::cout << "Opening file..." << std::endl;

    TTTrackReader reader;
    reader.init(src);

    // Get number of events
    //Long64_t nentries = reader.getEntries();
    Long64_t nentries = 10;
    std::cout << "Number of events: " << nentries << std::endl;

    // Loop over events
    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        // Retrieve the event
        reader.getEntry(ievt);
    }
}

// _____________________________________________________________________________
void RoadMerging::mergeRoads() const {

}

void RoadMerging::mergeRoad() const {

}
