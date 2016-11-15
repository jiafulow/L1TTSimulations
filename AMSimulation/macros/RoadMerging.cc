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

  // ___________________________________________________________________________
  // For reading pattern bank

  PatternBankReader pbreader(verbose_);
  pbreader.init(bank);

  pbreader.getEntry_toMerged();
  assert(pbreader.pb_indToMerged != 0);

  unsigned npatterns = pbreader.getEntries();


  // ___________________________________________________________________________
  // For reading

  TTTrackReader reader(verbose_);
  reader.init(src);

  // Get number of events
  //unsigned nentries = TChain::kBigNumber;
  unsigned nentries = 10;


  // ___________________________________________________________________________
  // Load patterns

  std::vector<Pattern> patterns;

  std::vector<unsigned> indToMerged;
  std::vector<std::vector<unsigned> > indFromMerged;

  for (unsigned ipatt = 0; ipatt < npatterns; ++ipatt) {
    pbreader.getPattern(ipatt);
    pbreader.getPatternAttr(ipatt);

    if (verbose_ && ipatt%200000 == 0)  std::cout << ".. Loading pattern: " << ipatt << std::endl;

    // Prepare temp pattern
    Pattern tempPattern;
    tempPattern.superstripIds  = *(pbreader.pb_superstripIds);
    tempPattern.frequency      = pbreader.pb_frequency;
    tempPattern.invPt_mean     = pbreader.pb_invPt_mean;
    tempPattern.index          = ipatt; // index in original frequency-sorted list
    patterns.push_back(tempPattern);

    unsigned indToMergedTemp   = pbreader.pb_indToMerged->at(ipatt);
    indToMerged.push_back(indToMergedTemp);
  }

  assert(patterns.size() == npatterns);
  assert(indToMerged.size() == npatterns);

  unsigned nmpatterns = pbreader.getTree_fromMerged()->GetEntries();

  for (unsigned jpatt = 0; jpatt < nmpatterns; ++jpatt) {
    pbreader.getEntry_fromMerged(jpatt);

    if (verbose_ && jpatt%200000 == 0)  std::cout << ".. Loading merged pattern: " << jpatt << std::endl;

    std::vector<unsigned> indFromMergedTemp = *(pbreader.pb_indFromMerged);
    indFromMerged.push_back(indFromMergedTemp);
  }

  assert(indFromMerged.size() == nmpatterns);


  // ___________________________________________________________________________
  // Loop over events

  std::vector<TTRoad> roads;
  std::vector<TTRoad> merged_roads;

  for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
    if (reader.loadTree(ievt) < 0)  break;
    reader.getEntry(ievt);

    const unsigned nroads = reader.vr_patternRef->size();

    if (verbose_)  std::cout << ".. Processing event: " << ievt << " nroads: " << nroads << std::endl;

    roads.clear();
    merged_roads.clear();
    patterns.clear();

    for (unsigned iroad=0; iroad<nroads; ++iroad) {

      if (verbose_)  std::cout << ".... Processing road: " << iroad << std::endl;

      // Reconstruct the road
      TTRoad aroad;
      aroad.patternRef             = reader.vr_patternRef->at(iroad);
      aroad.tower                  = reader.vr_tower->at(iroad);
      aroad.nstubs                 = reader.vr_nstubs->at(iroad);
      aroad.patternInvPt           = reader.vr_patternInvPt->at(iroad);
      aroad.patternFreq            = reader.vr_patternFreq->at(iroad);
      aroad.superstripIds          = reader.vr_superstripIds->at(iroad);
      aroad.stubRefs               = reader.vr_stubRefs->at(iroad);
      aroad.superstripIdsBigLeague = std::vector<std::vector<unsigned> >();
      roads.push_back(aroad);
    }

    assert(roads.size() == nroads);

    mergeRoads(patterns, indToMerged, indFromMerged, roads, merged_roads);
  }
}

// _____________________________________________________________________________
void RoadMerging::mergeRoads(
    const std::vector<Pattern>& patterns,
    const std::vector<unsigned>& indToMerged,
    const std::vector<std::vector<unsigned> >& indFromMerged,
    const std::vector<TTRoad>& roads,
    std::vector<TTRoad>& merged_roads
) const {

  const unsigned nroads = roads.size();

  for (unsigned iroad=0; iroad<nroads; ++iroad) {
    const TTRoad& aroad = roads.at(iroad);

  }

}
