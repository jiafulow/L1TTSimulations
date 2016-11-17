#include "RoadMerging.h"

#include <cassert>
#include <iostream>

#include "PatternBankReader.h"
#include "TTTrackReader.h"
#include "TTTrackWriter.h"

using namespace slhcl1tt;

namespace slhcl1tt {
  // Print superstripIds
  std::ostream& operator<<(std::ostream& o, const std::vector<unsigned>& v) {
    for (const auto& x : v)  o << std::setw(6) << x;
    return o;
  }
}  // namespace slhcl1tt


RoadMerging::RoadMerging() : verbose_(1) {}

RoadMerging::~RoadMerging() {}

// _____________________________________________________________________________
void RoadMerging::process(TString bank, TString src, TString out) const {

  float targetCoverage = 0.95;

  // ___________________________________________________________________________
  // For reading pattern bank

  PatternBankReader pbreader(verbose_);
  pbreader.init(bank);

  // Get pattern bank info
  pbreader.getPatternInfo();

  unsigned npatterns   = pbreader.getEntries();
  unsigned statistics  = pbreader.pb_count;
  float    coverage    = pbreader.pb_coverage;

  pbreader.getEntry_toMerged();
  assert(pbreader.pb_indToMerged != 0);

  std::cout << "Original coverage: " << coverage << " with " << npatterns << " patterns and " << statistics << " statistics" << std::endl;


  // ___________________________________________________________________________
  // For reading

  TTTrackReader reader(verbose_);
  reader.init(src);

  // Get number of events
  unsigned nentries = 4294967295;  // std::numeric_limits<unsigned>::max()


  // ___________________________________________________________________________
  // For writing

  TTTrackWriter writer(verbose_);
  writer.init(reader.getChain(), out);

  std::vector<unsigned> vr_patternRef;   // [road i]
  std::vector<unsigned> vr_tower;        // [road i]
  std::vector<unsigned> vr_nstubs;       // [road i]
  std::vector<float>    vr_patternInvPt; // [road i]
  std::vector<unsigned> vr_patternFreq;  // [road i]

  std::vector<std::vector<unsigned> > vr_superstripIds;  // [road i][layer j]
  std::vector<std::vector<std::vector<unsigned> > > vr_stubRefs;  // [road i][layer j][stub k]

  std::vector<std::vector<std::vector<unsigned> > > vr_superstripIdsUnited;  // [road i][layer j][superstrip k]

  writer.getTree()->Branch("AMTTRoads_superstripIdsUnited", &vr_superstripIdsUnited);  // extra branch


  // ___________________________________________________________________________
  // Load patterns

  std::vector<Pattern> patterns;
  std::vector<Pattern> merged_patterns;

  unsigned totalFrequency = 0; // accumulate sum frequency
  float newCoverage = 0.;      // running coverage
  unsigned stoppingPatternInd = 0;  // num of patterns to reach targetCoverage

  for (unsigned ipatt = 0; ipatt < npatterns; ++ipatt) {  // index: unmerged pattern
    pbreader.getPattern(ipatt);
    pbreader.getPatternAttr(ipatt);

    assert(pbreader.pb_superstripIds->size() == nLayers);

    // Accumulate frequency
    totalFrequency += pbreader.pb_frequency;
    newCoverage = coverage * totalFrequency / statistics;

    if (!(newCoverage >= targetCoverage)) {
      ++stoppingPatternInd;
    }

    if (verbose_ && ipatt%200000 == 0)  std::cout << ".. Loading pattern: " << ipatt << " freq: " << pbreader.pb_frequency << " cov: " << newCoverage << std::endl;

    // Prepare temp pattern
    Pattern tempPattern;
    tempPattern.superstripIds  = *(pbreader.pb_superstripIds);
    tempPattern.frequency      = pbreader.pb_frequency;
    tempPattern.invPt_mean     = pbreader.pb_invPt_mean;
    tempPattern.index          = ipatt; // index in original frequency-sorted list
    tempPattern.indToMerged    = pbreader.pb_indToMerged->at(ipatt);  // index in merged pattern bank
    patterns.push_back(tempPattern);
  }

  if (verbose_)  std::cout << "Unmerged pattern bank size: " << patterns.size() << " totalFrequency: " << totalFrequency << std::endl;
  if (verbose_)  std::cout << "Unmerged pattern bank cov: " << targetCoverage << " stops at " << stoppingPatternInd << std::endl;

  assert(patterns.size() == npatterns);
  assert(totalFrequency == statistics);

  unsigned nmpatterns = pbreader.getTree_fromMerged()->GetEntries();

  for (unsigned jpatt = 0; jpatt < nmpatterns; ++jpatt) {  // index: merged pattern
    pbreader.getEntry_fromMerged(jpatt);

    assert(!pbreader.pb_indFromMerged->empty());

    if (verbose_ && jpatt%200000 == 0)  std::cout << ".. Loading merged pattern: " << jpatt << std::endl;

    const std::vector<unsigned>& indFromMergedTemp = *(pbreader.pb_indFromMerged);

    // Clone a temp pattern, update only the following variables
    Pattern tempPattern;
    tempPattern.superstripIds  = patterns.at(indFromMergedTemp.front()).superstripIds;
    tempPattern.frequency      = 0;
    tempPattern.invPt_mean     = 0.;
    tempPattern.index          = jpatt;  // index in merged pattern bank
    tempPattern.indToMerged    = patterns.at(indFromMergedTemp.front()).indToMerged;  // index in merged pattern bank
    tempPattern.superstripIdsUnited = std::vector<std::vector<unsigned> >(6);
    tempPattern.indFromMerged  = indFromMergedTemp;

    assert(tempPattern.index == tempPattern.indToMerged);  // a ha!

    // Loop over siblings

    for (unsigned is = 0; is < indFromMergedTemp.size(); ++is) {
      unsigned kpatt = indFromMergedTemp.at(is);

      // Incremental update formula with weighted entries
      auto   x    = patterns.at(kpatt).invPt_mean;
      auto   w    = patterns.at(kpatt).frequency;
      auto&& n    = tempPattern.frequency;
      auto&& mean = tempPattern.invPt_mean;

      n += w;
      mean += ((x - mean) * w)/n;

      // Make a union of superstripIds
      const auto& ssids   = patterns.at(kpatt).superstripIds;
      auto&&      m_ssids = tempPattern.superstripIdsUnited;

      assert(m_ssids.size() == ssids.size());
      for (unsigned ilayer = 0; ilayer < ssids.size(); ++ilayer) {
        if (std::find(m_ssids[ilayer].begin(), m_ssids[ilayer].end(), ssids[ilayer]) == m_ssids[ilayer].end()) {
          m_ssids[ilayer].push_back(ssids[ilayer]);
        }
      }
    }
    merged_patterns.push_back(tempPattern);
  }

  if (verbose_)  std::cout << "Merged pattern bank size: " << merged_patterns.size() << std::endl;

  assert(merged_patterns.size() == nmpatterns);

  // Truncate at targetCoverage

  if (verbose_)  std::cout << "Sorting merged pattern bank by frequency ..." << std::endl;

  // Sort merged_patterns by frequency
  std::stable_sort(merged_patterns.begin(), merged_patterns.end(), [](const Pattern& lhs, const Pattern& rhs) {
    return (lhs.frequency >= rhs.frequency);  // higher frequency first
  });

  totalFrequency = 0; // accumulate sum frequency
  newCoverage = 0.;   // running coverage
  stoppingPatternInd = 0;  // num of patterns to reach targetCoverage

  for (unsigned jpatt = 0; jpatt < nmpatterns; ++jpatt) {  // index: merged pattern
    const Pattern& tempPattern = merged_patterns.at(jpatt);

    // Accumulate frequency
    totalFrequency += tempPattern.frequency;
    newCoverage = coverage * totalFrequency / statistics;

    if (verbose_ && jpatt%200000 == 0)  std::cout << ".. Checking merged pattern: " << jpatt << " freq: " << tempPattern.frequency << " cov: " << newCoverage << std::endl;

    if (!(newCoverage >= targetCoverage)) {
      ++stoppingPatternInd;
    }

    // Update the old mergedPatternRef to the new mergedPatternRef after sorting by frequency
    auto&& mergedPatternRef = merged_patterns.at(jpatt).index;
    const auto& indFromMergedTemp = merged_patterns.at(jpatt).indFromMerged;

    mergedPatternRef = jpatt;  // new index in merged pattern bank
    for (unsigned is = 0; is < indFromMergedTemp.size(); ++is) {
      unsigned kpatt = indFromMergedTemp.at(is);
      auto&& indToMerged = patterns.at(kpatt).indToMerged;
      indToMerged = jpatt;  // new index in merged pattern bank
    }
  }

  if (verbose_)  std::cout << "Merged pattern bank size: " << merged_patterns.size() << " totalFrequency: " << totalFrequency << std::endl;
  if (verbose_)  std::cout << "Merged pattern bank cov: " << targetCoverage << " stops at " << stoppingPatternInd << std::endl;

  assert(totalFrequency == statistics);

  //merged_patterns.resize(stoppingPatternInd);


  // ___________________________________________________________________________
  // Loop over events

  std::vector<TTRoad> roads;
  std::vector<TTRoad> merged_roads;

  for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
    if (reader.loadTree(ievt) < 0)  break;
    reader.getEntry(ievt);

    const unsigned nroads = reader.vr_patternRef->size();

    if (verbose_ && ievt%100 == 0)  std::cout << ".. Processing event: " << ievt << " nroads: " << nroads << std::endl;

    roads.clear();
    merged_roads.clear();

    for (unsigned iroad=0; iroad<nroads; ++iroad) {

      //if (verbose_)  std::cout << ".... Processing road: " << iroad << std::endl;

      // Reconstruct the road
      TTRoad aroad;
      aroad.patternRef          = reader.vr_patternRef->at(iroad);
      aroad.tower               = reader.vr_tower->at(iroad);
      aroad.nstubs              = reader.vr_nstubs->at(iroad);
      aroad.patternInvPt        = reader.vr_patternInvPt->at(iroad);
      aroad.patternFreq         = reader.vr_patternFreq->at(iroad);
      aroad.superstripIds       = reader.vr_superstripIds->at(iroad);
      aroad.stubRefs            = reader.vr_stubRefs->at(iroad);
      aroad.superstripIdsUnited = std::vector<std::vector<unsigned> >();
      roads.push_back(aroad);
    }

    assert(roads.size() == nroads);

    // Sort roads by pT
    //std::stable_sort(roads.begin(), roads.end(), [](const TTRoad& lhs, const TTRoad& rhs) {
    //  return (std::abs(lhs.patternInvPt) < std::abs(rhs.patternInvPt));  // higher pT first
    //});

    // The real work is done here
    mergeRoads(patterns, merged_patterns, roads, merged_roads);

    // Sort merged_roads by pT
    std::stable_sort(merged_roads.begin(), merged_roads.end(), [](const TTRoad& lhs, const TTRoad& rhs) {
      return (std::abs(lhs.patternInvPt) < std::abs(rhs.patternInvPt));  // higher pT first
    });

    //if (verbose_)  std::cout << ".. Num of unmerged roads: " << roads.size() << " least invPt: " << (roads.size() ? roads.front().patternInvPt : std::nan("")) << std::endl;
    //if (verbose_)  std::cout << ".. Num of merged roads: " << merged_roads.size() << " least invPt: " << (merged_roads.size() ? merged_roads.front().patternInvPt : std::nan("")) << std::endl;


    const unsigned nmroads = merged_roads.size();

    // Split TTRoad class into branches
    vr_patternRef         .clear();
    vr_tower              .clear();
    vr_nstubs             .clear();
    vr_patternInvPt       .clear();
    vr_patternFreq        .clear();
    vr_superstripIds      .clear();
    vr_stubRefs           .clear();
    vr_superstripIdsUnited.clear();

    for (unsigned jroad=0; jroad<nmroads; ++jroad) {
      const TTRoad& amroad = merged_roads.at(jroad);
      assert(amroad.superstripIdsUnited.size() > 0);

      vr_patternRef         .emplace_back(amroad.patternRef         );
      vr_tower              .emplace_back(amroad.tower              );
      vr_nstubs             .emplace_back(amroad.nstubs             );
      vr_patternInvPt       .emplace_back(amroad.patternInvPt       );
      vr_patternFreq        .emplace_back(amroad.patternFreq        );
      vr_superstripIds      .emplace_back(amroad.superstripIds      );
      vr_stubRefs           .emplace_back(amroad.stubRefs           );
      vr_superstripIdsUnited.emplace_back(amroad.superstripIdsUnited);
    }

    // Write out
    *(reader.vr_patternRef         ) = vr_patternRef         ;
    *(reader.vr_tower              ) = vr_tower              ;
    *(reader.vr_nstubs             ) = vr_nstubs             ;
    *(reader.vr_patternInvPt       ) = vr_patternInvPt       ;
    *(reader.vr_patternFreq        ) = vr_patternFreq        ;
    *(reader.vr_superstripIds      ) = vr_superstripIds      ;
    *(reader.vr_stubRefs           ) = vr_stubRefs           ;
    //*(reader.vr_superstripIdsUnited) = vr_superstripIdsUnited;

    writer.fill();

  }  // end loop over events

  // Write out
  writer.writeTree();
}

// _____________________________________________________________________________
void RoadMerging::mergeRoads(
    const std::vector<Pattern>& patterns,
    const std::vector<Pattern>& merged_patterns,
    const std::vector<TTRoad>& roads,
    std::vector<TTRoad>& merged_roads
) const {

  std::map<unsigned, TTRoad> awesome_map;

  const unsigned nroads = roads.size();

  //unsigned stoppingPatternInd = merged_patterns.size();

  for (unsigned iroad=0; iroad<nroads; ++iroad) {
    const TTRoad& aroad = roads.at(iroad);

    unsigned patternRef       = aroad.patternRef;
    unsigned mergedPatternRef = patterns.at(patternRef).indToMerged;

    //if (mergedPatternRef >= stoppingPatternInd)
    //  continue;

    std::pair<std::map<unsigned, TTRoad>::iterator, bool> ins = awesome_map.insert(std::make_pair(mergedPatternRef, aroad));

    if (ins.second) {  // insert success
      // the first road with this mergedPatternRef
      auto&& amroad = ins.first->second;
      amroad.patternRef          = mergedPatternRef;
      //amroad.tower               = aroad.tower;
      //amroad.nstubs              = aroad.nstubs;
      amroad.patternInvPt        = merged_patterns.at(mergedPatternRef).invPt_mean;
      amroad.patternFreq         = merged_patterns.at(mergedPatternRef).frequency;
      //amroad.superstripIds       = aroad.superstripIds;
      //amroad.stubRefs            = aroad.stubRefs;
      amroad.superstripIdsUnited = merged_patterns.at(mergedPatternRef).superstripIdsUnited;

    } else {           // insert fail
      // other roads with this mergedPatternRef
      auto&& amroad = ins.first->second;

      // Make a union of stubRefs
      const auto& stubrefs   = aroad.stubRefs;
      auto&&      m_stubrefs = amroad.stubRefs;
      auto&&      m_nstubs   = amroad.nstubs;

      assert(m_stubrefs.size() == stubrefs.size());
      for (unsigned ilayer = 0; ilayer < stubrefs.size(); ++ilayer) {
        for (unsigned jstub = 0; jstub < stubrefs[ilayer].size(); ++jstub) {
          if (std::find(m_stubrefs[ilayer].begin(), m_stubrefs[ilayer].end(), stubrefs[ilayer][jstub]) == m_stubrefs[ilayer].end()) {
            m_stubrefs[ilayer].push_back(stubrefs[ilayer][jstub]);
          }
        }

        const unsigned maxStubs = 4;
        if (m_stubrefs[ilayer].size() > maxStubs) {
          //m_stubrefs[ilayer].resize(maxStubs);  // keep first N stubs

          m_stubrefs[ilayer].erase(m_stubrefs[ilayer].begin(), m_stubrefs[ilayer].end() - maxStubs);  // keep last N stubs
          assert(m_stubrefs[ilayer].size() == maxStubs);
        }

        m_nstubs += m_stubrefs[ilayer].size();
      }
    }
  }

  std::map<unsigned, TTRoad>::const_iterator awesome_map_it  = awesome_map.begin();
  std::map<unsigned, TTRoad>::const_iterator awesome_map_end = awesome_map.end();

  for (; awesome_map_it != awesome_map_end; ++awesome_map_it) {
    merged_roads.push_back(awesome_map_it->second);
  }
}
