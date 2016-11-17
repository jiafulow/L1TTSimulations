#ifndef AMSimulation_PatternMerging_h_
#define AMSimulation_PatternMerging_h_

#include <vector>
#include <map>

#include "TString.h"


// Codes originally written by Roberto Rossin (Florida, now Padova)
// Modified for inclusion into AMSimulation

namespace slhcl1tt {

class RoadMerging {
public:

  // ___________________________________________________________________________
  // TTRoad

  struct TTRoad {
    unsigned patternRef;
    unsigned tower;
    unsigned nstubs;
    float    patternInvPt;
    unsigned patternFreq;

    std::vector<unsigned> superstripIds;
    std::vector<std::vector<unsigned> > stubRefs;  // [layer i][stub j]

    std::vector<std::vector<unsigned> > superstripIdsUnited;  // [layer i][superstrip j]
  };

  // ___________________________________________________________________________
  // Pattern

  struct Pattern {
    // pattern proper (vector of superstrips)
    std::vector<unsigned> superstripIds;

    // popularity
    unsigned frequency;

    // attributes
    float invPt_mean;
    //float invPt_sigma;
    //float cotTheta_mean;
    //float cotTheta_sigma;
    //float phi_mean;
    //float phi_sigma;
    //float z0_mean;
    //float z0_sigma;

    // index in original frequency-sorted list
    unsigned index;

    // index in merged pattern bank
    unsigned indToMerged;

    // pattern proper big-league (vector of vector of superstrips)
    std::vector<std::vector<unsigned> > superstripIdsUnited;

    // sibling indices in original frequency-sorted list
    std::vector<unsigned> indFromMerged;
  };

  // ___________________________________________________________________________
  // RoadMerging

  RoadMerging();
  ~RoadMerging();

  void process(TString bank, TString src, TString out) const;  //TODO: factor out this guy

  void mergeRoads(
      const std::vector<Pattern>& patterns,
      const std::vector<Pattern>& merged_patterns,
      const std::vector<TTRoad>& roads,
      std::vector<TTRoad>& merged_roads
  ) const;

private:
  static const int nLayers = 6;

  int verbose_;
};

}  // namespace slhcl1tt

#endif
