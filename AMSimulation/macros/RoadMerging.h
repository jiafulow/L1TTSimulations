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
  struct TTRoad {
    unsigned patternRef;
    unsigned tower;
    unsigned nstubs;
    float    patternInvPt;

    std::vector<unsigned> superstripIds;
    std::vector<std::vector<unsigned> > stubRefs;  // stubRefs[superstrip i][stub j]

    std::vector<std::vector<unsigned> > superstripIdsBigLeague;
  };

  // ___________________________________________________________________________
  // RoadMerging

  RoadMerging();
  ~RoadMerging();

  void process(TString src, TString bank) const;  //TODO: factor out this guy

  void mergeRoads() const;
  void mergeRoad() const;

private:
  int verbose_;
};

}  // namespace slhcl1tt

#endif
