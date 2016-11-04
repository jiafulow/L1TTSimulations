#ifndef AMSimulation_PatternMerging_h_
#define AMSimulation_PatternMerging_h_

#include <vector>
#include <map>
#include <unordered_set>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"


// Codes originally written by Luciano Ristori (FNAL) and developed by
// Roberto Rossin (Florida, now Padova)
// Modified for inclusion into AMSimulation

namespace slhcl1tt {

class PatternMerging {
public:

  // ___________________________________________________________________________
  // Sibling

  class Sibling {
  public:
    unsigned patternInd;
    unsigned siblingInd;
    int layer;
    int delta;
  };  // end class Sibling

  // ___________________________________________________________________________
  // Pattern

  class Pattern {
  public:
    // pattern proper (vector of superstrips)
    std::vector<unsigned> superstripIds;

    // popularity
    unsigned frequency;

    // attributes
    //float invPt_mean;
    //float invPt_sigma;
    //float cotTheta_mean;
    //float cotTheta_sigma;
    float phi_mean;
    float phi_sigma;
    //float z0_mean;
    //float z0_sigma;

    // index in original frequency-sorted list
    unsigned index;

    // Definition of Sibling
    // Requires two patterns to differ in one and only one layer
    // In case of success returns true and assigns layer number and delta
    // where delta is the distance of non-matching superstrip
    // layer and delta are undefined in case of false
    bool isSibling(const Pattern& p, int& layer, int& delta) const {
      int nDiff = 0;
      for (int i = 0; i < nLayers && nDiff <= 1; ++i) {
        if (superstripIds[i] != p.superstripIds[i]) {
          ++nDiff;
          layer = i;
          delta = (int) p.superstripIds[i] - (int) superstripIds[i];
        }
      }
      return (nDiff == 1); // one and only one layer
    }

  };  // end class Pattern

  // ___________________________________________________________________________
  // PatternMerging

  PatternMerging();
  ~PatternMerging();

  void mergePatterns(TString src, TString out, unsigned deltaN=36000, float targetCoverage=0.95) const;

  std::vector<unsigned> selectSiblings(
      unsigned patternInd,
      const std::vector<Sibling>& siblings,
      const std::vector<Pattern>& patternList,
      const std::vector<bool>& merged,
      const std::map<std::vector<unsigned>, unsigned>& patternMap
  ) const;

private:
  static const int nLayers = 6;

  int verbose_;
};

}  // namespace slhcl1tt

#endif
