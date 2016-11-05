#include "PatternMerging.h"

#include <cassert>
#include <iostream>
#include <iomanip>

#include "PatternBankReader.h"

using namespace slhcl1tt;

namespace slhcl1tt {
  // Print superstripIds
  std::ostream& operator<<(std::ostream& o, const std::vector<unsigned>& v) {
    for (const auto& x : v)  o << std::setw(6) << x;
    return o;
  }
}  // namespace slhcl1tt


PatternMerging::PatternMerging() : verbose_(1) {}

PatternMerging::~PatternMerging() {}

// _____________________________________________________________________________
void PatternMerging::mergePatterns(TString src, TString out, unsigned deltaN, float targetCoverage) const {

  //const unsigned deltaN         = 36000;  // search window (0 = all)
  //const float    targetCoverage = 0.95;   // get patterns until coverage >= targetCoverage
  const unsigned maxPatterns    = 0;      // maximum number of patterns to read from bank (0 = all)
  const unsigned maxTrials      = 0;      // number of patterns to process (0 = all)
  //const bool     randomize      = false;  // random sampling of patterns
  const int      nInfo          = 1000;   // interval for printed info

  // Book histograms

  if (verbose_) std::cout << "Book histograms" << std::endl;

  TH1::SetDefaultSumw2();
  TH1F * hNSiblings = new TH1F("NSiblings", "NSiblings", 41, -0.5, 40.5);
  TH1F * hNMerged   = new TH1F("NMerged", "NMerged", 8, 0.5, 8.5);
  TH1F * hLayer     = new TH1F("Layer", "Layer", 6, -0.5, 5.5);
  TH1F * hDeltaSS   = new TH1F("DeltaSS", "DeltaSS", 2001, -1000.5, +1000.5);
  TH1F * hDeltaN    = new TH1F("DeltaN", "DeltaN", 100, 0., -1.);  // automatic range
  TH1F * hDeltaPhi  = new TH1F("DeltaPhi", "DeltaPhi", 100, 0., -1.);  // automatic range


  // ___________________________________________________________________________
  // Open the pattern bank root file

  PatternBankReader pbreader(verbose_);
  pbreader.init(src);

  // Get pattern bank info
  pbreader.getPatternInfo();

  unsigned npatterns   = pbreader.getEntries();
  unsigned statistics  = pbreader.pb_count;
  float    coverage    = pbreader.pb_coverage;
  unsigned tower       = pbreader.pb_tower;
  //int      magicNumber = pbreader.pb_superstrip_nx;
  int      magicNumber = 243;  //FIXME

  if (verbose_) {
    std::cout << "Coverage: " << coverage;
    std::cout << " with " << npatterns << " patterns. ";
    std::cout << " Frequency sum: " << statistics;
    std::cout << " tower: " << tower;
    std::cout << " magicNumber: " << magicNumber;
    std::cout << std::endl;
  }


  // ___________________________________________________________________________
  // Load patterns into local data structure

  if (verbose_) std::cout << "Loading patterns..." << std::endl;

  if (maxPatterns) npatterns = maxPatterns;

  std::vector<Pattern> patternList; // save all patterns
  patternList.reserve(npatterns);

  unsigned totalFrequency = 0; // accumulate sum frequency
  float newCoverage = 0.0;     // running coverage

  // Loop to read all patterns and build vector patternList

  for (unsigned ipatt = 0; ipatt < npatterns; ++ipatt) {
    // Get pattern
    pbreader.getPattern(ipatt);
    pbreader.getPatternAttr(ipatt);

    assert(pbreader.pb_superstripIds->size() == nLayers);

    // Accumulate frequency
    totalFrequency += pbreader.pb_frequency;
    newCoverage = coverage * totalFrequency / statistics;

    if (verbose_ && (ipatt%1000000) == 0) {
      std::cout << std::setw(10) << ipatt << " patterns. ";
      std::cout << " Current coverage: " << newCoverage;
      std::cout << " frequency: " << pbreader.pb_frequency;
      std::cout << std::endl;
    }

    if (verbose_ >= 100) { // print pattern
      std::cout << ipatt << ": ";
      std::cout << std::setw(6)  << (*pbreader.pb_superstripIds);
      std::cout << std::setw(10) << pbreader.pb_phi_mean;
      std::cout << std::setw(6)  << pbreader.pb_frequency;
      std::cout << std::setw(15) << newCoverage;
      std::cout << std::endl;
    }

    // Prepare temp pattern

    Pattern tempPattern;
    tempPattern.superstripIds  = *(pbreader.pb_superstripIds);
    tempPattern.frequency      = pbreader.pb_frequency;
    //tempPattern.invPt_mean     = pbreader.pb_invPt_mean;
    //tempPattern.invPt_sigma    = pbreader.pb_invPt_sigma;
    //tempPattern.cotTheta_mean  = pbreader.pb_cotTheta_mean;
    //tempPattern.cotTheta_sigma = pbreader.pb_cotTheta_sigma;
    tempPattern.phi_mean       = pbreader.pb_phi_mean;
    tempPattern.phi_sigma      = pbreader.pb_phi_sigma;
    //tempPattern.z0_mean        = pbreader.pb_z0_mean;
    //tempPattern.z0_sigma       = pbreader.pb_z0_sigma;
    tempPattern.index          = ipatt; // index in original frequency-sorted list

    assert(tempPattern.superstripIds.size() == nLayers);
    assert(patternList.size() == ipatt);

    if (verbose_ >= 100) {  // print temp pattern
      std::cout << ipatt << ": ";
      std::cout << std::setw(6)  << tempPattern.superstripIds;
      std::cout << std::setw(10) << tempPattern.phi_mean;
      std::cout << std::setw(6)  << tempPattern.frequency;
      std::cout << std::setw(15) << newCoverage;
      std::cout << std::endl;
    }

    // append temp pattern to list
    patternList.push_back(tempPattern);

    // truncate at target coverage
    if (newCoverage >= targetCoverage) break;
  }  // end loop on patterns

  npatterns = patternList.size(); // number of patterns actually loaded for target coverage

  if (verbose_) {
    std::cout << "Coverage: " << newCoverage;
    std::cout << " with " << npatterns << " patterns. ";
    std::cout << " Frequency sum: " << totalFrequency;
    std::cout << std::endl;
  }

  if (!maxPatterns) assert(coverage < targetCoverage || newCoverage >= targetCoverage); // did we reach target coverage?

  // Create map from pattern to index in patternList

  if (verbose_) std::cout << "Setting up patternMap..." << std::endl;

  std::map<std::vector<unsigned>, unsigned> patternMap; // map pattern -> index

  for (const auto& tempPattern : patternList) {
    patternMap[tempPattern.superstripIds] = tempPattern.index;
  }


  // ___________________________________________________________________________
  // Clone and sort patterns by phi

  if (verbose_) std::cout << "Cloning patternList..." << std::endl;

  std::vector<Pattern> clonePatternList = patternList; // clone all patterns

  if (verbose_) std::cout << "Sorting clonePatternList by phi..." << std::endl;

  std::sort(clonePatternList.begin(), clonePatternList.end(), [](const Pattern& a, const Pattern& b) {
    return a.phi_mean < b.phi_mean;  // sort patterns by smaller phi first
  });

  // Set all indices to phi-sorted list in frequency-sorted list

  for (unsigned iq = 0; iq < npatterns; ++iq) {   // iq is index in phi-sorted list
    unsigned ip = clonePatternList.at(iq).index;  // ip is index in frequency-sorted list
    patternList.at(ip).index = iq;
    assert(clonePatternList.at(iq).frequency == patternList.at(ip).frequency);
  }

  if (verbose_) std::cout << "Finished sorting" << std::endl;


  // ___________________________________________________________________________
  //////////////////////////////////////////////////////////////
  // Merging starts here
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Vectors for pattern merging
  //////////////////////////////////////////////////////////////
  std::vector<unsigned>               indToMerged(npatterns, 0);
  std::vector<std::vector<unsigned> > indFromMerged;

  // flags for already merged patterns, all patterns start as not merged
  std::vector<bool>                   merged(npatterns, false);

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////


  // Adjust number of patterns to process

  unsigned nTrials = maxTrials;
  if (maxTrials == 0) nTrials = npatterns; // if maxTrials = 0, process all patterns in bank

  std::cout << "Begin merging..." << std::endl;

  for (unsigned iTrial = 0; iTrial < nTrials; ++iTrial) {
    unsigned ip = iTrial; // go through patterns sequentially - original list ordering

    // print some information at fixed intervals (nInfo)
    if (verbose_ && (iTrial%nInfo) == 0) {
      std::cout << std::setw(10) << iTrial << " patterns ";
      std::cout << std::setw(10) << indFromMerged.size() << " merged    ";
      std::cout << std::endl;
    }

    if (merged[ip]) continue; // skip this pattern if already merged

    if (verbose_ >= 10) { // print pattern number
      std::cout << std::endl << "-----" << std::setw(10) << ip << " " << patternList[ip].superstripIds << std::endl;
    }

    // Define range for sibling search
    // (deltaN = 0 means all patterns)
    const unsigned centralN = patternList[ip].index;
    const unsigned beginN   = (deltaN == 0 || centralN < deltaN) ? 0 : (centralN-deltaN);
    const unsigned endN     = (deltaN == 0 || (centralN+deltaN) > npatterns) ? npatterns : (centralN+deltaN);

    // find siblings for this pattern patternList[ip] in the search range

    std::vector<Sibling> siblings;

    // iq points into phi ordered list clonePatternList[iq]

    for (unsigned iq = beginN; iq < endN; ++iq) {
      if (iq == centralN) continue;  // skip if this pattern is itself

      if (merged[clonePatternList[iq].index]) continue; // skip this pattern if already merged

      int layer = 0, deltass = 0;

      // Is sibling?
      if (patternList[ip].isSibling(clonePatternList[iq], layer, deltass) ) {

        if (verbose_ >= 10) { // print sibling pattern
          std::cout << "     " << std::setw(10) << clonePatternList[iq].index << " " << clonePatternList[iq].superstripIds << std::endl;
        }

        // fill histograms
        float weight = patternList[ip].frequency;
        hLayer->Fill(layer, weight);
        hDeltaSS->Fill(deltass, weight);
        hDeltaN->Fill((float) iq - (float) centralN);
        hDeltaPhi->Fill(clonePatternList[iq].phi_mean - patternList[ip].phi_mean);

        // construct a sibling
        if (abs(deltass) == 1 || abs(deltass) == magicNumber) {
          Sibling sib;
          //sib.patternInd = ip;
          //sib.siblingInd = clonePatternList[iq].index;
          sib.index      = clonePatternList[iq].index;
          sib.frequency  = clonePatternList[iq].frequency;
          sib.layer      = layer;
          sib.delta      = deltass;
          siblings.push_back(sib);
        }
      }
    }  // end search range

    // fill histogram: number of siblings
    hNSiblings->Fill(siblings.size());


    // Interrogate selectSiblings() to pick patterns to merge

    merged[ip] = true; // mark this pattern as merged before selectSiblings() to avoid combining itself

    const std::vector<unsigned>& selectedSiblings = selectSiblings(ip, siblings, patternList, merged, patternMap);

    const int nm = selectedSiblings.size() + 1;
    assert(nm == 1 || nm == 2 || nm == 4 || nm == 8);

    // fill histogram: number of merged patterns
    hNMerged->Fill(nm);

    // write index of new indFromMerged element for original pattern
    indToMerged[ip] = indFromMerged.size();

    // temp vector to be appended to indFromMerged
    std::vector<unsigned> indFromMergedTemp;
    indFromMergedTemp.push_back(ip); // append this pattern

    // append all selected siblings and set index to merged for all of them
    for (unsigned is = 0; is < selectedSiblings.size(); ++is) {
      merged[selectedSiblings[is]] = true; // mark sibling as merged
      indToMerged[selectedSiblings[is]] = indFromMerged.size();
      indFromMergedTemp.push_back(selectedSiblings[is]);
    }

    // append temp vector to indFromMerged
    indFromMerged.push_back(indFromMergedTemp);
  }  // end loop on patterns


  if (verbose_ >= 10) {
     // Dump part of merged structure vectors for diagnostics
    for (unsigned i = 0; i < indFromMerged.size(); ++i) {
      if (i > 100) break;
      std::cout << i << ": ";
      for (unsigned j = 0; j < indFromMerged[i].size(); ++j) {
        std::cout << std::setw(10) << indFromMerged[i][j];
        std::cout << std::setw(10) << indToMerged[indFromMerged[i][j]];
        std::cout << "     ";
        for (int k = 0; k < nLayers; ++k) {
          std::cout << std::setw(4) << patternList[indFromMerged[i][j]].superstripIds[k];
        }
      }
      std::cout << std::endl;
    }
  }


  // ___________________________________________________________________________
  // Release memory

  {
    std::vector<Pattern>().swap(patternList);
    std::vector<Pattern>().swap(clonePatternList);
    std::map<std::vector<unsigned>, unsigned>().swap(patternMap);
  }


  // ___________________________________________________________________________
  // Consistency checks
  // 1. all indices in indFromMerged must have no repetition
  //    and cover the whole original pattern bank

  {
    std::cout << "Doing consistency checks..." << std::endl;

    std::unordered_set<unsigned> allIndices;
    unsigned NIndices = 0;

    for (unsigned i = 0; i < indFromMerged.size(); ++i) {
      for (unsigned j = 0; j < indFromMerged[i].size(); ++j) {
        allIndices.insert(indFromMerged[i][j]);
        ++NIndices;
      }
    }

    if (verbose_) {
      // Write final summary
      double gain = (double)(npatterns - indFromMerged.size())/npatterns;
      std::cout << allIndices.size() << " different indices out of " << NIndices << " in indFromMerged" << std::endl;
      std::cout << indFromMerged.size() << " merged patterns out of " << npatterns << " (- " << 100*gain << "%)" << std::endl;
    }
    assert(NIndices == allIndices.size());
  }

  // ___________________________________________________________________________
  // Add merged structure vectors to pattern bank file

  if (verbose_) std::cout << "Adding merged structure vectors to pattern file..." << std::endl;

  TFile* outfile = TFile::Open(out, "RECREATE");

  TTree * t1 = new TTree("toMerged", "toMerged");
  t1->Branch("indToMerged", &indToMerged);
  t1->Fill();
  t1->Write();

  TTree * t2 = new TTree("fromMerged", "fromMerged");
  std::vector<unsigned> indFromMergedTemp;
  t2->Branch("indFromMerged", &indFromMergedTemp);
  for (unsigned i = 0; i < indFromMerged.size(); ++i) {
    indFromMergedTemp = indFromMerged[i];
    t2->Fill();
  }
  t2->Write();

  pbreader.getTree()->CloneTree()->Write();
  pbreader.getInfoTree()->CloneTree()->Write();
  pbreader.getAttrTree()->CloneTree()->Write();

  outfile->mkdir("histograms")->cd();
  hNSiblings->Write();
  hNMerged->Write();
  hLayer->Write();
  hDeltaSS->Write();
  hDeltaN->Write();
  hDeltaPhi->Write();
  outfile->Close();

  if (verbose_) std::cout << "Wrote " << out << std::endl;

  std::cout << "The End" << std::endl;
}

// _____________________________________________________________________________
std::vector<unsigned> PatternMerging::selectSiblings(
    unsigned patternInd,
    const std::vector<Sibling>& siblings,
    const std::vector<Pattern>& patternList,
    const std::vector<bool>& merged,
    const std::map<std::vector<unsigned>, unsigned>& patternMap
) const {
  // list of selected siblings to be returned
  // indices point to patternList
  std::vector<unsigned> selectedSiblings;

  if (siblings.size() == 0) return selectedSiblings;

  unsigned maxTotalFrequency = 0, totalFrequency = 0;

  // Declare all the variables
  unsigned patternInd1 = 0, patternInd2 = 0, patternInd3 = 0;
  unsigned freq1 = 0, freq2 = 0, freq3 = 0;
  unsigned patternIndNew = 0, freqNew = 0;
  std::map<std::vector<unsigned>, unsigned>::const_iterator found;


  ///////////////////////////////////////////////////////////////////
  // try merge x8
  ///////////////////////////////////////////////////////////////////

  // build all triplets of siblings (is1,is2,is3)

  if (siblings.size() >= 3) {
    std::vector<std::vector<unsigned> > newPatternList;
    std::vector<unsigned> newPatternIndList;

    std::vector<unsigned> newPattern0, newPattern1, newPattern2, newPattern3;

    for (unsigned is1 = 0; is1+2 < siblings.size(); ++is1) {
      patternInd1 = siblings[is1].index;
      freq1 = siblings[is1].frequency;

      for (unsigned is2 = is1+1; is2+1 < siblings.size(); ++is2) {
        patternInd2 = siblings[is2].index;
        freq2 = siblings[is2].frequency;

        for (unsigned is3 = is2+1; is3 < siblings.size(); ++is3) {
          patternInd3 = siblings[is3].index;
          freq3 = siblings[is3].frequency;

          // for each triplet of siblings, build a vector with 4 combined patterns
          bool goodTriplet = true;

          newPatternList.clear();

          newPattern0 = patternList[patternInd].superstripIds;
          newPattern0[siblings[is1].layer] += siblings[is1].delta;
          newPattern0[siblings[is2].layer] += siblings[is2].delta;
          newPattern0[siblings[is3].layer] += siblings[is3].delta;
          newPatternList.push_back(newPattern0);

          newPattern1 = newPattern0;
          newPattern1[siblings[is1].layer] -= siblings[is1].delta;
          newPatternList.push_back(newPattern1);

          newPattern2 = newPattern0;
          newPattern2[siblings[is2].layer] -= siblings[is2].delta;
          newPatternList.push_back(newPattern2);

          newPattern3 = newPattern0;
          newPattern3[siblings[is3].layer] -= siblings[is3].delta;
          newPatternList.push_back(newPattern3);

          assert(newPatternList.size() == 4);

          // check if all the combined patterns in the list belong to the pattern bank
          // and build totalFrequency as we go

          newPatternIndList.clear();

          totalFrequency = freq1 + freq2 + freq3;

          // loop on all 4 new combined patterns

          for (unsigned iNew = 0; iNew < 4; ++iNew) {

            // find combined pattern in patternMap
            found = patternMap.find(newPatternList[iNew]);
            if (found == patternMap.end()) { // new pattern does not exist
              goodTriplet = false;
              break;
            }

            patternIndNew = found->second;
            if (merged[patternIndNew]) { // new pattern is already merged
              goodTriplet = false;
              break;
            }

            freqNew = patternList[patternIndNew].frequency;
            totalFrequency += freqNew;
            newPatternIndList.push_back(patternIndNew);
          } // end loop on combined patterns


          if (goodTriplet && totalFrequency > maxTotalFrequency) {
            assert(newPatternIndList.size() == 4);
            maxTotalFrequency = totalFrequency;
            selectedSiblings.clear();
            selectedSiblings.push_back(patternInd1);
            selectedSiblings.push_back(patternInd2);
            selectedSiblings.push_back(patternInd3);
            selectedSiblings.push_back(newPatternIndList[0]);
            selectedSiblings.push_back(newPatternIndList[1]);
            selectedSiblings.push_back(newPatternIndList[2]);
            selectedSiblings.push_back(newPatternIndList[3]);
          }
        } // end loop on is3
      } // end loop on is2
    } // end loop on is1
  } // end if having 3 or more siblings

  if (selectedSiblings.size() > 0) return selectedSiblings;


  ///////////////////////////////////////////////////////////////////
  // merge x8 did not work. Now try merge x4
  ///////////////////////////////////////////////////////////////////

  // build all pairs of siblings (is1,is2)

  maxTotalFrequency = 0;

  if (siblings.size() >= 2) {
    std::vector<unsigned> newPattern;

    for (unsigned is1 = 0; is1+1 < siblings.size(); ++is1) {
      patternInd1 = siblings[is1].index;
      freq1 = siblings[is1].frequency;

      for (unsigned is2 = is1+1; is2 < siblings.size(); ++is2) {
        patternInd2 = siblings[is2].index;
        freq2 = siblings[is2].frequency;

        // combine two siblings and make a new combined pattern
        newPattern = patternList[patternInd].superstripIds;
        newPattern[siblings[is1].layer] += siblings[is1].delta;
        newPattern[siblings[is2].layer] += siblings[is2].delta;

        // find combined pattern in patternMap
        found = patternMap.find(newPattern);
        if (found == patternMap.end()) { // new pattern does not exist
          continue;
        }

        patternIndNew = found->second;
        if (merged[patternIndNew]) { // new pattern is already merged
          continue;
        }

        freqNew = patternList[patternIndNew].frequency;
        totalFrequency = freq1 + freq2 + freqNew;

        if (totalFrequency > maxTotalFrequency) {
          maxTotalFrequency = totalFrequency;
          selectedSiblings.clear();
          selectedSiblings.push_back(patternInd1);
          selectedSiblings.push_back(patternInd2);
          selectedSiblings.push_back(patternIndNew);
        }
      } // end loop on is2
    } // end loop on is1
  } // end if having 2 or more siblings

  if (selectedSiblings.size() > 0) return selectedSiblings;


  ///////////////////////////////////////////////////////////////////
  // merge x4 did not work. Now try merge x2
  ///////////////////////////////////////////////////////////////////

  // look for max popularity sibling

  maxTotalFrequency = 0;

  for (unsigned is1 = 0; is1 < siblings.size(); ++is1) {
    patternInd1 = siblings[is1].index;
    freq1 = siblings[is1].frequency;
    totalFrequency = freq1;

    if (totalFrequency > maxTotalFrequency) {
      maxTotalFrequency = totalFrequency;
      selectedSiblings.clear();
      selectedSiblings.push_back(patternInd1);
    }
  } // end loop on is1

  return selectedSiblings;
}
