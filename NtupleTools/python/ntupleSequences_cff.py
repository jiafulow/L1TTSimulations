import FWCore.ParameterSet.Config as cms

from L1TTSimulations.NtupleTools.ntupleGen_cfi import *
#from L1TTSimulations.NtupleTools.ntupleGenExtra_cfi import *
#from L1TTSimulations.NtupleTools.ntupleSim_cfi import *
#from L1TTSimulations.NtupleTools.ntupleDigi_cfi import *
from L1TTSimulations.NtupleTools.ntupleL1TrackTrigger_cfi import *
from L1TTSimulations.NtupleTools.ntupleMaker_cfi import *

#ntupleGen += ntupleGenExtra

ntupleSequence = cms.Sequence(ntupleEventInfo * ntupleGen * ntupleL1TrackTrigger * ntupler)
#ntupleSequence = cms.Sequence(ntupleEventInfo * ntupleGen * ntupleSim * ntupleDigi * ntupleHLT * ntupleReco * ntupleL1TrackTrigger * ntupler)

#ntupleSequence_GENSIM = cms.Sequence(ntupleEventInfo * ntupleGen * ntupleSim * ntupler)
