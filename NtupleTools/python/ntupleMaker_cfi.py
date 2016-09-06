import FWCore.ParameterSet.Config as cms

ntupleEventInfo = cms.EDProducer("NtupleEventInfo",
    prefix = cms.string(''),
    suffix = cms.string(''),
)

ntupler = cms.EDAnalyzer("NtupleMaker",
    treeName = cms.string('tree'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
