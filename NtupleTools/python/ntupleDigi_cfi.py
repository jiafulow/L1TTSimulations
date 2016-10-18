import FWCore.ParameterSet.Config as cms

ntupleTrackingParticles = cms.EDProducer('NtupleTrackingParticles',
    inputTag = cms.InputTag('mix', 'MergedTrackTruth'),
    prefix = cms.string('trkParts@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleTrackingVertices = cms.EDProducer('NtupleTrackingVertices',
    inputTag = cms.InputTag('mix', 'MergedTrackTruth'),
    prefix = cms.string('trkVertices@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleMixedSimHits = cms.EDProducer('NtupleMixedSimHits',
    inputTag = cms.InputTag('mix'),
    trkPartTag = cms.InputTag('mix', 'MergedTrackTruth'),
    simHitCollections = cms.PSet(
        muon = cms.VInputTag(
            #cms.InputTag('g4SimHits','MuonDTHits'),
            cms.InputTag('g4SimHits','MuonCSCHits'),
            #cms.InputTag('g4SimHits','MuonRPCHits'),
        ),
    ),
    prefix = cms.string('mixedSimHits@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleDigi = cms.Sequence(ntupleTrackingParticles * ntupleTrackingVertices)
#ntupleDigi = cms.Sequence(ntupleTrackingParticles * ntupleTrackingVertices * ntupleMixedSimHits)

