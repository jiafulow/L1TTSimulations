import FWCore.ParameterSet.Config as cms

ntupleGenParticles = cms.EDProducer('NtupleGenParticles',
    inputTag = cms.InputTag('genParticles'),
    prefix = cms.string('genParts@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleGenJets = cms.EDProducer('NtupleGenJets',
    inputTag = cms.InputTag('ak4GenJets'),
    prefix = cms.string('genJets@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleGenMET = cms.EDProducer('NtupleGenMET',
    inputTag = cms.InputTag('genMetTrue'),
    prefix = cms.string('genMET@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

#from L1TTSimulations.NtupleTools.simBeamSpot_cfi import simBeamSpot

ntupleBeamSpot = cms.EDProducer('NtupleBeamSpot',
    inputTag = cms.InputTag('simBeamSpot', 'BeamSpot'),
    prefix = cms.string('beamSpot@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleGenEventInfo = cms.EDProducer('NtupleGenEventInfo',
    genEventInfo = cms.InputTag('generator'),
    pileupInfo = cms.InputTag('addPileupInfo'),
    pileupWeight = cms.InputTag(''),
    pdfWeight = cms.InputTag(''),
    randomSeed = cms.InputTag('randomEngineSeedKeeper'),
    prefix = cms.string('gen@'),
    suffix = cms.string(''),
)

ntupleGen = cms.Sequence(ntupleGenEventInfo * ntupleGenParticles)
#ntupleGen = cms.Sequence(ntupleGenEventInfo * ntupleGenParticles * (simBeamSpot * ntupleBeamSpot))
#ntupleGen = cms.Sequence(ntupleGenEventInfo * ntupleGenParticles * ntupleGenJets * ntupleGenMET * (simBeamSpot * ntupleBeamSpot))

