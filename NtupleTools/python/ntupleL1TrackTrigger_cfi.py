import FWCore.ParameterSet.Config as cms

ntuplePixelDigis = cms.EDProducer('NtuplePixelDigis',
    inputTag = cms.InputTag('simSiPixelDigis'),
    inputTagTP = cms.InputTag('mix', 'MergedTrackTruth'),
    simHitCollections = cms.PSet(
        pixel = cms.VInputTag(
            cms.InputTag('mix','g4SimHitsTrackerHitsPixelBarrelLowTof'),
            cms.InputTag('mix','g4SimHitsTrackerHitsPixelBarrelHighTof'),
            cms.InputTag('mix','g4SimHitsTrackerHitsPixelEndcapLowTof'),
            cms.InputTag('mix','g4SimHitsTrackerHitsPixelEndcapHighTof'),
        ),
    ),
    prefix = cms.string('pixelDigis@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleTTClusters = cms.EDProducer('NtupleTTClusters',
    inputTag = cms.InputTag('TTClustersFromPhase2TrackerDigis', 'ClusterInclusive'),
    inputTagMC = cms.InputTag('TTClusterAssociatorFromPixelDigis', 'ClusterInclusive'),
    inputTagDigi = cms.InputTag('mix', 'Tracker'),
    inputTagTP = cms.InputTag('mix', 'MergedTrackTruth'),
    prefix = cms.string('TTClusters@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleTTStubs = cms.EDProducer('NtupleTTStubs',
    inputTag = cms.InputTag('TTStubsFromPhase2TrackerDigis', 'StubAccepted'),
    inputTagMC = cms.InputTag('TTStubAssociatorFromPixelDigis', 'StubAccepted'),
    inputTagClus = ntupleTTClusters.inputTag,
    inputTagDigi = ntupleTTClusters.inputTagDigi,
    inputTagTP = ntupleTTClusters.inputTagTP,
    prefix = cms.string('TTStubs@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleTTTracks = cms.EDProducer('NtupleTTTracks',
    inputTag = cms.InputTag('TTTracksFromPhase2TrackerDigis', 'Level1TTTracks'),
    inputTagMC = cms.InputTag('TTTrackAssociatorFromPhase2TrackerDigis', 'Level1TTTracks'),
    inputTagStub = ntupleTTStubs.inputTag,
    nparameters = cms.int32(4),
    prefix = cms.string('TrackletTTTracks@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleTTClustersForTTI = ntupleTTClusters.clone(
    inputTag = cms.InputTag('TTStubsFromPhase2TrackerDigis', 'ClusterAccepted'),
    inputTagMC = cms.InputTag('TTClusterAssociatorFromPhase2TrackerDigis', 'ClusterAccepted'),
    inputTagDigi = cms.InputTag(''),
)

ntupleTTStubsForTTI = ntupleTTStubs.clone(
    inputTagClus = ntupleTTClustersForTTI.inputTag,
    inputTagDigi = cms.InputTag(''),
)

#ntupleL1TrackTrigger = cms.Sequence(ntuplePixelDigis * ntupleTTClusters * ntupleTTStubs * ntupleTTTracks)
#ntupleL1TrackTrigger_TTI = cms.Sequence(ntupleTTClustersForTTI * ntupleTTStubsForTTI * ntupleTTTracks)

ntupleL1TrackTrigger = cms.Sequence(ntupleTTStubs)

