import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C1)

process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D1_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.a1 = cms.EDAnalyzer("AnalyzerModuleVertices",
    csv = cms.string('module_vertices.csv'),
    verbosity = cms.int32(0),
)

process.p1 = cms.Path(process.a1)
process.schedule = cms.Schedule(process.p1)

# customisation of the process
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted
process = cust_2023tilted(process)

