import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
runOnMC = True

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()
# Modify the defaults
if not options.inputFiles:
    options.inputFiles = ['file:relval_tilted.root']
if options.outputFile == "output.root":
    options.outputFile = "test_ntuple.root"


## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

## Options and Output Report
process.options = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool( True ),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Geometry and Global Tags
process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D1_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

## Track trigger
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.L1TrackTrigger_step = cms.Path(process.TrackTriggerClustersStubs)
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.L1TTAssociator_step = cms.Path(process.TrackTriggerAssociatorClustersStubs)

## Schedule definition
process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1TTAssociator_step)

## HL-LHC upgrade customization
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

process = cust_2023tilted(process)

## Write the TTree
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.load('L1TTSimulations.NtupleTools.ntupleSequences_cff')
process.p = cms.Path(process.ntupleSequence)
process.schedule.append(process.p)

