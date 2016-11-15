import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')
runOnMC = True

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.setDefault('inputFiles', ['file:relval_tilted.root']
options.setDefault('outputFile', 'ntuple.root')
options.parseArguments()


## MessageLogger
process.load('FWCore.MessageService.MessageLogger_cfi')

## Options
process.options = cms.untracked.PSet(

)

## Input Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

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

## Make the ntuple
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)
process.load('L1TTSimulations.NtupleTools.ntupleSequences_cff')

process.p = cms.Path(process.ntupleSequence)
process.schedule.append(process.p)


# Configure framework report and summary
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

