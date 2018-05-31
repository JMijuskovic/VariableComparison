import FWCore.ParameterSet.Config as cms

process = cms.Process("VariableComparison")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("DQMServices.Core.DQM_cfg")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	 'file:/eos/cms/store/data/Run2018A/SingleMuon/USER/MuonPOGSkim-PromptReco-v2/000/316/877/00000/503AC7F6-DE61-E811-91DD-FA163E336D90.root',
         'file:/eos/cms/store/data/Run2018A/SingleMuon/USER/MuonPOGSkim-PromptReco-v2/000/316/877/00000/706491B8-F461-E811-81D5-FA163E510BE5.root',  
         'file:/eos/cms/store/data/Run2018A/SingleMuon/USER/MuonPOGSkim-PromptReco-v2/000/316/877/00000/7A9B8443-DB61-E811-BDC6-02163E01A083.root'
   )
)



process.VariableComparison = cms.EDAnalyzer('VariableComparison', 
       MuonCollection = cms.InputTag("slimmedMuons"),
       VertexLabel   = cms.InputTag("offlineSlimmedPrimaryVertices"),
       BeamSpotLabel = cms.InputTag("offlineBeamSpot"), 
       triggerObjects  = cms.InputTag("slimmedPatTrigger"),
       triggerResults = cms.InputTag("TriggerResults","","HLT"),

       pair_minInvMass  = cms.double(81),
       pair_maxInvMass  = cms.double(101),
       hlt_path=cms.string("HLT_IsoMu27_v16"),
       tag_hltFilter=cms.string("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"),
       tag_hltDrCut=cms.double(0.15),
       tag_minPt=cms.double(29),
       tag_isoCut=cms.double(0.15),
       tag_muonID = cms.string("TIGHT"),
       tagMuonID_selection = cms.string("pt>20"),
       probe_minPt=cms.double(20),
       probe_trackerIsoCut=cms.double(0.1),
       probe_muonIDs = cms.string("TIGHT"),
       noTrigger = cms.bool(False),
       etaMin = cms.vdouble(-2.4,-0.9,-2.1,+1.2,-2.4,+2.1),
       etaMax = cms.vdouble(+2.4,+0.9,-1.2,+2.1,-2.1,+2.4),
)

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('Plots.root')
)



process.p = cms.Path(process.VariableComparison)
