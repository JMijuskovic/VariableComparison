import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('VariableComparison'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
