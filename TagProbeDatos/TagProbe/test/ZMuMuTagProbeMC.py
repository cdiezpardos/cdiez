import FWCore.ParameterSet.Config as cms

process = cms.Process("MyProcess")

process.load("Configuration.EventContent.EventContent_cff")





process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('Quality'),
    categories = cms.untracked.vstring('FwkJob','FwkReport','FwkSummary','Root_NoDictionary'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
#        threshold = cms.untracked.string('INFO')
       #threshold = cms.untracked.string('DEBUG')
    ),
    destinations = cms.untracked.vstring('cout')
#    destinations = cms.untracked.vstring('TagProbeSummary'),
#    TagProbeSummary=cms.untracked.PSet(threshold = cms.untracked.string('INFO'))
)


process.AdaptorConfig = cms.Service("AdaptorConfig")




process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50000)
)


import os

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
#"rfio:/castor/cern.ch/user/g/goys/Wprime/3_8_X/filter/DY_M20_Fall10/Skim_10_1_Ek9.root"
     ) ,
     inputCommands = cms.untracked.vstring(
     'keep *', 'drop *_lumiProducer_*_*', 'drop *_MEtoEDMConverter_*_*', 'drop *_l1GtTriggerMenuLite_*_*' , 'drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT*','drop edmTriggerResults_TriggerResults__makeSD','drop *_*__JuntandoSkims','drop *_*__SuperMuSkim', 'drop edmTriggerResults_TriggerResults__EWKHighPtMuSkim'
     ) 
   
)

file_directory = "/ciet3b/data3/Fall10_All_MinimalAOD/ZmumuPOWHEG/"
for file in os.listdir(file_directory):
      process.source.fileNames.append("file://"+file_directory + "/" + file)
#      process.source.fileNames.append("dcap://dchadm.ciemat.es:22125"+file_directory + "/" + file)



process.Quality = cms.EDAnalyzer("ZMuMuTagProbe",
#-------------
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD"),
    hltL3FilterLabel = cms.InputTag("hltSingleMu21L3Filtered25::REDIGI311X"),#CHANGE TO ADAPT TO THE CHOSEN TRIGGER 
    hltL3FilterLabelMu11 = cms.InputTag("hltSingleMu11L3Filtered11::HLT"),
    hltL3FilterLabelMu15 = cms.InputTag("hltSingleMu15L3Filtered15::HLT"),
#------------
    # Main cuts -> 
    PtCut = cms.untracked.double(25.0),  
    EtaCut = cms.untracked.double(2.1),
    IsTotIso = cms.untracked.bool(True),    #----> This allows the change between track-only ("false") & combined ("true") isolation!
    IsoCut03 = cms.untracked.double(0.1),
    doReweight = cms.untracked.bool(False), 

#Sample 7TeV------------------------------
    
    UseMuonQualityCuts=cms.untracked.bool(True),

    SLOPE = cms.untracked.double(0.005),
    ABC = cms.untracked.double(0.66),
    ZPAR = cms.untracked.double(0.920),
    SLOPE2 = cms.untracked.double(-0.087),
    ABC2 = cms.untracked.double(11.6),
    ZPAR2 = cms.untracked.double(0.927)
)

process.Tracking = cms.EDAnalyzer("ZMuMuTracking",
#-------------
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD"),
#    hltL3FilterLabel = cms.InputTag("hltSingleMuNoIsoL3PreFiltered15::HLT"), 
    hltL3FilterLabel = cms.InputTag("hltSingleMu11L3Filtered11::REDIGI38XPU"),#CHANGE TO ADAPT TO THE CHOSEN TRIGGER
    hltL3FilterLabelMu11 = cms.InputTag("hltSingleMu11L3Filtered11::HLT"),
    hltL3FilterLabelMu15 = cms.InputTag("hltSingleMu15L3Filtered15::HLT"),
#------------
    # Main cuts -> 
    PtCut = cms.untracked.double(25.0),
    EtaCut = cms.untracked.double(2.1),
    IsTotIso = cms.untracked.bool(True),    #----> This allows the change between track-only ("false") & combined ("true") isolation!
    IsoCut03 = cms.untracked.double(0.1),
    UseMuonQualityCuts=cms.untracked.bool(True),

)




process.TFileService = cms.Service("TFileService", fileName = cms.string('MC.root'))

process.pReco = cms.Path(process.Quality)
process.pTrack = cms.Path(process.Tracking)


