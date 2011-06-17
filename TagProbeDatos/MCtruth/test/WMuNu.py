import FWCore.ParameterSet.Config as cms

process = cms.Process("MyProcess")

process.load("Configuration.EventContent.EventContent_cff")


process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('Quality','NoQuality'),
    categories = cms.untracked.vstring('FwkJob','FwkReport','FwkSummary','Root_NoDictionary'),
    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('ERROR')
        threshold = cms.untracked.string('INFO')
#       threshold = cms.untracked.string('DEBUG')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.AdaptorConfig = cms.Service("AdaptorConfig")




process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(False)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)

)
# For WMuNu -> 10pb-1 NLO stats ~= 97000 events  (LO ~= 88000)
# (Check in DBS for details)
# (1 pb-1 --> 9700)

#process.source = cms.Source("PoolSource",
#    debugVerbosity = cms.untracked.uint32(0),
#    debugFlag = cms.untracked.bool(True),
#    fileNames = cms.untracked.vstring("")
#)

process.source = cms.Source("PoolSource",

fileNames = cms.untracked.vstring(
# 'rfio:/castor/cern.ch/user/e/ecarrera/analyses/wprime_munu/wprime/Fall10MC/RECO/wprime_RECO_1000GeV_3_3_jXJ.root',
# 'rfio:/castor/cern.ch/user/e/ecarrera/analyses/wprime_munu/wprime/Fall10MC/RECO/wprime_RECO_1000GeV_4_6_7jq.root'

# 'rfio:/castor/cern.ch/user/e/ecarrera/analyses/wprime_munu/wprime/Fall10MC/RECO/wprime_RECO_1100GeV_1_2_bye.root',
# 'rfio:/castor/cern.ch/user/e/ecarrera/analyses/wprime_munu/wprime/Fall10MC/RECO/wprime_RECO_1100GeV_2_1_W9o.root'


)#'file:/data1/cdiez/EWK_SubSkim_Summer09_10TeV/ZMuMu_10TeV_10invpb_1.root')

)

import os
#dirname = "/ciet3b/data4/Summer10_10invpb_AODSIM/ZmumuPOWHEG/"
dirname = "/data3/Fall10_39X_Pileup_MinimalAOD/ZmumuPOWHEG/"
#dirname = "/pnfs/ciemat.es/data/cms/store/user/jalcaraz/Summer08_Zmumu_M20_IDEAL_V11_redigi_v1_AODRED/"
#basenamelist = os.listdir(dirname)
#for basename in basenamelist:
#    process.source.fileNames.append("file:" + dirname + basename) 
#print "Number of files to process is %s" % (len(process.source.fileNames))

process.NoQuality =  cms.EDAnalyzer("WMuNu",

 # Defaults follow...
    #
    # Input collections ->
    #MuonTag = cms.InputTag("muons"),
    METTag = cms.untracked.InputTag("corMetGlobalMuons"), # --> Check this line and the following one are in agreement!!
    MuonMETCor = cms.untracked.bool(True),  # --> For muon corrections (True means you do not need to correct from muonpt in the analyzer)
    UseMuonQualityCuts=cms.untracked.bool(False), #--> This applies quality cuts (from MuonID note) to the sample :-)
    #JetTag = cms.untracked.InputTag("iterativeCone5CaloJets"),
    #UseOnlyGlobalMuons = cms.untracked.bool(True), # --> By default, we only use global muons. This line is kinda obsolete after the quality cuts.. 
    #
    # Main cuts -> 
    #PtCut = cms.untracked.double(25.0),  
    #EtaCut = cms.untracked.double(2.0),
    #IsTotIso = cms.untracked.bool(False),    #----> This allows the change between track-only ("false") & combined ("true") isolation!
    #IsoCut03 = cms.untracked.double(0.1),
    #MassTMin = cms.untracked.double(50.0),
    #MassTMax = cms.untracked.double(200.0),
    #
    # To suppress Zmm ->
    #PtThrForZCount = cms.untracked.double(20.0),
    #
    # To suppress ttbar -> THIS CUTS HAVE BEEN REMOVED!!!
    #AcopCut = cms.untracked.double(99999.),
    #EJetMin = cms.untracked.double(10.),
    #NJetMax = cms.untracked.int32(999)

)

process.Quality = cms.EDAnalyzer("WMuNu",
    METTag = cms.untracked.InputTag("pfMet"),
    MuonMETCor = cms.untracked.bool(True),
    UseMuonQualityCuts=cms.untracked.bool(True),
    IsTotIso = cms.untracked.bool(True),
    IsoCut03 = cms.untracked.double(0.1),
#    PtCut = cms.untracked.double(25.0),  
#    EtaCut = cms.untracked.double(2.1),
#-------------
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD"),
#    hltL3FilterLabel = cms.InputTag("hltSingleMuNoIsoL3PreFiltered15::HLT"),
    hltL3FilterLabel = cms.InputTag("hltSingleMu11L3Filtered11::HLT"),
#------------


)


process.NoQuality_TotIso = cms.EDAnalyzer("WMuNu",
    METTag = cms.untracked.InputTag("pfMet"),
    MuonMETCor = cms.untracked.bool(True),
    UseMuonQualityCuts=cms.untracked.bool(False),
    IsTotIso = cms.untracked.bool(True),
)

process.Quality_TotIso = cms.EDAnalyzer("WMuNu",
    METTag = cms.untracked.InputTag("pfMet"),
    MuonMETCor = cms.untracked.bool(True),
    UseMuonQualityCuts=cms.untracked.bool(True),
    IsTotIso = cms.untracked.bool(True),
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('ZMuMu_hiMassGENycutPrueba.root'))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('ZmumuNotrigger.root'))

process.TFileService = cms.Service("TFileService", fileName = cms.string('WprimePRUEBA15noloop.root'))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('ZMMpt22_etaphi25_7TeV.root'))

#process.pReco = cms.Path(process.Quality+process.NoQuality+process.NoQuality_TotIso+process.Quality_TotIso)
process.pReco = cms.Path(process.Quality)


