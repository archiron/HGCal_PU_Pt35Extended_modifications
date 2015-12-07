# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --filein file:ElectronsPt35-2023HGCALUpg_SLHC25_step1.root --fileout file:ElectronsPt35-2023HGCALUpg_SLHC25.root 
#                                  --pileup_input dbs:/MinBias_TuneZ2star_14TeV-pythia6/TP2023HGCALGS-DES23_62_V1-v3/GEN-SIM --mc --eventcontent RECOSIM 
#                                  --pileup AVE_140_BX_25ns --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCal --datatier GEN-SIM-RECO 
#                                  --conditions PH2_1K_FB_V4::All --step RAW2DIGI,L1Reco,RECO --magField 38T_PostLS1 --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco 
#                                  --python_filename ElectronsPt35-2023HGCALUpg_SLHC25_step2_cfg.py --no_exec -n 100
import FWCore.ParameterSet.Config as cms
import os, sys

if len(sys.argv) > 1:
    print "step 2 - arg. 0 :", sys.argv[0]
    print "step 2 - arg. 1 :", sys.argv[1]
    print "step 2 - arg. 2 :", sys.argv[2]
    ind = int(sys.argv[2])
else:
    print "step 2 - rien"
    ind = 0
max_number = 50 # number of events
    
process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(max_number), # ne prend que max_number entrees
#    input = cms.untracked.int32(-1) # on prend tout
)

# Input source
max_skipped = ind * max_number
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(max_skipped), # c'est ici
    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring('file:ElectronsPt35-2023HGCALUpg_SLHC25_step1.root'),
#)
    fileNames = cms.untracked.vstring([
#        '/store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/00435638-29E2-E411-97CE-0026189438D7.root',
#        '/store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/022EA57F-40E3-E411-A6A5-0025904C68DE.root'
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/8C341D05-86E2-E411-B00F-002481E0DC66.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/E8AFE63C-DBE1-E411-BB0B-002590AC5082.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/AAF06B6B-C0E1-E411-9C6E-002590A2CD06.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/54D6CC80-C4E1-E411-9181-0025904B7C26.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/00399685-C9E1-E411-9311-0025905B85D8.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/008EB715-18E2-E411-A7AD-00266CFFBC64.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/8C341D05-86E2-E411-B00F-002481E0DC66.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/E8AFE63C-DBE1-E411-BB0B-002590AC5082.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/00B2B55C-8BE2-E411-B14D-003048947BB9.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0230E3F0-CCE1-E411-A4D6-0025904C637C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0406FD1F-A6E2-E411-B6B7-6CC2173BBA30.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/04191A4E-ADE1-E411-BBA8-0025907B4F78.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/045CA4E5-BBE1-E411-93A4-0025905B861C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/060C20B7-A3E1-E411-A051-0025905B8596.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/061754CD-B7E1-E411-8814-E0CB4E55364D.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/062DE76F-76E1-E411-8131-3417EBE2F3FA.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0647A551-BCE1-E411-93E3-20CF300E9EAF.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0689F94C-07E2-E411-AEE6-00266CFFBE5C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/08C1265D-E0E1-E411-B7D5-0025904C67B6.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/08F22FAA-72E1-E411-AFD3-20CF305B04F5.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/08FABA9F-3EE2-E411-8D8E-0025904CF93C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0C1C0DD1-E0E1-E411-879D-0025904C68D8.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0C79B2BE-9BE1-E411-AD01-008CFA051DA8.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0C97714F-92E1-E411-9CF1-001E673970C1.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0CA83C76-76E1-E411-9548-848F69FD284D.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0CFF6810-CCE1-E411-8141-C4346BC70B58.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/0E561972-06E2-E411-8B0C-002590A4FFE8.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/1040DDD5-FCE1-E411-8146-00259020083C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/127CB177-8AE1-E411-BE9E-0025904C66E4.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/12CEA9AC-B5E1-E411-A017-00266CFFCAC8.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/141F0C88-C5E1-E411-908D-002590A371AE.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/1601B21D-A1E1-E411-9105-0025904CF712.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/165DC246-01E2-E411-9225-C4346BBCB408.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/16DF8B15-0FE2-E411-9A77-6CC2173BC2E0.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/16EEA82D-8FE1-E411-87E4-0025905A6088.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/18004E93-C5E1-E411-8136-00266CFFBE14.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/185E2C90-91E1-E411-A61E-00266CFFC76C.root',
        ' root://llrxrd-redir.in2p3.fr//store/mc/TP2023HGCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-DIGI-RAW/HGCALMar26_PU140BX25_PH2_1K_FB_V6-v2/10000/189C6221-6EE1-E411-8E7B-002590D0B01E.root',
    ])
)
print "load OK"

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step2 nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('file:ElectronsPt35-2023HGCALUpg_SLHC25.root'),
    fileName = cms.untracked.string('file:/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/RecoEgamma/Examples/test/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola-2023HGCALUpg_SLHC25-' + '%003d'%ind + '-RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements ## non used in step2
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring([
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0010AE1F-6676-E411-8F16-002618943860.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0035CDEE-5C76-E411-8214-0023AEFDEEEC.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/004B2C7D-6876-E411-ABFA-002618943949.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/006DDC01-6276-E411-9E66-00259073E4E4.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/008F6D89-5976-E411-A05E-549F35AC7DEE.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02133DBD-6176-E411-967A-002590A8882A.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0253431B-4F76-E411-ABFE-0025904C66F4.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02758CA9-5F76-E411-A1D8-0015172C07E1.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02C7F040-7176-E411-B19E-0023AEFDEE68.root', 
    '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02DE880D-5576-E411-AE26-002590200A00.root'
    ])

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')
print "GlobalTag OK"

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
print "process schedule OK"

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCal 

#call to customisation function cust_2023HGCal imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCal(process)

# End of customisation functions
