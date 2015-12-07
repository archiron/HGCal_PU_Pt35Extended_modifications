# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein dbs:/RelValSingleElectronPt35Extended/CMSSW_6_2_0_SLHC24-PH2_1K_FB_V6_HGCalPandPU140-v1/GEN-SIM 
#                                  --fileout file:ElectronsPt35-2023HGCALUpg_SLHC25_step1.root 
#                                  --pileup_input dbs:/MinBias_TuneZ2star_14TeV-pythia6/TP2023HGCALGS-DES23_62_V1-v3/GEN-SIM --mc --eventcontent FEVTDEBUG 
#                                  --pileup AVE_140_BX_25ns --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCal 
#                                  --datatier GEN-SIM-DIGI-RAW --conditions PH2_1K_FB_V4::All --step DIGI:pdigi_valid,L1,DIGI2RAW --magField 38T_PostLS1 
#                                  --geometry Extended2023HGCalMuon,Extended2023HGCalMuonReco --python_filename ElectronsPt35-2023HGCALUpg_SLHC25_step1_cfg.py --no_exec -n 100
import FWCore.ParameterSet.Config as cms
import os, sys

if len(sys.argv) > 1:
    print "step 1 - arg. 0 :", sys.argv[0]
    print "step 1 - arg. 1 :", sys.argv[1]
    print "step 1 - arg. 2 :", sys.argv[2]
    ind = int(sys.argv[2])
else:
    print "step 1 - rien"
    ind = 0
max_number = 50 # number of events
    

process = cms.Process('DIGI2RAW')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(max_number), # ne prend que max_number entrees
)

# Input source
max_skipped = ind * max_number
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(max_skipped), # c'est ici
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
#        '/store/relval/CMSSW_6_2_0_SLHC24/RelValSingleElectronPt35Extended/GEN-SIM/PH2_1K_FB_V6_HGCalPandPU140-v1/00000/6E7DA459-52C4-E411-8AD3-0025905A48BC.root', 
#        '/store/relval/CMSSW_6_2_0_SLHC24/RelValSingleElectronPt35Extended/GEN-SIM/PH2_1K_FB_V6_HGCalPandPU140-v1/00000/CC22E259-52C4-E411-ACA4-0025905A6094.root'
#SIM
#        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/B64C033A-1820-E511-A27F-0025905964BC.root',
#        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/CA152C39-1820-E511-9469-0025905B8590.root',
#SIM-RECO
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/0E832677-8320-E511-9416-002618943916.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/14781706-7220-E511-BB30-0025905964A6.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/14CC1F94-7320-E511-A0EA-002590593876.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/1688D91F-7A20-E511-8278-003048D15DF0.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/306DEFC9-7420-E511-873A-0025905A48E4.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/32D9A603-7F20-E511-B971-0025905B8582.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/3CB8D7C6-7420-E511-8E35-0025905A6068.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/4AA02459-7220-E511-BE69-0026189437F8.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/58DDBDB1-7620-E511-A4AE-002618943916.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/607B198D-7120-E511-B486-00261894398B.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/643A7B2A-7D20-E511-BA8E-002354EF3BDF.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/645116CA-7120-E511-B455-00261894396B.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/6635CB7D-6E20-E511-B043-003048FFD76E.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/6C3A93E9-7A20-E511-A417-00261894391C.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/6EC35872-6F20-E511-943F-0026189438E8.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/76D0CA41-7320-E511-B490-00261894396B.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/924A2F20-7720-E511-A417-002618943916.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/927471BE-7120-E511-B5AC-002618943882.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/92E3862B-8A20-E511-A24C-0025905AA9CC.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/949DEF23-7720-E511-A845-003048D15DE4.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/A8EE3ED3-8120-E511-B678-0025905A609E.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/AAD209FF-7720-E511-8EF3-002618943896.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/B25EB269-8020-E511-8731-002354EF3BDC.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/B615F3D9-7820-E511-A262-0025905A608E.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/D00C0E9A-7320-E511-A28D-0025905964A6.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/D6537A57-7020-E511-972A-0025905A60DE.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DAAB0C68-8020-E511-A97E-0026189437F8.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DE4EF83E-7620-E511-AF0E-002618943916.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/E0DF4D26-7E20-E511-969B-0025905938A8.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/E6263331-7D20-E511-8C8C-002590593876.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/E650A3A2-7B20-E511-B336-0025905A60A0.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/F2B1D2D5-7920-E511-9348-0025905A608E.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/F2CCA736-7420-E511-A90F-0025905A6122.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/FA9FDEA1-7E20-E511-9AC5-0025905B8592.root',
        '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-RECO/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/FE231D03-7F20-E511-AB34-0025905A6066.root',
#RAW
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/0002BB08-5B20-E511-9ADA-003048FFCB74.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/0653DC79-5820-E511-8AF9-0025905A6118.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/0832D77C-5B20-E511-85FF-003048FFCC18.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/14E369A4-3420-E511-84ED-0025905964BA.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/185BAC47-5A20-E511-9915-0025905A60B6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/1A89DF69-5C20-E511-ADBF-0025905A48BC.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/1AF0D033-5820-E511-8755-0025905B85AE.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/22008578-5B20-E511-89FB-0025905A60B6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/24457127-5620-E511-A386-0025905A6068.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/28A3C2E9-5820-E511-A2E2-0025905A60D6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/2E342D09-5B20-E511-A541-0025905A60BC.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/2EB95267-3420-E511-BF4F-0025905A608A.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/32D6D90B-5B20-E511-8172-0025905A60CA.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/3679DD1C-5920-E511-8097-0025905B85AE.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/3AE0E9D6-5920-E511-91AE-0025905A612E.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/3E6816E8-5820-E511-910F-002618943930.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/42B2ACC1-5920-E511-8AB4-0025905AA9CC.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/44ED2DD8-5720-E511-8E0F-003048FF9AC6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/589DAB5F-3420-E511-8403-0025905A48E4.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/58DB0285-5620-E511-B0DB-00248C0BE014.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/5E736F7D-5820-E511-8465-0026189438EA.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/5E7FF178-5B20-E511-AD02-0025905A48BC.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/62B541DC-5920-E511-879F-0025905964A6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/66D57FE8-5920-E511-9900-0025905B855C.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/68D35AD8-5720-E511-8923-0025905B85AE.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/6A8DDE0D-5B20-E511-BC87-0025905A60A8.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/76D080D4-5720-E511-A947-002618943930.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/7CADC288-5620-E511-A0F8-0025905A612E.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/7E629633-5820-E511-A212-003048FF9AC6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/82A82BDF-5920-E511-AB3E-003048FF9AC6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/84593BDE-5620-E511-B307-00248C0BE014.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/8A12BA74-5A20-E511-991A-00248C0BE014.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/8CA81EE8-5820-E511-B76E-002618943854.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/8CE14DE8-5920-E511-97C5-0025905A60BC.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/96BD4348-5A20-E511-981B-0025905B855C.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/98699838-5820-E511-847F-002590593872.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/98E91137-5620-E511-ADDF-0025905A60D6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/9C070FD5-5720-E511-B65B-0026189438F6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/A00E08EE-5820-E511-9914-0025905B85EE.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/A4A57E71-5A20-E511-857E-0026189438B9.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/AA945D7F-5A20-E511-9A10-0025905A6070.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/AEBCDA72-3620-E511-A63B-002354EF3BE4.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/B408D907-5B20-E511-8F4C-0025905A6122.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/B4B37877-5A20-E511-A2E2-0025905B859E.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/BE4A3D44-5520-E511-89A5-003048FFCC18.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/C2B4E8DD-5620-E511-A17C-0026189438DF.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/C6D27809-5B20-E511-97E6-0025905B8598.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/CC8B4448-5A20-E511-9CB7-0025905A6118.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/CED74A91-5820-E511-A023-0025905A48D8.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/D4A57561-3420-E511-99F4-0025905B85F6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/D67035C3-5920-E511-82AF-0025905A60B4.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DA25F9C3-5920-E511-8A37-0025905A6118.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DC81A206-5B20-E511-BA68-00248C0BE014.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DCA7A17C-5B20-E511-854F-0025905A60A8.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/DE4965DF-5920-E511-9E96-0025905964C2.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/E4B4E61F-5920-E511-97B4-003048FFCB74.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/E644C669-5C20-E511-94BD-0025905B85D8.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/EC7303D9-5720-E511-9F29-002590593872.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/ECA7E1E2-5920-E511-982A-0025905A60D6.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/EE86D7D6-5720-E511-B965-0025905A608E.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/F64BDE09-5B20-E511-AB57-0025905A611C.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/F8DF2586-5620-E511-B2FC-0026189438DF.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/FC35121A-5E20-E511-A03A-003048FFCC18.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/FCB0C150-3920-E511-B529-0026189438EB.root',
#       '/store/relval/CMSSW_6_2_0_SLHC26_patch2/RelValSingleElectronPt35Extended/GEN-SIM-DIGI-RAW/PH2_1K_FB_V6_HGCalee28PU140-v1/00000/FE98A088-5620-E511-B459-0025905A612E.root',
        )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step1 nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/RecoEgamma/Examples/test/ElectronPt35-2023HGCALUpg_SLHC25-' + '%003d'%ind + '_step1.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0010AE1F-6676-E411-8F16-002618943860.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0035CDEE-5C76-E411-8214-0023AEFDEEEC.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/004B2C7D-6876-E411-ABFA-002618943949.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/006DDC01-6276-E411-9E66-00259073E4E4.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/008F6D89-5976-E411-A05E-549F35AC7DEE.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02133DBD-6176-E411-967A-002590A8882A.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0253431B-4F76-E411-ABFE-0025904C66F4.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02758CA9-5F76-E411-A1D8-0015172C07E1.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02C7F040-7176-E411-B19E-0023AEFDEE68.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02DE880D-5576-E411-AE26-002590200A00.root'])
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.endjob_step,process.FEVTDEBUGoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCal 

#call to customisation function cust_2023HGCal imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCal(process)

# End of customisation functions
