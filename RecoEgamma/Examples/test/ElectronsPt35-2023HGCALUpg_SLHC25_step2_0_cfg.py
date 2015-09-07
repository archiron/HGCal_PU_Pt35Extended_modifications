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
#    input = cms.untracked.int32(max_number), # ne prend que max_number entrees
    input = cms.untracked.int32(-1) # on prend tout
)

# Input source
max_skipped = ind * max_number
process.source = cms.Source("PoolSource",
#    skipEvents = cms.untracked.uint32(max_skipped), # c'est ici
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/RecoEgamma/Examples/test/ElectronPt35-2023HGCALUpg_SLHC25-' + '%003d'%ind + '_step1.root')
)

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
    fileName = cms.untracked.string('file:/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/RecoEgamma/Examples/test/ElectronPt35-2023HGCALUpg_SLHC25-' + '%003d'%ind + '-RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(140.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0010AE1F-6676-E411-8F16-002618943860.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0035CDEE-5C76-E411-8214-0023AEFDEEEC.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/004B2C7D-6876-E411-ABFA-002618943949.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/006DDC01-6276-E411-9E66-00259073E4E4.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/008F6D89-5976-E411-A05E-549F35AC7DEE.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02133DBD-6176-E411-967A-002590A8882A.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/0253431B-4F76-E411-ABFE-0025904C66F4.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02758CA9-5F76-E411-A1D8-0015172C07E1.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02C7F040-7176-E411-B19E-0023AEFDEE68.root', '/store/mc/TP2023HGCALGS/MinBias_TuneZ2star_14TeV-pythia6/GEN-SIM/DES23_62_V1-v3/00000/02DE880D-5576-E411-AE26-002590200A00.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V4::All', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023HGCal 

#call to customisation function cust_2023HGCal imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023HGCal(process)

# End of customisation functions
