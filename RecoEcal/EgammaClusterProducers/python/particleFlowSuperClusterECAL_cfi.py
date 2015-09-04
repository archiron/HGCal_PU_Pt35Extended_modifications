import FWCore.ParameterSet.Config as cms
import os

particleFlowSuperClusterECALBox = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("Box"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(False),
#CC@@ added PUEthardCut boolean
    usePUEtHardCut = cms.bool(False),
    
    #use shower start cut?
#CC@@ this one is also changed by the stamndard SLHC customization
#    useHGCEmPreID = cms.bool(False),
    useHGCEmPreID = cms.bool(True),

    #PFClusters collection
#CC@@ this one is also changed by the stamndard SLHC customization
#    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    PFClusters = cms.InputTag('particleFlowClusterHGCEE'),    
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),    
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    #rechit collections for lazytools
    ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
    ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowBasicClusterECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowSuperClusterECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowBasicClusterECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowSuperClusterECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBasicClusterECALPreshower"),
#CC@@ this one is also changed by the stamndard SLHC customization
#    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowSuperClusterECALEndcapWithPreshower"),                                          
    PFSuperClusterCollectionEndcapWithPreshower = cms.string(''),

    #use preshower ?
#CC@@ this one is also changed by the stamndard SLHC customization
#    use_preshower = cms.bool(True),
    use_preshower = cms.bool(False),

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),

    # regression setup
    useRegression = cms.bool(False), #regressions are mustache only
    regressionKeyEB = cms.string('pfecalsc_EBCorrection'),
    regressionKeyEE = cms.string('pfecalsc_EECorrection'),
    
    # fixed matrix energy determination for HGCAL
    useFixedMatrixHGCAL = cms.bool(True), 
    
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),    
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(3.0),
    thresh_PFClusterBarrel = cms.double(0.5),

    thresh_PFClusterSeedEndcap = cms.double(5.0),
    thresh_PFClusterEndcap = cms.double(0.5),

    # window width in ECAL
    phiwidth_SuperClusterBarrel = cms.double(0.28),
    etawidth_SuperClusterBarrel = cms.double(0.04),

    phiwidth_SuperClusterEndcap = cms.double(0.28),
#CC@@ new eta width value for the box
#    etawidth_SuperClusterEndcap = cms.double(0.04),
    etawidth_SuperClusterEndcap = cms.double(0.015),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),                                          

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.),                                             

    #corrections
    applyCrackCorrections = cms.bool(False)                                          
                                              
)

particleFlowSuperClusterECALMustache = cms.EDProducer(
    "PFECALSuperClusterProducer",
    # verbosity 
    verbose = cms.untracked.bool(False),
    #clustering type: "Box" or "Mustache"
    ClusteringType = cms.string("Mustache"),
    #energy weighting: "Raw", "CalibratedNoPS", "CalibratedTotal"
    EnergyWeight = cms.string("Raw"),

    #this overrides both dphi cuts below if true!
    useDynamicDPhiWindow = cms.bool(True), 
    usePUEtHardCut = cms.bool(False),
                                
    #use shower start cut?
#CC@@ new value (also applied in the custom config)
#    useHGCEmPreID = cms.bool(False),
    useHGCEmPreID = cms.bool(True),
              
    #PFClusters collection
#CC@@ (also applied in the custom config)
#    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    PFClusters = cms.InputTag('particleFlowClusterHGCEE'),    
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
                                              
    PFBasicClusterCollectionBarrel = cms.string("particleFlowBasicClusterECALBarrel"),                                       
    PFSuperClusterCollectionBarrel = cms.string("particleFlowSuperClusterECALBarrel"),
    PFBasicClusterCollectionEndcap = cms.string("particleFlowBasicClusterECALEndcap"),                                       
    PFSuperClusterCollectionEndcap = cms.string("particleFlowSuperClusterECALEndcap"),
    PFBasicClusterCollectionPreshower = cms.string("particleFlowBasicClusterECALPreshower"),
#CC@@ (also applied in the custom config)
#    PFSuperClusterCollectionEndcapWithPreshower = cms.string("particleFlowSuperClusterECALEndcapWithPreshower"),                                          
    PFSuperClusterCollectionEndcapWithPreshower = cms.string(''),

    #use preshower ?
#CC@@ (also applied in the custom config)
#    use_preshower = cms.bool(True),
    use_preshower = cms.bool(False),

    # are the seed thresholds Et or Energy?
    seedThresholdIsET = cms.bool(True),
    # regression setup
    useRegression = cms.bool(False),
    regressionConfig = cms.PSet(
       regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v1'),
       regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v1'),
       vertexCollection = cms.InputTag("offlinePrimaryVertices"),
       ecalRecHitsEB = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
       ecalRecHitsEE = cms.InputTag('ecalRecHit','EcalRecHitsEE')
       ),
#CC@@ (new)
    useFixedMatrixHGCAL = cms.bool(True), 
       
    #threshold for final SuperCluster Et
    thresh_SCEt = cms.double(4.0),
    
    # threshold in ECAL
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterBarrel = cms.double(0.0),

    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_PFClusterEndcap = cms.double(0.0),

    # window width in ECAL ( these don't mean anything for Mustache )
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    etawidth_SuperClusterBarrel = cms.double(0.04),

    phiwidth_SuperClusterEndcap = cms.double(0.6),
    etawidth_SuperClusterEndcap = cms.double(0.04),

    # threshold in preshower
    thresh_PFClusterES = cms.double(0.),           

    # turn on merging of the seed cluster to its nearest neighbors
    # that share a rechit
    doSatelliteClusterMerge = cms.bool(False),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    #thresh_PFClusterMustacheOutBarrel = cms.double(0.),
    #thresh_PFClusterMustacheOutEndcap = cms.double(0.), 

    #corrections
    applyCrackCorrections = cms.bool(False)
                                              
)

particleFlowSuperClusterECAL = particleFlowSuperClusterECALMustache.clone()
#particleFlowSuperClusterECAL = particleFlowSuperClusterECALBox.clone()
