#ifndef HGCALElectronClusterAnalyzer_h
#define HGCALElectronClusterAnalyzer_h

//
// Package:         RecoEgamma/Examples
// Class:           HGCALElectronClusterAnalyzer
//

//
// Original Author:  Ursula Berthon, Claude Charlot
//         Created:  Mon Mar 27 13:22:06 CEST 2006
// $Id: HGCALElectronClusterAnalyzer.h,v 1.21 2011/05/20 17:17:28 wmtan Exp $
//
//


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"


#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h" 
#include "FastSimulation/Utilities/interface/GammaFunctionGenerator.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h" 
#include "FastSimulation/ShowerDevelopment/interface/EMECALShowerParametrization.h"
#include "FastSimulation/ShowerDevelopment/interface/EMShower.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include <vector>

class MagneticField;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TH1I;
class TProfile;
class TTree;

class HGCALElectronClusterAnalyzer : public edm::EDAnalyzer
{
 public:

  explicit HGCALElectronClusterAnalyzer(const edm::ParameterSet& conf);

  virtual ~HGCALElectronClusterAnalyzer();

  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  void doTest(const HGCalGeometry&, ForwardSubdetector);
  GlobalVector rotateMomentum(const MagneticField *, GlobalVector, GlobalPoint, GlobalPoint, TrackCharge);

 private:

  void readParameters(const edm::ParameterSet&);   
  
  TrajectoryStateTransform transformer_;
  edm::ESHandle<TrackerGeometry> pDD;
  
  TFile *histfile_;
  TTree *tree_;
  float mcEnergy[10], mcEta[10], mcPhi[10], mcPt[10], mcQ[10];
  float superclusterEnergy[10], superclusterEta[10], superclusterPhi[10], superclusterEt[10];
  float seedMomentum[10], seedEta[10], seedPhi[10], seedPt[10], seedQ[10];

  TH1F *h_mcNum;
  TH1F *h_simEta;
  TH1F *h_simAbsEta;
  TH1F *h_simP;
  TH1F *h_simPt;
  TH1F *h_simPhi;
  TH1F *h_simZ;
  TH2F *h_simPtEta;

  TH1F *h_simZ_all;
  TH1F *h_simZ_electrons;
   
  // hgcal histos
  TH1F *h_hgcal_foundClusters_em;
  TH1F *h_hgcal_foundClustersVSeta_em;
  TH1F *h_hgcal_foundClustersVSetaEt5_em;
  TH1F *h_hgcal_clusterEnergy_em;
  TH2F *h_hgcal_clusterEtaVsPhi_em; 
  TH1F *h_hgcal_foundClusters_had;
  TH1F *h_hgcal_clusterEnergy_had;
  TH2F *h_hgcal_clusterEtaVsPhi_had; 
  int nevt;
  int ievent;

  std::vector<TH3F *> h_hgcal_shower_seed_em;
  std::vector<TH3F *> h_hgcal_shower_sc_em;
  std::vector<TH3F *> h_hgcal_shower_sc;
  std::vector<TH3F *> h_hgcal_shower_seed_rotated_em;
  std::vector<TH3F *> h_hgcal_allhits_em;
  std::vector<TH3F *> h_hgcal_allhits_weighted_em;
  std::vector<TH3F *> h_hgcal_clustershits_em;
  std::vector<TH3F *> h_hgcal_clustershits_weighted_em;
  std::vector<TH3F *> h_hgcal_clustershits_cutlayers_em;
  std::vector<TH3F *> h_hgcal_clustershits_cutlayers_weighted_em;
  std::vector<TH3F *> h_hgcal_clustershits_cutlayers_cutlength_em;
  std::vector<TH3F *> h_hgcal_clustershits_cutlayers_cutlength_weighted_em;
  std::vector<TH3F *> h_hgcal_superclustershits_em;
  std::vector<TH3F *> h_hgcal_superclustershits_weighted_em;
  std::vector<TH3F *> h_hgcal_superclustershits_cutlayers_em;
  std::vector<TH3F *> h_hgcal_superclustershits_cutlayers_weighted_em;
  std::vector<TH3F *> h_hgcal_superclustershits_cutlayers_cutlength_em;
  std::vector<TH3F *> h_hgcal_superclustershits_cutlayers_cutlength_weighted_em;

  TH1F *h_hgcal_sclusters_energy_em;
  TH1F *h_hgcal_sclusters_transverse_energy_em;
  TH1F *h_hgcal_sclusters_energy_pos_em;
  TH1F *h_hgcal_sclusters_energy_neg_em;
  TH2F *h_hgcal_sclusters_energyVSeta_em;
  TH2F *h_hgcal_sclusters_position_em;
  TH1F *h_hgcal_sclusters_etawidth_em;
  TH1F *h_hgcal_sclusters_phiwidth_em;
  TH2F *h_hgcal_sclusters_multiplicityVSeta_em;
  TH1F *h_hgcal_sclusters_newmultiplicity_em;
  TH2F *h_hgcal_sclusters_newmultiplicityVSeta_em;
  TH1F *h_hgcal_sclusters_multiplicity_em;
  TH1F *h_hgcal_sclusters_seedenergy_em;
  TH1F *h_hgcal_sclusters_newseedenergy_em;
  TH1F *h_hgcal_sclusters_seedfractions_em;
  TH2F *h_hgcal_sclusters_seedfractionVSeta_em;
  TH1F *h_hgcal_clusters_energy_em;
  TH2F *h_hgcal_clusters_position_em;
  TH1F *h_hgcal_clusters_multiplicity_em;
  TH2F *h_hgcal_clusters_multiplicityVSeta_em;
  TH1F *h_hgcal_clusters_rechitenergy_em;
  TH1F *h_hgcal_clusters_rechitenergy_12000_em;
  TH1F *h_hgcal_scclusters_ptoverpttrue_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_em ;
  TH2F *h_hgcal_scclusters_eoveretrueVSeta_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut00_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut04_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut1_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut2_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut4_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut10_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_cut20_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut4_em ;
  TH1F *h_hgcal_scclusters_eoveretrue_golden_em ;
  TH2F *h_hgcal_scclusters_detadphisubclusters_em;
  TH2F *h_hgcal_scclusters_detadphisubclusters_zoom_em;
  TH2F *h_hgcal_scclusters_detadphisubclusters_zoom_etagt2_em;
  TH2F *h_hgcal_scclusters_detadphisubclusters_zoom_etalt2_em;
  TH2F *h_hgcal_scclusters_detadphisubclusters_weighted_em;
  TH2F *h_hgcal_scclusters_detadphisubclusters_zoom_weighted_em;
  
  TH1F *h_hgcal_sclusters_etaPCAMinusEtaTrue;
  TH1F *h_hgcal_sclusters_phiPCAMinusPhiTrue;
  TH1F *h_hgcal_sclusters_etaSCMinusEtaTrue;
  TH1F *h_hgcal_sclusters_phiSCMinusPhiTrue;
  TH1F *h_hgcal_sclusters_etaPCAMinusEtaTrue_corr;
  TH1F *h_hgcal_sclusters_phiPCAMinusPhiTrue_corr;
  TH1F *h_hgcal_sclusters_etaSCMinusEtaTrue_corr;
  TH1F *h_hgcal_sclusters_phiSCMinusPhiTrue_corr;
  TH1F *h_hgcal_sclusters_etaPCAMinusEtaSC;
  TH1F *h_hgcal_sclusters_phiPCAMinusPhiSC;
  
  TH2F *h_hgcal_sclusters_etaSCMinusEtaTrueVsEta;
  TH2F *h_hgcal_sclusters_phiSCMinusPhiTrueVsPhi;  
  TH2F *h_hgcal_sclusters_etaSeedMinusEtaTrueVsEta;
  TH2F *h_hgcal_sclusters_etaPCAMinusEtaTrueVsEta;  
  TH2F *h_hgcal_sclusters_phiSCMinusPhiTrueVsEta;
  TH2F *h_hgcal_sclusters_phiPCAMinusPhiTrueVsEta;
  TH2F *h_hgcal_sclusters_phiPCAMinusPhiTrueVsPhi;
  
  TH1F *h_hgcal_sclusters_etaPCAMinusEtaTrue_golden;
  TH1F *h_hgcal_sclusters_phiPCAMinusPhiTrue_golden;
  TH1F *h_hgcal_sclusters_etaSCMinusEtaTrue_golden;
  TH1F *h_hgcal_sclusters_phiSCMinusPhiTrue_golden;
  TH1F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrue_golden;
  TH1F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrue_golden;
  TH1F *h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue_golden;

  TH1F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrue;
  TH1F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrue;
  TH1F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueUnsigned;
  TH1F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueNeg;
  TH1F *h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirVsEnergy;
  TH2F *h_hgcal_sclusters_thetadirPCAMinusThetaDirTrueVsEnergy;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEnergy;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsPhi;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaMC;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaPCA;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicity;
  TH2F *h_hgcal_sclusters_AbsetadirPCAMinusEtaDirTrueVsSeedmultiplicity;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedfraction;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicityPosZ;
  TH2F *h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSCmultiplicity;  
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEta;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiMC;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiPCA;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedmultiplicity;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedfraction;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSCmultiplicity;
  TH2F *h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsphiPCAMinusPhiTrue_corr;
  
  TH1F *h_hgcal_param_meanTvsy;
  TH1F *h_hgcal_param_meanAlphavsy ; 
  TH1F *h_hgcal_param_meanLnTvsy ;
  TH1F *h_hgcal_param_meanLnAlphavsy;
  TH1F *h_hgcal_param_sigmaLnTvsy;
  TH1F *h_hgcal_param_sigmaLnAlphavsy;
  TH1F *h_hgcal_param_corrAlphaTvsy;
  TH1F *h_hgcal_param_rC;
  TH1F *h_hgcal_param_rT;
  TH1F *h_hgcal_param_p;
    
  TH2F *h_hgcal_sclusters_predictedTvsy;
  TH2F *h_hgcal_sclusters_predictedAlphavsy ; 
  TH2F *h_hgcal_sclusters_predictedLnTvsy ;
  TH2F *h_hgcal_sclusters_predictedLnAlphavsy;
  TH2F *h_hgcal_sclusters_predictedsigmaLnTvsy;
  TH2F *h_hgcal_sclusters_predictedsigmaLnAlphavsy;
  TH2F *h_hgcal_sclusters_predictedcorrAlphaTvsy;
  TH2F *h_hgcal_sclusters_predictedrC;
  TH2F *h_hgcal_sclusters_predictedrT;
  TH2F *h_hgcal_sclusters_predictedp;    
  TH2F *h_hgcal_sclusters_predictedLength;
  TH2F *h_hgcal_sclusters_predictedLength_fullrange;
  
  TH2F *h_hgcal_sclusters_length;
  TH2F *h_hgcal_sclusters_energyVSlength;
  TH2F *h_hgcal_sclusters_length_fullrange;
  TH1F *h_hgcal_sclusters_entryposition;
  TH1F *h_hgcal_sclusters_entryposition_1mip;
  TH1F *h_hgcal_sclusters_entryposition_2mip;
  TH1F *h_hgcal_sclusters_entryposition_4mip;
  TH1F *h_hgcal_sclusters_entryposition_8mip;
  TH1F *h_hgcal_sclusters_entryposition_12mip;
  TH1F *h_hgcal_sclusters_entryposition_16mip;
  TH2F *h_hgcal_sclusters_energyVSlength_fullrange;
  TH2F *h_hgcal_sclusters_energyVSlength_cut_fullrange;
  TH2F *h_hgcal_sclusters_entryVSlength_fullrange;
  TH2F *h_hgcal_sclusters_entryVSlength_cut_fullrange;
  TH2F *h_hgcal_sclusters_energyVSlength_hasfirstlayer_fullrange;
  TH2F *h_hgcal_sclusters_energyVSlength_hasfirstlayer_cut_fullrange;
  TH2F *h_hgcal_sclusters_energyVSmeanradius_fullrange;
  TH2F *h_hgcal_sclusters_entryVSmeanradius_fullrange;
  TH2F *h_hgcal_sclusters_energyVSsigmaradius_fullrange;
  TH2F *h_hgcal_sclusters_entryVSsigmaradius_fullrange;
  TH2F *h_hgcal_sclusters_energyVSlongwidth_fullrange ; 
  TH2F *h_hgcal_sclusters_entryVSlongwidth_fullrange;  
  TH2F *h_hgcal_sclusters_energyVStranswidth_fullrange ;
  TH2F *h_hgcal_sclusters_entryVStranswidth_fullrange;
  TH2F *h_hgcal_sclusters_energyVSeigenratio_fullrange;
  TH2F *h_hgcal_sclusters_entryVSeigenratio_fullrange;
  TH2F *h_hgcal_sclusters_energyVSeigenratio_cut_fullrange;
  TH2F *h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_fullrange;
  TH2F *h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_cut_fullrange;
  TH1F *h_hgcal_sclusters_sigmaradius_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradius_fullrange;
  TH2F *h_hgcal_sclusters_sigmatransverseradiusVSeta_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradius_corr_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradiusaxis_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradiusaxismiddle_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradiusaxismiddle_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmatransverseradiusaxisVSeta_fullrange;
  TH2F *h_hgcal_sclusters_sigmatransverseradiusaxisVSphi_fullrange;
  TH2F *h_hgcal_sclusters_sigmatransverseradiusaxisVSlength_fullrange;
  TH1F *h_hgcal_sclusters_sigmatransverseradiusaxis_corr_fullrange;
  TH1F *h_hgcal_sclusters_sigmaradiusnorm_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetanorm20_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetanorm50_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetanorm100_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta20_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta50_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta100_fullrange;
  TH2F *h_hgcal_sclusters_sigmaradiusVSeta_fullrange;
  TH2F *h_hgcal_sclusters_sigmaradiusVSphi_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta_pu_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetanorm_fullrange;
  TH1F *h_hgcal_sclusters_sigmaphi_fullrange;
  TH1F *h_hgcal_sclusters_sigmaphinorm_fullrange;
  TH2F *h_hgcal_sclusters_sigmaetaVSeta_fullrange;
  TH2F *h_hgcal_sclusters_sigmaeta_puVSeta_fullrange;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_fullrange;
  TH2F *h_hgcal_sclusters_etaVSsigmatransverseradius_fullrange;
  TH2F *h_hgcal_sclusters_etaVSsigmatransverseradiusaxis_fullrange;
  TH2F *h_hgcal_sclusters_sigmaphiVSeta_fullrange;
  TH2F *h_hgcal_sclusters_sigmaphiVSphi_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta_corr_fullrange;
  TH1F *h_hgcal_sclusters_sigmaeta_pu_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaetaVSpt_pu_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaetaVSeta_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaeta_puVSeta_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaeta_puVSEt_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaeta_puVSseedfraction_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaeta10_puVSEt_corr_fullrange;
  TH1F *h_hgcal_sclusters_sigmaphi_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaphicorrVSphiminusphitrue;
  TH2F *h_hgcal_sclusters_sigmaphiVSeta_corr_fullrange;
  TH2F *h_hgcal_sclusters_sigmaphiVSEt_corr_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetaw_fullrange;
  TH2F *h_hgcal_sclusters_sigmaetawVSeta_fullrange;
  TH2F *h_hgcal_sclusters_sigmaetawnormVSeta_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetawnorm_fullrange;
  TH1F *h_hgcal_sclusters_sigmaetaw200_fullrange; 
  TH1F *h_hgcal_sclusters_deta_shower; 
  

  TH1F *h_hgcal_sclusters_sigmaeta_1;
  TH1F *h_hgcal_sclusters_sigmaeta_2;
  TH1F *h_hgcal_sclusters_sigmaeta_3;
  TH1F *h_hgcal_sclusters_sigmaeta_4;
  TH1F *h_hgcal_sclusters_sigmaeta_5;
  TH1F *h_hgcal_sclusters_sigmaeta_6;
  TH1F *h_hgcal_sclusters_sigmaeta_7;
  TH1F *h_hgcal_sclusters_sigmaeta_8;
  TH1F *h_hgcal_sclusters_sigmaeta_9;
  TH1F *h_hgcal_sclusters_sigmaeta_10;
  TH1F *h_hgcal_sclusters_sigmaeta_11;
  TH1F *h_hgcal_sclusters_sigmaeta_12;
  TH1F *h_hgcal_sclusters_sigmaeta_13;
  TH1F *h_hgcal_sclusters_sigmaeta_14;
  TH1F *h_hgcal_sclusters_sigmaeta_15;
  TH1F *h_hgcal_sclusters_sigmaeta_16;
  TH1F *h_hgcal_sclusters_sigmaeta_17;
  TH1F *h_hgcal_sclusters_sigmaeta_18;
  TH1F *h_hgcal_sclusters_sigmaeta_19;
  TH1F *h_hgcal_sclusters_sigmaeta_20;
  TH1F *h_hgcal_sclusters_sigmaeta_21;
  TH1F *h_hgcal_sclusters_sigmaeta_22;
  TH1F *h_hgcal_sclusters_sigmaeta_23;
  TH1F *h_hgcal_sclusters_sigmaeta_24;
  TH1F *h_hgcal_sclusters_sigmaeta_25;
  TH1F *h_hgcal_sclusters_sigmaeta_26;
  TH1F *h_hgcal_sclusters_sigmaeta_27;
  TH1F *h_hgcal_sclusters_sigmaeta_28;
  TH1F *h_hgcal_sclusters_sigmaeta_29;
  TH1F *h_hgcal_sclusters_sigmaeta_30;
  TH2F *h_hgcal_sclusters_sigmaetaVSlayer;
  TH2F *h_hgcal_sclusters_sigmaetaVSlayer_norm;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_1;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_2;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_3;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_4;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_5;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_6;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_7;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_8;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_9;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_10;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_11;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_12;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_13;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_14;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_15;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_16;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_17;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_18;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_19;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_20;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_21;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_22;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_23;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_24;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_25;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_26;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_27;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_28;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_29;
  TH1F *h_hgcal_sclusters_sigmaeta_cutlayer9_30;
  TH2F *h_hgcal_sclusters_sigmaetaVSlayer_cutlayer9;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_1;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_2;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_3;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_4;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_5;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_6;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_7;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_8;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_9;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_10;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_11;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_12;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_13;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_14;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_15;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_16;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_17;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_18;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_19;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_20;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_21;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_22;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_23;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_24;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_25;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_26;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_27;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_28;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_29;
  TH2F *h_hgcal_sclusters_etaVSsigmaeta_30;

  TH1F *h_hgcal_sclusters_longitudinal;
  TH1F *h_hgcal_sclusters_lengthCompatibility;
  TH1F *h_hgcal_sclusters_transversal;
  TH1F *h_hgcal_sclusters_transversalaxis;
  TH1F *h_hgcal_sclusters_transversalaxis_calib;
  TH1F *h_hgcal_sclusters_firsteigenvalue;
  TH1F *h_hgcal_sclusters_firsteigenvalue_nm1;
  TH1F *h_hgcal_sclusters_secondeigenvalue;
  TH1F *h_hgcal_sclusters_secondeigenvalue_nm1;
  TH1F *h_hgcal_sclusters_thirdeigenvalue;
  TH1F *h_hgcal_sclusters_thirdeigenvalue_nm1;
  TH2F *h_hgcal_sclusters_transversalVSeta;
  TH2F *h_hgcal_sclusters_firsteigenvalueVSeta;  
  TH2F *h_hgcal_sclusters_secondeigenvalueVSeta;  
  TH2F *h_hgcal_sclusters_firsteigenvalueVSsecond;  
  TH2F *h_hgcal_sclusters_transeigenvalueVSeta;  
  TH2F *h_hgcal_sclusters_eigenratioVSeta;  
  TH1F *h_hgcal_sclusters_firstsigma;
  TH1F *h_hgcal_sclusters_firstsigma_nm1;
  TH1F *h_hgcal_sclusters_secondsigma;
  TH1F *h_hgcal_sclusters_secondsigma_nm1;
  TH1F *h_hgcal_sclusters_thirdsigma;
  TH1F *h_hgcal_sclusters_thirdsigma_nm1;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinal;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinalscaled;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinalscaledmeasured;
  TH2F *h_hgcal_sclusters_transversalaxisVSlongitudinal;
  TH2F *h_hgcal_sclusters_transversalaxisVSlongitudinalscaled;
  TH2F *h_hgcal_sclusters_transversalaxisVSlongitudinalscaledmeasured;
  TH1F *h_hgcal_sclusters_layer;
  TH1F *h_hgcal_sclusters_longitudinal_cut20;
  TH1F *h_hgcal_sclusters_transversal_cut20;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinal_cut20;
  TH1F *h_hgcal_sclusters_longitudinal_cut50;
  TH1F *h_hgcal_sclusters_transversal_cut50;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinal_cut50;
  TH1F *h_hgcal_sclusters_longitudinal_cut100;
  TH1F *h_hgcal_sclusters_transversal_cut100;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinal_cut100;
  TH1F *h_hgcal_sclusters_longitudinal_cut200;
  TH1F *h_hgcal_sclusters_transversal_cut200;
  TH2F *h_hgcal_sclusters_transversalVSlongitudinal_cut200;
  
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutpos ;
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength ;
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength_cutsigmaeta ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength ;
  TH1F *h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength_cutsigmaeta ;

  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta;             
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem;             
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutkolmogorov;             
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos;             
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength;             
  TH1F *h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength_cutkolmogorov;             

  TH1F *h_hgcal_sclusters_hoverem;
  TH1F *h_hgcal_sclusters_hoverem_cone01;
  TH1F *h_hgcal_sclusters_hoverem_cutsigmaeta;
  TH2F *h_hgcal_sclusters_hoveremVSeta;
  TH2F *h_hgcal_sclusters_hoveremVSphi;
  TH2F *h_hgcal_sclusters_entryVShoverem;
  TH2F *h_hgcal_sclusters_expectedlengthVShoverem;
  TH1F *h_hgcal_sclusters_hoverem1;
  TH2F *h_hgcal_sclusters_hoverem1VSeta;
  TH1F *h_hgcal_sclusters_entry_cutsigmaeta_cuthadem;
  TH2F *h_hgcal_sclusters_entryVShoverem1;  
  TH2F *h_hgcal_sclusters_expectedlengthVShoverem1;
  TH1F *h_hgcal_sclusters_hoverem2;
  TH2F *h_hgcal_sclusters_hoverem2VSeta;
  TH2F *h_hgcal_sclusters_entryVShoverem2;
  TH2F *h_hgcal_sclusters_expectedlengthVShoverem2;
  TH2F *h_hgcal_sclusters_longitudinal_fit_alphaVSbeta;
  TH2F *h_hgcal_sclusters_longitudinal_fit_alphaVSinvbeta;
  TH2F *h_hgcal_sclusters_longitudinal_fit_alphaVSenergy;
  TH2F *h_hgcal_sclusters_longitudinal_fit_betaVSenergy;
  TH2F *h_hgcal_sclusters_longitudinal_fit_invbetaVSenergy;
  TH1F *h_hgcal_sclusters_longitudinal_fit_chi2;
  TH1F *h_hgcal_sclusters_longitudinal_fit_chi2_pnorm;
  TH1F *h_hgcal_sclusters_longitudinal_fit_chi2_bounded;
  TH1F *h_hgcal_sclusters_longitudinal_fit_chi2_nm1;
  TH2F *h_hgcal_sclusters_longitudinal_fit_chi2VSeta;
  TH1F *h_hgcal_sclusters_longitudinal_fit_leakage;
  TH2F *h_hgcal_sclusters_longitudinal_fit_leakageVShoverem;
  TH1F *h_hgcal_sclusters_longitudinal_fit_leakage_cutseedpos;
  TH2F *h_hgcal_sclusters_longitudinal_fit_leakageVShoverem_cutseedpos;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_prob;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_prob_best;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem_cutseedpos;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_nm1;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_best;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem_cutseedpos;
  TH1F *h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_nm1;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_prob;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_dist;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_prob_first;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_dist_first;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_prob_last;
  TH1F *h_hgcal_sclusters_transversal_kolmogorov_dist_last;
  TH2F *h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal;
  TH2F *h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_nm1;
  TH2F *h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cutsigmaeta_cuthadem;
  TH2F *h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem;
  TH2F *h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem_cutseedpos;
  TH1F *h_hgcal_sclusters_kolmogorov_prob_3D;
  TH1F *h_hgcal_sclusters_kolmogorov_dist_3D;
  TH1F *h_hgcal_sclusters_kolmogorov_dist_3D_px;
  TH1F *h_hgcal_sclusters_kolmogorov_dist_3D_py;
  TH2F *h_hgcal_sclusters_longitudinal_fit_chi2VSseedfraction;
  TH2F *h_hgcal_sclusters_longitudinal_fit_chi2VSeoveretrue;
  TH2F *h_hgcal_sclusters_longitudinal_fit_chi2VSdphidir;
  
  TH1F *h_hgcal_sclusters_e2530OverEtot;
  TH1F *h_hgcal_sclusters_e0110OverEtot;
  TH2F *h_hgcal_sclusters_e0110OverEtotVSe2530OverEtot;
  
  TH2F *h_hgcal_sclusters_subclustersVsEnergy;
  TH2F *h_hgcal_sclusters_seedEnergyVsEnergy;
  TH2F *h_hgcal_sclusters_seedEnergyVsEta;
  
  TH2F *h_hgcal_allclusters_predictedTvsy;
  TH2F *h_hgcal_allclusters_predictedLnTvsy ;

  TH1F *h_hgcal_allclusters_transverse_energy_em;
  TH2F *h_hgcal_allclusters_predictedLength;
  TH2F *h_hgcal_allclusters_length;
  TH2F *h_hgcal_allclusters_energyVSlength;
  TH2F *h_hgcal_allclusters_predictedLength_fullrange;
  TH2F *h_hgcal_allclusters_length_fullrange;
  TH2F *h_hgcal_allclusters_energyVSlength_fullrange;
  TH2F *h_hgcal_allclusters_energyVSlength_cut_fullrange;
  TH2F *h_hgcal_allclusters_energyVSpredictedLength_fullrange;
  TH2F *h_hgcal_allclusters_energyVSpredictedLength_cut_fullrange;
  TH2F *h_hgcal_allclusters_entryVSlength_fullrange;
  TH2F *h_hgcal_allclusters_entryVSlength_cut_fullrange;
  TH2F *h_hgcal_allclusters_entryVSpredictedLength_fullrange;
  TH2F *h_hgcal_allclusters_entryVSpredictedLength_cut_fullrange;
  TH2F *h_hgcal_allclusters_energyVSlength_hasfirstlayer_fullrange;
  TH2F *h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange;
  TH2F *h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange_weighted;

  TH1F *h_hgcal_puclusters_energy;
  TH2F *h_hgcal_puclusters_energyVSeta;
  TH1F *h_hgcal_puclusters_energy_eta16;
  TH1F *h_hgcal_puclusters_energy_eta20;
  TH1F *h_hgcal_puclusters_energy_eta25;
  TH1F *h_hgcal_puclusters_energy_eta29;

  TH1F *h_hgcal_allhits_energy_em;

  const HGCalGeometry *geometry_;
  const MagneticField *magField_;
  unsigned long long caloGeomCacheId_;  
  unsigned long long cacheIDMagField_ ;
  
  std::string outputFile_;
  edm::InputTag endcapRecHitCollection_;
  edm::InputTag endcapSuperClusterCollection_;
  edm::InputTag endcapClusterCollection_;
  edm::InputTag electronCollection_;
  edm::InputTag  mcTruthCollection_;
  bool readAOD_;
  bool withPileup_;

  double maxPt_;
  double maxAbsEta_;
  double deltaR_;
  std::vector<int> matchingIDs_;
  std::vector<int> matchingMotherIDs_;

  // histos limits and binning
  double etamin;
  double etamax;
  double phimin;
  double phimax;
  double ptmax;
  double pmax;
  double eopmax;
  double eopmaxsht;
  double detamin;
  double detamax;
  double dphimin;
  double dphimax;
  double detamatchmin;
  double detamatchmax;
  double dphimatchmin;
  double dphimatchmax;
  double fhitsmax;
  double lhitsmax;
  double poptruemin;
  double poptruemax;
  double meemin;
  double meemax;
  double hoemin;
  double hoemax;
  int nbineta;
  int nbinp;
  int nbinpt;
  int nbinpteff;
  int nbinphi;
  int nbinp2D;
  int nbinpt2D;
  int nbineta2D;
  int nbinphi2D;
  int nbineop;
  int nbineop2D;
  int nbinfhits;
  int nbinlhits;
  int nbinxyz;
  int nbindeta;
  int nbindphi;
  int nbindetamatch;
  int nbindphimatch;
  int nbindetamatch2D;
  int nbindphimatch2D;
  int nbinpoptrue;
  int nbinmee;
  int nbinhoe;
  
  //shower parametrisation
  CaloGeometryHelper* calohelper_;
  const RandomEngine* random_ ;
  GammaFunctionGenerator*  aGammaGenerator_;  
  EMECALShowerParametrization* showerparam_;
  std::vector<double> theCoreIntervals_;
  std::vector<double> theTailIntervals_;
  double RCFactor_;
  double RTFactor_;
  
 };

#endif



