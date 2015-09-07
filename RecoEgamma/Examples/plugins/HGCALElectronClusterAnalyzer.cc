
// -*- C++ -*-
//
// Package:    RecoEgamma/Examples
// Class:      HGCALElectronClusterAnalyzer
//
/**\class HGCALElectronClusterAnalyzer RecoEgamma/Examples/src/HGCALElectronClusterAnalyzer.cc

 Description: GsfElectrons analyzer using MC truth

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ursula Berthon
//         Created:  Mon Mar 27 13:22:06 CEST 2006
// $Id: HGCALElectronClusterAnalyzer.cc,v 1.51 2011/03/04 14:43:15 chamont Exp $
//
//

// user include files
#include "RecoEgamma/Examples/plugins/HGCALElectronClusterAnalyzer.h"
#include "RecoEgamma/Examples/interface/PCAShowerAnalysis.h"
#include "RecoEgamma/Examples/interface/HGCALShowerBasedEmIdentification.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "FastSimulation/Utilities/interface/RandomEngine.h"
#include "FastSimulation/Utilities/interface/GammaFunctionGenerator.h"
#include "FastSimulation/Utilities/interface/LandauFluctuationGenerator.h"
#include "FastSimulation/Event/interface/FSimTrack.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include <iostream>
#include <vector>
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"

#include <iostream>

using namespace reco;


// static Double_t gamma(Double_t *x, Double_t *par) {
//   if (par[1]<=0. || par[2]<=0.)return 0.;
//   else return par[0]*TMath::GammaDist(x[0],par[1],0.,par[2]);
// }
// 
// static Double_t transdetector(Double_t *x, Double_t *par) {
//   // par[4] is costheta where theta is the shower axis
//   // transversedetector = transverse / costheta
//   Double_t core = x[0]*par[2]*exp(-par[2]*par[4]*x[0]);
//   Double_t tail = x[0]*par[3]*exp(-par[3]*par[4]*x[0]);
//   //Double_t core = 2.*x[0]*par[2]*exp(-par[2]*x[0]);
//   //Double_t tail = 2.*x[0]*par[3]*exp(-par[3]*x[0]);
//   //Double_t core = 2.*x[0]*exp(-x[0]/par[2]);
//   //Double_t tail = 2.*x[0]*exp(-x[0]/par[3]);
//   //return par[0]*(par[1]*core + (1.-par[1])*tail);
//   return par[0]*(par[1]*core + (1.-par[1])*tail)/(par[4]*par[4]);
//   //return par[0]*(par[1]*core + (1.-par[1])*tail);
// }

/*
static Double_t shower3D(Double_t *x, Double_t *par) {
  // par[10] is costheta where theta is the shower axis
  Double_t longitudinal;
  if (par[1]<0. || par[2]<0.)return 0.;
  else longitudinal = TMath::GammaDist(x[0],par[1],0.,par[2]);
  Double_t reduced = x[0]/((par[1]-1.)*par[2]); // tmax = (alpha-1)/beta
  Double_t corefraction = par[3]+par[4]*reduced+par[5]*reduced*reduced;
  //corefraction *= 1.30;
  Double_t coreradius = (par[6]+par[7]*reduced)*2.27*0.80; // in cm
  Double_t tailradius = (par[8]+par[9]*reduced)*2.27*0.80; // in cm
  Double_t core = (x[1]*par[10]/coreradius)*std::exp(-x[1]*par[10]/coreradius);
  Double_t tail = (x[1]*par[10]/tailradius)*std::exp(-x[1]*par[10]/tailradius);
  Double_t transverse = (corefraction*core + (1.-corefraction)*tail)*par[10];
  return par[0]*longitudinal*transverse*2.*TMath::Pi();
}
*/
HGCALElectronClusterAnalyzer::HGCALElectronClusterAnalyzer(const edm::ParameterSet& conf)
 {
  outputFile_ = conf.getParameter<std::string>("outputFile");
  histfile_ = new TFile(outputFile_.c_str(),"RECREATE");
  electronCollection_=conf.getParameter<edm::InputTag>("electronCollection");
  endcapRecHitCollection_=conf.getParameter<edm::InputTag>("endcapRecHitCollection");
  endcapSuperClusterCollection_=conf.getParameter<edm::InputTag>("endcapSuperClusterCollection");
  endcapClusterCollection_=conf.getParameter<edm::InputTag>("endcapClusterCollection");
  mcTruthCollection_ = conf.getParameter<edm::InputTag>("mcTruthCollection");
  withPileup_ = conf.getParameter<bool>("withPileup");
  readAOD_ = conf.getParameter<bool>("readAOD");
  maxPt_ = conf.getParameter<double>("MaxPt");
  maxAbsEta_ = conf.getParameter<double>("MaxAbsEta");
  deltaR_ = conf.getParameter<double>("DeltaR");
  matchingIDs_ = conf.getParameter<std::vector<int> >("MatchingID");
  matchingMotherIDs_ = conf.getParameter<std::vector<int> >("MatchingMotherID");
  edm::ParameterSet pset =
   conf.getParameter<edm::ParameterSet>("HistosConfigurationMC") ;
  caloGeomCacheId_=0;
  cacheIDMagField_=0;
  nevt=-1;
  ievent=-1;
  
  etamin=pset.getParameter<double>("Etamin");
  etamax=pset.getParameter<double>("Etamax");
  phimin=pset.getParameter<double>("Phimin");
  phimax=pset.getParameter<double>("Phimax");
  ptmax=pset.getParameter<double>("Ptmax");
  pmax=pset.getParameter<double>("Pmax");
  eopmax=pset.getParameter<double>("Eopmax");
  eopmaxsht=pset.getParameter<double>("Eopmaxsht");
  detamin=pset.getParameter<double>("Detamin");
  detamax=pset.getParameter<double>("Detamax");
  dphimin=pset.getParameter<double>("Dphimin");
  dphimax=pset.getParameter<double>("Dphimax");
  detamatchmin=pset.getParameter<double>("Detamatchmin");
  detamatchmax=pset.getParameter<double>("Detamatchmax");
  dphimatchmin=pset.getParameter<double>("Dphimatchmin");
  dphimatchmax=pset.getParameter<double>("Dphimatchmax");
  fhitsmax=pset.getParameter<double>("Fhitsmax");
  lhitsmax=pset.getParameter<double>("Lhitsmax");
  nbineta=pset.getParameter<int>("Nbineta");
  nbineta2D=pset.getParameter<int>("Nbineta2D");
  nbinp=pset.getParameter<int>("Nbinp");
  nbinpt=pset.getParameter<int>("Nbinpt");
  nbinp2D=pset.getParameter<int>("Nbinp2D");
  nbinpt2D=pset.getParameter<int>("Nbinpt2D");
  nbinpteff=pset.getParameter<int>("Nbinpteff");
  nbinphi=pset.getParameter<int>("Nbinphi");
  nbinphi2D=pset.getParameter<int>("Nbinphi2D");
  nbineop=pset.getParameter<int>("Nbineop");
  nbineop2D=pset.getParameter<int>("Nbineop2D");
  nbinfhits=pset.getParameter<int>("Nbinfhits");
  nbinlhits=pset.getParameter<int>("Nbinlhits");
  nbinxyz=pset.getParameter<int>("Nbinxyz");
  nbindeta=pset.getParameter<int>("Nbindeta");
  nbindphi=pset.getParameter<int>("Nbindphi");
  nbindetamatch=pset.getParameter<int>("Nbindetamatch");
  nbindphimatch=pset.getParameter<int>("Nbindphimatch");
  nbindetamatch2D=pset.getParameter<int>("Nbindetamatch2D");
  nbindphimatch2D=pset.getParameter<int>("Nbindphimatch2D");
  nbinpoptrue= pset.getParameter<int>("Nbinpoptrue");
  poptruemin=pset.getParameter<double>("Poptruemin");
  poptruemax=pset.getParameter<double>("Poptruemax");
  nbinmee= pset.getParameter<int>("Nbinmee");
  meemin=pset.getParameter<double>("Meemin");
  meemax=pset.getParameter<double>("Meemax");
  nbinhoe= pset.getParameter<int>("Nbinhoe");
  hoemin=pset.getParameter<double>("Hoemin");
  hoemax=pset.getParameter<double>("Hoemax");
  
  // shower parametrization
  readParameters(conf.getParameter<edm::ParameterSet>("Calorimetry"));
  
  calohelper_ = new
   CaloGeometryHelper(conf.getParameter<edm::ParameterSet>("Calorimetry"));
  // random generators
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable() ) {
    throw cms::Exception("Configuration")
    << "FamosManager requires the RandomGeneratorService\n"
    "which is not present in the configuration file.\n"
    "You must add the service in the configuration file\n"
    "or remove the module that requires it";
  }
  random_ = new RandomEngine(&(*rng));  
  aGammaGenerator_ = new GammaFunctionGenerator(random_);  
  // The ECAL Properties
//   showerparam_ = new EMECALShowerParametrization(
//    calohelper_->ecalProperties(onEcal),
//    calohelper_->hcalProperties(onHcal),
//    calohelper_->layer1Properties(onLayer1),
//    calohelper_->layer2Properties(onLayer2),
//    theCoreIntervals_,
//    theTailIntervals_,
//    RCFactor_,
//    RTFactor_);
  showerparam_ = new EMECALShowerParametrization(
   calohelper_->ecalProperties(1),
   calohelper_->hcalProperties(0),
   calohelper_->layer1Properties(0),
   calohelper_->layer2Properties(0),
   theCoreIntervals_,
   theTailIntervals_,
   RCFactor_,
   RTFactor_);
   
 }

void HGCALElectronClusterAnalyzer::readParameters(const edm::ParameterSet& fastCalo) {

  edm::ParameterSet ECALparameters = fastCalo.getParameter<edm::ParameterSet>("ECAL");
  edm::ParameterSet CalorimeterProperties = fastCalo.getParameter<edm::ParameterSet>("CalorimeterProperties");
  edm::ParameterSet BarrelCalorimeterProperties = CalorimeterProperties.getParameter<edm::ParameterSet>("BarrelCalorimeterProperties");

//   evtsToDebug_ = fastCalo.getUntrackedParameter<std::vector<unsigned int> >("EvtsToDebug",std::vector<unsigned>());
//   debug_ = fastCalo.getUntrackedParameter<bool>("Debug");
//   useDQM_ = fastCalo.getUntrackedParameter<bool>("useDQM");
//  double bFixedLength = ECALparameters.getParameter<bool>("bFixedLength");
//  gridSize_ = ECALparameters.getParameter<int>("GridSize");
//  spotFraction_ = ECALparameters.getParameter<double>("SpotFraction");
  double bHom =  BarrelCalorimeterProperties.getParameter<bool>("bHom");
  std::cout << "CalorimeterProperties.BarrelCalorimeterProperties.bHom = " << bHom << std::endl;
//   pulledPadSurvivalProbability_ = ECALparameters.getParameter<double>("FrontLeakageProbability");
//   crackPadSurvivalProbability_ = ECALparameters.getParameter<double>("GapLossProbability");

  // we need those ones
  theCoreIntervals_ = ECALparameters.getParameter<std::vector<double> >("CoreIntervals");
  theTailIntervals_ = ECALparameters.getParameter<std::vector<double> >("TailIntervals");
  RCFactor_ = ECALparameters.getParameter<double>("RCFactor");
  RTFactor_ = ECALparameters.getParameter<double>("RTFactor");
/*
  //changed after tuning - Feb-July - Shilpi Jain
  // radiusFactor_ = ECALparameters.getParameter<double>("RadiusFactor");
  radiusFactorEE_ = ECALparameters.getParameter<double>("RadiusFactorEE");
  radiusFactorEB_ = ECALparameters.getParameter<double>("RadiusFactorEB");
  //(end of) changed after tuning - Feb-July - Shilpi Jain
  radiusPreshowerCorrections_ = ECALparameters.getParameter<std::vector<double> >("RadiusPreshowerCorrections");
  aTerm = 1.+radiusPreshowerCorrections_[1]*radiusPreshowerCorrections_[0];
  bTerm = radiusPreshowerCorrections_[0];
  mipValues_ = ECALparameters.getParameter<std::vector<double> >("MipsinGeV");
  simulatePreshower_ = ECALparameters.getParameter<bool>("SimulatePreshower");
  if(gridSize_ <1) gridSize_= 7;
  if(pulledPadSurvivalProbability_ <0. || pulledPadSurvivalProbability_>1 ) pulledPadSurvivalProbability_= 1.;
  if(crackPadSurvivalProbability_ <0. || crackPadSurvivalProbability_>1 ) crackPadSurvivalProbability_= 0.9;
  LogInfo("FastCalorimetry") << " Fast ECAL simulation parameters " << std::endl;
  LogInfo("FastCalorimetry") << " =============================== " << std::endl;
  if(simulatePreshower_)
    LogInfo("FastCalorimetry") << " The preshower is present " << std::endl;
  else
    LogInfo("FastCalorimetry") << " The preshower is NOT present " << std::endl;
  LogInfo("FastCalorimetry") << " Grid Size : " << gridSize_ << std::endl;
  if(spotFraction_>0.)
    LogInfo("FastCalorimetry") << " Spot Fraction : " << spotFraction_ << std::endl;
  else
  {
    LogInfo("FastCalorimetry") << " Core of the shower " << std::endl;
    for(unsigned ir=0; ir < theCoreIntervals_.size()/2;++ir)
    {
      LogInfo("FastCalorimetry") << " r < " << theCoreIntervals_[ir*2] << " R_M : " << theCoreIntervals_[ir*2+1] << " ";
    }
    LogInfo("FastCalorimetry") << std::endl;
    LogInfo("FastCalorimetry") << " Tail of the shower " << std::endl;
    for(unsigned ir=0; ir < theTailIntervals_.size()/2;++ir)
    {
      LogInfo("FastCalorimetry") << " r < " << theTailIntervals_[ir*2] << " R_M : " << theTailIntervals_[ir*2+1] << " ";
    }
    //changed after tuning - Feb-July - Shilpi Jain
    // LogInfo("FastCalorimetry") << "Radius correction factor " << radiusFactor_ << std::endl;
    LogInfo("FastCalorimetry") << "Radius correction factors: EB & EE " << radiusFactorEB_ << " : "<< radiusFactorEE_ << std::endl;
    //(end of) changed after tuning - Feb-July - Shilpi Jain
    LogInfo("FastCalorimetry") << std::endl;
    if(mipValues_.size()>2) {
      LogInfo("FastCalorimetry") << "Improper number of parameters for the preshower ; using 95keV" << std::endl;
      mipValues_.clear();
      mipValues_.resize(2,0.000095);
    }
  }
  LogInfo("FastCalorimetry") << " FrontLeakageProbability : " << pulledPadSurvivalProbability_ << std::endl;
  LogInfo("FastCalorimetry") << " GapLossProbability : " << crackPadSurvivalProbability_ << std::endl;
  // RespCorrP: p (momentum), ECAL and HCAL corrections = f(p)
  edm::ParameterSet CalorimeterParam = fastCalo.getParameter<edm::ParameterSet>("CalorimeterProperties");
  rsp = CalorimeterParam.getParameter<std::vector<double> >("RespCorrP");
  LogInfo("FastCalorimetry") << " RespCorrP (rsp) size " << rsp.size() << std::endl;
  if( rsp.size()%3 !=0 ) {
    LogInfo("FastCalorimetry")
    << " RespCorrP size is wrong -> no corrections applied !!!"
    << std::endl;
    p_knots.push_back(14000.);
    k_e.push_back (1.);
    k_h.push_back (1.);
  }
  else {
    for(unsigned i = 0; i < rsp.size(); i += 3) {
      LogInfo("FastCalorimetry") << "i = " << i/3 << " p = " << rsp [i]
      << " k_e(p) = " << rsp[i+1]
      << " k_e(p) = " << rsp[i+2] << std::endl;
      p_knots.push_back(rsp[i]);
      k_e.push_back (rsp[i+1]);
      k_h.push_back (rsp[i+2]);
    }
  }
  //FR
  edm::ParameterSet HCALparameters = fastCalo.getParameter<edm::ParameterSet>("HCAL");
  optionHDSim_ = HCALparameters.getParameter<int>("SimOption");
  hdGridSize_ = HCALparameters.getParameter<int>("GridSize");
  hdSimMethod_ = HCALparameters.getParameter<int>("SimMethod");
  //RF
  EcalDigitizer_ = ECALparameters.getUntrackedParameter<bool>("Digitizer",false);
  HcalDigitizer_ = HCALparameters.getUntrackedParameter<bool>("Digitizer",false);
  samplingHBHE_ = HCALparameters.getParameter< std::vector<double> >("samplingHBHE");
  samplingHF_ = HCALparameters.getParameter< std::vector<double> >("samplingHF");
  samplingHO_ = HCALparameters.getParameter< std::vector<double> >("samplingHO");
  ietaShiftHB_ = HCALparameters.getParameter< int >("ietaShiftHB");
  ietaShiftHE_ = HCALparameters.getParameter< int >("ietaShiftHE");
  ietaShiftHF_ = HCALparameters.getParameter< int >("ietaShiftHF");
  ietaShiftHO_ = HCALparameters.getParameter< int >("ietaShiftHO");
  timeShiftHB_ = HCALparameters.getParameter< std::vector<double> >("timeShiftHB");
  timeShiftHE_ = HCALparameters.getParameter< std::vector<double> >("timeShiftHE");
  timeShiftHF_ = HCALparameters.getParameter< std::vector<double> >("timeShiftHF");
  timeShiftHO_ = HCALparameters.getParameter< std::vector<double> >("timeShiftHO");
*/
}

void HGCALElectronClusterAnalyzer::beginJob(){

  histfile_->cd();

  // mc truth
  h_mcNum              = new TH1F( "h_mcNum",              "# mc particles",    nbinfhits,0.,fhitsmax );
  h_mcNum->Sumw2();
  h_simEta             = new TH1F( "h_mc_eta",             "gen #eta",           nbineta,etamin,etamax);
  h_simEta->Sumw2();
  h_simAbsEta             = new TH1F( "h_mc_abseta",             "gen |#eta|",           nbineta/2,0.,etamax);
  h_simAbsEta->Sumw2();
  h_simP               = new TH1F( "h_mc_P",               "gen p",              nbinp,0.,pmax);
  h_simP->Sumw2();
  h_simPt               = new TH1F( "h_mc_Pt",               "gen pt",            nbinpteff,5.,ptmax);
  h_simPt->Sumw2();
  h_simPhi               = new TH1F( "h_mc_phi",               "gen phi",        nbinphi,phimin,phimax);
  h_simPhi->Sumw2();
  h_simZ      = new TH1F( "h_mc_z",      "gen z ",    nbinxyz, -25, 25 );
  h_simZ->Sumw2();
  h_simPtEta           = new TH2F( "h_mc_pteta",   "gen pt vs #eta", nbineta2D,etamin,etamax, nbinpt2D,5.,ptmax );
  h_simPtEta->Sumw2();

  h_simZ_all      = new TH1F( "h_mc_z_all",      "gen z all particles",    nbinxyz, -25, 25 );
  h_simZ_electrons      = new TH1F( "h_mc_z_electrons",      "gen z electrons",    nbinxyz, -25, 25 );

  // HGCAL clusters histos

  // first in em part
  h_hgcal_clusterEnergy_em = new TH1F("h_hgcal_clusterEnergy_em","hgcal em cluster energy",100,0.,250.);
  h_hgcal_clusterEnergy_em              -> GetXaxis()-> SetTitle("Energy (GeV)");
  h_hgcal_clusterEnergy_em              -> GetYaxis()-> SetTitle("Events");
  h_hgcal_clusterEtaVsPhi_em = new TH2F( "h_hgcal_clusterEtaVsPhi_em","hgcal em clusters eta vs phi",nbineta2D,etamin,etamax,nbinphi2D,phimin,phimax );
  h_hgcal_clusterEtaVsPhi_em-> GetYaxis()-> SetTitle("#phi (rad)");
  h_hgcal_clusterEtaVsPhi_em-> GetXaxis()-> SetTitle("#eta");

  h_hgcal_foundClusters_em = new TH1F("h_hgcal_foundClusters_em","# hgcal em superclusters",50,0.,50.);
  h_hgcal_foundClusters_em              -> GetXaxis()-> SetTitle("N_{clusters}");
  h_hgcal_foundClusters_em              -> GetYaxis()-> SetTitle("Events");
  h_hgcal_foundClustersVSeta_em = new TH1F("h_hgcal_foundClustersVSeta_em","# hgcal em superclusters vs eta",100,-3.0,3.0);
  h_hgcal_foundClustersVSetaEt5_em = new TH1F("h_hgcal_foundClustersVSetaEt5_em","# hgcal em superclusters vs eta Et cut 5",100,-3.0,3.0);
  h_hgcal_sclusters_energy_em = new TH1F("h_hgcal_sclusters_energy_em","hgcal em supercluster energy",100,0.,500.);
  h_hgcal_sclusters_transverse_energy_em = new TH1F("h_hgcal_sclusters_transverse_energy_em","hgcal em supercluster transverse energy",100,0.,50.);
  //h_hgcal_sclusters_energy_em = new TH1F("h_hgcal_sclusters_energy_em","hgcal em supercluster energy",100,0.,5000.);
  h_hgcal_sclusters_energy_pos_em = new TH1F("h_hgcal_sclusters_energy_pos_em","hgcal em supercluster energy eta>0",100,0.,500.);
  h_hgcal_sclusters_energy_neg_em = new TH1F("h_hgcal_sclusters_energy_neg_em","hgcal em supercluster energy eta<0",100,0.,500.);
  h_hgcal_sclusters_energyVSeta_em = new TH2F("h_hgcal_sclusters_energyVSeta_em","hgcal em supercluster energy vs eta",100,0.,500.,100,-3.0,3.0);
  h_hgcal_sclusters_position_em = new TH2F("h_hgcal_sclusters_position_em","hgcal em supercluster eta vs phi",120,-3.0,3.0,100,-3.14,3.14);
  h_hgcal_sclusters_etawidth_em = new TH1F("h_hgcal_sclusters_etawidth_em","hgcal em supercluster eta width",100,0.,10.);
  h_hgcal_sclusters_phiwidth_em = new TH1F("h_hgcal_sclusters_phiwidth_em","hgcal em supercluster phi width",100,0.,10.);
  h_hgcal_sclusters_multiplicity_em = new TH1F("h_hgcal_sclusters_multiplicity_em","hgcal em supercluster multiplicity",50,0.,50.);
  h_hgcal_sclusters_newmultiplicity_em = new TH1F("h_hgcal_sclusters_newmultiplicity_em","hgcal em supercluster corrected multiplicity",50,0.,50.);
  h_hgcal_sclusters_multiplicityVSeta_em = new TH2F("h_hgcal_sclusters_multiplicityVSeta_em","hgcal em supercluster multiplicity vs eta",50,0.,50., 100.,-3.0,3.0);
  h_hgcal_sclusters_newmultiplicityVSeta_em = new TH2F("h_hgcal_sclusters_newmultiplicityVSeta_em","hgcal em supercluster corrected  multiplicity vs eta",50,0.,50., 100.,-3.0,3.0);
  h_hgcal_sclusters_seedenergy_em = new TH1F("h_hgcal_sclusters_seedenergy_em","hgcal em supercluster seed energy",100,0.,500.);
  h_hgcal_sclusters_newseedenergy_em = new TH1F("h_hgcal_sclusters_newseedenergy_em","hgcal em supercluster corrected seed energy",100,0.,500.);
  //h_hgcal_sclusters_seedenergy_em = new TH1F("h_hgcal_sclusters_seedenergy_em","hgcal em supercluster seed energy",100,0.,5000.);

  h_hgcal_clusters_energy_em = new TH1F("h_hgcal_clusters_energy_em","hgcal em subcluster energy",100,0.,500.);
  //h_hgcal_clusters_energy_em = new TH1F("h_hgcal_clusters_energy_em","hgcal em subcluster energy",100,0.,5000.);
  h_hgcal_clusters_position_em = new TH2F("h_hgcal_clusters_position_em","hgcal em subcluster eta vs phi",120,-3.0,3.0,100,-3.14,3.14);
  h_hgcal_clusters_multiplicity_em = new TH1F("h_hgcal_clusters_multiplicity_em","hgcal em subcluster multiplicity",100,0.,5000.);
  h_hgcal_clusters_multiplicityVSeta_em = new TH2F("h_hgcal_clusters_multiplicityVSeta_em","hgcal em subcluster multiplicity vs eta",100,0.,5000., 100, -3.0,3.0);
  h_hgcal_clusters_rechitenergy_em = new TH1F("h_hgcal_clusters_rechitenergy_em","hgcal em subcluster rechit energy",1200,0.,0.12);
  h_hgcal_clusters_rechitenergy_12000_em = new TH1F("h_hgcal_clusters_rechitenergy_12000_em","hgcal em subcluster rechit energy",12000,0.,0.12);

  h_hgcal_allhits_energy_em = new TH1F("h_hgcal_all_hits_energy_em","hgcal em rechit energy",12000,0.,0.12);
  
  // overlap
  h_hgcal_sclusters_seedfractions_em = new TH1F("h_hgcal_sclusters_seedfractions_em","hgcal em supercluster seed cluster fractions",100,0.,100.);
  h_hgcal_sclusters_seedfractionVSeta_em = new TH2F("h_hgcal_sclusters_seedfractionVseta_em","hgcal em supercluster seed energy fraction vs eta",100,0.,1.,100,-3.0,3.0);
  
  // then for hcal part
  h_hgcal_foundClusters_had = new TH1F("h_hgcal_foundClusters_had","# hgcal had clusters",50,0.,50.);
  h_hgcal_foundClusters_had              -> GetXaxis()-> SetTitle("N_{clusters}");
  h_hgcal_foundClusters_had              -> GetYaxis()-> SetTitle("Events");
  h_hgcal_clusterEnergy_had = new TH1F("h_hgcal_clusterEnergy_had","hgcal had cluster energy",100,0.,250.);
  h_hgcal_clusterEnergy_had              -> GetXaxis()-> SetTitle("Energy (GeV)");
  h_hgcal_clusterEnergy_had              -> GetYaxis()-> SetTitle("Events");
  h_hgcal_clusterEtaVsPhi_had = new TH2F( "h_hgcal_clusterEtaVsPhi_had","hgcal had clusters eta vs phi",nbineta2D,etamin,etamax,nbinphi2D,phimin,phimax );
  h_hgcal_clusterEtaVsPhi_had-> GetYaxis()-> SetTitle("#phi (rad)");
  h_hgcal_clusterEtaVsPhi_had-> GetXaxis()-> SetTitle("#eta");

  // supercluster energy reconstruction
  h_hgcal_scclusters_eoveretrue_em = new TH1F("h_hgcal_scclusters_eoveretrue_em","hgcal em supercluster E/Etrue",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrueVSeta_em = new TH2F("h_hgcal_scclusters_eoveretrueVSeta_em","hgcal em supercluster E/Etrue",130,0.0,1.3,60,1.5,3.0);
  h_hgcal_scclusters_eoveretrue_cut00_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut00_em","hgcal em supercluster E/Etrue, no cut",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut04_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut04_em","hgcal em supercluster E/Etrue, Et cut 0.4 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut1_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut1_em","hgcal em supercluster E/Etrue, E cut 1 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut2_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut2_em","hgcal em supercluster E/Etrue, E cut 2 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut4_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut4_em","hgcal em supercluster E/Etrue, E cut 4 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut10_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut10_em","hgcal em supercluster E/Etrue, E cut 10 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_cut20_em = new TH1F("h_hgcal_scclusters_eoveretrue_cut20_em","hgcal em supercluster E/Etrue, Et cut 20 mip",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_golden_em = new TH1F("h_hgcal_scclusters_eoveretrue_golden_em","hgcal em supercluster E/Etrue, golden",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut_em = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut_em","hgcal em supercluster E/Etrue, noise cut",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut4_em = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut4_em","hgcal em supercluster E/Etrue, noise cut",130,0.0,1.3);
  h_hgcal_scclusters_ptoverpttrue_em = new TH1F("h_hgcal_scclusters_ptoverpttrue_em","hgcal em supercluster pt/ptrue",130,0.0,1.3);
  h_hgcal_scclusters_detadphisubclusters_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_em","hgcal em subcluster deta vs dphi",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_zoom_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_zoom_em","hgcal em subcluster deta vs dphi",100,-0.05,0.05,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_zoom_etagt2_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_zoom_etagt2_em","hgcal em subcluster deta vs dphi eta gt 2",100,-0.05,0.05,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_zoom_etalt2_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_zoom_etalt2_em","hgcal em subcluster deta vs dphi eta lt 2",100,-0.05,0.05,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_weighted_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_weighted_em","hgcal em subcluster deta vs dphi, weighted",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_zoom_weighted_em = new TH2F("h_hgcal_scclusters_detadphisubclusters_zoom_weighted_em","hgcal em subcluster deta vs dphi, weigthed",100,-0.05,0.05,120,-0.3,0.3);
 
  // hits in HGCAL
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[0]","electron hits in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[1]","electron hits in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[2]","electron hits in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[3]","electron hits in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[4]","electron hits in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[5]","electron hits in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[6]","electron hits in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[7]","electron hits in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[8]","electron hits in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_em.push_back(new TH3F("h_hgcal_allhits_em[9]","electron hits in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits in HGCAL (weighted)
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[0]","electron hits in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[1]","electron hits in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[2]","electron hits in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[3]","electron hits in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[4]","electron hits in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[5]","electron hits in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[6]","electron hits in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[7]","electron hits in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[8]","electron hits in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_allhits_weighted_em.push_back(new TH3F("h_hgcal_allhits_weighted_em[9]","electron hits in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from pf clusters in HGCAL
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[0]","electron hits in pfclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[1]","electron hits in pfclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[2]","electron hits in pfclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[3]","electron hits in pfclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[4]","electron hits in pfclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[5]","electron hits in pfclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[6]","electron hits in pfclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[7]","electron hits in pfclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[8]","electron hits in pfclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_em.push_back(new TH3F("h_hgcal_clustershits_em[9]","electron hits in pfclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from pf clusters in HGCAL (weighted)
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[0]","electron hits in pfclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[1]","electron hits in pfclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[2]","electron hits in pfclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[3]","electron hits in pfclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[4]","electron hits in pfclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[5]","electron hits in pfclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[6]","electron hits in pfclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[7]","electron hits in pfclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[8]","electron hits in pfclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_weighted_em.push_back(new TH3F("h_hgcal_clustershits_weighted_em[9]","electron hits in pfclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from pf clusters in HGCAL after cut first layers
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[0]","electron hits in pfclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[1]","electron hits in pfclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[2]","electron hits in pfclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[3]","electron hits in pfclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[4]","electron hits in pfclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[5]","electron hits in pfclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[6]","electron hits in pfclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[7]","electron hits in pfclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[8]","electron hits in pfclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_em[9]","electron hits in pfclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from pf clusters in HGCAL (weighted) after cut first layers
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[0]","electron hits in pfclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[1]","electron hits in pfclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[2]","electron hits in pfclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[3]","electron hits in pfclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[4]","electron hits in pfclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[5]","electron hits in pfclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[6]","electron hits in pfclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[7]","electron hits in pfclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[8]","electron hits in pfclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_weighted_em[9]","electron hits in pfclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from pf clusters in HGCAL after cut first layers and shower length
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[0]","electron hits in pfclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[1]","electron hits in pfclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[2]","electron hits in pfclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[3]","electron hits in pfclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[4]","electron hits in pfclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[5]","electron hits in pfclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[6]","electron hits in pfclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[7]","electron hits in pfclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[8]","electron hits in pfclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_em[9]","electron hits in pfclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from pf clusters in HGCAL (weighted) after cut first layers and shower length
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[0]","electron hits in pfclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[1]","electron hits in pfclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[2]","electron hits in pfclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[3]","electron hits in pfclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[4]","electron hits in pfclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[5]","electron hits in pfclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[6]","electron hits in pfclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[7]","electron hits in pfclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[8]","electron hits in pfclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_clustershits_cutlayers_cutlength_weighted_em[9]","electron hits in pfclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from super clusters in HGCAL
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[0]","electron hits in superclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[1]","electron hits in superclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[2]","electron hits in superclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[3]","electron hits in superclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[4]","electron hits in superclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[5]","electron hits in superclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[6]","electron hits in superclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[7]","electron hits in superclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[8]","electron hits in superclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_em.push_back(new TH3F("h_hgcal_superclustershits_em[9]","electron hits in superclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from super clusters in HGCAL (weighted)
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[0]","electron hits in superclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[1]","electron hits in superclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[2]","electron hits in superclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[3]","electron hits in superclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[4]","electron hits in superclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[5]","electron hits in superclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[6]","electron hits in superclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[7]","electron hits in superclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[8]","electron hits in superclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_weighted_em[9]","electron hits in superclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from super clusters in HGCAL
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[0]","electron hits in superclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[1]","electron hits in superclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[2]","electron hits in superclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[3]","electron hits in superclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[4]","electron hits in superclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[5]","electron hits in superclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[6]","electron hits in superclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[7]","electron hits in superclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[8]","electron hits in superclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_em[9]","electron hits in superclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from super clusters in HGCAL (weighted)
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[0]","electron hits in superclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[1]","electron hits in superclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[2]","electron hits in superclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[3]","electron hits in superclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[4]","electron hits in superclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[5]","electron hits in superclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[6]","electron hits in superclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[7]","electron hits in superclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[8]","electron hits in superclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_weighted_em[9]","electron hits in superclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
   // hits from super clusters in HGCAL
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[0]","electron hits in superclusters in hgcal  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[1]","electron hits in superclusters in hgcal  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[2]","electron hits in superclusters in hgcal  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[3]","electron hits in superclusters in hgcal  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[4]","electron hits in superclusters in hgcal  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[5]","electron hits in superclusters in hgcal  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[6]","electron hits in superclusters in hgcal  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[7]","electron hits in superclusters in hgcal  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[8]","electron hits in superclusters in hgcal  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_em[9]","electron hits in superclusters in hgcal  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  // hits from super clusters in HGCAL (weighted)
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[0]","electron hits in superclusters in hgcal, weighted  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[1]","electron hits in superclusters in hgcal, weighted  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[2]","electron hits in superclusters in hgcal, weighted  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[3]","electron hits in superclusters in hgcal, weighted  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[4]","electron hits in superclusters in hgcal, weighted  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[5]","electron hits in superclusters in hgcal, weighted  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[6]","electron hits in superclusters in hgcal, weighted  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[7]","electron hits in superclusters in hgcal, weighted  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[8]","electron hits in superclusters in hgcal, weighted  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em.push_back(new TH3F("h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[9]","electron hits in superclusters in hgcal, weighted  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
 // few em showers
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[0]","electron shower seed  0", 120,0.,120.,120,0.,120.,350.,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[1]","electron shower seed  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[2]","electron shower seed  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[3]","electron shower seed  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[4]","electron shower seed  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[5]","electron shower seed  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[6]","electron shower seed  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[7]","electron shower seed  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[8]","electron shower seed  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_seed_em.push_back(new TH3F("h_hgcal_shower_seed_em[9]","electron shower seed  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[0]","electron shower sc  0", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[1]","electron shower sc  1", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[2]","electron shower sc  2", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[3]","electron shower sc  3", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[4]","electron shower sc  4", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[5]","electron shower sc  5", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[6]","electron shower sc  6", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[7]","electron shower sc  7", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[8]","electron shower sc  8", 120,0.,120.,120,0.,120.,350,315.,350.));	
  h_hgcal_shower_sc_em.push_back(new TH3F("h_hgcal_shower_sc_em[9]","electron shower sc  9", 120,0.,120.,120,0.,120.,350,315.,350.));	
   
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[0]","electron shower sc  0", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[1]","electron shower sc  1", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[2]","electron shower sc  2", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[3]","electron shower sc  3", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[4]","electron shower sc  4", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[5]","electron shower sc  5", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[6]","electron shower sc  6", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[7]","electron shower sc  7", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[8]","electron shower sc  8", 200,0.,200.,200,0.,200.,70,315.,385.));	
  h_hgcal_shower_sc.push_back(new TH3F("h_hgcal_shower_sc[9]","electron shower sc  9", 200,0.,200.,200,0.,200.,70,315.,385.));	
   
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[0]","electron shower seed rotated 0", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[1]","electron shower seed rotated 1", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[2]","electron shower seed rotated 2", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[3]","electron shower seed rotated 3", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[4]","electron shower seed rotated 4", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[5]","electron shower seed rotated 5", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[6]","electron shower seed rotated 6", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[7]","electron shower seed rotated 7", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[8]","electron shower seed rotated 8", 35,315.,350.,50,0.,50.,50,0.,50.));	
  h_hgcal_shower_seed_rotated_em.push_back(new TH3F("h_hgcal_shower_seed_rotated_em[9]","electron shower seed rotated 9", 35,315.,350.,50,0.,50.,50,0.,50.));	

  h_hgcal_sclusters_etaPCAMinusEtaTrue = new TH1F("h_hgcal_sclusters_etaPCAMinusEtaTrue","hgcal em supercluster eta position from PCA - eta true",200,-0.1,0.1);
  h_hgcal_sclusters_phiPCAMinusPhiTrue = new TH1F("h_hgcal_sclusters_phiPCAMinusPhiTrue","hgcal em supercluster phi position from PCA - phi true",200,-0.1,0.1);
  h_hgcal_sclusters_etaPCAMinusEtaTrue_corr = new TH1F("h_hgcal_sclusters_etaPCAMinusEtaTrue_corr","hgcal em supercluster eta position from PCA - eta true corrected",500,-0.1,0.1);
  h_hgcal_sclusters_phiPCAMinusPhiTrue_corr = new TH1F("h_hgcal_sclusters_phiPCAMinusPhiTrue_corr","hgcal em supercluster phi position from PCA - phi true corrected",500,-0.1,0.1);
  h_hgcal_sclusters_etaSCMinusEtaTrue = new TH1F("h_hgcal_sclusters_etaSCMinusEtaTrue","hgcal em supercluster eta position from PCA - eta true",200,-0.1,0.1);
  h_hgcal_sclusters_phiSCMinusPhiTrue = new TH1F("h_hgcal_sclusters_phiSCMinusPhiTrue","hgcal em supercluster phi position from PCA - phi true",200,-0.1,0.1);
  h_hgcal_sclusters_etaSCMinusEtaTrue_corr = new TH1F("h_hgcal_sclusters_etaSCMinusEtaTrue_corr","hgcal em supercluster eta position from PCA - eta true corrected",500,-0.1,0.1);
  h_hgcal_sclusters_phiSCMinusPhiTrue_corr = new TH1F("h_hgcal_sclusters_phiSCMinusPhiTrue_corr","hgcal em supercluster phi position from PCA - phi true corrected",500,-0.1,0.1);
  h_hgcal_sclusters_etaPCAMinusEtaSC = new TH1F("h_hgcal_sclusters_etaPCAMinusEtaSC","hgcal em supercluster eta position from PCA - eta sc position",200,-0.1,0.1);
  h_hgcal_sclusters_phiPCAMinusPhiSC = new TH1F("h_hgcal_sclusters_phiPCAMinusPhiSC","hgcal em supercluster phi position from PCA - phi sc position",200,-0.1,0.1);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrue = new TH1F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrue","hgcal em supercluster eta pdirection from PCA - eta ele direction",100,-0.5,0.5);
  h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue = new TH1F("thetadirPCAMinusThetaDirTrue","hgcal em supercluster theta pdirection from PCA - theeta ele direction",100,-0.1,0.1);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrue = new TH1F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrue","hgcal em supercluster phi position from PCA - phi ele direction",100,-1.0,1.0);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueUnsigned = new TH1F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueUnsigned","hgcal em supercluster phi position from PCA - phi ele direction, unsigned",100,-1.0,1.0);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueNeg = new TH1F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueNeg","hgcal em supercluster phi position from PCA - phi ele direction, phi<0.",100,-1.0,1.0);
  h_hgcal_sclusters_thetadirPCAMinusThetaDirTrueVsEnergy = new TH2F("thetadirPCAMinusThetaDirTrueVsEnergy","hgcal em supercluster theta pdirection from PCA - theta ele direction vs energy",100,100.,250.,100,-0.1,0.1);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy","hgcal em supercluster phi position from PCA - phi ele direction vs energy",100, 100.,250.,100,-1.0,1.0);

  h_hgcal_sclusters_etaSCMinusEtaTrueVsEta = new TH2F("h_hgcal_sclusters_etaSCMinusEtaTrueVsEta","hgcal em supercluster eta position from SC - eta true vs eta",200,-0.1,0.1,100,-3.0,3.0);
  h_hgcal_sclusters_etaSeedMinusEtaTrueVsEta = new TH2F("h_hgcal_sclusters_etaSeedMinusEtaTrueVsEta","hgcal em seedcluster eta position from SC - eta true vs eta",200,-0.1,0.1,100,-3.0,3.0);
  h_hgcal_sclusters_phiSCMinusPhiTrueVsPhi = new TH2F("h_hgcal_sclusters_phiSCMinusPhiTrueVsPhi","hgcal em supercluster phi position from SC - phi true vs phi",200,-0.1,0.1,60,-3.14,3.14);
  h_hgcal_sclusters_etaPCAMinusEtaTrueVsEta = new TH2F("h_hgcal_sclusters_etaPCAMinusEtaTrueVsEta","hgcal em supercluster eta position from PCA - eta true vs eta",200,-0.1,0.1,100,-3.0,3.0);
  h_hgcal_sclusters_phiPCAMinusPhiTrueVsEta = new TH2F("h_hgcal_sclusters_phiPCAMinusPhiTrueVsEta","hgcal em supercluster phi position from PCA - phi true vs eta",200,-0.1,0.1,100,-3.0,3.0);
  h_hgcal_sclusters_phiSCMinusPhiTrueVsEta = new TH2F("h_hgcal_sclusters_phiSCMinusPhiTrueVsEta","hgcal em supercluster phi position from SC - phi true vs eta",200,-0.1,0.1,100,-3.0,3.0);
  h_hgcal_sclusters_phiPCAMinusPhiTrueVsPhi = new TH2F("h_hgcal_sclusters_phiPCAMinusPhiTrueVsPhi","hgcal em supercluster phi position from PCA - phi true vs phi",200,-0.1,0.1,60,-3.14,3.14);

  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy","hgcal em supercluster phi position from PCA - phi ele direction vs energy",100,-1.0,1.0,100,0.,500.);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEta = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEta","hgcal em supercluster phi position from PCA - phi ele direction vs eta ",100,-1.0,1.0, 100,-3.0,3.0);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiMC = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiMC","hgcal em supercluster phi position from PCA - phi ele direction vs phi MC ",60,-3.14,3.14,60,-3.14,3.14);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiPCA = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiPCA","hgcal em supercluster phi position from PCA - phi ele direction vs phi PCA",60,-3.14,3.14,60,-3.14,3.14);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedmultiplicity = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedmultiplicity","hgcal em supercluster phi position from PCA - phi ele direction vs seed multiplicity",100,-1.,1.0,100,0.,1000.);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedfraction = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedfraction","hgcal em supercluster phi position from PCA - phi ele direction vs seed energy fraction",100,-1.,1.0,100,0.,1.);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSCmultiplicity = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSCmultiplicity","hgcal em supercluster phi position from PCA - phi ele direction vs SC multiplicity",100,-1.0,1.0,25,0.,25.);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsphiPCAMinusPhiTrue_corr = new TH2F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsphiPCAMinusPhiTrue_corr","hgcal em phi direction from PCA - phi ele direction vs phi position PCA - phi true",100,-1.0,1.0,100,-0.01,0.01);   
  
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEnergy = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEnergy","hgcal em supercluster eta position from PCA - eta ele direction vs energy",100,-1.0,1.0,100,0.,500.);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsPhi = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsPhi","hgcal em supercluster eta position from PCA - eta ele direction vs phi ",100,-1.0,1.0, 60,-3.14,3.14);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaMC = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaMC","hgcal em supercluster eta position from PCA - eta ele direction vs eta MC ",100,-1.0,1.0,100,-3.0,3.0);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaPCA = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaPCA","hgcal em supercluster eta position from PCA - eta ele direction vs eta PCA",100,-1.0,1.0,100,-3.0,3.0);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicity = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicity","hgcal em supercluster eta position from PCA - eta ele direction vs seed multiplicity",100,-1.,1.0,100,0.,1000.);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicityPosZ = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicityPosZ","hgcal em supercluster eta position from PCA - eta ele direction vs seed multiplicity",100,-1.,1.0,100,0.,1000.);
  h_hgcal_sclusters_AbsetadirPCAMinusEtaDirTrueVsSeedmultiplicity = new TH2F("h_hgcal_sclusters_AbsetadirPCAMinusEtaDirTrueVsSeedmultiplicity","hgcal em supercluster eta position from PCA - eta ele direction vs seed multiplicity, z>0",100,0.,1.0,100,0.,1000.);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedfraction = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedfraction","hgcal em supercluster eta position from PCA - eta ele direction vs seed multiplicity, z>0",100,-1.,1.0,100,0.,1.);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSCmultiplicity = new TH2F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSCmultiplicity","hgcal em supercluster eta position from PCA - eta ele direction vs SC multiplicity",100,-1.0,1.0,25,0.,25.);
  
  h_hgcal_sclusters_etaPCAMinusEtaTrue_golden = new TH1F("h_hgcal_sclusters_etaPCAMinusEtaTrue_golden","hgcal em supercluster eta position from PCA - eta true, golden",200,-0.1,0.1);
  h_hgcal_sclusters_phiPCAMinusPhiTrue_golden = new TH1F("h_hgcal_sclusters_phiPCAMinusPhiTrue_golden","hgcal em supercluster phi position from PCA - phi true, golden",200,-0.1,0.1);
  h_hgcal_sclusters_etaSCMinusEtaTrue_golden = new TH1F("h_hgcal_sclusters_etaSCMinusEtaTrue_golden","hgcal em supercluster eta position from PCA - eta true, golden",200,-0.1,0.1);
  h_hgcal_sclusters_phiSCMinusPhiTrue_golden = new TH1F("h_hgcal_sclusters_phiSCMinusPhiTrue_golden","hgcal em supercluster phi position from PCA - phi true, golden",200,-0.1,0.1);
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrue_golden = new TH1F("h_hgcal_sclusters_etadirPCAMinusEtaDirTrue_golden","hgcal em supercluster eta pdirection from PCA - eta ele direction, golden",100,-0.5,0.5);
  h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue_golden = new TH1F("thetadirPCAMinusThetaDirTrue_golden","hgcal em supercluster theta direction from PCA - theta ele direction, golden",100,-0.25,0.25);
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrue_golden = new TH1F("h_hgcal_sclusters_phidirPCAMinusPhiDirTrue_golden","hgcal em supercluster phi position from PCA - phi ele direction, golden",100,-1.0,1.0);
  
  h_hgcal_sclusters_subclustersVsEnergy = new TH2F("h_hgcal_sclusters_subclustersVsEnergy","hgcal em supercluster multiplicity vs energy",25,0.,25.,100,60.,350.);
  h_hgcal_sclusters_seedEnergyVsEnergy = new TH2F("h_hgcal_sclusters_seedEnergyVsEnergy","hgcal em seed cluster energy vs energy",100,60.,350.,100,60.,350.);
  h_hgcal_sclusters_seedEnergyVsEta= new TH2F("h_hgcal_sclusters_seedEnergyVsEta","hgcal em seed cluster energy vs eta",100,60.,350.,100,-3.0,3.0);

  
  // shower parametrization
  
  h_hgcal_param_meanTvsy = new TH1F("h_hgcal_param_meanTvsy","hgcal param mean T vs y",1000,20.,20000.);
  h_hgcal_param_meanAlphavsy = new TH1F("h_hgcal_param_meanAlphavsy","hgcal param mean alpha vs y",1000,20.,20000.);
  h_hgcal_param_meanLnTvsy = new TH1F("h_hgcal_param_meanLnTvsy","hgcal param mean lnT vs y",1000,20.,20000.);
  h_hgcal_param_meanLnAlphavsy = new TH1F("h_hgcal_param_meanLnAlphavsy","hgcal param mean lnAlpha vs y",1000,20.,20000.);
  h_hgcal_param_sigmaLnTvsy = new TH1F("h_hgcal_param_sigmaLnTvsy","hgcal param sigma lnT vs y",1000,20.,20000.);
  h_hgcal_param_sigmaLnAlphavsy = new TH1F("h_hgcal_param_sigmaLnAlphavsy","hgcal param sigma lnAlpha vs y",1000,20.,20000.);
  h_hgcal_param_corrAlphaTvsy = new TH1F("h_hgcal_param_corrAlphaTvsy","hgcal param corr alpha T vs y",1000,20.,20000.);
  h_hgcal_param_rC = new TH1F("h_hgcal_param_rC","hgcal param transverse rC vs tau, E=100 GeV",1000,0.,4.);
  h_hgcal_param_rT = new TH1F("h_hgcal_param_rT","hgcal param transverse rT vs tau, E=100 GeV",1000,0.,4.);
  h_hgcal_param_p = new TH1F("h_hgcal_param_p","hgcal param transverse p vs tau, E=100 GeV",1000,0.,4.);

  h_hgcal_sclusters_predictedTvsy = new TH2F("h_hgcal_sclusters_predictedTvsy","hgcal param shower max vs y",1000,20.,20000.,100,0.,10.);
  h_hgcal_sclusters_predictedLnTvsy = new TH2F("h_hgcal_sclusters_predictedLnTvsy","hgcal param shower max vs y",1000,20.,20000.,100,0.5,2.5);
  h_hgcal_sclusters_predictedLength = new TH2F("h_hgcal_sclusters_predictedLength","hgcal param shower length vs energy",100,5.,25.,100,0.,10.);
  h_hgcal_sclusters_predictedLength_fullrange = new TH2F("h_hgcal_sclusters_predictedLength_fullrange","hgcal param shower length vs energy",100,0.,25.,100,-10.,10.);
	
  // histogram measured shower quantities
  h_hgcal_sclusters_length = new TH2F("h_hgcal_sclusters_length","hgcal seed cluster shower length vs energy",40,5.,25.,40,4.,6.);
  h_hgcal_sclusters_energyVSlength = new TH2F("h_hgcal_sclusters_energyVSlength","hgcal seed cluster shower energy vs length",40,4.,6.,40,5.,25);
  h_hgcal_sclusters_length_fullrange = new TH2F("h_hgcal_sclusters_length_fullrange","hgcal seed cluster shower length vs energy",100,0.,25.,100,-10.,10.);
  h_hgcal_sclusters_energyVSlength_fullrange = new TH2F("h_hgcal_sclusters_energyVSlength_fullrange","hgcal seed cluster shower energy vs length",100,-10.,10.,100,0.,25.);
  h_hgcal_sclusters_energyVSlength_cut_fullrange = new TH2F("h_hgcal_sclusters_energyVSlength_cut_fullrange","hgcal seed cluster shower energy vs length, cut",100,-10.,10.,100,0.,25.);
  h_hgcal_sclusters_energyVSlength_hasfirstlayer_fullrange = new TH2F("h_hgcal_sclusters_energyVSlength_hasfirstlayer_fullrange","hgcal seed cluster shower energy vs length, first layer",100,-10.,10.,100,0.,25.);
  h_hgcal_sclusters_energyVSlength_hasfirstlayer_cut_fullrange = new TH2F("h_hgcal_sclusters_energyVSlength_hasfirstlayer_cut_fullrange","hgcal seed cluster shower energy vs length, first layer, cut",100,-10.,10.,100,0.,25.);
  h_hgcal_sclusters_entryVSlength_fullrange = new TH2F("h_hgcal_sclusters_entryVSlength_fullrange","hgcal seed cluster shower entry position vs length",350,320.,350.,100,0.,25.);
  h_hgcal_sclusters_entryVSlength_cut_fullrange = new TH2F("h_hgcal_sclusters_entryVSlength_cut_fullrange","hgcal seed cluster shower entry position vs length, cut",350,320.,350.,100,0.,25.);
  h_hgcal_sclusters_entry_cutsigmaeta_cuthadem = new TH1F("h_hgcal_sclusters_entry_cutsigmaeta_cuthadem","hgcal shower start position after cuts ",350,320.,350.);

  // for all clusters in HGCAL
  h_hgcal_allclusters_predictedTvsy = new TH2F("h_hgcal_allclusters_predictedTvsy","hgcal param shower max vs y",1000,20.,20000.,100,0.,10.);
  h_hgcal_allclusters_predictedLnTvsy = new TH2F("h_hgcal_allclusters_predictedLnTvsy","hgcal param shower max vs y",1000,20.,20000.,100,0.5,2.5);
  h_hgcal_allclusters_predictedLength = new TH2F("h_hgcal_allclusters_predictedLength","hgcal param shower length vs energy",100,5.,25.,100,0.,10.);
  h_hgcal_allclusters_predictedLength_fullrange = new TH2F("h_hgcal_allclusters_predictedLength_fullrange","hgcal param shower length vs energy",100,0.,25.,100,-10.,10.);
	
  // histogram measured shower quantities
  h_hgcal_allclusters_transverse_energy_em = new TH1F("h_hgcal_allclusters_transverse_energy","hgcal all cluster shower transverse energy",100,0.,10.);
  h_hgcal_allclusters_length_fullrange = new TH2F("h_hgcal_allclusters_length_fullrange","hgcal cluster shower length vs energy",100,0.,25.,100,-10.,10.);
  h_hgcal_allclusters_energyVSlength_fullrange = new TH2F("h_hgcal_allclusters_energyVSlength_fullrange","hgcal cluster shower energy vs length",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_energyVSlength_cut_fullrange = new TH2F("h_hgcal_allclusters_energyVSlength_cut_fullrange","hgcal cluster shower energy vs length, cut",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_entryVSlength_fullrange = new TH2F("h_hgcal_allclusters_entryVSlength_fullrange","hgcal cluster shower entry position vs length",350,320,350.,100,0.,25.);
  h_hgcal_allclusters_entryVSlength_cut_fullrange = new TH2F("h_hgcal_allclusters_entryVSlength_cut_fullrange","hgcal cluster shower entry position vs length, cut",350,320,350.,100,0.,25.);
  h_hgcal_allclusters_energyVSpredictedLength_fullrange = new TH2F("h_hgcal_allclusters_energyVSpredictedLength_fullrange","hgcal cluster shower energy vs PredictedLength",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_energyVSpredictedLength_cut_fullrange = new TH2F("h_hgcal_allclusters_energyVSpredictedLength_cut_fullrange","hgcal cluster shower energy vs PredictedLength, cut",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_entryVSpredictedLength_cut_fullrange = new TH2F("h_hgcal_allclusters_entryVSpredictedLength_cut_fullrange","hgcal cluster shower entry position vs PredictedLength, cut",350,320,350.,100,0.,25.);
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_fullrange = new TH2F("h_hgcal_allclusters_energyVSlength_hasfirstlayer_fullrange","hgcal cluster shower energy vs length, first layer",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange = new TH2F("h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange","hgcal cluster shower energy vs length, first layer, cut",100,-10.,10.,100,0.,25.);
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange_weighted = new TH2F("h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange_weighted","hgcal cluster shower energy vs length, first layer, cut, weighted",100,-10.,10.,100,0.,25.);
  
  // PU clusters
  h_hgcal_puclusters_energy = new TH1F("h_hgcal_pucluster_energy","hgcal pu cluster energy",100,0.,10.);
  h_hgcal_puclusters_energyVSeta = new TH2F("h_hgcal_puclusters_energyVSeta","hgcal pu cluster energy vs eta",100,0.0,10.,60,1.5,3.0);
  h_hgcal_puclusters_energy_eta16 = new TH1F("h_hgcal_pucluster_energy_eta16","hgcal pu cluster energy eta 1.6",100,0.,10.);
  h_hgcal_puclusters_energy_eta20 = new TH1F("h_hgcal_pucluster_energy_eta20","hgcal pu cluster energy eta 2.0",100,0.,10.);
  h_hgcal_puclusters_energy_eta25 = new TH1F("h_hgcal_pucluster_energy_eta25","hgcal pu cluster energy eta 2.5",100,0.,10.);
  h_hgcal_puclusters_energy_eta29 = new TH1F("h_hgcal_pucluster_energy_eta29","hgcal pu cluster energy eta 2.9",100,0.,10.);
  
  // transverse and PCA quantities
  h_hgcal_sclusters_energyVSmeanradius_fullrange = new TH2F("h_hgcal_sclusters_energyVSmeanradius_fullrange","hgcal seed cluster shower energy vs transverse mean",100,-10.,10.,100,0.,20.);
  h_hgcal_sclusters_entryVSmeanradius_fullrange = new TH2F("h_hgcal_sclusters_entryVSmeanradius_fullrange","hgcal seed cluster shower entry position vs transverse mean",350,320.,350.,100,0.,20.);
  h_hgcal_sclusters_energyVSsigmaradius_fullrange = new TH2F("h_hgcal_sclusters_energyVSsigmaradius_fullrange","hgcal seed cluster shower energy vs transverse sigma",100,-10.,10.,100,0.,5.);
  h_hgcal_sclusters_entryVSsigmaradius_fullrange = new TH2F("h_hgcal_sclusters_entryVSsigmaradius_fullrange","hgcal seed cluster shower entry position vs transverse sigma",350,320.,350.,100,0.,5.);
  h_hgcal_sclusters_energyVSlongwidth_fullrange = new TH2F("h_hgcal_sclusters_energyVSlongwidth_fullrange","hgcal seed cluster shower energy vs long eigenvalue PCA",100,-10.,10.,100,0.,1.);
  h_hgcal_sclusters_entryVSlongwidth_fullrange = new TH2F("h_hgcal_sclusters_entryVSlongwidth_fullrange","hgcal seed cluster shower entry position vs long eigenvalue PCA",350,320.,350.,100,0.,1.);
  h_hgcal_sclusters_energyVStranswidth_fullrange = new TH2F("h_hgcal_sclusters_energyVStrans_fullrange","hgcal seed cluster shower energy vs trans eigenvalue PCA",100,-10.,10.,100,0.,1.); 
  h_hgcal_sclusters_entryVStranswidth_fullrange = new TH2F("h_hgcal_sclusters_entryVStrans_fullrange","hgcal seed cluster shower entry position vs trans eigenvalue PCA",350,320.,350.,100,0.,1.);
  h_hgcal_sclusters_energyVSeigenratio_fullrange = new TH2F("h_hgcal_sclusters_energyVSratio_fullrange","hgcal seed cluster shower energy vs eigenvalue ratio PCA",100,-10.,10.,100,0.,10.); 
  h_hgcal_sclusters_entryVSeigenratio_fullrange = new TH2F("h_hgcal_sclusters_entryVSratio_fullrange","hgcal seed cluster shower entry position vs eigenvalue ratio PCA",350,320.,350.,100,0.,10.);
  h_hgcal_sclusters_energyVSeigenratio_cut_fullrange = new TH2F("h_hgcal_sclusters_energyVSratio_cut_fullrange","hgcal seed cluster shower energy vs eigenvalue ratio PCA, cut",100,-10.,10.,100,0.,10.); 
  h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_fullrange = new TH2F("h_hgcal_sclusters_energyVSratio_hasfirstlayer_fullrange","hgcal seed cluster shower energy vs eigenvalue ratio PCA, first layer",100,-10.,10.,100,0.,10.); 
  h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_cut_fullrange = new TH2F("h_hgcal_sclusters_energyVSratio_hasfirstlayer_cut_fullrange","hgcal seed cluster shower energy vs eigenvalue ratio PCA, first layer, cut",100,-10.,10.,100,0.,10.); 
  h_hgcal_sclusters_deta_shower = new TH1F("h_hgcal_sclusters_deta_shower","hgcal seed cluster shower eta hit - eta shower",100,0.,0.1);
  h_hgcal_sclusters_sigmaradius_fullrange = new TH1F("h_hgcal_sclusters_sigmaradius_fullrange","hgcal seed cluster shower transverse sigma",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradius_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradius_fullrange","hgcal seed cluster sigma rt",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradiusVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmatransverseradiusVSeta_fullrange","hgcal seed cluster sigma rt vs eta",100,0.,5.,60,1.5,3.);
  h_hgcal_sclusters_sigmatransverseradius_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradius_corr_fullrange","hgcal seed cluster sigma rt corr",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradiusaxis_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradiusaxis_fullrange","hgcal seed cluster sigma rt",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradiusaxismiddle_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradiusaxismiddle_fullrange","hgcal seed cluster sigma rt middle layers",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradiusaxismiddle_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradiusaxismiddle_corr_fullrange","hgcal seed cluster sigma rt middle layers",100,0.,5.);
  h_hgcal_sclusters_sigmatransverseradiusaxisVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmatransverseradiusaxisVSeta_fullrange","hgcal seed cluster sigma rt vs eta",100,0.,5.,60,1.5,3.);
  h_hgcal_sclusters_sigmatransverseradiusaxisVSphi_fullrange = new TH2F("h_hgcal_sclusters_sigmatransverseradiusaxisVSphi_fullrange","hgcal seed cluster sigma rt vs phi",100,0.,5.,60,-3.14,3.14);
  h_hgcal_sclusters_sigmatransverseradiusaxisVSlength_fullrange = new TH2F("h_hgcal_sclusters_sigmatransverseradiusaxisVSlength_fullrange","hgcal seed cluster sigma rt vs length",100,0.,5.,40,5.,25.);
  h_hgcal_sclusters_sigmatransverseradiusaxis_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmatransverseradiusaxis_corr_fullrange","hgcal seed cluster sigma rt corr",100,0.,5.);
  h_hgcal_sclusters_sigmaradiusnorm_fullrange = new TH1F("h_hgcal_sclusters_sigmaradiusnorm_fullrange","hgcal seed cluster shower transverse sigma",100,0.,5.);
  h_hgcal_sclusters_sigmaetanorm20_fullrange = new TH1F("h_hgcal_sclusters_sigmaetanorm20_fullrange","hgcal seed cluster shower transverse sigma20",100,0.,5.);
  h_hgcal_sclusters_sigmaetanorm50_fullrange = new TH1F("h_hgcal_sclusters_sigmaetanorm50_fullrange","hgcal seed cluster shower transverse sigma50",100,0.,5.);
  h_hgcal_sclusters_sigmaetanorm100_fullrange = new TH1F("h_hgcal_sclusters_sigmaetanorm100_fullrange","hgcal seed cluster shower transverse sigma100",100,0.,5.);
  h_hgcal_sclusters_sigmaeta20_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta20_fullrange","hgcal seed cluster shower transverse sigma20",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta50_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta50_fullrange","hgcal seed cluster shower transverse sigma50",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta100_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta100_fullrange","hgcal seed cluster shower transverse sigma100",100,0.,0.05);
  h_hgcal_sclusters_sigmaradiusVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaradiusVSeta_fullrange","hgcal seed cluster shower transverse sigma vs eta",100,0.,5.,60,1.5,3.);
  h_hgcal_sclusters_sigmaradiusVSphi_fullrange = new TH2F("h_hgcal_sclusters_sigmaradiusVSphi_fullrange","hgcal seed cluster shower transverse sigma vs eta",100,0.,5.,60,-3.14,3.14);
  h_hgcal_sclusters_sigmaeta_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta_fullrange","hgcal seed cluster shower eta sigma",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_pu_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta_pu_fullrange","hgcal seed cluster shower eta sigma",100,0.,0.05);
  h_hgcal_sclusters_sigmaetanorm_fullrange = new TH1F("h_hgcal_sclusters_sigmaetanorm_fullrange","hgcal seed cluster shower eta sigma",100,0.,0.05);
  h_hgcal_sclusters_sigmaetaVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaetaVSeta_fullrange","hgcal seed cluster shower sigma eta vs eta",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaeta_puVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaeta_puVSeta_fullrange","hgcal seed cluster shower sigma eta vs eta",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_etaVSsigmaeta_fullrange = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_fullrange","hgcal seed cluster shower sigma eta vs eta",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmatransverseradius_fullrange = new TH2F("h_hgcal_sclusters_etaVSsigmatransverseradius_fullrange","hgcal seed cluster shower sigma rt vs eta",60,1.5,3.,100,0.,5.);
  h_hgcal_sclusters_etaVSsigmatransverseradiusaxis_fullrange = new TH2F("h_hgcal_sclusters_etaVSsigmatransverseradiusaxis_fullrange","hgcal seed cluster shower sigma rt vs eta",60,1.5,3.,100,0.,5.);
  h_hgcal_sclusters_sigmaphi_fullrange = new TH1F("h_hgcal_sclusters_sigmaphi_fullrange","hgcal seed cluster shower phi sigma",100,0.,0.05);
  h_hgcal_sclusters_sigmaphinorm_fullrange = new TH1F("h_hgcal_sclusters_sigmaphinorm_fullrange","hgcal seed cluster shower phi sigma",100,0.,0.05);
  h_hgcal_sclusters_sigmaphiVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaphiVSeta_fullrange","hgcal seed cluster shower sigma phi vs eta",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaphiVSphi_fullrange = new TH2F("h_hgcal_sclusters_sigmaphiVSphi_fullrange","hgcal seed cluster shower sigma phi vs phi",100,0.,0.05,60,-3.14,3.14);
  h_hgcal_sclusters_sigmaeta_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta_corr_fullrange","hgcal seed cluster shower eta sigma corr",100,0.,0.025);
  h_hgcal_sclusters_sigmaeta_pu_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmaeta_pu_corr_fullrange","hgcal seed cluster shower eta sigma corr",100,0.,0.025);
  h_hgcal_sclusters_sigmaetaVSpt_pu_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaetaVSpt_pu_corr_fullrange","hgcal seed cluster shower eta sigma corr vs pt",100,0.,0.025,100,0.,50.);
  h_hgcal_sclusters_sigmaetaVSeta_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaetaVSeta_corr_fullrange","hgcal seed cluster shower sigma eta vs eta corr",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaeta_puVSeta_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaeta_puVSeta_corr_fullrange","hgcal seed cluster shower sigma eta vs eta corr",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaeta_puVSEt_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaeta_puVSEt_corr_fullrange","hgcal seed cluster shower sigma eta vs et, corr",100,0.,0.05,100,0.,50.);
  h_hgcal_sclusters_sigmaeta_puVSseedfraction_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaeta_puVSseedfraction_corr_fullrange","hgcal seed cluster shower sigma eta vs et, corr",100,0.,0.05,100,0.,1.);
  h_hgcal_sclusters_sigmaeta10_puVSEt_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaeta10_puVSEt_corr_fullrange","hgcal seed cluster shower sigma eta vs et, corr",100,0.,0.05,100,0.,50.);
  h_hgcal_sclusters_sigmaphi_corr_fullrange = new TH1F("h_hgcal_sclusters_sigmaphi_corr_fullrange","hgcal seed cluster shower eta sigma corr",100,0.,0.05);
  h_hgcal_sclusters_sigmaphicorrVSphiminusphitrue = new TH2F("h_hgcal_sclusters_sigmaphicorrVSphiminusphitrue","hgcal seed cluster shower phi sigma corr vs phi - phi true",100,0.,0.05, 100,-0.05,0.05 );
  h_hgcal_sclusters_sigmaphiVSeta_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaphiVSeta_corr_fullrange","hgcal seed cluster shower sigma eta vs eta corr",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaphiVSEt_corr_fullrange = new TH2F("h_hgcal_sclusters_sigmaphiVSEt_corr_fullrange","hgcal seed cluster shower sigma eta vs eta corr",100,0.,0.05,100,0.,50.);
  h_hgcal_sclusters_sigmaetaw_fullrange = new TH1F("h_hgcal_sclusters_sigmaetaw_fullrange","hgcal seed cluster shower eta sigma logweighted",100,0.,0.05);
  h_hgcal_sclusters_sigmaetawnorm_fullrange = new TH1F("h_hgcal_sclusters_sigmaetawnorm_fullrange","hgcal seed cluster shower eta sigma logweighted",100,0.,0.05);
  h_hgcal_sclusters_sigmaetawVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaetawVSeta_fullrange","hgcal seed cluster shower eta sigma logweighted vs eta",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaetawnormVSeta_fullrange = new TH2F("h_hgcal_sclusters_sigmaetawnormVSeta_fullrange","hgcal seed cluster shower eta sigma logweighted vs eta",100,0.,0.05,60,1.5,3.);
  h_hgcal_sclusters_sigmaetaw200_fullrange = new TH1F("h_hgcal_sclusters_sigmaetaw200_fullrange","hgcal seed cluster shower eta sigma logweighted, cut 200",100,0.,0.05);  
   
  h_hgcal_sclusters_sigmaetaVSlayer = new TH2F("h_hgcal_sclusters_sigmaetaVSlayer","hgcal seed cluster shower eta sigma logweighted per layer",30,1.,31.,100,0.,0.05);
  h_hgcal_sclusters_sigmaetaVSlayer_norm = new TH2F("h_hgcal_sclusters_sigmaetaVSlayer_norm","hgcal seed cluster shower eta sigma logweighted per layer",30,0.,4.,100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_1 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_1","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_2 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_2","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_3 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_3","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_4 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_4","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_5 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_5","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_6 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_6","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_7 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_7","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_8 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_8","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_9 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_9","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_10 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_10","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_11 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_11","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_12 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_12","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_13 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_13","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_14 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_14","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_15 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_15","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_16 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_16","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_17 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_17","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_18 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_18","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_19 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_19","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_20 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_20","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_21 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_21","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_22 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_22","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_23 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_23","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_24 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_24","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_25 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_25","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_26 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_26","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_27 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_27","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_28 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_28","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_29 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_29","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_30 = new TH1F("h_hgcal_sclusters_sigmaeta_layer_30","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaetaVSlayer_cutlayer9 = new TH2F("h_hgcal_sclusters_sigmaetaVSlayer_cutlayer9","hgcal seed cluster shower eta sigma logweighted per layer",30,1.,31.,100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_1 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_1","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_2 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_2","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_3 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_3","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_4 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_4","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_5 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_5","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_6 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_6","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_7 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_7","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_8 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_8","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_9 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_9","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_10 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_10","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_11 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_11","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_12 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_12","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_13 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_13","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_14 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_14","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_15 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_15","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_16 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_16","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_17 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_17","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_18 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_18","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_19 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_19","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_20 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_20","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_21 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_21","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_22 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_22","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_23 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_23","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_24 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_24","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_25 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_25","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_26 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_26","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_27 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_27","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_28 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_28","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_29 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_29","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_sigmaeta_cutlayer9_30 = new TH1F("h_hgcal_sclusters_sigmaeta_cutlayer9_layer_30","hgcal seed cluster shower eta sigma logweighted per layer",100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_1 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_1","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_2 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_2","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_3 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_3","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_4 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_4","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_5 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_5","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_6 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_6","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_7 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_7","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_8 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_8","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_9 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_9","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_10 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_10","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_11 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_11","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_12 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_12","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_13 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_13","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_14 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_14","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_15 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_15","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_16 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_16","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_17 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_17","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_18 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_18","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_19 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_19","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_20 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_20","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_21 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_21","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_22 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_22","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_23 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_23","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_24 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_24","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_25 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_25","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_26 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_26","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_27 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_27","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_28 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_28","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_29 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_29","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);
  h_hgcal_sclusters_etaVSsigmaeta_30 = new TH2F("h_hgcal_sclusters_etaVSsigmaeta_layer_30","hgcal seed cluster shower eta sigma logweighted per layer",60,1.5,3.,100,0.,0.05);

  h_hgcal_sclusters_entryposition = new TH1F("h_hgcal_sclusters_entryposition","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_1mip = new TH1F("h_hgcal_sclusters_entryposition_1mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_2mip = new TH1F("h_hgcal_sclusters_entryposition_2mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_4mip = new TH1F("h_hgcal_sclusters_entryposition_4mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_8mip = new TH1F("h_hgcal_sclusters_entryposition_8mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_12mip = new TH1F("h_hgcal_sclusters_entryposition_12mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_entryposition_16mip = new TH1F("h_hgcal_sclusters_entryposition_16mip","hgcal seed cluster shower start z position",300,320.,350.);
  h_hgcal_sclusters_longitudinal = new TH1F("h_hgcal_sclusters_longitudinal","hgcal seed cluster shower longitudinal profile",100,0.,30.);
  h_hgcal_sclusters_layer = new TH1F("h_hgcal_sclusters_layer","hgcal seed cluster shower layer profile",100,0.,30.);
  h_hgcal_sclusters_lengthCompatibility = new TH1F("h_hgcal_sclusters_lengthCompatibility","hgcal seed cluster shower lengthCompatibility",100,-10.,10.);
  h_hgcal_sclusters_transversal = new TH1F("h_hgcal_sclusters_transversal","hgcal seed cluster shower transversal profile",100,0.,10.);
  h_hgcal_sclusters_transversalaxis = new TH1F("h_hgcal_sclusters_transversalaxis","hgcal seed cluster shower transversal detector profile",100,0.,10.);
  h_hgcal_sclusters_transversalaxis_calib = new TH1F("h_hgcal_sclusters_transversalaxis_calib","hgcal seed cluster shower transversal detector profile",100,0.,10.);
  h_hgcal_sclusters_firsteigenvalue = new TH1F("h_hgcal_sclusters_firsteigenvalue","hgcal seed cluster shower first eigenvalue",100,0.,1.);
  h_hgcal_sclusters_firsteigenvalue_nm1 = new TH1F("h_hgcal_sclusters_firsteigenvalue_nm1","hgcal seed cluster shower first eigenvalue",100,0.,1.);
  h_hgcal_sclusters_secondeigenvalue = new TH1F("h_hgcal_sclusters_secondeigenvalue","hgcal seed cluster shower second eigenvalue",100,0.,1.);
  h_hgcal_sclusters_secondeigenvalue_nm1 = new TH1F("h_hgcal_sclusters_secondeigenvalue_nm1","hgcal seed cluster shower second eigenvalue",100,0.,1.);
  h_hgcal_sclusters_thirdeigenvalue = new TH1F("h_hgcal_sclusters_thirdeigenvalue","hgcal seed cluster shower third eigenvalue",100,0.,0.25);
  h_hgcal_sclusters_thirdeigenvalue_nm1 = new TH1F("h_hgcal_sclusters_thirdeigenvalue_nm1","hgcal seed cluster shower third eigenvalue nm1",100,0.,0.25);
  h_hgcal_sclusters_firstsigma = new TH1F("h_hgcal_sclusters_firstsigma","hgcal seed cluster shower first sigma",100,0.,5.);
  h_hgcal_sclusters_firstsigma_nm1 = new TH1F("h_hgcal_sclusters_firstsigma_nm1","hgcal seed cluster shower first sigma",100,0.,5.);
  h_hgcal_sclusters_secondsigma = new TH1F("h_hgcal_sclusters_secondsigma","hgcal seed cluster shower second sigma",100,0.,5.);
  h_hgcal_sclusters_secondsigma_nm1 = new TH1F("h_hgcal_sclusters_secondsigma_nm1","hgcal seed cluster shower second sigma",100,0.,5.);
  h_hgcal_sclusters_thirdsigma = new TH1F("h_hgcal_sclusters_thirdsigma","hgcal seed cluster shower third sigma",100,0.,10.);
  h_hgcal_sclusters_thirdsigma_nm1 = new TH1F("h_hgcal_sclusters_thirdsigma_nm1","hgcal seed cluster shower third sigma nm1",100,0.,10.);
  h_hgcal_sclusters_firsteigenvalueVSeta = new TH2F("h_hgcal_sclusters_firsteigenvalueVSeta","hgcal seed cluster shower first eigenvalue vs eta",100,0.,1., 60,1.5,3.0);
  h_hgcal_sclusters_secondeigenvalueVSeta = new TH2F("h_hgcal_sclusters_secondeigenvalueVSeta","hgcal seed cluster shower second eigenvalue vs eta",100,0.,1., 60,1.5,3.0);
  h_hgcal_sclusters_transeigenvalueVSeta = new TH2F("h_hgcal_sclusters_transeigenvalueVSeta","hgcal seed cluster shower transversal eigenvalue vs eta",100,0.,1., 60,1.5,3.0);
  h_hgcal_sclusters_eigenratioVSeta = new TH2F("h_hgcal_sclusters_eigenratioVSeta","hgcal seed cluster shower transversal eigenvalue ratio vs eta",100,0.,10., 60,1.5,3.0);
  h_hgcal_sclusters_firsteigenvalueVSsecond = new TH2F("h_hgcal_sclusters_firsteigenvalueVSfirsteigenvalueVSsecond","hgcal seed cluster shower first eigenvalue vs second",100,0.,1.,100,0.,1.);
  h_hgcal_sclusters_transversalVSlongitudinal = new TH2F("h_hgcal_sclusters_transversalVSlongitudinal","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinalscaled = new TH2F("h_hgcal_sclusters_transversalVSlongitudinalscaled","hgcal seed cluster shower transversal VS longitudinal profile scaled",100,0.,4., 100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinalscaledmeasured = new TH2F("h_hgcal_sclusters_transversalVSlongitudinalscaledmeasured","hgcal seed cluster shower transversal VS longitudinal profile scaled from measured length",100,0.,3., 100,0.,30.);
  h_hgcal_sclusters_transversalaxisVSlongitudinal = new TH2F("h_hgcal_sclusters_transversalaxisVSlongitudinal","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);
  h_hgcal_sclusters_transversalaxisVSlongitudinalscaled = new TH2F("h_hgcal_sclusters_transversalaxisVSlongitudinalscaled","hgcal seed cluster shower transversal VS longitudinal profile scaled",100,0.,4., 100,0.,30.);
  h_hgcal_sclusters_transversalaxisVSlongitudinalscaledmeasured = new TH2F("h_hgcal_sclusters_transversalaxisVSlongitudinalscaledmeasured","hgcal seed cluster shower transversal VS longitudinal profile scaled from measured length",100,0.,3., 100,0.,30.);
  h_hgcal_sclusters_transversalVSeta = new TH2F("h_hgcal_sclusters_transversalVSeta","hgcal eed cluster shower transversal VS MC eta",100,0.,30., 600,-3.,3.);
  h_hgcal_sclusters_longitudinal_cut20 = new TH1F("h_hgcal_sclusters_longitudinal_cut20","hgcal seed cluster shower longitudinal profile",100,0.,30.);
  h_hgcal_sclusters_transversal_cut20 = new TH1F("h_hgcal_sclusters_transversal_cut20","hgcal seed cluster shower transversal profile",100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinal_cut20 = new TH2F("h_hgcal_sclusters_transversalVSlongitudinal_cut20","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);
  h_hgcal_sclusters_longitudinal_cut50 = new TH1F("h_hgcal_sclusters_longitudinal_cut50","hgcal seed cluster shower longitudinal profile",100,0.,30.);
  h_hgcal_sclusters_transversal_cut50 = new TH1F("h_hgcal_sclusters_transversal_cut50","hgcal seed cluster shower transversal profile",100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinal_cut50 = new TH2F("h_hgcal_sclusters_transversalVSlongitudinal_cut50","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);
  h_hgcal_sclusters_longitudinal_cut100 = new TH1F("h_hgcal_sclusters_longitudinal_cut100","hgcal seed cluster shower longitudinal profile",100,0.,30.);
  h_hgcal_sclusters_transversal_cut100 = new TH1F("h_hgcal_sclusters_transversal_cut100","hgcal seed cluster shower transversal profile",100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinal_cut100 = new TH2F("h_hgcal_sclusters_transversalVSlongitudinal_cut100","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);
  h_hgcal_sclusters_longitudinal_cut200 = new TH1F("h_hgcal_sclusters_longitudinal_cut200","hgcal seed cluster shower longitudinal profile",100,0.,30.);
  h_hgcal_sclusters_transversal_cut200 = new TH1F("h_hgcal_sclusters_transversal_cut200","hgcal seed cluster shower transversal profile",100,0.,30.);
  h_hgcal_sclusters_transversalVSlongitudinal_cut200 = new TH2F("h_hgcal_sclusters_transversalVSlongitudinal_cut200","hgcal seed cluster shower transversal VS longitudinal profile",100,0.,30., 100,0.,30.);

  h_hgcal_sclusters_longitudinal_fit_chi2 = new TH1F("h_hgcal_sclusters_longitudinal_fit_chi2","hgcal seed cluster fit shower longitudinal profile chi2",200,0.,20.);
  h_hgcal_sclusters_longitudinal_fit_chi2_pnorm = new TH1F("h_hgcal_sclusters_longitudinal_fit_chi2_pnorm","hgcal seed cluster fit shower longitudinal profile chi2  normalization fixed to p",200,0.,20.);
  h_hgcal_sclusters_longitudinal_fit_chi2_bounded = new  TH1F("h_hgcal_sclusters_longitudinal_fit_chi2_bounded","hgcal seed cluster fit shower longitudinal profile chi2 bounded parameters",200,0.,20.);
  h_hgcal_sclusters_longitudinal_fit_chi2_nm1 = new TH1F("h_hgcal_sclusters_longitudinal_fit_chi2_nm1","hgcal seed cluster fit shower longitudinal profile chi2",200,0.,20.);
  h_hgcal_sclusters_longitudinal_fit_chi2VSeta = new TH2F("h_hgcal_sclusters_longitudinal_fit_chi2VSeta","hgcal seed cluster fit shower longitudinal profile chi2 vs eta",200,0.,20., 60,1.5,3.0);
  h_hgcal_sclusters_longitudinal_fit_alphaVSbeta = new TH2F("h_hgcal_sclusters_longitudinal_fit_alphaVSbeta","hgcal seed cluster fit shower longitudinal profile alpha vs beta",100,0.,15.,100,0.,1.);
  h_hgcal_sclusters_longitudinal_fit_alphaVSinvbeta = new TH2F("h_hgcal_sclusters_longitudinal_fit_alphaVSinvbeta","hgcal seed cluster fit shower longitudinal profile alpha vs invbeta",100,0.,15.,100,0.,4.);
  h_hgcal_sclusters_longitudinal_fit_alphaVSenergy = new TH2F("h_hgcal_sclusters_longitudinal_fit_alphaVSenergy","hgcal seed cluster fit shower longitudinal profile alpha vs energy",100,0.,10.,100,0.,15.);
  h_hgcal_sclusters_longitudinal_fit_betaVSenergy = new TH2F("h_hgcal_sclusters_longitudinal_fit_betaVSenergy","hgcal seed cluster fit shower longitudinal profile beta vs energy",100,0.,10.,100,0.,1.);
  h_hgcal_sclusters_longitudinal_fit_invbetaVSenergy = new TH2F("h_hgcal_sclusters_longitudinal_fit_invbetaVSenergy","hgcal seed cluster fit shower longitudinal profile invbeta vs energy",100,0.,10.,100,0.,4.);
  h_hgcal_sclusters_longitudinal_fit_leakage = new TH1F("h_hgcal_sclusters_longitudinal_fit_leakage","hgcal seed cluster fit shower longitudinal profile integrated leakage",200,0.,0.1);
  h_hgcal_sclusters_longitudinal_fit_leakageVShoverem = new TH2F("h_hgcal_sclusters_longitudinal_fit_leakageVShoverem","hgcal seed cluster fit shower longitudinal profile integrated leakage vs hovereem",200,0.,0.1,200,0.,0.1);
  h_hgcal_sclusters_longitudinal_fit_leakage_cutseedpos = new TH1F("h_hgcal_sclusters_longitudinal_fit_leakage_cutseedpos","hgcal seed cluster fit shower longitudinal profile integrated leakage cutseedpos",200,0.,0.1);
  h_hgcal_sclusters_longitudinal_fit_leakageVShoverem_cutseedpos = new TH2F("h_hgcal_sclusters_longitudinal_fit_leakageVShoverem_cutseedpos","hgcal seed cluster fit shower longitudinal profile integrated leakage vs hovereem cutseedpos",200,0.,0.1,200,0.,0.1);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist","hgcal seed cluster shower longitudinal profile kolmogorov dist",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem","hgcal seed cluster shower longitudinal profile kolmogorov dist cutsigmaeta_cuthadem",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem_cutseedpos = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem_cutseedpos","hgcal seed cluster shower longitudinal profile kolmogorov dist cutsigmaeta_cuthadem_cutseedpos",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_nm1 = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_nm1","hgcal seed cluster shower longitudinal profile kolmogorov dist N-1",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_best","hgcal seed cluster shower longitudinal profile kolmogorov dist",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem","hgcal seed cluster shower longitudinal profile kolmogorov dist cutsigmaeta_cuthadem",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem_cutseedpos = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem_cutseedpos","hgcal seed cluster shower longitudinal profile kolmogorov dist cutsigmaeta_cuthadem_cutseedpos",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_nm1 = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_nm1","hgcal seed cluster shower longitudinal profile kolmogorov dist n-1",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_prob = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_prob","hgcal seed cluster shower longitudinal profile kolmogorov proba",200,0.,1.);
  h_hgcal_sclusters_longitudinal_kolmogorov_prob_best = new TH1F("h_hgcal_sclusters_longitudinal_kolmogorov_prob_best","hgcal seed cluster shower longitudinal profile kolmogorov proba",200,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_prob = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_prob","hgcal seed cluster shower longitudinal profile kolmogorov proba transversal profile",400,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_dist = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_dist","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal profile",400,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_prob_first = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_prob_first","hgcal seed cluster shower longitudinal profile kolmogorov proba transversal profile_first",400,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_dist_first = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_dist_first","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal profile_first",400,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_prob_last = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_prob_last","hgcal seed cluster shower longitudinal profile kolmogorov proba transversal profile_last",400,0.,1.);
  h_hgcal_sclusters_transversal_kolmogorov_dist_last = new TH1F("h_hgcal_sclusters_transversal_kolmogorov_dist_last","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal profile_last",400,0.,1.);
  h_hgcal_sclusters_kolmogorov_prob_3D = new TH1F("h_hgcal_sclusters_kolmogorov_prob_3D","hgcal seed cluster shower longitudinal profile kolmogorov proba",200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_3D = new TH1F("h_hgcal_sclusters_kolmogorov_dist_3D","hgcal seed cluster shower longitudinal profile kolmogorov dist",200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_3D_px = new TH1F("h_hgcal_sclusters_kolmogorov_dist_3D_px","hgcal seed cluster shower longitudinal profile kolmogorov dist",200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_3D_py = new TH1F("h_hgcal_sclusters_kolmogorov_dist_3D_py","hgcal seed cluster shower longitudinal profile kolmogorov dist",200,0.,1.);
  h_hgcal_sclusters_longitudinal_fit_chi2VSseedfraction = new TH2F("h_hgcal_sclusters_longitudinal_fit_chi2VSseedfraction","hgcal seed cluster fit shower longitudinal profile chi2 vs seed energy fraction",100,0.,20., 100,0.,1.);
  h_hgcal_sclusters_longitudinal_fit_chi2VSdphidir = new TH2F("h_hgcal_sclusters_longitudinal_fit_chi2VSdphidir","hgcal seed cluster fit shower longitudinal profile chi2 vs phi direction",100,0.,20., 100,-1.,1.);
  h_hgcal_sclusters_longitudinal_fit_chi2VSeoveretrue = new TH2F("h_hgcal_sclusters_longitudinal_fit_eoveretrue","hgcal seed cluster fit shower longitudinal profile chi2 vs phi direction",100,0.,20., 130,0.,1.3);
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal = new TH2F("h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal VS longitudinal",200,0.,1.,200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_nm1 = new TH2F("h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_nm1","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal VS longitudinal",200,0.,1.,200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cutsigmaeta_cuthadem = new TH2F("h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cutsigmaeta_cuthadem","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal VS longitudinal cutsigmaeta_cuthadem",200,0.,1.,200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem = new TH2F("h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal VS longitudinal cuthadem",200,0.,1.,200,0.,1.);
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem_cutseedpos = new TH2F("h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem_cutseedpos","hgcal seed cluster shower longitudinal profile kolmogorov dist transversal VS longitudinal cuthadem_cutseedpos",200,0.,1.,200,0.,1.);

  h_hgcal_scclusters_eoveretrue_em_cutpos = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutpos","hgcal em supercluster E/Etrue cutpos",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength","hgcal em supercluster E/Etrue cutpos_cutlength",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta","hgcal em supercluster E/Etrue cutpos_cutlength_cutsigmaeta",130,0.0,1.3);

  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta","hgcal em supercluster E/Etrue cutsigmaeta",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem","hgcal em supercluster E/Etrue cutsigmaeta_cuthadem",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutkolmogorov = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutkolmogorov","hgcal em supercluster E/Etrue cutsigmaeta_cuthadem_cutkolmogorov",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos","hgcal em supercluster E/Etrue cutsigmaeta_cuthadem_cutpos",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength","hgcal em supercluster E/Etrue cutsigmaeta_cuthadem_cutpos_cutlength",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength_cutkolmogorov = new TH1F("h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength_cutkolmogorov","hgcal em supercluster E/Etrue cutsigmaeta_cutpos_cutlength_cutkolmogorov",130,0.0,1.3);
  
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos","hgcal em supercluster E/Etrue cutpos",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength","hgcal em supercluster E/Etrue cutpos_cutlength",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength_cutsigmaeta = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength_cutsigmaeta","hgcal em supercluster E/Etrue cutpos_cutlength_cutsigmaeta",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos","hgcal em supercluster E/Etrue cutpos",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength","hgcal em supercluster E/Etrue cutpos_cutlength",130,0.0,1.3);
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength_cutsigmaeta = new TH1F("h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength_cutsigmaeta","hgcal em supercluster E/Etrue cutpos_cutlength_cutsigmaeta",130,0.0,1.3);

  h_hgcal_sclusters_hoverem = new TH1F("h_hgcal_sclusters_hoverem","hgcal hadronic over em",100,0.0,0.5);
  h_hgcal_sclusters_hoverem_cone01 = new TH1F("h_hgcal_sclusters_hoverem_cone01","hgcal hadronic over em",100,0.0,0.5);
  h_hgcal_sclusters_hoverem_cutsigmaeta = new TH1F("h_hgcal_sclusters_hoverem_cutsigmaeta","hgcal hadronic over em_cutsigmaeta",1000,0.0,0.5);
  h_hgcal_sclusters_hoveremVSeta = new TH2F("h_hgcal_sclusters_hoveremVSeta","hgcal hadronic over em vs eta",100,0.0,0.5,60,1.5,3.0);
  h_hgcal_sclusters_hoveremVSphi = new TH2F("h_hgcal_sclusters_hoveremVSphi","hgcal hadronic over em vs phi",100,0.0,0.5,157,-3.14,3.14);
  h_hgcal_sclusters_entryVShoverem = new TH2F("h_hgcal_sclusters_entryVShoverem","hgcal entry pos vs hadronic over em",350,320.,350.,100,0.0,0.5);
  h_hgcal_sclusters_expectedlengthVShoverem = new TH2F("h_hgcal_sclusters_expectedlengthVShoverem","hgcal expectedlength vs hadronic over em",100,0.,25.,100,0.0,0.5);
  h_hgcal_sclusters_hoverem1 = new TH1F("h_hgcal_sclusters_hoverem1","hgcal hadronic over em1",100,0.0,0.5);
  h_hgcal_sclusters_hoverem1VSeta = new TH2F("h_hgcal_sclusters_hoverem1VSeta","hgcal hadronic over em1 vs eta",100,0.0,0.5,60,1.5,3.0);
  h_hgcal_sclusters_entryVShoverem1 = new TH2F("h_hgcal_sclusters_entryVShoverem1","hgcal entry pos vs hadronic over em1",350,320.,350.,100,0.0,0.5);
  h_hgcal_sclusters_expectedlengthVShoverem1 = new TH2F("h_hgcal_sclusters_expectedlengthVShoverem1","hgcal expectedlength vs hadronic over em1",100,0.,25.,100,0.0,0.5);
  h_hgcal_sclusters_hoverem2 = new TH1F("h_hgcal_sclusters_hoverem2","hgcal hadronic over em2",100,0.0,0.5);
  h_hgcal_sclusters_hoverem2VSeta = new TH2F("h_hgcal_sclusters_hoverem2VSeta","hgcal hadronic over em2 vs eta",100,0.0,0.5,60,1.5,3.0);
  h_hgcal_sclusters_entryVShoverem2 = new TH2F("h_hgcal_sclusters_entryVShoverem2","hgcal entry pos vs hadronic over em2",350,320.,350.,100,0.0,0.5);
  h_hgcal_sclusters_expectedlengthVShoverem2 = new TH2F("h_hgcal_sclusters_expectedlengthVShoverem2","hgcal expectedlength vs hadronic over em2",100,0.,25.,100,0.0,0.5);

  h_hgcal_sclusters_e2530OverEtot = new TH1F("h_hgcal_sclusters_e2530OverEtot","hgcal seed cluster longitudinal energy ratio",100,0.,1.);
  h_hgcal_sclusters_e0110OverEtot = new TH1F("h_hgcal_sclusters_e0110OverEtot","hgcal seed cluster longitudinal energy ratio",100,0.,1.);
  h_hgcal_sclusters_e0110OverEtotVSe2530OverEtot = new TH2F("h_hgcal_sclusters_e0110OverEtotVSe2530OverEtot","hgcal seed cluster longitudinal energy ratio",100,0.,1.,100,0.,1.);

  double y[9] = {50., 100., 200., 500., 1000., 2000., 5000., 10000., 19500.};
  double tau[9] = {0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, 3.9};
  // Sam should not be needed here but bHom=0 doesmn't seem to be taken into accout by
  // EMShowerParametrization ...
  for (int i=0; i<9; i++) {
    h_hgcal_param_meanTvsy->Fill(y[i],showerparam_->meanTSam(log(y[i])));
    h_hgcal_param_meanAlphavsy->Fill(y[i],showerparam_->meanAlphaSam(log(y[i])));
    h_hgcal_param_meanLnTvsy->Fill(y[i],showerparam_->meanLnTSam(log(y[i])));
    h_hgcal_param_meanLnAlphavsy->Fill(y[i],showerparam_->meanLnAlphaSam(log(y[i])));
    h_hgcal_param_sigmaLnTvsy->Fill(y[i],showerparam_->sigmaLnTSam(log(y[i])));
    h_hgcal_param_sigmaLnAlphavsy->Fill(y[i],showerparam_->sigmaLnAlphaSam(log(y[i])));
    h_hgcal_param_corrAlphaTvsy->Fill(y[i],showerparam_->correlationAlphaTSam(log(y[i])));
    h_hgcal_param_rC->Fill(tau[i],showerparam_->rCSam(tau[i],100.));
    h_hgcal_param_rT->Fill(tau[i],showerparam_->rTSam(tau[i],100.));
    h_hgcal_param_p->Fill(tau[i],showerparam_->pSam(tau[i],100.));
  }  
  

  // histos titles
  h_mcNum              -> GetXaxis()-> SetTitle("N_{gen}");
  h_mcNum              -> GetYaxis()-> SetTitle("Events");
  h_simEta             -> GetXaxis()-> SetTitle("#eta");
  h_simEta             -> GetYaxis()-> SetTitle("Events");
  h_simP               -> GetXaxis()-> SetTitle("p (GeV/c)");
  h_simP               -> GetYaxis()-> SetTitle("Events");

}

void
HGCALElectronClusterAnalyzer::endJob(){

  histfile_->cd();

  // mc truth
  h_mcNum->Write();

  // mc
  h_simEta->Write();
  h_simAbsEta->Write();
  h_simP->Write();
  h_simPt->Write();
  h_simZ->Write();
  h_simPhi->Write();
  h_simPtEta->Write();

  h_simZ_all->Write();
  h_simZ_electrons->Write();

  // HGCAL clusters histos

  h_hgcal_foundClusters_em              -> Write();
  h_hgcal_foundClustersVSeta_em-> Write();
  h_hgcal_foundClustersVSetaEt5_em-> Write();
  //h_hgcal_clusterEnergy_em              -> Write();
  //h_hgcal_clusterEtaVsPhi_em              -> Write();
  //h_hgcal_foundClusters_had              -> Write();
  //h_hgcal_clusterEnergy_had              -> Write();
  //h_hgcal_clusterEtaVsPhi_had              -> Write();

  h_hgcal_sclusters_energy_em->Write();
  h_hgcal_sclusters_transverse_energy_em->Write();
  h_hgcal_sclusters_energy_pos_em->Write();
  h_hgcal_sclusters_energy_neg_em->Write();
  h_hgcal_sclusters_energyVSeta_em->Write();
  h_hgcal_sclusters_position_em->Write();
  h_hgcal_sclusters_etawidth_em->Write();
  h_hgcal_sclusters_phiwidth_em->Write();
  h_hgcal_sclusters_multiplicity_em->Write();
  h_hgcal_sclusters_newmultiplicity_em->Write();
  h_hgcal_sclusters_multiplicityVSeta_em->Write();
  h_hgcal_sclusters_newmultiplicityVSeta_em->Write();
  h_hgcal_sclusters_seedenergy_em->Write();
  h_hgcal_sclusters_newseedenergy_em->Write();
  h_hgcal_sclusters_seedfractions_em->Write();
  h_hgcal_sclusters_seedfractionVSeta_em->Write();
  h_hgcal_clusters_energy_em->Write();
  h_hgcal_clusters_position_em->Write();
  h_hgcal_clusters_multiplicity_em->Write();
  h_hgcal_clusters_multiplicityVSeta_em->Write();
  h_hgcal_clusters_rechitenergy_em->Write();
  h_hgcal_clusters_rechitenergy_12000_em->Write();

  h_hgcal_scclusters_eoveretrue_em->Write();
  h_hgcal_scclusters_eoveretrueVSeta_em->Write();
  h_hgcal_scclusters_eoveretrue_cut00_em->Write();
  h_hgcal_scclusters_eoveretrue_cut04_em->Write();
  h_hgcal_scclusters_eoveretrue_cut1_em->Write();
  h_hgcal_scclusters_eoveretrue_cut2_em->Write();
  h_hgcal_scclusters_eoveretrue_cut4_em->Write();
  h_hgcal_scclusters_eoveretrue_cut10_em->Write();
  h_hgcal_scclusters_eoveretrue_cut20_em->Write();
  h_hgcal_scclusters_eoveretrue_golden_em->Write();
  h_hgcal_scclusters_eoveretrue_noisecut_em->Write();
  h_hgcal_scclusters_eoveretrue_noisecut4_em->Write();
  h_hgcal_scclusters_ptoverpttrue_em->Write();
  h_hgcal_scclusters_detadphisubclusters_em->Write();
  h_hgcal_scclusters_detadphisubclusters_zoom_em->Write();
  h_hgcal_scclusters_detadphisubclusters_zoom_etagt2_em->Write();
  h_hgcal_scclusters_detadphisubclusters_zoom_etalt2_em->Write();
  h_hgcal_scclusters_detadphisubclusters_weighted_em->Write();
  h_hgcal_scclusters_detadphisubclusters_zoom_weighted_em->Write();

  h_hgcal_allhits_em[0]->Write();
  h_hgcal_allhits_em[1]->Write();
  h_hgcal_allhits_em[2]->Write();
  h_hgcal_allhits_em[3]->Write();
  h_hgcal_allhits_em[4]->Write();
  h_hgcal_allhits_em[5]->Write();
  h_hgcal_allhits_em[6]->Write();
  h_hgcal_allhits_em[7]->Write();
  h_hgcal_allhits_em[8]->Write();
  h_hgcal_allhits_em[9]->Write();
  
  h_hgcal_allhits_weighted_em[0]->Write();
  h_hgcal_allhits_weighted_em[1]->Write();
  h_hgcal_allhits_weighted_em[2]->Write();
  h_hgcal_allhits_weighted_em[3]->Write();
  h_hgcal_allhits_weighted_em[4]->Write();
  h_hgcal_allhits_weighted_em[5]->Write();
  h_hgcal_allhits_weighted_em[6]->Write();
  h_hgcal_allhits_weighted_em[7]->Write();
  h_hgcal_allhits_weighted_em[8]->Write();
  h_hgcal_allhits_weighted_em[9]->Write();
  
  h_hgcal_clustershits_em[0]->Write();
  h_hgcal_clustershits_em[1]->Write();
  h_hgcal_clustershits_em[2]->Write();
  h_hgcal_clustershits_em[3]->Write();
  h_hgcal_clustershits_em[4]->Write();
  h_hgcal_clustershits_em[5]->Write();
  h_hgcal_clustershits_em[6]->Write();
  h_hgcal_clustershits_em[7]->Write();
  h_hgcal_clustershits_em[8]->Write();
  h_hgcal_clustershits_em[9]->Write();
  
  h_hgcal_clustershits_weighted_em[0]->Write();
  h_hgcal_clustershits_weighted_em[1]->Write();
  h_hgcal_clustershits_weighted_em[2]->Write();
  h_hgcal_clustershits_weighted_em[3]->Write();
  h_hgcal_clustershits_weighted_em[4]->Write();
  h_hgcal_clustershits_weighted_em[5]->Write();
  h_hgcal_clustershits_weighted_em[6]->Write();
  h_hgcal_clustershits_weighted_em[7]->Write();
  h_hgcal_clustershits_weighted_em[8]->Write();
  h_hgcal_clustershits_weighted_em[9]->Write();
  
  h_hgcal_clustershits_cutlayers_em[0]->Write();
  h_hgcal_clustershits_cutlayers_em[1]->Write();
  h_hgcal_clustershits_cutlayers_em[2]->Write();
  h_hgcal_clustershits_cutlayers_em[3]->Write();
  h_hgcal_clustershits_cutlayers_em[4]->Write();
  h_hgcal_clustershits_cutlayers_em[5]->Write();
  h_hgcal_clustershits_cutlayers_em[6]->Write();
  h_hgcal_clustershits_cutlayers_em[7]->Write();
  h_hgcal_clustershits_cutlayers_em[8]->Write();
  h_hgcal_clustershits_cutlayers_em[9]->Write();
  
  h_hgcal_clustershits_cutlayers_weighted_em[0]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[1]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[2]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[3]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[4]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[5]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[6]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[7]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[8]->Write();
  h_hgcal_clustershits_cutlayers_weighted_em[9]->Write();
  
  h_hgcal_clustershits_cutlayers_cutlength_em[0]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[1]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[2]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[3]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[4]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[5]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[6]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[7]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[8]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_em[9]->Write();
  
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[0]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[1]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[2]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[3]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[4]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[5]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[6]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[7]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[8]->Write();
  h_hgcal_clustershits_cutlayers_cutlength_weighted_em[9]->Write();
  
  h_hgcal_superclustershits_em[0]->Write();
  h_hgcal_superclustershits_em[1]->Write();
  h_hgcal_superclustershits_em[2]->Write();
  h_hgcal_superclustershits_em[3]->Write();
  h_hgcal_superclustershits_em[4]->Write();
  h_hgcal_superclustershits_em[5]->Write();
  h_hgcal_superclustershits_em[6]->Write();
  h_hgcal_superclustershits_em[7]->Write();
  h_hgcal_superclustershits_em[8]->Write();
  h_hgcal_superclustershits_em[9]->Write();
  
  h_hgcal_superclustershits_weighted_em[0]->Write();
  h_hgcal_superclustershits_weighted_em[1]->Write();
  h_hgcal_superclustershits_weighted_em[2]->Write();
  h_hgcal_superclustershits_weighted_em[3]->Write();
  h_hgcal_superclustershits_weighted_em[4]->Write();
  h_hgcal_superclustershits_weighted_em[5]->Write();
  h_hgcal_superclustershits_weighted_em[6]->Write();
  h_hgcal_superclustershits_weighted_em[7]->Write();
  h_hgcal_superclustershits_weighted_em[8]->Write();
  h_hgcal_superclustershits_weighted_em[9]->Write();
  
  h_hgcal_superclustershits_cutlayers_em[0]->Write();
  h_hgcal_superclustershits_cutlayers_em[1]->Write();
  h_hgcal_superclustershits_cutlayers_em[2]->Write();
  h_hgcal_superclustershits_cutlayers_em[3]->Write();
  h_hgcal_superclustershits_cutlayers_em[4]->Write();
  h_hgcal_superclustershits_cutlayers_em[5]->Write();
  h_hgcal_superclustershits_cutlayers_em[6]->Write();
  h_hgcal_superclustershits_cutlayers_em[7]->Write();
  h_hgcal_superclustershits_cutlayers_em[8]->Write();
  h_hgcal_superclustershits_cutlayers_em[9]->Write();
  
  h_hgcal_superclustershits_cutlayers_weighted_em[0]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[1]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[2]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[3]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[4]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[5]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[6]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[7]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[8]->Write();
  h_hgcal_superclustershits_cutlayers_weighted_em[9]->Write();
  
  h_hgcal_superclustershits_cutlayers_cutlength_em[0]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[1]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[2]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[3]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[4]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[5]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[6]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[7]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[8]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_em[9]->Write();
  
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[0]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[1]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[2]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[3]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[4]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[5]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[6]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[7]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[8]->Write();
  h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[9]->Write();

  h_hgcal_shower_seed_em[0]->Write();
  h_hgcal_shower_seed_em[1]->Write();
  h_hgcal_shower_seed_em[2]->Write();
  h_hgcal_shower_seed_em[3]->Write();
  h_hgcal_shower_seed_em[4]->Write();
  h_hgcal_shower_seed_em[5]->Write();
  h_hgcal_shower_seed_em[6]->Write();
  h_hgcal_shower_seed_em[7]->Write();
  h_hgcal_shower_seed_em[8]->Write();
  h_hgcal_shower_seed_em[9]->Write();
  
  h_hgcal_shower_sc_em[0]->Write();
  h_hgcal_shower_sc_em[1]->Write();
  h_hgcal_shower_sc_em[2]->Write();
  h_hgcal_shower_sc_em[3]->Write();
  h_hgcal_shower_sc_em[4]->Write();
  h_hgcal_shower_sc_em[5]->Write();
  h_hgcal_shower_sc_em[6]->Write();
  h_hgcal_shower_sc_em[7]->Write();
  h_hgcal_shower_sc_em[8]->Write();
  h_hgcal_shower_sc_em[9]->Write();
  
  h_hgcal_shower_sc[0]->Write();
  h_hgcal_shower_sc[1]->Write();
  h_hgcal_shower_sc[2]->Write();
  h_hgcal_shower_sc[3]->Write();
  h_hgcal_shower_sc[4]->Write();
  h_hgcal_shower_sc[5]->Write();
  h_hgcal_shower_sc[6]->Write();
  h_hgcal_shower_sc[7]->Write();
  h_hgcal_shower_sc[8]->Write();
  h_hgcal_shower_sc[9]->Write();
  
  h_hgcal_shower_seed_rotated_em[0]->Write();
  h_hgcal_shower_seed_rotated_em[1]->Write();
  h_hgcal_shower_seed_rotated_em[2]->Write();
  h_hgcal_shower_seed_rotated_em[3]->Write();
  h_hgcal_shower_seed_rotated_em[4]->Write();
  h_hgcal_shower_seed_rotated_em[5]->Write();
  h_hgcal_shower_seed_rotated_em[6]->Write();
  h_hgcal_shower_seed_rotated_em[7]->Write();
  h_hgcal_shower_seed_rotated_em[8]->Write();
  h_hgcal_shower_seed_rotated_em[9]->Write();
 
  h_hgcal_sclusters_etaPCAMinusEtaTrue->Write();
  h_hgcal_sclusters_phiPCAMinusPhiTrue->Write();
  h_hgcal_sclusters_etaPCAMinusEtaTrue_corr->Write();
  h_hgcal_sclusters_phiPCAMinusPhiTrue_corr->Write();
  h_hgcal_sclusters_etaPCAMinusEtaSC->Write();
  h_hgcal_sclusters_phiPCAMinusPhiSC->Write();
  h_hgcal_sclusters_etaSCMinusEtaTrue->Write();
  h_hgcal_sclusters_phiSCMinusPhiTrue->Write();
  h_hgcal_sclusters_etaSCMinusEtaTrue_corr->Write();
  h_hgcal_sclusters_phiSCMinusPhiTrue_corr->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrue->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrue->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueUnsigned->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueNeg->Write();
  h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue->Write();
  // pb here cannot write those histos otherwise root crashes when closing the files
  //h_hgcal_sclusters_phidirPCAMinusPhiDirVsEnergy->Write();
  //h_hgcal_sclusters_thetadirPCAMinusThetaDirTrueVsEnergy->Write();
  h_hgcal_sclusters_etaSCMinusEtaTrueVsEta->Write();  
  h_hgcal_sclusters_etaSeedMinusEtaTrueVsEta->Write();  
  h_hgcal_sclusters_phiSCMinusPhiTrueVsPhi->Write();  
  h_hgcal_sclusters_etaPCAMinusEtaTrueVsEta->Write();  
  h_hgcal_sclusters_phiPCAMinusPhiTrueVsEta->Write();
  h_hgcal_sclusters_phiSCMinusPhiTrueVsEta->Write();
  h_hgcal_sclusters_phiPCAMinusPhiTrueVsPhi->Write();
  
  h_hgcal_sclusters_etaPCAMinusEtaTrue_golden->Write();
  h_hgcal_sclusters_phiPCAMinusPhiTrue_golden->Write();
  h_hgcal_sclusters_etaSCMinusEtaTrue_golden->Write();
  h_hgcal_sclusters_phiSCMinusPhiTrue_golden->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrue_golden->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrue_golden->Write();
  h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue_golden->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEta->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiMC->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiPCA->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedmultiplicity->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedfraction->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSCmultiplicity->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEnergy->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsPhi->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaMC->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaPCA->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicity->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicityPosZ->Write();
  h_hgcal_sclusters_AbsetadirPCAMinusEtaDirTrueVsSeedmultiplicity->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedfraction->Write();
  h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSCmultiplicity->Write();
  h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsphiPCAMinusPhiTrue_corr->Write();   
  h_hgcal_sclusters_subclustersVsEnergy->Write();
  h_hgcal_sclusters_seedEnergyVsEnergy->Write();
  h_hgcal_sclusters_seedEnergyVsEta->Write();

  h_hgcal_param_meanTvsy->Write();
  h_hgcal_param_meanAlphavsy->Write();
  h_hgcal_param_meanLnTvsy->Write();
  h_hgcal_param_meanLnAlphavsy->Write();
  h_hgcal_param_sigmaLnTvsy->Write();
  h_hgcal_param_sigmaLnAlphavsy->Write();
  h_hgcal_param_corrAlphaTvsy->Write();
  h_hgcal_param_rC->Write();
  h_hgcal_param_rT->Write();
  h_hgcal_param_p->Write();

  h_hgcal_sclusters_predictedTvsy ->Write(); 
  h_hgcal_sclusters_predictedLnTvsy->Write();
  h_hgcal_sclusters_predictedLength->Write();
  h_hgcal_sclusters_predictedLength_fullrange->Write();
  	
  h_hgcal_sclusters_length->Write();
  h_hgcal_sclusters_energyVSlength->Write();
  h_hgcal_sclusters_length_fullrange->Write();
  h_hgcal_sclusters_energyVSlength_fullrange->Write();
  h_hgcal_sclusters_energyVSlength_cut_fullrange->Write();
  h_hgcal_sclusters_energyVSlength_hasfirstlayer_fullrange->Write();
  h_hgcal_sclusters_energyVSlength_hasfirstlayer_cut_fullrange->Write();
  h_hgcal_sclusters_entryVSlength_fullrange->Write();
  h_hgcal_sclusters_entryposition->Write();
  h_hgcal_sclusters_entryposition_1mip->Write();
  h_hgcal_sclusters_entryposition_2mip->Write();
  h_hgcal_sclusters_entryposition_4mip->Write();
  h_hgcal_sclusters_entryposition_8mip->Write();
  h_hgcal_sclusters_entryposition_12mip->Write();
  h_hgcal_sclusters_entryposition_16mip->Write();
  h_hgcal_sclusters_entry_cutsigmaeta_cuthadem->Write();
  h_hgcal_sclusters_entryVSlength_cut_fullrange->Write();
  h_hgcal_sclusters_energyVSmeanradius_fullrange->Write();
  h_hgcal_sclusters_entryVSmeanradius_fullrange->Write();
  h_hgcal_sclusters_energyVSsigmaradius_fullrange->Write();
  h_hgcal_sclusters_entryVSsigmaradius_fullrange->Write();
  h_hgcal_sclusters_energyVSlongwidth_fullrange->Write();  
  h_hgcal_sclusters_entryVSlongwidth_fullrange ->Write(); 
  h_hgcal_sclusters_energyVStranswidth_fullrange->Write();  
  h_hgcal_sclusters_entryVStranswidth_fullrange->Write();  
  h_hgcal_sclusters_energyVSeigenratio_fullrange->Write();  
  h_hgcal_sclusters_entryVSeigenratio_fullrange->Write();  
  h_hgcal_sclusters_energyVSeigenratio_cut_fullrange->Write();  
  h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_fullrange->Write();  
  h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_cut_fullrange->Write();  
  h_hgcal_sclusters_sigmaradius_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradius_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradius_corr_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxis_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxismiddle_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxismiddle_corr_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxisVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxisVSphi_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxisVSlength_fullrange->Write();
  h_hgcal_sclusters_sigmatransverseradiusaxis_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaradiusnorm_fullrange->Write();
  h_hgcal_sclusters_sigmaetanorm20_fullrange->Write();
  h_hgcal_sclusters_sigmaetanorm50_fullrange->Write();
  h_hgcal_sclusters_sigmaetanorm100_fullrange->Write();
  h_hgcal_sclusters_sigmaeta20_fullrange->Write();
  h_hgcal_sclusters_sigmaeta50_fullrange->Write();
  h_hgcal_sclusters_sigmaeta100_fullrange->Write();
  h_hgcal_sclusters_sigmaradiusVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmaradiusVSphi_fullrange->Write();
  h_hgcal_sclusters_sigmaphi_fullrange->Write();
  h_hgcal_sclusters_sigmaphinorm_fullrange->Write();
  h_hgcal_sclusters_sigmaphiVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmaphiVSphi_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_pu_fullrange->Write();
  h_hgcal_sclusters_sigmaetanorm_fullrange->Write();
  h_hgcal_sclusters_sigmaetaVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_puVSeta_fullrange->Write();
  h_hgcal_sclusters_etaVSsigmaeta_fullrange->Write();
  h_hgcal_sclusters_etaVSsigmatransverseradius_fullrange->Write();
  h_hgcal_sclusters_etaVSsigmatransverseradiusaxis_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_pu_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaetaVSpt_pu_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaetaVSeta_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_puVSeta_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_puVSEt_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaeta_puVSseedfraction_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaeta10_puVSEt_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaphi_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaphicorrVSphiminusphitrue->Write();
  h_hgcal_sclusters_sigmaphiVSeta_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaphiVSEt_corr_fullrange->Write();
  h_hgcal_sclusters_sigmaetaw_fullrange->Write();
  h_hgcal_sclusters_sigmaetawnorm_fullrange->Write();
  h_hgcal_sclusters_sigmaetawVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmaetawnormVSeta_fullrange->Write();
  h_hgcal_sclusters_sigmaetaw200_fullrange->Write();

  h_hgcal_sclusters_e2530OverEtot->Write();
  h_hgcal_sclusters_e0110OverEtot->Write();
  h_hgcal_sclusters_e0110OverEtotVSe2530OverEtot->Write();
  
  //std::vector<TProfile*> h_hgcal_sclusters_etaVSsigmaeta_pfx;
  h_hgcal_sclusters_sigmaetaVSlayer->Write();
  h_hgcal_sclusters_sigmaetaVSlayer_norm->Write();
  h_hgcal_sclusters_sigmaeta_1->Write();
  h_hgcal_sclusters_sigmaeta_2->Write();
  h_hgcal_sclusters_sigmaeta_3->Write();
  h_hgcal_sclusters_sigmaeta_4->Write();
  h_hgcal_sclusters_sigmaeta_5->Write();
  h_hgcal_sclusters_sigmaeta_6->Write();
  h_hgcal_sclusters_sigmaeta_7->Write();
  h_hgcal_sclusters_sigmaeta_8->Write();
  h_hgcal_sclusters_sigmaeta_9->Write();
  h_hgcal_sclusters_sigmaeta_10->Write();
  h_hgcal_sclusters_sigmaeta_11->Write();
  h_hgcal_sclusters_sigmaeta_12->Write();
  h_hgcal_sclusters_sigmaeta_13->Write();
  h_hgcal_sclusters_sigmaeta_14->Write();
  h_hgcal_sclusters_sigmaeta_15->Write();
  h_hgcal_sclusters_sigmaeta_16->Write();
  h_hgcal_sclusters_sigmaeta_17->Write();
  h_hgcal_sclusters_sigmaeta_18->Write();
  h_hgcal_sclusters_sigmaeta_19->Write();
  h_hgcal_sclusters_sigmaeta_20->Write();
  h_hgcal_sclusters_sigmaeta_21->Write();
  h_hgcal_sclusters_sigmaeta_22->Write();
  h_hgcal_sclusters_sigmaeta_23->Write();
  h_hgcal_sclusters_sigmaeta_24->Write();
  h_hgcal_sclusters_sigmaeta_25->Write();
  h_hgcal_sclusters_sigmaeta_26->Write();
  h_hgcal_sclusters_sigmaeta_27->Write();
  h_hgcal_sclusters_sigmaeta_28->Write();
  h_hgcal_sclusters_sigmaeta_29->Write();
  h_hgcal_sclusters_sigmaeta_30->Write();
  h_hgcal_sclusters_etaVSsigmaeta_1->Write();
  h_hgcal_sclusters_etaVSsigmaeta_2->Write();
  h_hgcal_sclusters_etaVSsigmaeta_3->Write();
  h_hgcal_sclusters_etaVSsigmaeta_4->Write();
  h_hgcal_sclusters_etaVSsigmaeta_5->Write();
  h_hgcal_sclusters_etaVSsigmaeta_6->Write();
  h_hgcal_sclusters_etaVSsigmaeta_7->Write();
  h_hgcal_sclusters_etaVSsigmaeta_8->Write();
  h_hgcal_sclusters_etaVSsigmaeta_9->Write();
  h_hgcal_sclusters_etaVSsigmaeta_10->Write();
  h_hgcal_sclusters_etaVSsigmaeta_11->Write();
  h_hgcal_sclusters_etaVSsigmaeta_12->Write();
  h_hgcal_sclusters_etaVSsigmaeta_13->Write();
  h_hgcal_sclusters_etaVSsigmaeta_14->Write();
  h_hgcal_sclusters_etaVSsigmaeta_15->Write();
  h_hgcal_sclusters_etaVSsigmaeta_16->Write();
  h_hgcal_sclusters_etaVSsigmaeta_17->Write();
  h_hgcal_sclusters_etaVSsigmaeta_18->Write();
  h_hgcal_sclusters_etaVSsigmaeta_19->Write();
  h_hgcal_sclusters_etaVSsigmaeta_20->Write();
  h_hgcal_sclusters_etaVSsigmaeta_21->Write();
  h_hgcal_sclusters_etaVSsigmaeta_22->Write();
  h_hgcal_sclusters_etaVSsigmaeta_23->Write();
  h_hgcal_sclusters_etaVSsigmaeta_24->Write();
  h_hgcal_sclusters_etaVSsigmaeta_25->Write();
  h_hgcal_sclusters_etaVSsigmaeta_26->Write();
  h_hgcal_sclusters_etaVSsigmaeta_27->Write();
  h_hgcal_sclusters_etaVSsigmaeta_28->Write();
  h_hgcal_sclusters_etaVSsigmaeta_29->Write();
  h_hgcal_sclusters_etaVSsigmaeta_30->Write();
  h_hgcal_sclusters_sigmaetaVSlayer_cutlayer9->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_1->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_2->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_3->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_4->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_5->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_6->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_7->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_8->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_9->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_10->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_11->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_12->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_13->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_14->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_15->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_16->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_17->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_18->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_19->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_20->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_21->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_22->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_23->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_24->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_25->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_26->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_27->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_28->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_29->Write();
  h_hgcal_sclusters_sigmaeta_cutlayer9_30->Write();
  h_hgcal_sclusters_longitudinal->Write();  
  h_hgcal_sclusters_layer->Write();  
  h_hgcal_sclusters_lengthCompatibility->Write();  
  h_hgcal_sclusters_transversal->Write();  
  h_hgcal_sclusters_transversalaxis->Write();  
  h_hgcal_sclusters_transversalaxis_calib->Write();  
  h_hgcal_sclusters_deta_shower->Write();  
  h_hgcal_sclusters_firsteigenvalue->Write();  
  h_hgcal_sclusters_firsteigenvalue_nm1->Write();  
  h_hgcal_sclusters_secondeigenvalue->Write();  
  h_hgcal_sclusters_secondeigenvalue_nm1->Write();  
  h_hgcal_sclusters_thirdeigenvalue->Write();  
  h_hgcal_sclusters_thirdeigenvalue_nm1->Write();  
  h_hgcal_sclusters_firstsigma->Write();  
  h_hgcal_sclusters_firstsigma_nm1->Write();  
  h_hgcal_sclusters_secondsigma->Write();  
  h_hgcal_sclusters_secondsigma_nm1->Write();  
  h_hgcal_sclusters_thirdsigma->Write();  
  h_hgcal_sclusters_thirdsigma_nm1->Write();  
  h_hgcal_sclusters_firsteigenvalueVSeta->Write();  
  h_hgcal_sclusters_secondeigenvalueVSeta->Write();  
  h_hgcal_sclusters_firsteigenvalueVSsecond->Write();  
  h_hgcal_sclusters_transeigenvalueVSeta->Write();  
  h_hgcal_sclusters_eigenratioVSeta->Write();  
  h_hgcal_sclusters_transversalVSlongitudinal->Write();  
  h_hgcal_sclusters_transversalVSlongitudinalscaled->Write();  
  h_hgcal_sclusters_transversalVSlongitudinalscaledmeasured->Write();  
  h_hgcal_sclusters_transversalaxisVSlongitudinal->Write();  
  h_hgcal_sclusters_transversalaxisVSlongitudinalscaled->Write();  
  h_hgcal_sclusters_transversalaxisVSlongitudinalscaledmeasured->Write();  
  h_hgcal_sclusters_transversalVSeta->Write();  
  h_hgcal_sclusters_longitudinal_cut20->Write();  
  h_hgcal_sclusters_transversal_cut20->Write();  
  h_hgcal_sclusters_transversalVSlongitudinal_cut20->Write();  
  h_hgcal_sclusters_longitudinal_cut50->Write();  
  h_hgcal_sclusters_transversal_cut50->Write();  
  h_hgcal_sclusters_transversalVSlongitudinal_cut50->Write();  
  h_hgcal_sclusters_longitudinal_cut100->Write();  
  h_hgcal_sclusters_transversal_cut100->Write();  
  h_hgcal_sclusters_transversalVSlongitudinal_cut100->Write();  
  h_hgcal_sclusters_longitudinal_cut200->Write();  
  h_hgcal_sclusters_transversal_cut200->Write();  
  h_hgcal_sclusters_transversalVSlongitudinal_cut200->Write();  
  
  h_hgcal_sclusters_longitudinal_fit_chi2->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2_pnorm->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2_bounded->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2_nm1->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2VSeta->Write();  
  h_hgcal_sclusters_longitudinal_fit_alphaVSbeta->Write();  
  h_hgcal_sclusters_longitudinal_fit_alphaVSenergy->Write();  
  h_hgcal_sclusters_longitudinal_fit_betaVSenergy->Write();  
  h_hgcal_sclusters_longitudinal_fit_invbetaVSenergy->Write();  
  h_hgcal_sclusters_longitudinal_fit_alphaVSinvbeta->Write();  
  h_hgcal_sclusters_longitudinal_fit_leakage->Write();  
  h_hgcal_sclusters_longitudinal_fit_leakageVShoverem->Write();  
  h_hgcal_sclusters_longitudinal_fit_leakage_cutseedpos->Write();  
  h_hgcal_sclusters_longitudinal_fit_leakageVShoverem_cutseedpos->Write();  
  h_hgcal_sclusters_longitudinal_kolmogorov_prob->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_prob_best->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem_cutseedpos->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_nm1->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem_cutseedpos->Write();   
  h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_nm1->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_prob->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_dist->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_prob_first->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_dist_first->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_prob_last->Write();   
  h_hgcal_sclusters_transversal_kolmogorov_dist_last->Write();   
  h_hgcal_sclusters_kolmogorov_prob_3D->Write();   
  h_hgcal_sclusters_kolmogorov_dist_3D->Write();   
  h_hgcal_sclusters_kolmogorov_dist_3D_px->Write();   
  h_hgcal_sclusters_kolmogorov_dist_3D_py->Write();   
  h_hgcal_sclusters_longitudinal_fit_chi2VSseedfraction->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2VSdphidir->Write();  
  h_hgcal_sclusters_longitudinal_fit_chi2VSeoveretrue->Write();  
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal->Write();    
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_nm1->Write();    
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cutsigmaeta_cuthadem->Write();    
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem->Write();    
  h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem_cutseedpos->Write();    
  
  h_hgcal_scclusters_eoveretrue_em_cutpos->Write();   ;
  h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength->Write();   ;
  h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta->Write();   ; 
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos->Write();   ;
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength->Write();   ;
  h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength_cutsigmaeta->Write();   ; 
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos->Write();   ;
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength->Write();   ;
  h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength_cutsigmaeta->Write();   ; 

  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta  ->Write();
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem  ->Write();
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutkolmogorov  ->Write();
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos->Write();
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength->Write();
  h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength_cutkolmogorov->Write();

  h_hgcal_sclusters_hoverem->Write();
  h_hgcal_sclusters_hoverem_cone01->Write();
  h_hgcal_sclusters_hoverem_cutsigmaeta->Write();
  h_hgcal_sclusters_hoveremVSeta  ->Write();
  h_hgcal_sclusters_hoveremVSphi  ->Write();
  h_hgcal_sclusters_entryVShoverem  ->Write();
  h_hgcal_sclusters_expectedlengthVShoverem  ->Write();
  h_hgcal_sclusters_hoverem1  ->Write();
  h_hgcal_sclusters_hoverem1VSeta   ->Write();
  h_hgcal_sclusters_entryVShoverem1  ->Write();
  h_hgcal_sclusters_expectedlengthVShoverem1  ->Write();
  h_hgcal_sclusters_hoverem2  ->Write();
  h_hgcal_sclusters_hoverem2VSeta  ->Write();
  h_hgcal_sclusters_entryVShoverem2->Write();
  h_hgcal_sclusters_expectedlengthVShoverem2  ->Write();
  
  h_hgcal_allclusters_predictedTvsy ->Write(); 
  h_hgcal_allclusters_predictedLnTvsy->Write();
  h_hgcal_allclusters_predictedLength_fullrange->Write();
  	
  h_hgcal_allclusters_transverse_energy_em->Write();
  h_hgcal_allclusters_length_fullrange->Write();
  h_hgcal_allclusters_energyVSlength_fullrange->Write();
  h_hgcal_allclusters_energyVSlength_cut_fullrange->Write();
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_fullrange->Write();
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange->Write();
  h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange_weighted->Write();
  h_hgcal_allclusters_energyVSpredictedLength_fullrange->Write();
  h_hgcal_allclusters_energyVSpredictedLength_cut_fullrange->Write();
  h_hgcal_allclusters_entryVSlength_fullrange->Write();
  h_hgcal_allclusters_entryVSlength_cut_fullrange->Write();

  h_hgcal_puclusters_energy ->Write();
  h_hgcal_puclusters_energyVSeta ->Write(); 
  h_hgcal_puclusters_energy_eta16 ->Write();  
  h_hgcal_puclusters_energy_eta20 ->Write();  
  h_hgcal_puclusters_energy_eta25 ->Write();  
  h_hgcal_puclusters_energy_eta29->Write();
   
  h_hgcal_allhits_energy_em->Write();
    
}

HGCALElectronClusterAnalyzer::~HGCALElectronClusterAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  histfile_->Write();
  histfile_->Close();

}

void HGCALElectronClusterAnalyzer::doTest(const HGCalGeometry& geom, ForwardSubdetector subdet) {

  const std::vector<DetId>& ids = geom.getValidDetIds();
  std::cout << ids.size() << " valid ids for " << geom.cellElement() << std::endl;
  int sectors[]= {1, 7, 13};
  int layers[] = {1, 5, 10};
  int zsides[] = {1, -1};
  int cells[] = {1, 101, 201};
  for (int iz = 0; iz < 2; ++iz) {
    int zside = zsides[iz];
    for (int is = 0; is < 3; ++is) {
      int sector = sectors[is];
      for (int il = 0; il < 3; ++il) {
        int layer = layers[il];
        for (int ic = 0; ic < 3; ++ic) {
	 int cell = cells[ic];
	 DetId id1= ((subdet == HGCEE) ? (DetId)(HGCEEDetId(subdet,zside,layer,sector,0,cell)) : (DetId)(HGCHEDetId(subdet,zside,layer,sector,0,cell)));
	 const CaloCellGeometry* icell1 = geom.getGeometry(id1);
	 GlobalPoint global1 = geom.getPosition(id1);
	 DetId idc1 = geom.getClosestCell(global1);
	 std::cout << "DetId (" << subdet << ":" << zside << ":" << layer
	 << ":" << sector << ":0:" << cell << ") Geom " << icell1
	 << " position (" << global1.x() << ", " << global1.y()
	 << ", " << global1.z() << ") rawId " << idc1.rawId() << " ids " << std::hex
	 << id1.rawId() << ":" << idc1.rawId() << std::dec
	 << " parameter[11] = " << icell1->param()[10] << ":"
	 << icell1->param()[11] << std::endl;
	 if (id1.rawId() != idc1.rawId()) std::cout << "***** ERROR *****\n";
	 DetId id2= ((subdet == HGCEE) ?
	 (DetId)(HGCEEDetId(subdet,zside,layer,sector,1,cell)) :
	 (DetId)(HGCHEDetId(subdet,zside,layer,sector,1,cell)));
	 const CaloCellGeometry* icell2 = geom.getGeometry(id2);
	 GlobalPoint global2 = geom.getPosition(id2);
	 DetId idc2 = geom.getClosestCell(global2);
	 std::cout << "DetId (" << subdet << ":" << zside << ":" << layer
	 << ":" << sector << ":1:" << cell << ") Geom " << icell2
	 << " position (" << global2.x() << ", " << global2.y()
	 << ", " << global2.z() << ") ids " << std::hex
	 << id2.rawId() << ":" << idc2.rawId() << std::dec
	 << " parameter[11] = " << icell2->param()[10] << ":"
	 << icell2->param()[11] << std::endl;
	 if (id2.rawId() != idc2.rawId()) std::cout << "***** ERROR *****\n";
        }
      }
    }
  }
  uint32_t probids[] = {1711603886, 1711603890, 1761408735, 1761411303, 1801744385, 1805447194};
  for (int k=0; k<6; ++k) {
    DetId id(probids[k]);
    if (id.det() == DetId::Forward && id.subdetId() == (int)(subdet)) {
      if (subdet == HGCEE) std::cout << "Test " << HGCEEDetId(id) << std::endl;
      else std::cout << "Test " << HGCHEDetId(id) << std::endl;
      const CaloCellGeometry* icell = geom.getGeometry(id);
      GlobalPoint global = geom.getPosition(id);
      std::cout << "Geom Cell: " << icell << " position (" << global.x() << ", " << global.y() << ", " << global.z() << ")"<< std::endl;
    }
  }

}

GlobalVector HGCALElectronClusterAnalyzer::rotateMomentum(const MagneticField *magField,GlobalVector momentum, GlobalPoint xmeas, GlobalPoint xvert, int charge) {

  //std::cout << "[HGCALElectronClusterAnalyzer::rotateMomentum] entering, magfield " << magField << std::endl;
  //double BInTesla = magField->inTesla(xmeas).z();
  double BInTesla = magField->inTesla(GlobalPoint(0.,0.,0.)).z();
  //std::cout << "[HGCALElectronClusterAnalyzer::rotateMomentum] B field " << BInTesla << std::endl;
  GlobalVector xdiff = xmeas - xvert;
  //double theta = xdiff.theta();
  //double phi= xdiff.phi();  
  double pt = momentum.perp();
  double pxOld = momentum.x();
  double pyOld = momentum.y();
  double pz = momentum.z();  
  double RadCurv = 100*pt/(BInTesla*0.29979);
  //std::cout << "[HGCALElectronClusterAnalyzer::rotateMomentum] argument of asin " << 0.5*xdiff.perp()/RadCurv << std::endl;
  double alpha = asin(0.5*xdiff.perp()/RadCurv);
  float ca = cos(charge*alpha);
  float sa = sin(charge*alpha);
  double pxNew = ca*pxOld + sa*pyOld;
  double pyNew = -sa*pxOld + ca*pyOld;
  
  return GlobalVector(pxNew, pyNew, pz);

}


// FreeTrajectoryState FTSFromVertexToPointFactory::get( MagneticField const & magField,
// GlobalPoint const & xmeas,
// GlobalPoint const & xvert,
// float momentum,
// TrackCharge charge )

//=========================================================================
// Main method
//=========================================================================

void
HGCALElectronClusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "analyzing new event " << std::endl;
  ievent++;
  
  // get electrons

  edm::Handle<GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons);
  std::cout << "=================> Treating event "<<iEvent.id()<<" Number of electrons "<<gsfElectrons.product()->size()<<std::endl;
  
  edm::Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(mcTruthCollection_, genParticles);

// //   simVertexs_g4simhits__SIM => Thevertex.Z()
//    edm::Handle<edm::SimVertexContainer> simVertices;
//    iEvent.getByLabel("g4SimHits", simVertices);
//    if (ievent<10) {
//      std::vector<SimVertex>::const_iterator iv;   
//      for (iv=simVertices->begin(); iv!=simVertices->end(); iv++) {
//      std::cout << "new sim vertex " << (*iv) << std::endl;
//      }
//    }
//   
  // hgcal super clusters 
  edm::Handle<SuperClusterCollection> superClusters ;
  iEvent.getByLabel(endcapSuperClusterCollection_,superClusters);
 
  edm::Handle<PFClusterCollection> clusters ;
  iEvent.getByLabel(endcapClusterCollection_,clusters);

  edm::Handle<PFClusterCollection> hcalClusters ;
  iEvent.getByLabel(edm::InputTag("particleFlowClusterHGCHEF"),hcalClusters);
  //const std::vector<reco::PFCluster>* hcalClusterColl = hcalClusters.product();
  
  edm::Handle<HGCRecHitCollection> recHits;
  iEvent.getByLabel(endcapRecHitCollection_,recHits);
  
  edm::ESHandle<HGCalGeometry> pGeometry;
  unsigned long long newCaloGeomCacheId= iSetup.get<IdealGeometryRecord>().cacheIdentifier() ;
  if (caloGeomCacheId_!=newCaloGeomCacheId) {
    caloGeomCacheId_ = newCaloGeomCacheId ;
//    iSetup.get<IdealGeometryRecord>().get(pGeometry) ;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",pGeometry) ;
    geometry_ = pGeometry.product() ;
  }  
  //std::cout << "geometry_ " << geometry_ << std::endl;
  
  edm::ESHandle<MagneticField> theMagField;
  unsigned long long newCacheIdMagField= iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier() ;
  if (cacheIDMagField_!=newCacheIdMagField) {
    cacheIDMagField_=newCacheIdMagField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
    magField_ = theMagField.product();   
  }
  
  // instantiate selection object
  HGCALShowerBasedEmIdentification
   hGCALShowerBasedEmIdentification(iEvent,iSetup,withPileup_);

/*  
  if (ievent<10) {
  
  // analyze all recHits
  HGCRecHitCollection::const_iterator ih;   
  for (ih=recHits->begin(); ih!=recHits->end(); ih++) {
    if ((*ih).detid().det()==DetId::Forward && (*ih).detid().subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((*ih).detid()) ;
      //GlobalPoint cellPos = geometry_->getGeometry(hgcid_)->getPosition();
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);
      //std::cout << "new HGCAL recHits with position " <<  cellPos << "and energy " << ih->energy() << std::endl;
      // only z>0 and x>0 and y>0
      if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.) {
       h_hgcal_allhits_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
       h_hgcal_allhits_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),ih->energy()*1000.);
      }
      h_hgcal_allhits_energy_em->Fill(ih->energy());
    }
  }
  
  }
*/

/* 
  // analyse PU clusters

  bool matchingIDpu, matchingMotherIDpu;
  reco::PFCluster bestCluster; 
  // association mc-reco
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {

    // select requested matching gen particle
    matchingIDpu=false;
    for (unsigned int i=0; i<matchingIDs_.size(); i++)
     if ( mcIter->pdgId() == matchingIDs_[i] ) matchingIDpu=true;
    //std::cout << "pdgID " << mcIter->pdgId() << std::endl;
    if (matchingIDpu) {

    // select requested mother matching gen particle
    // always include single particle with no mother
    const Candidate * mother = mcIter->mother();
    matchingMotherIDpu=false;
    for (unsigned int i=0; i<matchingMotherIDs_.size(); i++)
     if ((mother == 0) || ((mother != 0) &&  mother->pdgId() == matchingMotherIDs_[i]) ) matchingMotherIDpu=true;

    if (matchingMotherIDpu) {

      //if (mcIter->pt()> maxPt_ || std::abs(mcIter->eta())> maxAbsEta_) continue;

      // suppress the barrel
      if (std::abs(mcIter->eta()) < 1.5) continue;
  
      // looking for the best matching mc particle
//      bool okClusterFound = false;
      double gsfOkRatio = 999999.;

      //// find best matched electron
      //reco::GsfElectron bestGsfElectron;
      //for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
      // gsfIter!=gsfElectrons->end(); gsfIter++){
      // find best matched gen particle to reco supercluster
      for (reco::PFClusterCollection::const_iterator gsfIter=clusters->begin();
       gsfIter!=clusters->end(); gsfIter++){

	double dphi = gsfIter->position().phi()-mcIter->phi();
	if (std::abs(dphi)>CLHEP::pi)
	 dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	double deltaR = sqrt(std::pow((gsfIter->position().eta()-mcIter->eta()),2) + std::pow(dphi,2));
	if ( deltaR < deltaR_ ){
	  double tmpGsfRatio = gsfIter->energy()/mcIter->p();
	  if ( std::abs(tmpGsfRatio-1) < std::abs(gsfOkRatio-1) ) {
	    gsfOkRatio = tmpGsfRatio;
	    bestCluster=*gsfIter;
//	    okClusterFound = true;
	  }
	}
      } // loop over rec ele to look for the best one

    }
    
  } 
  
  }
  
//  }
   
  // here we have found the SC best matching the generated particle
  
  // analyse all pfclusters
  for (unsigned int i=0;i<clusters->size();++i) {

    const PFCluster & cl = (*clusters)[i] ;      
    
    // require the cluster to be away from the gen particle
	double dphi = bestCluster.position().phi()-cl.phi();
	if (std::abs(dphi)>CLHEP::pi)
	 dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	double deltaR = sqrt(std::pow((bestCluster.position().eta()-cl.eta()),2) + std::pow(dphi,2));
	if ( deltaR < deltaR_ ) continue;
      
    // Et cut
    double clet = cl.energy()*sin(2.*atan(exp(-cl.position().eta())));
    
    h_hgcal_puclusters_energy->Fill(cl.energy());
    h_hgcal_puclusters_energyVSeta->Fill(cl.energy(),fabs(cl.eta()));
    if (fabs(cl.eta())>1.55 && fabs(cl.eta())<1.65) h_hgcal_puclusters_energy_eta16->Fill(cl.energy());
    if (fabs(cl.eta())>1.95 && fabs(cl.eta())<2.05) h_hgcal_puclusters_energy_eta20->Fill(cl.energy());
    if (fabs(cl.eta())>2.45 && fabs(cl.eta())<2.55) h_hgcal_puclusters_energy_eta25->Fill(cl.energy());
    if (fabs(cl.eta())>2.85 && fabs(cl.eta())<2.95) h_hgcal_puclusters_energy_eta29->Fill(cl.energy());
    //if (clet < 0.1) continue;
    //if (clet < 0.5) continue;
    //if (clet < 5.0) continue;
    // eta cut 
    //if (fabs(cl.position().eta())<2.0) continue;
    h_hgcal_allclusters_transverse_energy_em->Fill(clet);
    
    // skip PCA for all clusters that takes time
    continue;
    
    PCAShowerAnalysis pcaShowerAnalysisAll(iEvent,iSetup);
    GlobalPoint firstPos;	
    double zmin = 400.; 
    if (ievent<20) std::cout << "new pf cluster in HGCAL with energy " << cl.energy() << " and position " << cl.position() << std::endl;

    for (unsigned int ih=0;ih<cl.hitsAndFractions().size();++ih) {
      const DetId & id_ = (cl.hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
      //if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
	GlobalPoint cellPos = geometry_->getPosition(hgcid_);
	if (fabs(cellPos.z())<zmin) {
	  firstPos = cellPos;
	  zmin = fabs(cellPos.z());
	}
//	if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.){
//          if (ievent<10) h_hgcal_clustershits_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
//          if (ievent<10) h_hgcal_clustershits_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);
//	}
      }    
    }
    
    // analyse position and direction with PCA 
    GlobalPoint pcaShowerPos;
    GlobalVector pcaShowerDir;
    pcaShowerAnalysisAll.showerParameters(&cl,pcaShowerPos,pcaShowerDir);
    if (ievent<20) std::cout << "*** Principal component analysis (all clusters) ****" << std::endl;
    if (ievent<20) std::cout << "shower average (x,y,z) = " << pcaShowerPos << std::endl;
    if (ievent<20) std::cout << "shower main axis (x,y,z) = " << pcaShowerDir << std::endl;
    double length =  (pcaShowerPos - firstPos).mag();
    if (ievent<20) std::cout << "shower length " << length << " first pos " << firstPos << std::endl;	 

    // now generate the fastsim shower for this candidate
    // change here for parametrisation from pileup pi0 in fullsim, assume no signal here
    GlobalVector dir = pcaShowerDir.unit();
    const XYZTLorentzVector momentum(dir.x()*cl.energy(),dir.y()*cl.energy(),dir.z()*cl.energy(),cl.energy());
    //RawParticle myEle(22,momentum);
    //std::vector<const RawParticle *> thePart;
    //if ( myEle.e() > 0.055 ) thePart.push_back(&myEle);
    //EMShower theShower(random_,aGammaGenerator_,showerparam_,&thePart,NULL,NULL,NULL,false); 
    ////double maxShower = theShower.getMaximumOfShower();	 
    ////double predictedLength = maxShower * calohelper_->ecalProperties(true)->radLenIncm();
    //double meanShower = theShower.getMeanDepth();	 
    double x0 = calohelper_->ecalProperties(true)->radLenIncm();
    //double predictedLength = meanShower * x0;
    double predictedLength = 3.6 + 1.383*log(cl.energy());
    
    // histogram parametrised shower quantities
    //h_hgcal_allclusters_predictedTvsy->Fill(meanShower,momentum.e()/calohelper_->ecalProperties(true)->criticalEnergy());
    //h_hgcal_allclusters_predictedLnTvsy->Fill(log(meanShower),momentum.e()/calohelper_->ecalProperties(true)->criticalEnergy());
    h_hgcal_allclusters_predictedLength_fullrange->Fill(predictedLength,log(momentum.e()));	  

    // histogram measured shower quantities
    h_hgcal_allclusters_length_fullrange->Fill(length,log(cl.energy()));
    h_hgcal_allclusters_energyVSlength_fullrange->Fill(log(cl.energy()),length);
    h_hgcal_allclusters_energyVSpredictedLength_fullrange->Fill(log(cl.energy()),predictedLength);
    h_hgcal_allclusters_entryVSlength_fullrange->Fill(fabs(firstPos.z()),length);
    // for photon and pi0 goes up to 5th layer
    if (fabs(firstPos.z())<322.50) h_hgcal_allclusters_energyVSlength_hasfirstlayer_fullrange->Fill(log(cl.energy()),length);
    // cut in length, inject here parametrization results
    double y = cl.energy()/0.00536;
    double sigma = predictedLength / (-2.506+1.245*log(y));
    bool cut = fabs(predictedLength-length)<4.*sigma/x0;
    if (ievent<20 && fabs(firstPos.z())>322.50) std::cout << "cluster rejected from firstpos cut" << std::endl;
    if (ievent<20 && !cut) std::cout << "cluster rejected from length cut, predicted length " << predictedLength << " sigma " << sigma/x0 << std::endl;
    if (cut) h_hgcal_allclusters_energyVSlength_cut_fullrange->Fill(log(cl.energy()),length);
    if (cut) h_hgcal_allclusters_energyVSpredictedLength_cut_fullrange->Fill(log(cl.energy()),predictedLength);
    if (fabs(firstPos.z())<322.50 && cut) h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange->Fill(log(cl.energy()),length);
    if (fabs(firstPos.z())<322.50 && cut) h_hgcal_allclusters_energyVSlength_hasfirstlayer_cut_fullrange_weighted->Fill(log(cl.energy()),length,cl.energy());
    if (cut) h_hgcal_allclusters_entryVSlength_cut_fullrange->Fill(fabs(firstPos.z()),length);

    // shower display after cuts
//    if (fabs(firstPos.z())<322.50) {
//      for (unsigned int ih=0;ih<cl.hitsAndFractions().size();++ih) {
//	const DetId & id_ = (cl.hitsAndFractions())[ih].first ;
//	HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
//	//if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
//	if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
//	  const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
//	  GlobalPoint cellPos = geometry_->getPosition(hgcid_);
//	  if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.){
//            if (ievent<10) h_hgcal_clustershits_cutlayers_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
//            if (ievent<10) h_hgcal_clustershits_cutlayers_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);
//            if (cut) {
//              if (ievent<10) h_hgcal_clustershits_cutlayers_cutlength_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
//              if (ievent<10) h_hgcal_clustershits_cutlayers_cutlength_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);	    
//	    }
//	  }
//	}    
//      }
//    }
           
  }
 
//  if (ievent<20) {

    // analyse all superclusters in HGCAL
    for (unsigned int i=0;i<superClusters->size();++i) {
      const SuperCluster & scl = (*superClusters)[i] ;      
      std::cout << "new em supercluster in HGCAL with energy " << scl.energy() << " and position " << scl.position() << std::endl;
 
     // eta cut 
     //if (fabs(scl.position().eta())<2.0) return;

     // access all sc cluster hits 
      for (CaloCluster_iterator itcl=scl.clustersBegin(); itcl!=scl.clustersEnd(); itcl++) {

	// Et cut
	//double clet = (*itcl)->energy()*sin(2.*atan(exp(-(*itcl)->position().eta())));
	//if (clet < 0.1) continue;
	//if (clet < 0.5) continue;
	//if (clet < 5.0) continue;

	PCAShowerAnalysis pcaShowerAnalysisSC(iEvent,iSetup);
	GlobalPoint firstPos;	
	double zmin = 400.; 

	for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
	  const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
	  HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
	  //if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
	  if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	    const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
	    GlobalPoint cellPos = geometry_->getPosition(hgcid_);
	if (fabs(cellPos.z())<zmin) {
	  firstPos = cellPos;
	  zmin = fabs(cellPos.z());
	}
	    if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.) {	
	      if (ievent<10) h_hgcal_superclustershits_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
	      if (ievent<10 )h_hgcal_superclustershits_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);
	    }
	  }
	}

	// analyse position and direction with PCA 
	GlobalPoint pcaShowerPos;
	GlobalVector pcaShowerDir;
	pcaShowerAnalysisSC.showerParameters(&(**itcl),pcaShowerPos,pcaShowerDir);
	if (ievent<1) std::cout << "*** Principal component analysis (SC clusters) ****" << std::endl;
	if (ievent<1) std::cout << "shower average (x,y,z) = " << pcaShowerPos << std::endl;
	if (ievent<1) std::cout << "shower main axis (x,y,z) = " << pcaShowerDir << std::endl;
	double length =  (pcaShowerPos - firstPos).mag();
	if (ievent<1) std::cout << "shower length " << length << " first pos " << firstPos << std::endl;	 

	// cut for all but seed cluster
	double x0 = calohelper_->ecalProperties(true)->radLenIncm();
	double predictedLength = 3.6 + 1.383*log((*itcl)->energy());
	double y = (*itcl)->energy()/0.00536;
	double sigma = predictedLength / (-2.506+1.245*log(y));
	bool cutlength = fabs(predictedLength-length)<4.*sigma/x0;
	bool cutpos = (fabs(firstPos.z())<322.50);

	// now generate the fastsim shower for the seed of this candidate
	if (itcl==scl.clustersBegin()) {
	  GlobalVector dir = pcaShowerDir.unit();
	  const XYZTLorentzVector momentum(dir.x()*(*itcl)->energy(),dir.y()*(*itcl)->energy(),dir.z()*(*itcl)->energy(),(*itcl)->energy());
	  RawParticle myEle(22,momentum);
	  std::vector<const RawParticle *> thePart;
	  if ( myEle.e() > 0.055 ) thePart.push_back(&myEle);
	  EMShower theShower(random_,aGammaGenerator_,showerparam_,&thePart,NULL,NULL,NULL,false); 
	  //double maxShower = theShower.getMaximumOfShower();	 
	  //double predictedLength = maxShower * calohelper_->ecalProperties(true)->radLenIncm();
	  double meanShower = theShower.getMeanDepth();	 
	  predictedLength = meanShower * x0;
	  // cut in length, inject here parametrization results
	  sigma = predictedLength / (-2.506+1.245*log(y));
	  cutlength = fabs(predictedLength-length)<4.*sigma/x0;
	  cutpos = (fabs(firstPos.z())<323.5);
	}

	// shower display after cuts
	bool cutall = (cutlength && cutpos);
	if (cutpos) {
	  for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
	    const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
	    HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
	    //if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
	    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	      const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
	      GlobalPoint cellPos = geometry_->getPosition(hgcid_);
	      if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.){
       	if (ievent<10) h_hgcal_superclustershits_cutlayers_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
       	if (ievent<10) h_hgcal_superclustershits_cutlayers_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);
       	if (cutall) {
       	  if (ievent<10) h_hgcal_superclustershits_cutlayers_cutlength_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
       	  if (ievent<10) h_hgcal_superclustershits_cutlayers_cutlength_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),theHit->energy()*1000.);	    
		}
	      }
	    }    
	  }
	} // end if cutpos

      } // end loop calclusters
      
    } // end loop superclusters 
*/ 
  // analyse superclusters in HGCAL
  for (unsigned int i=0;i<superClusters->size();++i) {
    const SuperCluster & scl = (*superClusters)[i] ;      
    //std::cout << "new em supercluster in HGCAL with energy " << scl.energy() << " and position " << scl.position() << std::endl;
     // access all sc cluster hits 
    h_hgcal_foundClustersVSeta_em ->Fill(scl.position().eta(),1.);
    if (scl.energy()*sin(2.*atan(exp(-scl.position().eta())))>5.) h_hgcal_foundClustersVSetaEt5_em ->Fill(scl.position().eta(),1.);
  }  
   
/*
  // test geometry
  std::string name;
  edm::ESHandle<HGCalGeometry> geom;
  name = "HGCalEESensitive";
  iSetup.get<IdealGeometryRecord>().get(name,geom);
  if (geom.isValid()) doTest(*geom, HGCEE);
  else std::cout << "Cannot get valid HGCalGeometry Object for " << name << std::endl;
  //name = "HGCalHESiliconSensitive";
  //iSetup.get<IdealGeometryRecord>().get(name,geom);
  //if (geom.isValid()) doTest(*geom, HGCHEF);
  //else std::cout << "Cannot get valid HGCalGeometry Object for " << name << std::endl;
  //name = "HGCalHEScintillatorSensitive";
  //iSetup.get<IdealGeometryRecord>().get(name,geom);
  //if (geom.isValid()) doTest(*geom,HGCHEB);
  //else std::cout << "Cannot get valid HGCalGeometry Object for " << name << std::endl;  
*/

  int mcNum=0;
  bool matchingID, matchingMotherID;

  h_hgcal_foundClusters_em->Fill(superClusters->size());
  
  // all generated particles
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {

    // all particles
    h_simZ_all   -> Fill( mcIter->vz() );
    
    // select requested matching gen particle
    matchingID=false;
    for (unsigned int i=0; i<matchingIDs_.size(); i++)
     if ( mcIter->pdgId() == matchingIDs_[i] ) matchingID=true;
    //std::cout << "pdgID " << mcIter->pdgId() << std::endl;
    if (matchingID) {

    // select requested mother matching gen particle
    // always include single particle with no mother
    const Candidate * mother = mcIter->mother();
    matchingMotherID=false;
    for (unsigned int i=0; i<matchingMotherIDs_.size(); i++)
     if ((mother == 0) || ((mother != 0) &&  mother->pdgId() == matchingMotherIDs_[i]) ) matchingMotherID=true;

    if (matchingMotherID) {
      h_simZ_electrons   -> Fill( mcIter->vz() );
    }
    
    }
    
  }   

  
  // association mc-reco
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {

    // number of mc particles
    mcNum++;

    // select requested matching gen particle
    matchingID=false;
    for (unsigned int i=0; i<matchingIDs_.size(); i++)
     if ( mcIter->pdgId() == matchingIDs_[i] ) matchingID=true;
    const Candidate * mother = mcIter->mother();
    //std::cout << "pdgID " << mcIter->pdgId() << " pt " << mcIter->pt() << " eta " << mcIter->eta() << " mother " << (mother == 0 ? 0. : mother->pdgId()) << std::endl;
    if (matchingID) {

    // select requested mother matching gen particle
    // always include single particle with no mother
    matchingMotherID=false;
    //std::cout << "matching ID " << mcIter->pdgId() << " mother ID " << (mother == 0 ? 0. : mother->pdgId()) << std::endl; 
    for (unsigned int i=0; i<matchingMotherIDs_.size(); i++)
     if ((mother == 0) || ((mother != 0) &&  mother->pdgId() == matchingMotherIDs_[i]) ) matchingMotherID=true;
    if (matchingMotherID) {

      if (mcIter->pt()> maxPt_ || std::abs(mcIter->eta())> maxAbsEta_) continue;
      //std::cout << "matching ID " << mcIter->pdgId() << " mother ID " << (mother == 0 ? 0. : mother->pdgId()) << std::endl; 

      // suppress the barrel
      if (std::abs(mcIter->eta()) < 1.5) continue;
      
      h_simEta -> Fill( mcIter->eta() );
      h_simAbsEta -> Fill( std::abs(mcIter->eta()) );
      h_simP   -> Fill( mcIter->p() );
      h_simPt   -> Fill( mcIter->pt() );
      h_simPhi   -> Fill( mcIter->phi() );
      h_simZ   -> Fill( mcIter->vz() );
      h_simPtEta   -> Fill( mcIter->eta(),mcIter->pt() );

      // looking for the best matching mc particle
      bool okSuperClusterFound = false;
      double gsfOkRatio = 999999.;

      //// find best matched electron
      //reco::GsfElectron bestGsfElectron;
      //for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
      // gsfIter!=gsfElectrons->end(); gsfIter++){
      // find best matched gen particle to reco supercluster
      reco::SuperCluster bestSuperCluster; 
      for (reco::SuperClusterCollection::const_iterator gsfIter=superClusters->begin();
       gsfIter!=superClusters->end(); gsfIter++){

        // correct for bending and z vertex position
	GlobalPoint posclu(gsfIter->position().x(),gsfIter->position().y(),gsfIter->position().z());
//	double dphi = posclu.phi()-mcIter->phi();
      // correct for bending and z vertex position
 	double dphi = posclu.phi() - 
 	 rotateMomentum(magField_,
 	                GlobalVector(mcIter->px(),mcIter->py(),mcIter->pz()),
 			posclu,
 			GlobalPoint(mcIter->vx(),mcIter->vy(),mcIter->vz()),mcIter->charge()).phi();
	if (std::abs(dphi)>CLHEP::pi)
	 dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	double posclutheta = atan2(posclu.perp(),posclu.z()-mcIter->vz());
	double posclueta = -log(tan(posclutheta/2.));	
//	double deltaR = sqrt(std::pow((posclueta-mcIter->eta()),2) + std::pow(dphi,2));
//	if ( deltaR < deltaR_ ){
        if (fabs(posclueta-mcIter->eta()) < 0.1 && fabs(dphi) <0.3) {
	  double tmpGsfRatio = gsfIter->energy()/mcIter->p();
	  if ( std::abs(tmpGsfRatio-1) < std::abs(gsfOkRatio-1) ) {
	    gsfOkRatio = tmpGsfRatio;
	    bestSuperCluster=*gsfIter;
	    okSuperClusterFound = true;
            //std::cout << "matched SC " <<bestSuperCluster << std::endl; 
	  }
	}
      } // loop over rec ele to look for the best one

      // analysis when the mc particle is found
     if (okSuperClusterFound){

     // analyze HGCAL clusters
     for (unsigned int i=0;i<superClusters->size();++i) {

       const SuperCluster & scl = (*superClusters)[i] ;       
       
       // select only the matched SC
       double scle = scl.energy() ;
       if (scle!=bestSuperCluster.energy()) continue;
              
       // Et cut
       double sclet = scle*sin(2.*atan(exp(-scl.position().eta())));
       if (sclet<5.) continue;
       
       // to simulate a non segmented HGCAL
       bool segmented = true;
       //bool segmented = false;
       
       // truth
       GlobalVector genMomentum(mcIter->px(),mcIter->py(),mcIter->pz());
       GlobalPoint genVertex(mcIter->vx(),mcIter->vy(),mcIter->vz());

       // PCA shower analysis standalone class
       PCAShowerAnalysis pcaShowerAnalysisLogwEle(iEvent,iSetup,segmented,true,false);
       // another one without logweighting for the eigenvalues and sigmas
       PCAShowerAnalysis pcaShowerAnalysisEle(iEvent,iSetup,segmented,false,false);
      
//        // redefine the seed cluster
//        double distmin = 10000.;
//        CaloCluster_iterator itseed = scl.clustersBegin();
//        for (CaloCluster_iterator itcl=scl.clustersBegin(); itcl!=scl.clustersEnd(); itcl++) {
//          double clet = (*itcl)->energy()*sin(2.*atan(exp(-(*itcl)->position().eta())));
// 	 if (clet<5.) continue;
// 	 GlobalPoint pcaShowerPos;
// 	 GlobalVector pcaShowerDir;
// 	 pcaShowerAnalysisLogwEle.showerParameters(&(**itcl),pcaShowerPos,pcaShowerDir);
// 	 double pcaShowerPosTheta = atan2(pcaShowerPos.perp(),pcaShowerPos.z()-mcIter->vz());
// 	 double pcaShowerPosEta = -log(tan(pcaShowerPosTheta/2.));
//   	 double deta = pcaShowerPosEta-mcIter->eta();
// 	 double dphi = pcaShowerDir.phi()-rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge()).phi();
// 	 //std::cout << "momentum rotation: old momentum " << genMomentum << " new momentum " <<
// 	 // rotateMomentum(magField_,genMomentum,pcaShowerPos,GlobalPoint(0.,0.,0.),mcIter->charge()) << std::endl;
// 	 if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
// 	 double dist = sqrt(deta*deta+dphi*dphi);
//          if (dist<distmin) {
// 	   itseed = itcl;
// 	   distmin = dist;
// 	 }
//        }      
//        // and reject cases where the seed is not the best matched
//        if ((*itseed)->energy()!=scl.seed()->energy()) continue;
       
       int detector = scl.seed()->hitsAndFractions()[0].first.subdetId() ;
       int component = scl.seed()->hitsAndFractions()[0].first.det() ;

       double newenergy = 0., newenergy000 = 0., newenergy04 = 0., newenergy1 = 0., newenergy2 = 0., newenergy4 = 0., newenergy10 = 0., newenergy20 = 0.;       
       double newseedenergy = 0.;       
       int newscmultiplicity = scl.clusters().size();       

       if (component==DetId::Forward && detector==HGCEE) {

	 std::cout << "new generated electron with energy " << mcIter->p() << " , transverse energy " << mcIter->pt() << " and position (eta,phi) : (" 
	 << mcIter->eta() << "," << mcIter->phi() << ")" << std::endl;
	 std::cout << "new supercluster in HGCAL with energy " << scle << " , transverse energy " << sclet << " and position (eta,phi) : (" 
	 << scl.position().eta() << "," << scl.position().phi() << ")" << std::endl;

	 nevt++;

	 // fill a few basic hgcal super cluster histos
	 h_hgcal_sclusters_energy_em->Fill(scle);
	 h_hgcal_sclusters_transverse_energy_em->Fill(sclet);
	 if (scl.position().eta()>0.) h_hgcal_sclusters_energy_pos_em->Fill(scle);
	 if (scl.position().eta()<0.) h_hgcal_sclusters_energy_neg_em->Fill(scle);
	 h_hgcal_sclusters_energyVSeta_em->Fill(scle,scl.position().eta());
	 h_hgcal_sclusters_position_em->Fill(scl.position().eta(),scl.position().phi());
	 h_hgcal_sclusters_etawidth_em->Fill(scl.etaWidth());
	 h_hgcal_sclusters_phiwidth_em->Fill(scl.phiWidth());
	 h_hgcal_sclusters_multiplicity_em->Fill(float(scl.clusters().size()));
	 h_hgcal_sclusters_multiplicityVSeta_em->Fill(float(scl.clusters().size()),mcIter->eta());
	 h_hgcal_sclusters_seedenergy_em->Fill(scl.seed()->energy());
	 h_hgcal_sclusters_seedfractionVSeta_em->Fill(scl.seed()->energy()/mcIter->p(),mcIter->eta());
	 if (scl.clusters().size()==1) h_hgcal_scclusters_eoveretrue_golden_em->Fill(scle/mcIter->p());
	 h_hgcal_scclusters_ptoverpttrue_em->Fill(sclet/mcIter->pt());
	 h_hgcal_scclusters_eoveretrue_em->Fill(scle/mcIter->p());
	 h_hgcal_scclusters_eoveretrueVSeta_em->Fill(scle/mcIter->p(),fabs(mcIter->eta()));

 	 // define few energy sums and fill seed shower 3D histos
	 double sumhits = 0., sumhits20 = 0., sumhits50 = 0., sumhits100 = 0., sumhits200 = 0.;
	 // access seed cluster hits 
	 for (unsigned int ih=0;ih<scl.seed()->hitsAndFractions().size();++ih) {
	   const DetId & id_ = (scl.seed()->hitsAndFractions())[ih].first ;
	   h_hgcal_sclusters_seedfractions_em->Fill((scl.seed()->hitsAndFractions())[ih].second);
	   HGCRecHitCollection::const_iterator theSeedHit = recHits->find(id_);    
	   if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
 	     const HGCEEDetId & hgcid_ = HGCEEDetId((scl.seed()->hitsAndFractions())[ih].first) ;
	     GlobalPoint cellPos = geometry_->getPosition(hgcid_);
	     GlobalPoint projPos(cellPos.x(),cellPos.y(),320.38);
	     // std::cout << "new cell subdet:zside:layer:sector,:subsector:cell position " << hgcid_.subdet() << ":" <<
	     // hgcid_.zside() << ":" << hgcid_.layer() << ":" << hgcid_.sector() << ":" << hgcid_.subsector() << ":" << hgcid_.cell() 
	     // << " " <<  cellPos << std::endl;
	     if (id_.subdetId()!=HGCEE) std::cout << "cell " << cellPos << " does not belong to HGCEE "<< std::endl;
	     sumhits += theSeedHit->energy();
	     if (theSeedHit->energy()>20.*0.000045)sumhits20 += theSeedHit->energy();
	     if (theSeedHit->energy()>50.*0.000045)sumhits50 += theSeedHit->energy();
	     if (theSeedHit->energy()>100.*0.000045)sumhits100 += theSeedHit->energy();
	     if (theSeedHit->energy()>200.*0.000045)sumhits200 += theSeedHit->energy();
	     if (segmented) {
	      if (nevt<10) h_hgcal_shower_seed_em[nevt]->Fill(fabs(cellPos.x()),fabs(cellPos.y()),fabs(cellPos.z()),theSeedHit->energy()*1000.);
            } else {
	      if (nevt<10) h_hgcal_shower_seed_em[nevt]->Fill(fabs(projPos.x()),fabs(projPos.y()),fabs(projPos.z()),theSeedHit->energy()*1000.);
	     }
	   }
	 }

	 // principal component analysis for the seed cluster
	 GlobalPoint pcaShowerPos;
	 GlobalVector pcaShowerDir;
	 pcaShowerAnalysisLogwEle.showerParameters(&scl,pcaShowerPos,pcaShowerDir);
	 if (nevt<10) std::cout << "*** Principal component analysis (standalone class) ****" << std::endl;
	 if (nevt<10) std::cout << "shower average (x,y,z) = " << pcaShowerPos << std::endl;
	 if (nevt<10) std::cout << "shower main axis (x,y,z) = " << pcaShowerDir << std::endl;
	 GlobalVector pcaShowerEigenValues = pcaShowerAnalysisEle.showerEigenValues(&scl);
	 if (nevt<10) std::cout << "shower eigen values (notlogweighted) (x,y,z) = " << pcaShowerEigenValues << std::endl;
 	 GlobalVector pcaShowerSigmas = pcaShowerAnalysisEle.showerSigmas(&scl);
	 if (nevt<10) std::cout << "shower sigmas(notlogweighted) (x,y,z) = " << pcaShowerSigmas << std::endl;
   
	 // initialize shower position and direction in identification object
	 hGCALShowerBasedEmIdentification.setShowerPosition(pcaShowerPos);
	 hGCALShowerBasedEmIdentification.setShowerDirection(pcaShowerDir);
	 // get first position, refined in x and y from the shower direction	 
	 GlobalPoint firstPos = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()));
	 h_hgcal_sclusters_entryposition->Fill(fabs(firstPos.z()));
	 GlobalPoint firstPos_1mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),1.);
	 h_hgcal_sclusters_entryposition_1mip->Fill(fabs(firstPos_1mip.z()));
	 GlobalPoint firstPos_2mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),2.);
	 h_hgcal_sclusters_entryposition_2mip->Fill(fabs(firstPos_2mip.z()));
	 GlobalPoint firstPos_4mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),4.);
	 h_hgcal_sclusters_entryposition_4mip->Fill(fabs(firstPos_4mip.z()));
	 GlobalPoint firstPos_8mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),8.);
	 h_hgcal_sclusters_entryposition_8mip->Fill(fabs(firstPos_8mip.z()));
	 GlobalPoint firstPos_12mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),12.);
	 h_hgcal_sclusters_entryposition_12mip->Fill(fabs(firstPos_12mip.z()));
	 GlobalPoint firstPos_16mip = hGCALShowerBasedEmIdentification.startPosition(&(*scl.seed()),16.);
	 h_hgcal_sclusters_entryposition_16mip->Fill(fabs(firstPos_16mip.z()));
	 
         // fill PCA analysis histograms
	 h_hgcal_sclusters_energyVSlongwidth_fullrange->Fill(log(scl.seed()->energy()),pcaShowerEigenValues.x());
	 h_hgcal_sclusters_entryVSlongwidth_fullrange->Fill(fabs(firstPos.z()),pcaShowerEigenValues.x());
	 double transverse = sqrt(pcaShowerEigenValues.y()*pcaShowerEigenValues.y()+pcaShowerEigenValues.z()*pcaShowerEigenValues.z());
         h_hgcal_sclusters_energyVStranswidth_fullrange->Fill(log(scl.seed()->energy()),transverse);
	 h_hgcal_sclusters_entryVStranswidth_fullrange->Fill(fabs(firstPos.z()),transverse);
         h_hgcal_sclusters_energyVSeigenratio_fullrange->Fill(log(scl.seed()->energy()),transverse/pcaShowerEigenValues.x());
	 h_hgcal_sclusters_entryVSeigenratio_fullrange->Fill(fabs(firstPos.z()),transverse/pcaShowerEigenValues.x());

	 h_hgcal_sclusters_firsteigenvalue->Fill(pcaShowerEigenValues.x());
	 h_hgcal_sclusters_secondeigenvalue->Fill(pcaShowerEigenValues.y());
 	 h_hgcal_sclusters_thirdeigenvalue->Fill(pcaShowerEigenValues.z());
         h_hgcal_sclusters_firsteigenvalueVSeta->Fill(pcaShowerEigenValues.x(),fabs(mcIter->eta()));
         h_hgcal_sclusters_secondeigenvalueVSeta->Fill(pcaShowerEigenValues.y(),fabs(mcIter->eta()));
	 h_hgcal_sclusters_firsteigenvalueVSsecond->Fill(pcaShowerEigenValues.x(),pcaShowerEigenValues.y());
         h_hgcal_sclusters_transeigenvalueVSeta->Fill(transverse,fabs(mcIter->eta()));
         h_hgcal_sclusters_eigenratioVSeta->Fill(transverse,fabs(mcIter->eta()));
	 h_hgcal_sclusters_firstsigma->Fill(pcaShowerSigmas.x());
	 h_hgcal_sclusters_secondsigma->Fill(pcaShowerSigmas.y());
 	 h_hgcal_sclusters_thirdsigma->Fill(pcaShowerSigmas.z());
        
	 double length =  (pcaShowerPos - firstPos).mag();
	 if (nevt<10) std::cout << "shower length " << length << " first pos " << firstPos << std::endl;	 

         // histogram reconstructed seed shower position 
  	 h_hgcal_sclusters_etaPCAMinusEtaTrue->Fill(pcaShowerPos.eta()-mcIter->eta());
         h_hgcal_sclusters_etaPCAMinusEtaTrueVsEta->Fill(pcaShowerPos.eta()-mcIter->eta(),pcaShowerPos.eta());
	 // correct theta/eta position for the z vertex
	 double pcaShowerPosTheta = atan2(pcaShowerPos.perp(),pcaShowerPos.z()-mcIter->vz());
	 double pcaShowerPosEta = -log(tan(pcaShowerPosTheta/2.));
  	 h_hgcal_sclusters_etaPCAMinusEtaTrue_corr->Fill(pcaShowerPosEta-mcIter->eta());
 	 if (scl.clusters().size()==1) h_hgcal_sclusters_etaPCAMinusEtaTrue_golden->Fill(pcaShowerPos.eta()-mcIter->eta());
 	 double deta, dphi, dtheta, dphicorr; 
 	 double mc_charge = mcIter->pdgId() > 0 ? -1. : 1. ;
	 dphi = pcaShowerPos.phi()-mcIter->phi();
	 if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	 h_hgcal_sclusters_phiPCAMinusPhiTrue->Fill(dphi*mc_charge);
         h_hgcal_sclusters_phiPCAMinusPhiTrueVsPhi->Fill(dphi*mc_charge,pcaShowerPos.phi());
	 // correct phi position for bending
	 dphicorr = pcaShowerPos.phi()-rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge()).phi();
	 if (std::abs(dphicorr)>CLHEP::pi) dphicorr = dphicorr < 0? (CLHEP::twopi) + dphicorr : dphicorr - CLHEP::twopi;
	 h_hgcal_sclusters_phiPCAMinusPhiTrue_corr->Fill(dphicorr*mc_charge);
	 if (scl.clusters().size()==1) h_hgcal_sclusters_phiPCAMinusPhiTrue_golden->Fill(dphi*mc_charge);
         h_hgcal_sclusters_etaPCAMinusEtaSC->Fill(pcaShowerPos.eta()-scl.position().eta());
	 dphi = pcaShowerPos.phi()-scl.position().phi();	 
	 if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
         h_hgcal_sclusters_phiPCAMinusPhiSC->Fill(dphi*mc_charge);
         h_hgcal_sclusters_etaSCMinusEtaTrue->Fill(scl.position().eta()-mcIter->eta());
	 // correct theta/eta position for the z vertex
	 GlobalPoint scShowerPos(scl.position().x(),scl.position().y(),scl.position().z());
	 double scShowerPosTheta = atan2(scShowerPos.perp(),scShowerPos.z()-mcIter->vz());
	 double scShowerPosEta = -log(tan(scShowerPosTheta/2.));
         h_hgcal_sclusters_etaSCMinusEtaTrue_corr->Fill(scShowerPosEta-mcIter->eta());
         h_hgcal_sclusters_etaSCMinusEtaTrueVsEta->Fill(scl.position().eta()-mcIter->eta(),scl.position().eta());
	 h_hgcal_sclusters_etaSeedMinusEtaTrueVsEta->Fill(scl.seed()->eta()-mcIter->eta(),mcIter->eta());
	 dphi = scl.position().phi()-mcIter->phi();	 	 
	 if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
         h_hgcal_sclusters_phiSCMinusPhiTrue->Fill(dphi*mc_charge);
         h_hgcal_sclusters_phiSCMinusPhiTrueVsPhi->Fill(dphi*mc_charge,scl.position().phi());
	 // correct phi position for bending;
	 dphicorr = scShowerPos.phi()-rotateMomentum(magField_,genMomentum,scShowerPos,genVertex,mcIter->charge()).phi();
	 if (std::abs(dphicorr)>CLHEP::pi) dphicorr = dphicorr < 0? (CLHEP::twopi) + dphicorr : dphicorr - CLHEP::twopi;
         h_hgcal_sclusters_phiSCMinusPhiTrue_corr->Fill(dphicorr*mc_charge);
         
	 // histogram reconstructed seed shower direction	 
	 std::cout << "PCA shower dir " << pcaShowerDir << "eta " << pcaShowerDir.eta() << " phi " << pcaShowerDir.phi() << std::endl;
	 std::cout << "MC ele dir eta " << mcIter->eta() << " phi " << mcIter->phi() << std::endl;
	 deta = pcaShowerDir.eta()-mcIter->eta();
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrue->Fill(deta);
         if (scl.clusters().size()==1) h_hgcal_sclusters_etadirPCAMinusEtaDirTrue_golden->Fill(deta);
	 dtheta = pcaShowerDir.theta()-mcIter->theta();	 	 
	 if (std::abs(dtheta)>CLHEP::pi) dtheta = dtheta < 0? (CLHEP::twopi) + dtheta : dtheta - CLHEP::twopi;	 
	 h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue->Fill(dtheta);
	 h_hgcal_sclusters_thetadirPCAMinusThetaDirTrueVsEnergy->Fill(mcIter->p(),dtheta);
	 if (scl.clusters().size()==1) h_hgcal_sclusters_thetadirPCAMinusThetaDirTrue_golden->Fill(dtheta);
         // correct for the rotation of the momentum
	 dphi = pcaShowerDir.phi()-rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge()).phi();
	 std::cout << "momentum rotation: old momentum " << genMomentum << " new momentum " <<
	  rotateMomentum(magField_,genMomentum,pcaShowerPos,GlobalPoint(0.,0.,0.),mcIter->charge()) << std::endl;
	 if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	 //std::cout << "deta is " << pcaShowerDir.eta()-mcIter->eta() << " dphi is " << dphi << std::endl;
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrue->Fill(dphi*mc_charge);
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy->Fill(mcIter->p(),dphi*mc_charge);
         if (scl.clusters().size()==1) h_hgcal_sclusters_phidirPCAMinusPhiDirTrue_golden->Fill(dphi*mc_charge);
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueUnsigned->Fill(dphi);
	 // pb diagonales phi>0
	 if (pcaShowerPos.phi()<0.) h_hgcal_sclusters_phidirPCAMinusPhiDirTrueNeg->Fill(dphi);
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEnergy->Fill(dphi,scl.seed()->energy());
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsEta->Fill(dphi,scl.seed()->eta());
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiMC->Fill(dphi,mcIter->phi());
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsPhiPCA->Fill(dphi,pcaShowerPos.phi());
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSCmultiplicity->Fill(dphi,float(scl.clusters().size()));
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedmultiplicity->Fill(dphi,float(scl.seed()->hitsAndFractions().size()));
         h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsSeedfraction->Fill(dphi*mc_charge,scl.seed()->energy()/mcIter->p());
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEnergy->Fill(deta,scl.seed()->energy());
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsPhi->Fill(deta,scl.seed()->phi());
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaMC->Fill(deta,mcIter->eta());
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsEtaPCA->Fill(deta,pcaShowerPos.eta());
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicity->Fill(deta,float(scl.seed()->hitsAndFractions().size()));
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedfraction->Fill(deta,scl.seed()->energy()/mcIter->p());
	 if (mcIter->eta()>0.) h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSeedmultiplicityPosZ->Fill(deta,float(scl.seed()->hitsAndFractions().size()));
	 h_hgcal_sclusters_AbsetadirPCAMinusEtaDirTrueVsSeedmultiplicity->Fill(fabs(deta),float(scl.seed()->hitsAndFractions().size()));
	 h_hgcal_sclusters_etadirPCAMinusEtaDirTrueVsSCmultiplicity->Fill(deta,float(scl.clusters().size()));
	 h_hgcal_sclusters_phidirPCAMinusPhiDirTrueVsphiPCAMinusPhiTrue_corr->Fill(dphi,dphicorr);

         // now generate the fastsim shower for this electron seed cluster
	 if (nevt<10) std::cout << "*** Fastsim shower simulation ****" << std::endl;
         if (nevt<10) std::cout << "HGCAL radiation length " << calohelper_->ecalProperties(true)->radLenIncm() << std::endl;
         if (nevt<10) std::cout << "HGCAL critical energy " << calohelper_->ecalProperties(true)->criticalEnergy() << std::endl;
         if (nevt<10) std::cout << "HGCAL Moliere radius " << calohelper_->ecalProperties(true)->moliereRadius() << std::endl;
	 GlobalVector dir = pcaShowerDir.unit();
         const XYZTLorentzVector momentum(dir.x()*scl.seed()->energy(),dir.y()*scl.seed()->energy(),dir.z()*scl.seed()->energy(),scl.seed()->energy());
	 RawParticle myEle(22,momentum);
	 std::vector<const RawParticle *> thePart;
	 if ( myEle.e() > 0.055 ) thePart.push_back(&myEle);
	 EMShower theShower(random_,aGammaGenerator_,showerparam_,&thePart,NULL,NULL,NULL,false); 
         double meanShower = theShower.getMeanDepth();    
         //double maxShower = theShower.getMaximumOfShower();    
	 double x0 = calohelper_->ecalProperties(true)->radLenIncm();	 
	 //double predictedMax = maxShower * x0;
	 double predictedLength = meanShower * x0;
	 // longitudinal scaling
	 double lny = scl.seed()->energy()/0.00536>1. ? std::log(scl.seed()->energy()/0.00536) : 0.;
	 double meant = showerparam_->meanT(lny)*x0;
	 //double meanalpha = showerparam_->meanAlpha(lny);
	 //double meanlnt = showerparam_->meanLnT(lny);
	 //double sigmalnt = showerparam_->sigmaLnT(lny);
	 //double meanlnalpha = showerparam_->meanLnAlpha(lny);
	 //double sigmalnalpha = showerparam_->sigmaLnAlpha(lny);
	 //double corralphat = showerparam_->correlationAlphaT(lny);
	 double scalelong = 1./meant;

         // histogram parametrised and measured shower quantities
	 h_hgcal_sclusters_predictedTvsy->Fill(meanShower,momentum.e()/calohelper_->ecalProperties(true)->criticalEnergy());
         h_hgcal_sclusters_predictedLnTvsy->Fill(log(meanShower),momentum.e()/calohelper_->ecalProperties(true)->criticalEnergy());
	 h_hgcal_sclusters_predictedLength->Fill(predictedLength,log(momentum.e()));	  
	 h_hgcal_sclusters_predictedLength_fullrange->Fill(predictedLength,log(momentum.e()));	  
	 h_hgcal_sclusters_length->Fill(length,log(scl.seed()->energy()));
	 h_hgcal_sclusters_energyVSlength->Fill(log(scl.seed()->energy()),length);
	 h_hgcal_sclusters_length_fullrange->Fill(length,log(scl.seed()->energy()));

         // define length compatibility cut
	 double lengthCompatibility = hGCALShowerBasedEmIdentification.lengthCompatibility(&(*scl.seed()));
	 h_hgcal_sclusters_lengthCompatibility->Fill(lengthCompatibility);
	 bool cutseedlength = (hGCALShowerBasedEmIdentification.cutLengthCompatibility(&(*scl.seed())));

	 h_hgcal_sclusters_energyVSlength_fullrange->Fill(log(scl.seed()->energy()),length);
	 if (cutseedlength) h_hgcal_sclusters_energyVSlength_cut_fullrange->Fill(log(scl.seed()->energy()),length);
	 h_hgcal_sclusters_entryVSlength_fullrange->Fill(fabs(firstPos.z()),length);
	 if (cutseedlength) h_hgcal_sclusters_entryVSlength_cut_fullrange->Fill(fabs(firstPos.z()),length);

	 // define cut in starting position
	 bool cutseedpos = hGCALShowerBasedEmIdentification.cutStartPosition(&(*scl.seed()));
	 //bool cutseedpostight = (fabs(firstPos.z())<321.0);
	 if (cutseedpos) h_hgcal_sclusters_energyVSlength_hasfirstlayer_fullrange->Fill(log(scl.seed()->energy()),length);
	 
	 // other longitudinal variables
	 double e2530OverEtot = hGCALShowerBasedEmIdentification.E2530OverEtot(&(*scl.seed()));
	 h_hgcal_sclusters_e2530OverEtot->Fill(e2530OverEtot);
	 double e0110OverEtot = hGCALShowerBasedEmIdentification.E0110OverEtot(&(*scl.seed()));
	 h_hgcal_sclusters_e0110OverEtot->Fill(e0110OverEtot);
	 h_hgcal_sclusters_e0110OverEtotVSe2530OverEtot->Fill(e0110OverEtot,e2530OverEtot);

	 // overal longitudinal cut on the seed cluster
	 bool cutseedall = cutseedlength && cutseedpos;
	 if (cutseedall) h_hgcal_sclusters_energyVSlength_hasfirstlayer_cut_fullrange->Fill(log(scl.seed()->energy()),length);

	 h_hgcal_sclusters_subclustersVsEnergy->Fill(float(scl.clusters().size()),mcIter->p());
	 h_hgcal_sclusters_seedEnergyVsEnergy->Fill(float(scl.seed()->energy()), mcIter->p());
	 h_hgcal_sclusters_seedEnergyVsEta->Fill(float(scl.seed()->energy()), mcIter->eta());
	 	 
         // now compute transverse variables
	 double sigmar = 0., sigmarxy = 0., meanr=0., sumnrj=0., sumw = 0., sumw200 = 0.;
	 double sigmaeta = 0., sigmaetaw = 0., sigmaetaw200 = 0.;
	 double sigmaeta20 = 0., sigmaeta50 = 0., sigmaeta100 = 0.;
	 double sumnrj20 = 0., sumnrj50 = 0., sumnrj100 = 0.;
	 double meaneta20 = 0., meaneta50 = 0., meaneta100 = 0.;
	 double sigmaeta_corr = 0., sigmaphi = 0., sigmaphi_corr = 0.;
	 double meaneta = 0., meanetaw = 0., meanphi = 0.;
	 double sigmaeta_pu = 0., sumnrj_pu = 0.; 
	 double sigmaeta_pu_corr = 0.;
	 double sigmart = 0., sigmartaxis = 0.;
	 double sigmartaxismiddle = 0., sumnrjmiddle = 0.;
	 //GlobalVector trueDir = rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge());
	 for (unsigned int ih=0;ih<scl.seed()->hitsAndFractions().size();++ih) {
	   const DetId & id_ = (scl.seed()->hitsAndFractions())[ih].first ;
	   h_hgcal_sclusters_seedfractions_em->Fill((scl.seed()->hitsAndFractions())[ih].second);
	   HGCRecHitCollection::const_iterator theSeedHit = recHits->find(id_);    
	   if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
 	     const HGCEEDetId & hgcid_ = HGCEEDetId((scl.seed()->hitsAndFractions())[ih].first) ;
	     GlobalPoint cellPos = geometry_->getPosition(hgcid_);	
	     int ilayer =  hgcid_.layer();    
	     GlobalPoint projPos(cellPos.x(),cellPos.y(),320.38);
	     GlobalPoint projShowerPos(pcaShowerPos.x(),pcaShowerPos.y(),320.38);
	     GlobalPoint axisPos;
	     GlobalVector radius; GlobalVector longitudinal; GlobalVector transverse;
	     if (segmented) {
	      radius = cellPos - firstPos;
	      // compute distance to shower axis
	      longitudinal =  (radius.dot(pcaShowerDir))*pcaShowerDir.unit()/pcaShowerDir.mag();
	      transverse = radius - longitudinal;
	     } else {
	      radius = projPos - projShowerPos;
	      // compute distance to shower axis
	      longitudinal = GlobalVector(0.,0.,0.);
	      transverse = projPos - projShowerPos;
	     }	     
	     // for all distributions apply 4mip threshold for pileup in SLHC21
	     if (!withPileup_ || theSeedHit->energy()>4.*0.0000551) {
	       if (theSeedHit->energy()>20*0.000051) { 
		sigmaeta20 += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy(); sumnrj20+= theSeedHit->energy();
		meaneta20 += (cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy();
	       }
	       if (theSeedHit->energy()>50*0.000051) { 
		sigmaeta50 += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy(); sumnrj50 += theSeedHit->energy();
		meaneta50 += (cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy();
	       }
	       if (theSeedHit->energy()>100*0.000051) { 
		sigmaeta100 += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy(); sumnrj100 += theSeedHit->energy();
		meaneta100 += (cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy();
	       }
	       if (segmented) {
	        sigmar += transverse.mag2() * theSeedHit->energy();
	        meanr += transverse.mag() * theSeedHit->energy();
	        sumnrj += theSeedHit->energy();
                double deta = (cellPos.eta()-pcaShowerPos.eta());
		meaneta += deta*theSeedHit->energy();
		sigmaeta += deta*deta*theSeedHit->energy();	     
		sigmart += (cellPos.perp()-pcaShowerPos.perp())*(cellPos.perp()-pcaShowerPos.perp()) * theSeedHit->energy();	     
		//sigmart += (cellPos-pcaShowerPos).perp() * (cellPos-pcaShowerPos).perp() * theSeedHit->energy();	     
		// another transverse radius, taken with respect to the shower axis rather than shower barycenter
	        double lambdatoaverage = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
	        axisPos = pcaShowerPos + lambdatoaverage*pcaShowerDir;	 
		sigmartaxis += (cellPos.perp()-axisPos.perp()) * (cellPos.perp()-axisPos.perp()) * theSeedHit->energy();	     				
		if (hgcid_.layer()>7 && hgcid_.layer()<19) {
		  sigmartaxismiddle += (cellPos.perp()-pcaShowerPos.perp())*(cellPos.perp()-pcaShowerPos.perp()) * theSeedHit->energy();	     
	          sumnrjmiddle += theSeedHit->energy();
		}
		//sigmartaxis += (cellPos-axisPos).perp() * (cellPos-axisPos).perp() * theSeedHit->energy();	     				
		//std::cout << "transverse axis pos: average position " << pcaShowerPos << " direction "	<< trueDir << std::endl;
		//std::cout << " new cell pos " << cellPos << " with energy " << theSeedHit->energy() << " axis position " << axisPos<< " trasnsverse axis distance " << (cellPos-axisPos).perp()<< std::endl;		
		// transversal cut, choose simple cut 1.5Rm, then can refine as function of depth
		double moliere = calohelper_->ecalProperties(true)->moliereRadius();
		double rmax = 1.5*moliere;
		//if (nevt<10) std::cout << "transverse.mag() " << transverse.mag() << " rmax " << rmax << std::endl;
		//if (longitudinal.mag() > lmin) {
		if (transverse.mag() < rmax) {
	          sigmaeta_pu += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * theSeedHit->energy();	     
	          sumnrj_pu += theSeedHit->energy();
		}
		double dphi = (cellPos.phi()-pcaShowerPos.phi());
		if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
		meanphi += dphi*theSeedHit->energy();
 		sigmaphi += dphi * dphi * theSeedHit->energy();
		if (theSeedHit->energy()>20*0.000045) { 
		 sigmaetaw += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * log(theSeedHit->energy()/20*0.000045);
		 sumw +=  log(theSeedHit->energy()/20*0.000045);
		} 
		if (theSeedHit->energy()>200*0.000045) { 
		 sigmaetaw200 += (cellPos.eta()-pcaShowerPos.eta())*(cellPos.eta()-pcaShowerPos.eta()) * log(theSeedHit->energy()/20*0.000045);
		 sumw200 +=  log(theSeedHit->energy()/20*0.000045);
		} 
	       } else {
		meaneta += (projPos.eta()-projShowerPos.eta())*theSeedHit->energy();
		sigmaeta += (projPos.eta()-projShowerPos.eta())*(projPos.eta()-projShowerPos.eta()) * theSeedHit->energy();	     
		double dphi = (projPos.phi()-projShowerPos.phi());
		if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
		meanphi += dphi*theSeedHit->energy();
 		sigmaphi += dphi * dphi * theSeedHit->energy();
		if (theSeedHit->energy()>20*0.000045) { 
		 sigmaetaw += (projPos.eta()-projShowerPos.eta())*(projPos.eta()-projShowerPos.eta()) * log(theSeedHit->energy()/20*0.000045);
		 sumw +=  log(theSeedHit->energy()/20*0.000045);
		} 
		if (theSeedHit->energy()>200*0.000045) { 
		 sigmaetaw200 += (projPos.eta()-projShowerPos.eta())*(projPos.eta()-projShowerPos.eta()) * log(theSeedHit->energy()/20*0.000045);
		 sumw200 +=  log(theSeedHit->energy()/20*0.000045);
		} 
	       }
	       // fill profile plots
	       h_hgcal_sclusters_longitudinal->Fill(longitudinal.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_layer->Fill(float(ilayer),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversal->Fill(transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalaxis->Fill(fabs(cellPos.perp()-axisPos.perp()),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_deta_shower->Fill(fabs(deta),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSeta->Fill(transverse.mag(),mcIter->eta(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinal->Fill(longitudinal.mag(),transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinalscaled->Fill(longitudinal.mag()*scalelong,transverse.mag(),theSeedHit->energy()/sumhits);	       
	       h_hgcal_sclusters_transversalVSlongitudinalscaledmeasured->Fill(longitudinal.mag()/length,transverse.mag(),theSeedHit->energy()/sumhits);	       
	       h_hgcal_sclusters_transversalaxisVSlongitudinal->Fill(longitudinal.mag(),fabs(cellPos.perp()-axisPos.perp()),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalaxisVSlongitudinalscaled->Fill(longitudinal.mag()*scalelong,fabs(cellPos.perp()-axisPos.perp()),theSeedHit->energy()/sumhits);	       
	       h_hgcal_sclusters_transversalaxisVSlongitudinalscaledmeasured->Fill(longitudinal.mag()/length,fabs(cellPos.perp()-axisPos.perp()),theSeedHit->energy()/sumhits);	       
	       // same with cuts
	       if (theSeedHit->energy()>20.*0.000045) {
	       h_hgcal_sclusters_longitudinal_cut20->Fill(longitudinal.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversal_cut20->Fill(transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinal_cut20->Fill(longitudinal.mag(),transverse.mag(),theSeedHit->energy()/sumhits20);
	       }
	       if (theSeedHit->energy()>50.*0.000045) {
	       h_hgcal_sclusters_longitudinal_cut50->Fill(longitudinal.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversal_cut50->Fill(transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinal_cut50->Fill(longitudinal.mag(),transverse.mag(),theSeedHit->energy()/sumhits50);
	       }
	       if (theSeedHit->energy()>100.*0.000045) {
	       h_hgcal_sclusters_longitudinal_cut100->Fill(longitudinal.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversal_cut100->Fill(transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinal_cut100->Fill(longitudinal.mag(),transverse.mag(),theSeedHit->energy()/sumhits100);
	       }
	       if (theSeedHit->energy()>200.*0.000045) {
	       h_hgcal_sclusters_longitudinal_cut200->Fill(longitudinal.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversal_cut200->Fill(transverse.mag(),theSeedHit->energy()/sumhits);
	       h_hgcal_sclusters_transversalVSlongitudinal_cut200->Fill(longitudinal.mag(),transverse.mag(),theSeedHit->energy()/sumhits200);
	       }
	     }
	   }
	 }

         sigmar /= sumnrj;
         sigmart /= sumnrj;
         sigmartaxis /= sumnrj;
         sigmartaxismiddle /= sumnrjmiddle;
	 sigmarxy /= sumnrj;
         sigmaeta20 /= sumnrj20;
         sigmaeta50 /= sumnrj50;
         sigmaeta100 /= sumnrj100;
         sigmaeta /= sumnrj;
         sigmaeta_pu /= sumnrj_pu;
         sigmaphi /= sumnrj;
         sigmaetaw /= sumw;
         sigmaetaw200 /= sumw200;
         meanr /= sumnrj;
         meaneta /= sumnrj;
         meanphi /= sumnrj;
         meanetaw /= sumw;
         meaneta20 /= sumnrj20;
         meaneta50 /= sumnrj50;
         meaneta100 /= sumnrj100;
         double sigmarnorm = sqrt(sigmar-meanr*meanr);
         double sigmaetanorm20 = sqrt(sigmaeta20-meaneta20*meaneta20);
         double sigmaetanorm50 = sqrt(sigmaeta50-meaneta50*meaneta50);
         double sigmaetanorm100 = sqrt(sigmaeta100-meaneta100*meaneta100);
	 double sigmaetawnorm = sqrt(sigmaetaw-meanetaw*meanetaw);
	 double sigmaetanorm = sqrt(sigmaeta-meaneta*meaneta);
	 double sigmaphinorm = sqrt(sigmaphi-meanphi*meanphi);
	 sigmart = sqrt(sigmart);
	 sigmartaxis = sqrt(sigmartaxis);
	 sigmartaxismiddle = sqrt(sigmartaxismiddle);
	 sigmaeta = sqrt(sigmaeta);
	 sigmaeta_pu = sqrt(sigmaeta_pu);
	 sigmaphi = sqrt(sigmaphi);
	 sigmaetaw = sqrt(sigmaetaw);
	 sigmaetaw200 = sqrt(sigmaetaw200);
	 
	 sigmaeta_corr = hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()));
	 //// same for pu_corr, the standalone class takes care of PU correction
	 sigmaeta_pu_corr = sigmaeta_corr;

 	 double sigmaeta_pu_corr10 = hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),10);
        
  	 sigmaphi_corr = hGCALShowerBasedEmIdentification.sigmaphiphi(&(*scl.seed()));
	 	 
         // now define sigmaetaeta cut
	 //bool cutsetaeta = (sigmaeta_corr<0.0055);
	 //if (withPileup_) cutsetaeta = (sigmaeta_pu_corr<0.00475);
         bool cutsetaeta = hGCALShowerBasedEmIdentification.cutSigmaetaeta(&(*scl.seed()));
	 
	 //double feta, feta_0, sigmart_corr;
	 double feta, feta_0, sigmart_corr;
	 feta = 4.94156-2.6929*fabs(scl.seed()->eta())+0.453452*fabs(scl.seed()->eta())*fabs(scl.seed()->eta());
	 feta_0 = 4.94156-2.6929*1.5+0.453452*1.5*1.5;
	 sigmart_corr = sigmart * feta_0 / feta;

	 double sigmartaxis_corr;
//          // to be removed with standalone class
// 	 if (fabs(scl.seed()->eta())<2.6)
//           feta = 4.13856 - 2.52796*fabs(scl.seed()->eta()) + 0.487409*fabs(scl.seed()->eta())*fabs(scl.seed()->eta());
//          else feta = 4.13856 - 2.52796*2.6 + 0.487409*2.6*2.6;
//          feta_0 = 4.13856 - 2.52796*1.5 + 0.487409*1.5*1.5;	 
// 	 sigmartaxis_corr = sigmartaxis * feta_0 / feta;
//         // to be removed with standalone class
	 
	 // standalone class
	 sigmartaxis_corr = hGCALShowerBasedEmIdentification.sigmartrt(&(*scl.seed()));
	 // standalone class
	 
	 // assume same correction for sigmartmiddle	 
	 double sigmartaxismiddle_corr;
	 sigmartaxismiddle_corr = sigmartaxismiddle * feta_0 / feta;
	 
 	 sigmaphi_corr = hGCALShowerBasedEmIdentification.sigmaphiphi(&(*scl.seed()));
	 	 
	 // transverse quantities and PCA eigenvalues
	 h_hgcal_sclusters_energyVSsigmaradius_fullrange->Fill(log(scl.seed()->energy()),sigmar);
	 h_hgcal_sclusters_entryVSsigmaradius_fullrange->Fill(fabs(firstPos.z()),sigmar);
	 h_hgcal_sclusters_energyVSmeanradius_fullrange->Fill(log(scl.seed()->energy()),meanr);
	 h_hgcal_sclusters_entryVSmeanradius_fullrange->Fill(fabs(firstPos.z()),meanr);
	 h_hgcal_sclusters_sigmaradiusnorm_fullrange->Fill(sigmarnorm);
	 h_hgcal_sclusters_sigmaetanorm20_fullrange->Fill(sigmaetanorm20);
	 h_hgcal_sclusters_sigmaetanorm50_fullrange->Fill(sigmaetanorm50);
	 h_hgcal_sclusters_sigmaetanorm100_fullrange->Fill(sigmaetanorm100);
	 h_hgcal_sclusters_sigmaradius_fullrange->Fill(sigmar);
	 h_hgcal_sclusters_sigmaeta20_fullrange->Fill(sigmaeta20);
	 h_hgcal_sclusters_sigmaeta50_fullrange->Fill(sigmaeta50);
	 h_hgcal_sclusters_sigmaeta100_fullrange->Fill(sigmaeta100);
	 h_hgcal_sclusters_sigmaradiusVSeta_fullrange->Fill(sigmar,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaradiusVSphi_fullrange->Fill(sigmar,mcIter->phi());
	 h_hgcal_sclusters_sigmatransverseradius_fullrange->Fill(sigmart);
	 h_hgcal_sclusters_sigmatransverseradiusVSeta_fullrange->Fill(sigmart,fabs(mcIter->eta()));
	 h_hgcal_sclusters_etaVSsigmatransverseradius_fullrange->Fill(fabs(mcIter->eta()),sigmart);
	 h_hgcal_sclusters_sigmatransverseradius_corr_fullrange->Fill(sigmart_corr);
	 h_hgcal_sclusters_sigmatransverseradiusaxis_fullrange->Fill(sigmartaxis);
	 h_hgcal_sclusters_sigmatransverseradiusaxismiddle_fullrange->Fill(sigmartaxismiddle);
	 h_hgcal_sclusters_sigmatransverseradiusaxismiddle_corr_fullrange->Fill(sigmartaxismiddle_corr);
	 h_hgcal_sclusters_sigmatransverseradiusaxisVSeta_fullrange->Fill(sigmartaxis,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmatransverseradiusaxisVSphi_fullrange->Fill(sigmartaxis,mcIter->phi());
	 h_hgcal_sclusters_sigmatransverseradiusaxisVSlength_fullrange->Fill(sigmartaxis,length);
	 h_hgcal_sclusters_etaVSsigmatransverseradiusaxis_fullrange->Fill(fabs(mcIter->eta()),sigmartaxis);
	 h_hgcal_sclusters_sigmatransverseradiusaxis_corr_fullrange->Fill(sigmartaxis_corr);
	 h_hgcal_sclusters_sigmaeta_fullrange->Fill(sigmaeta);
	 h_hgcal_sclusters_sigmaeta_pu_fullrange->Fill(sigmaeta_pu);
	 h_hgcal_sclusters_sigmaetanorm_fullrange->Fill(sigmaetanorm);
	 h_hgcal_sclusters_sigmaetaVSeta_fullrange->Fill(sigmaeta,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaeta_puVSeta_fullrange->Fill(sigmaeta_pu,fabs(mcIter->eta()));
	 h_hgcal_sclusters_etaVSsigmaeta_fullrange->Fill(fabs(mcIter->eta()),sigmaeta);
	 h_hgcal_sclusters_sigmaphi_fullrange->Fill(sigmaphi);
	 h_hgcal_sclusters_sigmaphinorm_fullrange->Fill(sigmaphinorm);
	 h_hgcal_sclusters_sigmaphiVSeta_fullrange->Fill(sigmaphi,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaphiVSphi_fullrange->Fill(sigmaphi,fabs(mcIter->phi()));
	 h_hgcal_sclusters_sigmaeta_corr_fullrange->Fill(sigmaeta_corr);
	 h_hgcal_sclusters_sigmaeta_pu_corr_fullrange->Fill(sigmaeta_pu_corr);
	 h_hgcal_sclusters_sigmaetaVSpt_pu_corr_fullrange->Fill(sigmaeta_pu_corr,mcIter->pt());
	 h_hgcal_sclusters_sigmaphi_corr_fullrange->Fill(sigmaphi_corr);
	 dphicorr = pcaShowerPos.phi()-rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge()).phi();
	 if (std::abs(dphicorr)>CLHEP::pi) dphicorr = dphicorr < 0? (CLHEP::twopi) + dphicorr : dphicorr - CLHEP::twopi;
	 h_hgcal_sclusters_sigmaphicorrVSphiminusphitrue->Fill(sigmaphi_corr,dphicorr*mcIter->charge());
	 h_hgcal_sclusters_sigmaetaVSeta_corr_fullrange->Fill(sigmaeta_corr,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaeta_puVSeta_corr_fullrange->Fill(sigmaeta_pu_corr,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaeta_puVSEt_corr_fullrange->Fill(sigmaeta_pu_corr,sclet);
	 h_hgcal_sclusters_sigmaeta_puVSseedfraction_corr_fullrange->Fill(sigmaeta_pu_corr,scl.seed()->energy()/scl.energy());
	 h_hgcal_sclusters_sigmaeta10_puVSEt_corr_fullrange->Fill(sigmaeta_pu_corr10,sclet);
	 h_hgcal_sclusters_sigmaphiVSeta_corr_fullrange->Fill(sigmaphi_corr,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaphiVSEt_corr_fullrange->Fill(sigmaphi_corr,sclet);
	 h_hgcal_sclusters_sigmaetaw_fullrange->Fill(sigmaetaw);
	 h_hgcal_sclusters_sigmaetawnorm_fullrange->Fill(sigmaetawnorm);
	 h_hgcal_sclusters_sigmaetawVSeta_fullrange->Fill(sigmaetaw,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaetawnormVSeta_fullrange->Fill(sigmaetawnorm,fabs(mcIter->eta()));
	 h_hgcal_sclusters_sigmaetaw200_fullrange->Fill(sigmaetaw200);
	         
	 for (int i=1; i<31; i++)
	  h_hgcal_sclusters_sigmaetaVSlayer->Fill(float(i),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),i));
	 for (int i=1; i<31; i++) {
	   double depth = hGCALShowerBasedEmIdentification.depth(&(*scl.seed()),i);
	   h_hgcal_sclusters_sigmaetaVSlayer_norm->Fill(depth/length,hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),i));
	 } 
	 h_hgcal_sclusters_sigmaeta_1->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),1));
	 h_hgcal_sclusters_sigmaeta_2->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),2));
	 h_hgcal_sclusters_sigmaeta_3->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),3));
	 h_hgcal_sclusters_sigmaeta_4->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),4));
	 h_hgcal_sclusters_sigmaeta_5->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),5));
	 h_hgcal_sclusters_sigmaeta_6->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),6));
	 h_hgcal_sclusters_sigmaeta_7->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),7));
	 h_hgcal_sclusters_sigmaeta_8->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),8));
	 h_hgcal_sclusters_sigmaeta_9->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),9));
	 h_hgcal_sclusters_sigmaeta_10->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),10));
	 h_hgcal_sclusters_sigmaeta_11->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),11));
	 h_hgcal_sclusters_sigmaeta_12->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),12));
	 h_hgcal_sclusters_sigmaeta_13->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),13));
	 h_hgcal_sclusters_sigmaeta_14->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),14));
	 h_hgcal_sclusters_sigmaeta_15->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),15));
	 h_hgcal_sclusters_sigmaeta_16->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),16));
	 h_hgcal_sclusters_sigmaeta_17->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),17));
	 h_hgcal_sclusters_sigmaeta_18->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),18));
	 h_hgcal_sclusters_sigmaeta_19->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),19));
	 h_hgcal_sclusters_sigmaeta_20->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),20));
	 h_hgcal_sclusters_sigmaeta_21->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),21));
	 h_hgcal_sclusters_sigmaeta_22->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),22));
	 h_hgcal_sclusters_sigmaeta_23->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),23));
	 h_hgcal_sclusters_sigmaeta_24->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),24));
	 h_hgcal_sclusters_sigmaeta_25->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),25));
	 h_hgcal_sclusters_sigmaeta_26->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),26));
	 h_hgcal_sclusters_sigmaeta_27->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),27));
	 h_hgcal_sclusters_sigmaeta_28->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),28));
	 h_hgcal_sclusters_sigmaeta_29->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),29));
	 h_hgcal_sclusters_sigmaeta_30->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),30));

	 double sigmaetaeta9 = hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),9);
	 if (sigmaetaeta9>0.0015 && sigmaetaeta9<0.0055) {
	   //std::cout << "hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),1)" << hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),1) << std::endl;
	   for (int i=1; i<31; i++)
	    h_hgcal_sclusters_sigmaetaVSlayer_cutlayer9->Fill(float(i),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),i));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_1->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),1));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_2->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),2));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_3->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),3));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_4->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),4));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_5->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),5));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_6->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),6));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_7->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),7));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_8->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),8));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_9->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),9));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_10->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),10));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_11->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),11));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_12->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),12));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_13->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),13));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_14->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),14));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_15->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),15));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_16->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),16));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_17->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),17));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_18->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),18));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_19->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),19));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_20->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),20));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_21->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),21));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_22->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),22));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_23->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),23));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_24->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),24));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_25->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),25));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_26->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),26));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_27->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),27));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_28->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),28));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_29->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),29));
	   h_hgcal_sclusters_sigmaeta_cutlayer9_30->Fill(hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),30));
	 }
	 h_hgcal_sclusters_etaVSsigmaeta_1->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),1));
	 h_hgcal_sclusters_etaVSsigmaeta_2->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),2));
	 h_hgcal_sclusters_etaVSsigmaeta_3->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),3));
	 h_hgcal_sclusters_etaVSsigmaeta_4->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),4));
	 h_hgcal_sclusters_etaVSsigmaeta_5->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),5));
	 h_hgcal_sclusters_etaVSsigmaeta_6->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),6));
	 h_hgcal_sclusters_etaVSsigmaeta_7->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),7));
	 h_hgcal_sclusters_etaVSsigmaeta_8->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),8));
	 h_hgcal_sclusters_etaVSsigmaeta_9->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),9));
	 h_hgcal_sclusters_etaVSsigmaeta_10->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),10));
	 h_hgcal_sclusters_etaVSsigmaeta_11->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),11));
	 h_hgcal_sclusters_etaVSsigmaeta_12->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),12));
	 h_hgcal_sclusters_etaVSsigmaeta_13->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),13));
	 h_hgcal_sclusters_etaVSsigmaeta_14->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),14));
	 h_hgcal_sclusters_etaVSsigmaeta_15->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),15));
	 h_hgcal_sclusters_etaVSsigmaeta_16->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),16));
	 h_hgcal_sclusters_etaVSsigmaeta_17->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),17));
	 h_hgcal_sclusters_etaVSsigmaeta_18->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),18));
	 h_hgcal_sclusters_etaVSsigmaeta_19->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),19));
	 h_hgcal_sclusters_etaVSsigmaeta_20->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),20));
	 h_hgcal_sclusters_etaVSsigmaeta_21->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),21));
	 h_hgcal_sclusters_etaVSsigmaeta_22->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),22));
	 h_hgcal_sclusters_etaVSsigmaeta_23->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),23));
	 h_hgcal_sclusters_etaVSsigmaeta_24->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),24));
	 h_hgcal_sclusters_etaVSsigmaeta_25->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),25));
	 h_hgcal_sclusters_etaVSsigmaeta_26->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),26));
	 h_hgcal_sclusters_etaVSsigmaeta_27->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),27));
	 h_hgcal_sclusters_etaVSsigmaeta_28->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),28));
	 h_hgcal_sclusters_etaVSsigmaeta_29->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),29));
	 h_hgcal_sclusters_etaVSsigmaeta_30->Fill(fabs(mcIter->eta()),hGCALShowerBasedEmIdentification.sigmaetaeta(&(*scl.seed()),30));
	 	 
	 if (cutseedlength) h_hgcal_sclusters_energyVSeigenratio_cut_fullrange->Fill(log(scl.seed()->energy()),transverse/pcaShowerEigenValues.x());
         if (cutseedpos) h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_fullrange->Fill(log(scl.seed()->energy()),transverse/pcaShowerEigenValues.x());
         if (cutseedall) h_hgcal_sclusters_energyVSeigenratio_hasfirstlayer_cut_fullrange->Fill(log(scl.seed()->energy()),transverse/pcaShowerEigenValues.x());
	 
 	 // H/E selection
         double hoverem = 0., hoverem1 = 0., hoverem2 = 0.;
         hoverem = hGCALShowerBasedEmIdentification.hadOverEm(&(*scl.seed()),"all");
         hoverem1 = hGCALShowerBasedEmIdentification.hadOverEm(&(*scl.seed()),"first");
         hoverem2 = hGCALShowerBasedEmIdentification.hadOverEm(&(*scl.seed()),"last");

 	 h_hgcal_sclusters_hoverem->Fill(hoverem); // defual cone size is 0.05
 	 h_hgcal_sclusters_hoverem_cone01->Fill(hoverem,0.1);
 	 h_hgcal_sclusters_hoveremVSeta->Fill(hoverem,fabs(mcIter->eta()));
 	 h_hgcal_sclusters_hoveremVSphi->Fill(hoverem,mcIter->phi());
 	 h_hgcal_sclusters_entryVShoverem->Fill(fabs(firstPos.z()),hoverem);
 	 h_hgcal_sclusters_expectedlengthVShoverem->Fill(predictedLength,hoverem);
 	 h_hgcal_sclusters_hoverem1->Fill(hoverem1);
 	 h_hgcal_sclusters_hoverem1VSeta->Fill(hoverem1,fabs(mcIter->eta()));
 	 h_hgcal_sclusters_entryVShoverem1->Fill(fabs(firstPos.z()),hoverem1);
 	 h_hgcal_sclusters_expectedlengthVShoverem1->Fill(predictedLength,hoverem1);
 	 h_hgcal_sclusters_hoverem2->Fill(hoverem2);
 	 h_hgcal_sclusters_hoverem2VSeta->Fill(hoverem2,fabs(mcIter->eta()));
 	 h_hgcal_sclusters_entryVShoverem2->Fill(fabs(firstPos.z()),hoverem2);
 	 h_hgcal_sclusters_expectedlengthVShoverem2->Fill(predictedLength,hoverem2);
	 
	 // now define H/E cut
	 bool cuthadem =  hGCALShowerBasedEmIdentification.cutHadOverEm(&(*scl.seed()));
	 
//          // to be removed with standalone class
// 	 // Now a longitudinal fit
// 	 double lmax = 30., rmax = 10.; int nbinz=20, nbinr=10;
// 	 double rtaxis = 0.;
// 	 TH1F *h_hgcal_sclusters_longitudinal_shower = new
// 	  TH1F("h_hgcal_sclusters_longitudinal_shower","hgcal seed cluster shower longitudinal profile",nbinz,0.,lmax);
// 	 h_hgcal_sclusters_longitudinal_shower->Sumw2();
// 	 TH1F *h_hgcal_sclusters_transversal_shower = new 
// 	  TH1F("h_hgcal_sclusters_transversal_shower","hgcal seed cluster shower transversal profile",nbinr,0.,rmax);
//          h_hgcal_sclusters_transversal_shower->Sumw2();
// 	 TH2F *h_hgcal_sclusters_3Dprofile_shower = new 
// 	  TH2F("h_hgcal_sclusters_3D_shower","hgcal seed cluster shower 3D profile",nbinz,0.,lmax,nbinr,0.,rmax);
// 	 h_hgcal_sclusters_3Dprofile_shower->Sumw2();
// 	 TH1F *h_hgcal_sclusters_transversal_shower_first = new 
// 	  TH1F("h_hgcal_sclusters_transversal_shower_first","hgcal seed cluster shower transversal profile first layers",nbinr,0.,rmax);
// 	 h_hgcal_sclusters_transversal_shower_first->Sumw2();
// 	 TH1F *h_hgcal_sclusters_transversal_shower_last = new 
// 	  TH1F("h_hgcal_sclusters_transversal_shower_last","hgcal seed cluster shower transversal profile last layers",nbinr,0.,rmax);
// 	 h_hgcal_sclusters_transversal_shower_last->Sumw2();
// 	 sumhits = 0;
// 	 double sumhitsfirst = 0., sumhitslast = 0.;
//          // new calib as of SLHC21
// 	 double clus_eta = scl.seed()->eta();
// 	 //const double _coef_a = 80.0837, _coef_b = -107.229, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
// 	 const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
//          const double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
// 	 double mip = 0.0000551;
// 	 double scale = mip*corr;
// 	 double weight[30] =
// //	    {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,1.19,1.19,1.19,1.19,1.19,1.19,1.19,1.19,1.19};
// 	  {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
// 	 // take initial pos as the intercept of the shower dir with z=320 rather than firstpos
//          //GlobalVector origMom = rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge());
// 	 double lambdaorig = (320.-pcaShowerPos.z())/pcaShowerDir.z();
// 	 //double lambdaorig = (320.-pcaShowerPos.z())/origMom.z();
// 	 if (pcaShowerPos.z()<0.) lambdaorig = (-320.-pcaShowerPos.z())/pcaShowerDir.z();
// 	 //if (pcaShowerPos.z()<0.) lambdaorig = (-320.-pcaShowerPos.z())/origMom.z();
// 	 GlobalPoint origPos = pcaShowerPos + lambdaorig*pcaShowerDir;	 
// 	 //GlobalPoint origPos = pcaShowerPos + lambdaorig*origMom;	 
// 	 for (unsigned int ih=0;ih<scl.seed()->hitsAndFractions().size();++ih) {
// 	   const DetId & id_ = (scl.seed()->hitsAndFractions())[ih].first ;
// 	   h_hgcal_sclusters_seedfractions_em->Fill((scl.seed()->hitsAndFractions())[ih].second);
// 	   HGCRecHitCollection::const_iterator theSeedHit = recHits->find(id_);    
// 	   if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
//  	     const HGCEEDetId & hgcid_ = HGCEEDetId((scl.seed()->hitsAndFractions())[ih].first) ;
// 	     GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
//              // recalibrate per layer
//              const int layer = hgcid_.layer();
// 	     GlobalVector radius; GlobalVector longitudinal; GlobalVector transverse;
// 	     radius = cellPos - origPos;
// 	     longitudinal =  (radius.dot(pcaShowerDir))*pcaShowerDir.unit()/pcaShowerDir.mag();
// 	     //longitudinal =  (radius.dot(origMom))*origMom.unit()/origMom.mag();
// 	     transverse = radius - longitudinal;
// 	     double lambdatoaverage = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
// 	     GlobalPoint axisPos = pcaShowerPos + lambdatoaverage*pcaShowerDir;	 
//              rtaxis = (cellPos.perp()-axisPos.perp()) ;	     				
//  	     // fill profile plots
// 	     if (!withPileup_ || theSeedHit->energy()>4.*0.0000551) {
// 	     h_hgcal_sclusters_longitudinal_shower->Fill(longitudinal.mag(),theSeedHit->energy()*weight[layer-1]/scale);
// 	     //h_hgcal_sclusters_transversal_shower->Fill(transverse.mag(),theSeedHit->energy()*weight[layer-1]/scale);
// 	     h_hgcal_sclusters_transversal_shower->Fill(rtaxis,theSeedHit->energy()*weight[layer-1]/scale);
// 	     h_hgcal_sclusters_transversalaxis_calib->Fill(rtaxis,theSeedHit->energy()*weight[layer-1]/scale);
// 	     h_hgcal_sclusters_3Dprofile_shower->Fill(longitudinal.mag(),rtaxis,theSeedHit->energy()*weight[layer-1]/scale);
// 	     if (longitudinal.mag()<lmax && transverse.mag()<rmax) sumhits += theSeedHit->energy()*weight[layer-1]/scale;
// 	     if (layer<10) {
// 	       h_hgcal_sclusters_transversal_shower_first->Fill(transverse.mag(),theSeedHit->energy()*weight[layer-1]/scale);
// 	       sumhitsfirst += theSeedHit->energy()*weight[layer-1]/scale;	       
// 	     }
// 	     if (layer>20) {
// 	       h_hgcal_sclusters_transversal_shower_last->Fill(transverse.mag(),theSeedHit->energy()*weight[layer-1]/scale);
// 	       sumhitslast += theSeedHit->energy()*weight[layer-1]/scale;	       
// 	     }
// 	     }
// 	   }
// 	 }
// 
// 	 TF1 *gammadist = new TF1("gammadist",gamma,0.,lmax,3);
// 	 // shower parametrization results
// 	 //double lny = scl.seed()->energy()/0.00536>1. ? std::log(scl.seed()->energy()/0.00536) : 0.;
// 	 //double alpha = -0.0433+0.54*lny; 
// 	 double alpha = meanalpha; 
// 	 //std::cout << "alpha " << alpha << " meanalpha " << meanalpha << std::endl;
// 	 //double tmax = -1.396+1.007*lny; 
// 	 double tmax = meant; 
// 	 //std::cout << "tmax " << tmax << " meant " << meant << std::endl;
// 	 //double sigmaalpha = alpha/(-0.8442+0.7904*lny);
// 	 double sigmaalpha = alpha*sigmalnalpha;
// 	 if (sigmaalpha<0.) sigmaalpha = 1.;
// 	 //double sigmatmax = tmax/(-2.506+1.245*lny);
// 	 double sigmatmax = tmax*sigmalnt;
// 	 if (sigmatmax<0.) sigmatmax = 1.;
// 	 double beta = (alpha -1.) / tmax;
// 	 // beta of the gammaDist is 1/beta !
// 	 double invbeta = 1./beta; 
// // 	 double sigmabeta = sigmaalpha*sigmaalpha/tmax + (alpha-1.)*(sigmatmax*sigmatmax)/(tmax*tmax);
//          double sigmainvbeta = pow(((1./(alpha-1.))*sigmaalpha),2) + pow(((tmax/pow((alpha-1.),2))*sigmaalpha),2);
//  	 sigmainvbeta = sqrt(sigmainvbeta);
// 	 //tmax*=x0; beta/=x0; sigmatmax*=x0; invbeta*=x0; sigmainvbeta*=x0;
//  	 if (nevt<10) std::cout << " *** Fitting the longitudinal profile *** " << std::endl;
//  	 if (nevt<10) std::cout << "alpha, invbeta " << alpha << " " << invbeta << std::endl;
//  	 if (nevt<10) std::cout << "sigmaalpha, sigmainvbeta " << sigmaalpha << " " << sigmainvbeta << std::endl;
// 	 //tmax*=x0; beta/=x0; sigmatmax*=x0; invbeta*=x0; sigmainvbeta*=x0;
// 
// 	 ///////////////////////////////////////
// 	 // now fit the longitudinal distribution
//          gammadist->SetParNames("Normalization","alpha","invbeta")	; 
//  	 gammadist->SetParameters(scl.seed()->energy(),alpha,invbeta); 
// 	 //h_hgcal_sclusters_longitudinal_shower->Fit("gammadist","","",0.,25.);
// 	 h_hgcal_sclusters_longitudinal_shower->Fit("gammadist","","",0.,30.);
// 	 ///////////////////////////////////////
// 	 	 
// 	 double fitparams[3] = {0.,0.,0.};
// 	 gammadist->GetParameters(fitparams);
// 	 double chi2 =gammadist->GetChisquare();
// 	 double ndf = gammadist->GetNDF();
// 	 double fittedalpha = fitparams[1];
// 	 double fittedinvbeta = fitparams[2];
//          // to be removed with standalone class

	 // standalone class
	 HGCALFitResults unbounded = hGCALShowerBasedEmIdentification.longitudinalFit(&(*scl.seed()),false, false);
	 double chi2 =  unbounded.chi2();
	 double ndf =  unbounded.ndf();
	 double fittedalpha = unbounded.alpha();
	 double fittedinvbeta = unbounded.invbeta();
	 //double fittednormalization = unbounded.normalization();
	 // standalone class
	 
// 	 if (nevt<10) std::cout << "Fitting the longitudinal profile chi2, ndf " << chi2 << ", " << ndf << std::endl;
// 	 if (nevt<10) {   
// 	   histfile_->cd();
//            h_hgcal_sclusters_transversal_shower->Write();
// 	   h_hgcal_sclusters_3Dprofile_shower->Write();
// 	 }

	 // plot fit result
	 h_hgcal_sclusters_longitudinal_fit_chi2->Fill(chi2/ndf);
	 h_hgcal_sclusters_longitudinal_fit_chi2VSeta->Fill(chi2/ndf, fabs(mcIter->eta()));
	 h_hgcal_sclusters_longitudinal_fit_alphaVSbeta->Fill(fittedalpha,1./fittedinvbeta);
	 h_hgcal_sclusters_longitudinal_fit_alphaVSinvbeta->Fill(fittedalpha,fittedinvbeta);
	 h_hgcal_sclusters_longitudinal_fit_alphaVSenergy->Fill(std::log(scl.seed()->energy()),fittedalpha);
	 h_hgcal_sclusters_longitudinal_fit_invbetaVSenergy->Fill(std::log(scl.seed()->energy()),fittedinvbeta);
	 h_hgcal_sclusters_longitudinal_fit_betaVSenergy->Fill(std::log(scl.seed()->energy()),1./fittedinvbeta);
	 h_hgcal_sclusters_longitudinal_fit_chi2VSseedfraction->Fill(chi2/ndf, scl.seed()->energy()/mcIter->p());
	 h_hgcal_sclusters_longitudinal_fit_chi2VSeoveretrue->Fill(chi2/ndf, scl.energy()/mcIter->p());
	 double dphidir = pcaShowerDir.phi()-rotateMomentum(magField_,genMomentum,pcaShowerPos,genVertex,mcIter->charge()).phi();
	 if (std::abs(dphidir)>CLHEP::pi) dphidir = dphidir < 0? (CLHEP::twopi) + dphidir : dphidir - CLHEP::twopi;
	 h_hgcal_sclusters_longitudinal_fit_chi2VSdphidir->Fill(chi2/ndf,dphi*mc_charge );
	 
// 	 // now integrate the fitted function from 30 to 50.7cm (4 first HCAL layers, to be refined)
// 	 std::cout << "fitparams[0],fitparams[1],fitparams[2] " << fittednormalization << " " << fittedalpha << " " <<fittedinvbeta << std::endl;
// 	 gammadist->SetParameters(fittednormalization,fittedalpha,fittedinvbeta); 
// 	 // from end of HGCEE rather: l = 30./costheta
// 	 double expectedh1 = gammadist->Integral(lmax/std::abs(std::cos(2.*atan(exp(-scl.seed()->eta())))),50.7);	 
//          // compute em integral to normalize
// 	 double expectedem = gammadist->Integral(0.,lmax/std::abs(std::cos(2.*atan(exp(-scl.seed()->eta())))));	 
// 	 h_hgcal_sclusters_longitudinal_fit_leakage->Fill(expectedh1*scl.seed()->energy()/expectedem);
// 	 h_hgcal_sclusters_longitudinal_fit_leakageVShoverem->Fill(expectedh1/expectedem, hoverem1);
// 	 if (cutseedpostight) h_hgcal_sclusters_longitudinal_fit_leakage_cutseedpos->Fill(expectedh1*scl.seed()->energy()/expectedem);
// 	 if (cutseedpostight) h_hgcal_sclusters_longitudinal_fit_leakageVShoverem_cutseedpos->Fill(expectedh1/expectedem, hoverem1);
	 
//          // to be removed with standalone class
// 	 // now build a discriminant based on statistical comparison with an histo
//  	 // inject here measured average length = alpha/beta
// 	 double alphascale = 1.03;
// 	 //double alphascale = 1.0;
// 	 //if (chi2/ndf>5.) { fittedlength = length; alphascale = 1.; }
// 	 gammadist->SetParameters(1.,(length/invbeta)*alphascale,invbeta); 
// 	 //// yet another protection
// 	 TH1F *h_hgcal_sclusters_expected_longitudinal_shower = new
// 	  TH1F("h_hgcal_sclusters_expected_longitudinal_shower","hgcal seed cluster shower longitudinal profile",nbinz,0.,lmax);
// 	 h_hgcal_sclusters_expected_longitudinal_shower->Sumw2();
// 	 if (nevt<10) std::cout << "Filling an histogram with the gamma function " << std::endl;
// 	 if (gammadist->Integral(0.,lmax)!=0.) h_hgcal_sclusters_expected_longitudinal_shower->FillRandom("gammadist");
// 	 // normalize prediction
// 	 if (h_hgcal_sclusters_expected_longitudinal_shower->GetEntries()!=0.) h_hgcal_sclusters_expected_longitudinal_shower->Scale(sumhits/h_hgcal_sclusters_expected_longitudinal_shower->GetEntries());
// 	 if (nevt<10) {
// 	   histfile_->cd();
// 	   h_hgcal_sclusters_expected_longitudinal_shower->Write();
// 	 }  
// 	 // get kolmogorov result
// 	 double kolmogorov = 0.; 
// 	 if (h_hgcal_sclusters_longitudinal_shower->Integral()!=0. && h_hgcal_sclusters_expected_longitudinal_shower->Integral()!=0.) 
// 	  kolmogorov = h_hgcal_sclusters_longitudinal_shower->KolmogorovTest(h_hgcal_sclusters_expected_longitudinal_shower);
// 	 double kolmogorov_dist = 10000.;
// 	 if (h_hgcal_sclusters_longitudinal_shower->Integral()!=0. && h_hgcal_sclusters_expected_longitudinal_shower->Integral()!=0.) 
// 	  kolmogorov_dist = h_hgcal_sclusters_longitudinal_shower->KolmogorovTest(h_hgcal_sclusters_expected_longitudinal_shower,"M");
//          // to be removed with standalone class

	 // standalone class
	 double kolmogorov = hGCALShowerBasedEmIdentification.longitudinalKolmogorov(&(*scl.seed()),false);
	 double kolmogorov_dist = hGCALShowerBasedEmIdentification.longitudinalKolmogorov(&(*scl.seed()),true);
	 // standalone class

	 // plot result
	 if (nevt<10) std::cout << "Kolmogorov test longitudinal profile " << kolmogorov << std::endl;
	 h_hgcal_sclusters_longitudinal_kolmogorov_prob->Fill(kolmogorov);
	 h_hgcal_sclusters_longitudinal_kolmogorov_dist->Fill(kolmogorov_dist);
	 // standalone class

	 // define cut from Kolmogorov dist
	 bool cutkolmogorov = (kolmogorov_dist<0.43);
	 if (nevt<10) std::cout << "shower longitudinal profile Kolmogorov dist cut " << cutkolmogorov << std::endl;
	  
//          // to be removed with standalone class
// 	 // now fit constraining the parameters alpha and beta
//  	 fitparams[0] = 0., fitparams[1] = 0, fitparams[2] = 0;
//  	 double alphamin = std::max(0.,alpha-3.*sigmaalpha), alphamax = std::max(0.,alpha+3*sigmaalpha);
//  	 double invbetamin = std::max(0.,invbeta-3.*sigmainvbeta), invbetamax = std::max(0.,invbeta+3*sigmainvbeta);
//  	 if (nevt<10) std::cout << "alphamin, alphamax " << alphamin << " " << alphamax << std::endl;
//  	 if (nevt<10) std::cout << "invbetamin, invbetamax " << invbetamin << " " << invbetamax << std::endl;
//  	 std::cout << "alphamin, alphamax " << alphamin << " " << alphamax << std::endl;
//  	 std::cout << "invbetamin, invbetamax " << invbetamin << " " << invbetamax << std::endl;
// 	 gammadist->SetParameters(scl.seed()->energy(),alpha,invbeta); 
//          gammadist->SetParLimits(1,alphamin, alphamax);
//          gammadist->SetParLimits(2,invbetamin,invbetamax);
// 	 h_hgcal_sclusters_longitudinal_shower->Fit("gammadist","","",0.,25.);
// 	 double chi2bounded =gammadist->GetChisquare();
// 	 double ndfbounded = gammadist->GetNDF();
//          // to be removed with standalone class
	 
	 // standalone class
	 HGCALFitResults bounded = hGCALShowerBasedEmIdentification.longitudinalFit(&(*scl.seed()),false, true);
	 double chi2bounded =  bounded.chi2();
	 double ndfbounded =  bounded.ndf();
	 // standalone class

	 h_hgcal_sclusters_longitudinal_fit_chi2_bounded->Fill(chi2bounded/ndfbounded);
	 
	 // define cut from bounded fit
	 bool cutfit = ((chi2bounded/ndfbounded)<3.);
	 if (nevt<10) std::cout << "shower fit cut " << cutfit << std::endl;
	 
//          // to be removed with standalone class
// 	 // another fit adding normalization constrained by incoming momentum
//  	 fitparams[0] = 0., fitparams[1] = 0, fitparams[2] = 0;
// 	 gammadist->SetParameters(mcIter->p(),alpha,invbeta); 
//  	 std::cout << "0.9*mcIter->p() " << 0.9*mcIter->p() << std::endl;
// 	 gammadist->SetParLimits(0,0.9*mcIter->p(),1.1*mcIter->p());
//          gammadist->SetParLimits(1,alphamin, alphamax);
//          gammadist->SetParLimits(2,invbetamin,invbetamax);
// 	 h_hgcal_sclusters_longitudinal_shower->Fit("gammadist","","",0.,25.);
// 	 double chi2boundedpnorm =gammadist->GetChisquare();
// 	 double ndfboundedpnorm = gammadist->GetNDF();
// //	 if (gammadist) delete gammadist;
// 	 if (h_hgcal_sclusters_longitudinal_shower) delete h_hgcal_sclusters_longitudinal_shower;
// //	 if (h_hgcal_sclusters_expected_longitudinal_shower) delete h_hgcal_sclusters_expected_longitudinal_shower;
//          // to be removed with standalone class

	 // standalone class
	 HGCALFitResults boundedpnorm = hGCALShowerBasedEmIdentification.longitudinalFit(&(*scl.seed()),true, true);
	 double chi2boundedpnorm =  boundedpnorm.chi2();
	 double ndfboundedpnorm =  boundedpnorm.ndf();
	 // standalone class
	 
	 h_hgcal_sclusters_longitudinal_fit_chi2_pnorm->Fill(chi2boundedpnorm/ndfboundedpnorm);

//          // to be removed with standalone class
// 	 // now define a Kolmogorv test for the transversal profile
// 	 //TF1 *transdist = new TF1("transdist",trans,0.,rmax,4);
// 	 TF1 *transdetectordist = new TF1("transdetectordist",transdetector,0.,rmax,5);
// 	 //TF1 *transdetectordist = new TF1("transdetectordist",transdetector,0.,rmax,4);
//  	 double corefraction = 0.8808;
// 	 double coreradius = 2.589*1.94;
// 	 double tailradius = 1.18*1.94;
// 	 // transverse rt parameters
// //  	 double corefraction = 9.96753e-01;
// // 	 double coreradius = 1.59421e-01; 
// // 	 double tailradius = 5.31567e-01; 
// 	 //transdist->SetParameters(1.,corefraction,coreradius); 
// 	 //transdetectordist->SetParameters(1.,corefraction,coreradius,tailradius); 
// 	 if (transdetectordist->Integral(0.,rmax) == 0.) transdetectordist->SetParameters(1.,0.,0.,1.); ;
// 	 transdetectordist->SetParameters(1.,corefraction,coreradius,tailradius,std::abs(std::cos(pcaShowerDir.theta()))); 
// 	 //transdetectordist->SetParameters(1.,corefraction,coreradius,tailradius); 
// 	 TH1F *h_hgcal_sclusters_expected_transversal_shower = new
// 	  TH1F("h_hgcal_sclusters_expected_transversal_shower","hgcal seed cluster shower transversal profile",nbinr,0.,rmax);
// 	 h_hgcal_sclusters_expected_transversal_shower->Sumw2();
// 	 if (nevt<10) std::cout << "Filling an histogram with the transversal profile function " << std::endl;
// 	 h_hgcal_sclusters_expected_transversal_shower->FillRandom("transdetectordist");
// 	 // normalize prediction
// 	 if (h_hgcal_sclusters_expected_transversal_shower->GetEntries()!=0.) h_hgcal_sclusters_expected_transversal_shower->Scale(sumhits/h_hgcal_sclusters_expected_transversal_shower->GetEntries());
// 	 if (nevt<10) { 
// 	   histfile_->cd();
// 	   h_hgcal_sclusters_expected_transversal_shower->Write();
// 	 }  
// 	 // get kolmogorov result	 
// 	 double kolmogorov_trans = 0.;
// 	 if (h_hgcal_sclusters_transversal_shower->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) 
// 	  kolmogorov_trans = h_hgcal_sclusters_transversal_shower->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower);
// 	 double kolmogorov_dist_trans = 10000.;
//          if (h_hgcal_sclusters_transversal_shower->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) 
// 	  kolmogorov_dist_trans = h_hgcal_sclusters_transversal_shower->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower,"M");
// 	 double kolmogorov_trans_first = 0., kolmogorov_trans_last = 0.;
// 	 double kolmogorov_dist_trans_first = 10000., kolmogorov_dist_trans_last = 10000.;
// 	 if (h_hgcal_sclusters_transversal_shower_first->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) kolmogorov_trans_first = h_hgcal_sclusters_transversal_shower_first->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower);
// 	 if (h_hgcal_sclusters_transversal_shower_first->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) kolmogorov_dist_trans_first = h_hgcal_sclusters_transversal_shower_first->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower,"M");
// 	 if (h_hgcal_sclusters_transversal_shower_last->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) kolmogorov_trans_last = h_hgcal_sclusters_transversal_shower_last->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower);
// 	 if (h_hgcal_sclusters_transversal_shower_last->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) kolmogorov_dist_trans_last = h_hgcal_sclusters_transversal_shower_last->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower,"M");
// 
// //	 if (transdist) delete transdist;
// 	 if (h_hgcal_sclusters_transversal_shower) delete h_hgcal_sclusters_transversal_shower;
// 	 if (h_hgcal_sclusters_transversal_shower_first) delete h_hgcal_sclusters_transversal_shower_first;
// 	 if (h_hgcal_sclusters_transversal_shower_last) delete h_hgcal_sclusters_transversal_shower_last;
// //	 if (h_hgcal_sclusters_expected_transversal_shower) delete h_hgcal_sclusters_expected_transversal_shower;
//          // to be removed with standalone class

	 // standalone class
	 double kolmogorov_trans = hGCALShowerBasedEmIdentification.transverseKolmogorov(&(*scl.seed()),false);
	 double kolmogorov_dist_trans = hGCALShowerBasedEmIdentification.transverseKolmogorov(&(*scl.seed()),true);
	 // standalone class

	 // plot result
	 if (nevt<10) std::cout << "Kolmogorov test transversal profile " << kolmogorov_trans << std::endl;
	 h_hgcal_sclusters_transversal_kolmogorov_prob->Fill(kolmogorov_trans);
	 h_hgcal_sclusters_transversal_kolmogorov_dist->Fill(kolmogorov_dist_trans);
         h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal->Fill(kolmogorov_dist,kolmogorov_dist_trans);	 
// 	 if (nevt<10) std::cout << "Kolmogorov test transversal profile first layers " << kolmogorov_trans_first << std::endl;
// 	 if (nevt<10) std::cout << "Kolmogorov test transversal profile last layers " << kolmogorov_trans_last << std::endl;
// 	 h_hgcal_sclusters_transversal_kolmogorov_prob_first->Fill(kolmogorov_trans_first);
// 	 h_hgcal_sclusters_transversal_kolmogorov_dist_first->Fill(kolmogorov_dist_trans_first);
// 	 h_hgcal_sclusters_transversal_kolmogorov_prob_last->Fill(kolmogorov_trans_last);
// 	 h_hgcal_sclusters_transversal_kolmogorov_dist_last->Fill(kolmogorov_dist_trans_last);

/*
         // ////////////////////////////
         // now a 3D fit/comparison
	 // static Double_t shower3D(Double_t *x, Double_t *par) {
	 //   Double_t longitidinal;
	 //   if (par[1]<0. || par[2]<0.)return 0.;
	 //   else longitudinal = TMath::GammaDist(x[0],par[1],0.,par[2]);
	 //   Double_t core = (x[1]/par[4])*std::exp(-x[1]/par[4]);
	 //   Double_t tail = (x[1]/par[5])*std::exp(-x[1]/par[5]);
	 //   Double_t transverse = par[3]*core + (1.-par[3])*tail;
	 //   return par[0]*longitudinal*transverse*2.*std::acos(-1.);
	 // }
	 // ////////////////////////////
	 
	 TF2 *shower3Ddist = new TF2("shower3Ddist",shower3D,0.,lmax,0.,rmax,11);
	 // longitudinal: alpha, invbeta already defined
         // transverse parametrization as a function of depth
	 double corefraction0 = 0.9733, corefraction1 = 0.02184, corefraction2 = -0.03102; 
	 double coreradius0 = 0.1483, coreradius1 = 0.02721; 
	 double tailradius0 = 0.3422, tailradius1 = 0.02829; 
 	 alphascale=1.20;
	 shower3Ddist->SetParameters(1.,(length/invbeta)*alphascale,invbeta,corefraction0,corefraction1, corefraction2, coreradius0,
	  coreradius1, tailradius0, tailradius1,std::abs(std::cos(pcaShowerDir.theta()))); 
	 TH2F *h_hgcal_sclusters_expected_3Dprofile_shower = new
	  TH2F("h_hgcal_sclusters_expected_3Dprofile_shower","hgcal seed cluster shower 3D profile",nbinz,0.,lmax,nbinr,0.,rmax);
	 h_hgcal_sclusters_expected_3Dprofile_shower->Sumw2();
	 h_hgcal_sclusters_expected_3Dprofile_shower->FillRandom("shower3Ddist");
	 // normalize prediction
	 h_hgcal_sclusters_expected_3Dprofile_shower->Scale(sumhits/h_hgcal_sclusters_expected_3Dprofile_shower->GetEntries());
	 if (nevt<10) h_hgcal_sclusters_expected_3Dprofile_shower->Write();
	 // get kolmogorov result
	 double kolmogorov_3D = h_hgcal_sclusters_3Dprofile_shower->KolmogorovTest(h_hgcal_sclusters_expected_3Dprofile_shower);
	 double kolmogorov_dist_3D = h_hgcal_sclusters_3Dprofile_shower->KolmogorovTest(h_hgcal_sclusters_expected_3Dprofile_shower,"M");
	 // plot result
	 if (nevt<10) std::cout << "Kolmogorov test 3D profile (prob) " << kolmogorov_3D << std::endl;
	 if (nevt<10) std::cout << "Kolmogorov test 3D profile (dist) " << kolmogorov_dist_3D << std::endl;
	 h_hgcal_sclusters_kolmogorov_prob_3D->Fill(kolmogorov_3D);
	 h_hgcal_sclusters_kolmogorov_dist_3D->Fill(kolmogorov_dist_3D);
	 // check projections
	 TH1D *h_hgcal_sclusters_expected_3Dprofile_shower_px = h_hgcal_sclusters_expected_3Dprofile_shower->ProjectionX();
	 TH1D *h_hgcal_sclusters_expected_3Dprofile_shower_py = h_hgcal_sclusters_expected_3Dprofile_shower->ProjectionY();
	 TH1D *h_hgcal_sclusters_3Dprofile_shower_px = h_hgcal_sclusters_3Dprofile_shower->ProjectionX();
	 TH1D *h_hgcal_sclusters_3Dprofile_shower_py = h_hgcal_sclusters_3Dprofile_shower->ProjectionY();
	 double kolmogorov_dist_3D_px = h_hgcal_sclusters_3Dprofile_shower_px->KolmogorovTest(h_hgcal_sclusters_expected_3Dprofile_shower_px,"M");
	 double kolmogorov_dist_3D_py = h_hgcal_sclusters_3Dprofile_shower_py->KolmogorovTest(h_hgcal_sclusters_expected_3Dprofile_shower_py,"M");
	 h_hgcal_sclusters_kolmogorov_dist_3D_px->Fill(kolmogorov_dist_3D_px);
	 h_hgcal_sclusters_kolmogorov_dist_3D_py->Fill(kolmogorov_dist_3D_py);

//	 if (shower3Ddist) delete shower3Ddist;
	 if (h_hgcal_sclusters_3Dprofile_shower) delete h_hgcal_sclusters_3Dprofile_shower;
//	 if (h_hgcal_sclusters_expected_3Dprofile_shower) delete h_hgcal_sclusters_expected_3Dprofile_shower;
	 if (h_hgcal_sclusters_expected_longitudinal_shower) delete h_hgcal_sclusters_expected_longitudinal_shower;
	 if (h_hgcal_sclusters_expected_transversal_shower) delete h_hgcal_sclusters_expected_transversal_shower;
	 if (h_hgcal_sclusters_expected_3Dprofile_shower) delete h_hgcal_sclusters_expected_3Dprofile_shower;
*/
	 ////////////////////////
	 // very loose selection
	 ////////////////////////
	 
	 if (cutseedpos) h_hgcal_scclusters_eoveretrue_em_cutpos->Fill(scle/mcIter->p());
	 if (cutseedpos && cutseedlength) h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength->Fill(scle/mcIter->p());
	 //if (cutseedpos && cutfit) h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength->Fill(scle/mcIter->p());
	 if (cutseedpos && cutseedlength && cutsetaeta) h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta->Fill(scle/mcIter->p());             
	 //if (cutseedpos && cutfit && cutsetaeta) h_hgcal_scclusters_eoveretrue_em_cutpos_cutlength_cutsigmaeta->Fill(scle/mcIter->p());                     
	 
	 // reorganized selection
	 if (cutsetaeta) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta->Fill(scle/mcIter->p());             
	 if (cutsetaeta && cuthadem) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem->Fill(scle/mcIter->p());             
	 if (cutsetaeta && cuthadem && cutkolmogorov) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutkolmogorov->Fill(scle/mcIter->p());             
	 if (cutsetaeta && cuthadem && cutseedpos) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos->Fill(scle/mcIter->p());             
	 if (cutsetaeta && cuthadem && cutseedpos && cutseedlength) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength->Fill(scle/mcIter->p());             

	 //if (cutsetaeta && cuthadem && cutseedpos && cutfit) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength->Fill(scle/mcIter->p());             
	 if (cutsetaeta && cuthadem && cutseedpos && cutseedlength && cutkolmogorov) h_hgcal_scclusters_eoveretrue_em_cutsigmaeta_cuthadem_cutpos_cutlength_cutkolmogorov->Fill(scle/mcIter->p());             
	 
	 // N-1 distributions
 	 if (cutsetaeta) h_hgcal_sclusters_hoverem_cutsigmaeta->Fill(hoverem);
	 if (cutsetaeta && cuthadem) h_hgcal_sclusters_entry_cutsigmaeta_cuthadem->Fill(fabs(firstPos.z()));
	 if (cutsetaeta && cuthadem) h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem->Fill(kolmogorov_dist);
	 if (cutsetaeta && cuthadem && cutseedpos) h_hgcal_sclusters_longitudinal_kolmogorov_dist_cutsigmaeta_cuthadem_cutseedpos->Fill(kolmogorov_dist);
         if (cutsetaeta && cuthadem) h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cutsigmaeta_cuthadem->Fill(kolmogorov_dist,kolmogorov_dist_trans);	 
         if (cuthadem) h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem->Fill(kolmogorov_dist,kolmogorov_dist_trans);	 
         if (cuthadem && cutseedpos) h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_cuthadem_cutseedpos->Fill(kolmogorov_dist,kolmogorov_dist_trans);
	 	 
//	 if (cutsetaeta && cuthadem) h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem->Fill(kolmogorov_dist_best);
//	 if (cutsetaeta && cuthadem && cutseedpos) h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_cutsigmaeta_cuthadem_cutseedpos->Fill(kolmogorov_dist_best);

	 if (cutsetaeta && cuthadem && cutseedpos && cutseedlength) {
	 //if (hGCALShowerBasedEmIdentification.isEm(&(*scl.seed())) && cutseedlength) {
           h_hgcal_sclusters_firsteigenvalue_nm1->Fill(pcaShowerEigenValues.x());
           h_hgcal_sclusters_secondeigenvalue_nm1->Fill(pcaShowerEigenValues.y());
           h_hgcal_sclusters_thirdeigenvalue_nm1->Fill(pcaShowerEigenValues.z());
           h_hgcal_sclusters_firstsigma_nm1->Fill(pcaShowerSigmas.x());
           h_hgcal_sclusters_secondsigma_nm1->Fill(pcaShowerSigmas.y());
           h_hgcal_sclusters_thirdsigma_nm1->Fill(pcaShowerSigmas.z());
	   h_hgcal_sclusters_longitudinal_fit_chi2_nm1->Fill(chi2/ndf);
	   h_hgcal_sclusters_longitudinal_kolmogorov_dist_nm1->Fill(kolmogorov_dist);	  
           h_hgcal_sclusters_kolmogorov_dist_transversalVSlongitudinal_nm1->Fill(kolmogorov_dist,kolmogorov_dist_trans);	 
//	   h_hgcal_sclusters_longitudinal_kolmogorov_dist_best_nm1->Fill(kolmogorov_dist_best);
	 }

	 // now loop on all subclusters
	 	 
	 for (CaloCluster_iterator itcl=scl.clustersBegin(); itcl!=scl.clustersEnd(); itcl++) {
	   //std::cout << "  new sub cluster in HGCAL with energy " << (*itcl)->energy() << 
	   //" , transverse energy " << (*itcl)->energy()*sin(2.*atan(exp(-(*itcl)->position().eta()))) << 
	   //" , position (x,y,z) " << (*itcl)->position() << " and position (eta,phi) : (" 
	   // << (*itcl)->position().eta() << "," << (*itcl)->position().phi() << std::endl;

	   PCAShowerAnalysis pcaShowerAnalysisSub(iEvent,iSetup);

           // check seed cluster is the leading cluster
	   if ((*itcl)->energy() > scl.seed()->energy()) std::cout << "   !!! sub cluster with energy higher than seed cluster !!! " << (*itcl)->energy() << " " << scl.seed()->energy() << std::endl;

	   // get first pos
	   GlobalPoint firstPos;	
	   double zmin = 400.; 

	   for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
	     const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
	     HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
	     //if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
	     if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	       const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
	       GlobalPoint cellPos = geometry_->getPosition(hgcid_);
	       if (fabs(cellPos.z())<zmin) {
		 firstPos = cellPos;
		 zmin = fabs(cellPos.z());
	       }
	     }    
	   }

	   // analyse position and direction with PCA 
	   GlobalPoint pcaShowerPos;
	   GlobalVector pcaShowerDir;
	   pcaShowerAnalysisSub.showerParameters(&(**itcl),pcaShowerPos,pcaShowerDir);
	   if (nevt<10) std::cout << "*** Principal component analysis (subclusters) ****" << std::endl;
	   if (nevt<10) std::cout << "shower average (x,y,z) = " << pcaShowerPos << std::endl;
	   if (nevt<10) std::cout << "shower main axis (x,y,z) = " << pcaShowerDir << std::endl;
	   double length =  (pcaShowerPos - firstPos).mag();
	   if (nevt<10) std::cout << "shower length " << length << " first pos " << firstPos << std::endl;	 

	   double x0 = calohelper_->ecalProperties(true)->radLenIncm();
	   // change here for parametrisation for pi0 extracted from full sim
	   //double predictedLength = meanShower * x0;
	   double predictedLength = 3.6 + 1.383*log((*itcl)->energy());

	   // cut in length, inject here parametrization results
	   double y = (*itcl)->energy()/0.00536;
	   double sigma = predictedLength / (-2.506+1.245*log(y));
	   bool cutsubcllength = fabs(predictedLength-length)<4.*sigma/x0;
	   if (nevt<10 && !cutsubcllength) std::cout << "subcluster rejected from length cut, predicted length " << predictedLength << " sigma " << sigma/x0 << std::endl;

	   // cut in position,cut adapted to pi0s and photons
	   bool cutsubclpos = (fabs(firstPos.z())<322.50);
	   if (nevt<10 && !cutsubcllength) std::cout << "subcluster rejected from position cut, inital z position" << firstPos.z() << std::endl;

	   // final cut on subclusters
	   bool cutsubclall = (cutsubcllength && cutsubclpos);

	   // apply cleaning cut here, do not cut the seed cluster
	   if (!cutsubclall && itcl!=scl.clustersBegin()) continue;
 	   if (!cutsubclall && nevt<10) std::cout << "subcluster rejected from position and length cuts " << std::endl;          
	   
	   // design new supercluster region
	   //bool superclusterphiroad = true;
           // first tighten the road in eta between 0.003 and 0.05
	   //double detamax = 0.003;
	   //if (std::abs(scl.seed()->eta())>2.) detamax = 0.04;
	   //double dphimax = 0.15;
	   //double detasub = (*itcl)->eta()-scl.seed()->eta();
	   //double dphisub = (*itcl)->phi()-scl.seed()->phi();
	   //if (std::abs(dphisub)>CLHEP::pi) dphisub = dphisub < 0? (CLHEP::twopi) + dphisub : dphisub - CLHEP::twopi;
	   
	   // remove subcluster not in new supercluster region
	   //superclusterphiroad = std::abs(detasub)<detamax && std::abs(dphisub)<dphimax;
	    
	   //if (!superclusterphiroad) continue;

	   newenergy += (*itcl)->energy();
	   
	   //GlobalVector dir, vect; GlobalPoint firstPos;
	   h_hgcal_clusters_energy_em->Fill((*itcl)->energy());
	   h_hgcal_clusters_position_em->Fill((*itcl)->eta(),(*itcl)->phi());
	   h_hgcal_clusters_multiplicity_em->Fill(float((*itcl)->hitsAndFractions().size()));
	   h_hgcal_clusters_multiplicityVSeta_em->Fill(float((*itcl)->hitsAndFractions().size()),mcIter->eta());
	   double deta = (*itcl)->eta() - scl.seed()->eta();
	   double dphi = (*itcl)->phi() - scl.seed()->phi();
	   if (std::abs(dphi)>CLHEP::pi)
            dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	   //double dx = (*itcl)->position().x() - scl.seed()->position().x(); 
	   //double dy = (*itcl)->position().y() - scl.seed()->position().y(); 
	   if (itcl!=scl.clustersBegin()) h_hgcal_scclusters_detadphisubclusters_em->Fill(deta,dphi);
	   if (itcl!=scl.clustersBegin()) h_hgcal_scclusters_detadphisubclusters_zoom_em->Fill(deta,dphi);
	   if (itcl!=scl.clustersBegin() && fabs(mcIter->eta())>2.) h_hgcal_scclusters_detadphisubclusters_zoom_etagt2_em->Fill(deta,dphi);
	   if (itcl!=scl.clustersBegin() && fabs(mcIter->eta())<2.) h_hgcal_scclusters_detadphisubclusters_zoom_etalt2_em->Fill(deta,dphi);
	   if (itcl!=scl.clustersBegin()) h_hgcal_scclusters_detadphisubclusters_weighted_em->Fill(deta,dphi,(*itcl)->energy());
	   if (itcl!=scl.clustersBegin()) h_hgcal_scclusters_detadphisubclusters_zoom_weighted_em->Fill(deta,dphi,(*itcl)->energy());
           double newenergy000cl = 0., newenergy04cl = 0., newenergy1cl = 0., newenergy2cl = 0., newenergy4cl = 0., newenergy10cl = 0., newenergy20cl = 0.;       
	   // eta correction
       	   const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
           const double _coef_b = -107.229,
  	   //double corr = 80.0837 / (1.0 + exp(-0.0472817 - 107.229*cosh((*itcl)->eta()))); 
           // new calib as of SLHC21
	   clus_eta = (*itcl)->eta();
           double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
	   double mip = 0.0000551;
	   double weight[30] =
	    {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
	   for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
             const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
             HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
	     //if (theHit->energy()>0) std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
	     if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	     const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
             h_hgcal_clusters_rechitenergy_em->Fill(theHit->energy());
             h_hgcal_clusters_rechitenergy_12000_em->Fill(theHit->energy());
	     //GlobalPoint cellPos = geometry_->getGeometry(id_)->getPosition();	     
	     GlobalPoint cellPos = geometry_->getPosition(hgcid_);
            if (nevt<10) {
	       h_hgcal_shower_sc_em[nevt]->Fill(fabs(cellPos.x()),fabs(cellPos.y()),fabs(cellPos.z()),theHit->energy()*1000.);
              h_hgcal_shower_sc[nevt]->Fill(fabs(cellPos.x()),fabs(cellPos.y()),fabs(cellPos.z()),theHit->energy()*1000.);
            }
	     //std::cout << "    cell position " << cellPos << " energy " << theSeedHit->energy() << std::endl;
	     // recompute calibrated SC energy
             const int layer = hgcid_.layer();
	     double scale = mip*corr;
	     newenergy000cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>0.4*0.0000551) newenergy04cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>1.*0.0000551) newenergy1cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>2.*0.0000551) newenergy2cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>4.*0.0000551) newenergy4cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>10.*0.0000551) newenergy10cl += theHit->energy()*weight[layer-1]/scale;
	     if (theHit->energy()>20.*0.0000551) newenergy20cl += theHit->energy()*weight[layer-1]/scale;
	   } 
	 } 
	 newenergy000cl = std::max(0., newenergy000cl-_coef_b/_coef_a);
	 newenergy04cl = std::max(0., newenergy04cl-_coef_b/_coef_a);
	 newenergy1cl = std::max(0., newenergy1cl-_coef_b/_coef_a);
	 newenergy2cl = std::max(0., newenergy2cl-_coef_b/_coef_a);
	 newenergy4cl = std::max(0., newenergy4cl-_coef_b/_coef_a);
	 newenergy10cl = std::max(0., newenergy10cl-_coef_b/_coef_a);
	 newenergy20cl = std::max(0., newenergy20cl-_coef_b/_coef_a);
	 if (itcl==scl.clustersBegin()) std::cout << "initial seed cluster energy " << (*itcl)->energy() << std::endl; 
	 if (itcl==scl.clustersBegin()) std::cout << "recomputed seed cluster energy " << newenergy000cl << std::endl; 
	 if (nevt<10) {
	   std::cout << "initial cluster energy " << (*itcl)->energy() << std::endl;  
	   std::cout << "recomputed cluster energy " << newenergy000cl << std::endl; 
	 }
	 newenergy000 += newenergy000cl;
	 newenergy04 += newenergy04cl;
	 newenergy1 += newenergy1cl;
	 newenergy2 += newenergy2cl;
	 newenergy4 += newenergy4cl;
	 newenergy10 += newenergy10cl;
	 newenergy20 += newenergy20cl;

       } // end loop on subclusters	  
         
       std::cout << "SC energy " << scl.energy() << std::endl;
       std::cout << "recomputed calibrated energy " << newenergy000 << std::endl;
       h_hgcal_scclusters_eoveretrue_cut00_em->Fill(newenergy000/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut04_em->Fill(newenergy04/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut1_em->Fill(newenergy1/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut2_em->Fill(newenergy2/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut4_em->Fill(newenergy4/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut10_em->Fill(newenergy10/mcIter->p());
       h_hgcal_scclusters_eoveretrue_cut20_em->Fill(newenergy20/mcIter->p());

       // plot eovertrue with selections for noise cut at 2 mip
       h_hgcal_scclusters_eoveretrue_noisecut_em->Fill(newenergy1/mcIter->p());
       if (cutseedpos) h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos->Fill(newenergy1/mcIter->p());
       if (cutseedpos && cutseedlength) h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength->Fill(newenergy1/mcIter->p());
       if (cutseedpos && cutseedlength && cutsetaeta) h_hgcal_scclusters_eoveretrue_noisecut_em_cutpos_cutlength_cutsigmaeta->Fill(newenergy1/mcIter->p());
       
       // plot eovertrue with selections for noise cut at 4 mip
       h_hgcal_scclusters_eoveretrue_noisecut4_em->Fill(newenergy4/mcIter->p());
       if (cutseedpos) h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos->Fill(newenergy4/mcIter->p());
       if (cutseedpos && cutseedlength) h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength->Fill(newenergy4/mcIter->p());
       if (cutseedpos && cutseedlength && cutsetaeta) h_hgcal_scclusters_eoveretrue_noisecut4_em_cutpos_cutlength_cutsigmaeta->Fill(newenergy4/mcIter->p());
       

     } // end if in HGCEE

     // here we are back to treatment for each supercluster best match to the generated electron
     h_hgcal_sclusters_newmultiplicity_em->Fill(float(newscmultiplicity));
     h_hgcal_sclusters_newmultiplicityVSeta_em->Fill(float(newscmultiplicity),mcIter->eta());
     h_hgcal_sclusters_newseedenergy_em->Fill(newseedenergy);

     } // end loop on SCs  
  
     } // gsf electron found

   } // mc particle found

   }

  } // loop over mc particle

  h_mcNum->Fill(mcNum);

}


