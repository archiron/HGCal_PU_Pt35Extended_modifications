// -*- C++ -*-
//
// Package:    RecoEgamma/Examples
// Class:      HGCALGsfElectronAnalyzer
//
/**\class HGCALGsfElectronAnalyzer RecoEgamma/Examples/src/HGCALGsfElectronAnalyzer.cc

 Description: GsfElectrons analyzer using MC truth

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ursula Berthon
//         Created:  Mon Mar 27 13:22:06 CEST 2006
// $Id: HGCALGsfElectronAnalyzer.cc,v 1.51 2011/03/04 14:43:15 chamont Exp $
//
//

// user include files
#include "RecoEgamma/Examples/plugins/HGCALGsfElectronAnalyzer.h"
#include "RecoEgamma/Examples/interface/PCAShowerAnalysis.h"
#include "RecoEgamma/Examples/interface/HGCALShowerBasedEmIdentification.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronClassification.h"

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
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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

HGCALGsfElectronAnalyzer::HGCALGsfElectronAnalyzer(const edm::ParameterSet& conf)
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

void HGCALGsfElectronAnalyzer::readParameters(const edm::ParameterSet& fastCalo) {

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

void HGCALGsfElectronAnalyzer::beginJob(){

  histfile_->cd();

  // mc truth
  h_mcNum              = new TH1F( "h_mcNum",              "# mc particles",    nbinfhits,0.,fhitsmax );
  h_mcNum->Sumw2();
  h_eleNum             = new TH1F( "h_mcNum_ele",             "# mc electrons",             nbinfhits,0.,fhitsmax);
  h_eleNum->Sumw2();
  h_neleinEE            = new TH1F( "h_num_ele_in_EE",             "# mc electrons",             nbinfhits,0.,fhitsmax);
  h_neleinEE->Sumw2();
  h_gamNum             = new TH1F( "h_mcNum_gam",             "# mc gammas",             nbinfhits,0.,fhitsmax);
  h_gamNum->Sumw2();

  // rec event
  histNum_= new TH1F("h_recEleNum","# rec electrons",20, 0.,20.);

  // mc
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

  // all electrons
  h_ele_EoverP_all       = new TH1F( "h_ele_EoverP_all",       "ele E/P_{vertex}, all reco electrons",  nbineop,0.,eopmax);
  h_ele_EoverP_all->Sumw2();
  h_ele_EoverP_all_barrel       = new TH1F( "h_ele_EoverP_all_barrel",       "ele E/P_{vertex}, all reco electrons, barrel",  nbineop,0.,eopmax);
  h_ele_EoverP_all_barrel->Sumw2();
  h_ele_EoverP_all_endcaps       = new TH1F( "h_ele_EoverP_all_endcaps",       "ele E/P_{vertex}, all reco electrons, endcaps",  nbineop,0.,eopmax);
  h_ele_EoverP_all_endcaps->Sumw2();
  h_ele_EseedOP_all            = new TH1F( "h_ele_EseedOP_all",            "ele E_{seed}/P_{vertex}, all reco electrons",        nbineop,0.,eopmax);
  h_ele_EseedOP_all->Sumw2();
  h_ele_EseedOP_all_barrel            = new TH1F( "h_ele_EseedOP_all_barrel",            "ele E_{seed}/P_{vertex}, all reco electrons, barrel",        nbineop,0.,eopmax);
  h_ele_EseedOP_all_barrel->Sumw2();
  h_ele_EseedOP_all_endcaps            = new TH1F( "h_ele_EseedOP_all_endcaps",            "ele E_{seed}/P_{vertex}, all reco electrons, endcaps",        nbineop,0.,eopmax);
  h_ele_EseedOP_all_endcaps->Sumw2();
  h_ele_EoPout_all         = new TH1F( "h_ele_EoPout_all",         "ele E_{seed}/P_{out}, all reco electrons",           nbineop,0.,eopmax);
  h_ele_EoPout_all->Sumw2();
  h_ele_EoPout_all_barrel         = new TH1F( "h_ele_EoPout_all_barrel",         "ele E_{seed}/P_{out}, all reco electrons barrel",           nbineop,0.,eopmax);
  h_ele_EoPout_all_barrel->Sumw2();
  h_ele_EoPout_all_endcaps         = new TH1F( "h_ele_EoPout_all_endcaps",         "ele E_{seed}/P_{out}, all reco electrons endcaps",           nbineop,0.,eopmax);
  h_ele_EoPout_all_endcaps->Sumw2();
  h_ele_EeleOPout_all         = new TH1F( "h_ele_EeleOPout_all",         "ele E_{ele}/P_{out}, all reco electrons",           nbineop,0.,eopmax);
  h_ele_EeleOPout_all->Sumw2();
  h_ele_EeleOPout_all_barrel         = new TH1F( "h_ele_EeleOPout_all_barrel",         "ele E_{ele}/P_{out}, all reco electrons barrel",           nbineop,0.,eopmax);
  h_ele_EeleOPout_all_barrel->Sumw2();
  h_ele_EeleOPout_all_endcaps         = new TH1F( "h_ele_EeleOPout_all_endcaps",         "ele E_{ele}/P_{out}, all reco electrons endcaps",           nbineop,0.,eopmax);
  h_ele_EeleOPout_all_endcaps->Sumw2();
  h_ele_dEtaSc_propVtx_all = new TH1F( "h_ele_dEtaSc_propVtx_all", "ele #eta_{sc} - #eta_{tr}, prop from vertex, all reco electrons",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_all->Sumw2();
  h_ele_dEtaSc_propVtx_all_barrel = new TH1F( "h_ele_dEtaSc_propVtx_all_barrel", "ele #eta_{sc} - #eta_{tr}, prop from vertex, all reco electrons barrel",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_all_barrel->Sumw2();
  h_ele_dEtaSc_propVtx_all_endcaps = new TH1F( "h_ele_dEtaSc_propVtx_all_endcaps", "ele #eta_{sc} - #eta_{tr}, prop from vertex, all reco electrons endcaps",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_all_endcaps->Sumw2();
  h_ele_dPhiSc_propVtx_all = new TH1F( "h_ele_dPhiSc_propVtx_all", "ele #phi_{sc} - #phi_{tr}, prop from vertex, all reco electrons",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_all->Sumw2();
  h_ele_dPhiSc_propVtx_all_barrel = new TH1F( "h_ele_dPhiSc_propVtx_all_barrel", "ele #phi_{sc} - #phi_{tr}, prop from vertex, all reco electrons barrel",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_all_barrel->Sumw2();
  h_ele_dPhiSc_propVtx_all_endcaps = new TH1F( "h_ele_dPhiSc_propVtx_all_endcaps", "ele #phi_{sc} - #phi_{tr}, prop from vertex, all reco electrons endcaps",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_all_endcaps->Sumw2();
  h_ele_dEtaCl_propOut_all = new TH1F( "h_ele_dEtaCl_propOut_all", "ele #eta_{cl} - #eta_{tr}, prop from outermost, all reco electrons",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_all->Sumw2();
  h_ele_dEtaCl_propOut_all_barrel = new TH1F( "h_ele_dEtaCl_propOut_all_barrel", "ele #eta_{cl} - #eta_{tr}, prop from outermost, all reco electrons barrel",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_all_barrel->Sumw2();
  h_ele_dEtaCl_propOut_all_endcaps = new TH1F( "h_ele_dEtaCl_propOut_all_endcaps", "ele #eta_{cl} - #eta_{tr}, prop from outermost, all reco electrons endcaps",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_all_endcaps->Sumw2();
  h_ele_dPhiCl_propOut_all = new TH1F( "h_ele_dPhiCl_propOut_all", "ele #phi_{cl} - #phi_{tr}, prop from outermost, all reco electrons",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_all->Sumw2();
  h_ele_dPhiCl_propOut_all_barrel = new TH1F( "h_ele_dPhiCl_propOut_all_barrel", "ele #phi_{cl} - #phi_{tr}, prop from outermost, all reco electrons barrel",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_all_endcaps = new TH1F( "h_ele_dPhiCl_propOut_all_endcaps", "ele #phi_{cl} - #phi_{tr}, prop from outermost, all reco electrons endcaps",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_all_barrel->Sumw2();
  h_ele_dPhiCl_propOut_all_endcaps->Sumw2();
  h_ele_HoE_all = new TH1F("h_ele_HoE_all", "ele hadronic energy / em energy, all reco electrons", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_all->Sumw2();
  h_ele_HoE_all_barrel = new TH1F("h_ele_HoE_all_barrel", "ele hadronic energy / em energy, all reco electrons barrel", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_all_barrel->Sumw2();
  h_ele_HoE_all_endcaps = new TH1F("h_ele_HoE_all_endcaps", "ele hadronic energy / em energy, all reco electrons endcaps", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_all_endcaps->Sumw2();
  h_ele_vertexPt_all       = new TH1F( "h_ele_vertexPt_all",       "ele p_{T}, all reco electrons",  nbinpteff,5.,ptmax);
  h_ele_vertexPt_all->Sumw2();
  h_ele_Et_all       = new TH1F( "h_ele_Et_all",       "ele SC E_{T}, all reco electrons",  nbinpteff,5.,ptmax);
  h_ele_Et_all->Sumw2();
  h_ele_vertexEta_all      = new TH1F( "h_ele_vertexEta_all",      "ele eta, all reco electrons",    nbineta,etamin,etamax);
  h_ele_vertexEta_all->Sumw2();
  h_ele_TIP_all       = new TH1F( "h_ele_TIP_all",       "ele vertex transverse radius, all reco electrons",  100,0.,0.2);
  h_ele_TIP_all->Sumw2();
  h_ele_TIP_all_barrel       = new TH1F( "h_ele_TIP_all_barrel",       "ele vertex transverse radius, all reco electrons barrel",  100,0.,0.2);
  h_ele_TIP_all_barrel->Sumw2();
  h_ele_TIP_all_endcaps       = new TH1F( "h_ele_TIP_all_endcaps",       "ele vertex transverse radius, all reco electrons endcaps",  100,0.,0.2);
  h_ele_TIP_all_endcaps->Sumw2();
  h_ele_mee_all      = new TH1F( "h_ele_mee_all", "ele pairs invariant mass, all reco electrons", nbinmee, meemin, meemax );
  h_ele_mee_all->Sumw2();
  h_ele_mee_seed_all      = new TH1F( "h_ele_mee_seed_all", "ele pairs invariant mass from seed, all reco electrons", nbinmee, meemin, meemax );
  h_ele_mee_seed_all->Sumw2();
  h_ele_mee_sc_all      = new TH1F( "h_ele_mee_sc_all", "ele pairs invariant mass from sc, all reco electrons", nbinmee, meemin, meemax );
  h_ele_mee_sc_all->Sumw2();
  h_ele_mee_newsc_all      = new TH1F( "h_ele_mee_newsc_all", "ele pairs invariant mass from sc, all reco electrons", nbinmee, meemin, meemax );
  h_ele_mee_newsc_all->Sumw2();
  h_ele_mee_best_all      = new TH1F( "h_ele_mee_best_all", "ele pairs invariant mass from sc, all reco electrons", nbinmee, meemin, meemax );
  h_ele_mee_best_all->Sumw2();
  h_ele_mee_seed_os      = new TH1F( "h_ele_mee_seed_os", "ele pairs invariant mass from seed, opp. sign", nbinmee, meemin, meemax );
  h_ele_mee_seed_os->Sumw2();
  h_ele_mee_best_os      = new TH1F( "h_ele_mee_best_os", "ele pairs invariant mass from seed, opp. sign", nbinmee, meemin, meemax );
  h_ele_mee_best_os->Sumw2();
  h_ele_mee_sc_os      = new TH1F( "h_ele_mee_sc_os", "ele pairs invariant mass from sc, opp. sign", nbinmee, meemin, meemax );
  h_ele_mee_sc_os->Sumw2();
  h_ele_mee_newsc_os      = new TH1F( "h_ele_mee_newsc_os", "ele pairs invariant mass from sc, opp. sign", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os->Sumw2();
  h_ele_mee_os      = new TH1F( "h_ele_mee_os", "ele pairs invariant mass, opp. sign", nbinmee, meemin, meemax );
  h_ele_mee_os->Sumw2();
  h_ele_mee_os_ebeb      = new TH1F( "h_ele_mee_os_ebeb", "ele pairs invariant mass, opp. sign, EB-EB", nbinmee, meemin, meemax );
  h_ele_mee_os_ebeb->Sumw2();
  h_ele_mee_best_os_ebeb      = new TH1F( "h_ele_mee_best_os_ebeb", "ele pairs invariant mass, opp. sign, EB-EB", nbinmee, meemin, meemax );
  h_ele_mee_best_os_ebeb->Sumw2();
  h_ele_mee_os_ebee      = new TH1F( "h_ele_mee_os_ebee", "ele pairs invariant mass, opp. sign, EB-EE", nbinmee, meemin, meemax );
  h_ele_mee_os_ebee->Sumw2();
  h_ele_mee_sc_os_ebee      = new TH1F( "h_ele_mee_sc_os_ebee", "ele pairs invariant mass, opp. sign, EB-EE", nbinmee, meemin, meemax );
  h_ele_mee_sc_os_ebee->Sumw2();
  h_ele_mee_newsc_os_ebee      = new TH1F( "h_ele_mee_newsc_os_ebee", "ele pairs invariant mass, opp. sign, EB-EE", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_ebee->Sumw2();
  h_ele_mee_seed_os_ebee      = new TH1F( "h_ele_mee_seed_os_ebee", "ele pairs invariant mass, opp. sign, EB-EE", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_ebee->Sumw2();
  h_ele_mee_best_os_ebee      = new TH1F( "h_ele_mee_best_os_ebee", "ele pairs invariant mass, opp. sign, EB-EE", nbinmee, meemin, meemax );
  h_ele_mee_best_os_ebee->Sumw2();
  h_ele_mee_os_eeee      = new TH1F( "h_ele_mee_os_eeee", "ele pairs invariant mass, opp. sign, EE-EE", nbinmee, meemin, meemax );
  h_ele_mee_os_eeee->Sumw2();
  h_ele_mee_sc_os_eeee      = new TH1F( "h_ele_mee_sc_os_eeee", "ele pairs invariant mass, opp. sign, EE-EE", nbinmee, meemin, meemax );
  h_ele_mee_sc_os_eeee->Sumw2();
  h_ele_mee_newsc_os_eeee      = new TH1F( "h_ele_mee_newsc_os_eeee", "ele pairs invariant mass, opp. sign, EE-EE", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_eeee->Sumw2();
  h_ele_mee_seed_os_eeee      = new TH1F( "h_ele_mee_seed_os_eeee", "ele pairs invariant mass, opp. sign, EE-EE", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_eeee->Sumw2();
  h_ele_mee_best_os_eeee      = new TH1F( "h_ele_mee_best_os_eeee", "ele pairs invariant mass, opp. sign, EE-EE", nbinmee, meemin, meemax );
  h_ele_mee_best_os_eeee->Sumw2();
  h_ele_mee_os_gg      = new TH1F( "h_ele_mee_os_gg", "ele pairs invariant mass, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_os_gg->Sumw2();
  h_ele_mee_os_gb      = new TH1F( "h_ele_mee_os_gb", "ele pairs invariant mass, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_os_gb->Sumw2();
  h_ele_mee_os_bb      = new TH1F( "h_ele_mee_os_bb", "ele pairs invariant mass, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_os_bb->Sumw2();
  h_ele_mee_seed_os_gg      = new TH1F( "h_ele_mee_seed_os_gg", "ele pairs invariant mass from seed, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_gg->Sumw2();
  h_ele_mee_seed_os_gg_ebeb      = new TH1F( "h_ele_mee_seed_os_gg_ebeb", "ele pairs invariant mass from seed, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_gg_ebeb->Sumw2();
  h_ele_mee_seed_os_gg_ebee      = new TH1F( "h_ele_mee_seed_os_gg_ebee", "ele pairs invariant mass from seed, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_gg_ebee->Sumw2();
  h_ele_mee_seed_os_gg_eeee      = new TH1F( "h_ele_mee_seed_os_gg_eeee", "ele pairs invariant mass from seed, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_gg_eeee->Sumw2();
  h_ele_mee_seed_os_gb      = new TH1F( "h_ele_mee_seed_os_gb", "ele pairs invariant mass from seed, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_gb->Sumw2();
  h_ele_mee_seed_os_bb      = new TH1F( "h_ele_mee_seed_os_bb", "ele pairs invariant mass from seed, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_seed_os_bb->Sumw2();
  h_ele_mee_sc_os_gg      = new TH1F( "h_ele_mee_sc_os_gg", "ele pairs invariant mass from sc, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_sc_os_gg->Sumw2();
  h_ele_mee_sc_os_gb      = new TH1F( "h_ele_mee_sc_os_gb", "ele pairs invariant mass from sc, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_sc_os_gb->Sumw2();
  h_ele_mee_sc_os_bb      = new TH1F( "h_ele_mee_sc_os_bb", "ele pairs invariant mass from sc, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_sc_os_bb->Sumw2();
  h_ele_mee_newsc_os_gg      = new TH1F( "h_ele_mee_newsc_os_gg", "ele pairs invariant mass from sc, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gg->Sumw2();
  h_ele_mee_newsc_os_gg_ebeb      = new TH1F( "h_ele_mee_newsc_os_gg_ebeb", "ele pairs invariant mass from sc, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gg_ebeb->Sumw2();
  h_ele_mee_newsc_os_gg_ebee      = new TH1F( "h_ele_mee_newsc_os_gg_ebee", "ele pairs invariant mass from sc, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gg_ebee->Sumw2();
  h_ele_mee_newsc_os_gg_eeee      = new TH1F( "h_ele_mee_newsc_os_gg_eeee", "ele pairs invariant mass from sc, opp. sign, good-good", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gg_eeee->Sumw2();
  h_ele_mee_newsc_os_gb      = new TH1F( "h_ele_mee_newsc_os_gb", "ele pairs invariant mass from sc, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gb->Sumw2();
  h_ele_mee_newsc_os_gb_ebee      = new TH1F( "h_ele_mee_newsc_os_gb_ebee", "ele pairs invariant mass from sc, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gb_ebee->Sumw2();
  h_ele_mee_newsc_os_gb_eeee      = new TH1F( "h_ele_mee_newsc_os_gb_eeee", "ele pairs invariant mass from sc, opp. sign, good-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_gb_eeee->Sumw2();
  h_ele_mee_newsc_os_bb      = new TH1F( "h_ele_mee_newsc_os_bb", "ele pairs invariant mass from sc, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_bb->Sumw2();
  h_ele_mee_newsc_os_bb_ebeb      = new TH1F( "h_ele_mee_newsc_os_bb_ebeb", "ele pairs invariant mass from sc, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_bb_ebeb->Sumw2();
  h_ele_mee_newsc_os_bb_ebee      = new TH1F( "h_ele_mee_newsc_os_bb_ebee", "ele pairs invariant mass from sc, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_bb_ebee->Sumw2();
  h_ele_mee_newsc_os_bb_eeee      = new TH1F( "h_ele_mee_newsc_os_bb_eeee", "ele pairs invariant mass from sc, opp. sign, bad-bad", nbinmee, meemin, meemax );
  h_ele_mee_newsc_os_bb_eeee->Sumw2();

  // duplicates
  h_ele_E2mnE1vsMee_all = new TH2F("h_ele_E2mnE1vsMee_all", "E2 - E1 vs ele pairs invariant mass, all electrons", nbinmee, meemin, meemax, 100, -50., 50. );
  h_ele_E2mnE1vsMee_egeg_all = new TH2F("h_ele_E2mnE1vsMee_egeg_all", "E2 - E1 vs ele pairs invariant mass, ecal driven pairs, all electrons", nbinmee, meemin, meemax, 100, -50., 50. );

  // charge ID
  h_ele_ChargeMnChargeTrue   = new TH1F( "h_ele_ChargeMnChargeTrue",   "ele charge - gen charge ",5,-1.,4.);
  h_ele_ChargeMnChargeTrue->Sumw2();
  h_ele_simEta_matched_qmisid             = new TH1F( "h_ele_eta_matched_qmisid","charge misid vs gen eta", nbineta,etamin,etamax);
  h_ele_simEta_matched_qmisid->Sumw2();
  h_ele_simAbsEta_matched_qmisid             = new TH1F( "h_ele_abseta_matched_qmisid", "charge misid vs gen |eta|", nbineta/2,0.,etamax);
  h_ele_simAbsEta_matched_qmisid->Sumw2();
  h_ele_simPt_matched_qmisid               = new TH1F( "h_ele_Pt_matched_qmisid", "charge misid vs gen transverse momentum", nbinpteff,5.,ptmax);
  h_ele_simPt_matched_qmisid->Sumw2();
  h_ele_simPhi_matched_qmisid               = new TH1F( "h_ele_phi_matched_qmisid","charge misid vs gen phi", nbinphi,phimin,phimax);
  h_ele_simPhi_matched_qmisid->Sumw2();
  h_ele_simZ_matched_qmisid      = new TH1F( "h_ele_z_matched_qmisid","charge misid vs gen z",nbinxyz, -25, 25 );
  h_ele_simZ_matched_qmisid->Sumw2();

  // matched electrons
  h_ele_charge         = new TH1F( "h_ele_charge",         "ele charge",             5,-2.,2.);
  h_ele_charge->Sumw2();
  h_ele_chargeVsEta    = new TH2F( "h_ele_chargeVsEta",         "ele charge vs eta", nbineta2D,etamin,etamax,5,-2.,2.);
  h_ele_chargeVsPhi    = new TH2F( "h_ele_chargeVsPhi",         "ele charge vs phi", nbinphi2D,phimin,phimax,5,-2.,2.);
  h_ele_chargeVsPt    = new TH2F( "h_ele_chargeVsPt",         "ele charge vs pt", nbinpt,0.,100.,5,-2.,2.);
  h_ele_vertexP        = new TH1F( "h_ele_vertexP",        "ele momentum",       nbinp,0.,pmax);
  h_ele_vertexP->Sumw2();
  h_ele_vertexPt       = new TH1F( "h_ele_vertexPt",       "ele transverse momentum",  nbinpt,0.,ptmax);
  h_ele_vertexPt->Sumw2();
  h_ele_Et       = new TH1F( "h_ele_Et",       "ele transverse energy",  nbinpt,0.,ptmax);
  h_ele_Et->Sumw2();
  h_ele_vertexPtVsEta   = new TH2F( "h_ele_vertexPtVsEta",       "ele transverse momentum vs eta",nbineta2D,etamin,etamax,nbinpt2D,0.,ptmax);
  h_ele_vertexPtVsPhi   = new TH2F( "h_ele_vertexPtVsPhi",       "ele transverse momentum vs phi",nbinphi2D,phimin,phimax,nbinpt2D,0.,ptmax);
  h_ele_simPt_matched       = new TH1F( "h_ele_simPt_matched",       "Efficiency vs gen transverse momentum",  nbinpteff,5.,ptmax);
  h_ele_vertexEta      = new TH1F( "h_ele_vertexEta",      "ele momentum eta",    nbineta,etamin,etamax);
  h_ele_vertexEta->Sumw2();
  h_ele_vertexEtaVsPhi  = new TH2F( "h_ele_vertexEtaVsPhi",      "ele momentum eta vs phi",nbineta2D,etamin,etamax,nbinphi2D,phimin,phimax );
  h_ele_simAbsEta_matched      = new TH1F( "h_ele_simAbsEta_matched",      "Efficiency vs gen |eta|",    nbineta/2,0.,etamax);
  h_ele_simAbsEta_matched->Sumw2();
  h_ele_simEta_matched      = new TH1F( "h_ele_simEta_matched",      "Efficiency vs gen eta",    nbineta,etamin,etamax);
  h_ele_simEta_matched->Sumw2();
  h_ele_simPtEta_matched           = new TH2F( "h_ele_simPtEta_matched",   "Efficiency vs pt #eta",  nbineta2D,etamin,etamax, nbinpt2D,5.,ptmax );
  h_ele_simPtEta_matched->Sumw2();
  h_ele_simPhi_matched               = new TH1F( "h_ele_simPhi_matched",               "Efficiency vs gen phi",        nbinphi,phimin,phimax);
  h_ele_simPhi_matched->Sumw2();
  h_ele_vertexPhi      = new TH1F( "h_ele_vertexPhi",      "ele  momentum #phi",    nbinphi,phimin,phimax);
  h_ele_vertexPhi->Sumw2();
  h_ele_vertexX      = new TH1F( "h_ele_vertexX",      "ele vertex x",    nbinxyz,-0.1,0.1 );
  h_ele_vertexX->Sumw2();
  h_ele_vertexY      = new TH1F( "h_ele_vertexY",      "ele vertex y",    nbinxyz,-0.1,0.1 );
  h_ele_vertexY->Sumw2();
  h_ele_vertexZ      = new TH1F( "h_ele_vertexZ",      "ele vertex z",    nbinxyz,-25, 25 );
  h_ele_vertexZ->Sumw2();
  h_ele_simZ_matched      = new TH1F( "h_ele_simZ_matched",      "Efficiency vs gen vertex z",    nbinxyz,-25,25);
  h_ele_simZ_matched->Sumw2();
  h_ele_vertexTIP      = new TH1F( "h_ele_vertexTIP",      "ele transverse impact parameter (wrt gen vtx)",    90,0.,0.15);
  h_ele_vertexTIP->Sumw2();
  h_ele_vertexTIPVsEta      = new TH2F( "h_ele_vertexTIPVsEta",      "ele transverse impact parameter (wrt gen vtx) vs eta", nbineta2D,etamin,etamax,45,0.,0.15);
  h_ele_vertexTIPVsPhi      = new TH2F( "h_ele_vertexTIPVsPhi",      "ele transverse impact parameter (wrt gen vtx) vs phi", nbinphi2D,phimin,phimax,45,0.,0.15);
  h_ele_vertexTIPVsPt      = new TH2F( "h_ele_vertexTIPVsPt",      "ele transverse impact parameter (wrt gen vtx) vs transverse momentum", nbinpt2D,0.,ptmax,45,0.,0.15);
  h_ele_PoPtrue        = new TH1F( "h_ele_PoPtrue",        "ele momentum / gen momentum", nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue->Sumw2();
  h_ele_PtoPttrue        = new TH1F( "h_ele_PtoPttrue",        "ele transverse momentum / gen transverse momentum", nbinpoptrue,poptruemin,poptruemax);
  h_ele_PtoPttrue->Sumw2();
  h_ele_PoPtrueVsEta   = new TH2F( "h_ele_PoPtrueVsEta",        "ele momentum / gen momentum vs eta", nbineta2D,etamin,etamax,50,poptruemin,poptruemax);
  h_ele_PoPtrueVsPhi   = new TH2F( "h_ele_PoPtrueVsPhi",        "ele momentum / gen momentum vs phi", nbinphi2D,phimin,phimax,50,poptruemin,poptruemax);
  h_ele_PoPtrueVsPt   = new TH2F( "h_ele_PoPtrueVsPt",        "ele momentum / gen momentum vs eta", nbinpt2D,0.,ptmax,50,poptruemin,poptruemax);
  h_ele_PoPtrue_barrel         = new TH1F( "h_ele_PoPtrue_barrel",        "ele momentum / gen momentum, barrel",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_barrel->Sumw2();
  h_ele_PoPtrue_endcaps        = new TH1F( "h_ele_PoPtrue_endcaps",        "ele momentum / gen momentum, endcaps",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_endcaps->Sumw2();
  h_ele_PoPtrue_golden_barrel         = new TH1F( "h_ele_PoPtrue_golden_barrel",        "ele momentum / gen momentum, golden, barrel",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_golden_barrel->Sumw2();
  h_ele_PoPtrue_golden_endcaps        = new TH1F( "h_ele_PoPtrue_golden_endcaps",        "ele momentum / gen momentum, golden, endcaps",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_golden_endcaps->Sumw2();
  h_ele_PoPtrue_showering_barrel         = new TH1F( "h_ele_PoPtrue_showering_barrel",        "ele momentum / gen momentum, showering, barrel",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_showering_barrel->Sumw2();
  h_ele_PoPtrue_showering_endcaps        = new TH1F( "h_ele_PoPtrue_showering_endcaps",        "ele momentum / gen momentum, showering, endcaps",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PoPtrue_showering_endcaps->Sumw2();
  h_ele_PtoPttrue_barrel         = new TH1F( "h_ele_PtoPttrue_barrel",        "ele transverse momentum / gen transverse momentum, barrel",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PtoPttrue_barrel->Sumw2();
  h_ele_PtoPttrue_endcaps        = new TH1F( "h_ele_PtoPttrue_endcaps",        "ele transverse momentum / gen transverse momentum, endcaps",nbinpoptrue,poptruemin,poptruemax);
  h_ele_PtoPttrue_endcaps->Sumw2();
  h_ele_EtaMnEtaTrue   = new TH1F( "h_ele_EtaMnEtaTrue",   "ele momentum  eta - gen  eta",nbindeta,detamin,detamax);
  h_ele_EtaMnEtaTrue->Sumw2();
  h_ele_EtaMnEtaTrue_barrel   = new TH1F( "h_ele_EtaMnEtaTrue_barrel",   "ele momentum  eta - gen  eta barrel",nbindeta,detamin,detamax);
  h_ele_EtaMnEtaTrue_barrel->Sumw2();
  h_ele_EtaMnEtaTrue_endcaps   = new TH1F( "h_ele_EtaMnEtaTrue_endcaps",   "ele momentum  eta - gen  eta endcaps",nbindeta,detamin,detamax);
  h_ele_EtaMnEtaTrue_endcaps->Sumw2();
  h_ele_EtaMnEtaTrueVsEta   = new TH2F( "h_ele_EtaMnEtaTrueVsEta",   "ele momentum  eta - gen  eta vs eta",nbineta2D,etamin,etamax,nbindeta/2,detamin,detamax);
  h_ele_EtaMnEtaTrueVsPhi   = new TH2F( "h_ele_EtaMnEtaTrueVsPhi",   "ele momentum  eta - gen  eta vs phi",nbinphi2D,phimin,phimax,nbindeta/2,detamin,detamax);
  h_ele_EtaMnEtaTrueVsPt   = new TH2F( "h_ele_EtaMnEtaTrueVsPt",   "ele momentum  eta - gen  eta vs pt",nbinpt,0.,ptmax,nbindeta/2,detamin,detamax);
  h_ele_PhiMnPhiTrue   = new TH1F( "h_ele_PhiMnPhiTrue",   "ele momentum  phi - gen  phi",nbindphi,dphimin,dphimax);
  h_ele_PhiMnPhiTrue->Sumw2();
  h_ele_PhiMnPhiTrue_barrel   = new TH1F( "h_ele_PhiMnPhiTrue_barrel",   "ele momentum  phi - gen  phi barrel",nbindphi,dphimin,dphimax);
  h_ele_PhiMnPhiTrue_barrel->Sumw2();
  h_ele_PhiMnPhiTrue_endcaps   = new TH1F( "h_ele_PhiMnPhiTrue_endcaps",   "ele momentum  phi - gen  phi endcaps",nbindphi,dphimin,dphimax);
  h_ele_PhiMnPhiTrue_endcaps->Sumw2();
  h_ele_PhiMnPhiTrue2   = new TH1F( "h_ele_PhiMnPhiTrue2",   "ele momentum  phi - gen  phi",nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_PhiMnPhiTrueVsEta   = new TH2F( "h_ele_PhiMnPhiTrueVsEta",   "ele momentum  phi - gen  phi vs eta",nbineta2D,etamin,etamax,nbindphi/2,dphimin,dphimax);
  h_ele_PhiMnPhiTrueVsPhi   = new TH2F( "h_ele_PhiMnPhiTrueVsPhi",   "ele momentum  phi - gen  phi vs phi",nbinphi2D,phimin,phimax,nbindphi/2,dphimin,dphimax);
  h_ele_PhiMnPhiTrueVsPt   = new TH2F( "h_ele_PhiMnPhiTrueVsPt",   "ele momentum  phi - gen  phi vs pt",nbinpt2D,0.,ptmax,nbindphi/2,dphimin,dphimax);

  // matched electron, superclusters
  histSclEn_ = new TH1F("h_scl_energy","ele supercluster energy",nbinp,0.,pmax);
  histSclEn_->Sumw2();
  histSclEoEtrue_barrel = new TH1F("h_scl_EoEtrue_barrel","ele supercluster energy / gen energy, barrel",50,0.2,1.2);
  histSclEoEtrue_barrel->Sumw2();
  histSclEoEtrue_barrel_eg = new TH1F("h_scl_EoEtrue_barrel_eg","ele supercluster energy / gen energy, barrel, ecal driven",50,0.2,1.2);
  histSclEoEtrue_barrel_eg->Sumw2();
  histSclEoEtrue_barrel_etagap = new TH1F("h_scl_EoEtrue_barrel_etagap","ele supercluster energy / gen energy, barrel, etagap",50,0.2,1.2);
  histSclEoEtrue_barrel_etagap->Sumw2();
  histSclEoEtrue_barrel_phigap = new TH1F("h_scl_EoEtrue_barrel_phigap","ele supercluster energy / gen energy, barrel, phigap",50,0.2,1.2);
  histSclEoEtrue_barrel_phigap->Sumw2();
  histSclEoEtrue_ebeegap = new TH1F("h_scl_EoEtrue_ebeegap","ele supercluster energy / gen energy, ebeegap",50,0.2,1.2);
  histSclEoEtrue_ebeegap->Sumw2();
  histSclEoEtrue_endcaps = new TH1F("h_scl_EoEtrue_endcaps","ele supercluster energy / gen energy, endcaps",50,0.2,1.2);
  histSclEoEtrue_endcaps->Sumw2();
  histSclEoEtrue_endcaps_golden = new TH1F("h_scl_EoEtrue_endcaps_golden","ele supercluster energy / gen energy, endcaps",50,0.2,1.2);
  histSclEoEtrue_endcaps_golden->Sumw2();
  histSclEoEtrue_endcaps_narrow = new TH1F("h_scl_EoEtrue_endcaps_narrow","ele supercluster energy / gen energy, endcaps",50,0.2,1.2);
  histSclEoEtrue_endcaps_narrow->Sumw2();
  histSclEoEtrue_endcaps_bigbrem = new TH1F("h_scl_EoEtrue_endcaps_bigbrem","ele supercluster energy / gen energy, endcaps",50,0.2,1.2);
  histSclEoEtrue_endcaps_bigbrem->Sumw2();
  histSclEoEtrue_endcaps_showering = new TH1F("h_scl_EoEtrue_endcaps_showering","ele supercluster energy / gen energy, endcaps",50,0.2,1.2);
  histSclEoEtrue_endcaps_showering->Sumw2();
  histSeedEoEtrue_endcaps_golden = new TH1F("h_seed_EoEtrue_endcaps_golden","ele seed energy / gen energy, endcaps",50,0.2,1.2);
  histSeedEoEtrue_endcaps_golden->Sumw2();
  histSeedEoEtrue_endcaps_narrow = new TH1F("h_seed_EoEtrue_endcaps_narrow","ele seed energy / gen energy, endcaps",50,0.2,1.2);
  histSeedEoEtrue_endcaps_narrow->Sumw2();
  histSeedEoEtrue_endcaps_bigbrem = new TH1F("h_seed_EoEtrue_endcaps_bigbrem","ele seed energy / gen energy, endcaps",50,0.2,1.2);
  histSeedEoEtrue_endcaps_bigbrem->Sumw2();
  histSeedEoEtrue_endcaps_showering = new TH1F("h_seed_EoEtrue_endcaps_showering","ele seed energy / gen energy, endcaps",50,0.2,1.2);
  histSeedEoEtrue_endcaps_showering->Sumw2();
  histSclEoEtrue_endcaps_eg = new TH1F("h_scl_EoEtrue_endcaps_eg","ele supercluster energy / gen energy, endcaps, ecal driven",50,0.2,1.2);
  histSclEoEtrue_endcaps_eg->Sumw2();
  histSclEoEtrue_endcaps_deegap = new TH1F("h_scl_EoEtrue_endcaps_deegap","ele supercluster energy / gen energy, endcaps, deegap",50,0.2,1.2);
  histSclEoEtrue_endcaps_deegap->Sumw2();
  histSclEoEtrue_endcaps_ringgap = new TH1F("h_scl_EoEtrue_endcaps_ringgap","ele supercluster energy / gen energy, endcaps, ringgap",50,0.2,1.2);
  histSclEoEtrue_endcaps_ringgap->Sumw2();
  histSclEoEtrue_barrel_new = new TH1F("h_scl_EoEtrue_barrel_new","ele supercluster energy / gen energy, barrel",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_barrel_new->Sumw2();
  histSclEoEtrue_barrel_eg_new = new TH1F("h_scl_EoEtrue_barrel_eg_new","ele supercluster energy / gen energy, barrel, ecal driven",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_barrel_eg_new->Sumw2();
  histSclEoEtrue_barrel_etagap_new = new TH1F("h_scl_EoEtrue_barrel_etagap_new","ele supercluster energy / gen energy, barrel, etagap",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_barrel_etagap_new->Sumw2();
  histSclEoEtrue_barrel_phigap_new = new TH1F("h_scl_EoEtrue_barrel_phigap_new","ele supercluster energy / gen energy, barrel, phigap",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_barrel_phigap_new->Sumw2();
  histSclEoEtrue_ebeegap_new = new TH1F("h_scl_EoEtrue_ebeegap_new","ele supercluster energy / gen energy, ebeegap",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_ebeegap_new->Sumw2();
  histSclEoEtrue_endcaps_new = new TH1F("h_scl_EoEtrue_endcaps_new","ele supercluster energy / gen energy, endcaps",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_endcaps_new->Sumw2();
  histSclEoEtrue_endcaps_eg_new = new TH1F("h_scl_EoEtrue_endcaps_eg_new","ele supercluster energy / gen energy, endcaps, ecal driven",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_endcaps_eg_new->Sumw2();
  histSclEoEtrue_endcaps_deegap_new = new TH1F("h_scl_EoEtrue_endcaps_deegap_new","ele supercluster energy / gen energy, endcaps, deegap",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_endcaps_deegap_new->Sumw2();
  histSclEoEtrue_endcaps_ringgap_new = new TH1F("h_scl_EoEtrue_endcaps_ringgap_new","ele supercluster energy / gen energy, endcaps, ringgap",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrue_endcaps_ringgap_new->Sumw2();
  histSclEt_ = new TH1F("h_scl_et","ele supercluster transverse energy",nbinpt,0.,ptmax);
  histSclEt_->Sumw2();
  histSclEtVsEta_ = new TH2F("h_scl_etVsEta","ele supercluster transverse energy vs eta",nbineta2D,etamin,etamax,nbinpt,0.,ptmax);
  histSclEtVsPhi_ = new TH2F("h_scl_etVsPhi","ele supercluster transverse energy vs phi",nbinphi2D,phimin,phimax,nbinpt,0.,ptmax);
  histSclEtaVsPhi_ = new TH2F("h_scl_etaVsPhi","ele supercluster eta vs phi",nbinphi2D,phimin,phimax,nbineta2D,etamin,etamax);
  histSclEta_ = new TH1F("h_scl_eta","ele supercluster eta",nbineta,etamin,etamax);
  histSclEta_->Sumw2();
  histSclPhi_ = new TH1F("h_scl_phi","ele supercluster phi",nbinphi,phimin,phimax);
  histSclPhi_->Sumw2();

  histSclSigEtaEta_ =  new TH1F("h_scl_sigetaeta","ele supercluster sigma eta eta",100,0.,0.05);
  histSclSigEtaEta_->Sumw2();
  histSclSigEtaEta_barrel_ =  new TH1F("h_scl_sigetaeta_barrel","ele supercluster sigma eta eta barrel",100,0.,0.05);
  histSclSigEtaEta_barrel_->Sumw2();
  histSclSigEtaEta_endcaps_ =  new TH1F("h_scl_sigetaeta_endcaps","ele supercluster sigma eta eta endcaps",100,0.,0.05);
  histSclSigEtaEta_endcaps_->Sumw2();
  histSclSigIEtaIEta_ =  new TH1F("h_scl_sigietaieta","ele supercluster sigma ieta ieta",100,0.,0.05);
  histSclSigIEtaIEta_->Sumw2();
  histSclSigIEtaIEta_barrel_ =  new TH1F("h_scl_sigietaieta_barrel","ele supercluster sigma ieta ieta, barrel",100,0.,0.05);
  histSclSigIEtaIEta_barrel_->Sumw2();
  histSclSigIEtaIEta_endcaps_ =  new TH1F("h_scl_sigietaieta_endcaps","ele supercluster sigma ieta ieta, endcaps",100,0.,0.05);
  histSclSigIEtaIEta_endcaps_->Sumw2();
  histSclE1x5_ =  new TH1F("h_scl_E1x5","ele supercluster energy in 1x5",nbinp,0., pmax);
  histSclE1x5_->Sumw2();
  histSclE1x5_barrel_ =  new TH1F("h_scl_E1x5_barrel","ele supercluster energy in 1x5 barrel",nbinp,0., pmax);
  histSclE1x5_barrel_->Sumw2();
  histSclE1x5_endcaps_ =  new TH1F("h_scl_E1x5_endcaps","ele supercluster energy in 1x5 endcaps",nbinp,0., pmax);
  histSclE1x5_endcaps_->Sumw2();
  histSclE2x5max_ =  new TH1F("h_scl_E2x5max","ele supercluster energy in 2x5 max",nbinp,0.,pmax);
  histSclE2x5max_->Sumw2();
  histSclE2x5max_barrel_ =  new TH1F("h_scl_E2x5max_barrel","ele supercluster energy in 2x5 max barrel",nbinp,0.,pmax);
  histSclE2x5max_barrel_->Sumw2();
  histSclE2x5max_endcaps_ =  new TH1F("h_scl_E2x5max_endcaps","ele supercluster energy in 2x5 max endcaps",nbinp,0.,pmax);
  histSclE2x5max_endcaps_->Sumw2();
  histSclE5x5_ =  new TH1F("h_scl_E5x5","ele supercluster energy in 5x5",nbinp,0.,pmax);
  histSclE5x5_->Sumw2();
  histSclE5x5_barrel_ =  new TH1F("h_scl_E5x5_barrel","ele supercluster energy in 5x5 barrel",nbinp,0.,pmax);
  histSclE5x5_barrel_->Sumw2();
  histSclE5x5_endcaps_ =  new TH1F("h_scl_E5x5_endcaps","ele supercluster energy in 5x5 endcaps",nbinp,0.,pmax);
  histSclE5x5_endcaps_->Sumw2();
  histSclSigEtaEta_eg_ =  new TH1F("h_scl_sigetaeta_eg","ele supercluster sigma eta eta, ecal driven",100,0.,0.05);
  histSclSigEtaEta_eg_->Sumw2();
  histSclSigEtaEta_eg_barrel_ =  new TH1F("h_scl_sigetaeta_eg_barrel","ele supercluster sigma eta eta, ecal driven barrel",100,0.,0.05);
  histSclSigEtaEta_eg_barrel_->Sumw2();
  histSclSigEtaEta_eg_endcaps_ =  new TH1F("h_scl_sigetaeta_eg_endcaps","ele supercluster sigma eta eta, ecal driven endcaps",100,0.,0.05);
  histSclSigEtaEta_eg_endcaps_->Sumw2();
  histSclSigIEtaIEta_eg_ =  new TH1F("h_scl_sigietaieta_eg","ele supercluster sigma ieta ieta, ecal driven",100,0.,0.05);
  histSclSigIEtaIEta_eg_->Sumw2();
  histSclSigIEtaIEta_eg_barrel_ =  new TH1F("h_scl_sigietaieta_barrel_eg","ele supercluster sigma ieta ieta, barrel, ecal driven",100,0.,0.05);
  histSclSigIEtaIEta_eg_barrel_->Sumw2();
  histSclSigIEtaIEta_eg_endcaps_ =  new TH1F("h_scl_sigietaieta_endcaps_eg","ele supercluster sigma ieta ieta, endcaps, ecal driven",100,0.,0.05);
  histSclSigIEtaIEta_eg_endcaps_->Sumw2();
  histSclE1x5_eg_ =  new TH1F("h_scl_E1x5_eg","ele supercluster energy in 1x5, ecal driven",nbinp,0., pmax);
  histSclE1x5_eg_->Sumw2();
  histSclE1x5_eg_barrel_ =  new TH1F("h_scl_E1x5_eg_barrel","ele supercluster energy in 1x5, ecal driven barrel",nbinp,0., pmax);
  histSclE1x5_eg_barrel_->Sumw2();
  histSclE1x5_eg_endcaps_ =  new TH1F("h_scl_E1x5_eg_endcaps","ele supercluster energy in 1x5, ecal driven endcaps",nbinp,0., pmax);
  histSclE1x5_eg_endcaps_->Sumw2();
  histSclE2x5max_eg_ =  new TH1F("h_scl_E2x5max_eg","ele supercluster energy in 2x5 max, ecal driven",nbinp,0.,pmax);
  histSclE2x5max_eg_->Sumw2();
  histSclE2x5max_eg_barrel_ =  new TH1F("h_scl_E2x5max_eg_barrel","ele supercluster energy in 2x5 max, ecal driven barrel",nbinp,0.,pmax);
  histSclE2x5max_eg_barrel_->Sumw2();
  histSclE2x5max_eg_endcaps_ =  new TH1F("h_scl_E2x5max_eg_endcaps","ele supercluster energy in 2x5 max, ecal driven endcaps",nbinp,0.,pmax);
  histSclE2x5max_eg_endcaps_->Sumw2();
  histSclE5x5_eg_ =  new TH1F("h_scl_E5x5_eg","ele supercluster energy in 5x5, ecal driven",nbinp,0.,pmax);
  histSclE5x5_eg_->Sumw2();
  histSclE5x5_eg_barrel_ =  new TH1F("h_scl_E5x5_eg_barrel","ele supercluster energy in 5x5, ecal driven barrel",nbinp,0.,pmax);
  histSclE5x5_eg_barrel_->Sumw2();
  histSclE5x5_eg_endcaps_ =  new TH1F("h_scl_E5x5_eg_endcaps","ele supercluster energy in 5x5, ecal driven endcaps",nbinp,0.,pmax);
  histSclE5x5_eg_endcaps_->Sumw2();

  histSclEoEtruePfVsEg = new TH2F("h_scl_EoEtruePfVsEg","ele supercluster energy / gen energy pflow vs eg",75,-0.1,1.4, 75, -0.1, 1.4);

  // matched electron, gsf tracks
  h_ele_ambiguousTracks      = new TH1F( "h_ele_ambiguousTracks", "ele # ambiguous tracks",  5,0.,5.);
  h_ele_ambiguousTracks->Sumw2();
  h_ele_ambiguousTracksVsEta      = new TH2F( "h_ele_ambiguousTracksVsEta","ele # ambiguous tracks  vs eta",  nbineta2D,etamin,etamax,5,0.,5.);
  h_ele_ambiguousTracksVsPhi      = new TH2F( "h_ele_ambiguousTracksVsPhi", "ele # ambiguous tracks  vs phi",  nbinphi2D,phimin,phimax,5,0.,5.);
  h_ele_ambiguousTracksVsPt      = new TH2F( "h_ele_ambiguousTracksVsPt", "ele # ambiguous tracks vs pt",  nbinpt2D,0.,ptmax,5,0.,5.);
  h_ele_foundHits      = new TH1F( "h_ele_foundHits",      "ele track # found hits",      nbinfhits,0.,fhitsmax);
  h_ele_foundHits->Sumw2();
  h_ele_foundHits_barrel      = new TH1F( "h_ele_foundHits_barrel",      "ele track # found hits, barrel",      nbinfhits,0.,fhitsmax);
  h_ele_foundHits_barrel->Sumw2();
  h_ele_foundHits_endcaps      = new TH1F( "h_ele_foundHits_endcaps",      "ele track # found hits, endcaps",      nbinfhits,0.,fhitsmax);
  h_ele_foundHits_endcaps->Sumw2();
  h_ele_foundHitsVsEta      = new TH2F( "h_ele_foundHitsVsEta",      "ele track # found hits vs eta",  nbineta2D,etamin,etamax,nbinfhits,0.,fhitsmax);
  h_ele_foundHitsVsPhi      = new TH2F( "h_ele_foundHitsVsPhi",      "ele track # found hits vs phi",  nbinphi2D,phimin,phimax,nbinfhits,0.,fhitsmax);
  h_ele_foundHitsVsPt      = new TH2F( "h_ele_foundHitsVsPt",      "ele track # found hits vs pt",  nbinpt2D,0.,ptmax,nbinfhits,0.,fhitsmax);
  h_ele_lostHits       = new TH1F( "h_ele_lostHits",       "ele track # lost hits",       5,0.,5.);
  h_ele_lostHits->Sumw2();
  h_ele_lostHits_barrel       = new TH1F( "h_ele_lostHits_barrel",       "ele track # lost hits, barrel",       5,0.,5.);
  h_ele_lostHits_barrel->Sumw2();
  h_ele_lostHits_endcaps       = new TH1F( "h_ele_lostHits_endcaps",       "ele track # lost hits, endcaps",       5,0.,5.);
  h_ele_lostHits_endcaps->Sumw2();
  h_ele_lostHitsVsEta       = new TH2F( "h_ele_lostHitsVsEta",       "ele track # lost hits vs eta",   nbineta2D,etamin,etamax,nbinlhits,0.,lhitsmax);
  h_ele_lostHitsVsPhi       = new TH2F( "h_ele_lostHitsVsPhi",       "ele track # lost hits vs eta",   nbinphi2D,phimin,phimax,nbinlhits,0.,lhitsmax);
  h_ele_lostHitsVsPt       = new TH2F( "h_ele_lostHitsVsPt",       "ele track # lost hits vs eta",   nbinpt2D,0.,ptmax,nbinlhits,0.,lhitsmax);
  h_ele_chi2           = new TH1F( "h_ele_chi2",           "ele track #chi^{2}",         100,0.,15.);
  h_ele_chi2->Sumw2();
  h_ele_chi2_barrel           = new TH1F( "h_ele_chi2_barrel",           "ele track #chi^{2}, barrel",         100,0.,15.);
  h_ele_chi2_barrel->Sumw2();
  h_ele_chi2_endcaps           = new TH1F( "h_ele_chi2_endcaps",           "ele track #chi^{2}, endcaps",         100,0.,15.);
  h_ele_chi2_endcaps->Sumw2();
  h_ele_chi2VsEta           = new TH2F( "h_ele_chi2VsEta",           "ele track #chi^{2} vs eta",  nbineta2D,etamin,etamax,50,0.,15.);
  h_ele_chi2VsPhi           = new TH2F( "h_ele_chi2VsPhi",           "ele track #chi^{2} vs phi",  nbinphi2D,phimin,phimax,50,0.,15.);
  h_ele_chi2VsPt           = new TH2F( "h_ele_chi2VsPt",           "ele track #chi^{2} vs pt",  nbinpt2D,0.,ptmax,50,0.,15.);
  h_ele_PinMnPout      = new TH1F( "h_ele_PinMnPout",      "ele track inner p - outer p, mean of GSF components"   ,nbinp,0.,200.);
  h_ele_PinMnPout->Sumw2();
  h_ele_PinMnPout_mode      = new TH1F( "h_ele_PinMnPout_mode",      "ele track inner p - outer p, mode of GSF components"   ,nbinp,0.,100.);
  h_ele_PinMnPout_mode->Sumw2();
  h_ele_PinMnPoutVsEta_mode = new TH2F( "h_ele_PinMnPoutVsEta_mode",      "ele track inner p - outer p vs eta, mode of GSF components" ,nbineta2D, etamin,etamax,nbinp2D,0.,100.);
  h_ele_PinMnPoutVsPhi_mode = new TH2F( "h_ele_PinMnPoutVsPhi_mode",      "ele track inner p - outer p vs phi, mode of GSF components" ,nbinphi2D, phimin,phimax,nbinp2D,0.,100.);
  h_ele_PinMnPoutVsPt_mode = new TH2F( "h_ele_PinMnPoutVsPt_mode",      "ele track inner p - outer p vs pt, mode of GSF components" ,nbinpt2D, 0.,ptmax,nbinp2D,0.,100.);
  h_ele_PinMnPoutVsE_mode = new TH2F( "h_ele_PinMnPoutVsE_mode",      "ele track inner p - outer p vs E, mode of GSF components" ,nbinp2D, 0.,200.,nbinp2D,0.,100.);
  h_ele_PinMnPoutVsChi2_mode = new TH2F( "h_ele_PinMnPoutVsChi2_mode",      "ele track inner p - outer p vs track chi2, mode of GSF components" ,50, 0.,20.,nbinp2D,0.,100.);
  h_ele_outerP         = new TH1F( "h_ele_outerP",         "ele track outer p, mean of GSF components",          nbinp,0.,pmax);
  h_ele_outerP->Sumw2();
  h_ele_outerP_mode         = new TH1F( "h_ele_outerP_mode",         "ele track outer p, mode of GSF components",          nbinp,0.,pmax);
  h_ele_outerP_mode->Sumw2();
  h_ele_outerPVsEta_mode         = new TH2F( "h_ele_outerPVsEta_mode",         "ele track outer p vs eta mode", nbineta2D,etamin,etamax,50,0.,pmax);
  h_ele_outerPt        = new TH1F( "h_ele_outerPt",        "ele track outer p_{T}, mean of GSF components",      nbinpt,0.,ptmax);
  h_ele_outerPt->Sumw2();
  h_ele_outerPt_mode        = new TH1F( "h_ele_outerPt_mode",        "ele track outer p_{T}, mode of GSF components",      nbinpt,0.,ptmax);
  h_ele_outerPt_mode->Sumw2();
  h_ele_outerPtVsEta_mode        = new TH2F( "h_ele_outerPtVsEta_mode", "ele track outer p_{T} vs eta, mode of GSF components", nbineta2D,etamin,etamax,nbinpt2D,0.,ptmax);
  h_ele_outerPtVsPhi_mode        = new TH2F( "h_ele_outerPtVsPhi_mode", "ele track outer p_{T} vs phi, mode of GSF components", nbinphi2D,phimin,phimax,nbinpt2D,0.,ptmax);
  h_ele_outerPtVsPt_mode        = new TH2F( "h_ele_outerPtVsPt_mode", "ele track outer p_{T} vs pt, mode of GSF components", nbinpt2D,0.,100.,nbinpt2D,0.,ptmax);

  // matched electrons, matching
  h_ele_EoP            = new TH1F( "h_ele_EoP",            "ele E/P_{vertex}",        nbineop,0.,eopmax);
  h_ele_EoP->Sumw2();
  h_ele_EoP_eg            = new TH1F( "h_ele_EoP_eg",            "ele E/P_{vertex}, ecal driven",        nbineop,0.,eopmax);
  h_ele_EoP_eg->Sumw2();
  h_ele_EoP_barrel            = new TH1F( "h_ele_EoP_barrel",            "ele E/P_{vertex} barrel",        nbineop,0.,eopmax);
  h_ele_EoP_barrel->Sumw2();
  h_ele_EoP_eg_barrel            = new TH1F( "h_ele_EoP_eg_barrel",            "ele E/P_{vertex}, ecal driven barrel",        nbineop,0.,eopmax);
  h_ele_EoP_eg_barrel->Sumw2();
  h_ele_EoP_endcaps            = new TH1F( "h_ele_EoP_endcaps",            "ele E/P_{vertex} endcaps",        nbineop,0.,eopmax);
  h_ele_EoP_endcaps->Sumw2();
  h_ele_EoP_eg_endcaps            = new TH1F( "h_ele_EoP_eg_endcaps",            "ele E/P_{vertex}, ecal driven endcaps",        nbineop,0.,eopmax);
  h_ele_EoP_eg_endcaps->Sumw2();
  h_ele_EoPVsEta            = new TH2F( "h_ele_EoPVsEta",            "ele E/P_{vertex} vs eta",  nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPVsPhi            = new TH2F( "h_ele_EoPVsPhi",            "ele E/P_{vertex} vs phi",  nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPVsE            = new TH2F( "h_ele_EoPVsE",            "ele E/P_{vertex} vs E",  50,0.,pmax ,50,0.,5.);
  h_ele_EseedOP            = new TH1F( "h_ele_EseedOP",            "ele E_{seed}/P_{vertex}",        nbineop,0.,eopmax);
  h_ele_EseedOP->Sumw2();
  h_ele_EseedOP_eg            = new TH1F( "h_ele_EseedOP_eg",            "ele E_{seed}/P_{vertex}, ecal driven",        nbineop,0.,eopmax);
  h_ele_EseedOP_eg->Sumw2();
  h_ele_EseedOP_barrel            = new TH1F( "h_ele_EseedOP_barrel",            "ele E_{seed}/P_{vertex} barrel",        nbineop,0.,eopmax);
  h_ele_EseedOP_barrel->Sumw2();
  h_ele_EseedOP_eg_barrel            = new TH1F( "h_ele_EseedOP_eg_barrel",            "ele E_{seed}/P_{vertex}, ecal driven barrel",        nbineop,0.,eopmax);
  h_ele_EseedOP_eg_barrel->Sumw2();
  h_ele_EseedOP_endcaps            = new TH1F( "h_ele_EseedOP_endcaps",            "ele E_{seed}/P_{vertex} endcaps",        nbineop,0.,eopmax);
  h_ele_EseedOP_endcaps->Sumw2();
  h_ele_EseedOP_eg_endcaps            = new TH1F( "h_ele_EseedOP_eg_endcaps",            "ele E_{seed}/P_{vertex}, ecal driven, endcaps",        nbineop,0.,eopmax);
  h_ele_EseedOP_eg_endcaps->Sumw2();
  h_ele_EseedOPVsEta            = new TH2F( "h_ele_EseedOPVsEta",            "ele E_{seed}/P_{vertex} vs eta",  nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EseedOPVsPhi            = new TH2F( "h_ele_EseedOPVsPhi",            "ele E_{seed}/P_{vertex} vs phi",  nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EseedOPVsE            = new TH2F( "h_ele_EseedOPVsE",            "ele E_{seed}/P_{vertex} vs E",  50,0.,pmax ,50,0.,5.);
  h_ele_EoPout         = new TH1F( "h_ele_EoPout",         "ele E_{seed}/P_{out}",           nbineop,0.,eopmax);
  h_ele_EoPout->Sumw2();
  h_ele_EoPout_eg         = new TH1F( "h_ele_EoPout_eg",         "ele E_{seed}/P_{out}, ecal driven",           nbineop,0.,eopmax);
  h_ele_EoPout_eg->Sumw2();
  h_ele_EoPout_barrel         = new TH1F( "h_ele_EoPout_barrel",         "ele E_{seed}/P_{out} barrel",           nbineop,0.,eopmax);
  h_ele_EoPout_barrel->Sumw2();
  h_ele_EoPout_eg_barrel         = new TH1F( "h_ele_EoPout_eg_barrel",         "ele E_{seed}/P_{out}, ecal driven, barrel",           nbineop,0.,eopmax);
  h_ele_EoPout_eg_barrel->Sumw2();
  h_ele_EoPout_endcaps         = new TH1F( "h_ele_EoPout_endcaps",         "ele E_{seed}/P_{out} endcaps",           nbineop,0.,eopmax);
  h_ele_EoPout_endcaps->Sumw2();
  h_ele_EoPout_endcaps_golden         = new TH1F( "h_ele_EoPout_endcaps_golden",         "ele E_{seed}/P_{out} endcaps",           nbineop,0.,eopmax);
  h_ele_EoPout_endcaps_golden->Sumw2();
  h_ele_EoPout_endcaps_showering         = new TH1F( "h_ele_EoPout_endcaps_showering",         "ele E_{seed}/P_{out} endcaps",           nbineop,0.,eopmax);
  h_ele_EoPout_endcaps_showering->Sumw2();
  h_ele_EoPout_eg_endcaps         = new TH1F( "h_ele_EoPout_eg_endcaps",         "ele E_{seed}/P_{out}, ecal driven, endcaps",           nbineop,0.,eopmax);
  h_ele_EoPout_eg_endcaps->Sumw2();
  h_ele_EoPoutVsEta         = new TH2F( "h_ele_EoPoutVsEta",         "ele E_{seed}/P_{out} vs eta",    nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsPhi         = new TH2F( "h_ele_EoPoutVsPhi",         "ele E_{seed}/P_{out} vs phi",    nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsE         = new TH2F( "h_ele_EoPoutVsE",         "ele E_{seed}/P_{out} vs E",    nbinp2D,0.,pmax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsEta_golden         = new TH2F( "h_ele_EoPoutVsEta_golden",         "ele E_{seed}/P_{out} vs eta",    nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsPhi_golden         = new TH2F( "h_ele_EoPoutVsPhi_golden",         "ele E_{seed}/P_{out} vs phi",    nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsE_golden         = new TH2F( "h_ele_EoPoutVsE_golden",         "ele E_{seed}/P_{out} vs E",    nbinp2D,0.,pmax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsEta_showering         = new TH2F( "h_ele_EoPoutVsEta_showering",         "ele E_{seed}/P_{out} vs eta",    nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsPhi_showering         = new TH2F( "h_ele_EoPoutVsPhi_showering",         "ele E_{seed}/P_{out} vs phi",    nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EoPoutVsE_showering         = new TH2F( "h_ele_EoPoutVsE_showering",         "ele E_{seed}/P_{out} vs E",    nbinp2D,0.,pmax,nbineop2D,0.,eopmaxsht);
  h_ele_EeleOPout         = new TH1F( "h_ele_EeleOPout",         "ele E_{ele}/P_{out}",           nbineop,0.,eopmax);
  h_ele_EeleOPout->Sumw2();
  h_ele_EeleOPout_eg         = new TH1F( "h_ele_EeleOPout_eg",         "ele E_{ele}/P_{out}, ecal driven",           nbineop,0.,eopmax);
  h_ele_EeleOPout_eg->Sumw2();
  h_ele_EeleOPout_barrel         = new TH1F( "h_ele_EeleOPout_barrel",         "ele E_{ele}/P_{out} barrel",           nbineop,0.,eopmax);
  h_ele_EeleOPout_barrel->Sumw2();
  h_ele_EeleOPout_eg_barrel         = new TH1F( "h_ele_EeleOPout_eg_barrel",         "ele E_{ele}/P_{out}, ecal driven, barrel",           nbineop,0.,eopmax);
  h_ele_EeleOPout_eg_barrel->Sumw2();
  h_ele_EeleOPout_endcaps         = new TH1F( "h_ele_EeleOPout_endcaps",         "ele E_{ele}/P_{out} endcaps",           nbineop,0.,eopmax);
  h_ele_EeleOPout_endcaps->Sumw2();
  h_ele_EeleOPout_eg_endcaps         = new TH1F( "h_ele_EeleOPout_eg_endcaps",         "ele E_{ele}/P_{out}, ecal driven, endcaps",           nbineop,0.,eopmax);
  h_ele_EeleOPout_eg_endcaps->Sumw2();
  h_ele_EeleOPoutVsEta         = new TH2F( "h_ele_EeleOPoutVsEta",         "ele E_{ele}/P_{out} vs eta",    nbineta2D,etamin,etamax,nbineop2D,0.,eopmaxsht);
  h_ele_EeleOPoutVsPhi         = new TH2F( "h_ele_EeleOPoutVsPhi",         "ele E_{ele}/P_{out} vs phi",    nbinphi2D,phimin,phimax,nbineop2D,0.,eopmaxsht);
  h_ele_EeleOPoutVsE         = new TH2F( "h_ele_EeleOPoutVsE",         "ele E_{ele}/P_{out} vs E",    nbinp2D,0.,pmax,nbineop2D,0.,eopmaxsht);
  h_ele_dEtaSc_propVtx = new TH1F( "h_ele_dEtaSc_propVtx", "ele #eta_{sc} - #eta_{tr}, prop from vertex",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx->Sumw2();
  h_ele_dEtaSc_propVtx_eg = new TH1F( "h_ele_dEtaSc_propVtx_eg", "ele #eta_{sc} - #eta_{tr}, prop from vertex, ecal driven",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_eg->Sumw2();
  h_ele_dEtaSc_propVtx_barrel = new TH1F( "h_ele_dEtaSc_propVtx_barrel", "ele #eta_{sc} - #eta_{tr}, prop from vertex, barrel",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_barrel->Sumw2();
  h_ele_dEtaSc_propVtx_eg_barrel = new TH1F( "h_ele_dEtaSc_propVtx_eg_barrel", "ele #eta_{sc} - #eta_{tr}, prop from vertex, ecal driven, barrel",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_eg_barrel->Sumw2();
  h_ele_dEtaSc_propVtx_endcaps = new TH1F( "h_ele_dEtaSc_propVtx_endcaps", "ele #eta_{sc} - #eta_{tr}, prop from vertex, endcaps",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_endcaps->Sumw2();
  h_ele_dEtaSc_propVtx_eg_endcaps = new TH1F( "h_ele_dEtaSc_propVtx_eg_endcaps", "ele #eta_{sc} - #eta_{tr}, prop from vertex, ecal driven, endcaps",      nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaSc_propVtx_eg_endcaps->Sumw2();
  h_ele_dEtaScVsEta_propVtx = new TH2F( "h_ele_dEtaScVsEta_propVtx", "ele #eta_{sc} - #eta_{tr} vs eta, prop from vertex", nbineta2D,etamin,etamax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaScVsPhi_propVtx = new TH2F( "h_ele_dEtaScVsPhi_propVtx", "ele #eta_{sc} - #eta_{tr} vs phi, prop from vertex", nbinphi2D,phimin,phimax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaScVsPt_propVtx = new TH2F( "h_ele_dEtaScVsPt_propVtx", "ele #eta_{sc} - #eta_{tr} vs pt, prop from vertex", nbinpt2D,0.,ptmax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dPhiSc_propVtx = new TH1F( "h_ele_dPhiSc_propVtx", "ele #phi_{sc} - #phi_{tr}, prop from vertex",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx->Sumw2();
  h_ele_dPhiSc_propVtx_eg = new TH1F( "h_ele_dPhiSc_propVtx_eg", "ele #phi_{sc} - #phi_{tr}, prop from vertex, ecal driven",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_eg->Sumw2();
  h_ele_dPhiSc_propVtx_barrel = new TH1F( "h_ele_dPhiSc_propVtx_barrel", "ele #phi_{sc} - #phi_{tr}, prop from vertex, barrel",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_barrel->Sumw2();
  h_ele_dPhiSc_propVtx_eg_barrel = new TH1F( "h_ele_dPhiSc_propVtx_eg_barrel", "ele #phi_{sc} - #phi_{tr}, prop from vertex, ecal driven, barrel",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_eg_barrel->Sumw2();
  h_ele_dPhiSc_propVtx_endcaps = new TH1F( "h_ele_dPhiSc_propVtx_endcaps", "ele #phi_{sc} - #phi_{tr}, prop from vertex, endcaps",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_endcaps->Sumw2();
  h_ele_dPhiSc_propVtx_eg_endcaps = new TH1F( "h_ele_dPhiSc_propVtx_eg_endcaps", "ele #phi_{sc} - #phi_{tr}, prop from vertex, ecal driven, endcaps",      nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiSc_propVtx_eg_endcaps->Sumw2();
  h_ele_dPhiScVsEta_propVtx = new TH2F( "h_ele_dPhiScVsEta_propVtx", "ele #phi_{sc} - #phi_{tr} vs eta, prop from vertex", nbineta2D,etamin,etamax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiScVsPhi_propVtx = new TH2F( "h_ele_dPhiScVsPhi_propVtx", "ele #phi_{sc} - #phi_{tr} vs phi, prop from vertex", nbinphi2D,phimin,phimax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiScVsPt_propVtx = new TH2F( "h_ele_dPhiScVsPt_propVtx", "ele #phi_{sc} - #phi_{tr} vs pt, prop from vertex", nbinpt2D,0.,ptmax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dEtaCl_propOut = new TH1F( "h_ele_dEtaCl_propOut", "ele #eta_{cl} - #eta_{tr}, prop from outermost",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut->Sumw2();
  h_ele_dEtaCl_propOut_eg = new TH1F( "h_ele_dEtaCl_propOut_eg", "ele #eta_{cl} - #eta_{tr}, prop from outermost, ecal driven",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_eg->Sumw2();
  h_ele_dEtaCl_propOut_barrel = new TH1F( "h_ele_dEtaCl_propOut_barrel", "ele #eta_{cl} - #eta_{tr}, prop from outermost, barrel",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_barrel->Sumw2();
  h_ele_dEtaCl_propOut_eg_barrel = new TH1F( "h_ele_dEtaCl_propOut_eg_barrel", "ele #eta_{cl} - #eta_{tr}, prop from outermost, ecal driven, barrel",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_eg_barrel->Sumw2();
  h_ele_dEtaCl_propOut_endcaps = new TH1F( "h_ele_dEtaCl_propOut_endcaps", "ele #eta_{cl} - #eta_{tr}, prop from outermost, endcaps",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_endcaps->Sumw2();
  h_ele_dEtaCl_propOut_eg_endcaps = new TH1F( "h_ele_dEtaCl_propOut_eg_endcaps", "ele #eta_{cl} - #eta_{tr}, prop from outermost, ecal driven, endcaps",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaCl_propOut_eg_endcaps->Sumw2();
  h_ele_dEtaClVsEta_propOut = new TH2F( "h_ele_dEtaClVsEta_propOut", "ele #eta_{cl} - #eta_{tr} vs eta, prop from out", nbineta2D,etamin,etamax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaClVsPhi_propOut = new TH2F( "h_ele_dEtaClVsPhi_propOut", "ele #eta_{cl} - #eta_{tr} vs phi, prop from out", nbinphi2D,phimin,phimax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaClVsPt_propOut = new TH2F( "h_ele_dEtaScVsPt_propOut", "ele #eta_{cl} - #eta_{tr} vs pt, prop from out", nbinpt2D,0.,ptmax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dPhiCl_propOut = new TH1F( "h_ele_dPhiCl_propOut", "ele #phi_{cl} - #phi_{tr}, prop from outermost",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut->Sumw2();
  h_ele_dPhiCl_propOut_eg = new TH1F( "h_ele_dPhiCl_propOut_eg", "ele #phi_{cl} - #phi_{tr}, prop from outermost, ecal driven",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_eg->Sumw2();
  h_ele_dPhiCl_propOut_barrel = new TH1F( "h_ele_dPhiCl_propOut_barrel", "ele #phi_{cl} - #phi_{tr}, prop from outermost, barrel",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_barrel->Sumw2();
  h_ele_dPhiCl_propOut_eg_barrel = new TH1F( "h_ele_dPhiCl_propOut_eg_barrel", "ele #phi_{cl} - #phi_{tr}, prop from outermost, ecal driven, barrel",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_eg_barrel->Sumw2();
  h_ele_dPhiCl_propOut_endcaps = new TH1F( "h_ele_dPhiCl_propOut_endcaps", "ele #phi_{cl} - #phi_{tr}, prop from outermost, endcaps",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_endcaps->Sumw2();
  h_ele_dPhiCl_propOut_eg_endcaps = new TH1F( "h_ele_dPhiCl_propOut_eg_endcaps", "ele #phi_{cl} - #phi_{tr}, prop from outermost, ecal driven, endcaps",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiCl_propOut_eg_endcaps->Sumw2();
  h_ele_dPhiClVsEta_propOut = new TH2F( "h_ele_dPhiClVsEta_propOut", "ele #phi_{cl} - #phi_{tr} vs eta, prop from out", nbineta2D,etamin,etamax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiClVsPhi_propOut = new TH2F( "h_ele_dPhiClVsPhi_propOut", "ele #phi_{cl} - #phi_{tr} vs phi, prop from out", nbinphi2D,phimin,phimax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiClVsPt_propOut = new TH2F( "h_ele_dPhiSClsPt_propOut", "ele #phi_{cl} - #phi_{tr} vs pt, prop from out", nbinpt2D,0.,ptmax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dEtaEleCl_propOut = new TH1F( "h_ele_dEtaEleCl_propOut", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut->Sumw2();
  h_ele_dEtaEleCl_propOut_eg = new TH1F( "h_ele_dEtaEleCl_propOut_eg", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost, ecal driven",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut_eg->Sumw2();
  h_ele_dEtaEleCl_propOut_barrel = new TH1F( "h_ele_dEtaEleCl_propOut_barrel", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost, barrel",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut_barrel->Sumw2();
  h_ele_dEtaEleCl_propOut_eg_barrel = new TH1F( "h_ele_dEtaEleCl_propOut_eg_barrel", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost, ecal driven, barrel",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut_eg_barrel->Sumw2();
  h_ele_dEtaEleCl_propOut_endcaps = new TH1F( "h_ele_dEtaEleCl_propOut_endcaps", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost, endcaps",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut_endcaps->Sumw2();
  h_ele_dEtaEleCl_propOut_eg_endcaps = new TH1F( "h_ele_dEtaEleCl_propOut_eg_endcaps", "ele #eta_{EleCl} - #eta_{tr}, prop from outermost, ecal driven, endcaps",   nbindetamatch,detamatchmin,detamatchmax);
  h_ele_dEtaEleCl_propOut_eg_endcaps->Sumw2();
  h_ele_dEtaEleClVsEta_propOut = new TH2F( "h_ele_dEtaEleClVsEta_propOut", "ele #eta_{EleCl} - #eta_{tr} vs eta, prop from out", nbineta2D,etamin,etamax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaEleClVsPhi_propOut = new TH2F( "h_ele_dEtaEleClVsPhi_propOut", "ele #eta_{EleCl} - #eta_{tr} vs phi, prop from out", nbinphi2D,phimin,phimax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dEtaEleClVsPt_propOut = new TH2F( "h_ele_dEtaScVsPt_propOut", "ele #eta_{EleCl} - #eta_{tr} vs pt, prop from out", nbinpt2D,0.,ptmax,nbindetamatch2D,detamatchmin,detamatchmax);
  h_ele_dPhiEleCl_propOut = new TH1F( "h_ele_dPhiEleCl_propOut", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut->Sumw2();
  h_ele_dPhiEleCl_propOut_eg = new TH1F( "h_ele_dPhiEleCl_propOut_eg", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost, ecal driven",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut_eg->Sumw2();
  h_ele_dPhiEleCl_propOut_barrel = new TH1F( "h_ele_dPhiEleCl_propOut_barrel", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost, barrel",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut_barrel->Sumw2();
  h_ele_dPhiEleCl_propOut_eg_barrel = new TH1F( "h_ele_dPhiEleCl_propOut_eg_barrel", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost, ecal driven, barrel",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut_eg_barrel->Sumw2();
  h_ele_dPhiEleCl_propOut_endcaps = new TH1F( "h_ele_dPhiEleCl_propOut_endcaps", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost, endcaps",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut_endcaps->Sumw2();
  h_ele_dPhiEleCl_propOut_eg_endcaps = new TH1F( "h_ele_dPhiEleCl_propOut_eg_endcaps", "ele #phi_{EleCl} - #phi_{tr}, prop from outermost, ecal driven, endcaps",   nbindphimatch,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleCl_propOut_eg_endcaps->Sumw2();
  h_ele_dPhiEleClVsEta_propOut = new TH2F( "h_ele_dPhiEleClVsEta_propOut", "ele #phi_{EleCl} - #phi_{tr} vs eta, prop from out", nbineta2D,etamin,etamax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleClVsPhi_propOut = new TH2F( "h_ele_dPhiEleClVsPhi_propOut", "ele #phi_{EleCl} - #phi_{tr} vs phi, prop from out", nbinphi2D,phimin,phimax,nbindphimatch2D,dphimatchmin,dphimatchmax);
  h_ele_dPhiEleClVsPt_propOut = new TH2F( "h_ele_dPhiSEleClsPt_propOut", "ele #phi_{EleCl} - #phi_{tr} vs pt, prop from out", nbinpt2D,0.,ptmax,nbindphimatch2D,dphimatchmin,dphimatchmax);

  h_ele_HoE = new TH1F("h_ele_HoE", "ele hadronic energy / em energy", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE->Sumw2();
  h_ele_HoE_eg = new TH1F("h_ele_HoE_eg", "ele hadronic energy / em energy, ecal driven", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_eg->Sumw2();
  h_ele_HoE_barrel = new TH1F("h_ele_HoE_barrel", "ele hadronic energy / em energy, barrel", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_barrel->Sumw2();
  h_ele_HoE_eg_barrel = new TH1F("h_ele_HoE_eg_barrel", "ele hadronic energy / em energy, ecal driven, barrel", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_eg_barrel->Sumw2();
  h_ele_HoE_endcaps = new TH1F("h_ele_HoE_endcaps", "ele hadronic energy / em energy, endcaps", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_endcaps->Sumw2();
  h_ele_HoE_eg_endcaps = new TH1F("h_ele_HoE_eg_endcaps", "ele hadronic energy / em energy, ecal driven, endcaps", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_eg_endcaps->Sumw2();
  h_ele_HoE_fiducial = new TH1F("h_ele_HoE_fiducial", "ele hadronic energy / em energy, fiducial region", nbinhoe, hoemin, hoemax) ;
  h_ele_HoE_fiducial->Sumw2();
  h_ele_HoEVsEta = new TH2F("h_ele_HoEVsEta", "ele hadronic energy / em energy vs eta", nbineta,etamin,etamax,nbinhoe, hoemin, hoemax) ;
  h_ele_HoEVsPhi = new TH2F("h_ele_HoEVsPhi", "ele hadronic energy / em energy vs phi", nbinphi2D,phimin,phimax,nbinhoe, hoemin, hoemax) ;
  h_ele_HoEVsE = new TH2F("h_ele_HoEVsE", "ele hadronic energy / em energy vs E", nbinp, 0.,300.,nbinhoe, hoemin, hoemax) ;

  h_ele_seed_dphi2_ = new TH1F("h_ele_seedDphi2", "ele seed dphi 2nd layer", 50,-0.003,+0.003) ;
  h_ele_seed_dphi2_->Sumw2();
  h_ele_seed_dphi2VsEta_ = new TH2F("h_ele_seedDphi2VsEta", "ele seed dphi 2nd layer vs eta", nbineta2D,etamin,etamax,50,-0.003,+0.003) ;
  h_ele_seed_dphi2VsPt_ = new TH2F("h_ele_seedDphi2VsPt", "ele seed dphi 2nd layer vs pt", nbinpt2D,0.,ptmax,50,-0.003,+0.003) ;
  h_ele_seed_drz2_ = new TH1F("h_ele_seedDrz2", "ele seed dr (dz) 2nd layer", 50,-0.03,+0.03) ;
  h_ele_seed_drz2_->Sumw2();
  h_ele_seed_drz2VsEta_ = new TH2F("h_ele_seedDrz2VsEta", "ele seed dr/dz 2nd layer vs eta", nbineta2D,etamin,etamax,50,-0.03,+0.03) ;
  h_ele_seed_drz2VsPt_ = new TH2F("h_ele_seedDrz2VsPt", "ele seed dr/dz 2nd layer vs pt", nbinpt2D,0.,ptmax,50,-0.03,+0.03) ;
  h_ele_seed_subdet2_ = new TH1F("h_ele_seedSubdet2", "ele seed subdet 2nd layer", 10,0.,10.) ;
  h_ele_seed_subdet2_->Sumw2();

  // classes
  h_ele_classes = new TH1F( "h_ele_classes", "ele new classes",      20,0.0,20.);
  h_ele_classes->Sumw2();
  h_ele_old_classes = new TH1F( "h_ele_old_classes", "ele old classes",      20,0.0,20.);
  h_ele_old_classes->Sumw2();
  h_ele_eta = new TH1F( "h_ele_eta", "ele electron eta",  nbineta/2,0.0,etamax);
  h_ele_eta->Sumw2();
  h_ele_eta_golden = new TH1F( "h_ele_eta_golden", "ele electron eta golden",  nbineta/2,0.0,etamax);
  h_ele_eta_golden->Sumw2();
  h_ele_eta_bbrem = new TH1F( "h_ele_eta_bbrem", "ele electron eta bbrem",  nbineta/2,0.0,etamax);
  h_ele_eta_bbrem->Sumw2();
  h_ele_eta_narrow = new TH1F( "h_ele_eta_narrow", "ele electron eta narrow",  nbineta/2,0.0,etamax);
  h_ele_eta_narrow->Sumw2();
  h_ele_eta_shower = new TH1F( "h_ele_eta_show", "ele electron eta showering",  nbineta/2,0.0,etamax);
  h_ele_eta_shower->Sumw2();
  h_ele_PinVsPoutGolden_mode = new TH2F( "h_ele_PinVsPoutGolden_mode",      "ele track inner p vs outer p vs eta, golden, mode of GSF components" ,nbinp2D,0.,pmax,50,0.,pmax);
  h_ele_PinVsPoutShowering_mode = new TH2F( "h_ele_PinVsPoutShowering_mode",      "ele track inner p vs outer p vs eta, showering, mode of GSF components" ,nbinp2D,0.,pmax,50,0.,pmax);
  h_ele_PinVsPoutGolden_mean = new TH2F( "h_ele_PinVsPoutGolden_mean",      "ele track inner p vs outer p vs eta, golden, mean of GSF components" ,nbinp2D,0.,pmax,50,0.,pmax);
  h_ele_PinVsPoutShowering_mean = new TH2F( "h_ele_PinVsPoutShowering_mean",      "ele track inner p vs outer p vs eta, showering, mean of GSF components" ,nbinp2D,0.,pmax,50,0.,pmax);
  h_ele_PtinVsPtoutGolden_mode = new TH2F( "h_ele_PtinVsPtoutGolden_mode",      "ele track inner pt vs outer pt vs eta, golden, mode of GSF components" ,nbinpt2D,0.,ptmax,50,0.,ptmax);
  h_ele_PtinVsPtoutShowering_mode = new TH2F( "h_ele_PtinVsPtoutShowering_mode",      "ele track inner pt vs outer pt vs eta, showering, mode of GSF components" ,nbinpt2D,0.,ptmax,50,0.,ptmax);
  h_ele_PtinVsPtoutGolden_mean = new TH2F( "h_ele_PtinVsPtoutGolden_mean",      "ele track inner pt vs outer pt vs eta, golden, mean of GSF components" ,nbinpt2D,0.,ptmax,50,0.,ptmax);
  h_ele_PtinVsPtoutShowering_mean = new TH2F( "h_ele_PtinVsPtoutShowering_mean",      "ele track inner pt vs outer pt vs eta, showering, mean of GSF components" ,nbinpt2D,0.,ptmax,50,0.,ptmax);
  histSclEoEtrueGolden_barrel = new TH1F("h_scl_EoEtrue_golden_barrel","ele supercluster energy / gen energy, golden, barrel",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrueGolden_barrel->Sumw2();
  histSclEoEtrueGolden_endcaps = new TH1F("h_scl_EoEtrue_golden_endcaps","ele supercluster energy / gen energy, golden, endcaps",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrueGolden_endcaps->Sumw2();
  histSclEoEtrueShowering_barrel = new TH1F("h_scl_EoEtrue_showering_barrel","ele supercluster energy / gen energy, showering, barrel",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrueShowering_barrel->Sumw2();
  histSclEoEtrueShowering_endcaps = new TH1F("h_scl_EoEtrue_showering_endcaps","ele supercluster energy / gen energy, showering, endcaps",nbinpoptrue,poptruemin,poptruemax);
  histSclEoEtrueShowering_endcaps->Sumw2();

  // isolation
  h_ele_tkSumPt_dr03 = new TH1F("h_ele_tkSumPt_dr03","tk isolation sum, dR=0.3",100,0.0,20.);
  h_ele_tkSumPt_dr03->Sumw2();
  h_ele_ecalRecHitSumEt_dr03= new TH1F("h_ele_ecalRecHitSumEt_dr03","ecal isolation sum, dR=0.3",100,0.0,20.);
  h_ele_ecalRecHitSumEt_dr03->Sumw2();
  h_ele_hcalDepth1TowerSumEt_dr03= new TH1F("h_ele_hcalDepth1TowerSumEt_dr03","hcal depth1 isolation sum, dR=0.3",100,0.0,20.);
  h_ele_hcalDepth1TowerSumEt_dr03->Sumw2();
  h_ele_hcalDepth2TowerSumEt_dr03= new TH1F("h_ele_hcalDepth2TowerSumEt_dr03","hcal depth2 isolation sum, dR=0.3",100,0.0,20.);
  h_ele_hcalDepth2TowerSumEt_dr03->Sumw2();
  h_ele_tkSumPt_dr04= new TH1F("h_ele_tkSumPt_dr04","tk isolation sum, dR=0.4",100,0.0,20.);
  h_ele_tkSumPt_dr04->Sumw2();
  h_ele_ecalRecHitSumEt_dr04= new TH1F("h_ele_ecalRecHitSumEt_dr04","ecal isolation sum, dR=0.4",100,0.0,20.);
  h_ele_ecalRecHitSumEt_dr04->Sumw2();
  h_ele_hcalDepth1TowerSumEt_dr04= new TH1F("h_ele_hcalDepth1TowerSumEt_dr04","hcal depth1 isolation sum, dR=0.4",100,0.0,20.);
  h_ele_hcalDepth1TowerSumEt_dr04->Sumw2();
  h_ele_hcalDepth2TowerSumEt_dr04= new TH1F("h_ele_hcalDepth2TowerSumEt_dr04","hcal depth2 isolation sum, dR=0.4",100,0.0,20.);
  h_ele_hcalDepth2TowerSumEt_dr04->Sumw2();

  // fbrem
  h_ele_fbrem = new TH1F( "h_ele_fbrem","ele brem fraction, mode of GSF components",100,0.,1.);
  h_ele_fbrem->Sumw2();
  h_ele_fbrem_eg = new TH1F( "h_ele_fbrem_eg","ele brem fraction, mode of GSF components, ecal driven",100,0.,1.);
  h_ele_fbrem_eg->Sumw2();
  h_ele_fbremVsEta_mode = new TProfile( "h_ele_fbremvsEtamode","mean ele brem fraction vs eta, mode of GSF components",nbineta2D,etamin,etamax,0.,1.);
  h_ele_fbremVsEta_mean = new TProfile( "h_ele_fbremvsEtamean","mean ele brem fraction vs eta, mean of GSF components",nbineta2D,etamin,etamax,0.,1.);

  // e/g et pflow electrons
  h_ele_mva = new TH1F( "h_ele_mva","ele identification mva",100,-1.,1.);
  h_ele_mva->Sumw2();
  h_ele_mva_eg = new TH1F( "h_ele_mva_eg","ele identification mva, ecal driven",100,-1.,1.);
  h_ele_mva_eg->Sumw2();
  h_ele_provenance = new TH1F( "h_ele_provenance","ele provenance",5,-2.,3.);
  h_ele_provenance->Sumw2();

  // histos titles
  h_mcNum              -> GetXaxis()-> SetTitle("N_{gen}");
  h_mcNum              -> GetYaxis()-> SetTitle("Events");
  h_eleNum             -> GetXaxis()-> SetTitle("# gen ele");
  h_eleNum             -> GetYaxis()-> SetTitle("Events");
  h_gamNum             -> GetXaxis()-> SetTitle("N_{gen #gamma}");
  h_gamNum             -> GetYaxis()-> SetTitle("Events");
  h_simEta             -> GetXaxis()-> SetTitle("#eta");
  h_simEta             -> GetYaxis()-> SetTitle("Events");
  h_simP               -> GetXaxis()-> SetTitle("p (GeV/c)");
  h_simP               -> GetYaxis()-> SetTitle("Events");
  h_ele_foundHits      -> GetXaxis()-> SetTitle("N_{hits}");
  h_ele_foundHits      -> GetYaxis()-> SetTitle("Events");
  h_ele_foundHits_barrel      -> GetXaxis()-> SetTitle("N_{hits}");
  h_ele_foundHits_barrel      -> GetYaxis()-> SetTitle("Events");
  h_ele_foundHits_endcaps      -> GetXaxis()-> SetTitle("N_{hits}");
  h_ele_foundHits_endcaps      -> GetYaxis()-> SetTitle("Events");
  h_ele_ambiguousTracks      -> GetXaxis()-> SetTitle("N_{ambiguous tracks}");
  h_ele_ambiguousTracks      -> GetYaxis()-> SetTitle("Events");
  h_ele_lostHits       -> GetXaxis()-> SetTitle("N_{lost hits}");
  h_ele_lostHits       -> GetYaxis()-> SetTitle("Events");
  h_ele_lostHits_barrel       -> GetXaxis()-> SetTitle("N_{lost hits}");
  h_ele_lostHits_barrel       -> GetYaxis()-> SetTitle("Events");
  h_ele_lostHits_endcaps       -> GetXaxis()-> SetTitle("N_{lost hits}");
  h_ele_lostHits_endcaps       -> GetYaxis()-> SetTitle("Events");
  h_ele_chi2           -> GetXaxis()-> SetTitle("#Chi^{2}");
  h_ele_chi2           -> GetYaxis()-> SetTitle("Events");
  h_ele_chi2_barrel           -> GetXaxis()-> SetTitle("#Chi^{2}");
  h_ele_chi2_barrel           -> GetYaxis()-> SetTitle("Events");
  h_ele_chi2_endcaps           -> GetXaxis()-> SetTitle("#Chi^{2}");
  h_ele_chi2_endcaps           -> GetYaxis()-> SetTitle("Events");
  h_ele_charge         -> GetXaxis()-> SetTitle("charge");
  h_ele_charge         -> GetYaxis()-> SetTitle("Events");
  h_ele_vertexP        -> GetXaxis()-> SetTitle("p_{vertex} (GeV/c)");
  h_ele_vertexP        -> GetYaxis()-> SetTitle("Events");
  h_ele_vertexPt       -> GetXaxis()-> SetTitle("p_{T vertex} (GeV/c)");
  h_ele_vertexPt       -> GetYaxis()-> SetTitle("Events");
  h_ele_Et       -> GetXaxis()-> SetTitle("E_{T} (GeV)");
  h_ele_Et       -> GetYaxis()-> SetTitle("Events");
  h_ele_Et_all       -> GetXaxis()-> SetTitle("E_{T} (GeV)");
  h_ele_Et_all       -> GetYaxis()-> SetTitle("Events");
  h_ele_vertexEta      -> GetXaxis()-> SetTitle("#eta");
  h_ele_vertexEta      -> GetYaxis()-> SetTitle("Events");
  h_ele_vertexPhi      -> GetXaxis()-> SetTitle("#phi (rad)");
  h_ele_vertexPhi      -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_barrel        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_barrel        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_endcaps        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_endcaps        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_golden_barrel        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_golden_barrel        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_showering_barrel        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_showering_barrel        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_golden_endcaps        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_golden_endcaps        -> GetYaxis()-> SetTitle("Events");
  h_ele_PoPtrue_showering_endcaps        -> GetXaxis()-> SetTitle("P/P_{gen}");
  h_ele_PoPtrue_showering_endcaps        -> GetYaxis()-> SetTitle("Events");
  h_ele_PtoPttrue        -> GetXaxis()-> SetTitle("P_{T}/P_{T}^{gen}");
  h_ele_PtoPttrue        -> GetYaxis()-> SetTitle("Events");
  h_ele_PtoPttrue_barrel        -> GetXaxis()-> SetTitle("P_{T}/P_{T}^{gen}");
  h_ele_PtoPttrue_barrel        -> GetYaxis()-> SetTitle("Events");
  h_ele_PtoPttrue_endcaps        -> GetXaxis()-> SetTitle("P_{T}/P_{T}^{gen}");
  h_ele_PtoPttrue_endcaps        -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_barrel -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_barrel -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps_golden -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps_golden -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps_bigbrem -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps_bigbrem -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps_narrow -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps_narrow -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps_showering -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps_showering -> GetYaxis()-> SetTitle("Events") ;
  histSeedEoEtrue_endcaps_golden -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSeedEoEtrue_endcaps_golden -> GetYaxis()-> SetTitle("Events") ;
  histSeedEoEtrue_endcaps_bigbrem -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSeedEoEtrue_endcaps_bigbrem -> GetYaxis()-> SetTitle("Events") ;
  histSeedEoEtrue_endcaps_narrow -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSeedEoEtrue_endcaps_narrow -> GetYaxis()-> SetTitle("Events") ;
  histSeedEoEtrue_endcaps_showering -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSeedEoEtrue_endcaps_showering -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrueGolden_barrel -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrueGolden_barrel -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrueShowering_barrel -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrueShowering_barrel -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrueGolden_endcaps -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrueGolden_endcaps -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrueShowering_endcaps -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrueShowering_endcaps -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_barrel_etagap -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_barrel_etagap -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_barrel_phigap -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_barrel_phigap -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_ebeegap -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_ebeegap -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_endcaps_deegap -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_endcaps_deegap -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_endcaps_ringgap -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_endcaps_ringgap -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_barrel_new -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_barrel_new -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_endcaps_new -> GetXaxis()-> SetTitle("E/E_{gen}") ;
  histSclEoEtrue_endcaps_new -> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtrue_barrel_etagap_new -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_barrel_etagap_new -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_barrel_phigap_new -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_barrel_phigap_new -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_ebeegap_new -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_ebeegap_new -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_endcaps_deegap_new -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_endcaps_deegap_new -> GetYaxis()-> SetTitle("Events");
  histSclEoEtrue_endcaps_ringgap_new -> GetXaxis()-> SetTitle("E/E_{gen}");
  histSclEoEtrue_endcaps_ringgap_new -> GetYaxis()-> SetTitle("Events");
  histSclSigEtaEta_-> GetXaxis()-> SetTitle("#sigma_{#eta #eta}") ;
  histSclSigEtaEta_-> GetYaxis()-> SetTitle("Events") ;
  histSclSigEtaEta_barrel_-> GetXaxis()-> SetTitle("#sigma_{#eta #eta}") ;
  histSclSigEtaEta_barrel_-> GetYaxis()-> SetTitle("Events") ;
  histSclSigEtaEta_endcaps_-> GetXaxis()-> SetTitle("#sigma_{#eta #eta}") ;
  histSclSigEtaEta_endcaps_-> GetYaxis()-> SetTitle("Events") ;
  histSclSigIEtaIEta_-> GetXaxis()-> SetTitle("#sigma_{i#eta i#eta}") ;
  histSclSigIEtaIEta_-> GetYaxis()-> SetTitle("Events") ;
  histSclSigIEtaIEta_barrel_-> GetXaxis()-> SetTitle("#sigma_{i#eta i#eta}") ;
  histSclSigIEtaIEta_barrel_-> GetYaxis()-> SetTitle("Events") ;
  histSclSigIEtaIEta_endcaps_-> GetXaxis()-> SetTitle("#sigma_{i#eta i#eta}") ;
  histSclSigIEtaIEta_endcaps_-> GetYaxis()-> SetTitle("Events") ;
  histSclE1x5_-> GetXaxis()-> SetTitle("E1x5 (GeV)") ;
  histSclE1x5_-> GetYaxis()-> SetTitle("Events") ;
  histSclE1x5_barrel_-> GetXaxis()-> SetTitle("E1x5 (GeV)") ;
  histSclE1x5_barrel_-> GetYaxis()-> SetTitle("Events") ;
  histSclE1x5_endcaps_-> GetXaxis()-> SetTitle("E1x5 (GeV)") ;
  histSclE1x5_endcaps_-> GetYaxis()-> SetTitle("Events") ;
  histSclE2x5max_-> GetXaxis()-> SetTitle("E2x5 (GeV)") ;
  histSclE2x5max_-> GetYaxis()-> SetTitle("Events") ;
  histSclE2x5max_barrel_-> GetXaxis()-> SetTitle("E2x5 (GeV)") ;
  histSclE2x5max_barrel_-> GetYaxis()-> SetTitle("Events") ;
  histSclE2x5max_endcaps_-> GetXaxis()-> SetTitle("E2x5 (GeV)") ;
  histSclE2x5max_endcaps_-> GetYaxis()-> SetTitle("Events") ;
  histSclE5x5_-> GetXaxis()-> SetTitle("E5x5 (GeV)") ;
  histSclE5x5_-> GetYaxis()-> SetTitle("Events") ;
  histSclE5x5_barrel_-> GetXaxis()-> SetTitle("E5x5 (GeV)") ;
  histSclE5x5_barrel_-> GetYaxis()-> SetTitle("Events") ;
  histSclE5x5_endcaps_-> GetXaxis()-> SetTitle("E5x5 (GeV)") ;
  histSclE5x5_endcaps_-> GetYaxis()-> SetTitle("Events") ;
  histSclEoEtruePfVsEg->GetXaxis()->SetTitle("E/E_{gen} (e/g)") ;
  histSclEoEtruePfVsEg->GetYaxis()->SetTitle("E/E_{gen} (pflow)") ;
  h_ele_ChargeMnChargeTrue   -> GetXaxis()-> SetTitle("q_{rec} - q_{gen}");
  h_ele_ChargeMnChargeTrue   -> GetYaxis()-> SetTitle("Events");
  h_ele_EtaMnEtaTrue   -> GetXaxis()-> SetTitle("#eta_{rec} - #eta_{gen}");
  h_ele_EtaMnEtaTrue   -> GetYaxis()-> SetTitle("Events");
  h_ele_EtaMnEtaTrue_barrel   -> GetXaxis()-> SetTitle("#eta_{rec} - #eta_{gen}");
  h_ele_EtaMnEtaTrue_barrel   -> GetYaxis()-> SetTitle("Events");
  h_ele_EtaMnEtaTrue_endcaps   -> GetXaxis()-> SetTitle("#eta_{rec} - #eta_{gen}");
  h_ele_EtaMnEtaTrue_endcaps   -> GetYaxis()-> SetTitle("Events");
  h_ele_PhiMnPhiTrue   -> GetXaxis()-> SetTitle("#phi_{rec} - #phi_{gen} (rad)");
  h_ele_PhiMnPhiTrue   -> GetYaxis()-> SetTitle("Events");
  h_ele_PhiMnPhiTrue_barrel   -> GetXaxis()-> SetTitle("#phi_{rec} - #phi_{gen} (rad)");
  h_ele_PhiMnPhiTrue_barrel   -> GetYaxis()-> SetTitle("Events");
  h_ele_PhiMnPhiTrue_endcaps   -> GetXaxis()-> SetTitle("#phi_{rec} - #phi_{gen} (rad)");
  h_ele_PhiMnPhiTrue_endcaps   -> GetYaxis()-> SetTitle("Events");
  h_ele_PinMnPout      -> GetXaxis()-> SetTitle("P_{vertex} - P_{out} (GeV/c)");
  h_ele_PinMnPout      -> GetYaxis()-> SetTitle("Events");
  h_ele_PinMnPout_mode      -> GetXaxis()-> SetTitle("P_{vertex} - P_{out}, mode of GSF components (GeV/c)");
  h_ele_PinMnPout_mode      -> GetYaxis()-> SetTitle("Events");
  h_ele_outerP         -> GetXaxis()-> SetTitle("P_{out} (GeV/c)");
  h_ele_outerP         -> GetYaxis()-> SetTitle("Events");
  h_ele_outerP_mode         -> GetXaxis()-> SetTitle("P_{out} (GeV/c)");
  h_ele_outerP_mode         -> GetYaxis()-> SetTitle("Events");
  h_ele_outerPt        -> GetXaxis()-> SetTitle("P_{T out} (GeV/c)");
  h_ele_outerPt        -> GetYaxis()-> SetTitle("Events");
  h_ele_outerPt_mode        -> GetXaxis()-> SetTitle("P_{T out} (GeV/c)");
  h_ele_outerPt_mode        -> GetYaxis()-> SetTitle("Events");
  h_ele_EoP            -> GetXaxis()-> SetTitle("E/P_{vertex}");
  h_ele_EoP            -> GetYaxis()-> SetTitle("Events");
  h_ele_EseedOP            -> GetXaxis()-> SetTitle("E_{seed}/P_{vertex}");
  h_ele_EseedOP            -> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout         -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout         -> GetYaxis()-> SetTitle("Events");
  h_ele_EeleOPout         -> GetXaxis()-> SetTitle("E_{ele}/P_{out}");
  h_ele_EeleOPout         -> GetYaxis()-> SetTitle("Events");
  h_ele_EoP_barrel            -> GetXaxis()-> SetTitle("E/P_{vertex}");
  h_ele_EoP_barrel            -> GetYaxis()-> SetTitle("Events");
  h_ele_EseedOP_barrel            -> GetXaxis()-> SetTitle("E_{seed}/P_{vertex}");
  h_ele_EseedOP_barrel            -> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout_barrel         -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout_barrel         -> GetYaxis()-> SetTitle("Events");
  h_ele_EeleOPout_barrel         -> GetXaxis()-> SetTitle("E_{ele}/P_{out}");
  h_ele_EeleOPout_barrel         -> GetYaxis()-> SetTitle("Events");
  h_ele_EoP_endcaps            -> GetXaxis()-> SetTitle("E/P_{vertex}");
  h_ele_EoP_endcaps            -> GetYaxis()-> SetTitle("Events");
  h_ele_EseedOP_endcaps            -> GetXaxis()-> SetTitle("E_{seed}/P_{vertex}");
  h_ele_EseedOP_endcaps            -> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout_endcaps         -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout_endcaps         -> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout_endcaps_golden         -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout_endcaps_golden         -> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout_endcaps_showering         -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout_endcaps_showering       -> GetYaxis()-> SetTitle("Events");
  h_ele_EeleOPout_endcaps         -> GetXaxis()-> SetTitle("E_{ele}/P_{out}");
  h_ele_EeleOPout_endcaps         -> GetYaxis()-> SetTitle("Events");
  h_ele_vertexX-> GetXaxis()-> SetTitle("x (cm)");
  h_ele_vertexX-> GetYaxis()-> SetTitle("Events");
  h_ele_vertexY-> GetXaxis()-> SetTitle("y (cm)");
  h_ele_vertexY-> GetYaxis()-> SetTitle("Events");
  h_ele_vertexZ-> GetXaxis()-> SetTitle("z (cm)");
  h_ele_vertexZ-> GetYaxis()-> SetTitle("Events");
  h_ele_vertexTIP-> GetXaxis()-> SetTitle("TIP (cm)");
  h_ele_vertexTIP-> GetYaxis()-> SetTitle("Events");
  h_ele_TIP_all-> GetXaxis()-> SetTitle("r_{T} (cm)");
  h_ele_TIP_all-> GetYaxis()-> SetTitle("Events");
  h_ele_vertexTIPVsEta-> GetYaxis()-> SetTitle("TIP (cm)");
  h_ele_vertexTIPVsEta-> GetXaxis()-> SetTitle("#eta");
  h_ele_vertexTIPVsPhi-> GetYaxis()-> SetTitle("TIP (cm)");
  h_ele_vertexTIPVsPhi-> GetXaxis()-> SetTitle("#phi (rad)");
  h_ele_vertexTIPVsPt-> GetYaxis()-> SetTitle("TIP (cm)");
  h_ele_vertexTIPVsPt-> GetXaxis()-> SetTitle("p_{T} (GeV/c)");
  h_ele_dEtaSc_propVtx-> GetXaxis()-> SetTitle("#eta_{sc} - #eta_{tr}");
  h_ele_dEtaSc_propVtx-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaCl_propOut-> GetXaxis()-> SetTitle("#eta_{seedcl} - #eta_{tr}");
  h_ele_dEtaCl_propOut-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaEleCl_propOut-> GetXaxis()-> SetTitle("#eta_{elecl} - #eta_{tr}");
  h_ele_dEtaEleCl_propOut-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiSc_propVtx-> GetXaxis()-> SetTitle("#phi_{sc} - #phi_{tr} (rad)");
  h_ele_dPhiSc_propVtx-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiCl_propOut-> GetXaxis()-> SetTitle("#phi_{seedcl} - #phi_{tr} (rad)");
  h_ele_dPhiCl_propOut-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiEleCl_propOut-> GetXaxis()-> SetTitle("#phi_{elecl} - #phi_{tr} (rad)");
  h_ele_dPhiEleCl_propOut-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaSc_propVtx_barrel-> GetXaxis()-> SetTitle("#eta_{sc} - #eta_{tr}");
  h_ele_dEtaSc_propVtx_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaCl_propOut_barrel-> GetXaxis()-> SetTitle("#eta_{seedcl} - #eta_{tr}");
  h_ele_dEtaCl_propOut_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaEleCl_propOut_barrel-> GetXaxis()-> SetTitle("#eta_{elecl} - #eta_{tr}");
  h_ele_dEtaEleCl_propOut_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiSc_propVtx_barrel-> GetXaxis()-> SetTitle("#phi_{sc} - #phi_{tr} (rad)");
  h_ele_dPhiSc_propVtx_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiCl_propOut_barrel-> GetXaxis()-> SetTitle("#phi_{seedcl} - #phi_{tr} (rad)");
  h_ele_dPhiCl_propOut_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiEleCl_propOut_barrel-> GetXaxis()-> SetTitle("#phi_{elecl} - #phi_{tr} (rad)");
  h_ele_dPhiEleCl_propOut_barrel-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaSc_propVtx_endcaps-> GetXaxis()-> SetTitle("#eta_{sc} - #eta_{tr}");
  h_ele_dEtaSc_propVtx_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaCl_propOut_endcaps-> GetXaxis()-> SetTitle("#eta_{seedcl} - #eta_{tr}");
  h_ele_dEtaCl_propOut_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaEleCl_propOut_endcaps-> GetXaxis()-> SetTitle("#eta_{elecl} - #eta_{tr}");
  h_ele_dEtaEleCl_propOut_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiSc_propVtx_endcaps-> GetXaxis()-> SetTitle("#phi_{sc} - #phi_{tr} (rad)");
  h_ele_dPhiSc_propVtx_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiCl_propOut_endcaps-> GetXaxis()-> SetTitle("#phi_{seedcl} - #phi_{tr} (rad)");
  h_ele_dPhiCl_propOut_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiEleCl_propOut_endcaps-> GetXaxis()-> SetTitle("#phi_{elecl} - #phi_{tr} (rad)");
  h_ele_dPhiEleCl_propOut_endcaps-> GetYaxis()-> SetTitle("Events");
  h_ele_HoE-> GetXaxis()-> SetTitle("H/E") ;
  h_ele_HoE-> GetYaxis()-> SetTitle("Events") ;
  h_ele_HoE_barrel-> GetXaxis()-> SetTitle("H/E") ;
  h_ele_HoE_barrel-> GetYaxis()-> SetTitle("Events") ;
  h_ele_HoE_endcaps-> GetXaxis()-> SetTitle("H/E") ;
  h_ele_HoE_endcaps-> GetYaxis()-> SetTitle("Events") ;
  h_ele_HoE_fiducial-> GetXaxis()-> SetTitle("H/E") ;
  h_ele_HoE_fiducial-> GetYaxis()-> SetTitle("Events") ;
  h_ele_fbrem-> GetXaxis()-> SetTitle("P_{in} - P_{out} / P_{in}");
  h_ele_fbrem-> GetYaxis()-> SetTitle("Events");
  h_ele_seed_dphi2_-> GetXaxis()-> SetTitle("#phi_{hit}-#phi_{pred} (rad)") ;
  h_ele_seed_dphi2_-> GetYaxis()-> SetTitle("Events") ;
  h_ele_seed_drz2_-> GetXaxis()-> SetTitle("r(z)_{hit}-r(z)_{pred} (cm)") ;
  h_ele_seed_drz2_-> GetYaxis()-> SetTitle("Events") ;
  h_ele_seed_subdet2_-> GetXaxis()-> SetTitle("2nd hit subdet Id") ;
  h_ele_seed_subdet2_-> GetYaxis()-> SetTitle("Events") ;
  h_ele_classes-> GetXaxis()-> SetTitle("class Id") ;
  h_ele_classes-> GetYaxis()-> SetTitle("Events") ;
  h_ele_old_classes-> GetXaxis()-> SetTitle("class Id") ;
  h_ele_old_classes-> GetYaxis()-> SetTitle("Events") ;
  h_ele_EoverP_all-> GetXaxis()-> SetTitle("E/P_{vertex}");
  h_ele_EoverP_all-> GetYaxis()-> SetTitle("Events");
  h_ele_EseedOP_all-> GetXaxis()-> SetTitle("E_{seed}/P_{vertex}");
  h_ele_EseedOP_all-> GetYaxis()-> SetTitle("Events");
  h_ele_EoPout_all -> GetXaxis()-> SetTitle("E_{seed}/P_{out}");
  h_ele_EoPout_all-> GetYaxis()-> SetTitle("Events");
  h_ele_EeleOPout_all-> GetXaxis()-> SetTitle("E_{ele}/P_{out}");
  h_ele_EeleOPout_all-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaSc_propVtx_all-> GetXaxis()-> SetTitle("#eta_{sc} - #eta_{tr}");
  h_ele_dEtaSc_propVtx_all-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiSc_propVtx_all-> GetXaxis()-> SetTitle("#phi_{sc} - #phi_{tr} (rad)");
  h_ele_dPhiSc_propVtx_all-> GetYaxis()-> SetTitle("Events");
  h_ele_dEtaCl_propOut_all-> GetXaxis()-> SetTitle("#eta_{sc} - #eta_{tr}");
  h_ele_dEtaCl_propOut_all-> GetYaxis()-> SetTitle("Events");
  h_ele_dPhiCl_propOut_all-> GetXaxis()-> SetTitle("#phi_{sc} - #phi_{tr} (rad)");
  h_ele_dPhiCl_propOut_all-> GetYaxis()-> SetTitle("Events");
  h_ele_HoE_all-> GetXaxis()-> SetTitle("H/E") ;
  h_ele_HoE_all-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_all-> GetXaxis()-> SetTitle("m_{ee} (GeV/c^{2})");
  h_ele_mee_all-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_seed_all-> GetXaxis()-> SetTitle("m_{ee} (GeV/c^{2})");
  h_ele_mee_seed_all-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_sc_all-> GetXaxis()-> SetTitle("m_{ee} (GeV/c^{2})");
  h_ele_mee_sc_all-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_seed_os-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_seed_os-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_sc_os-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_sc_os-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_ebeb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_ebeb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_ebee-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_ebee-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_eeee-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_eeee-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_gg-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_gg-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_gb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_gb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_os_bb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_os_bb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_seed_os_gg-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_seed_os_gg-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_seed_os_gb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_seed_os_gb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_seed_os_bb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_seed_os_bb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_sc_os_gg-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_sc_os_gg-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_sc_os_gb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_sc_os_gb-> GetYaxis()-> SetTitle("Events");
  h_ele_mee_sc_os_bb-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_mee_sc_os_bb-> GetYaxis()-> SetTitle("Events");
  h_ele_E2mnE1vsMee_all-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_E2mnE1vsMee_all-> GetYaxis()-> SetTitle("E2 - E1 (GeV)");
  h_ele_E2mnE1vsMee_egeg_all-> GetXaxis()-> SetTitle("m_{e^{+}e^{-}} (GeV/c^{2})");
  h_ele_E2mnE1vsMee_egeg_all-> GetYaxis()-> SetTitle("E2 - E1 (GeV)");
  histNum_-> GetXaxis()-> SetTitle("N_{ele}");
  histNum_-> GetYaxis()-> SetTitle("Events");
  h_ele_fbremVsEta_mode-> GetXaxis()-> SetTitle("#eta");
  h_ele_fbremVsEta_mean-> GetXaxis()-> SetTitle("#eta");

  h_hgcal_sclusters_multiplicity = new TH1F("h_hgcal_sclusters_multiplicity","hgcal em supercluster multiplicity",50,0.,50.);
  h_hgcal_sclusters_Etsubclusters = new TH1F("h_hgcal_sclusters__Etsubclusters","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning = new TH1F("h_hgcal_sclusters__Etsubclusters_etcut_detacut_cleaning","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning_long = new TH1F("h_hgcal_sclusters__Etsubclusters_etcut_detacut_cleaning_long","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_scclusters_detadphisubclusters = new TH2F("h_hgcal_scclusters_detadphisubclusters","hgcal em subcluster deta vs dphi",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_signed = new TH2F("h_hgcal_scclusters_detadphisubclusters_signed","hgcal em subcluster deta vs signed dphi",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_pca = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca","hgcal em subcluster deta vs dphi pca",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_pca_signed = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca_signed","hgcal em subcluster deta vs signed dphi pca",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_pca_weighted = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca_weighted","hgcal em subcluster deta vs dphi pca weighted",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted","hgcal em subcluster deta vs signed dphi pca weighted",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_bremkinematic= new TH2F("h_hgcal_scclusters_bremkinematic","hgcal em subcluster relative dphi pca vs cluster energy fraction ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_weighted= new TH2F("h_hgcal_scclusters_bremkinematic_weighted","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_etcut= new TH2F("h_hgcal_scclusters_bremkinematic_etcut","hgcal em subcluster relative dphi pca vs cluster energy fraction ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_weighted_etcut= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_etcut","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_etcut_detacut_long= new TH2F("h_hgcal_scclusters_bremkinematic_etcut_detacut_long","hgcal em subcluster relative dphi pca vs cluster energy fraction ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_long= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_long","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematicnew= new TH2F("h_hgcal_scclusters_bremkinematicnew","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_weighted= new TH2F("h_hgcal_scclusters_bremkinematicnew_weighted","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_etcut= new TH2F("h_hgcal_scclusters_bremkinematicnew_etcut","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_weighted_etcut= new TH2F("h_hgcal_scclusters_bremkinematicnew_weighted_etcut","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned_long= new TH2F("h_hgcal_scclusters_bremkinematicn_etcut_detacut_cleaned_long","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned_long= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned_long","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematicnew_weighted_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematicnew_weighted_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematic_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematic_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction ",200,-0.5,0.5,200,0.,1.0);
  h_hgcal_scclusters_bremkinematic_weighted_detacut_cleaned= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_detacut_cleaned","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",200,-0.5,0.5,200,0.,1.0);

  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_pt15 = new TH1F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_pt15","hgcal em supercluster energyfraction",200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_pt15 = new TH1F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_pt15","hgcal em supercluster energyfraction",200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_golden_pt15 = new TH1F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_golden_pt15","hgcal em supercluster energyfraction",200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_supergolden_pt15 = new TH1F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_supergolden_pt15","hgcal em supercluster energyfraction",200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeta_pt15 = new TH2F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeta_pt15","hgcal em supercluster energyfraction vs eta",150,1.5,3.0,200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVStruevertices_pt15 = new TH2F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVStruevertices_pt15","hgcal em supercluster energyfraction vs truevertices",200,0.,600.,200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSfbrem_pt15 = new TH2F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSfbrem_pt15","hgcal em supercluster energyfraction vs fbrem",200,0.,1.,200,0.,1.);
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeoveretrue_pt15  = new TH2F("h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeoveretrue_pt15","hgcal em supercluster energyfraction vs eoveretrue",150,0.,1.5,200,0.,1.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long = new TH1F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long","hgcal em supercluster cellmultiplicity",200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_golden = new TH1F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_golden","hgcal em supercluster cellmultiplicity",200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_supergolden = new TH1F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_supergolden","hgcal em supercluster cellmultiplicity",200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeta = new TH2F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeta","hgcal em supercluster cellmultiplicity vs eta",150,1.5,3.0,200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVStruevertices = new TH2F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVStruevertices","hgcal em supercluster cellmultiplicity vs truevertices",80,100.,180.,200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue = new TH2F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue","hgcal em supercluster cellmultiplicity vs eoveretrue",150,0.,1.5,200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_showering = new TH2F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_showering","hgcal em supercluster cellmultiplicity vs eoveretrue showering",150,0.,1.5,200,0.,2000.);
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_golden = new TH2F("h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_golden","hgcal em supercluster cellmultiplicity vs eoveretrue golden",150,0.,1.5,200,0.,2000.);
  h_hgcal_seedcluster_cellmultiplicity_golden = new TH1F("h_hgcal_seedcluster_cellmultiplicity_golden","hgcal em cluster cellmultiplicity golden",200,0.,1000.);
  h_hgcal_seedcluster_cellenergy_golden = new TH1F("h_hgcal_seedcluster_cellmenergy_golden","hgcal em cluster cellenergy golden",2000,0.00001,0.2);
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning = new TH1F("h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning","hgcal em supercluster multiplicity",50,0.,50.);
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long = new TH1F("h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long","hgcal em supercluster multiplicity",50,0.,50.);
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long_golden = new TH1F("h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long_golden","hgcal em supercluster multiplicity",50,0.,50.);
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVSeta = new TH2F("h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVSeta","hgcal em supercluster multiplicity vs eta",150,1.5,3.0,50,0.,50.);
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVStruevertices = new TH2F("h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVStruevertices","hgcal em supercluster multiplicity vs truevertices",200,0.,600.,50,0.,50.);
  h_hgcal_sclusters_Etsubclusters_long = new TH1F("h_hgcal_sclusters__Etsubclusters_long","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning = new TH1F("h_hgcal_sclusters__Etsubclusters_etcut_detacut_cleaning","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning_long = new TH1F("h_hgcal_sclusters__Etsubclusters_etcut_detacut_cleaning_long","hgcal em supercluster subclustrs transverse energy",100,0.,10.);
  h_hgcal_scclusters_detadphisubclusters_pca_signed_long = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca_signed_long","hgcal em subcluster deta vs signed dphi pca",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted_long = new TH2F("h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted_long","hgcal em subcluster deta vs signed dphi pca weighted",120,-0.3,0.3,120,-0.3,0.3);
  h_hgcal_scclusters_bremkinematic_long= new TH2F("h_hgcal_scclusters_bremkinematic_long","hgcal em subcluster relative dphi pca vs cluster energy fraction ",120,-0.3,0.3,100,0.,0.5);
  h_hgcal_scclusters_bremkinematic_weighted_long= new TH2F("h_hgcal_scclusters_bremkinematic_weighted_long","hgcal em subcluster relative dphi pca vs cluster energy fraction weighted ",120,-0.3,0.3,100,0.,0.5);
  
  h_hgcal_sclusters_scenergy = new TH1F("h_hgcal_sclusters_scenergy_","hgcal em supercluster energy",150,0.,1.5);
  h_hgcal_sclusters_scenergy_pt15 = new TH1F("h_hgcal_sclusters_scenergy_pt15","hgcal em supercluster energy",150,0.,1.5);
  h_hgcal_sclusters_scenergy_golden_pt15 = new TH1F("h_hgcal_sclusters_scenergy_golden_pt15","hgcal em supercluster energy",150,0.,1.5);
  h_hgcal_sclusters_scenergy_showering_pt15 = new TH1F("h_hgcal_sclusters_scenergy_showering_pt15","hgcal em supercluster energy",150,0.,1.5);
  h_hgcal_sclusters_scenergy_supergolden_pt15 = new TH1F("h_hgcal_sclusters_scenergy_supergolden_pt15","hgcal em supercluster energy",150,0.,1.5);
  h_hgcal_sclusters_scenergyVSvertices = new TH2F("h_hgcal_sclusters_scenergyVSvertices_","hgcal em supercluster energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_scenergyVStruevertices = new TH2F("h_hgcal_sclusters_scenergyVStruevertices_","hgcal em supercluster energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_seedenergy_goldenVSdeltaphiElePin = new TH2F("h_hgcal_sclusters_seedenergy_goldenVSdeltaphiElePin","hgcal em seedcluster energy golden vs delphi(ele,pin)",150,-0.15,0.15,150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden = new TH1F("h_hgcal_sclusters_seedenergy_golden","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_3x3 = new TH1F("h_hgcal_sclusters_seedenergy_golden_3x3","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_5x5 = new TH1F("h_hgcal_sclusters_seedenergy_golden_5x5","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_10x10 = new TH1F("h_hgcal_sclusters_seedenergy_golden_10x10","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden = new TH1F("h_hgcal_sclusters_seedenergy_supergolden","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt15 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt15","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt30 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt30","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt1020 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt1020","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt2030 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt2030","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt3040 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt3040","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt4050 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt4050","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt50 = new TH1F("h_hgcal_sclusters_seedenergy_supergolden_5x5_pt50","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt15 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt15","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt30 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt30","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt1020 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt1020","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt2030 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt2030","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt3040 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt3040","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt4050 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt4050","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_seedenergy_golden_pt50 = new TH1F("h_hgcal_sclusters_seedenergy_golden_pt50","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt15","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt30 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt30","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt1020 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt1020","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt2030 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt2030","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt3040 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt3040","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt4050 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt4050","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt50 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt50","hgcal em seedcluster energy supergolden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_golden_pt15 = new TH1F("h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_golden_pt15","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt15","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt30 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt30","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt1020 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt1020","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt2030 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt2030","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt3040 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt3040","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt4050 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt4050","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt50 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt50","hgcal em seedcluster energy golden",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_showering_pt15 = new TH1F("h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_showering_pt15","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt15","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt30 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt30","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt1020 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt1020","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt2030 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt2030","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt3040 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt3040","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt4050 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt4050","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt50 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt50","hgcal em seedcluster energy showering",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_detacut = new TH1F("h_hgcal_sclusters_newscenergy_detacut","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut = new TH1F("h_hgcal_sclusters_newscenergy_etcut","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_dynetcut = new TH1F("h_hgcal_sclusters_newscenergy_dynetcut","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_long = new TH1F("h_hgcal_sclusters_newscenergy_long","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_long = new TH1F("h_hgcal_sclusters_newscenergy_etcut_long","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_pt15 = new TH1F("h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_pt15","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt30 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt30","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1530 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1530","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1020 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1020","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt2030 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt2030","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt3040 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt3040","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt4050 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt4050","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt50 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt50","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_pt_pt1020 = new TH1F("h_hgcal_pt_pt1020","hgcal em supercluster new energy",100,0.,100.);
  h_hgcal_pt_pt2030 = new TH1F("h_hgcal_pt_pt2030","hgcal em supercluster new energy",100,0.,100.);
  h_hgcal_pt_pt3040 = new TH1F("h_hgcal_pt_pt3040","hgcal em supercluster new energy",100,0.,100.);
  h_hgcal_pt_pt4050 = new TH1F("h_hgcal_pt_pt4050","hgcal em supercluster new energy",100,0.,100.);
  h_hgcal_pt_pt50 = new TH1F("h_hgcal_pt_pt50","hgcal em supercluster new energy",100,0.,100.);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt40inf = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt40inf","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_frachigh_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_frachigh_pt15","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fracmiddle_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fracmiddle_pt15","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fraclow_pt15 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fralow_pt15","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_long = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_long","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSeta = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSeta","hgcal em supercluster new energy vs vertices",150,1.5,3.0,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi","hgcal em supercluster new energy vs vertices",1440,-180.,180.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi_aftercut = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi_aftercut","hgcal em supercluster new energy vs vertices",1440,-180.,180.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_longVSvertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_longVSvertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_longVStruevertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_longVStruevertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSvertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSvertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVStruevertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVStruevertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVSvertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVSvertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVStruevertices = new TH2F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVStruevertices","hgcal em supercluster new energy vs vertices",200,0.,600.,150,0.,1.5);
  h_hgcal_sclusters_newscenergy_detacut_cleaning = new TH1F("h_hgcal_sclusters_newscenergy_detacut_cleaning","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to120 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to120","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu120to130 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu120to130","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu130to140 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu130to140","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to150 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to150","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu150to160 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu150to160","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu160to170 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu160to170","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta16to17 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta16to17","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta17to18 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta17to18","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta18to19 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta18to19","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta19to20 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta19to20","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta20to21 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta20to21","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta21to22 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta21to22","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta22to23 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta22to23","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta23to24 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta23to24","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta24to25 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta24to25","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta25to26 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta25to26","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta26to27 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta26to27","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta27to28 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta27to28","hgcal em supercluster new energy",150,0.,1.5);
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta28to29 = new TH1F("h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta28to29","hgcal em supercluster new energy",150,0.,1.5);

  h_eb_sclusters_scenergy_pt15_pu110to125 = new TH1F("h_eb_sclusters_scenergy_pt15_pu110to125","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_scenergy_pt15_pu125to140 = new TH1F("h_eb_sclusters_scenergy_pt15_pu125to140","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_scenergy_pt15_pu140to155 = new TH1F("h_eb_sclusters_scenergy_pt15_pu140to155","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_scenergy_pt15_pu155to170 = new TH1F("h_eb_sclusters_scenergy_pt15_pu155to170","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to125 = new TH1F("h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to125","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu125to140 = new TH1F("h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu125to140","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to155 = new TH1F("h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to155","eb em supercluster new energy",150,0.,1.5);
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu155to170 = new TH1F("h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu155to170","eb em supercluster new energy",150,0.,1.5);

  h_ele_hgcal_hoverem = new TH1F("h_ele_hgcal_hoverem","hgcal hadronic over em",100,0.0,0.5);
  h_ele_hgcal_sietaieta = new TH1F("h_ele_hgcal_sietaieta","hgcal seed cluster shower eta sigma corr",100,0.,0.025);
  
}

void
HGCALGsfElectronAnalyzer::endJob(){

  histfile_->cd();

  std::cout << "[HGCALGsfElectronAnalyzer] efficiency calculation " << std::endl;
  // efficiency vs eta
  TH1F *h_ele_etaEff = (TH1F*)h_ele_simEta_matched->Clone("h_ele_etaEff");
  h_ele_etaEff->Reset();
  h_ele_etaEff->Divide(h_ele_simEta_matched,h_simEta,1,1,"b");
  h_ele_etaEff->Print();
  h_ele_etaEff->GetXaxis()->SetTitle("#eta");
  h_ele_etaEff->GetYaxis()->SetTitle("Efficiency");

   // efficiency vs z
   TH1F *h_ele_zEff = (TH1F*)h_ele_simZ_matched->Clone("h_ele_zEff");
   h_ele_zEff->Reset();
   h_ele_zEff->Divide(h_ele_simZ_matched,h_simZ,1,1,"b");
   h_ele_zEff->Print();
   h_ele_zEff->GetXaxis()->SetTitle("z (cm)");
   h_ele_zEff->GetYaxis()->SetTitle("Efficiency");
 
   // efficiency vs |eta|
   TH1F *h_ele_absetaEff = (TH1F*)h_ele_simAbsEta_matched->Clone("h_ele_absetaEff");
   h_ele_absetaEff->Reset();
   h_ele_absetaEff->Divide(h_ele_simAbsEta_matched,h_simAbsEta,1,1,"b");
   h_ele_absetaEff->GetXaxis()->SetTitle("|#eta|");
   h_ele_absetaEff->GetYaxis()->SetTitle("Efficiency");
 
   // efficiency vs pt
   TH1F *h_ele_ptEff = (TH1F*)h_ele_simPt_matched->Clone("h_ele_ptEff");
   h_ele_ptEff->Reset();
   h_ele_ptEff->Divide(h_ele_simPt_matched,h_simPt,1,1,"b");
   h_ele_ptEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_ele_ptEff->GetYaxis()->SetTitle("Efficiency");
 
   // efficiency vs phi
   TH1F *h_ele_phiEff = (TH1F*)h_ele_simPhi_matched->Clone("h_ele_phiEff");
   h_ele_phiEff->Reset();
   h_ele_phiEff->Divide(h_ele_simPhi_matched,h_simPhi,1,1,"b");
   h_ele_phiEff->GetXaxis()->SetTitle("#phi (rad)");
   h_ele_phiEff->GetYaxis()->SetTitle("Efficiency");
 
   // efficiency vs pt eta
   TH2F *h_ele_ptEtaEff = (TH2F*)h_ele_simPtEta_matched->Clone("h_ele_ptEtaEff");
   h_ele_ptEtaEff->Reset();
   h_ele_ptEtaEff->Divide(h_ele_simPtEta_matched,h_simPtEta,1,1,"b");
   h_ele_ptEtaEff->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   h_ele_ptEtaEff->GetXaxis()->SetTitle("#eta");
 
   std::cout << "[HGCALGsfElectronAnalyzer] q-misid calculation " << std::endl;
   // misid vs eta
   TH1F *h_ele_etaQmisid = (TH1F*)h_ele_simEta_matched_qmisid->Clone("h_ele_etaQmisid");
   h_ele_etaQmisid->Reset();
   h_ele_etaQmisid->Divide(h_ele_simEta_matched_qmisid,h_simEta,1,1,"b");
   h_ele_etaQmisid->Print();
   h_ele_etaQmisid->GetXaxis()->SetTitle("#eta");
   h_ele_etaQmisid->GetYaxis()->SetTitle("q misId");
 
   // misid vs z
   TH1F *h_ele_zQmisid = (TH1F*)h_ele_simZ_matched_qmisid->Clone("h_ele_zQmisid");
   h_ele_zQmisid->Reset();
   h_ele_zQmisid->Divide(h_ele_simZ_matched_qmisid,h_simZ,1,1,"b");
   h_ele_zQmisid->Print();
   h_ele_zQmisid->GetXaxis()->SetTitle("z (cm)");
   h_ele_zQmisid->GetYaxis()->SetTitle("q misId");
 
   // misid vs |eta|
   TH1F *h_ele_absetaQmisid = (TH1F*)h_ele_simAbsEta_matched_qmisid->Clone("h_ele_absetaQmisid");
   h_ele_absetaQmisid->Reset();
   h_ele_absetaQmisid->Divide(h_ele_simAbsEta_matched_qmisid,h_simAbsEta,1,1,"b");
   h_ele_absetaQmisid->GetXaxis()->SetTitle("|#eta|");
   h_ele_absetaQmisid->GetYaxis()->SetTitle("q misId");
 
   // misid vs pt
   TH1F *h_ele_ptQmisid = (TH1F*)h_ele_simPt_matched_qmisid->Clone("h_ele_ptQmisid");
   h_ele_ptQmisid->Reset();
   h_ele_ptQmisid->Divide(h_ele_simPt_matched_qmisid,h_simPt,1,1,"b");
   h_ele_ptQmisid->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_ele_ptQmisid->GetYaxis()->SetTitle("q misId");
 
   std::cout << "[HGCALGsfElectronAnalyzer] all reco electrons " << std::endl;
   // rec/gen all electrons
   TH1F *h_ele_etaEff_all = (TH1F*)h_ele_vertexEta_all->Clone("h_ele_etaEff_all");
   h_ele_etaEff_all->Reset();
   h_ele_etaEff_all->Divide(h_ele_vertexEta_all,h_simEta,1,1,"b");
   h_ele_etaEff_all->Print();
   h_ele_etaEff_all->GetXaxis()->SetTitle("#eta");
   h_ele_etaEff_all->GetYaxis()->SetTitle("N_{rec}/N_{gen}");
   TH1F *h_ele_ptEff_all = (TH1F*)h_ele_vertexPt_all->Clone("h_ele_ptEff_all");
  h_ele_ptEff_all->Reset();
  h_ele_ptEff_all->Divide(h_ele_vertexPt_all,h_simPt,1,1,"b");
  h_ele_ptEff_all->Print();
  h_ele_ptEff_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_ele_ptEff_all->GetYaxis()->SetTitle("N_{rec}/N_{gen}");

  // classes
   TH1F *h_ele_eta_goldenFrac = (TH1F*)h_ele_eta_golden->Clone("h_ele_eta_goldenFrac");
  h_ele_eta_goldenFrac->Reset();
  h_ele_eta_goldenFrac->Divide(h_ele_eta_golden,h_ele_eta,1,1,"b");
  h_ele_eta_goldenFrac->GetXaxis()->SetTitle("|#eta|");
  h_ele_eta_goldenFrac->GetYaxis()->SetTitle("Fraction of electrons");
  h_ele_eta_goldenFrac->SetTitle("fraction of golden electrons vs eta");
   TH1F *h_ele_eta_bbremFrac = (TH1F*)h_ele_eta_bbrem->Clone("h_ele_eta_bbremFrac");
  h_ele_eta_bbremFrac->Reset();
  h_ele_eta_bbremFrac->GetXaxis()->SetTitle("|#eta|");
  h_ele_eta_bbremFrac->GetYaxis()->SetTitle("Fraction of electrons");
  h_ele_eta_bbremFrac->Divide(h_ele_eta_bbrem,h_ele_eta,1,1,"b");
  h_ele_eta_bbremFrac->SetTitle("fraction of big brem electrons vs eta");
   TH1F *h_ele_eta_narrowFrac = (TH1F*)h_ele_eta_narrow->Clone("h_ele_eta_narrowFrac");
  h_ele_eta_narrowFrac->Reset();
  h_ele_eta_narrowFrac->Divide(h_ele_eta_narrow,h_ele_eta,1,1,"b");
  h_ele_eta_narrowFrac->GetXaxis()->SetTitle("|#eta|");
  h_ele_eta_narrowFrac->GetYaxis()->SetTitle("Fraction of electrons");
  h_ele_eta_narrowFrac->SetTitle("fraction of narrow electrons vs eta");
   TH1F *h_ele_eta_showerFrac = (TH1F*)h_ele_eta_shower->Clone("h_ele_eta_showerFrac");
  h_ele_eta_showerFrac->Reset();
  h_ele_eta_showerFrac->Divide(h_ele_eta_shower,h_ele_eta,1,1,"b");
  h_ele_eta_showerFrac->GetXaxis()->SetTitle("|#eta|");
  h_ele_eta_showerFrac->GetYaxis()->SetTitle("Fraction of electrons");
  h_ele_eta_showerFrac->SetTitle("fraction of showering electrons vs eta");

  // fbrem
  TH1F *h_ele_xOverX0VsEta = new TH1F( "h_ele_xOverx0VsEta","mean X/X_0 vs eta",nbineta/2,0.0,etamax);
  for (int ibin=1;ibin<h_ele_fbremVsEta_mean->GetNbinsX()+1;ibin++) {
    double xOverX0 = 0.;
    if (h_ele_fbremVsEta_mean->GetBinContent(ibin)>0.) xOverX0 = -log(h_ele_fbremVsEta_mean->GetBinContent(ibin));
    h_ele_xOverX0VsEta->SetBinContent(ibin,xOverX0);
  }

  //profiles from 2D histos
  TProfile *p_ele_PoPtrueVsEta = h_ele_PoPtrueVsEta->ProfileX();
  p_ele_PoPtrueVsEta->SetTitle("mean ele momentum / gen momentum vs eta");
  p_ele_PoPtrueVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_PoPtrueVsEta->GetYaxis()->SetTitle("<P/P_{gen}>");
  p_ele_PoPtrueVsEta->Write();
  TProfile *p_ele_PoPtrueVsPhi = h_ele_PoPtrueVsPhi->ProfileX();
  p_ele_PoPtrueVsPhi->SetTitle("mean ele momentum / gen momentum vs phi");
  p_ele_PoPtrueVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_PoPtrueVsPhi->GetYaxis()->SetTitle("<P/P_{gen}>");
  p_ele_PoPtrueVsPhi->Write();
  TProfile *p_ele_EoEtruePfVsEg_x = histSclEoEtruePfVsEg->ProfileX();
  p_ele_EoEtruePfVsEg_x->SetTitle("mean pflow sc energy / true energy vs e/g sc energy");
  p_ele_EoEtruePfVsEg_x->GetXaxis()->SetTitle("E/E_{gen} (e/g)") ;
  p_ele_EoEtruePfVsEg_x->GetYaxis()->SetTitle("<E/E_{gen}> (pflow)") ;
  p_ele_EoEtruePfVsEg_x->Write();
  TProfile *p_ele_EoEtruePfVsEg_y = histSclEoEtruePfVsEg->ProfileY();
  p_ele_EoEtruePfVsEg_y->SetTitle("mean e/g sc energy / true energy vs pflow sc energy");
  p_ele_EoEtruePfVsEg_y->GetXaxis()->SetTitle("E/E_{gen} (pflow)") ;
  p_ele_EoEtruePfVsEg_y->GetYaxis()->SetTitle("<E/E_{gen}> (eg)") ;
  p_ele_EoEtruePfVsEg_y->Write();
  TProfile *p_ele_EtaMnEtaTrueVsEta = h_ele_EtaMnEtaTrueVsEta->ProfileX();
  p_ele_EtaMnEtaTrueVsEta->SetTitle("mean ele eta - gen eta vs eta");
  p_ele_EtaMnEtaTrueVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_EtaMnEtaTrueVsEta->GetYaxis()->SetTitle("<#eta_{rec} - #eta_{gen}>");
  p_ele_EtaMnEtaTrueVsEta->Write();
  TProfile *p_ele_EtaMnEtaTrueVsPhi = h_ele_EtaMnEtaTrueVsPhi->ProfileX();
  p_ele_EtaMnEtaTrueVsPhi->SetTitle("mean ele eta - gen eta vs phi");
  p_ele_EtaMnEtaTrueVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_EtaMnEtaTrueVsPhi->GetYaxis()->SetTitle("<#eta_{rec} - #eta_{gen}>");
  p_ele_EtaMnEtaTrueVsPhi->Write();
  TProfile *p_ele_PhiMnPhiTrueVsEta = h_ele_PhiMnPhiTrueVsEta->ProfileX();
  p_ele_PhiMnPhiTrueVsEta->SetTitle("mean ele phi - gen phi vs eta");
  p_ele_PhiMnPhiTrueVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_PhiMnPhiTrueVsEta->GetYaxis()->SetTitle("<#phi_{rec} - #phi_{gen}> (rad)");
  p_ele_PhiMnPhiTrueVsEta->Write();
  TProfile *p_ele_PhiMnPhiTrueVsPhi = h_ele_PhiMnPhiTrueVsPhi->ProfileX();
  p_ele_PhiMnPhiTrueVsPhi->SetTitle("mean ele phi - gen phi vs phi");
  p_ele_PhiMnPhiTrueVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_PhiMnPhiTrueVsPhi->Write();
  TProfile *p_ele_vertexPtVsEta = h_ele_vertexPtVsEta->ProfileX();
  p_ele_vertexPtVsEta->SetTitle("mean ele transverse momentum vs eta");
  p_ele_vertexPtVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_vertexPtVsEta->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
  p_ele_vertexPtVsEta->Write();
  TProfile *p_ele_vertexPtVsPhi = h_ele_vertexPtVsPhi->ProfileX();
  p_ele_vertexPtVsPhi->SetTitle("mean ele transverse momentum vs phi");
  p_ele_vertexPtVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_vertexPtVsPhi->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
  p_ele_vertexPtVsPhi->Write();
  TProfile *p_ele_EoPVsEta = h_ele_EoPVsEta->ProfileX();
  p_ele_EoPVsEta->SetTitle("mean ele E/p vs eta");
  p_ele_EoPVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_EoPVsEta->GetYaxis()->SetTitle("<E/P_{vertex}>");
  p_ele_EoPVsEta->Write();
  TProfile *p_ele_EoPVsPhi = h_ele_EoPVsPhi->ProfileX();
  p_ele_EoPVsPhi->SetTitle("mean ele E/p vs phi");
  p_ele_EoPVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_EoPVsPhi->GetYaxis()->SetTitle("<E/P_{vertex}>");
  p_ele_EoPVsPhi->Write();
  TProfile *p_ele_EoPoutVsEta = h_ele_EoPoutVsEta->ProfileX();
  p_ele_EoPoutVsEta->SetTitle("mean ele E/pout vs eta");
  p_ele_EoPoutVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_EoPoutVsEta->GetYaxis()->SetTitle("<E_{seed}/P_{out}>");
  p_ele_EoPoutVsEta->Write();
  TProfile *p_ele_EoPoutVsPhi = h_ele_EoPoutVsPhi->ProfileX();
  p_ele_EoPoutVsPhi->SetTitle("mean ele E/pout vs phi");
  p_ele_EoPoutVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_EoPoutVsPhi->GetYaxis()->SetTitle("<E_{seed}/P_{out}>");
  p_ele_EoPoutVsPhi->Write();
  TProfile *p_ele_EeleOPoutVsEta = h_ele_EeleOPoutVsEta->ProfileX();
  p_ele_EeleOPoutVsEta->SetTitle("mean ele Eele/pout vs eta");
  p_ele_EeleOPoutVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_EeleOPoutVsEta->GetYaxis()->SetTitle("<E_{ele}/P_{out}>");
  p_ele_EeleOPoutVsEta->Write();
  TProfile *p_ele_EeleOPoutVsPhi = h_ele_EeleOPoutVsPhi->ProfileX();
  p_ele_EeleOPoutVsPhi->SetTitle("mean ele Eele/pout vs phi");
  p_ele_EeleOPoutVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_EeleOPoutVsPhi->GetYaxis()->SetTitle("<E_{ele}/P_{out}>");
  p_ele_EeleOPoutVsPhi->Write();
  TProfile *p_ele_HoEVsEta = h_ele_HoEVsEta->ProfileX();
  p_ele_HoEVsEta->SetTitle("mean ele H/E vs eta");
  p_ele_HoEVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_HoEVsEta->GetYaxis()->SetTitle("<H/E>");
  p_ele_HoEVsEta->Write();
  TProfile *p_ele_HoEVsPhi = h_ele_HoEVsPhi->ProfileX();
  p_ele_HoEVsPhi->SetTitle("mean ele H/E vs phi");
  p_ele_HoEVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_HoEVsPhi->GetYaxis()->SetTitle("<H/E>");
  p_ele_HoEVsPhi->Write();
  TProfile *p_ele_chi2VsEta = h_ele_chi2VsEta->ProfileX();
  p_ele_chi2VsEta->SetTitle("mean ele track chi2 vs eta");
  p_ele_chi2VsEta->GetXaxis()->SetTitle("#eta");
  p_ele_chi2VsEta->GetYaxis()->SetTitle("<#Chi^{2}>");
  p_ele_chi2VsEta->Write();
  TProfile *p_ele_chi2VsPhi = h_ele_chi2VsPhi->ProfileX();
  p_ele_chi2VsPhi->SetTitle("mean ele track chi2 vs phi");
  p_ele_chi2VsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_chi2VsPhi->GetYaxis()->SetTitle("<#Chi^{2}>");
  p_ele_chi2VsPhi->Write();
  TProfile *p_ele_foundHitsVsEta = h_ele_foundHitsVsEta->ProfileX();
  p_ele_foundHitsVsEta->SetTitle("mean ele track # found hits vs eta");
  p_ele_foundHitsVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_foundHitsVsEta->GetYaxis()->SetTitle("<N_{hits}>");
  p_ele_foundHitsVsEta->Write();
  TProfile *p_ele_foundHitsVsPhi = h_ele_foundHitsVsPhi->ProfileX();
  p_ele_foundHitsVsPhi->SetTitle("mean ele track # found hits vs phi");
  p_ele_foundHitsVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_foundHitsVsPhi->GetYaxis()->SetTitle("<N_{hits}>");
  p_ele_foundHitsVsPhi->Write();
  TProfile *p_ele_lostHitsVsEta = h_ele_lostHitsVsEta->ProfileX();
  p_ele_lostHitsVsEta->SetTitle("mean ele track # lost hits vs eta");
  p_ele_lostHitsVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_lostHitsVsEta->GetYaxis()->SetTitle("<N_{hits}>");
  p_ele_lostHitsVsEta->Write();
  TProfile *p_ele_lostHitsVsPhi = h_ele_lostHitsVsPhi->ProfileX();
  p_ele_lostHitsVsPhi->SetTitle("mean ele track # lost hits vs phi");
  p_ele_lostHitsVsPhi->GetXaxis()->SetTitle("#phi (rad)");
  p_ele_lostHitsVsPhi->GetYaxis()->SetTitle("<N_{hits}>");
  p_ele_lostHitsVsPhi->Write();
  TProfile *p_ele_vertexTIPVsEta = h_ele_vertexTIPVsEta->ProfileX();
  p_ele_vertexTIPVsEta->SetTitle("mean tip (wrt gen vtx) vs eta");
  p_ele_vertexTIPVsEta->GetXaxis()->SetTitle("#eta");
  p_ele_vertexTIPVsEta->GetYaxis()->SetTitle("<TIP> (cm)");
  p_ele_vertexTIPVsEta->Write();
  TProfile *p_ele_vertexTIPVsPhi = h_ele_vertexTIPVsPhi->ProfileX();
  p_ele_vertexTIPVsPhi->SetTitle("mean tip (wrt gen vtx) vs phi");
  p_ele_vertexTIPVsPhi->GetXaxis()->SetTitle("#phi");
  p_ele_vertexTIPVsPhi->GetYaxis()->SetTitle("<TIP> (cm)");
  p_ele_vertexTIPVsPhi->Write();
  TProfile *p_ele_vertexTIPVsPt = h_ele_vertexTIPVsPt->ProfileX();
  p_ele_vertexTIPVsPt->SetTitle("mean tip (wrt gen vtx) vs phi");
  p_ele_vertexTIPVsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_ele_vertexTIPVsPt->GetYaxis()->SetTitle("<TIP> (cm)");
  p_ele_vertexTIPVsPt->Write();

  // mc truth
  h_mcNum->Write();
  h_eleNum->Write();
  h_neleinEE->Write();
  h_gamNum->Write();

  // rec event
  histNum_->Write();

  // mc
  h_simEta->Write();
  h_simAbsEta->Write();
  h_simP->Write();
  h_simPt->Write();
  h_simZ->Write();
  h_simPhi->Write();
  h_simPtEta->Write();

  // all electrons
  h_ele_EoverP_all->Write();
  h_ele_EseedOP_all->Write();
  h_ele_EoPout_all->Write();
  h_ele_EeleOPout_all ->Write();
  h_ele_dEtaSc_propVtx_all->Write();
  h_ele_dPhiSc_propVtx_all->Write();
  h_ele_dEtaCl_propOut_all ->Write();
  h_ele_dPhiCl_propOut_all->Write();
  h_ele_HoE_all->Write();
  h_ele_TIP_all->Write();
  h_ele_vertexPt_all->Write();
  h_ele_Et_all->Write();
  h_ele_vertexEta_all->Write();
  h_ele_mee_all->Write();
  h_ele_mee_best_all->Write();
  h_ele_mee_os->Write();
  h_ele_mee_best_os->Write();
  h_ele_mee_seed_all->Write();
  h_ele_mee_seed_os->Write();
  h_ele_mee_sc_all->Write();
  h_ele_mee_sc_os->Write();
  h_ele_mee_newsc_all->Write();
  h_ele_mee_newsc_os->Write();
  h_ele_mee_os_ebeb->Write();
  h_ele_mee_os_ebee->Write();
  h_ele_mee_sc_os_ebee->Write();
  h_ele_mee_newsc_os_ebee->Write();
  h_ele_mee_seed_os_ebee->Write();
  h_ele_mee_os_eeee->Write();
  h_ele_mee_best_os_ebeb->Write();
  h_ele_mee_best_os_ebee->Write();
  h_ele_mee_best_os_eeee->Write();
  h_ele_mee_sc_os_eeee->Write();
  h_ele_mee_newsc_os_eeee->Write();
  h_ele_mee_seed_os_eeee->Write();
  h_ele_mee_os_gg->Write();
  h_ele_mee_os_gb->Write();
  h_ele_mee_os_bb->Write();
  h_ele_mee_seed_os_gg->Write();
  h_ele_mee_seed_os_gg_ebeb->Write();
  h_ele_mee_seed_os_gg_ebee->Write();
  h_ele_mee_seed_os_gg_eeee->Write();
  h_ele_mee_seed_os_gb->Write();
  h_ele_mee_seed_os_bb->Write();
  h_ele_mee_sc_os_gg->Write();
  h_ele_mee_sc_os_gb->Write();
  h_ele_mee_sc_os_bb->Write();
  h_ele_mee_newsc_os_gg->Write();
  h_ele_mee_newsc_os_gg_ebeb->Write();
  h_ele_mee_newsc_os_gg_ebee->Write();
  h_ele_mee_newsc_os_gg_eeee->Write();
  h_ele_mee_newsc_os_gb->Write();
  h_ele_mee_newsc_os_gb_ebee->Write();
  h_ele_mee_newsc_os_gb_eeee->Write();
  h_ele_mee_newsc_os_bb->Write();
  h_ele_mee_newsc_os_bb_ebeb->Write();
  h_ele_mee_newsc_os_bb_ebee->Write();
  h_ele_mee_newsc_os_bb_eeee->Write();
  h_ele_E2mnE1vsMee_all ->Write();
  h_ele_E2mnE1vsMee_egeg_all->Write();

  // charge ID
  h_ele_charge->Write();
  h_ele_simEta_matched_qmisid->Write();
  h_ele_simAbsEta_matched_qmisid->Write();
  h_ele_simPt_matched_qmisid->Write();
  h_ele_simPhi_matched_qmisid->Write();
  h_ele_simZ_matched_qmisid->Write();

  // matched electrons
  h_ele_vertexP->Write();
  h_ele_vertexPt->Write();
  h_ele_Et->Write();
  h_ele_vertexPtVsEta->Write();
  h_ele_vertexPtVsPhi->Write();
  h_ele_simPt_matched->Write();
  h_ele_vertexEta->Write();
  h_ele_vertexEtaVsPhi->Write();
  h_ele_simAbsEta_matched->Write();
  h_ele_simEta_matched->Write();
  h_ele_simPhi_matched->Write();
  h_ele_simPtEta_matched->Write();
  h_ele_vertexPhi->Write();
  h_ele_vertexX->Write();
  h_ele_vertexY ->Write();
  h_ele_vertexZ->Write();
  h_ele_vertexTIP->Write();
  h_ele_simZ_matched->Write();
  h_ele_vertexTIPVsEta->Write();
  h_ele_vertexTIPVsPhi->Write();
  h_ele_vertexTIPVsPt->Write();
  h_ele_PoPtrue->Write();
  h_ele_PoPtrueVsEta ->Write();
  h_ele_PoPtrueVsPhi->Write();
  h_ele_PoPtrueVsPt->Write();
  h_ele_PoPtrue_barrel ->Write();
  h_ele_PoPtrue_endcaps->Write();
  h_ele_PoPtrue_golden_barrel ->Write();
  h_ele_PoPtrue_golden_endcaps->Write();
  h_ele_PoPtrue_showering_barrel ->Write();
  h_ele_PoPtrue_showering_endcaps->Write();
  h_ele_PtoPttrue->Write();
  h_ele_PtoPttrue_barrel ->Write();
  h_ele_PtoPttrue_endcaps->Write();
  h_ele_ChargeMnChargeTrue->Write();
  h_ele_EtaMnEtaTrue->Write();
  h_ele_EtaMnEtaTrue_barrel->Write();
  h_ele_EtaMnEtaTrue_endcaps->Write();
  h_ele_EtaMnEtaTrueVsEta ->Write();
  h_ele_EtaMnEtaTrueVsPhi->Write();
  h_ele_EtaMnEtaTrueVsPt->Write();
  h_ele_PhiMnPhiTrue ->Write();
  h_ele_PhiMnPhiTrue_barrel ->Write();
  h_ele_PhiMnPhiTrue_endcaps ->Write();
  h_ele_PhiMnPhiTrue2 ->Write();
  h_ele_PhiMnPhiTrueVsEta->Write();
  h_ele_PhiMnPhiTrueVsPhi->Write();
  h_ele_PhiMnPhiTrueVsPt->Write();

  // matched electron, superclusters
  histSclEn_->Write();
  histSclEoEtrue_barrel->Write();
  histSclEoEtrue_endcaps->Write();
  histSclEoEtrue_barrel_eg->Write();
  histSclEoEtrue_endcaps_eg->Write();
  histSclEoEtrue_barrel_etagap->Write();
  histSclEoEtrue_barrel_phigap->Write();
  histSclEoEtrue_ebeegap->Write();
  histSclEoEtrue_endcaps->Write();
  histSclEoEtrue_endcaps_golden->Write();
  histSclEoEtrue_endcaps_narrow->Write();
  histSclEoEtrue_endcaps_bigbrem->Write();
  histSclEoEtrue_endcaps_showering->Write();
  histSeedEoEtrue_endcaps_golden->Write();
  histSeedEoEtrue_endcaps_narrow->Write();
  histSeedEoEtrue_endcaps_bigbrem->Write();
  histSeedEoEtrue_endcaps_showering->Write();
  histSclEoEtrue_endcaps_deegap->Write();
  histSclEoEtrue_endcaps_ringgap->Write();
  histSclEoEtruePfVsEg->Write();
  histSclEoEtrue_barrel_new->Write();
  histSclEoEtrue_endcaps_new->Write();
  histSclEoEtrue_barrel_eg_new->Write();
  histSclEoEtrue_endcaps_eg_new->Write();
  histSclEoEtrue_barrel_etagap_new->Write();
  histSclEoEtrue_barrel_phigap_new->Write();
  histSclEoEtrue_ebeegap_new->Write();
  histSclEoEtrue_endcaps_new->Write();
  histSclEoEtrue_endcaps_deegap_new->Write();
  histSclEoEtrue_endcaps_ringgap_new->Write();
  histSclEoEtruePfVsEg->Write();
  histSclEt_->Write();
  histSclEtVsEta_->Write();
  histSclEtVsPhi_->Write();
  histSclEtaVsPhi_ ->Write();
  histSclEta_->Write();
  histSclPhi_->Write();
  histSclSigEtaEta_->Write();
  histSclSigEtaEta_barrel_->Write();
  histSclSigEtaEta_endcaps_->Write();
  histSclSigIEtaIEta_->Write();
  histSclSigIEtaIEta_barrel_->Write();
  histSclSigIEtaIEta_endcaps_->Write();
  histSclE1x5_->Write();
  histSclE1x5_barrel_->Write();
  histSclE1x5_endcaps_->Write();
  histSclE2x5max_->Write();
  histSclE2x5max_barrel_->Write();
  histSclE2x5max_endcaps_->Write();
  histSclE5x5_->Write();
  histSclE5x5_barrel_->Write();
  histSclE5x5_endcaps_->Write();
  histSclSigEtaEta_eg_->Write();
  histSclSigEtaEta_eg_barrel_->Write();
  histSclSigEtaEta_eg_endcaps_->Write();
  histSclSigIEtaIEta_eg_->Write();
  histSclSigIEtaIEta_eg_barrel_->Write();
  histSclSigIEtaIEta_eg_endcaps_->Write();
  histSclE1x5_eg_->Write();
  histSclE1x5_eg_barrel_->Write();
  histSclE1x5_eg_endcaps_->Write();
  histSclE2x5max_eg_->Write();
  histSclE2x5max_eg_barrel_->Write();
  histSclE2x5max_eg_endcaps_->Write();
  histSclE5x5_eg_->Write();
  histSclE5x5_eg_barrel_->Write();
  histSclE5x5_eg_endcaps_->Write();

  // matched electron, gsf tracks
  h_ele_ambiguousTracks->Write();
  h_ele_ambiguousTracksVsEta->Write();
  h_ele_ambiguousTracksVsPhi->Write();
  h_ele_ambiguousTracksVsPt->Write();
  h_ele_foundHits->Write();
  h_ele_foundHits_barrel->Write();
  h_ele_foundHits_endcaps->Write();
  h_ele_foundHitsVsEta->Write();
  h_ele_foundHitsVsPhi->Write();
  h_ele_foundHitsVsPt->Write();
  h_ele_lostHits->Write();
  h_ele_lostHits_barrel->Write();
  h_ele_lostHits_endcaps->Write();
  h_ele_lostHitsVsEta->Write();
  h_ele_lostHitsVsPhi->Write();
  h_ele_lostHitsVsPt->Write();
  h_ele_chi2 ->Write();
  h_ele_chi2_barrel ->Write();
  h_ele_chi2_endcaps ->Write();
  h_ele_chi2VsEta ->Write();
  h_ele_chi2VsPhi ->Write();
  h_ele_chi2VsPt->Write();
  h_ele_PinMnPout->Write();
  h_ele_PinMnPout_mode->Write();
  h_ele_PinMnPoutVsEta_mode->Write();
  h_ele_PinMnPoutVsPhi_mode->Write();
  h_ele_PinMnPoutVsPt_mode->Write();
  h_ele_PinMnPoutVsE_mode->Write();
  h_ele_PinMnPoutVsChi2_mode->Write();
  h_ele_outerP ->Write();
  h_ele_outerP_mode->Write();
  h_ele_outerPVsEta_mode->Write();
  h_ele_outerPt->Write();
  h_ele_outerPt_mode ->Write();
  h_ele_outerPtVsEta_mode->Write();
  h_ele_outerPtVsPhi_mode->Write();
  h_ele_outerPtVsPt_mode->Write();

  // matched electrons, matching
  h_ele_EoP ->Write();
  h_ele_EoP_eg ->Write();
  h_ele_EoP_barrel ->Write();
  h_ele_EoP_eg_barrel ->Write();
  h_ele_EoP_endcaps ->Write();
  h_ele_EoP_eg_endcaps ->Write();
  h_ele_EoPVsEta ->Write();
  h_ele_EoPVsPhi->Write();
  h_ele_EoPVsE->Write();
  h_ele_EseedOP ->Write();
  h_ele_EseedOP_eg ->Write();
  h_ele_EseedOP_barrel ->Write();
  h_ele_EseedOP_eg_barrel ->Write();
  h_ele_EseedOP_endcaps ->Write();
  h_ele_EseedOP_eg_endcaps ->Write();
  h_ele_EseedOPVsEta ->Write();
  h_ele_EseedOPVsPhi->Write();
  h_ele_EseedOPVsE->Write();
  h_ele_EoPout->Write();
  h_ele_EoPout_eg->Write();
  h_ele_EoPout_barrel->Write();
  h_ele_EoPout_eg_barrel->Write();
  h_ele_EoPout_endcaps->Write();
  h_ele_EoPout_eg_endcaps->Write();
  h_ele_EoPoutVsEta->Write();
  h_ele_EoPoutVsPhi->Write();
  h_ele_EoPoutVsE ->Write();
  h_ele_EoPoutVsEta_golden->Write();
  h_ele_EoPoutVsPhi_golden->Write();
  h_ele_EoPoutVsE_golden ->Write();
  h_ele_EoPoutVsEta_showering->Write();
  h_ele_EoPoutVsPhi_showering->Write();
  h_ele_EoPoutVsE_showering ->Write();
  h_ele_EeleOPout->Write();
  h_ele_EeleOPout_eg->Write();
  h_ele_EeleOPout_barrel->Write();
  h_ele_EeleOPout_eg_barrel->Write();
  h_ele_EeleOPout_endcaps->Write();
  h_ele_EoPout_endcaps_golden->Write();
  h_ele_EoPout_endcaps_showering->Write();
  h_ele_EeleOPout_eg_endcaps->Write();
  h_ele_EeleOPoutVsEta->Write();
  h_ele_EeleOPoutVsPhi->Write();
  h_ele_EeleOPoutVsE ->Write();
  h_ele_dEtaSc_propVtx->Write();
  h_ele_dEtaSc_propVtx_eg->Write();
  h_ele_dEtaSc_propVtx_barrel->Write();
  h_ele_dEtaSc_propVtx_eg_barrel->Write();
  h_ele_dEtaSc_propVtx_endcaps->Write();
  h_ele_dEtaSc_propVtx_eg_endcaps->Write();
  h_ele_dEtaScVsEta_propVtx->Write();
  h_ele_dEtaScVsPhi_propVtx->Write();
  h_ele_dEtaScVsPt_propVtx ->Write();
  h_ele_dPhiSc_propVtx->Write();
  h_ele_dPhiSc_propVtx_eg->Write();
  h_ele_dPhiSc_propVtx_barrel->Write();
  h_ele_dPhiSc_propVtx_eg_barrel->Write();
  h_ele_dPhiSc_propVtx_endcaps->Write();
  h_ele_dPhiSc_propVtx_eg_endcaps->Write();
  h_ele_dPhiScVsEta_propVtx ->Write();
  h_ele_dPhiScVsPhi_propVtx->Write();
  h_ele_dPhiScVsPt_propVtx->Write();
  h_ele_dEtaCl_propOut->Write();
  h_ele_dEtaCl_propOut_eg->Write();
  h_ele_dEtaCl_propOut_barrel->Write();
  h_ele_dEtaCl_propOut_eg_barrel->Write();
  h_ele_dEtaCl_propOut_endcaps->Write();
  h_ele_dEtaCl_propOut_eg_endcaps->Write();
  h_ele_dEtaClVsEta_propOut->Write();
  h_ele_dEtaClVsPhi_propOut->Write();
  h_ele_dEtaClVsPt_propOut->Write();
  h_ele_dPhiCl_propOut->Write();
  h_ele_dPhiCl_propOut_eg->Write();
  h_ele_dPhiCl_propOut_barrel->Write();
  h_ele_dPhiCl_propOut_eg_barrel->Write();
  h_ele_dPhiCl_propOut_endcaps->Write();
  h_ele_dPhiCl_propOut_eg_endcaps->Write();
  h_ele_dPhiClVsEta_propOut->Write();
  h_ele_dPhiClVsPhi_propOut->Write();
  h_ele_dPhiClVsPt_propOut->Write();
  h_ele_dEtaEleCl_propOut->Write();
  h_ele_dEtaEleCl_propOut_eg->Write();
  h_ele_dEtaEleCl_propOut_barrel->Write();
  h_ele_dEtaEleCl_propOut_eg_barrel->Write();
  h_ele_dEtaEleCl_propOut_endcaps->Write();
  h_ele_dEtaEleCl_propOut_eg_endcaps->Write();
  h_ele_dEtaEleClVsEta_propOut->Write();
  h_ele_dEtaEleClVsPhi_propOut->Write();
  h_ele_dEtaEleClVsPt_propOut->Write();
  h_ele_dPhiEleCl_propOut->Write();
  h_ele_dPhiEleCl_propOut_eg->Write();
  h_ele_dPhiEleCl_propOut_barrel->Write();
  h_ele_dPhiEleCl_propOut_eg_barrel->Write();
  h_ele_dPhiEleCl_propOut_endcaps->Write();
  h_ele_dPhiEleCl_propOut_eg_endcaps->Write();
  h_ele_dPhiEleClVsEta_propOut->Write();
  h_ele_dPhiEleClVsPhi_propOut->Write();
  h_ele_dPhiEleClVsPt_propOut->Write();
  h_ele_HoE->Write();
  h_ele_HoE_eg->Write();
  h_ele_HoE_barrel->Write();
  h_ele_HoE_eg_barrel->Write();
  h_ele_HoE_endcaps->Write();
  h_ele_HoE_eg_endcaps->Write();
  h_ele_HoE_fiducial->Write();
  h_ele_HoEVsEta->Write();
  h_ele_HoEVsPhi->Write();
  h_ele_HoEVsE->Write();

  h_ele_seed_dphi2_->Write();
  h_ele_seed_subdet2_->Write();
  TProfile *p_ele_seed_dphi2VsEta_ = h_ele_seed_dphi2VsEta_->ProfileX();
  p_ele_seed_dphi2VsEta_->SetTitle("mean ele seed dphi 2nd layer vs eta");
  p_ele_seed_dphi2VsEta_->GetXaxis()->SetTitle("#eta");
  p_ele_seed_dphi2VsEta_->GetYaxis()->SetTitle("<#phi_{pred} - #phi_{hit}, 2nd layer> (rad)");
  p_ele_seed_dphi2VsEta_->SetMinimum(-0.004);
  p_ele_seed_dphi2VsEta_->SetMaximum(0.004);
  p_ele_seed_dphi2VsEta_->Write();
  TProfile *p_ele_seed_dphi2VsPt_ = h_ele_seed_dphi2VsPt_->ProfileX();
  p_ele_seed_dphi2VsPt_->SetTitle("mean ele seed dphi 2nd layer vs pt");
  p_ele_seed_dphi2VsPt_->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_ele_seed_dphi2VsPt_->GetYaxis()->SetTitle("<#phi_{pred} - #phi_{hit}, 2nd layer> (rad)");
  p_ele_seed_dphi2VsPt_->SetMinimum(-0.004);
  p_ele_seed_dphi2VsPt_->SetMaximum(0.004);
  p_ele_seed_dphi2VsPt_->Write();
  h_ele_seed_drz2_->Write();
  TProfile *p_ele_seed_drz2VsEta_ = h_ele_seed_drz2VsEta_->ProfileX();
  p_ele_seed_drz2VsEta_->SetTitle("mean ele seed dr(dz) 2nd layer vs eta");
  p_ele_seed_drz2VsEta_->GetXaxis()->SetTitle("#eta");
  p_ele_seed_drz2VsEta_->GetYaxis()->SetTitle("<r(z)_{pred} - r(z)_{hit}, 2nd layer> (cm)");
  p_ele_seed_drz2VsEta_->SetMinimum(-0.15);
  p_ele_seed_drz2VsEta_->SetMaximum(0.15);
  p_ele_seed_drz2VsEta_->Write();
  TProfile *p_ele_seed_drz2VsPt_ = h_ele_seed_drz2VsPt_->ProfileX();
  p_ele_seed_drz2VsPt_->SetTitle("mean ele seed dr(dz) 2nd layer vs pt");
  p_ele_seed_drz2VsPt_->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_ele_seed_drz2VsPt_->GetYaxis()->SetTitle("<r(z)_{pred} - r(z)_{hit}, 2nd layer> (cm)");
  p_ele_seed_drz2VsPt_->SetMinimum(-0.15);
  p_ele_seed_drz2VsPt_->SetMaximum(0.15);
  p_ele_seed_drz2VsPt_->Write();

  // classes
  h_ele_classes->Write();
  h_ele_old_classes->Write();
  h_ele_eta->Write();
  h_ele_eta_golden->Write();
  h_ele_eta_bbrem->Write();
  h_ele_eta_narrow->Write();
  h_ele_eta_shower->Write();
  h_ele_PinVsPoutGolden_mode->Write();
  h_ele_PinVsPoutShowering_mode->Write();
  h_ele_PinVsPoutGolden_mean->Write();
  h_ele_PinVsPoutShowering_mean->Write();
  h_ele_PtinVsPtoutGolden_mode->Write();
  h_ele_PtinVsPtoutShowering_mode->Write();
  h_ele_PtinVsPtoutGolden_mean->Write();
  h_ele_PtinVsPtoutShowering_mean->Write();
  histSclEoEtrueGolden_barrel->Write();
  histSclEoEtrueGolden_endcaps->Write();
  histSclEoEtrueShowering_barrel->Write();
  histSclEoEtrueShowering_endcaps->Write();

  // fbrem
  h_ele_fbrem->Write();
  h_ele_fbrem_eg->Write();
  h_ele_fbremVsEta_mode->GetXaxis()->SetTitle("#eta");
  h_ele_fbremVsEta_mode->GetYaxis()->SetTitle("<P_{in} - P_{out} / P_{in}>");
  h_ele_fbremVsEta_mode->Write();
  h_ele_fbremVsEta_mean->GetXaxis()->SetTitle("#eta");
  h_ele_fbremVsEta_mean->GetYaxis()->SetTitle("<P_{in} - P_{out} / P_{in}>");
  h_ele_fbremVsEta_mean->Write();
  h_ele_eta_goldenFrac->Write();
  h_ele_eta_bbremFrac->Write();
  h_ele_eta_narrowFrac->Write();
  h_ele_eta_showerFrac->Write();
  h_ele_xOverX0VsEta->Write();

  // efficiencies
  h_ele_etaEff->Write();
  h_ele_zEff->Write();
  h_ele_phiEff->Write();
  h_ele_absetaEff->Write();
  h_ele_ptEff->Write();
  h_ele_ptEtaEff->Write();
  h_ele_etaEff_all->Write();
  h_ele_ptEff_all->Write();

  // q misid
  h_ele_etaQmisid->Write();
  h_ele_zQmisid->Write();
  h_ele_absetaQmisid->Write();
  h_ele_ptQmisid->Write();

  // e/g et pflow electrons
  h_ele_mva->Write();
  h_ele_mva_eg->Write();
  h_ele_provenance->Write();

  // isolation
  h_ele_tkSumPt_dr03->GetXaxis()->SetTitle("TkIsoSum, cone 0.3 (GeV/c)");
  h_ele_tkSumPt_dr03->GetYaxis()->SetTitle("Events");
  h_ele_tkSumPt_dr03->Write();
  h_ele_ecalRecHitSumEt_dr03->GetXaxis()->SetTitle("EcalIsoSum, cone 0.3 (GeV)");
  h_ele_ecalRecHitSumEt_dr03->GetYaxis()->SetTitle("Events");
  h_ele_ecalRecHitSumEt_dr03->Write();
  h_ele_hcalDepth1TowerSumEt_dr03->GetXaxis()->SetTitle("Hcal1IsoSum, cone 0.3 (GeV)");
  h_ele_hcalDepth1TowerSumEt_dr03->GetYaxis()->SetTitle("Events");
  h_ele_hcalDepth1TowerSumEt_dr03->Write();
  h_ele_hcalDepth2TowerSumEt_dr03->GetXaxis()->SetTitle("Hcal2IsoSum, cone 0.3 (GeV)");
  h_ele_hcalDepth2TowerSumEt_dr03->GetYaxis()->SetTitle("Events");
  h_ele_hcalDepth2TowerSumEt_dr03->Write();
  h_ele_tkSumPt_dr04->GetXaxis()->SetTitle("TkIsoSum, cone 0.4 (GeV/c)");
  h_ele_tkSumPt_dr04->GetYaxis()->SetTitle("Events");
  h_ele_tkSumPt_dr04->Write();
  h_ele_ecalRecHitSumEt_dr04->GetXaxis()->SetTitle("EcalIsoSum, cone 0.4 (GeV)");
  h_ele_ecalRecHitSumEt_dr04->GetYaxis()->SetTitle("Events");
  h_ele_ecalRecHitSumEt_dr04->Write();
  h_ele_hcalDepth1TowerSumEt_dr04->GetXaxis()->SetTitle("Hcal1IsoSum, cone 0.4 (GeV)");
  h_ele_hcalDepth1TowerSumEt_dr04->GetYaxis()->SetTitle("Events");
  h_ele_hcalDepth1TowerSumEt_dr04->Write();
  h_ele_hcalDepth2TowerSumEt_dr04->GetXaxis()->SetTitle("Hcal2IsoSum, cone 0.4 (GeV)");
  h_ele_hcalDepth2TowerSumEt_dr04->GetYaxis()->SetTitle("Events");
  h_ele_hcalDepth2TowerSumEt_dr04->Write();

  h_hgcal_sclusters_multiplicity->Write();
  h_hgcal_sclusters_Etsubclusters->Write();
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning->Write();
  h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning_long->Write();
  h_hgcal_scclusters_detadphisubclusters->Write();
  h_hgcal_scclusters_detadphisubclusters_signed->Write();
  h_hgcal_scclusters_detadphisubclusters_pca->Write();
  h_hgcal_scclusters_detadphisubclusters_pca_signed->Write();
  h_hgcal_scclusters_detadphisubclusters_pca_weighted->Write();
  h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted->Write();
  h_hgcal_scclusters_bremkinematic->Write();
  h_hgcal_scclusters_bremkinematic_weighted->Write();
  h_hgcal_scclusters_bremkinematic_etcut->Write();
  h_hgcal_scclusters_bremkinematic_weighted_etcut->Write();
  h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematic_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematic_weighted_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematicnew->Write();
  h_hgcal_scclusters_bremkinematicnew_weighted->Write();
  h_hgcal_scclusters_bremkinematicnew_etcut->Write();
  h_hgcal_scclusters_bremkinematicnew_weighted_etcut->Write();
  h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematicnew_detacut_cleaned->Write();
  h_hgcal_scclusters_bremkinematicnew_weighted_detacut_cleaned->Write();

  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_golden_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_supergolden_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeta_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSfbrem_pt15->Write(); 
  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeoveretrue_pt15->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_golden->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_supergolden->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeta->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVStruevertices->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_showering->Write(); 
  h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_golden->Write(); 
  h_hgcal_seedcluster_cellmultiplicity_golden->Write();   
  h_hgcal_seedcluster_cellenergy_golden->Write(); 
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning->Write(); 
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long->Write(); 
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long_golden->Write(); 
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVSeta->Write(); 
  h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVStruevertices->Write(); 
  h_hgcal_sclusters_Etsubclusters_long->Write();
  h_hgcal_scclusters_detadphisubclusters_pca_signed_long->Write();
  h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted_long->Write();
  h_hgcal_scclusters_bremkinematic_long->Write();
  h_hgcal_scclusters_bremkinematic_weighted_long->Write();
  h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned_long->Write();
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned_long->Write();
  h_hgcal_scclusters_bremkinematic_etcut_detacut_long->Write();
  h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_long->Write();
  
  h_hgcal_sclusters_scenergy->Write();
  h_hgcal_sclusters_scenergy_pt15->Write();
  h_hgcal_sclusters_scenergy_golden_pt15->Write();
  h_hgcal_sclusters_scenergy_showering_pt15->Write();
  h_hgcal_sclusters_scenergy_supergolden_pt15->Write();
  h_hgcal_sclusters_scenergyVSvertices->Write();
  h_hgcal_sclusters_scenergyVStruevertices->Write();
  h_hgcal_sclusters_seedenergy_goldenVSdeltaphiElePin->Write();
  h_hgcal_sclusters_seedenergy_golden->Write();
  h_hgcal_sclusters_seedenergy_golden_3x3->Write();
  h_hgcal_sclusters_seedenergy_golden_5x5->Write();
  h_hgcal_sclusters_seedenergy_golden_10x10->Write();
  h_hgcal_sclusters_seedenergy_supergolden->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt15->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt30->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt1020->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt2030->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt3040->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt4050->Write();
  h_hgcal_sclusters_seedenergy_supergolden_5x5_pt50->Write();
  h_hgcal_sclusters_seedenergy_golden_pt15->Write();
  h_hgcal_sclusters_seedenergy_golden_pt30->Write();
  h_hgcal_sclusters_seedenergy_golden_pt1020->Write();
  h_hgcal_sclusters_seedenergy_golden_pt2030->Write();
  h_hgcal_sclusters_seedenergy_golden_pt3040->Write();
  h_hgcal_sclusters_seedenergy_golden_pt4050->Write();
  h_hgcal_sclusters_seedenergy_golden_pt50->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt30->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt1020->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt2030->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt3040->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt4050->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt50->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden->Write();
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_golden_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt30->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt1020->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt2030->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt3040->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt4050->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt50->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering->Write();
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_showering_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt30->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt1020->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt2030->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt3040->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt4050->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt50->Write();
  h_hgcal_sclusters_newscenergy_etcut->Write();
  h_hgcal_sclusters_newscenergy_dynetcut->Write();
  h_hgcal_sclusters_newscenergy_detacut->Write();
  h_hgcal_sclusters_newscenergy_long->Write();
  h_hgcal_sclusters_newscenergy_etcut_long->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long->Write();
  h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt30->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1530->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1020->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt2030->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt3040->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt4050->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt50->Write();
  h_hgcal_pt_pt1020->Write();
  h_hgcal_pt_pt2030->Write();
  h_hgcal_pt_pt3040->Write();
  h_hgcal_pt_pt4050->Write();
  h_hgcal_pt_pt50->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt40inf->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_frachigh_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fracmiddle_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fraclow_pt15->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_longVSvertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_longVStruevertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVSvertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVStruevertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSvertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVStruevertices->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSeta->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi_aftercut->Write();
  h_hgcal_sclusters_newscenergy_detacut_cleaning->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to120->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu120to130->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu130to140->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to150->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu150to160->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu160to170->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta16to17->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta17to18->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta18to19->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta19to20->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta20to21->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta21to22->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta22to23->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta23to24->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta24to25->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta25to26->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta26to27->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta27to28->Write();
  h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta28to29->Write();

  h_eb_sclusters_scenergy_pt15_pu110to125->Write();
  h_eb_sclusters_scenergy_pt15_pu125to140->Write();
  h_eb_sclusters_scenergy_pt15_pu140to155->Write();
  h_eb_sclusters_scenergy_pt15_pu155to170->Write();
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to125->Write();
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu125to140->Write();
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to155->Write();
  h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu155to170->Write();

  h_ele_hgcal_hoverem->Write();
  h_ele_hgcal_sietaieta->Write();

}

HGCALGsfElectronAnalyzer::~HGCALGsfElectronAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  histfile_->Write();
  histfile_->Close();

}

GlobalVector HGCALGsfElectronAnalyzer::rotateMomentum(const MagneticField *magField,GlobalVector momentum, GlobalPoint xmeas, GlobalPoint xvert, int charge) {

  //std::cout << "[HGCALGsfElectronAnalyzer::rotateMomentum] entering, magfield " << magField << std::endl;
  //double BInTesla = magField->inTesla(xmeas).z();
  double BInTesla = magField->inTesla(GlobalPoint(0.,0.,0.)).z();
  //std::cout << "[HGCALGsfElectronAnalyzer::rotateMomentum] B field " << BInTesla << std::endl;
  GlobalVector xdiff = xmeas - xvert;
  //double theta = xdiff.theta();
  //double phi= xdiff.phi();  
  double pt = momentum.perp();
  double pxOld = momentum.x();
  double pyOld = momentum.y();
  double pz = momentum.z();  
  double RadCurv = 100*pt/(BInTesla*0.29979);
  //std::cout << "[HGCALGsfElectronAnalyzer::rotateMomentum] argument of asin " << 0.5*xdiff.perp()/RadCurv << std::endl;
  double alpha = asin(0.5*xdiff.perp()/RadCurv);
  float ca = cos(charge*alpha);
  float sa = sin(charge*alpha);
  double pxNew = ca*pxOld + sa*pyOld;
  double pyNew = -sa*pxOld + ca*pyOld;
  
  return GlobalVector(pxNew, pyNew, pz);

}

//=========================================================================
// Main method
//=========================================================================

static int oldclassification(const GsfElectron & ele) {

  int iclass;
  
  // old classification
  // iclass = 0 (GOLDEN), 1 (BIGBREM), 2 (BADTRACK), 3 (SHOWERING), 4 (CRACK)
  if (std::abs(ele.superCluster()->eta()-1.5) < 0.1 || std::abs(ele.superCluster()->eta()-3.0) < 0.1)
  {
    iclass = 4 ;
    return iclass;
  }
//   if (ele.isEB()) {
//     if (ele.isEBEEGap() || ele.isEBEtaGap()) {iclass = 4; return iclass;}
//   } else if (ele.isEE()) {        
//     if (std::abs(std::abs(ele.superCluster()->eta())-1.5) < 0.1 ||
//      std::abs(std::abs(ele.superCluster()->eta())-3.0) < 0.1) {iclass = 4; return iclass;}
//   }  
   
  float fbrem = ele.trackFbrem() ;
  int nbrem = ele.numberOfBrems() ;
  if (nbrem == 0 && fbrem < 0.5) // part (pin - scEnergy)/pin < 0.1 removed - M.D.
  { iclass = 0 ; }
  else if (nbrem == 0 && fbrem >= 0.5) // part (pin - scEnergy)/pin < 0.1 removed - M.D.
  { iclass = 1 ; }
  else
  { iclass = 3 ; }

  float pf_fbrem =( ele.superCluster()->energy() - ele.eEleClusterOverPout()*ele.trackMomentumOut().R()) / ele.superCluster()->energy();
//   std::cout << "classification, superClusterFbrem " << ele.superClusterFbrem() << std::endl;
   if (( pf_fbrem-ele.trackFbrem())>=0.15)
   { iclass = 2 ; }    

//  if (ele.isEE()) iclass += 10;
  
  return iclass;
   	  
}	
      
static int classification(const GsfElectron & ele) {

  int iclass;
  
  // new classification
  // redefine golden, narrow and showering 
  // iclass = 0 (GOLDEN), 1 (BIGBREM), 2 (NARROW), 3 (SHOWERING), 4 (CRACK)
  double fbremcut = 0.1;
  if (ele.isEE()) fbremcut = 0.2;
  
  if (ele.isEB()) {
    if (ele.isEBEEGap() || ele.isEBEtaGap()) {iclass = 4; return iclass;}
  } else if (ele.isEE()) {        
    if (std::abs(std::abs(ele.superCluster()->eta())-1.5) < 0.05 ||
     std::abs(std::abs(ele.superCluster()->eta())-3.0) < 0.05) {iclass = 4; return iclass;}
  }  

  if (std::abs(ele.fbrem()) < fbremcut && std::abs(ele.deltaPhiEleClusterTrackAtCalo()) < 0.008) {
    iclass = 0;
  } else if (std::abs(ele.fbrem()) < 0.4) {
    iclass = 2;
  } else if (std::abs(ele.fbrem()) > 0.8) {
    iclass = 1;
//   } else if (std::abs(std::abs(ele.superCluster()->eta())-1.5) < 0.1 ||
//    std::abs(std::abs(ele.superCluster()->eta())-3.0) < 0.1) {
//     iclass = 4;
//     return iclass;
  } else {
    iclass = 3; 
  }	
  
//  if (ele.isEE()) iclass += 10;
  
  return iclass;
   	  
}	

static bool isBad(const GsfElectron & ele) {

  // bad electrons defined as shower+bigbrem+cracks
  bool isbad=false;
  if (classification(ele)==3 || classification(ele)==1 || classification(ele)==4) isbad=true;
  
  return isbad;
  
}


static CaloCluster_iterator seedCluster(const GsfElectron & ele) {

  CaloCluster_iterator itseed = ele.superCluster()->clustersBegin();
  for (CaloCluster_iterator itcl=ele.superCluster()->clustersBegin(); itcl!=ele.superCluster()->clustersEnd(); itcl++) {
    if ((*itcl)->energy() == ele.superCluster()->seed()->energy()) {
      itseed = itcl;
      break;
    }
  }
  
  return itseed;
  
}

static CaloCluster_iterator electronCluster(const GsfElectron & ele) {

  double enrjEle = ele.eEleClusterOverPout()*ele.trackMomentumOut().R();
  CaloCluster_iterator itseed = seedCluster(ele);
  CaloCluster_iterator itelecl = itseed;
  double eps = 1.E-05;

  for (CaloCluster_iterator itcl=ele.superCluster()->clustersBegin(); itcl!=ele.superCluster()->clustersEnd(); itcl++) {

    if (std::abs((*itcl)->energy() - enrjEle)<eps) {
      itelecl = itcl;
      break;
    }

  }
  
  //if ((*itelecl)->energy()*sin(2.*atan(exp(-(*itelecl)->position().eta())))<2.) itelecl = itseed;
  if ((*itelecl)->energy()/std::cosh((*itelecl)->eta())<2.) itelecl = itseed;
  //std::cout << "selected ele cluster energy " << (*itelecl)->energy() << std::endl;
  
  return itelecl;
  
}

static int detector(CaloCluster_iterator itcl) {

  DetId seedXtalId = (*itcl)->hitsAndFractions()[0].first ; 
  int det = seedXtalId.subdetId() ;

//   if (det==EcalBarrel) {
//     std::cout << "new cluster at eta = " << (*itcl)->position().eta() << " cluster is in EcalBarrel " << std::endl;
//    } else if (det == HGCEE) {
//     std::cout << "new cluster at eta = " << (*itcl)->position().eta() << " detector is in HGCEE " << std::endl;
//    } else {
//     std::cout << "new cluster at eta = " << (*itcl)->position().eta() << " unknown detector!! " << std::endl;
//   }
  
  return det;
  
}

static GlobalPoint position(CaloCluster_iterator itcl, const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int det = detector(itcl);
  if (det==EcalBarrel) {

    return GlobalPoint((*itcl)->position().x(),(*itcl)->position().y(),(*itcl)->position().z());

   } else if (det==HGCEE) {

    PCAShowerAnalysis pcaShowerAnalysisEle(iEvent,iSetup);
    return pcaShowerAnalysisEle.showerBarycenter(&(**itcl));	     

  } else {
    
    std::cout << "unknown detector !!" << det << std::endl;
    return GlobalPoint(0.,0.,0.);

  }

}

static double transversal_leakage(int ipart, double energy, GlobalPoint entryPos, GlobalVector pcaShowerDir, 
 const RandomEngine *random_, GammaFunctionGenerator *aGammaGenerator_, EMECALShowerParametrization *showerparam_,
 int matrix_size)
{
  
  // evaluate the fraction of energy leaking out the 5x5 matrix in each layer
  double correction = 1.;

  if (ipart!=11 && ipart!=22) return correction;

  GlobalVector dir = pcaShowerDir.unit();
  const XYZTLorentzVector momentum(dir.x()*energy,dir.y()*energy,dir.z()*energy,energy);
  RawParticle myPart(ipart,momentum);
  std::vector<const RawParticle *> thePart;
  if ( myPart.e() > 0.055 ) thePart.push_back(&myPart);
  EMShower theShower(random_,aGammaGenerator_,showerparam_,&thePart,NULL,NULL,NULL,false); 
  
  return 1;
  
}

static double energycl(CaloCluster_iterator itcl, int ipart, const edm::Event& iEvent, const
 edm::EventSetup& iSetup, edm::Handle<HGCRecHitCollection> recHits, const HGCalGeometry
 *geometry, const RandomEngine *random_, GammaFunctionGenerator *aGammaGenerator_, 
 EMECALShowerParametrization *showerparam_,int matrix=5) {

   //std::cout << "new cluster with initial energy " << (*itcl)->energy() << std::endl;  
  int det = detector(itcl);
  if (det==EcalBarrel) {
    
    return (*itcl)->energy();
  
   } else if (det==HGCEE) {
  
    double energy = 0.;
    GlobalPoint pcaShowerPos;
    GlobalVector pcaShowerDir;

    PCAShowerAnalysis pcaShowerAnalysisSub(iEvent,iSetup); 
    pcaShowerAnalysisSub.showerParameters(&(**itcl),pcaShowerPos,pcaShowerDir);	

//     // for low energy clusters determine the barycenter in each layer
//     if ((*itcl)->energy()<20.) {
//     double x[30], y[30], nrj[30]; 
//     for (int i=0;i<30;i++) {
//       x[i]=0.; y[i]=0.; nrj[i]=0.;
//     }   
//     for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
//       const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
//       HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
//       //std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
//       if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
// 	const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
// 	GlobalPoint cellPos = geometry->getPosition(hgcid_);
// 	// recompute calibrated SC energy
// 	const int layer = hgcid_.layer();
// 	x[layer-1] += cellPos.x()*theHit->energy();
// 	y[layer-1] += cellPos.y()*theHit->energy();
// 	nrj[layer-1] += theHit->energy();
//       }
//     }
//     }
          
    // eta correction
    const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
    // new calib as of SLHC21
    double clus_eta = (*itcl)->eta();
    double corr = _coef_a*std::abs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*std::abs(clus_eta)+_coef_e));
    double mip = 0.0000551;
    double weight[30] =
     {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};
    for (unsigned int ih=0;ih<(*itcl)->hitsAndFractions().size();++ih) {
      const DetId & id_ = ((*itcl)->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
      //std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
	GlobalPoint cellPos = geometry->getPosition(hgcid_);
        //std::cout << "new hit with energy: " << theHit->energy() << " and position " << cellPos<< std::endl; 
 	// recompute calibrated SC energy
	const int layer = hgcid_.layer();
	double scale = mip*corr;
 	// energy as matrix sum around cells intercepted by the shower axis
	//energy += theHit->energy()*weight[layer-1]/scale;
	// rather use the barycenter in each layer
	double lambda = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
	GlobalPoint interceptPos = pcaShowerPos + lambda*pcaShowerDir;	 
	//std::cout << "intercept of shower axis is  " << interceptPos << std::endl;  
	//GlobalPoint barycenter(x[layer-1]/nrj[layer-1],y[layer-1]/nrj[layer-1],cellPos.z());	
	//std::cout << "barycenter in that layer is " << barycenter << std::endl;  
	int matrix_layer = matrix;
	//if (layer<6) matrix_layer = matrix-2;
	if (std::abs(cellPos.x()-interceptPos.x())<0.95*matrix_layer/2. && 
	  (std::abs(cellPos.y()-interceptPos.y())<0.95*matrix_layer/2.)) {
           energy += theHit->energy()*weight[layer-1]/scale;
	  //std::cout << "adding new cell at " <<  cellPos << " with calibrated energy " << theHit->energy()*weight[layer-1]/scale << std::endl;  
        }
      } 
    } 

    //std::cout << "cluster reevaluated energy " << energy << std::endl;  
    // correct for energy leaking outside window
    GlobalPoint entryPos=pcaShowerPos; // to be modified and passed as argument from ecalMomentum 
    double correction = transversal_leakage(ipart,energy,entryPos,pcaShowerDir,random_,aGammaGenerator_,
     showerparam_,matrix);
    return energy/correction;
  
  } else {
    
    std::cout << "unknown detector !!" << det << std::endl;
    return 0.;

  }
   
}

static math::XYZTLorentzVector ecalMomentum(const GsfElectron & ele, const edm::Event& iEvent, const edm::EventSetup& iSetup, 
 edm::Handle<HGCRecHitCollection> recHits, const HGCalGeometry *geometry, bool withPileup_, const RandomEngine *random_,
 GammaFunctionGenerator *aGammaGenerator_, EMECALShowerParametrization *showerparam_, int energy_option=0) {

  // option allows to steer the energy determination
  // option < 10: energy from ele cluster 7x7 matrix for supergolden, as with 20 gfor other electrons
  // option = 20: energy from ele cluster 7x7 matrix and other subclusters 5x5 matrix for all electrons
  
  // This implementation applies the same cleaning treatment in EB appart from the longitudinal
  // In particular the ele cluster only is used there for the definition of the energy for goldens
  // This can probably be improved for |eta|<1 using back the subcluster multiplicity
  
  math::XYZTLorentzVector momentum;
  CaloCluster_iterator itelecl = electronCluster(ele);
  int ncl =0;
  for (CaloCluster_iterator itcl=(ele.superCluster())->clustersBegin(); itcl!=(ele.superCluster())->clustersEnd(); itcl++) ncl++;
  GlobalPoint eleclPos = position(itelecl,iEvent,iSetup);;
  double newenergy = (*itelecl)->energy();
  if (energy_option==20) newenergy =
   energycl(itelecl,11,iEvent,iSetup,recHits,geometry,random_,aGammaGenerator_,showerparam_,7);
  double newx = eleclPos.x()*newenergy;
  double newy = eleclPos.y()*newenergy;
  double newz = eleclPos.z()*newenergy;
  double radius;

  // for the goldens in HGCAL with only one cluster, use the ele cluster with fixed matrix energy
  if (classification(ele)==0 && ncl==1 && energy_option<10) {

    newenergy = energycl(itelecl,11,iEvent,iSetup,recHits,geometry,
     random_,aGammaGenerator_,showerparam_,7); // this function returns cl->energy() in EB
    newx = eleclPos.x()*newenergy;
    newy = eleclPos.y()*newenergy;
    newz = eleclPos.z()*newenergy;
    radius=eleclPos.mag();
    return math::XYZTLorentzVector(newx/radius,newy/radius,newz/radius,newenergy);

   } else { // multicluster patterns

    // instantiate HGCALEmShowerID class
    HGCALShowerBasedEmIdentification
     hGCALShowerBasedEmIdentification(iEvent,iSetup,withPileup_);
    
    for (CaloCluster_iterator itcl=(ele.superCluster())->clustersBegin(); itcl!=(ele.superCluster())->clustersEnd(); itcl++) {

      // skip the ele cluster
      if (itcl==itelecl) continue;
      double clet = (*itcl)->energy()/std::cosh((*itcl)->eta());;
      GlobalPoint clPos = position(itcl,iEvent,iSetup);
      double dphipca = clPos.phi() - eleclPos.phi();
      if (std::abs(dphipca)>CLHEP::pi)
       dphipca = dphipca < 0? (CLHEP::twopi) + dphipca : dphipca - CLHEP::twopi;
      double detapca = std::abs(clPos.eta() - eleclPos.eta());
      // apply PU cleaning cuts
      double etcut = 0.25;
      if (clet>etcut) {
	if (detapca<0.015) {
	  // thight qdphi
	  if (ele.charge()*dphipca>-0.02 || clet>0.5) {
	    if (ele.isEE()) {
              PCAShowerAnalysis pcaShowerAnalysisSub(iEvent,iSetup); 
	      GlobalVector clDir = pcaShowerAnalysisSub.showerAxis(&(**itcl));	
	      hGCALShowerBasedEmIdentification.setShowerPosition(clPos);
	      hGCALShowerBasedEmIdentification.setShowerDirection(clDir); 	 
 	      // cut in length compatibility
	      double x0 = 0.968;
	      // here use parametrisation for pi0 extracted from full sim
	      double predictedLength = 3.6 + 1.383*log((*itcl)->energy());
	      double y = (*itcl)->energy()/0.00536;
	      double sigma = predictedLength / (-2.506+1.245*log(y));
	      GlobalPoint startPos = hGCALShowerBasedEmIdentification.startPosition(&(**itcl));
	      double length = (clPos-startPos).mag();
	      // loose cut
	      bool cutsubcllength = std::abs(predictedLength-length)<4.*sigma/x0;
	      // cut in start position, adapted to pi0s and photons	  
	      bool cutsubclpos = (std::abs(startPos.z())<322.50);
              //!!!!
	      // here remove the longitudinal cuts
	      cutsubclpos = true;
	      cutsubcllength = true;
	      //!!!!
	      if (cutsubcllength && cutsubclpos) {
		double ecl = (*itcl)->energy();
		if (energy_option==20) ecl =
	         energycl(itcl,22,iEvent,iSetup,recHits,geometry,random_,aGammaGenerator_,showerparam_,5);
		newenergy += ecl;
		newx +=clPos.x()*ecl;
		newy +=clPos.y()*ecl;
		newz +=clPos.z()*ecl;
  	      } 
	    } else {
	      double ecl = (*itcl)->energy();
	      //if (energy_option==20) ecl =
	      // energycl(itcl,22,iEvent,iSetup,recHits,geometry,random_,aGammaGenerator_,showerparam_,5);
	      newenergy += ecl;
	      newx +=clPos.x()*ecl;
	      newy +=clPos.y()*ecl;
	      newz +=clPos.z()*ecl;	  
            }  
	  } // end dphi cut
	} // end deta cut
      } // end Et cut

    } // end loop on subclusters 

//     // now correct for the eta variation
//     double feta;
//     feta = 0.988175;
//     if (std::abs(ele.eta())>2.2) feta = 0.86488392 + 0.0560414*std::abs(ele.eta());
//     newenergy  = newenergy / feta;
    
    newx /= newenergy;	 
    newy /= newenergy;	 
    newz /= newenergy;	 
    radius = TMath::Sqrt(newx*newx + newy*newy + newz*newz);
    return math::XYZTLorentzVector(newenergy*newx/radius,newenergy*newy/radius,newenergy*newz/radius,newenergy);

  }  
  
}

static math::XYZTLorentzVector momentum(const GsfElectron & ele, const edm::Event& iEvent, const edm::EventSetup& iSetup, 
 edm::Handle<HGCRecHitCollection> recHits, const HGCalGeometry *geometry, bool withPileup_, 
 const RandomEngine *random_, GammaFunctionGenerator *aGammaGenerator_, EMECALShowerParametrization *showerparam_,
 int energy_option=0) {

  // this implementation uses the track direction rather than ECAL position
  
  // first retreive the (cleaned) energy from the ECAL only momentum
  double energy = ecalMomentum(ele,iEvent,iSetup,recHits,geometry,withPileup_,random_,aGammaGenerator_,
   showerparam_,energy_option).energy(); 

  // use the track for the direction as in the GsfElectron
  math::XYZTLorentzVector momentum = ele.p4();
  momentum *= energy/momentum.e();

  return momentum;
  
}  

static GlobalVector rotateMomentum(const MagneticField *magField,GlobalVector momentum, GlobalPoint xmeas, GlobalPoint xvert, int charge) {

  double BInTesla = magField->inTesla(GlobalPoint(0.,0.,0.)).z();
  //std::cout << "[HGCALElectronClusterAnalyzer::rotateMomentum] B field " << BInTesla << std::endl;
  GlobalVector xdiff = xmeas - xvert;
  double pt = momentum.perp();
  double pxOld = momentum.x();
  double pyOld = momentum.y();
  double pz = momentum.z();  
  double RadCurv = 100*pt/(BInTesla*0.29979);
  double alpha = asin(0.5*xdiff.perp()/RadCurv);
  float ca = cos(charge*alpha);
  float sa = sin(charge*alpha);
  double pxNew = ca*pxOld + sa*pyOld;
  double pyNew = -sa*pxOld + ca*pyOld;
  
  return GlobalVector(pxNew, pyNew, pz);

}

void
HGCALGsfElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "analyzing new event " << std::endl;
  // get electrons
  
//  bool display = true;
//  ievent++;
  
  edm::Handle<GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel(electronCollection_,gsfElectrons);
  edm::LogInfo("")<<"\n\n =================> Treating event "<<iEvent.id()<<" Number of electrons "<<gsfElectrons.product()->size();

  edm::Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(mcTruthCollection_, genParticles);

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

  edm::Handle<VertexCollection> vertices;
  iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"),vertices);

  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  iEvent.getByLabel("addPileupInfo", puInfo);   
  int npu = 0;
  for (auto pvi = puInfo->begin(); pvi!=puInfo->end(); ++pvi) {
    int BX = pvi->getBunchCrossing();
    if(BX == 0) {
    npu = pvi->getPU_NumInteractions();
    }
  }   

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
  HGCALShowerBasedEmIdentification
   hGCALShowerBasedEmIdentificationelId(iEvent,iSetup,withPileup_);

  histNum_->Fill((*gsfElectrons).size());

  // association mc-reco for Zee
  reco::GsfElectronCollection::const_iterator bestRecoElectron;
  //double gsfOkRatio = 10000.;
  double drmin = 0.1;
  double ptcut = 15.;
  double found = false;
  
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {
     
    if ( mcIter->pdgId() == 11) {      
    
      if (mcIter->mother()!=0 && mcIter->mother()->pdgId() == 23) { // matching mohter to Zs
      
	// preselect electrons
	if (mcIter->pt() < ptcut || std::abs(mcIter->eta())> maxAbsEta_) continue;

 	// then remove electrons going to cracks in phi 
	double phideg = std::abs(180.*mcIter->phi()/CLHEP::pi);     
	if ((phideg>8.&&phideg<12.)||(phideg>28.&&phideg<32.)||(phideg>48.&&phideg<52.)||(phideg>68.&&phideg<72.)||
            (phideg>88.&&phideg<92.)||(phideg>108.&&phideg<112.)||(phideg>128.&&phideg<132.)||(phideg>148.&&phideg<152.)||
            (phideg>168.&&phideg<172.)) continue;

        std::cout << "New electron: pdgID " << mcIter->pdgId() << " pt " << mcIter->pt() << " eta " << mcIter->eta() << " mother " << (mcIter->mother() == 0 ? 0. : mcIter->mother()->pdgId()) << std::endl;

	for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
         gsfIter!=gsfElectrons->end(); gsfIter++){

          double dphi = gsfIter->phi()-mcIter->phi();
          if (std::abs(dphi)>CLHEP::pi)
           dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
    	  double deltaR = sqrt(std::pow((gsfIter->eta()-mcIter->eta()),2) + std::pow(dphi,2));
	   if ( deltaR < 0.1 ) {
 	     if ( gsfIter->charge() < 0.) {
	      //double tmpGsfRatio = gsfIter->p()/mcIter->p();
	      //if ( std::abs(tmpGsfRatio-1) < std::abs(gsfOkRatio-1) ) {
	      if ( deltaR < drmin ) {
		//gsfOkRatio = tmpGsfRatio;
		drmin = deltaR;
		bestRecoElectron=gsfIter;
		found = true;
	      }
	    }	  
	  }
	} // end loop on reco electrons
      	
      }
     
    }    

  } // end loop on MC particles
  if (found) std:: cout << "best matching reco electron : " << bestRecoElectron->p4() << std::endl; 
  
  reco::GsfElectronCollection::const_iterator bestRecoPositron;
  //gsfOkRatio = 10000.;
  drmin = 0.1;
  found = false;
  
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {
     
    if ( mcIter->pdgId() == -11) { 
    
      if (mcIter->mother()!=0 && mcIter->mother()->pdgId() == 23) { // matching mohter to Zs
      
	// preselect electrons
	if (mcIter->pt() < ptcut || std::abs(mcIter->eta())> maxAbsEta_) continue;

 	// then remove electrons going to cracks in phi 
	double phideg = std::abs(180.*mcIter->phi()/CLHEP::pi);     
	if ((phideg>8.&&phideg<12.)||(phideg>28.&&phideg<32.)||(phideg>48.&&phideg<52.)||(phideg>68.&&phideg<72.)||
            (phideg>88.&&phideg<92.)||(phideg>108.&&phideg<112.)||(phideg>128.&&phideg<132.)||(phideg>148.&&phideg<152.)||
            (phideg>168.&&phideg<172.)) continue;

        std::cout << "New positron: pdgID " << mcIter->pdgId() << " pt " << mcIter->pt() << " eta " << mcIter->eta() << " mother " << (mcIter->mother() == 0 ? 0. : mcIter->mother()->pdgId()) << std::endl;

	for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
         gsfIter!=gsfElectrons->end(); gsfIter++){

          double dphi = gsfIter->phi()-mcIter->phi();
          if (std::abs(dphi)>CLHEP::pi)
           dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
    	  double deltaR = sqrt(std::pow((gsfIter->eta()-mcIter->eta()),2) + std::pow(dphi,2));
	   if ( deltaR < 0.1 ) {
 	     if ( gsfIter->charge() > 0.) {
	      //double tmpGsfRatio = gsfIter->p()/mcIter->p();
	      //if ( std::abs(tmpGsfRatio-1) < std::abs(gsfOkRatio-1) ) {
	      if ( deltaR < drmin ) {
		//gsfOkRatio = tmpGsfRatio;
		drmin = deltaR;
		bestRecoPositron=gsfIter;
		found = true;
	      }
	    }	  
	  }
	} // end loop on reco electrons
      	
      }
     
    }    

  } // end loop on MC particles
  if (found) std:: cout << "best matching reco positron : " << bestRecoPositron->p4() << std::endl; 

  // all rec electrons
  for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
   gsfIter!=gsfElectrons->end(); gsfIter++){

    // select the matched MC electron
    if (gsfIter!=bestRecoElectron) continue;
    
    // preselect electrons
    //if (std::abs(gsfIter->eta())>maxAbsEta_) continue;
    double scet = gsfIter->superCluster()->energy()*sin(2.*atan(exp(-gsfIter->superCluster()->position().eta())));
    //if (scet<15.) continue;
    std::cout << "new electron with transverse energy " << scet << std::endl;

    double EeleOpout_new = (*electronCluster(*gsfIter))->energy()/gsfIter->trackMomentumOut().R();
    h_ele_EoverP_all     -> Fill( gsfIter->eSuperClusterOverP() );
    h_ele_EseedOP_all     -> Fill( gsfIter->eSeedClusterOverP() );
    h_ele_EoPout_all     -> Fill( gsfIter->eSeedClusterOverPout() );
    //h_ele_EeleOPout_all     -> Fill( gsfIter->eEleClusterOverPout() );
    h_ele_EeleOPout_all     -> Fill( EeleOpout_new );
    h_ele_dEtaSc_propVtx_all -> Fill(gsfIter->deltaEtaSuperClusterTrackAtVtx());
    h_ele_dPhiSc_propVtx_all -> Fill(gsfIter->deltaPhiSuperClusterTrackAtVtx());
    h_ele_dEtaCl_propOut_all -> Fill(gsfIter->deltaEtaSeedClusterTrackAtCalo());
    h_ele_dPhiCl_propOut_all -> Fill(gsfIter->deltaPhiSeedClusterTrackAtCalo());
    h_ele_HoE_all     -> Fill( gsfIter->hadronicOverEm() );
    double d = gsfIter->vertex().x()*gsfIter->vertex().x()+gsfIter->vertex().y()*gsfIter->vertex().y();
    h_ele_TIP_all     -> Fill( sqrt(d) );
    h_ele_vertexEta_all     -> Fill( gsfIter->eta() );
    h_ele_vertexPt_all      -> Fill( gsfIter->pt() );
    h_ele_Et_all      -> Fill( gsfIter->superCluster()->energy()/cosh(gsfIter->superCluster()->eta()));

    // momentum from original SC
    double enrjsc1=gsfIter->superCluster()->energy();
    //double r1sc=TMath::Sqrt(gsfIter->superCluster()->x()*gsfIter->superCluster()->x() + 
    // gsfIter->superCluster()->y()*gsfIter->superCluster()->y() +gsfIter->superCluster()->z()*gsfIter->superCluster()->z());
    //math::XYZTLorentzVector esc1(enrjsc1*gsfIter->superCluster()->x()/r1sc,
    //                             enrjsc1*gsfIter->superCluster()->y()/r1sc,
    // 			         enrjsc1*gsfIter->superCluster()->z()/r1sc,enrjsc1);
    math::XYZTLorentzVector esc1 = gsfIter->p4();
    esc1 *= enrjsc1/esc1.e();
    // momnetum from ele cluster
    CaloCluster_iterator itelecl = electronCluster(*gsfIter);
    double enrj1=(*itelecl)->energy();
    //GlobalPoint eleclPos = position(itelecl, iEvent, iSetup);
    //std::cout << "ele cluster (cmssw) position " << (*itelecl)->position() << 
    // " and PCA position " << eleclPos << std::endl;
    //double r1=eleclPos.mag();
    //math::XYZTLorentzVector e1(enrj1*eleclPos.x()/r1,
    //                           enrj1*eleclPos.y()/r1,
    //			       enrj1*eleclPos.z()/r1,enrj1);
    math::XYZTLorentzVector e1 = gsfIter->p4();
    e1 *= enrj1/e1.e();
    // momentum from ECAL after cleaning with best energy assignment
    math::XYZTLorentzVector enewsc1 = momentum(*gsfIter,iEvent,iSetup,recHits,geometry_,withPileup_,
     random_,aGammaGenerator_,showerparam_,20);
    // momentum from ECAL after cleaning with best energy assignment
    math::XYZTLorentzVector e1best = momentum(*gsfIter,iEvent,iSetup,recHits,geometry_,withPileup_,
     random_,aGammaGenerator_,showerparam_,20);
 
    // mee
    //for (reco::GsfElectronCollection::const_iterator gsfIter2=gsfIter+1;
    // gsfIter2!=gsfElectrons->end(); gsfIter2++){
    for (reco::GsfElectronCollection::const_iterator gsfIter2=gsfElectrons->begin();
     gsfIter2!=gsfElectrons->end(); gsfIter2++){

        // select the matched MC electron
        if (gsfIter2!=bestRecoPositron) continue;

	// preselect electrons
	//if (std::abs(gsfIter2->eta())>maxAbsEta_) continue;
	double scet2 = gsfIter2->superCluster()->energy()*sin(2.*atan(exp(-gsfIter2->superCluster()->position().eta())));
        //if (scet2<15.) continue;
        std::cout << "new positron with transverse energy " << scet2 << std::endl;
	
	// keep only EBEE and EEEE combinations
	//if (gsfIter->isEB() && gsfIter2->isEB()) continue;

        // momentum from original combination    
        math::XYZTLorentzVector p12 = (*gsfIter).p4()+(*gsfIter2).p4();
        float mee2 = p12.Dot(p12);     
	mee2 = std::sqrt(std::abs(mee2));
        // momentum from old SC
        double enrjsc2=gsfIter2->superCluster()->energy();
	//double r2sc=TMath::Sqrt(gsfIter2->superCluster()->x()*gsfIter2->superCluster()->x() + 
        // gsfIter2->superCluster()->y()*gsfIter2->superCluster()->y()
	// +gsfIter2->superCluster()->z()*gsfIter2->superCluster()->z());
        //math::XYZTLorentzVector esc2(enrjsc2*gsfIter2->superCluster()->x()/r2sc,
	// enrjsc2*gsfIter2->superCluster()->y()/r2sc,
	// enrjsc2*gsfIter2->superCluster()->z()/r2sc,enrjsc2);
        math::XYZTLorentzVector esc2 = gsfIter2->p4();
        esc2 *= enrjsc2/esc2.e();
        // momentum from ele cluster
	CaloCluster_iterator itelecl2 = electronCluster(*gsfIter2);
        double enrj2=(*itelecl2)->energy();
	//GlobalPoint eleclPos2 = position(itelecl2, iEvent, iSetup);
	//std::cout << "ele cluster (cmssw) position " << (*itelecl)->position() << 
	// " and PCA position " << eleclPos << std::endl;
	//double r2=eleclPos2.mag();
	//math::XYZTLorentzVector e2(enrj2*eleclPos2.x()/r2,
        //                	   enrj2*eleclPos2.y()/r2,
	//			   enrj2*eleclPos2.z()/r2,enrj2);
        math::XYZTLorentzVector e2 = gsfIter2->p4();
        e2 *= enrj2/e2.e();
        // momentum from ECAL after cleaning with best energy assignment
        math::XYZTLorentzVector enewsc2 = momentum(*gsfIter2,iEvent,iSetup,recHits,geometry_,withPileup_,
	 random_,aGammaGenerator_,showerparam_,20);
        // momentum from ECAL after cleaning with best energy assignment
        math::XYZTLorentzVector e2best = momentum(*gsfIter2,iEvent,iSetup,recHits,geometry_,withPileup_,
	 random_,aGammaGenerator_,showerparam_,20);
        math::XYZTLorentzVector e12 = e1+e2;
	double mee2seed = e12.Dot(e12); 
        math::XYZTLorentzVector e12sc = esc1+esc2;
	double mee2sc = e12sc.Dot(e12sc); 
	math::XYZTLorentzVector e12newsc = enewsc1+enewsc2;
	double mee2newsc = e12newsc.Dot(e12newsc); 
//	// best inclusive energy measurement
//        math::XYZTLorentzVector e1best = enewsc1;
//	if (classification(*gsfIter)==0) e1best = e1;
// 	math::XYZTLorentzVector e2best = enewsc2;
//	if (classification(*gsfIter2)==0) e2best = e2;
	math::XYZTLorentzVector e12best = e1best+e2best;
	double mee2best = e12best.Dot(e12best); 	
       //std::cout << "meeele " << sqrt(mee2) << " mee SC " << sqrt(mee2sc) << " mee new SC "
	//<< sqrt(mee2newsc) <<  " mee seed " << sqrt(mee2seed) <<std::endl;
	h_ele_mee_all -> Fill(sqrt(mee2));
	h_ele_mee_sc_all -> Fill(sqrt(mee2sc));
	h_ele_mee_newsc_all -> Fill(sqrt(mee2newsc));
	h_ele_mee_seed_all -> Fill(sqrt(mee2seed));
	h_ele_mee_best_all -> Fill(sqrt(mee2best));
        h_ele_E2mnE1vsMee_all->Fill(sqrt(mee2),enrj2-enrj1);
        if (gsfIter->ecalDrivenSeed() && gsfIter2->ecalDrivenSeed()) h_ele_E2mnE1vsMee_egeg_all->Fill(sqrt(mee2),enrj2-enrj1);
	if (gsfIter->charge()*gsfIter2->charge()<0.) {
	  h_ele_mee_os -> Fill(sqrt(mee2));
	  h_ele_mee_sc_os -> Fill(sqrt(mee2sc));
	  h_ele_mee_newsc_os -> Fill(sqrt(mee2newsc));
	  h_ele_mee_seed_os -> Fill(sqrt(mee2seed));
	  h_ele_mee_best_os -> Fill(sqrt(mee2best));
	  if (gsfIter->isEE() && gsfIter2->isEE()) h_ele_mee_os_eeee -> Fill(sqrt(mee2));
	  if (gsfIter->isEE() && gsfIter2->isEE()) h_ele_mee_sc_os_eeee -> Fill(sqrt(mee2sc));
	  if (gsfIter->isEE() && gsfIter2->isEE()) h_ele_mee_newsc_os_eeee -> Fill(sqrt(mee2newsc));
	  if (gsfIter->isEE() && gsfIter2->isEE()) h_ele_mee_seed_os_eeee -> Fill(sqrt(mee2seed));
	  if (gsfIter->isEE() && gsfIter2->isEE()) h_ele_mee_best_os_eeee -> Fill(sqrt(mee2best));
	  if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) h_ele_mee_os_ebee -> Fill(sqrt(mee2));
	  if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) h_ele_mee_sc_os_ebee -> Fill(sqrt(mee2sc));
	  if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) h_ele_mee_newsc_os_ebee -> Fill(sqrt(mee2newsc));
	  if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) h_ele_mee_seed_os_ebee -> Fill(sqrt(mee2seed));
	  if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) h_ele_mee_best_os_ebee -> Fill(sqrt(mee2best));
	  if (gsfIter->isEB() && gsfIter2->isEB()) h_ele_mee_best_os_ebeb -> Fill(sqrt(mee2best));
	  if ((classification(*gsfIter)==0 && classification(*gsfIter2)==0))
	   { h_ele_mee_os_gg -> Fill(sqrt(mee2));
//	       h_ele_mee_seed_os_gg -> Fill(sqrt(mee2seed));
	     h_ele_mee_seed_os_gg -> Fill(sqrt(mee2best));
	     h_ele_mee_sc_os_gg -> Fill(sqrt(mee2sc));h_ele_mee_newsc_os_gg -> Fill(sqrt(mee2newsc));
	     if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) {
//	       h_ele_mee_seed_os_gg_ebee -> Fill(sqrt(mee2seed));
	       h_ele_mee_seed_os_gg_ebee -> Fill(sqrt(mee2best));h_ele_mee_newsc_os_gg_ebee -> Fill(sqrt(mee2newsc));
	     } else if (gsfIter->isEE() && gsfIter2->isEE()) {	     
//	       h_ele_mee_seed_os_gg_eeee -> Fill(sqrt(mee2seed));
	       h_ele_mee_seed_os_gg_eeee -> Fill(sqrt(mee2best));h_ele_mee_newsc_os_gg_eeee -> Fill(sqrt(mee2newsc));
	     } else if (gsfIter->isEB() && gsfIter2->isEB()) {	
//	       h_ele_mee_seed_os_gg_ebeb -> Fill(sqrt(mee2seed));
	       h_ele_mee_seed_os_gg_ebeb -> Fill(sqrt(mee2best));h_ele_mee_newsc_os_gg_ebeb -> Fill(sqrt(mee2newsc));
	     }
	   }
	  else if (
// 	     (classification(*gsfIter)==3 && classification(*gsfIter2)==3) ||
// 	     (classification(*gsfIter)==3 && classification(*gsfIter2)==4) ||
// 	     (classification(*gsfIter)==4 && classification(*gsfIter2)==3) ||
// 	     (classification(*gsfIter)==4 && classification(*gsfIter2)==4))
             isBad(*gsfIter) && isBad(*gsfIter2))
	   { h_ele_mee_os_bb -> Fill(sqrt(mee2));h_ele_mee_seed_os_bb -> Fill(sqrt(mee2seed));
	     h_ele_mee_sc_os_bb -> Fill(sqrt(mee2sc));h_ele_mee_newsc_os_bb -> Fill(sqrt(mee2newsc));
	     if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) {
	       h_ele_mee_newsc_os_bb_ebee -> Fill(sqrt(mee2newsc));
	     } else if (gsfIter->isEE() && gsfIter2->isEE()) {	     
	       h_ele_mee_newsc_os_bb_eeee -> Fill(sqrt(mee2newsc));
	     } else if (gsfIter->isEB() && gsfIter2->isEB()) {	     
	       h_ele_mee_newsc_os_bb_ebeb -> Fill(sqrt(mee2newsc));
	     }	     
	   }
	  else if (
//	     (classification(*gsfIter)==0 && classification(*gsfIter2)==3) ||
//	     (classification(*gsfIter)==0 && classification(*gsfIter2)==4) ||
// 	     (classification(*gsfIter)==3 && classification(*gsfIter2)==0) ||
// 	     (classification(*gsfIter)==4 && classification(*gsfIter2)==0))	     
	     (isBad(*gsfIter) && (classification(*gsfIter2)%10)==0) ||
	     (isBad(*gsfIter2) && (classification(*gsfIter)%10)==0))
	   { h_ele_mee_os_gb -> Fill(sqrt(mee2)); h_ele_mee_seed_os_gb -> Fill(sqrt(mee2seed));
	     h_ele_mee_sc_os_gb -> Fill(sqrt(mee2sc));h_ele_mee_newsc_os_gb -> Fill(sqrt(mee2best));
	     if ((gsfIter->isEE() && gsfIter2->isEB()) || (gsfIter->isEB() && gsfIter2->isEE())) {
	       h_ele_mee_newsc_os_gb_ebee -> Fill(sqrt(mee2best));
	     } else if (gsfIter->isEE() && gsfIter2->isEE()) {	     
	       h_ele_mee_newsc_os_gb_eeee -> Fill(sqrt(mee2best));
	     }	     	     
	   }
        }
    }
  }

  int mcNum=0, gamNum=0, eleNum=0;
  bool matchingID, matchingMotherID;
  int neleinEE = 0;
  
  // association mc-reco
  for (reco::GenParticleCollection::const_iterator mcIter=genParticles->begin(); mcIter != genParticles->end(); mcIter++ ) {

    // number of mc particles
    mcNum++;

    // counts photons
    if (mcIter->pdgId() == 22 ){ gamNum++; }

      // select requested matching gen particle
      matchingID=false;
      for (unsigned int i=0; i<matchingIDs_.size(); i++)
       if ( mcIter->pdgId() == matchingIDs_[i] ) matchingID=true;
     //std::cout << "pdgID " << mcIter->pdgId() << " pt " << mcIter->pt() << " eta " << mcIter->eta() << " mother " << (mother == 0 ? 0. : mother->pdgId()) << std::endl;
 
      if (matchingID) {

      // select requested mother matching gen particle
      // always include single particle with no mother
      const Candidate * mother = mcIter->mother();
      matchingMotherID=false;
      for (unsigned int i=0; i<matchingMotherIDs_.size(); i++)
       if ((mother == 0) || ((mother != 0) &&  mother->pdgId() == matchingMotherIDs_[i]) ) matchingMotherID=true;

      if (matchingMotherID) {

      if (mcIter->pt()> maxPt_ || std::abs(mcIter->eta())> maxAbsEta_) continue;
      //std::cout << "matching ID " << mcIter->pdgId() << " mother ID " << (mother == 0 ? 0. : mother->pdgId()) << std::endl; 

      // suppress the barrel
//      if (std::abs(mcIter->eta()) < 1.5) continue;
      if (std::abs(mcIter->eta()) > 1.5 && std::abs(mcIter->eta())<3.0) {
        neleinEE++;
      }
      
      // select central z
      //if ( std::abs(mcIter->production_vertex()->position().z())>50.) continue;
      
      eleNum++;
      h_simEta -> Fill( mcIter->eta() );
      h_simAbsEta -> Fill( std::abs(mcIter->eta()) );
      h_simP   -> Fill( mcIter->p() );
      h_simPt   -> Fill( mcIter->pt() );
      h_simPhi   -> Fill( mcIter->phi() );
      h_simZ   -> Fill( mcIter->vz() );
      h_simPtEta   -> Fill( mcIter->eta(),mcIter->pt() );

      // looking for the best matching gsf electron
      bool okGsfFound = false;
      double gsfOkRatio = 999999.;

      // find best matched electron
      reco::GsfElectron bestGsfElectron;
      for (reco::GsfElectronCollection::const_iterator gsfIter=gsfElectrons->begin();
       gsfIter!=gsfElectrons->end(); gsfIter++){

        double dphi = gsfIter->phi()-mcIter->phi();
        if (std::abs(dphi)>CLHEP::pi)
         dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
    	double deltaR = sqrt(std::pow((gsfIter->eta()-mcIter->eta()),2) + std::pow(dphi,2));
	if ( deltaR < deltaR_ ){
//	  if ( ( (mcIter->pdgId() == 11) && (gsfIter->charge() < 0.) ) ||
//	       ( (mcIter->pdgId() == -11) && (gsfIter->charge() > 0.) ) )
//	   {
	    double tmpGsfRatio = gsfIter->p()/mcIter->p();
	    if ( std::abs(tmpGsfRatio-1) < std::abs(gsfOkRatio-1) ) {
	      gsfOkRatio = tmpGsfRatio;
	      bestGsfElectron=*gsfIter;
	      okGsfFound = true;
	    }
//	  }
	}
      } // loop over rec ele to look for the best one

      // analysis when the mc track is found
     if (okGsfFound){

	// histos for electrons in EB or in EE

	//int eleClass = bestGsfElectron.classification();
	int eleClass = classification(bestGsfElectron);
	int eleOldClass = oldclassification(bestGsfElectron);
	h_ele_classes ->Fill(eleClass);
	h_ele_old_classes ->Fill(eleOldClass);

	// histos for EB only
	double newscenergyeb = ecalMomentum(bestGsfElectron,iEvent,iSetup,recHits,geometry_,withPileup_,
	 random_,aGammaGenerator_,showerparam_,10).e();
        double fnpu = float(npu);
	if (mcIter->pt()>15.) {
 	    if ((bestGsfElectron.isEB() && !bestGsfElectron.isEBEEGap() && !bestGsfElectron.isEBEtaGap())) {
	    if (fnpu>=110.&&fnpu<125.) h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to125->Fill(newscenergyeb/mcIter->p());
	    if (fnpu>=125.&&fnpu<140.) h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu125to140->Fill(newscenergyeb/mcIter->p());
	    if (fnpu>=140.&&fnpu<155.) h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to155->Fill(newscenergyeb/mcIter->p());
	    if (fnpu>=155.&&fnpu<170.) h_eb_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu155to170->Fill(newscenergyeb/mcIter->p());
	    if (fnpu>=110.&&fnpu<125.) h_eb_sclusters_scenergy_pt15_pu110to125->Fill(bestGsfElectron.superCluster()->energy()/mcIter->p());
	    if (fnpu>=125.&&fnpu<140.) h_eb_sclusters_scenergy_pt15_pu125to140->Fill(bestGsfElectron.superCluster()->energy()/mcIter->p());
	    if (fnpu>=140.&&fnpu<155.) h_eb_sclusters_scenergy_pt15_pu140to155->Fill(bestGsfElectron.superCluster()->energy()/mcIter->p());
	    if (fnpu>=155.&&fnpu<170.) h_eb_sclusters_scenergy_pt15_pu155to170->Fill(bestGsfElectron.superCluster()->energy()/mcIter->p());
	  }
	}
	
        // fiducial region in HGCAL EE
	
	// remove electrons going to cracks in eta and cracks in phi
	bool etacrack=false;
	//double etamax=2.7;
	double etamax=2.9;
        if (std::abs(mcIter->eta()) < 1.6 || std::abs(mcIter->eta())>etamax) etacrack=true;		
	// then remove electrons going to cracks in phi 
	bool phicrack=false;
	double phideg = std::abs(180.*mcIter->phi()/CLHEP::pi);     
	if ((phideg>8.&&phideg<12.)||(phideg>28.&&phideg<32.)||(phideg>48.&&phideg<52.)||(phideg>68.&&phideg<72.)||
            (phideg>88.&&phideg<92.)||(phideg>108.&&phideg<112.)||(phideg>128.&&phideg<132.)||(phideg>148.&&phideg<152.)||
            (phideg>168.&&phideg<172.)) phicrack=true;

	if (bestGsfElectron.isEE() && !phicrack && !etacrack) {
	
	// electron related distributions
	h_ele_charge        -> Fill( bestGsfElectron.charge() );
	h_ele_chargeVsEta        -> Fill( bestGsfElectron.eta(),bestGsfElectron.charge() );
	h_ele_chargeVsPhi        -> Fill( bestGsfElectron.phi(),bestGsfElectron.charge() );
	h_ele_chargeVsPt        -> Fill( bestGsfElectron.pt(),bestGsfElectron.charge() );
	h_ele_vertexP       -> Fill( bestGsfElectron.p() );
	h_ele_vertexPt      -> Fill( bestGsfElectron.pt() );
        h_ele_Et      -> Fill( bestGsfElectron.superCluster()->energy()/cosh(bestGsfElectron.superCluster()->eta()));
	h_ele_vertexPtVsEta      -> Fill(  bestGsfElectron.eta(),bestGsfElectron.pt() );
	h_ele_vertexPtVsPhi      -> Fill(  bestGsfElectron.phi(),bestGsfElectron.pt() );
	h_ele_vertexEta     -> Fill( bestGsfElectron.eta() );
	// generated distributions for matched electrons
	h_ele_simPt_matched      -> Fill( mcIter->pt() );
        h_ele_simPhi_matched   -> Fill( mcIter->phi() );
	h_ele_simAbsEta_matched     -> Fill( std::abs(mcIter->eta()) );
	h_ele_simEta_matched     -> Fill( mcIter->eta() );
	h_ele_simPtEta_matched      -> Fill(  mcIter->eta(),mcIter->pt() );
	h_ele_vertexEtaVsPhi     -> Fill(  bestGsfElectron.phi(),bestGsfElectron.eta() );
	h_ele_vertexPhi     -> Fill( bestGsfElectron.phi() );
	h_ele_vertexX     -> Fill( bestGsfElectron.vertex().x() );
	h_ele_vertexY     -> Fill( bestGsfElectron.vertex().y() );
	h_ele_vertexZ     -> Fill( bestGsfElectron.vertex().z() );
        h_ele_simZ_matched   -> Fill( mcIter->vz() );
	double d = (bestGsfElectron.vertex().x()-mcIter->vx())
	          *(bestGsfElectron.vertex().x()-mcIter->vx())+
	           (bestGsfElectron.vertex().y()-mcIter->vy())
		  *(bestGsfElectron.vertex().y()-mcIter->vy());
	d = sqrt(d);
	h_ele_vertexTIP     -> Fill( d );
	h_ele_vertexTIPVsEta     -> Fill(  bestGsfElectron.eta(), d );
	h_ele_vertexTIPVsPhi     -> Fill(  bestGsfElectron.phi(), d );
	h_ele_vertexTIPVsPt     -> Fill(  bestGsfElectron.pt(), d );
	h_ele_EtaMnEtaTrue  -> Fill( bestGsfElectron.eta()-mcIter->eta());
	if (bestGsfElectron.isEB()) h_ele_EtaMnEtaTrue_barrel  -> Fill( bestGsfElectron.eta()-mcIter->eta());
	if (bestGsfElectron.isEE()) h_ele_EtaMnEtaTrue_endcaps  -> Fill( bestGsfElectron.eta()-mcIter->eta());
	h_ele_EtaMnEtaTrueVsEta  -> Fill( bestGsfElectron.eta(), bestGsfElectron.eta()-mcIter->eta());
	h_ele_EtaMnEtaTrueVsPhi  -> Fill( bestGsfElectron.phi(), bestGsfElectron.eta()-mcIter->eta());
	h_ele_EtaMnEtaTrueVsPt  -> Fill( bestGsfElectron.pt(), bestGsfElectron.eta()-mcIter->eta());
	h_ele_PhiMnPhiTrue  -> Fill( bestGsfElectron.phi()-mcIter->phi());
	if (bestGsfElectron.isEB()) h_ele_PhiMnPhiTrue_barrel  -> Fill( bestGsfElectron.phi()-mcIter->phi());
	if (bestGsfElectron.isEE()) h_ele_PhiMnPhiTrue_endcaps  -> Fill( bestGsfElectron.phi()-mcIter->phi());
	h_ele_PhiMnPhiTrue2  -> Fill( bestGsfElectron.phi()-mcIter->phi());
	h_ele_PhiMnPhiTrueVsEta  -> Fill( bestGsfElectron.eta(), bestGsfElectron.phi()-mcIter->phi());
	h_ele_PhiMnPhiTrueVsPhi  -> Fill( bestGsfElectron.phi(), bestGsfElectron.phi()-mcIter->phi());
	h_ele_PhiMnPhiTrueVsPt  -> Fill( bestGsfElectron.pt(), bestGsfElectron.phi()-mcIter->phi());
	h_ele_PoPtrue       -> Fill( bestGsfElectron.p()/mcIter->p());
	h_ele_PtoPttrue       -> Fill( bestGsfElectron.pt()/mcIter->pt());
	h_ele_PoPtrueVsEta       -> Fill( bestGsfElectron.eta(), bestGsfElectron.p()/mcIter->p());
	h_ele_PoPtrueVsPhi       -> Fill( bestGsfElectron.phi(), bestGsfElectron.p()/mcIter->p());
	h_ele_PoPtrueVsPt       -> Fill( bestGsfElectron.py(), bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEB()) h_ele_PoPtrue_barrel       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEE()) h_ele_PoPtrue_endcaps       -> Fill( bestGsfElectron.p()/mcIter->p());
//	if (bestGsfElectron.isEB() && bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PoPtrue_golden_barrel       -> Fill( bestGsfElectron.p()/mcIter->p());
//	if (bestGsfElectron.isEE() && bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PoPtrue_golden_endcaps       -> Fill( bestGsfElectron.p()/mcIter->p());
//	if (bestGsfElectron.isEB() && bestGsfElectron.classification() == GsfElectron::SHOWERING) h_ele_PoPtrue_showering_barrel       -> Fill( bestGsfElectron.p()/mcIter->p());
//	if (bestGsfElectron.isEE() && bestGsfElectron.classification() == GsfElectron::SHOWERING) h_ele_PoPtrue_showering_endcaps       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEB() && classification(bestGsfElectron)==0) h_ele_PoPtrue_golden_barrel       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEE() && classification(bestGsfElectron)==0) h_ele_PoPtrue_golden_endcaps       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEB() && classification(bestGsfElectron)==3) h_ele_PoPtrue_showering_barrel       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEE() && classification(bestGsfElectron)==3) h_ele_PoPtrue_showering_endcaps       -> Fill( bestGsfElectron.p()/mcIter->p());
	if (bestGsfElectron.isEB()) h_ele_PtoPttrue_barrel       -> Fill( bestGsfElectron.pt()/mcIter->pt());
	if (bestGsfElectron.isEE()) h_ele_PtoPttrue_endcaps       -> Fill( bestGsfElectron.pt()/mcIter->pt());

	// supercluster related distributions
	reco::SuperClusterRef sclRef = bestGsfElectron.superCluster();
	if (!bestGsfElectron.ecalDrivenSeed()&&bestGsfElectron.trackerDrivenSeed()) sclRef = bestGsfElectron.pflowSuperCluster();
        histSclEn_->Fill(sclRef->energy());
        double R=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
        double Rt=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
        histSclEt_->Fill(sclRef->energy()*(Rt/R));
        histSclEtVsEta_->Fill(sclRef->eta(),sclRef->energy()*(Rt/R));
        histSclEtVsPhi_->Fill(sclRef->phi(),sclRef->energy()*(Rt/R));
        if (bestGsfElectron.isEB())  histSclEoEtrue_barrel->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE())  histSclEoEtrue_endcaps->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.ecalDrivenSeed())  histSclEoEtrue_barrel_eg->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && bestGsfElectron.ecalDrivenSeed())  histSclEoEtrue_endcaps_eg->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.isEBEtaGap())  histSclEoEtrue_barrel_etagap->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.isEBPhiGap())  histSclEoEtrue_barrel_phigap->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEBEEGap())  histSclEoEtrue_ebeegap->Fill(sclRef->energy()/mcIter->p());
        //if (bestGsfElectron.isEE())  histSclEoEtrue_endcaps->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && bestGsfElectron.isEEDeeGap())  histSclEoEtrue_endcaps_deegap->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==4)  histSclEoEtrue_endcaps_ringgap->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB())  histSclEoEtrue_barrel_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE())  histSclEoEtrue_endcaps_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.ecalDrivenSeed())  histSclEoEtrue_barrel_eg_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && bestGsfElectron.ecalDrivenSeed())  histSclEoEtrue_endcaps_eg_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.isEBEtaGap())  histSclEoEtrue_barrel_etagap_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEB() && bestGsfElectron.isEBPhiGap())  histSclEoEtrue_barrel_phigap_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEBEEGap())  histSclEoEtrue_ebeegap_new->Fill(sclRef->energy()/mcIter->p());
        //if (bestGsfElectron.isEE())  histSclEoEtrue_endcaps_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && bestGsfElectron.isEEDeeGap())  histSclEoEtrue_endcaps_deegap_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==4)  histSclEoEtrue_endcaps_ringgap_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==4)  histSclEoEtrue_endcaps_ringgap_new->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==0) histSclEoEtrue_endcaps_golden->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==2) histSclEoEtrue_endcaps_narrow->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==1) histSclEoEtrue_endcaps_bigbrem->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==3) histSclEoEtrue_endcaps_showering->Fill(sclRef->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==0) histSeedEoEtrue_endcaps_golden->Fill(sclRef->seed()->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==2) histSeedEoEtrue_endcaps_narrow->Fill(sclRef->seed()->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==1) histSeedEoEtrue_endcaps_bigbrem->Fill(sclRef->seed()->energy()/mcIter->p());
        if (bestGsfElectron.isEE() && classification(bestGsfElectron)==3) histSeedEoEtrue_endcaps_showering->Fill(sclRef->seed()->energy()/mcIter->p());
        histSclEta_->Fill(sclRef->eta());
        histSclEtaVsPhi_->Fill(sclRef->phi(),sclRef->eta());
        histSclPhi_->Fill(sclRef->phi());
        histSclSigEtaEta_->Fill(bestGsfElectron.scSigmaEtaEta());
        if (bestGsfElectron.isEB()) histSclSigEtaEta_barrel_->Fill(bestGsfElectron.scSigmaEtaEta());
        if (bestGsfElectron.isEE()) histSclSigEtaEta_endcaps_->Fill(bestGsfElectron.scSigmaEtaEta());
        histSclSigIEtaIEta_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        if (bestGsfElectron.isEB()) histSclSigIEtaIEta_barrel_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        if (bestGsfElectron.isEE()) histSclSigIEtaIEta_endcaps_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        histSclE1x5_->Fill(bestGsfElectron.scE1x5());
        if (bestGsfElectron.isEB()) histSclE1x5_barrel_->Fill(bestGsfElectron.scE1x5());
        if (bestGsfElectron.isEE()) histSclE1x5_endcaps_->Fill(bestGsfElectron.scE1x5());
        histSclE2x5max_->Fill(bestGsfElectron.scE2x5Max());
        if (bestGsfElectron.isEB()) histSclE2x5max_barrel_->Fill(bestGsfElectron.scE2x5Max());
        if (bestGsfElectron.isEE()) histSclE2x5max_endcaps_->Fill(bestGsfElectron.scE2x5Max());
        histSclE5x5_->Fill(bestGsfElectron.scE5x5());
        if (bestGsfElectron.isEB()) histSclE5x5_barrel_->Fill(bestGsfElectron.scE5x5());
        if (bestGsfElectron.isEE()) histSclE5x5_endcaps_->Fill(bestGsfElectron.scE5x5());
        if (bestGsfElectron.ecalDrivenSeed()) histSclSigIEtaIEta_eg_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) histSclSigIEtaIEta_eg_barrel_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) histSclSigIEtaIEta_eg_endcaps_->Fill(bestGsfElectron.scSigmaIEtaIEta());
        if (bestGsfElectron.ecalDrivenSeed())histSclE1x5_eg_->Fill(bestGsfElectron.scE1x5());
        if (bestGsfElectron.isEB() && bestGsfElectron.ecalDrivenSeed())histSclE1x5_eg_barrel_->Fill(bestGsfElectron.scE1x5());
        if (bestGsfElectron.isEE() && bestGsfElectron.ecalDrivenSeed())histSclE1x5_eg_endcaps_->Fill(bestGsfElectron.scE1x5());
        if (bestGsfElectron.ecalDrivenSeed())histSclE2x5max_eg_->Fill(bestGsfElectron.scE2x5Max());
        if (bestGsfElectron.isEB() && bestGsfElectron.ecalDrivenSeed())histSclE2x5max_eg_barrel_->Fill(bestGsfElectron.scE2x5Max());
        if (bestGsfElectron.isEE() && bestGsfElectron.ecalDrivenSeed())histSclE2x5max_eg_endcaps_->Fill(bestGsfElectron.scE2x5Max());
        if (bestGsfElectron.ecalDrivenSeed())histSclE5x5_eg_->Fill(bestGsfElectron.scE5x5());
        if (bestGsfElectron.isEB() && bestGsfElectron.ecalDrivenSeed())histSclE5x5_eg_barrel_->Fill(bestGsfElectron.scE5x5());
        if (bestGsfElectron.isEE() && bestGsfElectron.ecalDrivenSeed())histSclE5x5_eg_endcaps_->Fill(bestGsfElectron.scE5x5());
        float pfEnergy=0., egEnergy=0.;
	if (!bestGsfElectron.superCluster().isNull()) egEnergy = bestGsfElectron.superCluster()->energy();
	if (!bestGsfElectron.pflowSuperCluster().isNull()) pfEnergy = bestGsfElectron.pflowSuperCluster()->energy();
	histSclEoEtruePfVsEg->Fill(egEnergy/mcIter->p(),pfEnergy/mcIter->p());

	// suclusters in the SC
        //int detector = sclRef.seed()->hitsAndFractions()[0].first.subdetId() ;
        //int component = sclRef.seed()->hitsAndFractions()[0].first.det() ;
        //if (component==DetId::Forward && detector==HGCEE) {

	// first identify the ele cluster
	CaloCluster_iterator itelecl = electronCluster(bestGsfElectron);
	// get its PCA position
	GlobalPoint eleclPos = position(itelecl,iEvent,iSetup);	  
	double eleclet = (*itelecl)->energy()*sin(2.*atan(exp(-eleclPos.eta())));
	std::cout << "  new ele cluster in HGCAL with energy " << (*itelecl)->energy() << 
	  " , transverse energy " << eleclet << 
	  " , position (x,y,z) " << (*itelecl)->position() << " and position (eta,phi) : (" 
	   << (*itelecl)->position().eta() << "," << (*itelecl)->position().phi() << std::endl;
	std::cout << "  PCA position " << eleclPos << std::endl;
	
	// produce HGCAL elID variables for Et>10
	PCAShowerAnalysis pcaShowerAnalysisElId(iEvent,iSetup);
	GlobalPoint clPos = pcaShowerAnalysisElId.showerBarycenter(&(*sclRef->seed()));	
	GlobalVector clDir = pcaShowerAnalysisElId.showerAxis(&(*sclRef->seed()));	  	
	hGCALShowerBasedEmIdentificationelId.setShowerPosition(clPos);
	hGCALShowerBasedEmIdentificationelId.setShowerDirection(clDir); 	 
        if (sclRef->energy()*(Rt/R)>10.) {
	double hoverem = hGCALShowerBasedEmIdentificationelId.hadOverEm(&(*sclRef->seed()),"all");
	h_ele_hgcal_hoverem->Fill(hoverem);
	double sietaieta = hGCALShowerBasedEmIdentificationelId.sigmaetaeta(&(*sclRef->seed())); 
	h_ele_hgcal_sietaieta->Fill(sietaieta);
        }
	
	// now loop on all subclusters
	int nsubclusters_etcut_detacut_cleaning=1;
	int nsubclusters_etcut_detacut_cleaning_long=1;
	int ncells_etcut_detacut_cleaning_long=(*itelecl)->hitsAndFractions().size();
	double newenergy_etcut = (*itelecl)->energy();
	double newenergy_dynetcut = (*itelecl)->energy();
	double newenergy_detacut = (*itelecl)->energy();
	double newenergy_long = (*itelecl)->energy();
	double newenergy_etcut_long = (*itelecl)->energy();
	double newenergy_etcut_detacut_cleaning = (*itelecl)->energy();
	double newenergy_etcut_detacut_cleaning_long = (*itelecl)->energy();
	double newenergy_etcut_detacut_long = (*itelecl)->energy();
	double newenergy_detacut_cleaning = (*itelecl)->energy();
	double etcut = 0.25;
	double detacut = 0.015;
	for (CaloCluster_iterator itcl=sclRef->clustersBegin(); itcl!=sclRef->clustersEnd(); itcl++) {
          // skip the electron cluster 	  
	  if (itcl==itelecl) continue;
	  //std::cout << "  new sub cluster in HGCAL with energy " << (*itcl)->energy() << 
	  //" , transverse energy " << clet << 
	  //" , position (x,y,z) " << (*itcl)->position() << " and position (eta,phi) : (" 
	  // << (*itcl)->position().eta() << "," << (*itcl)->position().phi() << std::endl;	    
	  // here call the PCA directly instead of the static function to avoid filling the
	  // PCA twice
	  PCAShowerAnalysis pcaShowerAnalysisSub(iEvent,iSetup);
	  GlobalPoint clPos = pcaShowerAnalysisSub.showerBarycenter(&(**itcl));	
	  GlobalVector clDir = pcaShowerAnalysisSub.showerAxis(&(**itcl));	  
	  std::cout << "  ele cluster direction " << clDir <<  std::endl;
	  //std::cout << "  PCA position " << clPos << std::endl;
	  double clet = (*itcl)->energy()*sin(2.*atan(exp(-clPos.eta())));
	  h_hgcal_sclusters_Etsubclusters->Fill(clet);
	  double deta = (*itelecl)->eta() - (*itcl)->eta();
	  double dphi = (*itelecl)->phi() - (*itcl)->phi();
	  if (std::abs(dphi)>CLHEP::pi)
           dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	  double detapca = clPos.eta() - eleclPos.eta();
	  double dphipca = clPos.phi() - eleclPos.phi();
	  if (std::abs(dphipca)>CLHEP::pi)
           dphipca = dphipca < 0? (CLHEP::twopi) + dphipca : dphipca - CLHEP::twopi;
	  h_hgcal_scclusters_detadphisubclusters->Fill(deta,dphi);
	  h_hgcal_scclusters_detadphisubclusters_pca->Fill(detapca,dphipca);	  
	  h_hgcal_scclusters_detadphisubclusters_pca_weighted->Fill(detapca,dphipca,(*itcl)->energy());	  
	  h_hgcal_scclusters_detadphisubclusters_signed->Fill(deta,bestGsfElectron.charge()*dphi);
	  h_hgcal_scclusters_detadphisubclusters_pca_signed->Fill(detapca,bestGsfElectron.charge()*dphipca);
	  h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted->Fill(detapca,bestGsfElectron.charge()*dphipca,(*itcl)->energy());
	  double dphigbar = (*itcl)->phi() - sclRef->phi();
	  if (std::abs(dphigbar)>CLHEP::pi)
           dphigbar = dphigbar < 0? (CLHEP::twopi) + dphigbar : dphigbar - CLHEP::twopi;
	  double dphigele = ((*itcl)->phi() - (*itelecl)->phi());
	  if (std::abs(dphigele)>CLHEP::pi)
           dphigele = dphigele < 0? (CLHEP::twopi) + dphigele : dphigele - CLHEP::twopi;
	  //double dphiratio = dphigbar / dphigele;
	  double dphiratio = dphigele;
	  //double energyratio = (*itcl)->energy()/((*itelecl)->energy()+(*itcl)->energy());
	  double energyratio = (*itcl)->energy()/sclRef->energy();
// 	  h_hgcal_scclusters_bremkinematic->Fill(bestGsfElectron.charge()*dphigbar,energyratio);
// 	  h_hgcal_scclusters_bremkinematic_weighted->Fill(bestGsfElectron.charge()*dphigbar,energyratio,(*itcl)->energy());
// 	  h_hgcal_scclusters_bremkinematicnew->Fill(bestGsfElectron.charge()*dphiratio,energyratio);
// 	  h_hgcal_scclusters_bremkinematicnew_weighted->Fill(bestGsfElectron.charge()*dphiratio,energyratio,(*itcl)->energy());
	  h_hgcal_scclusters_bremkinematic->Fill(bestGsfElectron.charge()*dphigele,energyratio);
	  h_hgcal_scclusters_bremkinematic_weighted->Fill(bestGsfElectron.charge()*dphigele,energyratio,(*itcl)->energy());
	  h_hgcal_scclusters_bremkinematicnew->Fill(bestGsfElectron.charge()*dphiratio,energyratio);
	  h_hgcal_scclusters_bremkinematicnew_weighted->Fill(bestGsfElectron.charge()*dphiratio,energyratio,(*itcl)->energy());
	  // apply em loose ID
	  hGCALShowerBasedEmIdentification.setShowerPosition(clPos);
	  hGCALShowerBasedEmIdentification.setShowerDirection(clDir); 	 
          //bool cutsetaeta = hGCALShowerBasedEmIdentification.cutSigmaetaeta(&(**itelecl));
	  //bool cuthadem =  hGCALShowerBasedEmIdentification.cutHadOverEm(&(**itelecl));
 	  //bool cutseedpos = hGCALShowerBasedEmIdentification.cutStartPosition(&(**itelecl));
	  //bool cutseedlength = hGCALShowerBasedEmIdentification.cutLengthCompatibility(&(**itelecl));
          //if (cutsetaeta && cuthadem && cutseedpos && cutseedlength) {
	  // cut in length compatibility
	  double x0 = 0.968;
	  // here use parametrisation for pi0 extracted from full sim
	  double predictedLength = 3.6 + 1.383*log((*itcl)->energy());
	  double y = (*itcl)->energy()/0.00536;
	  double sigma = predictedLength / (-2.506+1.245*log(y));
	  GlobalPoint startPos = hGCALShowerBasedEmIdentification.startPosition(&(**itcl));
	  double length = (clPos-startPos).mag();
	  // tighter cut
	  bool cutsubcllength = std::abs(predictedLength-length)<4.*sigma/x0;
	  //bool cutsubcllength = std::abs(predictedLength-length)<3.*sigma/x0;
	  //bool cutsubcllength = std::abs(predictedLength-length)<2.*sigma/x0;
	  // cut in start position, adapted to pi0s and photons	  
	  bool cutsubclpos = (std::abs(startPos.z())<322.50);
          //!!!!!
	  // here remove the longitudinal cuts 
	  cutsubclpos = true;
	  cutsubcllength = true;
	  //!!!!!
	  //if (cuthadem) {
          if (cutsubcllength && cutsubclpos) {
	    newenergy_long += (*itcl)->energy();
	    if (clet>etcut) newenergy_etcut_long += (*itcl)->energy();
	    h_hgcal_sclusters_Etsubclusters_long->Fill(clet);
	    h_hgcal_scclusters_detadphisubclusters_pca_signed_long->Fill(detapca,bestGsfElectron.charge()*dphipca);	  
	    h_hgcal_scclusters_detadphisubclusters_pca_signed_weighted_long->Fill(detapca,bestGsfElectron.charge()*dphipca,(*itcl)->energy());
	    // also apply qdphi cut for the long plots
	    if (bestGsfElectron.charge()*dphipca>-0.02 || clet>0.5) {
// 	    h_hgcal_scclusters_bremkinematic_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy());
// 	    h_hgcal_scclusters_bremkinematic_weighted_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());
	    h_hgcal_scclusters_bremkinematic_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy());
	    h_hgcal_scclusters_bremkinematic_weighted_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());
	    }
	  }
	  // apply cleaning deta and dphi cut
	  if (std::abs(detapca)<detacut) {
	    newenergy_detacut += (*itcl)->energy();
	    if (bestGsfElectron.charge()*dphipca>-0.02 || clet>0.5) {
	      newenergy_detacut_cleaning += (*itcl)->energy();
// 	      h_hgcal_scclusters_bremkinematic_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigbar,energyratio);
// 	      h_hgcal_scclusters_bremkinematic_weighted_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigbar,energyratio,(*itcl)->energy());	    
// 	      h_hgcal_scclusters_bremkinematicnew_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,energyratio);
// 	      h_hgcal_scclusters_bremkinematicnew_weighted_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,energyratio,(*itcl)->energy());	    
	      h_hgcal_scclusters_bremkinematic_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigele,energyratio);
	      h_hgcal_scclusters_bremkinematic_weighted_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigele,energyratio,(*itcl)->energy());	    
	      h_hgcal_scclusters_bremkinematicnew_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,energyratio);
	      h_hgcal_scclusters_bremkinematicnew_weighted_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,energyratio,(*itcl)->energy());	    
	    }
	  }
          // test dynamic Et threshold
          double dynetcut = std::min(0.5+(1.5-0.5)*(std::abs(dphipca)-0.15)/(0.3-0.15),1.5);  
	  if (clet>dynetcut) {
	    newenergy_dynetcut += (*itcl)->energy();	  
	  }
	  // apply all cut
	  if (clet>etcut) {
	    newenergy_etcut += (*itcl)->energy();
	    h_hgcal_scclusters_bremkinematic_etcut->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy());
	    h_hgcal_scclusters_bremkinematic_weighted_etcut->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());
	    if (std::abs(detapca)<detacut) {
	      if (cutsubcllength && cutsubclpos) {
	        newenergy_etcut_detacut_long += (*itcl)->energy();
//	        h_hgcal_scclusters_bremkinematic_etcut_detacut_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy());
//	        h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
	        h_hgcal_scclusters_bremkinematic_etcut_detacut_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy());
	        h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
	      }
	      if (bestGsfElectron.charge()*dphipca>-0.02 || clet>0.5) {
	        nsubclusters_etcut_detacut_cleaning++;
		newenergy_etcut_detacut_cleaning += (*itcl)->energy();
	        h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning->Fill(clet);
// 		h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy());
// 		h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
// 		h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,(*itcl)->energy()/sclRef->energy());
// 		h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
		h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy());
		h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
		h_hgcal_scclusters_bremkinematicnew_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,(*itcl)->energy()/sclRef->energy());
		h_hgcal_scclusters_bremkinematicnew_weighted_etcut_detacut_cleaned->Fill(bestGsfElectron.charge()*dphiratio,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
		if (cutsubcllength && cutsubclpos) {
	          nsubclusters_etcut_detacut_cleaning_long++;
	          ncells_etcut_detacut_cleaning_long += (*itcl)->hitsAndFractions().size();
	          newenergy_etcut_detacut_cleaning_long += (*itcl)->energy();
	          h_hgcal_sclusters_Etsubclusters_etcut_detacut_cleaning_long->Fill(clet);
// 	          h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy());
// 	          h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned_long->Fill(bestGsfElectron.charge()*dphigbar,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
	          h_hgcal_scclusters_bremkinematic_etcut_detacut_cleaned_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy());
	          h_hgcal_scclusters_bremkinematic_weighted_etcut_detacut_cleaned_long->Fill(bestGsfElectron.charge()*dphigele,(*itcl)->energy()/sclRef->energy(),(*itcl)->energy());	    
		}
	      }
	    }
	  }
	}

        // do the energy response vs phi
	double newscenergy = ecalMomentum(bestGsfElectron,iEvent,iSetup,recHits,geometry_,withPileup_,
	 random_,aGammaGenerator_,showerparam_,15).e();
	double newscenergy5x5 = ecalMomentum(bestGsfElectron,iEvent,iSetup,recHits,geometry_,withPileup_,
	 random_,aGammaGenerator_,showerparam_,20).e();
	// here force usage of newscenergy5x5 in all hitsos
	newscenergy = newscenergy5x5;
	h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi->Fill(180.*mcIter->phi()/CLHEP::pi,newscenergy/mcIter->p());

// 	// remove electrons going to cracks in eta
// 	// note: all this part has fiducial requirement 1.5<|eta|<2.7
// 	bool etacrack=false;
//         if (std::abs(mcIter->eta()) < 1.6 || std::abs(mcIter->eta())>2.9) etacrack=true;
// 		
// 	// then remove electrons going to cracks in phi 
// 	bool phicrack=false;
// 	double phideg = std::abs(180.*mcIter->phi()/CLHEP::pi);     
// 	if ((phideg>8.&&phideg<12.)||(phideg>28.&&phideg<32.)||(phideg>48.&&phideg<52.)||(phideg>68.&&phideg<72.)||
//             (phideg>88.&&phideg<92.)||(phideg>108.&&phideg<112.)||(phideg>128.&&phideg<132.)||(phideg>148.&&phideg<152.)||
//             (phideg>168.&&phideg<172.)) phicrack=true;
// 
// 	if (!phicrack && !etacrack) {
// 	
	h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSphi_aftercut->Fill(180.*mcIter->phi()/CLHEP::pi,newscenergy/mcIter->p());
      
	h_hgcal_sclusters_scenergy->Fill(sclRef->energy()/mcIter->p());
	if (mcIter->pt()>15.) {
	  h_hgcal_sclusters_scenergy_pt15->Fill(sclRef->energy()/mcIter->p());
	  if (classification(bestGsfElectron)==0) h_hgcal_sclusters_scenergy_golden_pt15->Fill(sclRef->energy()/mcIter->p());
	  if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1)) h_hgcal_sclusters_scenergy_supergolden_pt15->Fill(sclRef->energy()/mcIter->p());
	  if (isBad(bestGsfElectron)) h_hgcal_sclusters_scenergy_showering_pt15->Fill(sclRef->energy()/mcIter->p());
	}
	
	h_hgcal_sclusters_scenergyVSvertices->Fill(vertices->size(),sclRef->energy()/mcIter->p());
	h_hgcal_sclusters_scenergyVStruevertices->Fill(float(npu),sclRef->energy()/mcIter->p());
	h_hgcal_sclusters_multiplicity->Fill(float(sclRef->clusters().size()));
	
	h_hgcal_sclusters_newscenergy_etcut->Fill(newenergy_etcut/mcIter->p());
	h_hgcal_sclusters_newscenergy_dynetcut->Fill(newenergy_dynetcut/mcIter->p());
	h_hgcal_sclusters_newscenergy_detacut->Fill(newenergy_detacut/mcIter->p());
	h_hgcal_sclusters_newscenergy_long->Fill(newenergy_long/mcIter->p());
	h_hgcal_sclusters_newscenergy_etcut_long->Fill(newenergy_etcut_long/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_long->Fill(newenergy_etcut_detacut_long/mcIter->p());
	h_hgcal_sclusters_newscenergy_detacut_cleaning->Fill(newenergy_detacut_cleaning/mcIter->p());
	h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning->Fill(newenergy_etcut_detacut_cleaning/mcIter->p());
        
//  	// test new classification implementation in ElectronClassification
// 	std::cout << "electron old classification is " <<  bestGsfElectron.classification() << std::endl; 
// 	std::cout << "electron new classification is " <<  classification(bestGsfElectron) << std::endl; 
// 	ElectronClassification standaloneClassification;
// 	standaloneClassification.classify(bestGsfElectron);
// 	std::cout << "electron new standalone classification is " <<  bestGsfElectron.classification() << std::endl; 
//         if (classification(bestGsfElectron) !=  bestGsfElectron.classification()) 
// 	 std::cout << "ALERT!! new standalone classification is giving a different result!!! " << std::endl; 
//  	
	// new fixed matrix energy assignment for goldens
        double ecl = (*itelecl)->energy();
	double ecl3x3 = energycl(itelecl,11,iEvent,iSetup,recHits,geometry_,random_,aGammaGenerator_,showerparam_,3);
        double ecl5x5 = energycl(itelecl,11,iEvent,iSetup,recHits,geometry_,random_,aGammaGenerator_,showerparam_,5);
        double ecl10x10 = energycl(itelecl,11,iEvent,iSetup,recHits,geometry_,random_,aGammaGenerator_,showerparam_,10);	
 	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_seedenergy_golden->Fill(ecl/mcIter->p());
 	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_seedenergy_golden_3x3->Fill(ecl3x3/mcIter->p());
 	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_seedenergy_golden_5x5->Fill(ecl5x5/mcIter->p());
 	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_seedenergy_golden_10x10->Fill(ecl10x10/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>15.) h_hgcal_sclusters_seedenergy_golden_pt15->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>30.) h_hgcal_sclusters_seedenergy_golden_pt30->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>10. && mcIter->pt()<20.) h_hgcal_sclusters_seedenergy_golden_pt1020->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>20. && mcIter->pt()<30.) h_hgcal_sclusters_seedenergy_golden_pt2030->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>30. && mcIter->pt()<40.) h_hgcal_sclusters_seedenergy_golden_pt3040->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>40. && mcIter->pt()<50.) h_hgcal_sclusters_seedenergy_golden_pt4050->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>50. ) h_hgcal_sclusters_seedenergy_golden_pt50->Fill(ecl5x5/mcIter->p());
	
	// supergolden
 	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1)) h_hgcal_sclusters_seedenergy_supergolden->Fill(ecl/mcIter->p());
 	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1)) h_hgcal_sclusters_seedenergy_supergolden_5x5->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>15.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt15->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>30.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt30->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>10. && mcIter->pt()<20.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt1020->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>20. && mcIter->pt()<30.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt2030->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>30. && mcIter->pt()<40.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt3040->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>40. && mcIter->pt()<50.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt4050->Fill(ecl5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>50.) h_hgcal_sclusters_seedenergy_supergolden_5x5_pt50->Fill(ecl5x5/mcIter->p());

	// cleaned sc energy for all electrons (standalone function)
	//double newscenergy = newenergy_etcut_detacut_cleaning_long;
 	h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>15.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>15.) h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_pt15->Fill(newscenergy5x5/mcIter->p());
        if (mcIter->pt()>30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt30->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>15. && mcIter->pt()<30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1530->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>10. && mcIter->pt()<20.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt1020->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>20. && mcIter->pt()<30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt2030->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>30. && mcIter->pt()<40.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt3040->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>40. && mcIter->pt()<50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt4050->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt50->Fill(newscenergy/mcIter->p());
        if (mcIter->pt()>10. && mcIter->pt()<20.) h_hgcal_pt_pt1020->Fill(mcIter->pt());
        if (mcIter->pt()>20. && mcIter->pt()<30.) h_hgcal_pt_pt2030->Fill(mcIter->pt());
        if (mcIter->pt()>30. && mcIter->pt()<40.) h_hgcal_pt_pt3040->Fill(mcIter->pt());
        if (mcIter->pt()>40. && mcIter->pt()<50.) h_hgcal_pt_pt4050->Fill(mcIter->pt());
        if (mcIter->pt()>50.) h_hgcal_pt_pt50->Fill(mcIter->pt());
        if (mcIter->pt()>40.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt40inf->Fill(newscenergy5x5/mcIter->p());
        //double fnpu = float(npu);
	if (mcIter->pt()>15.) {
          if (std::abs(bestGsfElectron.superCluster()->eta())>1.6 && std::abs(bestGsfElectron.superCluster()->eta())<2.9) {
	    if (fnpu>=110.&&fnpu<120.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu110to120->Fill(newscenergy/mcIter->p());
	    if (fnpu>=120.&&fnpu<130.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu120to130->Fill(newscenergy/mcIter->p());
	    if (fnpu>=130.&&fnpu<140.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu130to140->Fill(newscenergy/mcIter->p());
	    if (fnpu>=140.&&fnpu<150.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu140to150->Fill(newscenergy/mcIter->p());
	    if (fnpu>=150.&&fnpu<160.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu150to160->Fill(newscenergy/mcIter->p());
	    if (fnpu>=160.&&fnpu<170.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_pu160to170->Fill(newscenergy/mcIter->p());
	  }
	  if (std::abs(mcIter->eta())>1.6 && std::abs(mcIter->eta())<1.7) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta16to17->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>1.7 && std::abs(mcIter->eta())<1.8) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta17to18->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>1.8 && std::abs(mcIter->eta())<1.9) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta18to19->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>1.9 && std::abs(mcIter->eta())<2.0) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta19to20->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.0 && std::abs(mcIter->eta())<2.1) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta20to21->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.1 && std::abs(mcIter->eta())<2.2) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta21to22->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.2 && std::abs(mcIter->eta())<2.3) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta22to23->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.3 && std::abs(mcIter->eta())<2.4) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta23to24->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.4 && std::abs(mcIter->eta())<2.5) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta24to25->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.5 && std::abs(mcIter->eta())<2.6) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta25to26->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.6 && std::abs(mcIter->eta())<2.7) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta26to27->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.7 && std::abs(mcIter->eta())<2.8) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta27to28->Fill(newscenergy/mcIter->p());
	  if (std::abs(mcIter->eta())>2.8 && std::abs(mcIter->eta())<2.9) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_pt15_eta28to29->Fill(newscenergy/mcIter->p());
	}
	
	// cleaned sc energy for bad electrons
	if (isBad(bestGsfElectron)) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>15.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt15->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>15.) h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_showering_pt15->Fill(newscenergy5x5/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt30->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>10.&& mcIter->pt()<20.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt1020->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>20.&& mcIter->pt()<30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt2030->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>30.&& mcIter->pt()<40.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt3040->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>40.&& mcIter->pt()<50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt4050->Fill(newscenergy/mcIter->p());
	if (isBad(bestGsfElectron) && mcIter->pt()>50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_showering_pt50->Fill(newscenergy/mcIter->p());
        
	// cleaned sc energy for golden electrons
	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>15.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt15->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>15.) h_hgcal_sclusters_newscenergy5x5_etcut_detacut_cleaning_long_golden_pt15->Fill(newscenergy5x5/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt30->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>10.&& mcIter->pt()<20.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt1020->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>20.&& mcIter->pt()<30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt2030->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>30.&& mcIter->pt()<40.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt3040->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>40.&& mcIter->pt()<50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt4050->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && mcIter->pt()>50.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_golden_pt50->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1)) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>15.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt15->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>30.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt30->Fill(newscenergy/mcIter->p());
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1) && mcIter->pt()>30.&& mcIter->pt()<40.) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_supergolden_pt3040->Fill(newscenergy/mcIter->p());
	
	// cleaned sc energy vs pileup
        h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSvertices->Fill(vertices->size(),newscenergy/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVStruevertices->Fill(float(npu),newscenergy/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_longVSeta->Fill(std::abs(mcIter->eta()),newscenergy/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVSvertices->Fill(vertices->size(),newscenergy/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_cleaningVStruevertices->Fill(float(npu),newscenergy/mcIter->p());

        // other cleaning tests
	h_hgcal_sclusters_newscenergy_etcut_detacut_longVSvertices->Fill(vertices->size(),newenergy_etcut_detacut_long/mcIter->p());
        h_hgcal_sclusters_newscenergy_etcut_detacut_longVStruevertices->Fill(float(npu),newenergy_etcut_detacut_long/mcIter->p());

	// energy fraction
	if (mcIter->pt()>15.) {
	  double fraction = (*itelecl)->energy()/newscenergy;
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_pt15->Fill(fraction);
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_pt15->Fill(fraction);
	  if (classification(bestGsfElectron)==0) h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_golden_pt15->Fill(fraction);
	  if (classification(bestGsfElectron)==0 && nsubclusters_etcut_detacut_cleaning_long==1)
	   h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_long_supergolden_pt15->Fill(fraction);
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeta_pt15->Fill(std::abs(mcIter->eta()),fraction);
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVStruevertices_pt15->Fill(float(npu),fraction);
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSfbrem_pt15->Fill(bestGsfElectron.fbrem(),fraction);
	  h_hgcal_sclusters_energyfraction_etcut_detacut_cleaning_longVSeoveretrue_pt15->Fill(newscenergy/mcIter->p(),fraction);
          if (fraction>0.9) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_frachigh_pt15->Fill(newscenergy/mcIter->p());
          if (fraction>0.7 && fraction<0.9) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fracmiddle_pt15->Fill(newscenergy/mcIter->p());
          if (fraction<0.7) h_hgcal_sclusters_newscenergy_etcut_detacut_cleaning_long_fraclow_pt15->Fill(newscenergy/mcIter->p());
	}
	
	// multiplicity
	h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning->Fill(float(nsubclusters_etcut_detacut_cleaning));
	h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long->Fill(float(nsubclusters_etcut_detacut_cleaning_long));
	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_long_golden->Fill(float(nsubclusters_etcut_detacut_cleaning_long));
	h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVSeta->Fill(std::abs(mcIter->eta()),float(nsubclusters_etcut_detacut_cleaning_long));
	h_hgcal_sclusters_multiplicity_etcut_detacut_cleaning_longVStruevertices->Fill(float(npu),float(nsubclusters_etcut_detacut_cleaning_long));
	
	// cell multiplicity
	h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long->Fill(float(ncells_etcut_detacut_cleaning_long));
	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_golden->Fill(float(ncells_etcut_detacut_cleaning_long));
	if (classification(bestGsfElectron)==0 && (nsubclusters_etcut_detacut_cleaning_long==1)) h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_long_supergolden->Fill(float(ncells_etcut_detacut_cleaning_long));
	h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeta->Fill(std::abs(mcIter->eta()),float(ncells_etcut_detacut_cleaning_long));
	h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVStruevertices->Fill(float(npu),float(ncells_etcut_detacut_cleaning_long));
	h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue->Fill(newscenergy/mcIter->p(),float(ncells_etcut_detacut_cleaning_long));
	if (classification(bestGsfElectron)==0) h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_golden->Fill(newscenergy/mcIter->p(),float(ncells_etcut_detacut_cleaning_long));
	if (classification(bestGsfElectron)==3) h_hgcal_sclusters_cellmultiplicity_etcut_detacut_cleaning_longVSeoveretrue_showering->Fill(newscenergy/mcIter->p(),float(ncells_etcut_detacut_cleaning_long));
	
	// electron cluster (golden)
	if (classification(bestGsfElectron)==0) h_hgcal_seedcluster_cellmultiplicity_golden->Fill(float((*itelecl)->hitsAndFractions().size()));
	for (unsigned int ih=0;ih<(*itelecl)->hitsAndFractions().size();++ih) {
	  const DetId & id_ = ((*itelecl)->hitsAndFractions())[ih].first ;
	  HGCRecHitCollection::const_iterator theHit = recHits->find(id_);    
	  //std::cout << "new hit with energy: " << theHit->energy() << std::endl; 
	  if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	    if (classification(bestGsfElectron)==0)
	    h_hgcal_seedcluster_cellenergy_golden->Fill(theHit->energy());
	  }
	}  
	
	if (classification(bestGsfElectron)==0) {
 	  GlobalVector direction(bestGsfElectron.trackMomentumAtVtx().x(),bestGsfElectron.trackMomentumAtVtx().y(),
	                         bestGsfElectron.trackMomentumAtVtx().z());
 	  GlobalPoint origin(bestGsfElectron.vertex().x(),bestGsfElectron.vertex().y(),
	                         bestGsfElectron.vertex().z());
	  double dphi =(*itelecl)->phi() - rotateMomentum(&(*magField_),direction,eleclPos,origin,bestGsfElectron.charge()).phi();
	  if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	  h_hgcal_sclusters_seedenergy_goldenVSdeltaphiElePin->Fill(dphi,(*itelecl)->energy()/mcIter->p());
	}
	
//	}
	
// 	// shower displays
// 	if (display) {
// 
// 	  if (ievent<10) {
// 
// 	  // analyze all recHits
// 	  HGCRecHitCollection::const_iterator ih;   
// 	  for (ih=recHits->begin(); ih!=recHits->end(); ih++) {
// 	    if ((*ih).detid().det()==DetId::Forward && (*ih).detid().subdetId()==HGCEE) {
// 	      const HGCEEDetId & hgcid_ = HGCEEDetId((*ih).detid()) ;
// 	      //GlobalPoint cellPos = geometry_->getGeometry(hgcid_)->getPosition();
// 	      GlobalPoint cellPos = geometry_->getPosition(hgcid_);
// 	      //std::cout << "new HGCAL recHits with position " <<  cellPos << "and energy " << ih->energy() << std::endl;
// 	      // only z>0 and x>0 and y>0
// 	      if (cellPos.z()>0. && cellPos.x()>0. && cellPos.y()>0.) {
// 	       h_hgcal_allhits_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z());
// 	       h_hgcal_allhits_weighted_em[ievent]->Fill(cellPos.x(),cellPos.y(),cellPos.z(),ih->energy()*1000.);
// 	      }
// 	    }
// 	  }
// 
// 	  }
// 
// 	}
// 	
	// track related distributions
	h_ele_ambiguousTracks     -> Fill( bestGsfElectron.ambiguousGsfTracksSize() );
	h_ele_ambiguousTracksVsEta     -> Fill( bestGsfElectron.eta(), bestGsfElectron.ambiguousGsfTracksSize() );
	h_ele_ambiguousTracksVsPhi     -> Fill( bestGsfElectron.phi(), bestGsfElectron.ambiguousGsfTracksSize() );
	h_ele_ambiguousTracksVsPt     -> Fill( bestGsfElectron.pt(), bestGsfElectron.ambiguousGsfTracksSize() );
	if (!readAOD_) { // track extra does not exist in AOD
	  h_ele_foundHits     -> Fill( bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  if (bestGsfElectron.isEB()) h_ele_foundHits_barrel     -> Fill( bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  if (bestGsfElectron.isEE()) h_ele_foundHits_endcaps     -> Fill( bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  h_ele_foundHitsVsEta     -> Fill( bestGsfElectron.eta(), bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  h_ele_foundHitsVsPhi     -> Fill( bestGsfElectron.phi(), bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  h_ele_foundHitsVsPt     -> Fill( bestGsfElectron.pt(), bestGsfElectron.gsfTrack()->numberOfValidHits() );
	  h_ele_lostHits      -> Fill( bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  if (bestGsfElectron.isEB()) h_ele_lostHits_barrel      -> Fill( bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  if (bestGsfElectron.isEE()) h_ele_lostHits_endcaps      -> Fill( bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  h_ele_lostHitsVsEta      -> Fill( bestGsfElectron.eta(), bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  h_ele_lostHitsVsPhi      -> Fill( bestGsfElectron.phi(), bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  h_ele_lostHitsVsPt      -> Fill( bestGsfElectron.pt(), bestGsfElectron.gsfTrack()->numberOfLostHits() );
	  h_ele_chi2          -> Fill( bestGsfElectron.gsfTrack()->normalizedChi2() );
	  if (bestGsfElectron.isEB()) h_ele_chi2_barrel          -> Fill( bestGsfElectron.gsfTrack()->normalizedChi2() );
	  if (bestGsfElectron.isEE()) h_ele_chi2_endcaps          -> Fill( bestGsfElectron.gsfTrack()->normalizedChi2() );
	  h_ele_chi2VsEta          -> Fill( bestGsfElectron.eta(), bestGsfElectron.gsfTrack()->normalizedChi2() );
	  h_ele_chi2VsPhi          -> Fill( bestGsfElectron.phi(), bestGsfElectron.gsfTrack()->normalizedChi2() );
	  h_ele_chi2VsPt          -> Fill( bestGsfElectron.pt(), bestGsfElectron.gsfTrack()->normalizedChi2() );
	}
	// from gsf track interface, hence using mean
	if (!readAOD_) { // track extra does not exist in AOD
	  h_ele_PinMnPout     -> Fill( bestGsfElectron.gsfTrack()->innerMomentum().R() - bestGsfElectron.gsfTrack()->outerMomentum().R() );
	  h_ele_outerP        -> Fill( bestGsfElectron.gsfTrack()->outerMomentum().R() );
	  h_ele_outerPt       -> Fill( bestGsfElectron.gsfTrack()->outerMomentum().Rho() );
        }
	// from electron interface, hence using mode
	h_ele_PinMnPout_mode     -> Fill( bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	h_ele_PinMnPoutVsEta_mode     -> Fill(  bestGsfElectron.eta(), bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	h_ele_PinMnPoutVsPhi_mode     -> Fill(  bestGsfElectron.phi(), bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	h_ele_PinMnPoutVsPt_mode     -> Fill(  bestGsfElectron.pt(), bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	h_ele_PinMnPoutVsE_mode     -> Fill(  bestGsfElectron.caloEnergy(), bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	if (!readAOD_)  // track extra does not exist in AOD
 	 h_ele_PinMnPoutVsChi2_mode     -> Fill(  bestGsfElectron.gsfTrack()->normalizedChi2(), bestGsfElectron.trackMomentumAtVtx().R() - bestGsfElectron.trackMomentumOut().R() );
	h_ele_outerP_mode        -> Fill( bestGsfElectron.trackMomentumOut().R() );
	h_ele_outerPVsEta_mode        -> Fill(bestGsfElectron.eta(),  bestGsfElectron.trackMomentumOut().R() );
	h_ele_outerPt_mode       -> Fill( bestGsfElectron.trackMomentumOut().Rho() );
	h_ele_outerPtVsEta_mode       -> Fill(bestGsfElectron.eta(),  bestGsfElectron.trackMomentumOut().Rho() );
	h_ele_outerPtVsPhi_mode       -> Fill(bestGsfElectron.phi(),  bestGsfElectron.trackMomentumOut().Rho() );
	h_ele_outerPtVsPt_mode       -> Fill(bestGsfElectron.pt(),  bestGsfElectron.trackMomentumOut().Rho() );

	if (!readAOD_) { // track extra does not exist in AOD
          edm::RefToBase<TrajectorySeed> seed = bestGsfElectron.gsfTrack()->extra()->seedRef();
	  ElectronSeedRef elseed=seed.castTo<ElectronSeedRef>();
	  h_ele_seed_dphi2_-> Fill(elseed->dPhi2());
          h_ele_seed_dphi2VsEta_-> Fill(bestGsfElectron.eta(), elseed->dPhi2());
          h_ele_seed_dphi2VsPt_-> Fill(bestGsfElectron.pt(), elseed->dPhi2()) ;
          h_ele_seed_drz2_-> Fill(elseed->dRz2());
          h_ele_seed_drz2VsEta_-> Fill(bestGsfElectron.eta(), elseed->dRz2());
          h_ele_seed_drz2VsPt_-> Fill(bestGsfElectron.pt(), elseed->dRz2());
          h_ele_seed_subdet2_-> Fill(elseed->subDet2());
        }
	// match distributions
	h_ele_EoP    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_EoP_eg    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	if (bestGsfElectron.isEB()) h_ele_EoP_barrel    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EoP_eg_barrel    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	if (bestGsfElectron.isEE()) h_ele_EoP_endcaps    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EoP_eg_endcaps    -> Fill( bestGsfElectron.eSuperClusterOverP() );
	h_ele_EoPVsEta    -> Fill(bestGsfElectron.eta(),  bestGsfElectron.eSuperClusterOverP() );
	h_ele_EoPVsPhi    -> Fill(bestGsfElectron.phi(),  bestGsfElectron.eSuperClusterOverP() );
	h_ele_EoPVsE    -> Fill(bestGsfElectron.caloEnergy(),  bestGsfElectron.eSuperClusterOverP() );
	h_ele_EseedOP    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_EseedOP_eg    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	if (bestGsfElectron.isEB()) h_ele_EseedOP_barrel    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EseedOP_eg_barrel    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	if (bestGsfElectron.isEE()) h_ele_EseedOP_endcaps    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EseedOP_eg_endcaps    -> Fill( bestGsfElectron.eSeedClusterOverP() );
	h_ele_EseedOPVsEta    -> Fill(bestGsfElectron.eta(),  bestGsfElectron.eSeedClusterOverP() );
	h_ele_EseedOPVsPhi    -> Fill(bestGsfElectron.phi(),  bestGsfElectron.eSeedClusterOverP() );
	h_ele_EseedOPVsE    -> Fill(bestGsfElectron.caloEnergy(),  bestGsfElectron.eSeedClusterOverP() );
	h_ele_EoPout -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_EoPout_eg -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.isEB()) h_ele_EoPout_barrel -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EoPout_eg_barrel -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.isEE()) h_ele_EoPout_endcaps -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EoPout_eg_endcaps -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsEta -> Fill( bestGsfElectron.eta(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsPhi -> Fill( bestGsfElectron.phi(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsE -> Fill( bestGsfElectron.caloEnergy(), bestGsfElectron.eSeedClusterOverPout() );
	if (bestGsfElectron.isEE() && classification(bestGsfElectron)==0) {
	if (bestGsfElectron.isEE()) h_ele_EoPout_endcaps_golden -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsEta_golden -> Fill( bestGsfElectron.eta(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsPhi_golden -> Fill( bestGsfElectron.phi(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsE_golden -> Fill( bestGsfElectron.caloEnergy(), bestGsfElectron.eSeedClusterOverPout() );
	} else if (bestGsfElectron.isEE() && classification(bestGsfElectron)==3) {
	if (bestGsfElectron.isEE()) h_ele_EoPout_endcaps_showering -> Fill( bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsEta_showering -> Fill( bestGsfElectron.eta(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsPhi_showering -> Fill( bestGsfElectron.phi(), bestGsfElectron.eSeedClusterOverPout() );
	h_ele_EoPoutVsE_showering -> Fill( bestGsfElectron.caloEnergy(), bestGsfElectron.eSeedClusterOverPout() );
	}
// 	h_ele_EeleOPout -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	if (bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	if (bestGsfElectron.isEB()) h_ele_EeleOPout_barrel -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg_barrel -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	if (bestGsfElectron.isEE()) h_ele_EeleOPout_endcaps -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg_endcaps -> Fill( bestGsfElectron.eEleClusterOverPout() );
// 	h_ele_EeleOPoutVsEta -> Fill( bestGsfElectron.eta(), bestGsfElectron.eEleClusterOverPout() );
// 	h_ele_EeleOPoutVsPhi -> Fill( bestGsfElectron.phi(), bestGsfElectron.eEleClusterOverPout() );
// 	h_ele_EeleOPoutVsE -> Fill( bestGsfElectron.caloEnergy(), bestGsfElectron.eEleClusterOverPout() );
        double EeleOpout_new = (*electronCluster(bestGsfElectron))->energy()/bestGsfElectron.trackMomentumOut().R();
	h_ele_EeleOPout -> Fill( EeleOpout_new );
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg -> Fill( EeleOpout_new );
	if (bestGsfElectron.isEB()) h_ele_EeleOPout_barrel -> Fill( EeleOpout_new );
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg_barrel -> Fill( EeleOpout_new );
	if (bestGsfElectron.isEE()) h_ele_EeleOPout_endcaps -> Fill( EeleOpout_new );
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_EeleOPout_eg_endcaps -> Fill(EeleOpout_new );
	h_ele_EeleOPoutVsEta -> Fill( bestGsfElectron.eta(), EeleOpout_new );
	h_ele_EeleOPoutVsPhi -> Fill( bestGsfElectron.phi(), EeleOpout_new );
	h_ele_EeleOPoutVsE -> Fill( bestGsfElectron.caloEnergy(), EeleOpout_new );
	h_ele_dEtaSc_propVtx -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaSc_propVtx_eg -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEB()) h_ele_dEtaSc_propVtx_barrel -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaSc_propVtx_eg_barrel -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEE())h_ele_dEtaSc_propVtx_endcaps -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaSc_propVtx_eg_endcaps -> Fill(bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	h_ele_dEtaScVsEta_propVtx -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	h_ele_dEtaScVsPhi_propVtx -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	h_ele_dEtaScVsPt_propVtx -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaEtaSuperClusterTrackAtVtx());
	h_ele_dPhiSc_propVtx -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiSc_propVtx_eg -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEB()) h_ele_dPhiSc_propVtx_barrel -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiSc_propVtx_eg_barrel -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEE())h_ele_dPhiSc_propVtx_endcaps -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiSc_propVtx_eg_endcaps -> Fill(bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	h_ele_dPhiScVsEta_propVtx -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	h_ele_dPhiScVsPhi_propVtx -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	h_ele_dPhiScVsPt_propVtx -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaPhiSuperClusterTrackAtVtx());
	h_ele_dEtaCl_propOut -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaCl_propOut_eg -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEB()) h_ele_dEtaCl_propOut_barrel -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaCl_propOut_eg_barrel -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEE()) h_ele_dEtaCl_propOut_endcaps -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaCl_propOut_eg_endcaps -> Fill(bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	h_ele_dEtaClVsEta_propOut -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	h_ele_dEtaClVsPhi_propOut -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	h_ele_dEtaClVsPt_propOut -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaEtaSeedClusterTrackAtCalo());
	h_ele_dPhiCl_propOut -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiCl_propOut_eg -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEB()) h_ele_dPhiCl_propOut_barrel -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiCl_propOut_eg_barrel -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEE()) h_ele_dPhiCl_propOut_endcaps -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiCl_propOut_eg_endcaps -> Fill(bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	h_ele_dPhiClVsEta_propOut -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	h_ele_dPhiClVsPhi_propOut -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	h_ele_dPhiClVsPt_propOut -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaPhiSeedClusterTrackAtCalo());
	h_ele_dEtaEleCl_propOut -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaEleCl_propOut_eg -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	if (bestGsfElectron.isEB()) h_ele_dEtaEleCl_propOut_barrel -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaEleCl_propOut_eg_barrel -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	if (bestGsfElectron.isEE()) h_ele_dEtaEleCl_propOut_endcaps -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dEtaEleCl_propOut_eg_endcaps -> Fill(bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	h_ele_dEtaEleClVsEta_propOut -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	h_ele_dEtaEleClVsPhi_propOut -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	h_ele_dEtaEleClVsPt_propOut -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaEtaEleClusterTrackAtCalo());
	h_ele_dPhiEleCl_propOut -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiEleCl_propOut_eg -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	if (bestGsfElectron.isEB()) h_ele_dPhiEleCl_propOut_barrel -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiEleCl_propOut_eg_barrel -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	if (bestGsfElectron.isEE()) h_ele_dPhiEleCl_propOut_endcaps -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_dPhiEleCl_propOut_eg_endcaps -> Fill(bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	h_ele_dPhiEleClVsEta_propOut -> Fill( bestGsfElectron.eta(),bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	h_ele_dPhiEleClVsPhi_propOut -> Fill(bestGsfElectron.phi(),bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	h_ele_dPhiEleClVsPt_propOut -> Fill(bestGsfElectron.pt(),bestGsfElectron.deltaPhiEleClusterTrackAtCalo());
	h_ele_HoE -> Fill(bestGsfElectron.hadronicOverEm());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_HoE_eg -> Fill(bestGsfElectron.hadronicOverEm());
	if (bestGsfElectron.isEB()) h_ele_HoE_barrel -> Fill(bestGsfElectron.hadronicOverEm());
	if (bestGsfElectron.isEB()&&bestGsfElectron.ecalDrivenSeed()) h_ele_HoE_eg_barrel -> Fill(bestGsfElectron.hadronicOverEm());
	if (bestGsfElectron.isEE()) h_ele_HoE_endcaps -> Fill(bestGsfElectron.hadronicOverEm());
	if (bestGsfElectron.isEE()&&bestGsfElectron.ecalDrivenSeed()) h_ele_HoE_eg_endcaps -> Fill(bestGsfElectron.hadronicOverEm());
	if (!bestGsfElectron.isEBEtaGap() && !bestGsfElectron.isEBPhiGap()&& !bestGsfElectron.isEBEEGap() &&
	    !bestGsfElectron.isEERingGap() && !bestGsfElectron.isEEDeeGap()) h_ele_HoE_fiducial -> Fill(bestGsfElectron.hadronicOverEm());
	h_ele_HoEVsEta -> Fill( bestGsfElectron.eta(),bestGsfElectron.hadronicOverEm());
	h_ele_HoEVsPhi -> Fill(bestGsfElectron.phi(),bestGsfElectron.hadronicOverEm());
	h_ele_HoEVsE -> Fill(bestGsfElectron.caloEnergy(),bestGsfElectron.hadronicOverEm());
	 
//        if (bestGsfElectron.classification() == GsfElectron::GOLDEN && bestGsfElectron.isEB())  histSclEoEtrueGolden_barrel->Fill(sclRef->energy()/mcIter->p());
//        if (bestGsfElectron.classification() == GsfElectron::GOLDEN && bestGsfElectron.isEE())  histSclEoEtrueGolden_endcaps->Fill(sclRef->energy()/mcIter->p());
//        if (bestGsfElectron.classification() == GsfElectron::SHOWERING && bestGsfElectron.isEB())  histSclEoEtrueShowering_barrel->Fill(sclRef->energy()/mcIter->p());
//        if (bestGsfElectron.classification() == GsfElectron::SHOWERING && bestGsfElectron.isEE())  histSclEoEtrueShowering_endcaps->Fill(sclRef->energy()/mcIter->p());
        if (classification(bestGsfElectron)==0 && bestGsfElectron.isEB())  histSclEoEtrueGolden_barrel->Fill(sclRef->energy()/mcIter->p());
        if (classification(bestGsfElectron)==0 && bestGsfElectron.isEE())  histSclEoEtrueGolden_endcaps->Fill(sclRef->energy()/mcIter->p());
        if (classification(bestGsfElectron)==3 && bestGsfElectron.isEB())  histSclEoEtrueShowering_barrel->Fill(sclRef->energy()/mcIter->p());
        if (classification(bestGsfElectron)==3 && bestGsfElectron.isEE())  histSclEoEtrueShowering_endcaps->Fill(sclRef->energy()/mcIter->p());

	//eleClass = eleClass%100; // get rid of barrel/endcap distinction
  h_ele_eta->Fill(std::abs(bestGsfElectron.eta()));
//  if (bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_eta_golden ->Fill(std::abs(bestGsfElectron.eta()));
//  if (bestGsfElectron.classification() == GsfElectron::BIGBREM) h_ele_eta_bbrem ->Fill(std::abs(bestGsfElectron.eta()));
//  //if (bestGsfElectron.classification() == GsfElectron::NARROW) h_ele_eta_narrow ->Fill(std::abs(bestGsfElectron.eta()));
//  if (bestGsfElectron.classification() == GsfElectron::SHOWERING) h_ele_eta_shower ->Fill(std::abs(bestGsfElectron.eta()));

  if (classification(bestGsfElectron)==0) h_ele_eta_golden ->Fill(std::abs(bestGsfElectron.eta()));
  if (classification(bestGsfElectron)==1) h_ele_eta_bbrem ->Fill(std::abs(bestGsfElectron.eta()));
  if (classification(bestGsfElectron)==2) h_ele_eta_narrow ->Fill(std::abs(bestGsfElectron.eta()));
  if (classification(bestGsfElectron)==3) h_ele_eta_shower ->Fill(std::abs(bestGsfElectron.eta()));

	//fbrem
	double fbrem_mean=0.;
	if (!readAOD_) // track extra does not exist in AOD
	 fbrem_mean =  1. - bestGsfElectron.gsfTrack()->outerMomentum().R()/bestGsfElectron.gsfTrack()->innerMomentum().R();
	double fbrem_mode =  bestGsfElectron.fbrem();
	h_ele_fbrem->Fill(fbrem_mode);
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_fbrem_eg->Fill(fbrem_mode);
	h_ele_fbremVsEta_mode->Fill(bestGsfElectron.eta(),fbrem_mode);
	if (!readAOD_) // track extra does not exist in AOD
	 h_ele_fbremVsEta_mean->Fill(bestGsfElectron.eta(),fbrem_mean);

//        if (bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PinVsPoutGolden_mode -> Fill(bestGsfElectron.trackMomentumOut().R(), bestGsfElectron.trackMomentumAtVtx().R());
//        if (bestGsfElectron.classification() == GsfElectron::SHOWERING)
        if (classification(bestGsfElectron)==0) h_ele_PinVsPoutGolden_mode -> Fill(bestGsfElectron.trackMomentumOut().R(), bestGsfElectron.trackMomentumAtVtx().R());
        if (classification(bestGsfElectron)==3)
	 h_ele_PinVsPoutShowering_mode -> Fill(bestGsfElectron.trackMomentumOut().R(), bestGsfElectron.trackMomentumAtVtx().R());
        if (!readAOD_) { // track extra not available in AOD
//          if (bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PinVsPoutGolden_mean -> Fill(bestGsfElectron.gsfTrack()->outerMomentum().R(), bestGsfElectron.gsfTrack()->innerMomentum().R());
//          if (bestGsfElectron.classification() == GsfElectron::SHOWERING)
          if (classification(bestGsfElectron)==0) h_ele_PinVsPoutGolden_mean -> Fill(bestGsfElectron.gsfTrack()->outerMomentum().R(), bestGsfElectron.gsfTrack()->innerMomentum().R());
          if (classification(bestGsfElectron)==3)
	   h_ele_PinVsPoutShowering_mean ->  Fill(bestGsfElectron.gsfTrack()->outerMomentum().R(), bestGsfElectron.gsfTrack()->innerMomentum().R());
        }
//	if (bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PtinVsPtoutGolden_mode -> Fill(bestGsfElectron.trackMomentumOut().Rho(), bestGsfElectron.trackMomentumAtVtx().Rho());
//        if (bestGsfElectron.classification() == GsfElectron::SHOWERING)
	if (classification(bestGsfElectron)==0) h_ele_PtinVsPtoutGolden_mode -> Fill(bestGsfElectron.trackMomentumOut().Rho(), bestGsfElectron.trackMomentumAtVtx().Rho());
        if (classification(bestGsfElectron)==3)
	 h_ele_PtinVsPtoutShowering_mode -> Fill(bestGsfElectron.trackMomentumOut().Rho(), bestGsfElectron.trackMomentumAtVtx().Rho());
        if (!readAOD_) { // track extra not available in AOD
//	  if (bestGsfElectron.classification() == GsfElectron::GOLDEN) h_ele_PtinVsPtoutGolden_mean -> Fill(bestGsfElectron.gsfTrack()->outerMomentum().Rho(), bestGsfElectron.gsfTrack()->innerMomentum().Rho());
//          if (bestGsfElectron.classification() == GsfElectron::SHOWERING)
	  if (classification(bestGsfElectron)==0) h_ele_PtinVsPtoutGolden_mean -> Fill(bestGsfElectron.gsfTrack()->outerMomentum().Rho(), bestGsfElectron.gsfTrack()->innerMomentum().Rho());
          if (classification(bestGsfElectron)==3)
	   h_ele_PtinVsPtoutShowering_mean ->  Fill(bestGsfElectron.gsfTrack()->outerMomentum().Rho(), bestGsfElectron.gsfTrack()->innerMomentum().Rho());
        }

        h_ele_mva->Fill(bestGsfElectron.mva());
        if (bestGsfElectron.ecalDrivenSeed()) h_ele_mva_eg->Fill(bestGsfElectron.mva());
	if (bestGsfElectron.ecalDrivenSeed()) h_ele_provenance->Fill(1.);
	if (bestGsfElectron.trackerDrivenSeed()) h_ele_provenance->Fill(-1.);
	if (bestGsfElectron.trackerDrivenSeed()||bestGsfElectron.ecalDrivenSeed()) h_ele_provenance->Fill(0.);
	if (bestGsfElectron.trackerDrivenSeed()&&!bestGsfElectron.ecalDrivenSeed()) h_ele_provenance->Fill(-2.);
	if (!bestGsfElectron.trackerDrivenSeed()&&bestGsfElectron.ecalDrivenSeed()) h_ele_provenance->Fill(2.);

        h_ele_tkSumPt_dr03->Fill(bestGsfElectron.dr03TkSumPt());
        h_ele_ecalRecHitSumEt_dr03->Fill(bestGsfElectron.dr03EcalRecHitSumEt());
        h_ele_hcalDepth1TowerSumEt_dr03->Fill(bestGsfElectron.dr03HcalDepth1TowerSumEt());
        h_ele_hcalDepth2TowerSumEt_dr03->Fill(bestGsfElectron.dr03HcalDepth2TowerSumEt());
        h_ele_tkSumPt_dr04->Fill(bestGsfElectron.dr04TkSumPt());
        h_ele_ecalRecHitSumEt_dr04->Fill(bestGsfElectron.dr04EcalRecHitSumEt());
        h_ele_hcalDepth1TowerSumEt_dr04->Fill(bestGsfElectron.dr04HcalDepth1TowerSumEt());
        h_ele_hcalDepth2TowerSumEt_dr04->Fill(bestGsfElectron.dr04HcalDepth2TowerSumEt());

      } // gsf electron found

      } // gen electrons in EE
      
    } // mc particle found

    }

  } // loop over mc particle

  h_mcNum->Fill(mcNum);
  h_eleNum->Fill(eleNum);
  h_neleinEE->Fill(neleinEE);

}


