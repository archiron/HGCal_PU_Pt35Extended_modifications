#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronClassification.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

//===================================================================
// Purpose: electron classification for CMS phase II upgrade
// Author: Claude Charlot - LLR- Ecole Polytechnique
// 02/2015
//===================================================================


using namespace reco;

void ElectronClassification::classify( GsfElectron & electron )
 {
 
  if ((!electron.isEB())&&(!electron.isEE()))
   {
    edm::LogWarning("")
      << "ElectronClassification::init(): Undefined electron, eta = "
      << electron.eta() << "!!!!" ;
    electron.setClassification(GsfElectron::UNKNOWN) ;
    return ;
   }

   // barrel gaps
   if ( electron.isEB() && ( electron.isEBEEGap() || electron.isEBEtaGap()) )
   {
    electron.setClassification(GsfElectron::GAP) ;
    return ;
   }
   
  // endcap gaps
   
  int det = electron.superCluster()->seed()->hitsAndFractions()[0].first.subdetId() ; 
  int component = electron.superCluster()->seed()->hitsAndFractions()[0].first.det();  

  if ( electron.isEE() && det == CaloID::DET_ECAL_ENDCAP ) { // not HGCAL   
  
  if ( electron.isEBEEGap() )
   {
    electron.setClassification(GsfElectron::GAP) ;
    return ;
   }
  
  } else if  ( electron.isEE() && component == DetId::Forward && det==HGCEE ) { // HGCAL

  // eta gaps in HGCAL
  if (std::fabs(std::fabs(electron.superCluster()->eta())-1.5) < 0.05 ||
   std::fabs(std::fabs(electron.superCluster()->eta())-3.0) < 0.05) 
   {
    electron.setClassification(GsfElectron::GAP) ;
    return ;
   }
   
  // phi gaps in HGCAL
  double phideg = std::fabs(180.*electron.phi()/CLHEP::pi); 
  if ((phideg>8.&&phideg<12.)||(phideg>28.&&phideg<32.)||(phideg>48.&&phideg<52.)||(phideg>68.&&phideg<72.)||
   (phideg>88.&&phideg<92.)||(phideg>108.&&phideg<112.)||(phideg>128.&&phideg<132.)||(phideg>148.&&phideg<152.)||
   (phideg>168.&&phideg<172.)) 
   {
    electron.setClassification(GsfElectron::GAP) ;
    return ;
   } 
   
  } // end crack/gaps depending on det
   
  double fbremcut = 0.1;
  if (electron.isEE()) fbremcut = 0.2;
  
  if (std::fabs(electron.trackFbrem()) < fbremcut && std::fabs(electron.deltaPhiEleClusterTrackAtCalo()) < 0.008) {
    electron.setClassification(GsfElectron::GOLDEN) ; 
  } else if (std::fabs(electron.trackFbrem()) < 0.4) {
    // beware this is back the narrow class
    // keeep here the BADTRAK enum to avoid modifying the dataformat
    electron.setClassification(GsfElectron::BADTRACK) ;
  } else if (std::fabs(electron.trackFbrem()) > 0.8) {
    electron.setClassification(GsfElectron::BIGBREM) ; 
  } else {
    electron.setClassification(GsfElectron::SHOWERING); 
  }	
  
 }

void ElectronClassification::refineWithPflow( GsfElectron & electron )
 {
  // no badtrack at the moment
  return;
 }

