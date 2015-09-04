
#ifndef ElectronEnergyCorrector_H
#define ElectronEnergyCorrector_H

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "RecoEgamma/EgammaTools/interface/HGCALSCFixedMatrixEnergy.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

// all this below is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

class EcalClusterFunctionBaseClass ;

class ElectronEnergyCorrector
 {
  public:

    ElectronEnergyCorrector( EcalClusterFunctionBaseClass * crackCorrectionFunction ) ;

    void update(const edm::EventSetup&);
  // the line below is only temporary
  // not needed anymore when axis will be stored inside the CaloCluster   
    void setEvent(edm::Event *);
  
    void classBasedParameterizationEnergy( reco::GsfElectron &, const reco::BeamSpot & bs ) ;
    void classBasedParameterizationUncertainty( reco::GsfElectron & ) ;
    void simpleParameterizationUncertainty( reco::GsfElectron & ) ;
    void hgcalSuperClusterCleanedEnergy( reco::GsfElectron &) ;

  private:

    double fEtaBarrelBad( double scEta ) const ;
    double fEtaBarrelGood( double scEta ) const ;
    double fEtaEndcapBad( double scEta ) const ;
    double fEtaEndcapGood( double scEta ) const ;

    // new corrections (N. Chanon et al.)
    float fEta  (float energy, float eta, int algorithm) const ;
    //float fBrem (float e,  float eta, int algorithm) const ;
    //float fEtEta(float et, float eta, int algorithm) const ;
    float fBremEta(float sigmaPhiSigmaEta, float eta, int algorithm, reco::GsfElectron::Classification cl ) const ;
    float fEt(float et, int algorithm, reco::GsfElectron::Classification cl ) const ;
    float fEnergy(float e, int algorithm, reco::GsfElectron::Classification cl ) const ;

    EcalClusterFunctionBaseClass * crackCorrectionFunction_ ;
  
    // HGCAL
    std::unique_ptr<HGCALSCFixedMatrixEnergy> fixedMatrix_;

// this is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
   const reco::PFClusterCollection *clusters_; 
   const reco::PFCluster *getPFCluster(const reco::CaloCluster *) const;
   edm::Event * event_;

 } ;

#endif




