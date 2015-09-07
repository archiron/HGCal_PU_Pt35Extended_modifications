#ifndef HGCALShowerBasedEmIdentification_h
#define HGCALShowerBasedEmIdentification_h

//===================================================================
// Purpose: HGCAL em shower ID variables
// Author: Claude Charlot - LLR- Ecole Polytechnique
// 02/2015
//===================================================================

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoEgamma/Examples/interface/PCAShowerAnalysis.h"
#include "RecoEgamma/Examples/interface/HGCALEmParam.h"
#include "RecoEgamma/Examples/interface/HGCALFitResults.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "TPrincipal.h"

class HGCALShowerBasedEmIdentification
{

  public:

  HGCALShowerBasedEmIdentification(const edm::Event&, const edm::EventSetup&, bool
   withPileup, bool debug=false) ;
  
  bool isEm(const reco::CaloCluster* clu); 
  bool cutStartPosition(const reco::CaloCluster* clu);
  bool cutSigmaetaeta(const reco::CaloCluster* clu);
  bool cutHadOverEm(const reco::CaloCluster* clu);
  bool cutLengthCompatibility(const reco::CaloCluster* clu);

  // shower start position in HGCAL
  GlobalPoint startPosition(const reco::CaloCluster* clu, double cut = 0.);
  // shower entry position in HGCAL
  GlobalPoint entryPosition(const reco::CaloCluster* clu);
  // shower length in cm
  double length(const reco::CaloCluster* clu);
  // shower length compatibility
  double lengthCompatibility(const reco::CaloCluster* clu);
  // layer of shower max 
  int showerMaximum(const reco::CaloCluster* clu);
  // depth at layer i in cm
  double depth(const reco::CaloCluster* clu, int ilayer);
  // energy longitudinal ratioE2530/Etot
  double E2530OverEtot (const reco::CaloCluster* clu, double cut = 0.);
  // energy longitudinal ratioE0110/Etot
  double E0110OverEtot (const reco::CaloCluster* clu, double cut = 0.);
  // shower transverse width in eta
  double sigmaetaeta(const reco::CaloCluster* clu);
  // shower transverse width in phi
  double sigmaphiphi(const reco::CaloCluster* clu);
  // shower transverse with in cm projected on the layer plane
  double sigmartrt(const reco::CaloCluster* clu);
  // shower transverse width in eta per layer
  double sigmaetaeta(const reco::CaloCluster* clu, int ilayer, bool logweight=false);  
  // shower hadronic over em ratio: algo can be "all" (all HCAL layers), "first" (first 4 HCAL layers)
  //  or "last" (last 4 HCAL layers)
  double hadOverEm(const reco::CaloCluster* clu, std::string algo="all");
  // fit and Kolmogorov tests
  HGCALFitResults longitudinalFit(const reco::CaloCluster* clu, bool normalized=false, bool
   bounded=true);
  double longitudinalKolmogorov(const reco::CaloCluster* clu, bool dist=true);
  double transverseKolmogorov(const reco::CaloCluster* clu, bool dist=true);
  
  void setShowerPosition(const GlobalPoint & pos);
  void setShowerDirection(const GlobalVector & dir);

  ~HGCALShowerBasedEmIdentification();
  
private:

  HGCALEmParam hgcalEmParam_;
   
  edm::Handle<HGCRecHitCollection> recHits_;
  edm::Handle<reco::PFClusterCollection> hcalClusters_;
  const HGCalGeometry *geometry_;
  unsigned long long caloGeomCacheId_;  
  
  GlobalPoint showerPos_;
  GlobalVector showerDir_;
  
  bool showerPosIsSet_;
  bool showerDirIsSet_;
   
  bool withPileup_;
  bool debug_;

  // parameters and cut values, to be moved as configurable values
  double mip_;
  double minenergy_;
  double rmax_;
  double hovereConesize_;
  double cutStartPosition_;
  double cutSigmaetaeta_;
  double cutHoverem_;
  double cutLengthCompatibility_;
  
  // longitudinal parametrisation
  double criticalEnergy_;
  double radiationLength_;
  double meant0_;
  double meant1_;
  double meanalpha0_;
  double meanalpha1_;
  double sigmalnt0_;
  double sigmalnt1_;
  double sigmalnalpha0_;
  double sigmalnalpha1_;
  double corrlnalphalnt0_;
  double corrlnalphalnt1_;

};
#endif
