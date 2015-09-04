#ifndef PCAShowerAnalysis_h
#define PCAShowerAnalysis_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "TPrincipal.h"

class PCAShowerAnalysis
{

  public:

  PCAShowerAnalysis( const edm::Event&, const edm::EventSetup&, bool segmented=true, bool logweighting=true, bool debug=false ) ;
  
  void showerParameters( const reco::SuperCluster* , GlobalPoint&, GlobalVector&) ;
  GlobalPoint showerBarycenter( const reco::SuperCluster* ) ;
  GlobalVector showerAxis( const reco::SuperCluster* ) ;
  GlobalVector showerEigenValues( const reco::SuperCluster* ) ;
  GlobalVector showerSigmas( const reco::SuperCluster* ) ;

  void showerParameters( const reco::CaloCluster* , GlobalPoint&, GlobalVector& ) ;
  GlobalPoint showerBarycenter( const reco::CaloCluster* ) ;
  GlobalVector showerAxis( const reco::CaloCluster* ) ;
  GlobalVector showerEigenValues( const reco::CaloCluster* ) ;
  GlobalVector showerSigmas( const reco::CaloCluster* ) ;

  ~PCAShowerAnalysis();
  
private:

  void showerParameters( const reco::CaloClusterPtr , GlobalPoint&, GlobalVector& ) ;
  GlobalPoint showerBarycenter( const reco::CaloClusterPtr ) ;
  GlobalVector showerAxis( const reco::CaloClusterPtr ) ;
  GlobalVector showerEigenValues( const reco::CaloClusterPtr ) ;
  GlobalVector showerSigmas( const reco::CaloClusterPtr ) ;
  
  edm::Handle<HGCRecHitCollection> recHits_;
  const HGCalGeometry *geometry_;
  unsigned long long caloGeomCacheId_;  
  TPrincipal *principal_;
  
  double mip_;
  double entryz_;
  
  bool logweighting_;
  bool segmented_;
  
  bool alreadyfilled_;
  bool debug_;
  
};
#endif
