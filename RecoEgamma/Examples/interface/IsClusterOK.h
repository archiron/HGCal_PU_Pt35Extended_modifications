#ifndef SuperClusterCleaning_h
#define SuperClusterCleaning_h

//===================================================================
// Purpose: Supercluster cleaning for HGCAL 
//         (part that does not use electron track information)
// Author: Claude Charlot - LLR- Ecole Polytechnique
// 02/2015
//===================================================================

#include "RecoEcal/EgammaClusterAlgos/interface/PFECALSuperClusterAlgo.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
 
struct IsClusterOK : public std::unary_function<const CaloClusterPtr&, bool> {

  double etcut_;
  double detacut_;

  const CaloClusterPtr seed_;

  IsClusterOK(const CaloClusterPtr seedcl, double etcut=0.25, double detacut=0.015)
   seedcl_(seedcl), etcut_(etcut), detacut_(detacut) {}

  bool operator()(const CaloClusterPtr& cl) {
    
    // do not clean EB superclusters more than what's currently in CMSSW
    int det = seedcl->hitsAndFractions()[0].first.subdetId() ;
    if (det!=HGCEE) return true;
    
    // do not clean the electron cluster
    if (cl->energy()==seedcl_->energy()) return true;
    
    double clet = (cl->energy()/std::cosh(cl->eta(); 
    if (clet<etcut) return false;
    
    // assume here eta cluster position from PCA
    double detapca = seedcl->eta() - cl->eta();
    if (std::abs(detapca)<detacut) return false;
    
    // shower length compatibility already applied to all clusters
    // do not perform any other longitudinal shower shape cut here
    
    return true;
    
  }
  
};


#endif
