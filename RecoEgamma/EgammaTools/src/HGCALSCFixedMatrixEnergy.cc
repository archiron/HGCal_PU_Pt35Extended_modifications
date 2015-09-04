#include "RecoEgamma/EgammaTools/interface/HGCALSCFixedMatrixEnergy.h"

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include <vector>
#include <memory>

using namespace reco;

HGCALSCFixedMatrixEnergy::HGCALSCFixedMatrixEnergy(const edm::ParameterSet&) 
 : caloGeomCacheId_(0), geometry_(0), clusters_(0)
 { 
 }
 
HGCALSCFixedMatrixEnergy::HGCALSCFixedMatrixEnergy() 
 : caloGeomCacheId_(0), geometry_(0), clusters_(0)
 {  
 }

void HGCALSCFixedMatrixEnergy::update(const edm::EventSetup& es) {

  edm::ESHandle<HGCalGeometry> pGeometry;
  unsigned long long newCaloGeomCacheId= es.get<IdealGeometryRecord>().cacheIdentifier() ;
  if (caloGeomCacheId_!=newCaloGeomCacheId) {
    caloGeomCacheId_ = newCaloGeomCacheId ;
    es.get<IdealGeometryRecord>().get("HGCalEESensitive",pGeometry) ;
    geometry_ = pGeometry.product() ;
  }  
  
}

// this function is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
void HGCALSCFixedMatrixEnergy::setEvent(const edm::Event& iEvent) {

  edm::Handle<PFClusterCollection> pfclustersHandle;
  iEvent.getByLabel( "particleFlowClusterHGCEE", pfclustersHandle );  
  clusters_ = pfclustersHandle.product() ;
  
}
   
// this function is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
const PFCluster* HGCALSCFixedMatrixEnergy::getPFCluster(const CaloCluster *cl) const {

  // find the PFCluster corresponding to the given CaloCluster
  const PFCluster *res = 0;  
  for (PFClusterCollection::const_iterator itcl=clusters_->begin(); itcl!=clusters_->end();
   itcl++) {  
    if (itcl->energy()==cl->energy()) {
      res = &(*itcl);
      break;
    } 
  }  
  //std::cout << "[HGCALSCFixedMatrixEnergy::getPFCluster] found corresponding PFCluster " 
  // << cl->energy() << " " << res->energy() << std::endl;
  return res;

}

double HGCALSCFixedMatrixEnergy::getEnergy(const reco::SuperCluster& sc, int matrix_size_seed, 
 int matrix_size_sub) const {

  double energysc = 0.;
  int matrix_size = matrix_size_sub;
   
  int det = sc.seed()->hitsAndFractions()[0].first.subdetId() ;
  if (det==EcalBarrel) {    
    throw cms::Exception("HGCALSCFixedMatrixEnergy::getEnergy")
      << "Trying to run HGCALSCFixedMatrixEnergy on a non HGCAL cluster!"
      << det << std::endl;    
    return sc.energy();  
  }
   
  // here we are in HGCAL
  //std::cout << "[HGCALSCFixedMatrixEnergy::getEnergy] new supercluster with energy " 
  // << sc.energy() << std::endl;
   
  for (CaloCluster_iterator itcl=sc.clustersBegin(); itcl!=sc.clustersEnd(); itcl++) {

    //const PFCluster *cl = dynamic_cast<const PFCluster*>(&(**itcl));
    const PFCluster *cl = getPFCluster(&(**itcl));
    if (cl==0) continue;
    
    // change the matrix size if the considered cluster is the seed
    if ((*itcl)->energy() == sc.seed()->energy()) matrix_size = matrix_size_seed;
 
    // SC energy from the sum of the basic clusters fixed matrix energies
    //std::cout << "[HGCALSCFixedMatrixEnergy::getEnergy] adding new subcluster with energy " 
    // << (*itcl)->energy() << " and fixed matrix energy " << getEnergy(cl,matrix_size) << std::endl;
    energysc += getEnergy(cl,matrix_size);
    
  }

  //std::cout << "[HGCALSCFixedMatrixEnergy::getEnergy] supercluster fixed matrix energy " << energysc << std::endl;  
  // now apply overal calibration correction for the energy lost outside the window
  // only a global factor is applied
  // value taken from https://indico.cern.ch/event/391727/session/3/contribution/14/material/slides/0.pdf
  double calib = 0.963;
  energysc /= calib;
  //std::cout << "[HGCALSCFixedMatrixEnergy::getEnergy] supercluster fixed matrix energy (calibrated) " << energysc << std::endl;    
  
  return energysc;
     
}

double HGCALSCFixedMatrixEnergy::getEnergy(const PFCluster *cl, int matrix_size) const {
  
  double energycl = 0.;
  
  GlobalPoint pcaShowerPos(cl->position().x(),cl->position().y(),cl->position().z());
  GlobalVector pcaShowerDir(cl->axis().x(),cl->axis().y(),cl->axis().z());

  // Uncalibrated recHits, need here to hard-code the calibration factors
  // Quite bad, an acess should exist to calibrated recHits (or calibration factors)

  // eta correction
  const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
  // new calib as of SLHC21
  double clus_eta = cl->eta();
  double corr = _coef_a*std::abs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*std::abs(clus_eta)+_coef_e));
  double mip = 0.0000551;
  double weight[30] =
   {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};

  for (unsigned int ih=0;ih<cl->hitsAndFractions().size();++ih) {
    const DetId & id_ = (cl->hitsAndFractions())[ih].first ;
    const auto& refhit = cl->recHitFractions()[ih].recHitRef();
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);
      // recompute calibrated rechit energy
      const int layer = hgcid_.layer();
      double scale = mip*corr;
      // energy as matrix sum around cells intercepted by the shower axis
      double lambda = (cellPos.z()-pcaShowerPos.z())/pcaShowerDir.z();
      GlobalPoint interceptPos = pcaShowerPos + lambda*pcaShowerDir;	 
      // hard coded cell size, would be better to retreicve here the cell dimension
      if (std::abs(cellPos.x()-interceptPos.x())<0.95*matrix_size/2. && 
	(std::abs(cellPos.y()-interceptPos.y())<0.95*matrix_size/2.)) {
         energycl += refhit->energy()*weight[layer-1]/scale;
	//std::cout << "adding new cell at " <<  cellPos << " with calibrated energy " << theHit->energy()*weight[layer-1]/scale << std::endl;  
      }
    } 
  }
  
  return energycl; 

}
