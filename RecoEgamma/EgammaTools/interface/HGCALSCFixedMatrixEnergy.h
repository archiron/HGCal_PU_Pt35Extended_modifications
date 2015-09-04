#ifndef __HGCALSCFixedMatrixEnergy_H__
#define __HGCALSCFixedMatrixEnergy_H__

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"

// this is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <vector>
#include <memory>

class HGCALSCFixedMatrixEnergy {

 public:

  HGCALSCFixedMatrixEnergy();
  HGCALSCFixedMatrixEnergy(const edm::ParameterSet&);
  void update(const edm::EventSetup&);
// setEvent function is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
  void setEvent(const edm::Event&);

  double getEnergy(const reco::SuperCluster&, int matrixe_size_seed=7, int matrix_size_sub=5) const;
  double getEnergy(const reco::PFCluster *, int matrixe_size) const;

 private:

  unsigned long long caloGeomCacheId_=0;  
  const HGCalGeometry* geometry_;
  //const IdealGeometryRecord* geom_record;  
  
// this is only needed temporarily waiting for the axis to be stored in 
// the CaloCluster
   const reco::PFClusterCollection *clusters_; 
   const reco::PFCluster *getPFCluster(const reco::CaloCluster *) const;
   
};

#endif
