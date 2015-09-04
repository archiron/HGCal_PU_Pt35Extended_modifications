#include "RecoEgamma/EgammaElectronAlgos/interface/HGCALShowerBasedEmIdentification.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/PCAShowerAnalysis.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "TPrincipal.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"


static Double_t gamma(Double_t *x, Double_t *par) {
  Double_t resu = 0.;
  if (par[1]<=0. || par[2]<=0.) resu = 0.;
  else resu = par[0]*TMath::GammaDist(x[0],par[1],0.,par[2]);
  return resu;
}

static Double_t doublexp(Double_t *x, Double_t *par) {
  // par[4] is costheta where theta is the shower axis
  // transversedetector = transverse / costheta
  Double_t resu = 0.;
  Double_t core = x[0]*par[2]*exp(-par[2]*par[4]*x[0]);
  Double_t tail = x[0]*par[3]*exp(-par[3]*par[4]*x[0]);
  resu = par[0]*(par[1]*core + (1.-par[1])*tail)/(par[4]*par[4]);
  return resu;
}

HGCALShowerBasedEmIdentification::HGCALShowerBasedEmIdentification (const edm::Event&
 iEvent, const edm::EventSetup& iSetup, bool withPileup, bool debug): 
 hgcalEmParam_(0.968, 2.270, 5.36E-3), withPileup_(withPileup), debug_(debug)
{

  iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);
  iEvent.getByLabel(edm::InputTag("particleFlowClusterHGCHEF"),hcalClusters_);
  
  edm::ESHandle<HGCalGeometry> pGeometry;
  unsigned long long newCaloGeomCacheId= iSetup.get<IdealGeometryRecord>().cacheIdentifier() ;
  if (caloGeomCacheId_!=newCaloGeomCacheId) {
    caloGeomCacheId_ = newCaloGeomCacheId ;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",pGeometry) ;
    geometry_ = pGeometry.product() ;
  }  
  
  // initialize showerPos and showerDir
  showerPos_ = GlobalPoint(0.,0.,0.);
  showerDir_ = GlobalVector(0.,0.,0.);
  showerPosIsSet_ = false;
  showerDirIsSet_ = false;
    
  // parameters
  mip_ = 0.0000551;
  minenergy_ = 4.;
  //minenergy_ = 10.;
  rmax_ = 100.; // no transverse limitation for no PU case
  if (withPileup_) rmax_ = 1.5*2.27;
  hovereConesize_ = 0.05;
    
  // cut values, to be moved as configurable parameters
  cutStartPosition_ = 322.5;
  cutSigmaetaeta_ = 0.0055;
  if (withPileup_) cutSigmaetaeta_ = 0.00480;
  cutHoverem_ = 0.003;
  if (withPileup_) cutHoverem_ = 0.065;
  cutLengthCompatibility_ = 4.0;
  
  std::cout << "*** HGCAL ShowerBased EmIdentification ***" << std::endl;
  std::cout << "- max transverse radius: " << rmax_ << std::endl;
  std::cout << "- hovere cone size: " << hovereConesize_ << std::endl;
  std::cout << "- cut sigmaetaeta: " << cutSigmaetaeta_ << std::endl;
  std::cout << "- cut hoverem: " << cutHoverem_ << std::endl;
  std::cout << "- cut start position: " << cutStartPosition_ << std::endl;
  std::cout << "- cut length compatibility (in sigmas): " << cutLengthCompatibility_ << std::endl;
  
}

HGCALShowerBasedEmIdentification::~HGCALShowerBasedEmIdentification ()
{
//  delete pcaShowerAnalysis_;
}

void HGCALShowerBasedEmIdentification::setShowerPosition(const GlobalPoint &pos)
{
  showerPos_ = pos;
  showerPosIsSet_ = true;
}

void HGCALShowerBasedEmIdentification::setShowerDirection(const GlobalVector &dir)
{
  showerDir_ = dir;
  showerDirIsSet_ = true;
}

bool HGCALShowerBasedEmIdentification::isEm(const reco::CaloCluster* clu)
{
  return (cutStartPosition(clu) &&  cutSigmaetaeta(clu) && cutHadOverEm(clu));
} 

GlobalPoint HGCALShowerBasedEmIdentification::startPosition(const reco::CaloCluster*
 clu, double cut)
{

  // the position in z is the first z position of the cells within the cluster 
  // in x and y take the max energy in this first layer

  // first determine the energy in each layer
  double energy[30];
  for (int i=0; i<30; i++) energy[i]=0.;
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      int layer = hgcid_.layer();
      energy[layer-1] += theHit->energy();
    }
  }
  
  GlobalPoint firstPos;
  double zmin = 10000., maxfirstenergy=0.; 
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);
      // add a cut on energy for the pileup
      if (energy[hgcid_.layer()-1]>cut*mip_) {
	if (fabs(cellPos.z())<zmin || ((fabs(cellPos.z())==zmin)&&(theHit->energy()>maxfirstenergy))) {
	  firstPos = cellPos;
	  zmin = fabs(cellPos.z());
	  maxfirstenergy = theHit->energy();
	}
      }
    }
  }

  // finally refine firstPos x and y using the meaured direction 
  if (!showerPosIsSet_ || !showerDirIsSet_) return firstPos;
 
  double lambda = (firstPos-showerPos_).z()/showerDir_.z();
  GlobalPoint extraPos = showerPos_ + lambda*showerDir_;	 
  firstPos = extraPos;
  if (debug_) {
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] showerPos " << showerPos_ << std::endl;
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] showerDir " << showerDir_ << std::endl;
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] firstPos " << firstPos << std::endl;
  }
 
  return firstPos;

}

GlobalPoint HGCALShowerBasedEmIdentification::entryPosition(const reco::CaloCluster* clu)
{
 
  // first check that showerPos has been set
  if (!showerPosIsSet_ || !showerDirIsSet_) {
   std::cout << "[HGCALShowerBasedEmIdentification::entryPosition] error, showwer position not set " << std::endl;
   return GlobalPoint(0.,0.,0.);
  }

  double lambda = (320.-showerPos_.z())/showerDir_.z();
  if (showerPos_.z()<0.) lambda = (-320.-showerPos_.z())/showerDir_.z();
  GlobalPoint entryPos = showerPos_ + lambda*showerDir_;	 
  if (debug_) {
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] showerPos " << showerPos_ << std::endl;
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] showerDir " << showerDir_ << std::endl;
   std::cout << "[HGCALShowerBasedEmIdentification::startPosition] firstPos " << entryPos << std::endl;
  }
 
  return entryPos;

}

int HGCALShowerBasedEmIdentification::showerMaximum(const reco::CaloCluster* clu)
{

  // the layer which contains the cell of maximimum energy within the cluster 

  int layermax = 10000., maxfirstenergy=0.; 
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theSeedHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      if (theSeedHit->energy()>maxfirstenergy) {
	layermax = HGCEEDetId(id_).layer();
	maxfirstenergy = theSeedHit->energy();
      }
    }
  }

  return layermax;

}

double HGCALShowerBasedEmIdentification::depth(const reco::CaloCluster* clu, int ilayer)
{

  // compute the depth from the shower start at the layer i extrapolating 
  // along the shower axis
  GlobalPoint maxPos;
  double zlayer=0.;
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId(id_) ;
      if (hgcid_.layer() == ilayer) {
	zlayer = geometry_->getPosition(hgcid_).z();
	break;
      }
    }
  }
  
  double lambda = (zlayer-showerPos_.z())/showerDir_.z();
  GlobalPoint extraPos = showerPos_ + lambda*showerDir_;	 

  return (extraPos - startPosition(clu)).mag();

}

double HGCALShowerBasedEmIdentification::E2530OverEtot(const reco::CaloCluster*
 clu, double cut)
{

  // first determine the energy in each layer
  double energy[30]; double etot=0.; // not calibrated
  for (int i=0; i<30; i++) energy[i]=0.;
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      int layer = hgcid_.layer();
      // add a cut on energy for the pileup
      if (theHit->energy()>cut*mip_) {
        energy[layer-1] += theHit->energy();
        etot += theHit->energy();
      }	
    }
  }

  // then makes sums
  double e2530 = 0.;
  for (int layer=24; layer<30; layer++) e2530 += energy[layer];
  
  return e2530/etot;
  
}

double HGCALShowerBasedEmIdentification::E0110OverEtot(const reco::CaloCluster*
 clu, double cut)
{

  // first determine the energy in each layer
  double energy[30]; double etot=0.; // not calibrated
  for (int i=0; i<30; i++) energy[i]=0.;
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      int layer = hgcid_.layer();
      // add a cut on energy for the pileup
      if (theHit->energy()>cut*mip_) {
        energy[layer-1] += theHit->energy();
        etot += theHit->energy();
      }	
    }
  }

  // then makes sums
  double e0110 = 0.;
  for (int layer=0; layer<10; layer++) e0110 += energy[layer];
  
  return e0110/etot;
  
}

double HGCALShowerBasedEmIdentification::length(const reco::CaloCluster* clu)
{

  // first check that showerPos has been set
  if (!showerPosIsSet_) {
   std::cout << "[HGCALShowerBasedEmIdentification::length] error, showwer position not set " << std::endl;
   std::cout << "[HGCALShowerBasedEmIdentification::length] error, please invoke setShowerPos and setShowerDir before invoking this function " << std::endl;
   return 0.;
  }

  // shower length	 
  double length = (showerPos_ - startPosition(clu)).mag();

  return length;
  
}

double HGCALShowerBasedEmIdentification::lengthCompatibility(const reco::CaloCluster* clu)
{

  double lny = clu->energy()/hgcalEmParam_.getCriticalEnergy()>1. ? std::log(clu->energy()/hgcalEmParam_.getCriticalEnergy()) : 0.;

  // inject here parametrization results
  double meantmax = hgcalEmParam_.meanT(lny);
  double meanalpha = hgcalEmParam_.meanAlpha(lny);
  double sigmalntmax = hgcalEmParam_.sigmaLnT(lny);
  double sigmalnalpha = hgcalEmParam_.sigmaLnAlpha(lny);
  double corrlnalphalntmax = hgcalEmParam_.correlationAlphaT(lny);
  
  double invbeta = meantmax/(meanalpha-1.);
  double predictedLength = meanalpha*invbeta;
  predictedLength *= hgcalEmParam_.getRadiationLength();
  
  double sigmaalpha = meanalpha*sigmalnalpha;
  if (sigmaalpha<0.) sigmaalpha = 1.;
  double sigmatmax = meantmax*sigmalntmax;
  if (sigmatmax<0.) sigmatmax = 1.;
  
  double predictedSigma = sigmalnalpha*sigmalnalpha/((meanalpha-1.)*(meanalpha-1.));
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] sigma/m squared 1st term " << predictedSigma << std::endl;
  predictedSigma += sigmalntmax*sigmalntmax;
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] sigma/m squared with 2nd term " << predictedSigma << std::endl;
  predictedSigma -= 2*sigmalnalpha*sigmalntmax*corrlnalphalntmax/(meanalpha-1.);
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] sigma/m squared with correlation " << predictedSigma << std::endl;
  predictedSigma = predictedLength*sqrt(predictedSigma);
  
  double showerLength = length(clu);
  double lengthCompatibility = (showerLength-predictedLength)/predictedSigma;
  
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] predictedLength-showerLength " << predictedLength-showerLength << std::endl;
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] sigmalntmax " << sigmalntmax << std::endl;
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] sigmalnalpha " << sigmalnalpha << std::endl;
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] predictedSigma " << predictedSigma << std::endl;
  if (debug_) std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] lengthCompatibility " << lengthCompatibility << std::endl;
  
//  // ancien calcul  
//  double sigma = predictedLength / (-2.506+1.245*lny);
//  double oldlengthCompatibility = fabs(predictedLength-length)*radiationLength_/sigma;
//  std::cout << "[HGCALShowerBasedEmIdentification::lengthCompatibility] previous lengthCompatibility " << oldlengthCompatibility << std::endl;
 
  return lengthCompatibility;
  
}

double HGCALShowerBasedEmIdentification::sigmaetaeta(const reco::CaloCluster* clu)
{

  double sigmaetaeta=0., sumnrj=0.;
  GlobalPoint firstPos = startPosition(clu);
    
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      GlobalVector radius, longitudinal, transverse;
      radius = cellPos - firstPos;
      // distances in local coordinates
      longitudinal =  (radius.dot(showerDir_))*showerDir_.unit()/showerDir_.mag();
      transverse = radius - longitudinal;
      // apply energy cut cut
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
	// simple transversal cut, later can refine as function of depth
	if (!withPileup_ || transverse.mag() < rmax_) {
	  sigmaetaeta += (cellPos.eta()-showerPos_.eta())*(cellPos.eta()-showerPos_.eta()) * theHit->energy();	     
	  sumnrj += theHit->energy();
	}
      }
    }
  }

  sigmaetaeta /= sumnrj;
  sigmaetaeta = sqrt(sigmaetaeta);

  // now correct the eta dependency
  double feta; double feta_0; 
  if (!withPileup_) {
  // NoPU corrections
    feta = 0.00964148 - 0.0107843*fabs(clu->eta()) + 0.00495703*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00964148 - 0.01078431*1.5 + 0.00495703*1.5*1.5;
  } else {
    feta = 0.0146784 - 0.0176412*fabs(clu->eta()) + 0.00722875*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.014678 - 0.0176412*1.5 + 0.00722875*1.5*1.5;
  }
  sigmaetaeta *= feta_0 / feta ;

  return sigmaetaeta;

}

double HGCALShowerBasedEmIdentification::sigmaphiphi(const reco::CaloCluster* clu)
{

  double sigmaphiphi=0., sumnrj=0.;
  GlobalPoint firstPos = startPosition(clu);
    
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      GlobalVector radius, longitudinal, transverse;
      radius = cellPos - firstPos;
      // distances in local coordinates
      longitudinal =  (radius.dot(showerDir_))*showerDir_.unit()/showerDir_.mag();
      transverse = radius - longitudinal;
      // apply energy cut cut
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
	// simple transversal cut, later can refine as function of depth
	if (!withPileup_ || transverse.mag() < rmax_) {
	  double dphi = (cellPos.phi()-showerPos_.phi());
	  if (std::abs(dphi)>CLHEP::pi) dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
 	  sigmaphiphi += dphi * dphi * theHit->energy();
	  sumnrj += theHit->energy();
	}
      }
    }
  }

  sigmaphiphi /= sumnrj;
  sigmaphiphi = sqrt(sigmaphiphi);

  // now correct the eta dependency
  double fphi; double fphi_0; 
  if (!withPileup_) {
  // NoPU corrections
    fphi = -1.33123e-04 - 1.47429e-03*fabs(clu->eta()) + 2.17692e-03*fabs(clu->eta())*fabs(clu->eta());
    fphi_0 = -1.33123e-04 - 1.47429e-03*1.5 + 2.17692e-03*1.5*1.5;
  } else {
    fphi = 0.01351 - 0.014881*fabs(clu->eta()) + 0.00665377*fabs(clu->eta())*fabs(clu->eta());
    fphi_0 = 0.01351 - 0.014881*1.5 + 0.00665377*1.5*1.5;
  }
  sigmaphiphi *= fphi_0 / fphi ;

  return sigmaphiphi;

}

double HGCALShowerBasedEmIdentification::sigmartrt(const reco::CaloCluster* clu)
{

  double sigmartrt=0., sumnrj=0.;
  GlobalPoint firstPos = startPosition(clu);
    
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      GlobalPoint axisPos;
      GlobalVector radius, longitudinal, transverse;
      radius = cellPos - firstPos;
      // distances in local coordinates
      longitudinal =  (radius.dot(showerDir_))*showerDir_.unit()/showerDir_.mag();
      transverse = radius - longitudinal;
      // apply energy cut cut
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
	// simple transversal cut, later can refine as function of depth
	if (!withPileup_ || transverse.mag() < rmax_) {
	  double lambdatoaverage = (cellPos.z()-showerPos_.z())/showerDir_.z();
	  axisPos = showerPos_ + lambdatoaverage*showerDir_;	 
          sigmartrt += (cellPos.perp()-axisPos.perp()) * (cellPos.perp()-axisPos.perp()) * theHit->energy();	     				
	  sumnrj += theHit->energy();
	}
      }
    }
  }

  sigmartrt /= sumnrj;
  sigmartrt = sqrt(sigmartrt);

//   // now correct the eta dependency
//   double feta; double feta_0; 
//   // correction derived with 140PU
//   if (fabs(clu->eta())<2.6)
//    feta = 4.13856 - 2.52796*fabs(clu->eta()) + 0.487409*fabs(clu->eta())*fabs(clu->eta());
//   else feta = 4.13856 - 2.52796*2.6 + 0.487409*2.6*2.6;
//   feta_0 = 4.13856 - 2.52796*1.5 + 0.487409*1.5*1.5;
//  
//  sigmartrt *= feta_0 / feta ;

  return sigmartrt;

}

double HGCALShowerBasedEmIdentification::sigmaetaeta(const reco::CaloCluster* clu, 
 int ilayer, bool logweight)
{

  double sigmaetaeta=0., sumnrj=0., sumw=0.;
  TPrincipal principal(2,"D"); 
  GlobalPoint position, cellPos;
  
  double variables[2] = {0.,0.};

  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      if (HGCEEDetId(id_).layer()==ilayer) {
      cellPos = geometry_->getPosition(HGCEEDetId(id_));
      variables[0] = cellPos.x(); variables[1] = cellPos.y(); 
      sumnrj += theHit->energy();
      //// energy weighting
      //for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
      //} else {
      // a log-weighting, energy not in fraction of total
      double w0 = -log(4.); // threshold, reduced here
      double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                   //  for the highest hit of ~0.1 GeV
      int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
      if (nhit<0) nhit=1;
      for (int i=0; i<nhit; i++) principal.AddRow(variables);	
      }	     
    }
  }
    
  //std::cout << "[HGCALShowerBasedEmIdentification::sigmaetaeta], layer " << ilayer
  //<< " sumnrj " << sumnrj << std::endl;

  principal.MakePrincipals();  
  position = GlobalPoint((*principal.GetMeanValues())[0],(*principal.GetMeanValues())[1],cellPos.z());	 

  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      if (HGCEEDetId(id_).layer()==ilayer) {
      GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));	     
      double weight = theHit->energy();
      if (logweight) {
        double w0 = 2.; 
        weight = std::max(0.,w0 + log(theHit->energy()/sumnrj));
      }	
      sigmaetaeta += (cellPos.eta()-position.eta())*(cellPos.eta()-position.eta()) * weight;
      sumw += weight;
      }
    }
  }

  //std::cout << "[HGCALShowerBasedEmIdentification::sigmaetaeta], layer " << ilayer
  //<< " position " << position << " sigmaeta " << sigmaetaeta << std::endl;

  if (sumw==0.) return 0.;

  sigmaetaeta /= sumw;
  sigmaetaeta = sqrt(sigmaetaeta);
  //std::cout << "[HGCALShowerBasedEmIdentification::sigmaetaeta], layer " << ilayer
  //<< " position " << position << " sigmaeta " << sigmaetaeta << std::endl;
   
  // now correct the eta dependency
  double feta=1.; double feta_0=1.; 

  if (!withPileup_) {
  // NoPU corrections
  if (ilayer==1) {
    feta = 0.0101201 - 0.0102206*fabs(clu->eta()) + 0.00362097*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0101201 - 0.0102206*1.5 + 0.00362097*1.5*1.5;
  } else if (ilayer==2) {
    feta = 0.00766728 - 0.00791464*fabs(clu->eta()) + 0.00303245*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00766728 - 0.00791464*1.5 + 0.00303245*1.5*1.5;
  } else if (ilayer==3) {
    feta = 0.00634266 - 0.00656176*fabs(clu->eta()) + 0.00261383*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00634266 - 0.00656176*1.5 + 0.00261383*1.5*1.5;
  } else if (ilayer==4) {
    feta = 0.0078056 - 0.00834319*fabs(clu->eta()) + 0.00335996*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0078056 - 0.00834319*1.5 + 0.00335996*1.5*1.5;
  } else if (ilayer==5) {
    feta = 0.00749733 - 0.007857*fabs(clu->eta()) + 0.00311549*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00749733 - 0.007857*1.5 + 0.00311549*1.5*1.5;
  } else if (ilayer==6) {
    feta =  0.00913178 - 0.00956428*fabs(clu->eta()) + 0.00381994*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00913178 - 0.00956428*1.5 + 0.00381994*1.5*1.5;
  } else if (ilayer==7) {
    feta =  0.00699386 - 0.0074963*fabs(clu->eta()) + 0.0032499*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00699386 - 0.0074963*1.5 + 0.0032499*1.5*1.5;
  } else if (ilayer==8) {
    feta =  0.00816775 - 0.00871664*fabs(clu->eta()) + 0.00382981*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00816775 - 0.00871664*1.5 + 0.00382981*1.5*1.5;
  } else if (ilayer==9) {
    feta =  0.00617746 - 0.00665314*fabs(clu->eta()) + 0.00323033*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00617746 - 0.00665314*1.5 + 0.00323033*1.5*1.5;
  } else if (ilayer==10) {
    feta =  0.00872935 - 0.00928504*fabs(clu->eta()) + 0.00416514*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00872935 - 0.00928504*1.5 + 0.00416514*1.5*1.5;
  } else if (ilayer==11) {
    feta =  0.00766227 - 0.00822116*fabs(clu->eta()) + 0.00383486*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00766227 - 0.00822116*1.5 + 0.00383486*1.5*1.5;
  } else if (ilayer==12) {
    feta =  0.0088128 - 0.00946105 *fabs(clu->eta()) + 0.00437748*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.0088128 -0.00946105 *1.5 + 0.00437748*1.5*1.5;
  } else if (ilayer==13) {
    feta =  0.00759962 - 0.00827032 *fabs(clu->eta()) + 0.00403462*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00759962 - 0.00827032 *1.5 + 0.00403462*1.5*1.5;
  } else if (ilayer==14) {
    feta = 0.00975987 - 0.0103078 *fabs(clu->eta()) + 0.00477958*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  0.00975987 - 0.0103078 *1.5 + 0.00477958*1.5*1.5;
  } else if (ilayer==15) {
    feta  = 0.00805003 - 0.00864848 *fabs(clu->eta()) + 0.00437212*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00805003 - 0.00864848 *1.5 + 0.00437212*1.5*1.5;
  } else if (ilayer==16) {
    feta  = 0.0114562 - 0.0121809 *fabs(clu->eta()) + 0.00553499*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0114562 - 0.0121809 *1.5 + 0.00553499*1.5*1.5;
  } else if (ilayer==17) {
    feta  = 0.00902615 - 0.00966597 *fabs(clu->eta()) + 0.00486833*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00902615 - 0.00966597 *1.5 + 0.00486833*1.5*1.5;
  } else if (ilayer==18) {
    feta  = 0.0132086 - 0.0140273 *fabs(clu->eta()) + 0.00627827*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0132086 - 0.0140273 *1.5 + 0.00627827*1.5*1.5;
  } else if (ilayer==19) {
    feta  = 0.00839401 - 0.00930494 *fabs(clu->eta()) + 0.00511968*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00839401 - 0.00930494 *1.5 + 0.00511968*1.5*1.5;
  } else if (ilayer==20) {
    feta  = 0.009421 - 0.0103944 *fabs(clu->eta()) + 0.00568878*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.009421 - 0.0103944 *1.5 + 0.00568878*1.5*1.5;
  } else if (ilayer==21) {
    feta  = 0.00607798 - 0.0081057 *fabs(clu->eta()) + 0.00530424*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00607798 - 0.0081057 *1.5 + 0.00530424*1.5*1.5;
  } else if (ilayer==22) {
    feta  = 0.00243842 - 0.0048886 *fabs(clu->eta()) + 0.00482695*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00243842 - 0.0048886 *1.5 + 0.00482695*1.5*1.5;
  } else if (ilayer==23) {
    feta  = 0.00233566 - 0.00745343 *fabs(clu->eta()) + 0.00604013*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00233566 - 0.00745343 *1.5 + 0.00604013*1.5*1.5;
  } else if (ilayer==24) {
    feta  = 0.00233566 - 0.00745343 *fabs(clu->eta()) + 0.00604013*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00233566 - 0.00745343 *1.5 + 0.00604013*1.5*1.5;
  } else if (ilayer==25) {
    feta  = 0.0124781 - 0.0180787 *fabs(clu->eta()) + 0.00848182*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0124781 - 0.0180787 *1.5 + 0.00848182*1.5*1.5;
  } else if (ilayer==26) {
    feta  = 0.0178317 - 0.0247529 *fabs(clu->eta()) + 0.010354*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0178317 - 0.0247529 *1.5 + 0.010354*1.5*1.5;
  } else if (ilayer==27) {
    feta  = 0.0236844 - 0.0301237 *fabs(clu->eta()) + 0.0111237*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0236844 - 0.0301237 *1.5 + 0.0111237*1.5*1.5;
  } else if (ilayer==28) {
    feta  = 0.0114311 - 0.0179117 *fabs(clu->eta()) + 0.00802573*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0114311 - 0.0179117 *1.5 + 0.00802573*1.5*1.5;
  } else if (ilayer==29) {
    feta  = 0.0158825 - 0.0203313 *fabs(clu->eta()) + 0.00774348*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0158825 - 0.0203313 *1.5 + 0.00774348*1.5*1.5;
  } else if (ilayer==30) {
    feta  = 0.0232539 - 0.0272789 *fabs(clu->eta()) + 0.00920109*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.0232539 - 0.0272789 *1.5 + 0.00920109*1.5*1.5;
  }
  } else {
  if (ilayer==1) {
    feta = 1.33362e-02 - 1.48793e-02*fabs(clu->eta()) + 5.44801e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.33362e-02 - 1.48793e-02*1.5 + 5.44801e-03*1.5*1.5;
  } else if (ilayer==2) {
    feta = 2.54405e-02 - 2.83450e-02*fabs(clu->eta()) + 9.02912e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 2.54405e-02 - 2.83450e-02*1.5 + 9.02912e-03*1.5*1.5;
  } else if (ilayer==3) {
    feta = 2.20685e-02 - 2.53305e-02*fabs(clu->eta()) + 8.32861e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 2.20685e-02 - 2.53305e-02*1.5 + 8.32861e-03*1.5*1.5;
  } else if (ilayer==4) {
    feta = 2.19241e-02 - 2.50043e-02*fabs(clu->eta()) + 8.39974e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 2.19241e-02 - 2.50043e-02*1.5 + 8.39974e-03*1.5*1.5;
  } else if (ilayer==5) {
    feta =  1.89229e-02 - 2.15890e-02*fabs(clu->eta()) + 7.56937e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.89229e-02 - 2.15890e-02*1.5 + 7.56937e-03*1.5*1.5;
  } else if (ilayer==6) {
    feta =  2.10086e-02 - 2.40040e-02*fabs(clu->eta()) + 8.01320e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  2.10086e-02 - 2.40040e-02*1.5 + 8.01320e-03*1.5*1.5;
  } else if (ilayer==7) {
    feta =  1.83148e-02 - 2.07314e-02*fabs(clu->eta()) + 7.19385e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.83148e-02 - 2.07314e-02*1.5 + 7.19385e-03*1.5*1.5;
  } else if (ilayer==8) {
    feta =  1.70924e-02 - 1.93284e-02*fabs(clu->eta()) + 7.03540e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.70924e-02 - 1.93284e-02*1.5 + 7.03540e-03*1.5*1.5;
  } else if (ilayer==9) {
    feta =  1.66959e-02 - 1.88356e-02*fabs(clu->eta()) + 6.75873e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.66959e-02 - 1.88356e-02*1.5 + 6.75873e-03*1.5*1.5;
  } else if (ilayer==10) {
    feta =  1.24162e-02 - 1.40265e-02*fabs(clu->eta()) + 5.69077e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.24162e-02 - 1.40265e-02*1.5 + 5.69077e-03*1.5*1.5;
  } else if (ilayer==11) {
    feta =  1.17059e-02 - 1.32282e-02*fabs(clu->eta()) + 5.40480e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.17059e-02 - 1.32282e-02*1.5 + 5.40480e-03*1.5*1.5;
  } else if (ilayer==12) {
    feta =  1.16689e-02 - 1.30371e-02 *fabs(clu->eta()) + 5.49795e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.16689e-02 -1.30371e-02 *1.5 + 5.49795e-03*1.5*1.5;
  } else if (ilayer==13) {
    feta =  1.02646e-02 - 1.17211e-02 *fabs(clu->eta()) + 5.14711e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.02646e-02 - 1.17211e-02 *1.5 + 5.14711e-03*1.5*1.5;
  } else if (ilayer==14) {
    feta = 1.26116e-02 - 1.36569e-02 *fabs(clu->eta()) + 5.75172e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 =  1.26116e-02 - 1.36569e-02 *1.5 + 5.75172e-03*1.5*1.5;
  } else if (ilayer==15) {
    feta  = 8.28188e-03 - 9.39396e-03 *fabs(clu->eta()) + 4.68145e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 8.28188e-03 - 9.39396e-03 *1.5 + 4.68145e-03*1.5*1.5;
  } else if (ilayer==16) {
    feta  = 1.14035e-02 - 1.23078e-02 *fabs(clu->eta()) + 5.60658e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.14035e-02 - 1.23078e-02 *1.5 + 5.60658e-03*1.5*1.5;
  } else if (ilayer==17) {
    feta  = 1.15703e-02 - 1.26016e-02 *fabs(clu->eta()) + 5.73669e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.15703e-02 - 1.26016e-02 *1.5 + 5.73669e-03*1.5*1.5;
  } else if (ilayer==18) {
    feta  = 1.16629e-02 - 1.22541e-02 *fabs(clu->eta()) + 5.89914e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.16629e-02 - 1.22541e-02 *1.5 + 5.89914e-03*1.5*1.5;
  } else if (ilayer==19) {
    feta  = 1.29627e-02 - 1.42037e-02 *fabs(clu->eta()) + 6.56636e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.29627e-02 - 1.42037e-02 *1.5 + 6.56636e-03*1.5*1.5;
  } else if (ilayer==20) {
    feta  = 7.21102e-03 - 7.96553e-03 *fabs(clu->eta()) + 5.29562e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 7.21102e-03 - 7.96553e-03 *1.5 + 5.29562e-03*1.5*1.5;
  } else if (ilayer==21) {
    feta  = 0.00607798 - 0.0081057 *fabs(clu->eta()) + 0.00530424*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 0.00607798 - 0.0081057 *1.5 + 0.00530424*1.5*1.5;
  } else if (ilayer==22) {
    feta  = -3.58073e-03 + 2.20787e-03 *fabs(clu->eta()) + 3.27484e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -3.58073e-03 + 2.20787e-03 *1.5 + 3.27484e-03*1.5*1.5;
  } else if (ilayer==23) {
    feta  = -4.88647e-03 + 2.74604e-03 *fabs(clu->eta()) + 3.39906e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -4.88647e-03 + 2.74604e-03 *1.5 + 3.39906e-03*1.5*1.5;
  } else if (ilayer==24) {
    feta  = -1.02614e-02 + 7.83279e-03 *fabs(clu->eta()) + 2.49495e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -1.02614e-02 + 7.83279e-03 *1.5 + 2.49495e-03*1.5*1.5;
  } else if (ilayer==25) {
    feta  = -1.30804e-02 + 9.54999e-03 *fabs(clu->eta()) + 2.25360e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -1.30804e-02 + 9.54999e-03 *1.5 + 2.25360e-03*1.5*1.5;
  } else if (ilayer==26) {
    feta  = -1.63714e-02 + 1.34650e-02 *fabs(clu->eta()) + 1.37441e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -1.63714e-02 + 1.34650e-02 *1.5 + 1.37441e-03*1.5*1.5;
  } else if (ilayer==27) {
    feta  = -8.82053e-03 + 4.55760e-03 *fabs(clu->eta()) + 3.65359e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -8.82053e-03 + 4.55760e-03 *1.5 + 3.65359e-03*1.5*1.5;
  } else if (ilayer==28) {
    feta  = 3.03863e-03 - 7.43808e-03 *fabs(clu->eta()) +  6.53450e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 3.03863e-03 - 7.43808e-03 *1.5 +  6.53450e-03*1.5*1.5;
  } else if (ilayer==29) {
    feta  = 1.09964e-02 - 1.53346e-02 *fabs(clu->eta()) + 8.16445e-03*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = 1.09964e-02 - 1.53346e-02 *1.5 + 8.16445e-03*1.5*1.5;
  } else if (ilayer==30) {
    feta  = -2.10798e-02 + 1.59415e-02 *fabs(clu->eta()) + 7.63047e-04*fabs(clu->eta())*fabs(clu->eta());
    feta_0 = -2.10798e-02 + 1.59415e-02 *1.5 + 7.63047e-04*1.5*1.5;
  }
  }
  sigmaetaeta *= feta_0 / feta ;

  return sigmaetaeta;

}

double HGCALShowerBasedEmIdentification::hadOverEm(const reco::CaloCluster* clu, std::string algo)
{

  // check algo, algo can be :
  // - all for summing all HCAL clusters in the cone
  // - first for summing only HCAL clusters starting in the first 4 layers
  // - last for summing only HCAL clusters starting in the last 4 layers
  
  if (algo!="all" && algo!="first" && algo!="last") {
   std::cout << "[HGCALShowerBasedEmIdentification::hadOverEm] error, wrong algo name " << std::endl;
   return 0.;
  }
   
  // first check that showerPos has been set
  if (!showerPosIsSet_) {
   std::cout << "[HGCALShowerBasedEmIdentification::hadOverEm] error, showwer position not set " << std::endl;
   return 0.;
  }

  // H/E 
  math::XYZVector vectorSC(showerPos_.x(),showerPos_.y(),showerPos_.z());
  std::vector<reco::PFCluster>::const_iterator trItr = hcalClusters_->begin();
  std::vector<reco::PFCluster>::const_iterator trItrEnd = hcalClusters_->end();
  double hoverem=0.;  
  for( ;  trItr != trItrEnd ; ++trItr){
    math::XYZVector vectorHgcalHFECluster(trItr->position().x(),trItr->position().y(),trItr->position().z());
    double dR = ROOT::Math::VectorUtil::DeltaR(vectorSC,vectorHgcalHFECluster);
    if (dR<hovereConesize_) {
      const DetId & id_ = (trItr->hitsAndFractions())[0].first ;
      HGCHEDetId HGCHEid_(id_);
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCHEF) {
        // HCAL cluster starts in the first four layers
	if (algo=="all") hoverem += trItr->energy();
	else if (algo=="first" && HGCHEid_.layer()<5) hoverem += trItr->energy();
	else if (algo=="last" && HGCHEid_.layer()>8) hoverem += trItr->energy();
      }		
    }  
  }  
  
  // here use the seed cluster
  //hoverem /= scle;
  hoverem /= clu->energy();

  return hoverem;
  
}
  
HGCALFitResults HGCALShowerBasedEmIdentification::longitudinalFit(const reco::CaloCluster* clu, bool norm, bool bounded)
{

  double lmax = 30.; int nbinz = 20;
  TH1F *h_hgcal_sclusters_longitudinal_shower = new
   TH1F("h_hgcal_sclusters_longitudinal_shower","hgcal seed cluster shower longitudinal profile",nbinz,0.,lmax);
  h_hgcal_sclusters_longitudinal_shower->Sumw2();

  // calib as of SLHC21
  double clus_eta = clu->eta();
  const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
  const double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
  double scale = mip_*corr;
  double weight[30] =
   {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};

  GlobalPoint origPos = entryPosition(clu);

  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      const int layer = hgcid_.layer();
      GlobalVector radius; GlobalVector longitudinal; 
      radius = cellPos - origPos;
      longitudinal =  (radius.dot(showerDir_))*showerDir_.unit()/showerDir_.mag();
      // fill profile plots
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
        if (longitudinal.mag()<lmax) {
	  h_hgcal_sclusters_longitudinal_shower->Fill(longitudinal.mag(),theHit->energy()*weight[layer-1]/scale);
        }
      }
    }
  }

  TF1 *pgamma = new TF1("gammafunc",gamma,0.,lmax,3);
  // inject here parametrization results
  double lny = clu->energy()/hgcalEmParam_.getCriticalEnergy()>1. ? std::log(clu->energy()/hgcalEmParam_.getCriticalEnergy()) : 0.;
  double alpha = hgcalEmParam_.meanAlpha(lny);
  double tmax = hgcalEmParam_.meanT(lny);
  double sigmaalpha = alpha* hgcalEmParam_.sigmaLnAlpha(lny);
  if (sigmaalpha<0.) sigmaalpha = 1.;
  double sigmatmax = tmax*hgcalEmParam_.sigmaLnT(lny);
  if (sigmatmax<0.) sigmatmax = 1.;
  double beta = (alpha -1.) / tmax;  
  double invbeta = 1./beta; // beta of the gammaDist is 1/beta !
  double sigmainvbeta = pow(((1./(alpha-1.))*sigmaalpha),2) + pow(((tmax/pow((alpha-1.),2))*sigmaalpha),2);
  sigmainvbeta = sqrt(sigmainvbeta);

  // set the parameters
  double normalization = 1.;
  if (norm) normalization = clu->energy();
  pgamma->SetParameters(normalization,alpha,invbeta); 
  pgamma->SetParNames("Normalization","alpha","invbeta")	; 

  if (debug_) std::cout << " *** Fitting the longitudinal profile *** " << std::endl;
  if (debug_) std::cout << "alpha, invbeta " << alpha << " " << invbeta << std::endl;
  if (debug_) std::cout << "sigmaalpha, sigmainvbeta " << sigmaalpha << " " << sigmainvbeta << std::endl;

  // fix normalization if requested  
  if (norm) {
    pgamma->SetParLimits(0,normalization,normalization);
  }

  // set parameters boundaries if requested  
  if (bounded) {
    double alphamin = std::max(0.,alpha-3.*sigmaalpha), alphamax = std::max(0.,alpha+3*sigmaalpha);
    double invbetamin = std::max(0.,invbeta-3.*sigmainvbeta), invbetamax = std::max(0.,invbeta+3*sigmainvbeta);
    pgamma->SetParLimits(1,alphamin, alphamax);
    pgamma->SetParLimits(2,invbetamin,invbetamax); 
    if (debug_) std::cout << "alphamin, alphamax " << alphamin << " " << alphamax << std::endl;
    if (debug_) std::cout << "invbetamin, invbetamax " << invbetamin << " " << invbetamax << std::endl;
  }
  
  ///////////////////////////////////////
  // now fit the longitudinal distribution
  h_hgcal_sclusters_longitudinal_shower->Fit("gammafunc","","",0.,lmax);
  ///////////////////////////////////////

  // get fit results
  double fitparams[3] = {0.,0.,0.};
  pgamma->GetParameters(fitparams);
  double chi2 =pgamma->GetChisquare();
  int ndf = pgamma->GetNDF();
    
  delete h_hgcal_sclusters_longitudinal_shower; 
  delete pgamma;
  
  return HGCALFitResults(chi2,ndf,fitparams[0],fitparams[1],fitparams[2]);
  
}

double HGCALShowerBasedEmIdentification::longitudinalKolmogorov(const reco::CaloCluster* clu, bool dist)
{

  double lmax = 30.; int nbinz=20;
  TH1F *h_hgcal_sclusters_longitudinal_shower = new
   TH1F("h_hgcal_sclusters_longitudinal_shower","hgcal seed cluster shower longitudinal profile",nbinz,0.,lmax);
  h_hgcal_sclusters_longitudinal_shower->Sumw2();
  double sumhits = 0;

  double clus_eta = clu->eta();
  const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
  const double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
  double scale = mip_*corr;
  double weight[30] =
   {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};

  GlobalPoint origPos = entryPosition(clu);	 
  
  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      const int layer = hgcid_.layer();
      GlobalVector radius; GlobalVector longitudinal; 
      radius = cellPos - origPos;
      longitudinal =  (radius.dot(showerDir_))*showerDir_.unit()/showerDir_.mag();
      // fill profile plots
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
        if (longitudinal.mag()<lmax) {
          h_hgcal_sclusters_longitudinal_shower->Fill(longitudinal.mag(),theHit->energy()*weight[layer-1]/scale);
          sumhits += theHit->energy()*weight[layer-1]/scale;
	}  
      }
    }
  }

  TF1 *pgamma = new TF1("gammafunc",gamma,0.,lmax,3);
  // inject here parametrization results
  double lny = clu->energy()/hgcalEmParam_.getCriticalEnergy()>1. ? std::log(clu->energy()/hgcalEmParam_.getCriticalEnergy()) : 0.;
  double alpha = hgcalEmParam_.meanAlpha(lny);
  double tmax = hgcalEmParam_.meanT(lny);
  double beta = (alpha -1.) / tmax;  
  double invbeta = 1./beta; // beta of the gammaDist is 1/beta !

  // use the measured length to constrain the parameters
  double alphascale = 1.03; // adhoc factor ..
  pgamma->SetParameters(1.,(length(clu)/invbeta)*alphascale,invbeta); 
  TH1F *h_hgcal_sclusters_expected_longitudinal_shower = new
   TH1F("h_hgcal_sclusters_expected_longitudinal_shower","hgcal seed cluster shower longitudinal profile",nbinz,0.,lmax);
  h_hgcal_sclusters_expected_longitudinal_shower->Sumw2();
  if (debug_) std::cout << "Filling an histogram with the gamma function " << std::endl;
  if (pgamma->Integral(0.,lmax)!=0.) h_hgcal_sclusters_expected_longitudinal_shower->FillRandom("gammafunc");
  // normalize prediction
  if (h_hgcal_sclusters_expected_longitudinal_shower->GetEntries()!=0.) h_hgcal_sclusters_expected_longitudinal_shower->Scale(sumhits/h_hgcal_sclusters_expected_longitudinal_shower->GetEntries());

  // now get the kolmogorov result
  double kolmogorov;
  std::string option = "M";
  if (!dist) {
    kolmogorov = 0.;
    // if dist is false, get the Kolmogorov probability
    if (h_hgcal_sclusters_longitudinal_shower->Integral()!=0. && h_hgcal_sclusters_expected_longitudinal_shower->Integral()!=0.) 
     kolmogorov = h_hgcal_sclusters_longitudinal_shower->KolmogorovTest(h_hgcal_sclusters_expected_longitudinal_shower);
  } else {
    // if dist is true, get the Kolmogorov distance  
    kolmogorov = 10000.;
    if (h_hgcal_sclusters_longitudinal_shower->Integral()!=0. && h_hgcal_sclusters_expected_longitudinal_shower->Integral()!=0.) 
     kolmogorov = h_hgcal_sclusters_longitudinal_shower->KolmogorovTest(h_hgcal_sclusters_expected_longitudinal_shower,"M");
  }
  
  delete h_hgcal_sclusters_longitudinal_shower; 
  delete h_hgcal_sclusters_expected_longitudinal_shower; 
  delete pgamma;
  
  // result
  if (debug_) std::cout << "Kolmogorov test longitudinal profile " << kolmogorov << std::endl;

  return kolmogorov;

}


double HGCALShowerBasedEmIdentification::transverseKolmogorov(const reco::CaloCluster* clu, bool dist)
{

  double rmax = 10.; int nbinr=10;
  TH1F *h_hgcal_sclusters_transversal_shower = new 
   TH1F("h_hgcal_sclusters_transversal_shower","hgcal seed cluster shower transversal profile",nbinr,0.,rmax);
  h_hgcal_sclusters_transversal_shower->Sumw2();
  double sumhits = 0;

  double clus_eta = clu->eta();
  const double _coef_a = 80.0837, _coef_c = 0.0472817, _coef_d = -0.266294, _coef_e = 0.34684;
  const double corr = _coef_a*fabs(std::tanh(clus_eta))/(1.-(_coef_c*pow(clus_eta,2)+_coef_d*fabs(clus_eta)+_coef_e));
  double scale = mip_*corr;
  double weight[30] =
   {0.080,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.62,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,0.809,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239,1.239};

  for (unsigned int ih=0;ih<clu->hitsAndFractions().size();++ih) {
    const DetId & id_ = (clu->hitsAndFractions())[ih].first ;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
    if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
      const HGCEEDetId & hgcid_ = HGCEEDetId((clu->hitsAndFractions())[ih].first) ;
      GlobalPoint cellPos = geometry_->getPosition(hgcid_);	     
      const int layer = hgcid_.layer();
      double lambdatoaverage = (cellPos.z()-showerPos_.z())/showerDir_.z();
      GlobalPoint axisPos = showerPos_ + lambdatoaverage*showerDir_;	 
      double rtaxis = (cellPos.perp()-axisPos.perp()) ;	     				
      // fill profile plots
      if (!withPileup_ || theHit->energy()>minenergy_*mip_) {
        if (rtaxis<rmax) {	
          h_hgcal_sclusters_transversal_shower->Fill(rtaxis,theHit->energy()*weight[layer-1]/scale);
	  sumhits += theHit->energy()*weight[layer-1]/scale;
        }
      }
    }
  }

  TF1 *ptrans = new TF1("transfunc",doublexp,0.,rmax,5);
  double corefraction = 0.8808;
  double coreradius = 2.589*1.94;
  double tailradius = 1.18*1.94;

  ptrans->SetParameters(1.,corefraction,coreradius,tailradius,std::abs(std::cos(showerDir_.theta()))); 
  TH1F *h_hgcal_sclusters_expected_transversal_shower = new
   TH1F("h_hgcal_sclusters_expected_transversal_shower","hgcal seed cluster shower transversal profile",nbinr,0.,rmax);
  h_hgcal_sclusters_expected_transversal_shower->Sumw2();
  if (debug_) std::cout << "Filling an histogram with the transversal profile function " << std::endl;
  h_hgcal_sclusters_expected_transversal_shower->FillRandom("transfunc");
  // normalize prediction
  if (h_hgcal_sclusters_expected_transversal_shower->GetEntries()!=0.) h_hgcal_sclusters_expected_transversal_shower->Scale(sumhits/h_hgcal_sclusters_expected_transversal_shower->GetEntries());

  // now get the kolmogorov result
  double kolmogorov;
  std::string option = "M";
  if (!dist) {
    kolmogorov = 0.;
    // if dist is false, get the Kolmogorov probability
    if (h_hgcal_sclusters_transversal_shower->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) 
     kolmogorov = h_hgcal_sclusters_transversal_shower->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower);
  } else {
    // if dist is true, get the Kolmogorov distance  
    kolmogorov = 10000.;
    if (h_hgcal_sclusters_transversal_shower->Integral()!=0. && h_hgcal_sclusters_expected_transversal_shower->Integral()!=0.) 
     kolmogorov = h_hgcal_sclusters_transversal_shower->KolmogorovTest(h_hgcal_sclusters_expected_transversal_shower,"M");
  }
  
  delete h_hgcal_sclusters_transversal_shower; 
  delete h_hgcal_sclusters_expected_transversal_shower; 
  delete ptrans;

  // result
  if (debug_) std::cout << "Kolmogorov test trasnversal profile " << kolmogorov << std::endl;

  return kolmogorov;  

}

bool HGCALShowerBasedEmIdentification::cutSigmaetaeta(const reco::CaloCluster* clu)
{
  return (sigmaetaeta(clu)<cutSigmaetaeta_);
}

bool HGCALShowerBasedEmIdentification::cutStartPosition(const reco::CaloCluster* clu)
{
  return (std::abs(startPosition(clu).z())<cutStartPosition_);
}

bool HGCALShowerBasedEmIdentification::cutHadOverEm(const reco::CaloCluster* clu)
{
  return (hadOverEm(clu)<cutHoverem_); 
}

bool HGCALShowerBasedEmIdentification::cutLengthCompatibility(const reco::CaloCluster* clu)
{
  return (std::abs(lengthCompatibility(clu))<cutLengthCompatibility_); 
}

