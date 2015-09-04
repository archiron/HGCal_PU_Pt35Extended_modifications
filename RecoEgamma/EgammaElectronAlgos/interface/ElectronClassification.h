#ifndef ElectronClassification_H
#define ElectronClassification_H

//===================================================================
// Purpose: electron classification for CMS phase II upgrade
// Author: Claude Charlot - LLR- Ecole Polytechnique
// 02/2015
//===================================================================

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


class ElectronClassification
{
 public:

  ElectronClassification(){}

  void classify( reco::GsfElectron & ) ;
  void refineWithPflow( reco::GsfElectron & ) ;

 private:

//  void classify(const reco::GsfElectron &);

//  bool isInCrack(float eta) const;
//  bool isInEtaGaps(float eta) const;
//  bool isInPhiGaps(float phi) const;

//  reco::GsfElectron::Classification electronClass_;

};

#endif




