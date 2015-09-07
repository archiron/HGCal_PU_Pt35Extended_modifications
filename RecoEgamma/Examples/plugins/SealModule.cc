//#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Utilities/interface/typelookup.h"

#include "HGCALGsfElectronAnalyzer.h"
#include "HGCALElectronClusterAnalyzer.h"


#include "CommonTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

typedef Merger<reco::SuperClusterCollection> EgammaSuperClusterMerger;
DEFINE_FWK_MODULE(EgammaSuperClusterMerger);
DEFINE_FWK_MODULE(HGCALElectronClusterAnalyzer);
DEFINE_FWK_MODULE(HGCALGsfElectronAnalyzer);
