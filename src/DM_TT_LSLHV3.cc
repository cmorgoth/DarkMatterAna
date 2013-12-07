#include "DM_TT_LSLHV3.hh"
#include <iostream>
#include <iomanip>

TTJets::TTJets(){
  
};

TTJets::TTJets(int MetIndex ): BaseDM( MetIndex, "TTJets" ){
  
  T = new TChain("outTree");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag_sf.root");
  
  effT = new TChain("effTree");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag_sf.root");
  
};

TTJets::~TTJets(){
  delete T;
};

