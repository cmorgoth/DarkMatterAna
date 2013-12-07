#include "DM_WJetsHTBinsV3.hh"
#include <iostream>
#include <iomanip>

WJetsHTBins::WJetsHTBins(){
  
};

WJetsHTBins::WJetsHTBins(int MetIndex ): BaseDM( MetIndex, "WJetsHTBins" ){
  
  T = new TChain("outTree");
   T->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_150_HT_200_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_200_HT_250_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_250_HT_300_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_300_HT_400_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_400_HT_Inf_pu_mu_LooseBtag_sf.root");
  effT = new TChain("effTree");
   effT->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_150_HT_200_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_200_HT_250_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_250_HT_300_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_300_HT_400_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/WJets/BtagCorrMC/WJetsToLNu_400_HT_Inf_pu_mu_LooseBtag_sf.root");
  
};

WJetsHTBins::~WJetsHTBins(){
  delete T;
};

