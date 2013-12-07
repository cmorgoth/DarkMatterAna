#include "DM_ZJetsNuNuV3.hh"
#include <iostream>
#include <iomanip>

ZJetsNuNu::ZJetsNuNu(){
  
};

ZJetsNuNu::ZJetsNuNu(int MetIndex ): BaseDM( MetIndex, "ZJetsNuNu" ){
  
  T = new TChain("outTree");
  T->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_50_HT_100_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_100_HT_200_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_200_HT_400_pu_mu_LooseBtag_sf.root");
  T->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_400_HT_inf_pu_mu_LooseBtag_sf.root");
  
  effT = new TChain("effTree");
  effT->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_50_HT_100_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_100_HT_200_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_200_HT_400_pu_mu_LooseBtag_sf.root");
  effT->Add("/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_400_HT_inf_pu_mu_LooseBtag_sf.root");
};

ZJetsNuNu::~ZJetsNuNu(){
  delete T;
};

