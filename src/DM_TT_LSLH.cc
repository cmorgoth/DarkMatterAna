#include "DM_TT_LSLH.hh"
#include <iostream>

TTJets::TTJets(){ };

TTJets::TTJets(int metIndex){
  
  MRMin = 200.0;//nominal 150.0
  RSQMin = 0.5;//nominal 0.5
  this->metIndex = metIndex;

  double N_In = 0 , Ntot = 0;
  double Nexp = 0;
  
  //////////////////////////////////////////
  ////////////Trigger Efficiency///////////
  /////////// files and TEfficiency///////
  ///////////////////////////////////////
  TFile* file = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_DoubleMuonPD_Final.root");
  eff = (TEfficiency*)file->Get("Eff2d");

  TFile* file1 = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_SignleElePD_Final.root");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
  
  ///////////////////////////////////
  ////////Opening Ntuples///////////
  //////////////////////////////////
  
  //Leptonic
  //F = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJets_ILV_FullyLeptMGDecaysTauola_BtagCorr_TightMedLoose.root");
  //F = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJetsFullyLept_pu_mu_tot.root");
  F = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsFullyLept_pu_mu_LooseAndTightBtag.root");
  
  TTree* effT = (TTree*)F->Get("effTree");
  T = (TTree*)F->Get("outTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  weight0 = Lumi*sigma0/Ntot;
  Nexp += weight0*T->GetEntries(); 
  
  //SemiLeptonic
  //F1 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJetsSemiLept_BtagCorr_Loose.root");
  //F1 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJetsSemiLept_pu_mu_tot.root");
  F1 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsSemiLept_pu_mu_LooseAndTightBtag.root");

  effT = (TTree*)F1->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  T1 = (TTree*)F1->Get("outTree");
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  
  weight1 = Lumi*sigma1/Ntot;
  Nexp += weight1*T1->GetEntries(); 
  
  //Hadronic
  //F2 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJetsHad_BtagCorr_Loose.root");
  //F2 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/TTJetsHad_pu_mu_tot.root");
  F2 = TFile::Open("/media/data/cmorgoth/Data/DMData/TTJets/BtagCorrMC/TTJetsHad_pu_mu_LooseAndTightBtag.root");
  
  effT = (TTree*)F2->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  T2 = (TTree*)F2->Get("outTree");
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  weight2 = Lumi*sigma2/Ntot;
  Nexp += weight2*T2->GetEntries(); 
  
  
  std::cout << "Weight 0: " << weight0 << std::endl;
  std::cout << "Weight 1: " << weight1 << std::endl;
  std::cout << "Weight 2: " << weight2 << std::endl;
  
  
};

TTJets::TTJets(const char* FileName ){
  
  //F = new TFile(FileName);
  //T0 = (TTree*)F->Get("outTree");
  
};

TTJets::~TTJets(){
  delete T;
  delete T1;
  delete T2;
  
  delete F;
  delete F1;
  delete F2;
};


bool TTJets::PrintEvents(){
  
  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;

  //Leptonic
  TTree* effT = (TTree*)F->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In*weight0;
    Nt_PV += N_PV*weight0;
    Nt_2J += N_2J*weight0;
    Nt_0b += N_0b*weight0;
    Nt_LepVeto += N_LepVeto*weight0;
  }
  
  //SemiLeptonic
  effT = (TTree*)F1->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In*weight1;
    Nt_PV += N_PV*weight1;
    Nt_2J += N_2J*weight1;
    Nt_0b += N_0b*weight1;
    Nt_LepVeto += N_LepVeto*weight1;
  }
  
  //Hadronic
  effT = (TTree*)F2->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In*weight2;
    Nt_PV += N_PV*weight2;
    Nt_2J += N_2J*weight2;
    Nt_0b += N_0b*weight2;
    Nt_LepVeto += N_LepVeto*weight2;
  }

  std::cout << "========================================" << std::endl;
  std::cout << "================TTJets===================" << std::endl;
  std::cout << "========================================" << std::endl;

  std::cout << "TTJets weighted Nt_In: " << NtotGen << std::endl;
  std::cout << "TTJets weighted Nt_PV: " << Nt_PV << std::endl;
  std::cout << "TTJets weighted Nt_2J: " << Nt_2J << std::endl;
  std::cout << "TTJets weighted Nt_0b: " << Nt_0b << std::endl;
  std::cout << "TTJets weighted Nt_LepVeto: " << Nt_LepVeto << std::endl;

  std::cout << "TTJets weighted Nt_tree: " << T->GetEntries()*weight0 +T1->GetEntries()*weight1 +  T2->GetEntries()*weight2 << "\n\n" << std::endl;

  
  //After Selection cuts
  double RSQ[4], MR[4];
  int BOX, nBtag;
  
  double Nt_MR_RSQ_cut0BTag = 0.0, Nt_2muBox0BTag = 0.0,  Nt_1muBox0BTag = 0.0, Nt_0muBox0BTag = 0.0, \
    Nt_MR_RSQ_cut = 0.0, Nt_2muBox = 0.0,  Nt_1muBox = 0.0, Nt_0muBox = 0.0, Nt_Btags = 0.0;
  
  
  SetStatus();
  //Leptonic
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += weight0;
      
      if( BOX == 0 )Nt_0muBox += weight0;
      if( BOX == 1 )Nt_1muBox += weight0;
      if( BOX == 2 )Nt_2muBox += weight0;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag += weight0;
	
	if( BOX == 0 )Nt_0muBox0BTag += weight0;
	if( BOX == 1 )Nt_1muBox0BTag += weight0;
	if( BOX == 2 )Nt_2muBox0BTag += weight0;
	
      }
      
    }
    
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  //SemiLeptonic
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += weight1;
      
      if( BOX == 0 )Nt_0muBox += weight1;
      if( BOX == 1 )Nt_1muBox += weight1;
      if( BOX == 2 )Nt_2muBox += weight1;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag += weight1;
	
	if( BOX == 0 )Nt_0muBox0BTag += weight1;
	if( BOX == 1 )Nt_1muBox0BTag += weight1;
	if( BOX == 2 )Nt_2muBox0BTag += weight1;
	
      }
      
    }

  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  //Hadronic
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += weight2;
      
      if( BOX == 0 )Nt_0muBox += weight2;
      if( BOX == 1 )Nt_1muBox += weight2;
      if( BOX == 2 )Nt_2muBox += weight2;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag += weight2;
	
	if( BOX == 0 )Nt_0muBox0BTag += weight2;
	if( BOX == 1 )Nt_1muBox0BTag += weight2;
	if( BOX == 2 )Nt_2muBox0BTag += weight2;
	
      }
      
    }
    
  }
  T2->SetBranchStatus("*", 0);

  std::cout << "TTJets weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut << std::endl;
  std::cout << "TTJets weighted Nt_0muBox: " << Nt_0muBox << std::endl;
  std::cout << "TTJets weighted Nt_1muBox: " << Nt_1muBox << std::endl;
  std::cout << "TTJets weighted Nt_2muBox: " << Nt_2muBox << std::endl;

  std::cout << "TTJets weighted Nt_MR_RSQ_cut0Btag: " << Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "TTJets weighted Nt_0muBox0Btag: " << Nt_0muBox0BTag << std::endl;
  std::cout << "TTJets weighted Nt_1muBox0Btag: " << Nt_1muBox0BTag << std::endl;
  std::cout << "TTJets weighted Nt_2muBox0Btag: " << Nt_2muBox0BTag << std::endl;

  std::cout << "TTJets Btag EVENTS: " << Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
  
};


TH1F TTJets::PlotMR_2Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* MR2 = new TH1F("MR2", "MR2_BOX", MR_Bins, MR_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( /*N_Jets == 2 &&*/ BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR2->Fill(MR[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( /*N_Jets == 2 && */BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR2->Fill(MR[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( /*N_Jets == 2 && */BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR2->Fill(MR[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *MR2;

};


TH1F TTJets::PlotMR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* MR1 = new TH1F("MR1", "MR1BOX", MR_Bins, MR_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *MR1;
  
};

TH1F TTJets::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0; 
  int BOX, N_Jets, nBtag[2];
  TH1F* MR0 = new TH1F("MR0", "MR0BOX", MR_Bins, MR_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *MR0;
};


TH1F  TTJets::PlotRSQ_2Box(){

  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ2 = new TH1F("RSQ2", "RSQ2_BOX", RSQ_Bins, RSQ_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ2->Fill(RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ2->Fill(RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ2->Fill(RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ2;

};

TH1F  TTJets::PlotRSQ_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ1 = new TH1F("RSQ1", "RSQ1BOX", RSQ_Bins, RSQ_BinArr);

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ1;
  
};


TH1F  TTJets::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ0 = new TH1F("RSQ0", "RSQ0BOX", RSQ_Bins, RSQ_BinArr);
  
  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ0;
  
};


TH2F TTJets::PlotRSQ_vs_MR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH2F* RSQ_MR_0BOX = new TH2F("RSQ_MR_0BOX", "RSQ_VS_MR_0BOX", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  RSQ_MR_0BOX->Sumw2();

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F TTJets::PlotRSQ_vs_MR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH2F* RSQ_MR_1BOX = new TH2F("RSQ_MR_1BOX", "RSQ_VS_MR_1BOX", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  RSQ_MR_1BOX->Sumw2();

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ_MR_1BOX;
  
};

TH2F TTJets::PlotRSQ_vs_MR_2Box(){

  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH2F* RSQ_MR_2BOX = new TH2F("RSQ_MR_2BOX", "RSQ_VS_MR_2BOX", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  RSQ_MR_2BOX->Sumw2();

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_2BOX->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_2BOX->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_2BOX->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
    }
  }
  T2->SetBranchStatus("*", 0);

  return *RSQ_MR_2BOX;
  
};


std::vector<TH1F*> TTJets::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30], pu_w, mu_w;
  int BOX, nBtag[2], N_Jets;
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
   
  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, 0, 1000 );
    }
  }
  
  MET[9] = new TH1F( "NJETS0_Z", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( "NJETS1_Z", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( "NJETS2_Z", "NJETS 2 BOX", 9, 1, 10);
  
  SetMetStatus();
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  T->SetBranchAddress("N_Jets", &N_Jets);

  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  
  float metmag = 0.;
  float metmagcorr = 0.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
	MET[0]->Fill(metmag, weight0*hltWeight*pu_w*mu_w);
        MET[1]->Fill(metmagcorr, weight0*hltWeight*pu_w*mu_w);
        MET[2]->Fill(metmagcorr-metmag, weight0*hltWeight*pu_w*mu_w);
	MET[9]->Fill(N_Jets,weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight0*hltWeight*pu_w*mu_w);
        MET[4]->Fill(metmagcorr, weight0*hltWeight*pu_w*mu_w);
        MET[5]->Fill(metmagcorr-metmag, weight0*hltWeight*pu_w*mu_w);
	MET[10]->Fill(N_Jets,weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight0*hltWeight*pu_w*mu_w);
        MET[7]->Fill(metmagcorr, weight0*hltWeight*pu_w*mu_w);
        MET[8]->Fill(metmagcorr-metmag, weight0*hltWeight*pu_w*mu_w);
	MET[11]->Fill(N_Jets,weight0*hltWeight*pu_w*mu_w);
      }
    }
    
  }
  T->SetBranchStatus("*", 0);

  SetMetStatus1();
  T1->SetBranchAddress("pu_w", &pu_w);
  T1->SetBranchAddress("mu_w", &mu_w);
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  
  //T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  T1->SetBranchAddress("ht", &ht);
  T1->SetBranchAddress("mht", &mht[0]);
  T1->SetBranchAddress("metX", metX);
  T1->SetBranchAddress("metCorrX", metcorrX);
  T1->SetBranchAddress("metY", metY);
  T1->SetBranchAddress("metCorrY", metcorrY);
  T1->SetBranchAddress("N_Jets", &N_Jets);

  T1->SetBranchAddress("pTHem1", &pTHem1);
  T1->SetBranchAddress("pTHem2", &pTHem2);
  T1->SetBranchAddress("etaHem1", &etaHem1);
  T1->SetBranchAddress("etaHem2", &etaHem2);
  T1->SetBranchAddress("phiHem1", &phiHem1);
  T1->SetBranchAddress("phiHem2", &phiHem2);
  
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[0]->Fill(metmag, weight1*hltWeight*pu_w*mu_w);
        MET[1]->Fill(metmagcorr, weight1*hltWeight*pu_w*mu_w);
        MET[2]->Fill(metmagcorr-metmag, weight1*hltWeight*pu_w*mu_w);
	MET[9]->Fill(N_Jets,weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight1*hltWeight*pu_w*mu_w);
        MET[4]->Fill(metmagcorr, weight1*hltWeight*pu_w*mu_w);
        MET[5]->Fill(metmagcorr-metmag, weight1*hltWeight*pu_w*mu_w);
	MET[10]->Fill(N_Jets,weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight1*hltWeight*pu_w*mu_w);
        MET[7]->Fill(metmagcorr, weight1*hltWeight*pu_w*mu_w);
        MET[8]->Fill(metmagcorr-metmag, weight1*hltWeight*pu_w*mu_w);
	MET[11]->Fill(N_Jets,weight1*hltWeight*pu_w*mu_w);
      }
    }
    
  }
  T1->SetBranchStatus("*", 0);

  SetMetStatus2();
  T2->SetBranchAddress("pu_w", &pu_w);
  T2->SetBranchAddress("mu_w", &mu_w);
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  
  //T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  T2->SetBranchAddress("ht", &ht);
  T2->SetBranchAddress("mht", &mht[0]);
  T2->SetBranchAddress("metX", metX);
  T2->SetBranchAddress("metCorrX", metcorrX);
  T2->SetBranchAddress("metY", metY);
  T2->SetBranchAddress("metCorrY", metcorrY);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  
  T2->SetBranchAddress("pTHem1", &pTHem1);
  T2->SetBranchAddress("pTHem2", &pTHem2);
  T2->SetBranchAddress("etaHem1", &etaHem1);
  T2->SetBranchAddress("etaHem2", &etaHem2);
  T2->SetBranchAddress("phiHem1", &phiHem1);
  T2->SetBranchAddress("phiHem2", &phiHem2);

  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
	MET[0]->Fill(metmag, weight2*hltWeight*pu_w*mu_w);
        MET[1]->Fill(metmagcorr, weight2*hltWeight*pu_w*mu_w);
        MET[2]->Fill(metmagcorr-metmag, weight2*hltWeight*pu_w*mu_w);
	MET[9]->Fill(N_Jets,weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight2*hltWeight*pu_w*mu_w);
        MET[4]->Fill(metmagcorr, weight2*hltWeight*pu_w*mu_w);
        MET[5]->Fill(metmagcorr-metmag, weight2*hltWeight*pu_w*mu_w);
	MET[10]->Fill(N_Jets,weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight2*hltWeight*pu_w*mu_w);
        MET[7]->Fill(metmagcorr, weight2*hltWeight*pu_w*mu_w);
	MET[8]->Fill(metmagcorr-metmag, weight2*hltWeight*pu_w*mu_w);
	MET[11]->Fill(N_Jets,weight2*hltWeight*pu_w*mu_w);
      }
    }
  }
  T2->SetBranchStatus("*", 0);

  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec;
  
};

bool TTJets::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int TTJets::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH1F*> TTJets::DoubleMuBoxPlots(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30], Mu_E[2], Mu_Px[2], Mu_Py[2], Mu_Pz[2], pu_w, mu_w;
  int BOX, N_Jets, nBtag[2];
  double hltWeight;
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;

  std::vector<TH1F*> vec_plot;
  TH1F* plot_2mu[6];

  plot_2mu[0] = new TH1F( "Mass", "Mass", 15, .0, 500.);
  plot_2mu[1] = new TH1F( "Angle", "Angle", 15, .0, 2*3.1416);
  plot_2mu[2] = new TH1F( "Pt1", "Pt1", 15, .0, 500);
  plot_2mu[3] = new TH1F( "Pt2", "Pt2", 15, .0, 500);
  plot_2mu[4] = new TH1F( "Eta1", "Eta1", 15, -3.0, 3.0);
  plot_2mu[5] = new TH1F( "Eta2", "Eta2", 15, -3.0, 3.0);

  SetMetStatus();
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);

  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);

    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      TLorentzVector mu1(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
      TLorentzVector mu2(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
      TLorentzVector sum_mu;
      sum_mu = mu1 + mu2;
      double Mass = sum_mu.M();
      double angle = mu1.Angle(mu2.Vect());
      double pt1 = mu1.Pt();
      double pt2 = mu2.Pt();
      double eta1 = asinh(Mu_Pz[0]/pt1);
      double eta2 = asinh(Mu_Pz[1]/pt2);
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      if( hltWeight == 0.0 )hltWeight = 1.0;
      if(Mass > .0 && Mass < 1800.){
        plot_2mu[0]->Fill(Mass,weight0*hltWeight*pu_w*mu_w);
        plot_2mu[1]->Fill(angle,weight0*hltWeight*pu_w*mu_w);
        plot_2mu[2]->Fill(pt1,weight0*hltWeight*pu_w*mu_w);
        plot_2mu[3]->Fill(pt2,weight0*hltWeight*pu_w*mu_w);
        plot_2mu[4]->Fill(eta1,weight0*hltWeight*pu_w*mu_w);
        plot_2mu[5]->Fill(eta2,weight0*hltWeight*pu_w*mu_w);
      }
    }
  }

  SetMetStatus1();
  T1->SetBranchAddress("pu_w", &pu_w);
  T1->SetBranchAddress("mu_w", &mu_w);
  T1->SetBranchStatus("Mu_E",1);
  T1->SetBranchStatus("Mu_Px",1);
  T1->SetBranchStatus("Mu_Py",1);
  T1->SetBranchStatus("Mu_Pz",1);
  T1->SetBranchAddress("Mu_E", Mu_E);
  T1->SetBranchAddress("Mu_Px", Mu_Px);
  T1->SetBranchAddress("Mu_Py", Mu_Py);
  T1->SetBranchAddress("Mu_Pz", Mu_Pz);
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  
  //T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  T1->SetBranchAddress("ht", &ht);
  T1->SetBranchAddress("mht", &mht[0]);
  T1->SetBranchAddress("metX", metX);
  T1->SetBranchAddress("metCorrX", metcorrX);
  T1->SetBranchAddress("metY", metY);

  T1->SetBranchAddress("pTHem1", &pTHem1);
  T1->SetBranchAddress("pTHem2", &pTHem2);
  T1->SetBranchAddress("etaHem1", &etaHem1);
  T1->SetBranchAddress("etaHem2", &etaHem2);
  T1->SetBranchAddress("phiHem1", &phiHem1);
  T1->SetBranchAddress("phiHem2", &phiHem2);
  
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);

    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      TLorentzVector mu1(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
      TLorentzVector mu2(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
      TLorentzVector sum_mu;
      sum_mu = mu1 + mu2;
      double Mass = sum_mu.M();
      double angle = mu1.Angle(mu2.Vect());
      double pt1 = mu1.Pt();
      double pt2 = mu2.Pt();
      double eta1 = asinh(Mu_Pz[0]/pt1);
      double eta2 = asinh(Mu_Pz[1]/pt2);
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      if(Mass > .0 && Mass < 1800.){
        plot_2mu[0]->Fill(Mass,weight1*hltWeight*pu_w*mu_w);
        plot_2mu[1]->Fill(angle,weight1*hltWeight*pu_w*mu_w);
        plot_2mu[2]->Fill(pt1,weight1*hltWeight*pu_w*mu_w);
        plot_2mu[3]->Fill(pt2,weight1*hltWeight*pu_w*mu_w);
        plot_2mu[4]->Fill(eta1,weight1*hltWeight*pu_w*mu_w);
        plot_2mu[5]->Fill(eta2,weight1*hltWeight*pu_w*mu_w);
      }
    }
  }
  
  SetMetStatus2();
  T2->SetBranchAddress("pu_w", &pu_w);
  T2->SetBranchAddress("mu_w", &mu_w);
  T2->SetBranchStatus("Mu_E",1);
  T2->SetBranchStatus("Mu_Px",1);
  T2->SetBranchStatus("Mu_Py",1);
  T2->SetBranchStatus("Mu_Pz",1);
  T2->SetBranchAddress("Mu_E", Mu_E);
  T2->SetBranchAddress("Mu_Px", Mu_Px);
  T2->SetBranchAddress("Mu_Py", Mu_Py);
  T2->SetBranchAddress("Mu_Pz", Mu_Pz);
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  
  //T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  T2->SetBranchAddress("ht", &ht);
  T2->SetBranchAddress("mht", &mht[0]);
  T2->SetBranchAddress("metX", metX);
  T2->SetBranchAddress("metCorrX", metcorrX);
  T2->SetBranchAddress("metY", metY);

  T2->SetBranchAddress("pTHem1", &pTHem1);
  T2->SetBranchAddress("pTHem2", &pTHem2);
  T2->SetBranchAddress("etaHem1", &etaHem1);
  T2->SetBranchAddress("etaHem2", &etaHem2);
  T2->SetBranchAddress("phiHem1", &phiHem1);
  T2->SetBranchAddress("phiHem2", &phiHem2);

  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      TLorentzVector mu1(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
      TLorentzVector mu2(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
      TLorentzVector sum_mu;
      sum_mu = mu1 + mu2;
      double Mass = sum_mu.M();
      double angle = mu1.Angle(mu2.Vect());
      double pt1 = mu1.Pt();
      double pt2 = mu2.Pt();
      double eta1 = asinh(Mu_Pz[0]/pt1);
      double eta2 = asinh(Mu_Pz[1]/pt2);
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      if(Mass > .0 && Mass < 1800.){
        plot_2mu[0]->Fill(Mass,weight2*hltWeight*pu_w*mu_w);
        plot_2mu[1]->Fill(angle,weight2*hltWeight*pu_w*mu_w);
        plot_2mu[2]->Fill(pt1,weight2*hltWeight*pu_w*mu_w);
        plot_2mu[3]->Fill(pt2,weight2*hltWeight*pu_w*mu_w);
        plot_2mu[4]->Fill(eta1,weight2*hltWeight*pu_w*mu_w);
        plot_2mu[5]->Fill(eta2,weight2*hltWeight*pu_w*mu_w);
      }
    }
  }
    
  for(int j = 0; j < 6; j++){
    vec_plot.push_back(plot_2mu[j]);
  }
  
  return vec_plot;

};


std::vector<TH1F*> TTJets::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w;
  int BOX, N_Jets, nBtag[2], nBtagMed;

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 6; l++ ){
    name = TString(Form("MR_1D_TT_%dmu_Box",l));
    name1 = TString(Form("R2_1D_TT_%dmu_Box",l));
    Razor1D[2*l] = new TH1F( name, name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F( name1, name1, RSQ_Bins, RSQ_BinArr);
    if(l < 3){
      Razor1D[2*l]->Sumw2();
      Razor1D[2*l+1]->Sumw2();
    }
  }

  SetStatus();
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);    
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight0*hltWeight*pu_w*mu_w);
        Razor1D[1]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
	Razor1D[6]->Fill(MR[metIndex], weight0*hltWeight*pu_w*mu_w);
	Razor1D[7]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex],  weight0*hltWeight*pu_w*mu_w);
        Razor1D[3]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
	Razor1D[8]->Fill(MR[metIndex],  weight0*hltWeight*pu_w*mu_w);
        Razor1D[9]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight0*hltWeight*pu_w*mu_w);
        Razor1D[5]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
	Razor1D[10]->Fill(MR[metIndex], weight0*hltWeight*pu_w*mu_w);
	Razor1D[11]->Fill(RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("pu_w", &pu_w);
  T1->SetBranchAddress("mu_w", &mu_w);
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  //T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  T1->SetBranchAddress("pTHem1", &pTHem1);
  T1->SetBranchAddress("pTHem2", &pTHem2);
  T1->SetBranchAddress("etaHem1", &etaHem1);
  T1->SetBranchAddress("etaHem2", &etaHem2);
  T1->SetBranchAddress("phiHem1", &phiHem1);
  T1->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    
    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
        Razor1D[1]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
	Razor1D[6]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
        Razor1D[7]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
        Razor1D[3]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
	Razor1D[8]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
	Razor1D[9]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
        Razor1D[5]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
	Razor1D[10]->Fill(MR[metIndex], weight1*hltWeight*pu_w*mu_w);
        Razor1D[11]->Fill(RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("pu_w", &pu_w);
  T2->SetBranchAddress("mu_w", &mu_w);
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  //T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("nBtagTCorr", &nBtag[1]);

  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  T2->SetBranchAddress("pTHem1", &pTHem1);
  T2->SetBranchAddress("pTHem2", &pTHem2);
  T2->SetBranchAddress("etaHem1", &etaHem1);
  T2->SetBranchAddress("etaHem2", &etaHem2);
  T2->SetBranchAddress("phiHem1", &phiHem1);
  T2->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);

    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[1]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
	Razor1D[6]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[7]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[3]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
	Razor1D[8]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[9]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[5]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
	Razor1D[10]->Fill(MR[metIndex], weight2*hltWeight*pu_w*mu_w);
        Razor1D[11]->Fill(RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }
    }
  }
  T2->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }
  
  return Razor1DVec;
  
};


std::vector<TH2F*> TTJets::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w;
  int BOX, N_Jets, nBtag[2], nBtagMed;

  std::vector< TH2F* > Razor2DVec;
  TH2F* Razor2D[3];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("MR_2D_TT_%dmu_Box",l));
    name1 = TString(Form("R2_2D_TT_%dmu_Box",l));
    Razor2D[l] = new TH2F( name, name, MR_Bins, MR_BinArr,  RSQ_Bins, RSQ_BinArr);
    Razor2D[l]->Sumw2();
  }

  SetStatus();
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    TLorentzVector j1;
    TLorentzVector j2;  
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight*pu_w*mu_w);
      }
    }
  }
  T->SetBranchStatus("*", 0);

  SetStatus1();
  T1->SetBranchAddress("pu_w", &pu_w);
  T1->SetBranchAddress("mu_w", &mu_w);
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  //T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  T1->SetBranchAddress("pTHem1", &pTHem1);
  T1->SetBranchAddress("pTHem2", &pTHem2);
  T1->SetBranchAddress("etaHem1", &etaHem1);
  T1->SetBranchAddress("etaHem2", &etaHem2);
  T1->SetBranchAddress("phiHem1", &phiHem1);
  T1->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    
    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight*pu_w*mu_w);
      }
    }
  }
  T1->SetBranchStatus("*", 0);

  SetStatus2();
  T2->SetBranchAddress("pu_w", &pu_w);
  T2->SetBranchAddress("mu_w", &mu_w);
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  //T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  //T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("nBtagTCorr", &nBtag[1]);
  
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  T2->SetBranchAddress("pTHem1", &pTHem1);
  T2->SetBranchAddress("pTHem2", &pTHem2);
  T2->SetBranchAddress("etaHem1", &etaHem1);
  T2->SetBranchAddress("etaHem2", &etaHem2);
  T2->SetBranchAddress("phiHem1", &phiHem1);
  T2->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);

    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    
    nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight*pu_w*mu_w);
      }
    }
  }
  T2->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;
  
};

bool TTJets::SetStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets", 1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("CSV", 1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  //T->SetBranchStatus("nBtagMCorr", 1);
  //T->SetBranchStatus("nBtagMCorrUp", 1);
  //T->SetBranchStatus("nBtagMCorrDown", 1);
  T->SetBranchStatus("nBtagTCorr", 1);
  T->SetBranchStatus("nBtagTCorrUp", 1);
  T->SetBranchStatus("nBtagTCorrDown", 1);
};

bool TTJets::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches                                                                                                         
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets", 1);
  T->SetBranchStatus("CSV", 1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  //T->SetBranchStatus("nBtagMCorr", 1);
  //T->SetBranchStatus("nBtagMCorrUp", 1);
  //T->SetBranchStatus("nBtagMCorrDown", 1);
  T->SetBranchStatus("nBtagTCorr", 1);
  T->SetBranchStatus("nBtagTCorrUp", 1);
  T->SetBranchStatus("nBtagTCorrDown", 1);
};

bool TTJets::SetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches
  T1->SetBranchStatus("pu_w",1);
  T1->SetBranchStatus("mu_w",1);
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("N_Jets", 1);
  T1->SetBranchStatus("CSV", 1);
  T1->SetBranchStatus("pTHem1", 1);
  T1->SetBranchStatus("pTHem2", 1);
  T1->SetBranchStatus("etaHem1", 1);
  T1->SetBranchStatus("etaHem2", 1);
  T1->SetBranchStatus("phiHem1", 1);
  T1->SetBranchStatus("phiHem2", 1);
  T1->SetBranchStatus("nBtagLCorr", 1);
  T1->SetBranchStatus("nBtagLCorrUp", 1);
  T1->SetBranchStatus("nBtagLCorrDown", 1);
  //T1->SetBranchStatus("nBtagMCorr", 1);
  //T1->SetBranchStatus("nBtagMCorrUp", 1);
  //T1->SetBranchStatus("nBtagMCorrDown", 1);
  T1->SetBranchStatus("nBtagTCorr", 1);
  T1->SetBranchStatus("nBtagTCorrUp", 1);
  T1->SetBranchStatus("nBtagTCorrDown", 1);
};

bool TTJets::SetMetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches
  T1->SetBranchStatus("pu_w",1);
  T1->SetBranchStatus("mu_w",1);
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("N_Jets", 1);
  T1->SetBranchStatus("CSV", 1);
  T1->SetBranchStatus("ht",1);
  T1->SetBranchStatus("metX",1);
  T1->SetBranchStatus("metY",1);
  T1->SetBranchStatus("metCorrX",1);
  T1->SetBranchStatus("metCorrY",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("nBtagLCorr", 1);
  T1->SetBranchStatus("nBtagLCorrUp", 1);
  T1->SetBranchStatus("nBtagLCorrDown", 1);
  //T1->SetBranchStatus("nBtagMCorr", 1);
  //T1->SetBranchStatus("nBtagMCorrUp", 1);
  //T1->SetBranchStatus("nBtagMCorrDown", 1);
  T1->SetBranchStatus("nBtagTCorr", 1);
  T1->SetBranchStatus("nBtagTCorrUp", 1);
  T1->SetBranchStatus("nBtagTCorrDown", 1);
};

bool TTJets::SetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches
  T2->SetBranchStatus("pu_w",1);
  T2->SetBranchStatus("mu_w",1);
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("N_Jets", 1);
  T2->SetBranchStatus("CSV", 1);
  T2->SetBranchStatus("pTHem1", 1);
  T2->SetBranchStatus("pTHem2", 1);
  T2->SetBranchStatus("etaHem1", 1);
  T2->SetBranchStatus("etaHem2", 1);
  T2->SetBranchStatus("phiHem1", 1);
  T2->SetBranchStatus("phiHem2", 1);
  T2->SetBranchStatus("nBtagLCorr", 1);
  T2->SetBranchStatus("nBtagLCorrUp", 1);
  T2->SetBranchStatus("nBtagLCorrDown", 1);
  //T2->SetBranchStatus("nBtagMCorr", 1);
  //T2->SetBranchStatus("nBtagMCorrUp", 1);
  //T2->SetBranchStatus("nBtagMCorrDown", 1);
  T2->SetBranchStatus("nBtagTCorr", 1);
  T2->SetBranchStatus("nBtagTCorrUp", 1);
  T2->SetBranchStatus("nBtagTCorrDown", 1);

};

bool TTJets::SetMetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches                                                                                                        
  T2->SetBranchStatus("pu_w",1);
  T2->SetBranchStatus("mu_w",1);
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("N_Jets", 1);
  T2->SetBranchStatus("CSV", 1);
  T2->SetBranchStatus("ht",1);
  T2->SetBranchStatus("metX",1);
  T2->SetBranchStatus("metY",1);
  T2->SetBranchStatus("metCorrX",1);
  T2->SetBranchStatus("metCorrY",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("nBtagLCorr", 1);
  T2->SetBranchStatus("nBtagLCorrUp", 1);
  T2->SetBranchStatus("nBtagLCorrDown", 1);
  //T2->SetBranchStatus("nBtagMCorr", 1);
  //T2->SetBranchStatus("nBtagMCorrUp", 1);
  //T2->SetBranchStatus("nBtagMCorrDown", 1);
  T2->SetBranchStatus("nBtagTCorr", 1);
  T2->SetBranchStatus("nBtagTCorrUp", 1);
  T2->SetBranchStatus("nBtagTCorrDown", 1);
};

double TTJets::HLTscale(double MR, double R2){

  int MRbin = -1;
  int R2bin = -1;
  
  const double R2A[] = {0.3, 0.4, 0.5, 0.6, 2.5};
  const double MRA[] = {200., 300., 400., 3500.};

  int Nbins = 3;
  int NbinsR2 = 4;
  
  for(int j = 0; j <= NbinsR2; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
        break;
      }
    }
  }

  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }

  return eff->GetEfficiency(eff->GetGlobalBin(MRbin, R2bin , 0));
};

double TTJets::HLTscaleEle(double MR, double R2){
  
  int MRbin = -1;
  int R2bin = -1;
  
  const double R2A[] = {0.3, 0.4, 0.5, 0.6, 2.5};
  const double MRA[] = {200., 300., 400., 3500.};

  int Nbins = 3;
  int NbinsR2 = 4;
  
  for(int j = 0; j <= NbinsR2; j++){
    if( R2 > R2A[j]){
      if(R2 < R2A[j + 1]){
        R2bin = j+1;
        break;
      }
    }
  }

  for(int j = 0; j <= Nbins; j++){
    if( MR > MRA[j]){
      if(MR < MRA[j + 1]){
        MRbin = j+1;
        break;
      }
    }
  }
  
  return eff_ele->GetEfficiency(eff_ele->GetGlobalBin(MRbin, R2bin , 0));
  
};

