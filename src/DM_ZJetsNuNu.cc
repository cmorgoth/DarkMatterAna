#include "DM_ZJetsNuNu.hh"
#include <iostream>
#include "TLorentzVector.h"

ZJetsNuNu::ZJetsNuNu(){ };

ZJetsNuNu::ZJetsNuNu(int metIndex){
  
  MRMin = 200.0;//nominal 150.0
  RSQMin = 0.5;//nominal 0.5
  this->metIndex = metIndex;
  //std::cout << "+++++++++++++++++++++++++++" << std::endl;
  double N_In = 0 , Ntot = 0;
  
  double Nexp = 0;
  
  TFile* file = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_DoubleMuonPD_Final.root");
  eff = (TEfficiency*)file->Get("Eff2d");
  
  TFile* file1 = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_SignleElePD_Final.root");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
  
  //////////////////////////////////////////////////
  ////////////////MC Files/////////////////////////
  ////////////////////////////////////////////////
  F = TFile::Open("/media/data/cmorgoth/Data/DMData/ZJetsToNuNu_ILV_50_HT_100.root");
  
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
  
  F1 = TFile::Open("/media/data/cmorgoth/Data/DMData/ZJetsToNuNu_ILV_100_HT_200.root");
  
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
  
  F2 = TFile::Open("/media/data/cmorgoth/Data/DMData/ZJetsToNuNu_ILV_200_HT_400.root");
  
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
  
  F3 = TFile::Open("/media/data/cmorgoth/Data/DMData/ZJetsToNuNu_ILV_400_HT_inf.root");
  
  effT = (TTree*)F3->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  T3 = (TTree*)F3->Get("outTree");
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  weight3 = Lumi*sigma3/Ntot;
  Nexp += weight3*T3->GetEntries(); 
  //delete F;
  //delete effT;
  //std::cout << "Expected Number of events ZJetsNuNu: " << Nexp << std::endl;
  std::cout << "Weight 0: " << weight0 << std::endl;
  std::cout << "Weight 1: " << weight1 << std::endl;
  std::cout << "Weight 2: " << weight2 << std::endl;
  std::cout << "Weight 3: " << weight3 << std::endl;
  
  ////////////////////////////////////////////////////////////////
  ////////// Including Btag capability///////////////////////////
  ///////////////////////////////////////////////////////////////

  /*
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
  */
};

ZJetsNuNu::ZJetsNuNu(const char* FileName ){
  
  //F = new TFile(FileName);
  //T0 = (TTree*)F->Get("outTree");
  
};

ZJetsNuNu::~ZJetsNuNu(){
  delete T;
  delete T1;
  delete T2;
  delete T3;
  delete F;
  delete F1;
  delete F2;
  delete F3;
};


bool ZJetsNuNu::PrintEvents(){
  
  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;
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
  
  //second HTbin
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
  
  //Thrid HTbin
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

  //Fourth HTbin
  effT = (TTree*)F3->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In*weight3;
    Nt_PV += N_PV*weight3;
    Nt_2J += N_2J*weight3;
    Nt_0b += N_0b*weight3;
    Nt_LepVeto += N_LepVeto*weight3;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "================ZJets===================" << std::endl;
  std::cout << "========================================" << std::endl;

  std::cout << "Znunu weighted Nt_In: " << NtotGen << std::endl;
  std::cout << "Znunu weighted Nt_PV: " << Nt_PV << std::endl;
  std::cout << "Znunu weighted Nt_2J: " << Nt_2J << std::endl;
  std::cout << "Znunu weighted Nt_0b: " << Nt_0b << std::endl;
  std::cout << "Znunu weighted Nt_LepVeto: " << Nt_LepVeto << std::endl;

  std::cout << "Znunu weighted Nt_tree: " << T->GetEntries()*weight0 +T1->GetEntries()*weight1 +  T2->GetEntries()*weight2 + T3->GetEntries()*weight3 << "\n\n" << std::endl;

  //std::cout << "========================================" << std::endl;

  //After Selection cuts
  double RSQ[4], MR[4];
  int BOX, nBtag;
  
  double Nt_MR_RSQ_cut0BTag = 0.0, Nt_2muBox0BTag = 0.0,  Nt_1muBox0BTag = 0.0, Nt_0muBox0BTag = 0.0,\
    Nt_MR_RSQ_cut = 0.0, Nt_2muBox = 0.0,  Nt_1muBox = 0.0, Nt_0muBox = 0.0, Nt_Btags = 0.0;

  //std::cout << "weights: " << weight0 << " " << weight1 << " " << weight2 << " " << weight3 << std::endl; 
  
  SetStatus();
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += weight3;
      
      if( BOX == 0 )Nt_0muBox += weight3;
      if( BOX == 1 )Nt_1muBox += weight3;
      if( BOX == 2 )Nt_2muBox += weight3;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag += weight3;
	
	if( BOX == 0 )Nt_0muBox0BTag += weight3;
	if( BOX == 1 )Nt_1muBox0BTag += weight3;
	if( BOX == 2 )Nt_2muBox0BTag += weight3;
	
      }
      
    }
    
  }
  T3->SetBranchStatus("*", 0);

  std::cout << "Znunu weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut << std::endl;
  std::cout << "Znunu weighted Nt_0muBox: " << Nt_0muBox << std::endl;
  std::cout << "Znunu weighted Nt_1muBox: " << Nt_1muBox << std::endl;
  std::cout << "Znunu weighted Nt_2muBox: " << Nt_2muBox << std::endl;

  std::cout << "Znunu weighted Nt_MR_RSQ_cut0Btag: " << Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "Znunu weighted Nt_0muBox0Btag: " << Nt_0muBox0BTag << std::endl;
  std::cout << "Znunu weighted Nt_1muBox0Btag: " << Nt_1muBox0BTag << std::endl;
  std::cout << "Znunu weighted Nt_2muBox0Btag: " << Nt_2muBox0BTag << std::endl;

  std::cout << "Znunu Btag EVENTS: " << Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
  
};

TH1F ZJetsNuNu::PlotMR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  TH1F* MR1 = new TH1F("MR1", "MR1BOX", MR_Bins, MR_BinArr);
  MR1->Sumw2();
  
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
      if(RSQ[metIndex]<0.5)std::cout << "HLT: " << hltWeight << std::endl;
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
  std::cout << "debug0" << std::endl;
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
  std::cout << "debug1"<< std::endl;
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
  std::cout << "debug2"<< std::endl;
  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);
  std::cout << "debug3"<< std::endl;
  return *MR1;
};

TH1F ZJetsNuNu::PlotMR_0Box(){
  
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);

  return *MR0;
};


TH1F  ZJetsNuNu::PlotRSQ_1Box(){
  
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);

  return *RSQ1;
  
};


TH1F  ZJetsNuNu::PlotRSQ_0Box(){
  
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);

  return *RSQ0;
  
};


TH2F ZJetsNuNu::PlotRSQ_vs_MR_0Box(){
  
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F ZJetsNuNu::PlotRSQ_vs_MR_1Box(){
  
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

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], weight3*hltWeight);
    }
  }
  T3->SetBranchStatus("*", 0);

  return *RSQ_MR_1BOX;
  
};

std::vector<TH1F*> ZJetsNuNu::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, nBtag[2], N_Jets;

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
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", &mht[0]);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  
  float metmag = 0.;
  float metmagcorr = 0.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[0]->Fill(metmag, weight0*hltWeight);
        MET[1]->Fill(metmagcorr, weight0*hltWeight);
        MET[2]->Fill(metmagcorr-metmag, weight0*hltWeight);
	MET[9]->Fill(N_Jets,weight0*hltWeight);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight0*hltWeight);
        MET[4]->Fill(metmagcorr, weight0*hltWeight);
        MET[5]->Fill(metmagcorr-metmag, weight0*hltWeight);
	MET[10]->Fill(N_Jets,weight0*hltWeight);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight0*hltWeight);
        MET[7]->Fill(metmagcorr, weight0*hltWeight);
        MET[8]->Fill(metmagcorr-metmag, weight0*hltWeight);
	MET[11]->Fill(N_Jets,weight0*hltWeight);
      }
    }
    
  }
  T->SetBranchStatus("*", 0);

  SetMetStatus1();
  T1->SetBranchAddress("RSQ", RSQ);
  T1->SetBranchAddress("MR", MR);
  T1->SetBranchAddress("BOX_NUM", &BOX);
  T1->SetBranchAddress("nBtag", &nBtag[0]);
  T1->SetBranchAddress("nBtagTight", &nBtag[1]);
  T1->SetBranchAddress("N_Jets", &N_Jets);
  T1->SetBranchAddress("CSV", CSV);
  T1->SetBranchAddress("ht", &ht);
  T1->SetBranchAddress("mht", &mht[0]);
  T1->SetBranchAddress("metX", metX);
  T1->SetBranchAddress("metCorrX", metcorrX);
  T1->SetBranchAddress("metY", metY);
  T1->SetBranchAddress("metCorrY", metcorrY);
  
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[0]->Fill(metmag, weight1*hltWeight);
        MET[1]->Fill(metmagcorr, weight1*hltWeight);
        MET[2]->Fill(metmagcorr-metmag, weight1*hltWeight);
	MET[9]->Fill(N_Jets,weight1*hltWeight);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight1*hltWeight);
        MET[4]->Fill(metmagcorr, weight1*hltWeight);
        MET[5]->Fill(metmagcorr-metmag, weight1*hltWeight);
	MET[10]->Fill(N_Jets,weight1*hltWeight);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight1*hltWeight);
        MET[7]->Fill(metmagcorr, weight1*hltWeight);
        MET[8]->Fill(metmagcorr-metmag, weight1*hltWeight);
	MET[11]->Fill(N_Jets,weight1*hltWeight);
      }
    }
    
  }
  T1->SetBranchStatus("*", 0);

  SetMetStatus2();
  T2->SetBranchAddress("RSQ", RSQ);
  T2->SetBranchAddress("MR", MR);
  T2->SetBranchAddress("BOX_NUM", &BOX);
  T2->SetBranchAddress("nBtag", &nBtag[0]);
  T2->SetBranchAddress("nBtagTight", &nBtag[1]);
  T2->SetBranchAddress("N_Jets", &N_Jets);
  T2->SetBranchAddress("CSV", CSV);
  T2->SetBranchAddress("ht", &ht);
  T2->SetBranchAddress("mht", &mht[0]);
  T2->SetBranchAddress("metX", metX);
  T2->SetBranchAddress("metCorrX", metcorrX);
  T2->SetBranchAddress("metY", metY);
  T2->SetBranchAddress("metCorrY", metcorrY);
    
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[0]->Fill(metmag, weight2*hltWeight);
        MET[1]->Fill(metmagcorr, weight2*hltWeight);
        MET[2]->Fill(metmagcorr-metmag, weight2*hltWeight);
	MET[9]->Fill(N_Jets,weight2*hltWeight);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight2*hltWeight);
        MET[4]->Fill(metmagcorr, weight2*hltWeight);
        MET[5]->Fill(metmagcorr-metmag, weight2*hltWeight);
	MET[10]->Fill(N_Jets,weight2*hltWeight);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight2*hltWeight);
        MET[7]->Fill(metmagcorr, weight2*hltWeight);
	MET[8]->Fill(metmagcorr-metmag, weight2*hltWeight);
	MET[11]->Fill(N_Jets,weight2*hltWeight);
      }
    }
  }
  T2->SetBranchStatus("*", 0);

  SetMetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  T3->SetBranchAddress("ht", &ht);
  T3->SetBranchAddress("mht", &mht[0]);
  T3->SetBranchAddress("metX", metX);
  T3->SetBranchAddress("metCorrX", metcorrX);
  T3->SetBranchAddress("metY", metY);
  T3->SetBranchAddress("metCorrY", metcorrY);
    
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[0]->Fill(metmag, weight3*hltWeight);
        MET[1]->Fill(metmagcorr, weight3*hltWeight);
        MET[2]->Fill(metmagcorr-metmag, weight3*hltWeight);
	MET[9]->Fill(N_Jets,weight3*hltWeight);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[3]->Fill(metmag, weight3*hltWeight);
        MET[4]->Fill(metmagcorr, weight3*hltWeight);
        MET[5]->Fill(metmagcorr-metmag, weight3*hltWeight);
	MET[10]->Fill(N_Jets,weight3*hltWeight);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        MET[6]->Fill(metmag, weight3*hltWeight);
        MET[7]->Fill(metmagcorr, weight3*hltWeight);
        MET[8]->Fill(metmagcorr-metmag, weight3*hltWeight);
	MET[11]->Fill(N_Jets,weight3*hltWeight);
      }
    }
    
  }
  T3->SetBranchStatus("*", 0);

  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec;
  
};


bool ZJetsNuNu::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int ZJetsNuNu::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH2F*> ZJetsNuNu::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int BOX, N_Jets, nBtag[2];
  
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
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight0*hltWeight);
      }
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5 ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight1*hltWeight);
      }
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight2*hltWeight);
      }
    }
  }
  T2->SetBranchStatus("*", 0);

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  T3->SetBranchAddress("pTHem1", &pTHem1);
  T3->SetBranchAddress("pTHem2", &pTHem2);
  T3->SetBranchAddress("etaHem1", &etaHem1);
  T3->SetBranchAddress("etaHem2", &etaHem2);
  T3->SetBranchAddress("phiHem1", &phiHem1);
  T3->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);

    TLorentzVector j1;
    TLorentzVector j2;
    
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight3*hltWeight);
      }
    }
  }
  T3->SetBranchStatus("*", 0);

  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;

};


std::vector<TH1F*> ZJetsNuNu::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int BOX, N_Jets,nBtag[2];

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 6; l++ ){
    name = TString(Form("MR_1D_Z_%dmu_Box",l));
    name1 = TString(Form("R2_1D_Z_%dmu_Box",l));
    Razor1D[2*l] = new TH1F( name, name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F( name1, name1, RSQ_Bins, RSQ_BinArr);
    if( l < 3 ){
      Razor1D[2*l]->Sumw2();
      Razor1D[2*l+1]->Sumw2();
    }
  }

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight0*hltWeight);
        Razor1D[1]->Fill(RSQ[metIndex], weight0*hltWeight);
	Razor1D[6]->Fill(MR[metIndex], weight0*hltWeight);
        Razor1D[7]->Fill(RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex],  weight0*hltWeight);
        Razor1D[3]->Fill(RSQ[metIndex], weight0*hltWeight);
	Razor1D[8]->Fill(MR[metIndex], weight0*hltWeight);
        Razor1D[9]->Fill(RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight0*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight0*hltWeight);
	Razor1D[10]->Fill(MR[metIndex], weight0*hltWeight);
        Razor1D[11]->Fill(RSQ[metIndex], weight0*hltWeight);
      }
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight1*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight1*hltWeight);
	Razor1D[6]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[7]->Fill(RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[3]->Fill(RSQ[metIndex], weight1*hltWeight);
	Razor1D[8]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[9]->Fill(RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight1*hltWeight);
	Razor1D[10]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[11]->Fill(RSQ[metIndex], weight1*hltWeight);
      }
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
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[1]->Fill(RSQ[metIndex], weight2*hltWeight);
	Razor1D[6]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[7]->Fill(RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[3]->Fill(RSQ[metIndex], weight2*hltWeight);
	Razor1D[8]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[9]->Fill(RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight2*hltWeight);
	Razor1D[10]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[11]->Fill(RSQ[metIndex], weight2*hltWeight);
      }
    }
  }
  T2->SetBranchStatus("*", 0);

  SetStatus3();
  T3->SetBranchAddress("RSQ", RSQ);
  T3->SetBranchAddress("MR", MR);
  T3->SetBranchAddress("BOX_NUM", &BOX);
  T3->SetBranchAddress("nBtag", &nBtag[0]);
  T3->SetBranchAddress("nBtagTight", &nBtag[1]);
  T3->SetBranchAddress("N_Jets", &N_Jets);
  T3->SetBranchAddress("CSV", CSV);
  T3->SetBranchAddress("pTHem1", &pTHem1);
  T3->SetBranchAddress("pTHem2", &pTHem2);
  T3->SetBranchAddress("etaHem1", &etaHem1);
  T3->SetBranchAddress("etaHem2", &etaHem2);
  T3->SetBranchAddress("phiHem1", &phiHem1);
  T3->SetBranchAddress("phiHem2", &phiHem2);
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    TLorentzVector j1;
    TLorentzVector j2;

    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && Dphi < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[1]->Fill(RSQ[metIndex], weight3*hltWeight);
	Razor1D[6]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[7]->Fill(RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[3]->Fill(RSQ[metIndex], weight3*hltWeight);
	Razor1D[8]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[9]->Fill(RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight3*hltWeight);
	Razor1D[10]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[11]->Fill(RSQ[metIndex], weight3*hltWeight);
      }
    }
  }
  
  T3->SetBranchStatus("*", 0);
  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }
  
  return Razor1DVec;
  
};

bool ZJetsNuNu::SetStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem2", 1);
};

bool ZJetsNuNu::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
};
  
bool ZJetsNuNu::SetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("CSV",1);
  T1->SetBranchStatus("pTHem1", 1);
  T1->SetBranchStatus("pTHem2", 1);
  T1->SetBranchStatus("etaHem1", 1);
  T1->SetBranchStatus("etaHem2", 1);
  T1->SetBranchStatus("phiHem1", 1);
  T1->SetBranchStatus("phiHem2", 1);
};

bool ZJetsNuNu::SetMetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("ht",1);
  T1->SetBranchStatus("metX",1);
  T1->SetBranchStatus("metY",1);
  T1->SetBranchStatus("metCorrX",1);
  T1->SetBranchStatus("metCorrY",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("CSV",1);
};

bool ZJetsNuNu::SetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("CSV",1);
  T2->SetBranchStatus("pTHem1", 1);
  T2->SetBranchStatus("pTHem2", 1);
  T2->SetBranchStatus("etaHem1", 1);
  T2->SetBranchStatus("etaHem2", 1);
  T2->SetBranchStatus("phiHem1", 1);
  T2->SetBranchStatus("phiHem2", 1);
};

bool ZJetsNuNu::SetMetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("ht",1);
  T2->SetBranchStatus("metX",1);
  T2->SetBranchStatus("metY",1);
  T2->SetBranchStatus("metCorrX",1);
  T2->SetBranchStatus("metCorrY",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("CSV",1);
};

bool ZJetsNuNu::SetStatus3(){
  T3->SetBranchStatus("*",0); //disable all branches
  T3->SetBranchStatus("RSQ",1);
  T3->SetBranchStatus("MR",1);
  T3->SetBranchStatus("BOX_NUM",1);
  T3->SetBranchStatus("nBtag",1);
  T3->SetBranchStatus("nBtagTight",1);
  T3->SetBranchStatus("N_Jets",1);
  T3->SetBranchStatus("CSV",1);
  T3->SetBranchStatus("pTHem1", 1);
  T3->SetBranchStatus("pTHem2", 1);
  T3->SetBranchStatus("etaHem1", 1);
  T3->SetBranchStatus("etaHem2", 1);
  T3->SetBranchStatus("phiHem1", 1);
  T3->SetBranchStatus("phiHem2", 1);
};

bool ZJetsNuNu::SetMetStatus3(){
  T3->SetBranchStatus("*",0); //disable all branches
  T3->SetBranchStatus("RSQ",1);
  T3->SetBranchStatus("MR",1);
  T3->SetBranchStatus("BOX_NUM",1);
  T3->SetBranchStatus("nBtag",1);
  T3->SetBranchStatus("nBtagTight",1);
  T3->SetBranchStatus("ht",1);
  T3->SetBranchStatus("metX",1);
  T3->SetBranchStatus("metY",1);
  T3->SetBranchStatus("metCorrX",1);
  T3->SetBranchStatus("metCorrY",1);
  T3->SetBranchStatus("N_Jets",1);
  T3->SetBranchStatus("CSV",1);
};

double ZJetsNuNu::HLTscale(double MR, double R2){

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

double ZJetsNuNu::HLTscaleEle(double MR, double R2){
  
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


