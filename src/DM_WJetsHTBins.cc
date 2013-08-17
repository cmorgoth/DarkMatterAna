#include "DM_WJetsHTBins.hh"
#include <iostream>

WJetsHTBins::WJetsHTBins(){ };

WJetsHTBins::WJetsHTBins(int metIndex){
  
  MRMin = 200.0;//nominal 150.0
  RSQMin = 0.5;//nominal 0.5
  this->metIndex = metIndex;
  
  TFile* file = new TFile("/media/data/cmorgoth/Data/DMData/hlt_eff_mr200_MoreBin_ABCD_PT80v2Muon.root");
  //hlt = (TH2F*)file->Get("h");
  eff = (TEfficiency*)file->Get("Eff2d");
  TFile* file1 = new TFile("/media/data/cmorgoth/Data/DMData/hlt_eff_mr200_MoreBin_PT80.root");
  //hlt_ele = (TH2F*)file1->Get("h");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
  //std::cout << hlt->GetNbinsX() << std::endl;
  //std::cout << hlt_ele->GetNbinsX() << std::endl;
  
  double N_In = 0 , Ntot = 0;
  double Nexp = 0;
  
  F = TFile::Open("/media/data/cmorgoth/Data/DMData/WJetsToLNu_ILV_150_HT_200.root");

  //F = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/WJetsToLNu_150_HT_200_NoTauVeto.root");
  
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
  
  F1 = TFile::Open("/media/data/cmorgoth/Data/DMData/WJetsToLNu_ILV_200_HT_250.root");
  
  //F1 = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/WJetsToLNu_200_HT_250_NoTauVeto.root");
  
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
  
  F2 = TFile::Open("/media/data/cmorgoth/Data/DMData/WJetsToLNu_ILV_250_HT_300.root");
  
  //F2 =  TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/WJetsToLNu_250_HT_300_NoTauVeto.root");
  
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

  F3 = TFile::Open("/media/data/cmorgoth/Data/DMData/WJetsToLNu_ILV_300_HT_400.root");
    
  //F3 = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/WJetsToLNu_300_HT_400_Total_NoTauVeto.root");

  effT = (TTree*)F3->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  T3 = (TTree*)F3->Get("outTree");
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
    //std::cout << "N_IN3: " << N_In << std::endl;
  }
  //std::cout << "Ntot HT400-inf: " << Ntot << std::endl;
  weight3 = Lumi*sigma3/Ntot;
  Nexp += weight3*T3->GetEntries(); 
  //delete F;
  //delete effT;
  //std::cout << "Expected Number of events WJetsHTBins: " << Nexp << std::endl;

  F4 = TFile::Open("/media/data/cmorgoth/Data/DMData/WJetsToLNu_ILV_400_HT_Inf.root");
  
  //F4 = TFile::Open("root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/WJetsToLNu_400_HT_Inf_NoTauVeto.root");
  
  effT = (TTree*)F4->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  T4 = (TTree*)F4->Get("outTree");
  Ntot = 0;
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
    //std::cout << "N_IN3: " << N_In << std::endl;                                                                 
  }
  //std::cout << "Ntot HT400-inf: " << Ntot << std::endl;                                                          
  weight4 = Lumi*sigma4/Ntot;
  Nexp += weight4*T4->GetEntries();
  
  
  std::cout << "Weight 0: " << weight0 << std::endl;
  std::cout << "Weight 1: " << weight1 << std::endl;
  std::cout << "Weight 2: " << weight2 << std::endl;
  std::cout << "Weight 3: " << weight3 << std::endl;
  std::cout << "Weight 4: " << weight4 << std::endl;
  
  //////////////////////////////////////////////////////////////// 
  ////////// Including Btag capability///////////////////////////                                   
  ///////////////////////////////////////////////////////////////

  
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
  std::cout << "------------------BranchIndex: " << this->btagIndex << std::endl;

};


WJetsHTBins::WJetsHTBins(const char* FileName ){
  
  //F = new TFile(FileName);
  //T0 = (TTree*)F->Get("outTree");
  
};

WJetsHTBins::~WJetsHTBins(){
  delete T;
  delete T1;
  delete T2;
  delete T3;
  delete T4;
  delete F;
  delete F1;
  delete F2;
  delete F3;
  delete F4;
};


bool WJetsHTBins::PrintEvents(){

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
  
  //Fifth HTbin                                                                                                   
  effT = (TTree*)F4->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In*weight4;
    Nt_PV += N_PV*weight4;
    Nt_2J += N_2J*weight4;
    Nt_0b += N_0b*weight4;
    Nt_LepVeto += N_LepVeto*weight4;
  }
  

  std::cout << "========================================" << std::endl;
  std::cout << "================ W + Jets ==============" << std::endl;
  std::cout << "========================================" << std::endl;

  std::cout << "W+jets weighted Nt_In: " << NtotGen << std::endl;
  std::cout << "W+jets weighted Nt_PV: " << Nt_PV << std::endl;
  std::cout << "W+jets weighted Nt_2J: " << Nt_2J << std::endl;
  std::cout << "W+jets weighted Nt_0b: " << Nt_0b << std::endl;
  std::cout << "W+jets weighted Nt_LepVeto: " << Nt_LepVeto << std::endl;

  std::cout << "W+jets weighted Nt_tree: " << T->GetEntries()*weight0 +T1->GetEntries()*weight1 +  T2->GetEntries()*weight2 + T3->GetEntries()*weight3 + T4->GetEntries()*weight4 << "\n\n" << std::endl;

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
  T->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

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
  T1->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance
  
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
  T2->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

  //fourth HT bin
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
  T3->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

  //fifth HT bin 
  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag);
  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += weight4;
      
      if( BOX == 0 )Nt_0muBox += weight4;
      if( BOX == 1 )Nt_1muBox += weight4;
      if( BOX == 2 )Nt_2muBox += weight4;
      
      if( nBtag == 0 ){
	
        Nt_MR_RSQ_cut0BTag += weight4;
	
        if( BOX == 0 )Nt_0muBox0BTag += weight4;
        if( BOX == 1 )Nt_1muBox0BTag += weight4;
        if( BOX == 2 )Nt_2muBox0BTag += weight4;
	
      }
      
    }
    
  }

  T4->SetBranchStatus("*", 0);

  //std::cout << "========================================" << std::endl;
  std::cout << "WJets weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut << std::endl;
  std::cout << "WJets weighted Nt_0muBox: " << Nt_0muBox << std::endl;
  std::cout << "WJets weighted Nt_1muBox: " << Nt_1muBox << std::endl;
  std::cout << "WJets weighted Nt_2muBox: " << Nt_2muBox << std::endl;

  std::cout << "WJets weighted Nt_MR_RSQ_cut0Btag: " << Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "WJets weighted Nt_0muBox0Btag: " << Nt_0muBox0BTag << std::endl;
  std::cout << "WJets weighted Nt_1muBox0Btag: " << Nt_1muBox0BTag << std::endl;
  std::cout << "WJets weighted Nt_2muBox0Btag: " << Nt_2muBox0BTag << std::endl;

  std::cout << "WJets Btag EVENTS: " << Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
  
};

TH1F WJetsHTBins::PlotMR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0;
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
    fBtag[4] = (nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight0*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
  //////////////////////////
  
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
  ///////////////////////////
  
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
  ////////////////////////////
  
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
  ////////////////////////////
  
  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR1->Fill(MR[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);
  
  return *MR1;
};

TH1F WJetsHTBins::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0;
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
  ////
  
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      MR0->Fill(MR[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);
  
  
  return *MR0;
};


TH1F  WJetsHTBins::PlotRSQ_1Box(){
  
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ1->Fill(RSQ[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);
  
  return *RSQ1;
    
};


TH1F  WJetsHTBins::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  //TH1F* RSQ0 = new TH1F("RSQ0", "RSQ0BOX", RSQ_Bins, RSQ_Low, RSQ_High);
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ0->Fill(RSQ[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);

  return *RSQ0;
  
};


TH2F WJetsHTBins::PlotRSQ_vs_MR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* RSQ_MR_0BOX = new TH2F("RSQ_MR_0BOX_W", "RSQ_VS_MR_0BOX_W", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F WJetsHTBins::PlotRSQ_vs_MR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* RSQ_MR_1BOX = new TH2F("RSQ_MR_1BOX_W", "RSQ_VS_MR_1BOX_W", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex]){
      hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
      if( hltWeight == 0.0 )hltWeight = 1.0;
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], weight4*hltWeight);
    }
  }
  T4->SetBranchStatus("*", 0);

  return *RSQ_MR_1BOX;
  
};


std::vector<TH1F*> WJetsHTBins::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, nBtag[2], N_Jets;

  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, 0, 1000 );
    }
  }
  
  MET[9] = new TH1F( "NJETS0_W", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( "NJETS1_W", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( "NJETS2_W", "NJETS 2 BOX", 9, 1, 10);
  
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
  T->SetBranchAddress("N_Jets", &N_Jets);
  
  float metmag = .0;
  float metmagcorr = .0;
  double wt = 1.;
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
      if( BOX == 0 ){
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, weight0*wt);
	MET[1]->Fill(metmagcorr, weight0*wt);
	MET[2]->Fill(metmagcorr-metmag, weight0*wt);
	MET[9]->Fill(N_Jets,weight0*wt);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, weight0*wt);
	MET[4]->Fill(metmagcorr, weight0*wt);
	MET[5]->Fill(metmagcorr-metmag, weight0*wt);
	MET[10]->Fill(N_Jets, weight0*wt);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, weight0*wt);
	MET[7]->Fill(metmagcorr, weight0*wt);
	MET[8]->Fill(metmagcorr-metmag, weight0*wt);
	MET[11]->Fill(N_Jets, weight0*wt);
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
  T1->SetBranchAddress("N_Jets", &N_Jets);
  
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
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, weight1*wt);
        MET[1]->Fill(metmagcorr, weight1*wt);
        MET[2]->Fill(metmagcorr-metmag, weight1*wt);
	MET[9]->Fill(N_Jets,weight1*wt);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, weight1*wt);
        MET[4]->Fill(metmagcorr, weight1*wt);
        MET[5]->Fill(metmagcorr-metmag, weight1*wt);
	MET[10]->Fill(N_Jets,weight1*wt);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, weight1*wt);
        MET[7]->Fill(metmagcorr, weight1*wt);
        MET[8]->Fill(metmagcorr-metmag, weight1*wt);
	MET[11]->Fill(N_Jets,weight1*wt);
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
  T2->SetBranchAddress("N_Jets", &N_Jets);

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
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, weight2*wt);
        MET[1]->Fill(metmagcorr, weight2*wt);
        MET[2]->Fill(metmagcorr-metmag, weight2*wt);
	MET[9]->Fill(N_Jets,weight2*wt);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, weight2*wt);
        MET[4]->Fill(metmagcorr, weight2*wt);
        MET[5]->Fill(metmagcorr-metmag, weight2*wt);
	MET[10]->Fill(N_Jets,weight2*wt);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, weight2*wt);
        MET[7]->Fill(metmagcorr, weight2*wt);
	MET[8]->Fill(metmagcorr-metmag, weight2*wt);
	MET[11]->Fill(N_Jets,weight2*wt);
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
  T3->SetBranchAddress("N_Jets", &N_Jets);

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
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, weight3*wt);
        MET[1]->Fill(metmagcorr, weight3*wt);
        MET[2]->Fill(metmagcorr-metmag, weight3*wt);
	MET[9]->Fill(N_Jets,weight3*wt);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, weight3*wt);
        MET[4]->Fill(metmagcorr, weight3*wt);
        MET[5]->Fill(metmagcorr-metmag, weight3*wt);
	MET[10]->Fill(N_Jets,weight3*wt);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, weight3*wt);
        MET[7]->Fill(metmagcorr, weight3*wt);
	MET[8]->Fill(metmagcorr-metmag, weight3*wt);
	MET[11]->Fill(N_Jets,weight3*wt);
      }
    }

  }
  T3->SetBranchStatus("*", 0);

  SetMetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);
  T4->SetBranchAddress("ht", &ht);
  T4->SetBranchAddress("mht", &mht[0]);
  T4->SetBranchAddress("metX", metX);
  T4->SetBranchAddress("metCorrX", metcorrX);
  T4->SetBranchAddress("metY", metY);
  T4->SetBranchAddress("metCorrY", metcorrY);
  T4->SetBranchAddress("N_Jets", &N_Jets);

  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, weight4*wt);
        MET[1]->Fill(metmagcorr, weight4*wt);
        MET[2]->Fill(metmagcorr-metmag, weight4*wt);
	MET[9]->Fill(N_Jets,weight4*wt);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, weight4*wt);
        MET[4]->Fill(metmagcorr, weight4*wt);
        MET[5]->Fill(metmagcorr-metmag, weight4*wt);
	MET[10]->Fill(N_Jets,weight4*wt);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, weight4*wt);
        MET[7]->Fill(metmagcorr, weight4*wt);
	MET[8]->Fill(metmagcorr-metmag, weight4*wt);
	MET[11]->Fill(N_Jets,weight4*wt);
      }
    }
    
  }
  T4->SetBranchStatus("*", 0);

  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }

  return metvec;
};

std::vector<TH2F*> WJetsHTBins::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  int BOX, N_Jets, nBtag[2];
  
  std::vector< TH2F* > Razor2DVec;
  TH2F* Razor2D[3];
  TString name;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("Razor2D_Wj_%dmu_Box",l));
    Razor2D[l] = new TH2F( name, name, 200, 100., 1500., 200, 0.0, 1.5);
  }

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

    if( /*RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  &&*/ fBtag[btagIndex] ){
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
  
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( /*RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && */fBtag[btagIndex] ){
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

  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( /*RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && */fBtag[btagIndex] ){
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
  
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( /*RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && */fBtag[btagIndex] ){
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

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);
  
  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( /*RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && */fBtag[btagIndex] ){
      if( BOX == 0){
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], weight4*hltWeight);
      }else if( BOX == 1 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], weight4*hltWeight);
      }else if( BOX == 2 ){
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], weight4*hltWeight);
      }
    }

  }
  T4->SetBranchStatus("*", 0);

  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;

};

bool WJetsHTBins::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int WJetsHTBins::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH1F*> WJetsHTBins::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[10];
  int BOX, N_Jets, nBtag[2];

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[6];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("MR_1D_Wj_%dmu_Box",l));
    name1 = TString(Form("R2_1D_Wj_%dmu_Box",l));
    Razor1D[2*l] = new TH1F( name, name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F( name1, name1, RSQ_Bins, RSQ_BinArr);
    Razor1D[2*l]->Sumw2();
    Razor1D[2*l+1]->Sumw2();
  }

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

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight0*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex],  weight0*hltWeight);
	Razor1D[3]->Fill(RSQ[metIndex], weight0*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight0*hltWeight);
	Razor1D[5]->Fill(RSQ[metIndex], weight0*hltWeight);
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
  for(int i = 0; i < T1->GetEntries(); i++){
    T1->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight1*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight1*hltWeight);
	Razor1D[3]->Fill(RSQ[metIndex], weight1*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight1*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight1*hltWeight);
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
  for(int i = 0; i < T2->GetEntries(); i++){
    T2->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight2*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight2*hltWeight);
	Razor1D[3]->Fill(RSQ[metIndex], weight2*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight2*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight2*hltWeight);
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
  for(int i = 0; i < T3->GetEntries(); i++){
    T3->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight3*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight3*hltWeight);
	Razor1D[3]->Fill(RSQ[metIndex], weight3*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight3*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight3*hltWeight);
      }
    }
  }
  T3->SetBranchStatus("*", 0);

  SetStatus4();
  T4->SetBranchAddress("RSQ", RSQ);
  T4->SetBranchAddress("MR", MR);
  T4->SetBranchAddress("BOX_NUM", &BOX);
  T4->SetBranchAddress("nBtag", &nBtag[0]);
  T4->SetBranchAddress("nBtagTight", &nBtag[1]);
  T4->SetBranchAddress("N_Jets", &N_Jets);
  T4->SetBranchAddress("CSV", CSV);
  for(int i = 0; i < T4->GetEntries(); i++){
    T4->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], weight4*hltWeight);
	Razor1D[1]->Fill(RSQ[metIndex], weight4*hltWeight);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex], weight4*hltWeight);
	Razor1D[3]->Fill(RSQ[metIndex], weight4*hltWeight);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], weight4*hltWeight);
        Razor1D[5]->Fill(RSQ[metIndex], weight4*hltWeight);
      }
    }
  }
  T4->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 6; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }

  return Razor1DVec;
  
};


bool WJetsHTBins::SetStatus(){

  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
};

bool WJetsHTBins::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches                                                                                                                                         
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("idMC", 1);
  T->SetBranchStatus("statusMC", 1);
};

bool WJetsHTBins::SetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("CSV",1);
};

bool WJetsHTBins::SetMetStatus1(){
  T1->SetBranchStatus("*",0); //disable all branches                                                                                                                                         
  T1->SetBranchStatus("RSQ",1);
  T1->SetBranchStatus("MR",1);
  T1->SetBranchStatus("BOX_NUM",1);
  T1->SetBranchStatus("nBtag",1);
  T1->SetBranchStatus("nBtagTight",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("CSV",1);
  T1->SetBranchStatus("ht",1);
  T1->SetBranchStatus("metX",1);
  T1->SetBranchStatus("metY",1);
  T1->SetBranchStatus("metCorrX",1);
  T1->SetBranchStatus("metCorrY",1);
  T1->SetBranchStatus("N_Jets",1);
  T1->SetBranchStatus("idMC", 1);
  T1->SetBranchStatus("statusMC", 1);
};

bool WJetsHTBins::SetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("CSV",1);
};

bool WJetsHTBins::SetMetStatus2(){
  T2->SetBranchStatus("*",0); //disable all branches                                                                                                                                         
  T2->SetBranchStatus("RSQ",1);
  T2->SetBranchStatus("MR",1);
  T2->SetBranchStatus("BOX_NUM",1);
  T2->SetBranchStatus("nBtag",1);
  T2->SetBranchStatus("nBtagTight",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("CSV",1);
  T2->SetBranchStatus("ht",1);
  T2->SetBranchStatus("metX",1);
  T2->SetBranchStatus("metY",1);
  T2->SetBranchStatus("metCorrX",1);
  T2->SetBranchStatus("metCorrY",1);
  T2->SetBranchStatus("N_Jets",1);
  T2->SetBranchStatus("idMC", 1);
  T2->SetBranchStatus("statusMC", 1);
};

bool WJetsHTBins::SetStatus3(){

  T3->SetBranchStatus("*",0); //disable all branches
  T3->SetBranchStatus("RSQ",1);
  T3->SetBranchStatus("MR",1);
  T3->SetBranchStatus("BOX_NUM",1);
  T3->SetBranchStatus("nBtag",1);
  T3->SetBranchStatus("nBtagTight",1);
  T3->SetBranchStatus("N_Jets",1);
  T3->SetBranchStatus("CSV",1);
};

bool WJetsHTBins::SetMetStatus3(){
  T3->SetBranchStatus("*",0); //disable all branches                                                                                                                                         
  T3->SetBranchStatus("RSQ",1);
  T3->SetBranchStatus("MR",1);
  T3->SetBranchStatus("BOX_NUM",1);
  T3->SetBranchStatus("nBtag",1);
  T3->SetBranchStatus("nBtagTight",1);
  T3->SetBranchStatus("N_Jets",1);
  T3->SetBranchStatus("CSV",1);
  T3->SetBranchStatus("ht",1);
  T3->SetBranchStatus("metX",1);
  T3->SetBranchStatus("metY",1);
  T3->SetBranchStatus("metCorrX",1);
  T3->SetBranchStatus("metCorrY",1);
  T3->SetBranchStatus("N_Jets",1);
  T3->SetBranchStatus("idMC", 1);
  T3->SetBranchStatus("statusMC", 1);
};

bool WJetsHTBins::SetStatus4(){

  T4->SetBranchStatus("*",0); //disable all branches                                                          
  T4->SetBranchStatus("RSQ",1);
  T4->SetBranchStatus("MR",1);
  T4->SetBranchStatus("BOX_NUM",1);
  T4->SetBranchStatus("nBtag",1);
  T4->SetBranchStatus("nBtagTight",1);
  T4->SetBranchStatus("N_Jets",1);
  T4->SetBranchStatus("CSV",1);
};

bool WJetsHTBins::SetMetStatus4(){
  T4->SetBranchStatus("*",0); //disable all branches                                                                                                                                         
  T4->SetBranchStatus("RSQ",1);
  T4->SetBranchStatus("MR",1);
  T4->SetBranchStatus("BOX_NUM",1);
  T4->SetBranchStatus("nBtag",1);
  T4->SetBranchStatus("nBtagTight",1);
  T4->SetBranchStatus("N_Jets",1);
  T4->SetBranchStatus("CSV",1);
  T4->SetBranchStatus("ht",1);
  T4->SetBranchStatus("metX",1);
  T4->SetBranchStatus("metY",1);
  T4->SetBranchStatus("metCorrX",1);
  T4->SetBranchStatus("metCorrY",1);
  T4->SetBranchStatus("N_Jets",1);
  T4->SetBranchStatus("idMC", 1);
  T4->SetBranchStatus("statusMC", 1);
};


double WJetsHTBins::HLTscale(double MR, double R2){

  int MRbin = -1;
  int R2bin = -1;

  const double R2A[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MRA[] = {200., 300., 400. ,500., 2000.};
  
  int Nbins = 4;

  for(int j = 0; j <= Nbins; j++){
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
  
  //return hlt->GetBinContent( MRbin, R2bin );
  return eff->GetEfficiency(eff->GetGlobalBin(MRbin, R2bin , 0));

};

double WJetsHTBins::HLTscaleEle(double MR, double R2){

  int MRbin = -1;
  int R2bin = -1;

  const double R2A[] = {0.35, 0.5, 0.6, 0.8, 1.5};
  const double MRA[] = {200., 300., 400. ,500., 2000.};

  int Nbins = 4;

  for(int j = 0; j <= Nbins; j++){
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

  //return hlt_ele->GetBinContent( MRbin, R2bin );
  return eff_ele->GetEfficiency(eff_ele->GetGlobalBin(MRbin, R2bin , 0));
  
};



/*
double WJetsHTBins::HLTscale(double MR, double R2){
  double MRhigh = 2000.;
  double MRlow = 200.;
  int Nbins = 50;
  std::vector<double> MRlimit;
  double stepMR = (MRhigh - MRlow)/Nbins;

  double R2high = 1.5;
  double R2low = 0.5;
  std::vector<double> R2limit;
  double stepR2 = (R2high - R2low)/Nbins;

  for(int i = 0; i <= Nbins; i++){
    MRlimit.push_back(MRlow+i*stepMR);
    R2limit.push_back(R2low+i*stepR2);
  }

  int binH = Nbins;
  int binL = 0;
  int bin = (int)binH/2;
  int MRbin = -1;
  
  for(int j = 0; j < Nbins; j++){//Binary Search!!
    if( MR > MRlimit[bin]){
      if(MR < MRlimit[bin + 1]){
        MRbin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  MR < MRlimit[bin] ){
      if(MR > MRlimit[bin - 1]){
        MRbin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      MRbin = bin+1;
      break;
    }

  }
  
  binH = Nbins;
  binL = 0;
  bin = (int)binH/2;
  int R2bin = -1;
  for(int j = 0; j < Nbins; j++){//Binary Search!!                                                                
    if( R2 > R2limit[bin]){
      if(R2 < R2limit[bin + 1]){
        R2bin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  R2 < R2limit[bin] ){
      if(R2 > R2limit[bin - 1]){
        R2bin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      R2bin = bin+1;
      break;
    }
    
  }
  
  return hlt->GetBinContent( MRbin, R2bin );
  
};

double WJetsHTBins::HLTscaleEle(double MR, double R2){
  double MRhigh = 2000.;
  double MRlow = 200.;
  int Nbins = 50;
  std::vector<double> MRlimit;
  double stepMR = (MRhigh - MRlow)/Nbins;
  
  double R2high = 1.5;
  double R2low = 0.5;
  std::vector<double> R2limit;
  double stepR2 = (R2high - R2low)/Nbins;
  
  for(int i = 0; i <= Nbins; i++){
    MRlimit.push_back(MRlow+i*stepMR);
    R2limit.push_back(R2low+i*stepR2);
  }
  
  int binH = Nbins;
  int binL = 0;
  int bin = (int)binH/2;
  int MRbin = -1;
  
  for(int j = 0; j < Nbins; j++){//Binary Search!!                                                      
    if( MR > MRlimit[bin]){
      if(MR < MRlimit[bin + 1]){
        MRbin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  MR < MRlimit[bin] ){
      if(MR > MRlimit[bin - 1]){
        MRbin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{
      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      MRbin = bin+1;
      break;
    }

  }

  binH = Nbins;
  binL = 0;
  bin = (int)binH/2;
  int R2bin = -1;
  for(int j = 0; j < Nbins; j++){//Binary Search!!
    if( R2 > R2limit[bin]){
      if(R2 < R2limit[bin + 1]){
        R2bin = bin+1;
        break;
      }else{
        binL = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else if(  R2 < R2limit[bin] ){
      if(R2 > R2limit[bin - 1]){
        R2bin = bin;
        break;
      }else{
        binH = bin;
        bin = binL + (int)(binH-binL)/2;
        continue;
      }
    }else{

      std::cout << "number is:  " << MRlimit[bin] << std::endl;
      R2bin = bin+1;
      break;
    }
    
  }
  
  return hlt_ele->GetBinContent( MRbin, R2bin );
  
};
*/
