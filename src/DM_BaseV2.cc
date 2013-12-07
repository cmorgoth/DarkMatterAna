#include "DM_BaseV2.hh"
#include <iostream>
#include "TLorentzVector.h"

BaseDM::BaseDM(){ };

BaseDM::BaseDM(int metIndex, TString pName){
  
  MRMin = 200.0;//nominal 150.0
  RSQMin = 0.5;//nominal 0.5
  this->metIndex = metIndex;
  this->pName = pName;
  
  if(pName != "Data"){
    TFile* file = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_DoubleMuonPD_Final.root");
    eff = (TEfficiency*)file->Get("Eff2d");
    
    TFile* file1 = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_SignleElePD_Final.root");
    eff_ele = (TEfficiency*)file->Get("Eff2d");
  }else{
    std::cout << "======== IS DATA!==========" << std::endl;
  }
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


BaseDM::BaseDM(const char* FileName ){
  
  //F = new TFile(FileName);
  //T0 = (TTree*)F->Get("outTree");
  
};

BaseDM::~BaseDM(){
  delete T;
};


bool BaseDM::PrintEvents(){

  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;

  effT->SetBranchAddress("Npassed_In", &N_In);
  effT->SetBranchAddress("Npassed_PV", &N_PV);
  effT->SetBranchAddress("Npassed_2Jet", &N_2J);
  effT->SetBranchAddress("Npassed_0btag", &N_0b);
  effT->SetBranchAddress("Npassed_LepVeto", &N_LepVeto);
  
  
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    NtotGen += N_In;
    Nt_PV += N_PV;
    Nt_2J += N_2J;
    Nt_0b += N_0b;
    Nt_LepVeto += N_LepVeto;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "================ W + Jets ==============" << std::endl;
  std::cout << "========================================" << std::endl;
  
  std::cout << "W+jets weighted Nt_In: " << NtotGen << std::endl;
  std::cout << "W+jets weighted Nt_PV: " << Nt_PV << std::endl;
  std::cout << "W+jets weighted Nt_2J: " << Nt_2J << std::endl;
  std::cout << "W+jets weighted Nt_0b: " << Nt_0b << std::endl;
  std::cout << "W+jets weighted Nt_LepVeto: " << Nt_LepVeto << std::endl;
  
  //After Selection cuts
  double RSQ[4], MR[4];
  int BOX, nBtag;
  
  double Nt_MR_RSQ_cut0BTag = 0.0, Nt_2muBox0BTag = 0.0,  Nt_1muBox0BTag = 0.0, Nt_0muBox0BTag = 0.0,\
    Nt_MR_RSQ_cut = 0.0, Nt_2muBox = 0.0,  Nt_1muBox = 0.0, Nt_0muBox = 0.0, Nt_Btags = 0.0, sf_w;

  SetStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  T->SetBranchAddress("sf_w", &sf_w);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut += sf_w;
      
      if( BOX == 0 )Nt_0muBox += sf_w;
      if( BOX == 1 )Nt_1muBox += sf_w;
      if( BOX == 2 )Nt_2muBox += sf_w;
      
      if( nBtag == 0 ){
	
	Nt_MR_RSQ_cut0BTag += sf_w;
	
	if( BOX == 0 )Nt_0muBox0BTag += sf_w;
	if( BOX == 1 )Nt_1muBox0BTag += sf_w;
	if( BOX == 2 )Nt_2muBox0BTag += sf_w;
	
      }
      
    }
    
  }
  T->SetBranchStatus("*", 0);//Disable all branches, trying to gain some performance

  
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

TH1F BaseDM::PlotMR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
    
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
      MR1->Fill(MR[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR1;
};

TH1F BaseDM::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = .0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
  
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
      MR0->Fill(MR[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *MR0;
};


TH1F  BaseDM::PlotRSQ_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
  
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
      RSQ1->Fill(RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
    
  return *RSQ1;
    
};


TH1F  BaseDM::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
  
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
      RSQ0->Fill(RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ0;
  
};

TH2F BaseDM::PlotRSQ_vs_MR_0Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
  
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
      RSQ_MR_0BOX->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);

  return *RSQ_MR_0BOX;
  
};

TH2F BaseDM::PlotRSQ_vs_MR_1Box(){
  
  double RSQ[4], MR[4], CSV[30];
  double hltWeight = 0.0, sf_w;
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
  T->SetBranchAddress("sf_w", &sf_w);
  
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
      RSQ_MR_1BOX->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight);
    }
  }
  T->SetBranchStatus("*", 0);
  
  return *RSQ_MR_1BOX;
  
};


std::vector<TH1F*> BaseDM::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30], pu_w, mu_w, sf_w;
  int BOX, nBtag[2], N_Jets;
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
 
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
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
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
  
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("sf_w", &sf_w);
  
  float metmag = .0;
  float metmagcorr = .0;
  double wt = 1.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= nBtagCut[0]);
    fBtag[3] = (nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
    double Dphi = j1.DeltaPhi(j2);

    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex]  && fabs(Dphi) < 2.5){
      if( BOX == 0 ){
	wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[0]->Fill(metmag, sf_w*wt*pu_w*mu_w);
	MET[1]->Fill(metmagcorr, sf_w*wt*pu_w*mu_w);
	MET[2]->Fill(metmagcorr-metmag, sf_w*wt*pu_w*mu_w);
	MET[9]->Fill(N_Jets,sf_w*wt*pu_w*mu_w);
      }else if( BOX == 1 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[3]->Fill(metmag, sf_w*wt*pu_w*mu_w);
	MET[4]->Fill(metmagcorr, sf_w*wt*pu_w*mu_w);
	MET[5]->Fill(metmagcorr-metmag, sf_w*wt*pu_w*mu_w);
	MET[10]->Fill(N_Jets, sf_w*wt*pu_w*mu_w);
      }else if( BOX == 2 ){
	wt = HLTscale(MR[metIndex], RSQ[metIndex]);
	if( wt == 0.0)wt = 1.0;
	MET[6]->Fill(metmag, sf_w*wt*pu_w*mu_w);
	MET[7]->Fill(metmagcorr, sf_w*wt*pu_w*mu_w);
	MET[8]->Fill(metmagcorr-metmag, sf_w*wt*pu_w*mu_w);
	MET[11]->Fill(N_Jets, sf_w*wt*pu_w*mu_w);
      }
    }
  }
  T->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }

  return metvec;
};

std::vector<TH2F*> BaseDM::Plot_2DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w;
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
  T->SetBranchAddress("pu_w", &pu_w);
  T->SetBranchAddress("mu_w", &mu_w);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  //T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagLCorr", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("sf_w", &sf_w);
  
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
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }
    }
  }
  T->SetBranchStatus("*", 0);

  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;

};

bool BaseDM::pfJetPassCSVM(double btagOutput){
  if(btagOutput < 0.679)   return false;
  return true;
};

int BaseDM::pfJetPassCSVM(double* CSVM, int N_Jets){
  int nMBtag = 0;
  for(int i = 0; i < N_Jets; i++)if(CSVM[i] >= 0.679)nMBtag++;
  return nMBtag;
};

std::vector<TH1F*> BaseDM::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[10];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2, pu_w, mu_w, sf_w;
  int BOX, N_Jets, nBtag[2];

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[12];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 6; l++ ){
    name = TString(Form("MR_1D_Wj_%dmu_Box",l));
    name1 = TString(Form("R2_1D_Wj_%dmu_Box",l));
    Razor1D[2*l] = new TH1F( name, name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F( name1, name1, RSQ_Bins, RSQ_BinArr);
    if( l < 3 ){
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
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  T->SetBranchAddress("sf_w", &sf_w);
  
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
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
        //hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
	hltWeight = eff_ele->GetEfficiency(eff_ele->FindFixBin(MR[2], RSQ[2]));
        //if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[0]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[1]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[6]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[7]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }else if( BOX == 1 ){
        //hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	hltWeight = eff->GetEfficiency(eff->FindFixBin(MR[2], RSQ[2]));
        //if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[2]->Fill(MR[metIndex],  sf_w*hltWeight*pu_w*mu_w);
	Razor1D[3]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[8]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[9]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }else if( BOX == 2 ){
        hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        //if( hltWeight == 0.0 )hltWeight = 1.0;
        Razor1D[4]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[5]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[10]->Fill(MR[metIndex], sf_w*hltWeight*pu_w*mu_w);
	Razor1D[11]->Fill(RSQ[metIndex], sf_w*hltWeight*pu_w*mu_w);
      }
    }
  }  
  T->SetBranchStatus("*", 0);

  for(int j = 0; j < 12; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }

  return Razor1DVec;
  
};


bool BaseDM::SetStatus(){

  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);
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
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  T->SetBranchStatus("sf_w", 1);
};

bool BaseDM::SetMetStatus(){
  T->SetBranchStatus("*",0); //disable all branches                                                                                                     
  T->SetBranchStatus("pu_w",1);
  T->SetBranchStatus("mu_w",1);           
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
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("idMC", 1);
  T->SetBranchStatus("statusMC", 1);
  T->SetBranchStatus("nBtagLCorr", 1);
  T->SetBranchStatus("nBtagLCorrUp", 1);
  T->SetBranchStatus("nBtagLCorrDown", 1);
  T->SetBranchStatus("sf_w", 1);
};

double BaseDM::HLTscale(double MR, double R2){

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

double BaseDM::HLTscaleEle(double MR, double R2){

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


