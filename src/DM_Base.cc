#include "DM_Base.hh"
#include <iostream>
#include <iomanip>

BaseDM::BaseDM(){
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
  //Does nothing
};


BaseDM::BaseDM( const char* FileName, TString pName, float sig, int MetIndex ){
  
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");

  this->sigma = sig;//Sets cross sections from derived class namely (DYJets, WJets, TTJets, Data)
  this->processName = pName;
  this->metIndex = MetIndex;
  
  ///////////////////////////////////////
  ///////////Trigger Turn-On////////////
  /////////////////////////////////////
  TFile* file = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_DoubleMuonPD_Final.root");
  eff = (TEfficiency*)file->Get("Eff2d");

  TFile* file1 = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_SignleElePD_Final.root");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
  
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
  
  CalcWeight();
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
};

BaseDM::BaseDM( const char* FileName, TString pName , int MetIndex ){
  
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  
  this->processName = pName;
  this->weight = 1.0;
  this->sigma = -99.;
  this->metIndex = MetIndex;
  
  ///////////////////////////////////////   
  ///////////Trigger Turn-On////////////                                                                               
  /////////////////////////////////////                                                                                 
  
  TFile* file = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_DoubleMuonPD_Final.root");
  eff = (TEfficiency*)file->Get("Eff2d");
  
  TFile* file1 = new TFile("/media/data/cmorgoth/TriggerDM/hlt_eff_SignleElePD_Final.root");
  eff_ele = (TEfficiency*)file->Get("Eff2d");
     
  if( btagIndex == 0 || btagIndex == 1 ){
    this->BtagBranch = "nBtag";
  }else{
    this->BtagBranch = "nBtagTight";
  }
  
  std::cout << "----------------Branch Name:  " << BtagBranch << std::endl;
    
  MRMin = 200.0;//nominal 150.
  RSQMin = 0.5;//nomimal 0.5
  
};


BaseDM::~BaseDM(){
  //delete T;
};


/*bool BaseDM::FillVecBadHcalEvents(){
  
  ifstream file ;
  file.open("AllBadHCALLaser/AllBadHCALLaser.txt" );
  std::string buffer;
  
  if( file.is_open() ){
    while( file.good() ){
      file >> buffer;
      if( file.eof() )break;
      HcalBadRuns.insert( std::pair<std::string, bool>(buffer, true) );
      //std::cout << "=== buffer: " << buffer << std::endl;
    }
    
  }else{
    std::cout << "==========unable to open bad Hcal runs=======" << std::endl;
  }
  
  file.close();
  
};
*/

bool BaseDM::PrintEvents(){
  std::setprecision(15);
  std::cout << std::fixed;
  
  double NtotGen = 0, Nt_PV = 0, Nt_2J = 0, Nt_0b = 0, Nt_LepVeto = 0, N_In, N_PV, N_2J, N_0b, N_LepVeto;
  TTree* effT = (TTree*)F->Get("effTree");
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
  std::cout << "============== " << this->processName << " ==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "weighted Nt_In: " << NtotGen*weight << std::endl;
  std::cout << "weighted Nt_PV: " << Nt_PV*weight << std::endl;
  std::cout << "weighted Nt_2J: " << Nt_2J*weight << std::endl;
  std::cout << "weighted Nt_0b: " << Nt_0b*weight << std::endl;
  std::cout << "weighted Nt_LepVeto: " << Nt_LepVeto*weight << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries()*weight << "\n\n" << std::endl; 

  std::cout << "========================================" << std::endl;
  std::cout << "============== " << this->processName << " ==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "MC Nt_In: " << NtotGen << std::endl;
  std::cout << "MC Nt_PV: " << Nt_PV << std::endl;
  std::cout << "MC Nt_2J: " << Nt_2J << std::endl;
  std::cout << "MC Nt_0b: " << Nt_0b << std::endl;
  std::cout << "MC Nt_LepVeto: " << Nt_LepVeto << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries() << "\n\n" << std::endl;

  double RSQ[4], MR[4]/*, run, ls, evNum*/;
  int BOX, nBtag;
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  //T->SetBranchAddress("run", &run);
  //T->SetBranchAddress("ls", &ls);
  //T->SetBranchAddress("evNum", &evNum);
  
  long Nt_MR_RSQ_cut = .0,  Nt_2muBox = .0, Nt_1muBox = .0, Nt_0muBox = .0,	\
    Nt_MR_RSQ_cut0BTag = .0,  Nt_2muBox0BTag = .0, Nt_1muBox0BTag = .0, Nt_0muBox0BTag = .0;
   
  
     
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    /*if( isData ){
      ss << run << ":" << ls << ":" << evNum;
      
      if( HcalBadRuns.find(ss.str()) != HcalBadRuns.end() )std::cout << "======RJE=====" << std::endl;
      
      }*/
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin ){
      Nt_MR_RSQ_cut++;
      if( BOX == 0 )Nt_0muBox++;
      if( BOX == 1 )Nt_1muBox++;
      if( BOX == 2 )Nt_2muBox++;
      
      if( nBtag == 0 ){
	Nt_MR_RSQ_cut0BTag++;
	if( BOX == 0 )Nt_0muBox0BTag++;
	if( BOX == 1 )Nt_1muBox0BTag++;
	if( BOX == 2 )Nt_2muBox0BTag++;
	
      }
      
      
    }
  }

  std::cout << "weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut*weight << std::endl;
  std::cout << "weighted Nt_0muBox: " << Nt_0muBox*weight << std::endl;
  std::cout << "weighted Nt_1muBox: " << Nt_1muBox*weight << std::endl;
  std::cout << "weighted Nt_2muBox: " << Nt_2muBox*weight << std::endl;
  
  std::cout << "weighted Nt_MR_RSQ_cut0BTag: " << Nt_MR_RSQ_cut0BTag*weight << std::endl;
  std::cout << "weighted Nt_0muBox0BTag: " << Nt_0muBox0BTag*weight << std::endl;
  std::cout << "weighted Nt_1muBox0BTag: " << Nt_1muBox0BTag*weight << std::endl;
  std::cout << "weighted Nt_2muBox0BTag: " << Nt_2muBox0BTag*weight << std::endl;
  
  std::cout << this->processName << "Btag EVENTS: " << (Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag)*weight << std::endl;

  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;
  
  return true;
  
};


bool BaseDM::CalcWeight(){
  
  double Ntot = 0, N_In;
  TTree* effT = (TTree*)F->Get("effTree");
  effT->SetBranchAddress("Npassed_In", &N_In);
  for (int i = 0; i < effT->GetEntries(); i++){
    effT->GetEntry(i);
    Ntot += N_In;
  }
  
  weight = Lumi*sigma/Ntot;
  std::cout << "BaseDM weight Process " << this->processName << " is ----> " << weight << std::endl;
  return true;
  
};

TH1F BaseDM::PlotMR_1Box(){
  
  TString BoxCut = "BOX == 1 && ";
  double RSQ[4], MR[4], CSV[30], metX[4], metY[4];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* MR1 = new TH1F("MR1", "MR1BOX", MR_Bins, MR_BinArr);
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX==1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
	MR1->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
	if( hltWeight == 0.0 )hltWeight = 1.0;
	MR1->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *MR1;
};

TH1F BaseDM::PlotMR_0Box(){
  
  double RSQ[4], MR[4], CSV[30], metX[4], metY[4];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* MR0 = new TH1F("MR0", "MR0BOX", MR_Bins, MR_BinArr);
  
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );

    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        MR0->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	MR0->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  return *MR0;
};


TH1F  BaseDM::PlotRSQ_1Box(){
  
  double RSQ[4], MR[4], CSV[30], metX[4], metY[4];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* RSQ1 = new TH1F("RSQ1", "RSQ1BOX", RSQ_Bins, RSQ_BinArr);
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(RSQ[metIndex] > 1.5)std::cout << "out of limits" << std::endl; 
      if(this->processName == "Data"){
        RSQ1->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ1->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ1;
  
};


TH1F  BaseDM::PlotRSQ_0Box(){
  
  double RSQ[4], MR[4], CSV[30], metX[4], metY[4];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH1F* RSQ0 = new TH1F("RSQ0", "RSQ0BOX", RSQ_Bins, RSQ_BinArr);
 
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);

  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        RSQ0->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ0->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ0;
  
};


TH2F BaseDM::PlotRSQ_vs_MR_0Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_0BOX_data", "RSQ_VS_MR_0BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 0 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscaleEle( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
};

TH2F BaseDM::PlotRSQ_vs_MR_1Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_1BOX_data", "RSQ_VS_MR_1BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  std::cout << "Nentries:  " << T->GetEntries() << std::endl;
  for(long j = 0; j < T->GetEntries(); j++){
    T->GetEntry(j);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 1 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
};

TH2F BaseDM::PlotRSQ_vs_MR_2Box(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4];
  double mht[3], CSV[30];
  double hltWeight = 0.0;
  int BOX, N_Jets, nBtag[2];
  
  TH2F* h = new TH2F("RSQ_MR_2BOX_data", "RSQ_VS_MR_2BOX_data", MR_Bins, MR_BinArr, RSQ_Bins, RSQ_BinArr);
  h->Sumw2();
  
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("ht", &ht);
  T->SetBranchAddress("mht", mht);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metCorrX", metcorrX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchAddress("metCorrY", metcorrY);
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    
    if( BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        h->Fill(MR[metIndex], RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	h->Fill(MR[metIndex], RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *h;
  
};


TH1F  BaseDM::PlotRSQ_2Box(){
  
  double RSQ[4], MR[4], CSV[30], run, evNum, ls, metX[4], metY[4];
  double hltWeight = 0.0, Mu_E[2], Mu_Px[2], Mu_Py[2], Mu_Pz[2];
  int BOX, N_Jets, nBtag[2];
  TH1F* RSQ2 = new TH1F("RSQ2_DY", "RSQ2BOX_DY", RSQ_Bins, RSQ_BinArr);
  
  SetBrachStatus();
  //T->SetBranchAddress("run", &run);
  //T->SetBranchAddress("evNum", &evNum);
  //T->SetBranchAddress("ls", &ls);
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    //if((run == 191859 && ls == 96) || (run == 193621 && ls == 628) || (run == 194424 && ls == 537) || (run == 194151 && ls == 551) || (run == 195398 && ls == 58) || (run == 194897 && ls == 50))continue;
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    
    TLorentzVector mu1(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
    TLorentzVector mu2(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
    TLorentzVector sum_mu;
    sum_mu = mu1 + mu2;
    double Mass = sum_mu.M();
    
    if( /*Mass > 50.  && Mass < 180. N_Jets == 2 &&*/ BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        RSQ2->Fill(RSQ[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	RSQ2->Fill(RSQ[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *RSQ2;
  
};


TH1F BaseDM::PlotMR_2Box(){
  
  double RSQ[4], MR[4], CSV[30], metX[4], metY[4];
  double hltWeight = 0.0, Mu_E[4], Mu_Px[4], Mu_Py[4], Mu_Pz[4];
  int BOX, N_Jets, nBtag[2];
  TH1F* MR2 = new TH1F("MR2_DY", "MR2BOX_DY", MR_Bins, MR_BinArr);
  SetBrachStatus();
  T->SetBranchAddress("RSQ", RSQ);
  T->SetBranchAddress("MR", MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
  T->SetBranchAddress("N_Jets", &N_Jets);
  T->SetBranchAddress("CSV", CSV);
  T->SetBranchAddress("metX", metX);
  T->SetBranchAddress("metY", metY);
  T->SetBranchStatus("Mu_E",1);
  T->SetBranchStatus("Mu_Px",1);
  T->SetBranchStatus("Mu_Py",1);
  T->SetBranchStatus("Mu_Pz",1);
  T->SetBranchAddress("Mu_E", Mu_E);
  T->SetBranchAddress("Mu_Px", Mu_Px);
  T->SetBranchAddress("Mu_Py", Mu_Py);
  T->SetBranchAddress("Mu_Pz", Mu_Pz);

  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    double MET = sqrt(metX[2]*metX[2]+metY[2]*metY[2]);
    TLorentzVector mu1(Mu_Px[0], Mu_Py[0], Mu_Pz[0], Mu_E[0]);
    TLorentzVector mu2(Mu_Px[1], Mu_Py[1], Mu_Pz[1], Mu_E[1]);
    TLorentzVector sum_mu;
    sum_mu = mu1 + mu2;
    double Mass = sum_mu.M();

    if( /*Mass > 50. && Mass < 180. && N_Jets == 2 &&*/ BOX == 2 && RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin && fBtag[btagIndex] ){
      if(this->processName == "Data"){
        MR2->Fill(MR[metIndex], weight);
      }else{
	hltWeight = HLTscale( MR[metIndex], RSQ[metIndex]);
        if( hltWeight == 0.0 )hltWeight = 1.0;
	MR2->Fill(MR[metIndex], weight*hltWeight);
      }
    }
  }
  
  return *MR2;
};

std::vector<TH1F*> BaseDM::DoubleMuBoxPlots(){
  
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30], Mu_E[2], Mu_Px[2], Mu_Py[2], Mu_Pz[2];
  int BOX, N_Jets, nBtag[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;

  std::vector<TH1F*> vec_plot;
  TH1F* plot_2mu[6];
  
  plot_2mu[0] = new TH1F( "Mass", "Mass", 15, .0, 500.);
  plot_2mu[1] = new TH1F( "Angle", "Angle", 15, .0, 2*3.1416);
  plot_2mu[2] = new TH1F( "Pt1", "Pt1", 15, .0, 500);
  plot_2mu[3] = new TH1F( "Pt2", "Pt2", 15, .0, 500);
  plot_2mu[4] = new TH1F( "Eta1", "Eta1", 15, -3.0, 3.0);
  plot_2mu[5] = new TH1F( "Eta2", "Eta2", 15, -3.0, 3.0);
  
  SetBrachStatus();
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
  T->SetBranchAddress("nBtag", &nBtag[0]);
  T->SetBranchAddress("nBtagTight", &nBtag[1]);
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
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere                                         
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
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
      
      if(Mass > 0. && Mass < 1800.){
	plot_2mu[0]->Fill(Mass);
	plot_2mu[1]->Fill(angle);
	plot_2mu[2]->Fill(pt1);
	plot_2mu[3]->Fill(pt2);
	plot_2mu[4]->Fill(eta1);
	plot_2mu[5]->Fill(eta2);
      }
      
    }
    
  }
  
  for(int j = 0; j < 6; j++){
    vec_plot.push_back(plot_2mu[j]);
  }
  
  return vec_plot;
  
};

std::vector<TH1F*> BaseDM::PlotMETx(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  
  std::vector< TH1F* > metvec;
  
  TH1F* METx[9];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_metXBox%d_plotType%d",l,m));
      METx[3*l + m] = new TH1F( name, name, 50, -300, 300 );
    }
  }
  
  SetBrachStatus();
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
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	METx[0]->Fill(metX[metIndex], weight);
	METx[1]->Fill(metcorrX[metIndex], weight);
	METx[2]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }else if( BOX == 1 ){
	METx[3]->Fill(metX[metIndex], weight);
	METx[4]->Fill(metcorrX[metIndex], weight);
	METx[5]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }else if( BOX == 2 ){
	METx[6]->Fill(metX[metIndex], weight);
	METx[7]->Fill(metcorrX[metIndex], weight);
	METx[8]->Fill(metcorrX[metIndex]-metX[metIndex], weight);
      }
    }
    
  }
  
  for(int j = 0; j < 9; j++){
    metvec.push_back(METx[j]);
  }
    
  return metvec;
};


std::vector<TH1F*> BaseDM::PlotMETy(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  
  std::vector< TH1F* > metvec;
  TH1F* MET[9];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_metY_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, -300, 300 );
    }
  }
  
  SetBrachStatus();
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
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	MET[0]->Fill(metY[metIndex], weight);
	MET[1]->Fill(metcorrY[metIndex], weight);
	MET[2]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
      }else if( BOX == 1 ){
	MET[3]->Fill(metY[metIndex], weight);
	MET[4]->Fill(metcorrY[metIndex], weight);
	MET[5]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
      }else if( BOX == 2 ){
	MET[6]->Fill(metY[metIndex], weight);
	MET[7]->Fill(metcorrY[metIndex], weight);
	MET[8]->Fill(metcorrY[metIndex]-metY[metIndex], weight);
	//std::cout << "METY debug: " << "entry: "  <<  i << " RSQ: " << RSQ << " MR: " << MR \
	//	  << " nBtag:" << nBtag << std::endl;
      }
    }
    
  }
  
  for(int j = 0; j < 9; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec; 
};


std::vector<TH1F*> BaseDM::PlotMETmag(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4], run/*, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, nBtag[2], N_Jets;
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  
  std::vector< TH1F* > metvec;
  TH1F* MET[12];
  TString name;
  for(int l = 0; l < 3; l++ ){
    for( int m = 0; m < 3; m++ ){
      name = TString(Form("BaseDM_METmag_Box%d_plotType%d",l,m));
      MET[3*l + m] = new TH1F( name, name, 50, 0, 1000);
    }
  }
  
  MET[9] = new TH1F( this->processName + "_NJETS0_Z", "NJETS 0 BOX", 9, 1, 10);
  MET[10] = new TH1F( this->processName + "_NJETS1_Z", "NJETS 1 BOX", 9, 1, 10);
  MET[11] = new TH1F( this->processName + "_NJETS2_Z", "NJETS 2 BOX", 9, 1, 10);
  
  SetBrachStatus();
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
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);

  double metmag =0;
  double metmagcorr = 0;
  double wt = 1.;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);

    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere                                         
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    metmag = sqrt(metX[metIndex]*metX[metIndex]+metY[metIndex]*metY[metIndex]);
    metmagcorr = sqrt(metcorrX[metIndex]*metcorrX[metIndex]+metcorrY[metIndex]*metcorrY[metIndex]);
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0){
	if(this->processName == "Data"){
	  wt = 1.0;
	}else{
	  wt = HLTscaleEle(MR[metIndex], RSQ[metIndex]);
	  if( wt == 0.0)wt = 1.0;
	}
	MET[0]->Fill(metmag, weight*wt);
	MET[1]->Fill(metmagcorr, weight*wt);
	MET[2]->Fill(metmagcorr-metmag, weight*wt);
	MET[9]->Fill(N_Jets, weight*wt);
      }else if( BOX == 1 ){
	if(this->processName == "Data"){
          wt = 1.0;
        }else{
          wt = HLTscale(MR[metIndex], RSQ[metIndex]);
          if( wt == 0.0)wt = 1.0;
        }
	MET[3]->Fill(metmag, weight*wt);
	MET[4]->Fill(metmagcorr, weight*wt);
	MET[5]->Fill(metmagcorr-metmag, weight*wt);
	MET[10]->Fill(N_Jets,weight*wt);
      }else if( BOX == 2 ){
	if(this->processName == "Data"){
          wt = 1.0;
        }else{
          wt = HLTscale(MR[metIndex], RSQ[metIndex]);
          if( wt == 0.0)wt = 1.0;
        }
	MET[6]->Fill(metmag, weight*wt);
	MET[7]->Fill(metmagcorr, weight*wt);
	MET[8]->Fill(metmagcorr-metmag, weight*wt);
	MET[11]->Fill(N_Jets,weight*wt);
      }
    }
    
  }
  
  for(int j = 0; j < 12; j++){
    metvec.push_back(MET[j]);
  }
  
  return metvec; 
};

std::vector<TH1F*> BaseDM::PlotHT(){
  double metX[4], metcorrX[4], metY[4], metcorrY[4], ht, RSQ[4], MR[4]/*, run, ls, evNum*/;
  double mht[3], CSV[30];
  int BOX, N_Jets, nBtag[2];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  
  std::vector< TH1F* > metvec;
  
  TH1F* HT[3];
  TString name;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("BaseDM_HT_Box%d",l));
    HT[l] = new TH1F( name, name, 50, 135, 600 );
  }
  
  SetBrachStatus();
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
  T->SetBranchAddress("pTHem1", &pTHem1);
  T->SetBranchAddress("pTHem2", &pTHem2);
  T->SetBranchAddress("etaHem1", &etaHem1);
  T->SetBranchAddress("etaHem2", &etaHem2);
  T->SetBranchAddress("phiHem1", &phiHem1);
  T->SetBranchAddress("phiHem2", &phiHem2);
  
  //std::cout << "Tree entries: " << T->GetEntries() << std::endl;
    for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    TLorentzVector j1;
    TLorentzVector j2;
    j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
    j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere                                         
    double Dphi = j1.DeltaPhi(j2);
    
    fBtag[0] = (nBtag[0] == 0);
    fBtag[1] = fBtag[2] = (nBtag[0] >= 1);
    fBtag[3] = ( nBtag[1] >= nBtagCut[2] && nBtag[0] >= nBtagCut[0] );
    int nBtagMed = pfJetPassCSVM(CSV, N_Jets);
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1] );
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if(  RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5 ){
      if( BOX == 0 ){
	HT[0]->Fill(ht, weight);
      }else if( BOX == 1 ){
	HT[1]->Fill(ht, weight);
      }else if( BOX == 2 ){
	HT[2]->Fill(ht, weight);
      }
      
    }
    
  }
  
  for(int j = 0; j < 3; j++){
    metvec.push_back(HT[j]);
  }
  
  return metvec;
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

std::vector<TH2F*> BaseDM::Plot_2DRazor(){
  std::cout << "Weight: " << weight << std::endl;
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

  SetBrachStatus();
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
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);
    
    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
	Razor2D[0]->Fill(MR[metIndex], RSQ[metIndex]);
      }else if( BOX == 1 ){
	Razor2D[1]->Fill(MR[metIndex], RSQ[metIndex]);
      }else if( BOX == 2 ){
	Razor2D[2]->Fill(MR[metIndex], RSQ[metIndex]);
      }
    }
  }
  T->SetBranchStatus("*", 0);
  
  for(int j = 0; j < 3; j++){
    Razor2DVec.push_back(Razor2D[j]);
  }
  
  return Razor2DVec;

};

std::vector<TH1F*> BaseDM::Plot_1DRazor(){
  double RSQ[4], MR[4], CSV[30];
  double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
  int BOX, N_Jets, nBtag[2];

  std::vector< TH1F* > Razor1DVec;
  TH1F* Razor1D[6];
  TString name, name1;
  double hltWeight;
  for(int l = 0; l < 3; l++ ){
    name = TString(Form("MR_1D_TT_%dmu_Box",l));
    name1 = TString(Form("R2_1D_TT_%dmu_Box",l));
    Razor1D[2*l] = new TH1F( name, name, MR_Bins, MR_BinArr);
    Razor1D[2*l+1] = new TH1F( name1, name1, RSQ_Bins, RSQ_BinArr);
    Razor1D[2*l]->Sumw2();
    Razor1D[2*l+1]->Sumw2();
  }

  SetBrachStatus();
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
    //fBtag[4] = ( nBtag[1] >= nBtagCut[2] && nBtagMed >= nBtagCut[1]);
    fBtag[4] = (nBtag[1] >= nBtagCut[2]);

    if( RSQ[metIndex] > RSQMin && MR[metIndex] > MRMin  && fBtag[btagIndex] && fabs(Dphi) < 2.5){
      if( BOX == 0){
	Razor1D[0]->Fill(MR[metIndex]);
        Razor1D[1]->Fill(RSQ[metIndex]);
      }else if( BOX == 1 ){
	Razor1D[2]->Fill(MR[metIndex]);
        Razor1D[3]->Fill(RSQ[metIndex]);
      }else if( BOX == 2 ){
	Razor1D[4]->Fill(MR[metIndex]);
        Razor1D[5]->Fill(RSQ[metIndex]);
      }
    }
  }

  for(int j = 0; j < 6; j++){
    Razor1DVec.push_back(Razor1D[j]);
  }
  
  return Razor1DVec;
  
}

bool BaseDM::SetBrachStatus(){
  T->SetBranchStatus("*",0); //disable all branches
  T->SetBranchStatus("run",1);
  T->SetBranchStatus("evNum",1);
  T->SetBranchStatus("ls",1);
  T->SetBranchStatus("RSQ",1);
  T->SetBranchStatus("MR",1);
  T->SetBranchStatus("BOX_NUM",1);
  T->SetBranchStatus("nBtag",1);
  T->SetBranchStatus("nBtagTight",1);
  T->SetBranchStatus("N_Jets",1);
  T->SetBranchStatus("CSV",1);
  T->SetBranchStatus("ht",1);
  T->SetBranchStatus("mht",1);
  T->SetBranchStatus("metX",1);
  T->SetBranchStatus("metY",1);
  T->SetBranchStatus("metCorrX",1);
  T->SetBranchStatus("metCorrY",1);
  T->SetBranchStatus("pTHem1", 1);
  T->SetBranchStatus("pTHem2", 1);
  T->SetBranchStatus("etaHem1", 1);
  T->SetBranchStatus("etaHem2", 1);
  T->SetBranchStatus("phiHem1", 1);
  T->SetBranchStatus("phiHem2", 1);
  T->SetBranchStatus("N_Jets",1);
  //T->SetBranchStatus("run",1);
  //T->SetBranchStatus("ls",1);
  //T->SetBranchStatus("evNum",1);
  
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

