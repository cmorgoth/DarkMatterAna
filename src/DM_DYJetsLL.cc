#include "DM_DYJetsLL.hh"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

DYJets::DYJets(){
  
};

DYJets::DYJets(const char* FileName, int MetIndex ): BaseDM( FileName, "DY+Jets", sigma, MetIndex ){
  
  //F = new TFile(FileName);
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
    
};

DYJets::~DYJets(){
  delete T;
};

/*
bool DYJets::PrintEvents(){
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
  std::cout << "================DY JETS==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "DYLL weighted Nt_In: " << NtotGen*weight << std::endl;
  std::cout << "DYLL weighted Nt_PV: " << Nt_PV*weight << std::endl;
  std::cout << "DYLL weighted Nt_2J: " << Nt_2J*weight << std::endl;
  std::cout << "DYLL weighted Nt_0b: " << Nt_0b*weight << std::endl;
  std::cout << "DYLL weighted Nt_LepVeto: " << Nt_LepVeto*weight << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries()*weight << "\n\n" << std::endl; 
  //std::cout << "========================================" << std::endl;

  double RSQ, MR;
  int BOX, nBtag;
  SetBrachStatus();
  T->SetBranchAddress("RSQ", &RSQ);
  T->SetBranchAddress("MR", &MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  
  long Nt_MR_RSQ_cut = 0, Nt_2muBox = 0,  Nt_1muBox = 0, Nt_0muBox = 0,	\
    Nt_MR_RSQ_cut0BTag = 0, Nt_2muBox0BTag = 0,  Nt_1muBox0BTag = 0, Nt_0muBox0BTag = 0;
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    if( RSQ > RSQMin && MR > MRMin ){
      Nt_MR_RSQ_cut++;
      if( BOX == 0 )Nt_0muBox++;
      if( BOX == 1 )Nt_1muBox++;
      if( BOX == 2 )Nt_2muBox++;
      
      if( nBtag == 0){
	Nt_MR_RSQ_cut0BTag++;
	if( BOX == 0 )Nt_0muBox0BTag++;
	if( BOX == 1 )Nt_1muBox0BTag++;
	if( BOX == 2 )Nt_2muBox0BTag++;
      }
      
    }
  }
  
  //std::cout << "========================================" << std::endl;

  std::cout << "DYLL weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut*weight << std::endl;
  std::cout << "DYLL weighted Nt_0muBox: " << Nt_0muBox*weight << std::endl;
  std::cout << "DYLL  weighted Nt_1muBox: " << Nt_1muBox*weight << std::endl;
  std::cout << "DYLL  weighted Nt_2muBox: " << Nt_2muBox*weight << std::endl;

  std::cout << "DYLL weighted Nt_MR_RSQ_cut0BTag: " << Nt_MR_RSQ_cut0BTag*weight << std::endl;
  std::cout << "DYLL weighted Nt_0muBox0BTag: " << Nt_0muBox0BTag*weight << std::endl;
  std::cout << "DYLL  weighted Nt_1muBox0BTag: " << Nt_1muBox0BTag*weight << std::endl;
  std::cout << "DYLL  weighted Nt_2muBox0BTag: " << Nt_2muBox0BTag*weight << std::endl;
  
  std::cout << "DYLL Btag EVENTS: " << (Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag)*weight << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;

  return true;
  
};
*/


/*
  bool DYJets::PrintEventsMC(){
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
  std::cout << "================DY JETS==================" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "DYLL weighted Nt_In: " << NtotGen*weight << std::endl;
  std::cout << "DYLL weighted Nt_PV: " << Nt_PV*weight << std::endl;
  std::cout << "DYLL weighted Nt_2J: " << Nt_2J*weight << std::endl;
  std::cout << "DYLL weighted Nt_0b: " << Nt_0b*weight << std::endl;
  std::cout << "DYLL weighted Nt_LepVeto: " << Nt_LepVeto*weight << std::endl;

  std::cout << "NEntries weighted: " << T->GetEntries()*weight << "\n\n" << std::endl; 
  //std::cout << "========================================" << std::endl;
  
  double RSQ, MR;
  int BOX, nBtag;
  
  int arrMax = 150;
  int nMC, idMC[arrMax], statusMC[arrMax], mothMC[arrMax];
  float energyMC[arrMax], etaMC[arrMax];
  double Mu_E[2], Mu_Py[2], Mu_Px[2];
  
  
  T->SetBranchAddress("RSQ", &RSQ);
  T->SetBranchAddress("MR", &MR);
  T->SetBranchAddress("BOX_NUM", &BOX);
  T->SetBranchAddress("nBtag", &nBtag);
  T->SetBranchAddress("nMC", &nMC);
  T->SetBranchAddress("idMC", &idMC);
  T->SetBranchAddress("mothMC", &mothMC);
  T->SetBranchAddress("statusMC", &statusMC);
  T->SetBranchAddress("energyMC", &energyMC);
  T->SetBranchAddress("etaMC", &etaMC);
  T->SetBranchAddress("Mu_E", &Mu_E);
  T->SetBranchAddress("Mu_Py", &Mu_Py);
  T->SetBranchAddress("Mu_Px", &Mu_Px);
   
  
  long Nt_MR_RSQ_cut = 0, Nt_2muBox = 0,  Nt_1muBox = 0, Nt_0muBox = 0,	\
    Nt_MR_RSQ_cut0BTag = 0, Nt_2muBox0BTag = 0,  Nt_1muBox0BTag = 0, Nt_0muBox0BTag = 0;

  float diff1, diff2;
  int jmin1 = 0, jmin2 = 0, ctr = 0, flag = 0;
  
  for(int i = 0; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    diff1 = diff2 = 999999.;
    flag = 0;
    
    for( int j = 0; j < nMC; j++ ){
      if( abs(idMC[j]) == 13 && energyMC[j]/cosh( etaMC[j] ) > 25&& statusMC[j] == 1 ){
	if( fabs( Mu_E[0] - energyMC[j] ) < diff1 ){
	  diff1 = fabs( Mu_E[0] - energyMC[j] );
	  jmin1 = j;
	}
	if( fabs( Mu_E[1] - energyMC[j] ) < diff2 ){
	  diff2 = fabs( Mu_E[1] - energyMC[j] );
	  jmin2 = j;
	}
	if(flag == 0){
	  ctr++;
	  flag = 1;
	}
	//std::cout << "=====================PASSSSSSSSSS=========================  " << energyMC[j]/cosh( etaMC[j] ) << std::endl;
	//std::cout << "=================================================" << std::endl;
	//std::cout << "j: "  << j << "  energyMC[j]: " << energyMC[j] << \
	  //" etaMC[j]: " << etaMC[j] << " PtMC[j]: " << energyMC[j]/cosh( etaMC[j] ) << "  mothMC[j]:  " << mothMC[j] << std::endl;
	  
	  
	  //std::cout << "\nMu_E[0]: " << Mu_E[0] << "  Mu_Pt[0]: " <<  sqrt( Mu_Py[0]*Mu_Py[0] + Mu_Px[0]*Mu_Px[0] )<< std::endl;
	  //std::cout << "\nMu_E[1]: " << Mu_E[1] << "  Mu_Pt[1]: " <<  sqrt( Mu_Py[1]*Mu_Py[1] + Mu_Px[1]*Mu_Px[1] )<< std::endl;
	  //std::cout << "=================================================" << std::endl;
	  }
      //std::cout << "nMC: " << nMC << "  idMC[6]:  " << idMC[6] << "  mothMC[6]:  " << mothMC[6] << "  statusMC[6]: " << statusMC[6] << std::endl;
      //std::cout << "nMC: " << nMC << "  idMC[7]:  " << idMC[7] << "  mothMC[7]:  " << mothMC[7] << "  statusMC[7]: " << statusMC[7] << std::endl;
      //std::cout << "nMC: " << nMC << "  idMC[8]:  " << idMC[8] << "  mothMC[8]:  " << mothMC[8] << "  statusMC[8]: " << statusMC[8] << std::endl;
    }
    
    //std::cout << "====================CTRCTRCTR===========================  " << ctr << std::endl;
    //std::cout << "j: "  << jmin1 << "  energyMC[j]: " << energyMC[jmin1] << \
      ///  " etaMC[j]: " << etaMC[jmin1] << " PtMC[j]: " << energyMC[jmin1]/cosh( etaMC[jmin1] ) << "  mothMC[j]:  " << mothMC[jmin1] << std::endl;
    
    
      // std::cout << "\nMu_E[0]: " << Mu_E[0] << "  Mu_Pt[0]: " <<  sqrt( Mu_Py[0]*Mu_Py[0] + Mu_Px[0]*Mu_Px[0] )<< std::endl;

      //std::cout << "j: "  << jmin2 << "  energyMC[j]: " << energyMC[jmin2] <<	\
      //      " etaMC[j]: " << etaMC[jmin2] << " PtMC[j]: " << energyMC[jmin2]/cosh( etaMC[jmin2] ) << "  mothMC[j]:  " << mothMC[jmin2] << std::endl;
      // std::cout << "\nMu_E[1]: " << Mu_E[1] << "  Mu_Pt[1]: " <<  sqrt( Mu_Py[1]*Mu_Py[1] + Mu_Px[1]*Mu_Px[1] )<< std::endl;
    
    if( RSQ > RSQMin && MR > MRMin ){
      Nt_MR_RSQ_cut++;
      if( BOX == 0 )Nt_0muBox++;
      if( BOX == 1 )Nt_1muBox++;
      if( BOX == 2 )Nt_2muBox++;
      
      if( nBtag == 0){
	Nt_MR_RSQ_cut0BTag++;
	if( BOX == 0 )Nt_0muBox0BTag++;
	if( BOX == 1 )Nt_1muBox0BTag++;
	if( BOX == 2 )Nt_2muBox0BTag++;
      }
      
    }
  }
  std::cout << "====================CTRCTRCTR===========================  " << ctr*weight << std::endl;
  //std::cout << "========================================" << std::endl;

  std::cout << "DYLL weighted Nt_MR_RSQ_cut: " << Nt_MR_RSQ_cut*weight << std::endl;
  std::cout << "DYLL weighted Nt_0muBox: " << Nt_0muBox*weight << std::endl;
  std::cout << "DYLL  weighted Nt_1muBox: " << Nt_1muBox*weight << std::endl;
  std::cout << "DYLL  weighted Nt_2muBox: " << Nt_2muBox*weight << std::endl;

  std::cout << "DYLL weighted Nt_MR_RSQ_cut0BTag: " << Nt_MR_RSQ_cut0BTag*weight << std::endl;
  std::cout << "DYLL weighted Nt_0muBox0BTag: " << Nt_0muBox0BTag*weight << std::endl;
  std::cout << "DYLL  weighted Nt_1muBox0BTag: " << Nt_1muBox0BTag*weight << std::endl;
  std::cout << "DYLL  weighted Nt_2muBox0BTag: " << Nt_2muBox0BTag*weight << std::endl;
  
  std::cout << "DYLL Btag EVENTS: " << (Nt_MR_RSQ_cut - Nt_MR_RSQ_cut0BTag)*weight << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "========================================" << "\n\n" << std::endl;

  return true;
  
};
*/

