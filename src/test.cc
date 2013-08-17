#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DM_WJetsHTBins.hh"
#include "DM_ZJetsNuNu.hh"
#include "DM_DY_HTBins.hh"
#include "DM_TT_LSLH.hh"
#include "DM_METPlots.hh"
#include "DM_Data.hh"
#include "THStack.h"
#include "TString.h"
#include "DM_StackPlots.hh"
#include "DM_RatioPlots.hh"
#include "DM_2D_MR_RSQ_Dist.hh"
#include "THStack.h"
#include "DM_1DRatio.hh"

using namespace std;
/*
const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.50};
const float ZJetsNuNu::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 3500.};

const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.50};
const float WJetsHTBins::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 3500.};

const float DY::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.50};
const float DY::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 3500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.50};
const float TTJets::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 3500.};

const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 3500.};
*/

const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.7, 0.9, 1.1, 2.50};
const float ZJetsNuNu::MR_BinArr[] = {200., 466., 732., 1000., 3500.};

const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.7, 0.9, 1.1,  2.50};
const float WJetsHTBins::MR_BinArr[] = {200., 466., 732., 1000., 3500.};

const float DY::RSQ_BinArr[] = {0.5, 0.7, 0.9, 1.1, 2.50};
const float DY::MR_BinArr[] = {200., 466., 732., 1000., 3500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.7, 0.9, 1.1, 2.50};
const float TTJets::MR_BinArr[] = {200., 466., 732., 1000., 3500.};

const float BaseDM::RSQ_BinArr[] = {0.5, 0.7, 0.9, 1.1, 2.50};
const float BaseDM::MR_BinArr[] = {200., 466., 732., 1000., 3500.};

int main(){
  
  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;
  ///////////////////////////////
  ////////// tt + jets//////////
  //////////////////////////////
  
  TTJets* TT = new TTJets(2);
  TT->SetBtagCut(bL,bM,bT);
  
  std::vector<TH1F*> TTjets = TT->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (TTjets[i])->Integral() << std::endl;
  
  TH1F* MR_22_TT = new TH1F( *TTjets[4] );
  std::cout << "TTJets MR 2BOX: " << MR_22_TT->Integral() << std::endl;
  
  TH1F* MR_11_TT = new TH1F( *TTjets[2] );
  std::cout << "TTJets MR 1BOX: " << MR_11_TT->Integral() << std::endl;
  
  TH1F* MR_00_TT = new TH1F( *TTjets[0] );
  std::cout << "TTJets MR 0BOX: " << MR_00_TT->Integral() << std::endl;
  
  TH1F* RSQ_22_TT = new TH1F( *TTjets[5] );
  std::cout << "TTJets RSQ 2BOX: " << RSQ_22_TT->Integral() << std::endl;
  
  TH1F* RSQ_11_TT = new TH1F( *TTjets[3] );
  std::cout << "TTJets RSQ 1BOX: " << RSQ_11_TT->Integral() << std::endl;
  
  TH1F* RSQ_00_TT = new TH1F( *TTjets[1] );
  std::cout << "TTJets RSQ 0BOX: " << RSQ_00_TT->Integral() << std::endl;
  
  TH2F* TT_2D_0mu = new TH2F(TT->PlotRSQ_vs_MR_0Box());
  TH2F* TT_2D_1mu = new TH2F(TT->PlotRSQ_vs_MR_1Box());
  TH2F* TT_2D_2mu = new TH2F(TT->PlotRSQ_vs_MR_2Box());
  
  std::cout << "2d tt 0mu: " << TT_2D_0mu->Integral() << std::endl;
  std::cout << "2d tt 1mu: " << TT_2D_1mu->Integral() << std::endl;
  std::cout << "2d tt 2mu: " << TT_2D_2mu->Integral() << std::endl;
  
  std::vector<TH1F*> TTplots = TT->DoubleMuBoxPlots();

  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  TTplots[0]->Draw();
  C1->SaveAs("MassTT.pdf");
  C1->SaveAs("MassTT.png");
  TTplots[1]->Draw();
  C1->SaveAs("delta_thetaTT.pdf");
  C1->SaveAs("delta_thetaTT.png");
  TTplots[2]->Draw();
  C1->SaveAs("pt1TT.pdf");
  C1->SaveAs("pt1TT.png");
  TTplots[3]->Draw();
  C1->SaveAs("pt2TT.pdf");
  C1->SaveAs("pt2TT.png");
  TTplots[4]->Draw();
  C1->SaveAs("eta1TT.pdf");
  C1->SaveAs("eta1TT.png");
  TTplots[5]->Draw();
  C1->SaveAs("eta2TT.pdf");
  C1->SaveAs("eta2TT.png");
  
  ///////////////////////////////////
  ////////////Z(nunu)+Jets///////////
  //////////////////////////////////
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  TH2F*  MR_RSQ_0BOX_Z = new TH2F( Z->PlotRSQ_vs_MR_0Box() );
  MR_RSQ_0BOX_Z->Sumw2();
  std::cout << "Z(nunu) 0 box: " << MR_RSQ_0BOX_Z->Integral() << std::endl;
  
  std::vector<TH1F*> Zjets = Z->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (Zjets[i])->Integral() << std::endl;
  
  TH1F* MR_22_Z = new TH1F( *Zjets[4] );
  std::cout << "ZJets MR 2BOX: " << MR_22_Z->Integral() << std::endl;
  
  TH1F* MR_11_Z = new TH1F( *Zjets[2] );
  std::cout << "ZJets MR 1BOX: " << MR_11_Z->Integral() << std::endl;
  
  TH1F* MR_00_Z = new TH1F( *Zjets[0] );
  std::cout << "ZJets MR 0BOX: " << MR_00_Z->Integral() << std::endl;
  
  TH1F* RSQ_22_Z = new TH1F( *Zjets[5] );
  std::cout << "ZJets RSQ 2BOX: " << RSQ_22_Z->Integral() << std::endl;

  TH1F* RSQ_11_Z = new TH1F( *Zjets[3] );
  std::cout << "ZJets RSQ 1BOX: " << RSQ_11_Z->Integral() << std::endl;
  
  TH1F* RSQ_00_Z = new TH1F( *Zjets[1] );
  std::cout << "ZJets RSQ 0BOX: " << RSQ_00_Z->Integral() << std::endl;
  
  //////////////////////////
  /////////W+Jets//////////
  /////////////////////////
  WJetsHTBins* W = new WJetsHTBins( 2 );
  W->SetBtagCut(bL,bM,bT);
  
  std::vector<TH1F*> Wjets = W->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (Wjets[i])->Integral() << std::endl;
  
  TH1F* RSQ_0 = new TH1F( *Wjets[1] );
  std::cout << "WJets R2 0BOX: " << RSQ_0->Integral() << std::endl;
  
  TH1F* RSQ_1 = new TH1F( *Wjets[3] );
  std::cout << "WJets R2 1BOX: " << RSQ_1->Integral() << std::endl;
  
  TH1F* RSQ_2 = new TH1F( *Wjets[5] );
  std::cout << "WJets R2 2BOX: " << RSQ_2->Integral() << std::endl;
  
  TH1F* MR_0 = new TH1F( *Wjets[0] );
  std::cout << "WJets MR 0BOX: " << MR_0->Integral() << std::endl;

  TH1F* MR_1 = new TH1F( *Wjets[2] );
  std::cout << "WJets MR 1BOX: " << MR_1->Integral() << std::endl;
  
  TH1F* MR_2 = new TH1F( *Wjets[4] );
  std::cout << "WJets MR 2BOX: " << MR_2->Integral() << std::endl;
  
  TH2F*  MR_R2_1BOX = new TH2F( W->PlotRSQ_vs_MR_1Box() );
  TH2F*  MR_RSQ_0BOX = new TH2F( W->PlotRSQ_vs_MR_0Box() );
  

  /////////////////////////
  //////////Drell-Yan//////
  /////////////////////////
  DY* dy = new DY( 2 );
  dy->SetBtagCut(bL,bM,bT);

  
  std::vector<TH1F*> dy_jets = dy->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (dy_jets[i])->Integral() << std::endl;
  
  TH1F* MR_dy_22 = new TH1F( *dy_jets[4] );
  std::cout << "dy Jets MR 2BOX: " << MR_dy_22->Integral() << std::endl;
  
  TH1F* MR_dy_11 = new TH1F( *dy_jets[2] );
  std::cout << "dy Jets MR 1BOX: " << MR_dy_11->Integral() << std::endl;

  TH1F* MR_dy_00 = new TH1F( *dy_jets[0] );
  std::cout << "dy Jets MR 0BOX: " << MR_dy_00->Integral() << std::endl;

  TH1F* RSQ_dy_22 = new TH1F( *dy_jets[5] );
  std::cout << "dy Jets RSQ 2BOX: " << RSQ_dy_22->Integral() << std::endl;

  TH1F* RSQ_dy_11 = new TH1F( *dy_jets[3] );
  std::cout << "dy Jets RSQ 1BOX: " << RSQ_dy_11->Integral() << std::endl;

  TH1F* RSQ_dy_00 = new TH1F( *dy_jets[1] );
  std::cout << "dy Jets RSQ 0BOX: " << RSQ_dy_00->Integral() << std::endl;

  TH2F*  MR_RSQ_0BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_0Box() );
  MR_RSQ_0BOX_DY->Sumw2();
  TH2F*  MR_RSQ_1BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_1Box() );
  MR_RSQ_1BOX_DY->Sumw2();
  TH2F*  MR_RSQ_2BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_2Box() );
  MR_RSQ_2BOX_DY->Sumw2();

  std::cout << "DY0mu: " << MR_RSQ_0BOX_DY->Integral() << std::endl;
  std::cout << "DY1mu: " << MR_RSQ_1BOX_DY->Integral() << std::endl;
  std::cout << "DY2mu: " << MR_RSQ_2BOX_DY->Integral() << std::endl;
  

  std::vector<TH1F*> dyplots = dy->DoubleMuBoxPlots();

  //TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  //C1->cd();
  dyplots[0]->Draw();
  C1->SaveAs("MassDY.pdf");
  C1->SaveAs("MassDY.png");
  dyplots[1]->Draw();
  C1->SaveAs("delta_thetaDY.pdf");
  C1->SaveAs("delta_thetaDY.png");
  dyplots[2]->Draw();
  C1->SaveAs("pt1DY.pdf");
  C1->SaveAs("pt1DY.png");
  dyplots[3]->Draw();
  C1->SaveAs("pt2DY.pdf");
  C1->SaveAs("pt2DY.png");
  dyplots[4]->Draw();
  C1->SaveAs("eta1DY.pdf");
  C1->SaveAs("eta1DY.png");
  dyplots[5]->Draw();
  C1->SaveAs("eta2DY.pdf");
  C1->SaveAs("eta2DY.png");
  

  const char* data_file = "/afs/cern.ch/work/c/cpena/DarkMatter/CMSSW_5_2_3/src/VecbosApp/53X/HTMHT_Run2012A_ILV/out/HTMHT_ILV_Run2012AB.root";

  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  
  
  TH1F* MR_22_data = new TH1F( data->PlotMR_2Box() );
  std::cout << "Data MR 2BOX: " << MR_22_data->Integral() << std::endl;
  MR_22_data->Sumw2();
  
  TH1F* MR_11_data = new TH1F( data->PlotMR_1Box() );
  std::cout << "Data MR 1BOX: " << MR_11_data->Integral() << std::endl;
  MR_11_data->Sumw2();
  
  TH1F* MR_00_data = new TH1F( data->PlotMR_0Box() );
  std::cout << "Data MR 0BOX: " << MR_00_data->Integral() << std::endl;
  MR_00_data->Sumw2();

  TH1F* RSQ_22_data = new TH1F( data->PlotRSQ_2Box() );
  std::cout << "Data RSQ 2BOX: " << RSQ_22_data->Integral() << std::endl;
  RSQ_22_data->Sumw2();
  
  TH1F* RSQ_11_data = new TH1F( data->PlotRSQ_1Box() );
  std::cout << "Data RSQ 1BOX: " << RSQ_11_data->Integral() << std::endl;
  RSQ_11_data->Sumw2();
  
  TH1F* RSQ_00_data = new TH1F( data->PlotRSQ_0Box() );
  std::cout << "Data RSQ 0BOX: " << RSQ_00_data->Integral() << std::endl;
  RSQ_00_data->Sumw2();
  
  TH2F* data_2d_0mu = new TH2F(data->PlotRSQ_vs_MR_0Box());
  TH2F* data_2d_1mu = new TH2F(data->PlotRSQ_vs_MR_1Box());
  TH2F* data_2d_2mu = new TH2F(data->PlotRSQ_vs_MR_2Box());

  std::cout << "data 0mu: " << data_2d_0mu->Integral() << std::endl;
  std::cout << "data 1mu: " << data_2d_1mu->Integral() << std::endl;
  std::cout << "data 2mu: " << data_2d_2mu->Integral() << std::endl;
  
  
  std::vector<TH1F*> dplots = data->DoubleMuBoxPlots();
  dplots[0]->Draw();
  C1->SaveAs("Mass.pdf");
  C1->SaveAs("Mass.png");
  dplots[1]->Draw();
  C1->SaveAs("delta_theta.pdf");
  C1->SaveAs("delta_theta.png");
  dplots[2]->Draw();
  C1->SaveAs("pt1.pdf");
  C1->SaveAs("pt1.png");
  dplots[3]->Draw();
  C1->SaveAs("pt2.pdf");
  C1->SaveAs("pt2.png");
  dplots[4]->Draw();
  C1->SaveAs("eta1.pdf");
  C1->SaveAs("eta1.png");
  dplots[5]->Draw();
  C1->SaveAs("eta2.pdf");
  C1->SaveAs("eta2.png");

  THStack* stack;// = new THStack("stack", "");
  TLegend* leg;
  TString name;
  TH1F* aux2;
  for(int i = 0; i < 6; i++){
    TTplots[i]->SetFillColor(kPink+9);
    TTplots[i]->SetStats(0);
    TTplots[i]->SetTitle("");
    dyplots[i]->SetFillColor(kViolet+9);
    dyplots[i]->SetStats(0);
    dyplots[i]->SetTitle("");
    
    dplots[i]->Sumw2();
    dplots[i]->SetMarkerStyle(20);
    dplots[i]->SetLineColor(1);
    dplots[i]->SetMarkerSize(1.5);
    dplots[i]->SetStats(0);
    dplots[i]->SetTitle("");
    
    leg = new TLegend(0.7,0.7,0.9,0.92);
    
    leg->AddEntry(dyplots[i],"Z(ll) + jets","fl");
    leg->AddEntry(TTplots[i],"t #bar{t} + jets","fl");
    leg->AddEntry(dplots[i],"data","lep");
    
    stack = new THStack("stack", "");
    stack->Add(dyplots[i]);
    stack->Add(TTplots[i]);
    
    name = TString(Form("StackPlots/PlotN_%d.pdf",i));;
    aux2 = new TH1F(*dyplots[i]);
    aux2->Sumw2();
    aux2->Add(TTplots[i],1);
    TString type;
    
    if(i == 0){
      type = "Mass";
    }else if(i == 1){
      type = "Angle";
    }else if(i == 2 || i == 3){
      type = "PT";
    }else if(i == 4 || i == 5){
      type = "eta";
    }
    std::cout << "type: " << type << std::endl; 
    std::cout << "debug 1: " << std::endl;
    RatioPlotsV2(stack, aux2, dplots[i], "MC 2 #mu BOX", "Data 2 #mu BOX", name, type, leg);
    std::cout << "debug 2: " << std::endl;
    
    delete aux2;
    delete stack;
    std::cout << "debug 3: " << std::endl;
  }
  
  
  
  TFile* f1 = new TFile("Btag_1Tight_1Medium.root","RECREATE");

  RSQ_dy_00->Write("dy_R2_0mu");
  RSQ_dy_11->Write("dy_R2_1mu");
  RSQ_dy_22->Write("dy_R2_2mu");
  MR_dy_00->Write("dy_MR_0mu");
  MR_dy_11->Write("dy_MR_1mu");
  MR_dy_22->Write("dy_MR_2mu");

  RSQ_0->Write("W_R2_0mu");
  RSQ_1->Write("W_R2_1mu");
  RSQ_2->Write("W_R2_2mu");
  MR_0->Write("W_MR_0mu");
  MR_1->Write("W_MR_1mu");
  MR_2->Write("W_MR_2mu");

  MR_00_TT->Write("TT_MR_0mu");
  MR_11_TT->Write("TT_MR_1mu");
  MR_22_TT->Write("TT_MR_2mu");
  RSQ_00_TT->Write("TT_R2_0mu");
  RSQ_11_TT->Write("TT_R2_1mu");
  RSQ_22_TT->Write("TT_R2_2mu");
  TT_2D_0mu->Write("TT_2d_0mu");
  TT_2D_1mu->Write("TT_2d_1mu");
  TT_2D_2mu->Write("TT_2d_2mu");

  MR_22_data->Write("data_MR_2mu");
  MR_11_data->Write("data_MR_1mu");
  MR_00_data->Write("data_MR_0mu");

  MR_RSQ_0BOX_Z->Write("Z_2d_0mu");
  MR_00_Z->Write("Z_MR_0mu");
  MR_11_Z->Write("Z_MR_1mu");
  MR_22_Z->Write("Z_MR_2mu");
  RSQ_00_Z->Write("Z_R2_0mu");
  RSQ_11_Z->Write("Z_R2_1mu");
  RSQ_22_Z->Write("Z_R2_2mu");

  MR_RSQ_0BOX_DY->Write("dy_2d_0mu");
  MR_RSQ_1BOX_DY->Write("dy_2d_1mu");
  MR_RSQ_2BOX_DY->Write("dy_2d_2mu");
  MR_R2_1BOX->Write("W_2d_1mu");
  MR_RSQ_0BOX->Write("W_2d_0mu");
  
  RSQ_22_data->Write("data_R2_2mu");
  RSQ_11_data->Write("data_R2_1mu");
  RSQ_00_data->Write("data_R2_0mu");
  
  data_2d_0mu->Write("data_2d_0mu");
  data_2d_1mu->Write("data_2d_1mu");
  data_2d_2mu->Write("data_2d_2mu");

  f1->Close();
  
  return 0;
  
}  






