#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DM_WJetsHTBins.hh"
#include "DM_ZJetsNuNu_v2.hh"
#include "DM_DY_HTBins.hh"
#include "DM_TT_LSLH.hh"
#include "DM_METPlots.hh"
#include "DM_Data.hh"
#include "THStack.h"
#include "TString.h"
#include "DM_StackPlots.hh"
#include "DM_RatioPlots.hh"
#include "DM_2D_MR_RSQ_Dist.hh"

using namespace std;

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


int main(){
  

  ///////////////////////////////
  ////////// tt + jets//////////
  //////////////////////////////
  TTJets* TT = new TTJets(2);
  TT->SetnBtagCut(2,1);
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

  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  TH2F*  MR_RSQ_0BOX_Z = new TH2F( Z->PlotRSQ_vs_MR_0Box() );
  MR_RSQ_0BOX_Z->Sumw2();
  std::cout << "Z(nunu) 0 box: " << MR_RSQ_0BOX_Z->Integral() << std::endl;


  //////////////////////////
  /////////W+Jets//////////
  /////////////////////////
  WJetsHTBins* W = new WJetsHTBins( 2 );
  W->SetBtagCut(2,1);
  
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
  dy->SetBtagCut(2,2,1);//>= 2 loose, >= 2 med, >= 1 Tight
  
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


  const char* data_file = "/afs/cern.ch/work/c/cpena/DarkMatter/CMSSW_5_2_3/src/VecbosApp/53X/HTMHT_Run2012A_ILV/out/HTMHT_ILV_Run2012AB.root";

  Data* data = new Data(data_file, 2);
  data->SetBtagCut(2,2,1);
  //data->PrintEvents();                                                                           

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
  
  TFile* f1 = new TFile("VetoBtag.root","RECREATE");

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






