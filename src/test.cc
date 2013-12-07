#include <iostream>
#include <sstream>
#include <fstream>
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
/*5x5 v1
const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0, 2.50};
const float ZJetsNuNu::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};

const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0,  2.50};
const float WJetsHTBins::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};

const float DY::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0, 2.50};
const float DY::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0, 2.50};
const float TTJets::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};

const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.0, 2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 800., 3500.};
*/

//5x5 v2
const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};
const float ZJetsNuNu::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.};

const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1,  2.50};
const float WJetsHTBins::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.};

const float DY::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};
const float DY::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};                         
const float TTJets::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.};

const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};                       
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.}; 


int main(){
  
  int bL, bM, bT;
  bL = bM = 1;
  bT = 1;

  std::ofstream ofs("Yields_LO_TwoBtag_NewBinning_New.tex", std::ofstream::out);
  
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
  
  std::vector<TH2F*> TTjets2D = TT->Plot_2DRazor();
  for(int i = 0; i < 3; i++)std::cout << i << " " << (TTjets2D[i])->Integral() << std::endl;
  TH2F* TT_2D_0mu = new TH2F(*TTjets2D[0]);
  TH2F* TT_2D_1mu = new TH2F(*TTjets2D[1]);
  TH2F* TT_2D_2mu = new TH2F(*TTjets2D[2]);
  
  std::cout << "2d tt 0mu: " << TT_2D_0mu->Integral() << std::endl;
  std::cout << "2d tt 1mu: " << TT_2D_1mu->Integral() << std::endl;
  std::cout << "2d tt 2mu: " << TT_2D_2mu->Integral() << std::endl;
  
  std::vector<TH1F*> TTplots = TT->DoubleMuBoxPlots();

  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  /*
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
  */
  
  ///////////////////////////////////
  ////////////Z(nunu)+Jets///////////
  //////////////////////////////////
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  
  std::vector<TH2F*> Zjets2D = Z->Plot_2DRazor();
  for(int i = 0; i < 3; i++)std::cout << i << " " << (Zjets2D[i])->Integral() << std::endl;
  
  TH2F*  MR_RSQ_0BOX_Z = new TH2F( *Zjets2D[0] );
  TH2F*  MR_RSQ_1BOX_Z = new TH2F( *Zjets2D[1] );
  TH2F*  MR_RSQ_2BOX_Z = new TH2F( *Zjets2D[2] );
  //MR_RSQ_0BOX_Z->Sumw2();
  std::cout << "Z(nunu) 0 box: " << MR_RSQ_0BOX_Z->Integral() << std::endl;
  std::cout << "Z(nunu) 1 box: " << MR_RSQ_2BOX_Z->Integral() << std::endl;
  std::cout << "Z(nunu) 2 box: " << MR_RSQ_1BOX_Z->Integral() << std::endl;
  
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
  
  
  std::vector<TH2F*> Wjets2D = W->Plot_2DRazor();
  for(int i = 0; i < 3; i++)std::cout << i << " " << (Wjets2D[i])->Integral() << std::endl;
  
  TH2F*  MR_RSQ_0BOX_W = new TH2F( *Wjets2D[0] );
  TH2F*  MR_RSQ_1BOX_W = new TH2F( *Wjets2D[1] );
  TH2F*  MR_RSQ_2BOX_W = new TH2F( *Wjets2D[2] );

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

  std::vector<TH2F*> dy_jets2D = dy->Plot_2DRazor();
  for(int i = 0; i < 3; i++)std::cout << i << " " << (dy_jets2D[i])->Integral() << std::endl;

  TH2F*  MR_RSQ_0BOX_DY = new TH2F( *dy_jets2D[0] );
  MR_RSQ_0BOX_DY->Sumw2();
  TH2F*  MR_RSQ_1BOX_DY = new TH2F( *dy_jets2D[1] );
  MR_RSQ_1BOX_DY->Sumw2();
  TH2F*  MR_RSQ_2BOX_DY = new TH2F( *dy_jets2D[2] );
  MR_RSQ_2BOX_DY->Sumw2();

  std::cout << "DY0mu: " << MR_RSQ_0BOX_DY->Integral() << std::endl;
  std::cout << "DY1mu: " << MR_RSQ_1BOX_DY->Integral() << std::endl;
  std::cout << "DY2mu: " << MR_RSQ_2BOX_DY->Integral() << std::endl;
  
  //Creating Total Bkg Histograms
  TH1F* R2_bkg_0mu = new TH1F( *RSQ_0 );//Wjets
  R2_bkg_0mu->Add(RSQ_00_TT);//TTJets
  R2_bkg_0mu->Add(RSQ_00_Z);//Z(nunu)Jets
  R2_bkg_0mu->Add(RSQ_dy_00);//Z(ll)Jets
  
  TH1F* R2_bkg_1mu = new TH1F( *RSQ_1 );//Wjets                                                                    
  R2_bkg_1mu->Add(RSQ_11_TT);//TTJets                                                                              
  R2_bkg_1mu->Add(RSQ_11_Z);//Z(nunu)Jets                                                                          
  R2_bkg_1mu->Add(RSQ_dy_11);//Z(ll)Jets
  
  TH1F* R2_bkg_2mu = new TH1F( *RSQ_2 );//Wjets                                                                    
  R2_bkg_2mu->Add(RSQ_22_TT);//TTJets                                                                              
  R2_bkg_2mu->Add(RSQ_22_Z);//Z(nunu)Jets                                                                          
  R2_bkg_2mu->Add(RSQ_dy_22);//Z(ll)Jets

  //MR Total
  TH1F* MR_bkg_0mu = new TH1F( *MR_0 );//Wjets
  MR_bkg_0mu->Add(MR_00_TT);//TTJets
  MR_bkg_0mu->Add(MR_00_Z);//Z(nunu)Jets
  MR_bkg_0mu->Add(MR_dy_00);//Z(ll)Jets 
  
  TH1F* MR_bkg_1mu = new TH1F( *MR_1 );//Wjets                                                                     
  MR_bkg_1mu->Add(MR_11_TT);//TTJets                                                                               
  MR_bkg_1mu->Add(MR_11_Z);//Z(nunu)Jets                                                                           
  MR_bkg_1mu->Add(MR_dy_11);//Z(ll)Jets

  TH1F* MR_bkg_2mu = new TH1F( *MR_2 );//Wjets
  MR_bkg_2mu->Add(MR_22_TT);//TTJets                                                                     
  MR_bkg_2mu->Add(MR_22_Z);//Z(nunu)Jets                                                                           
  MR_bkg_2mu->Add(MR_dy_22);//Z(ll)Jets
  
  //2D Plots

  TH2F* bkg_2d_0mu = new TH2F( *MR_RSQ_0BOX_DY );
  bkg_2d_0mu->Add(MR_RSQ_0BOX_W);
  bkg_2d_0mu->Add(MR_RSQ_0BOX_Z);
  bkg_2d_0mu->Add(TT_2D_0mu);

  TH2F* bkg_2d_1mu = new TH2F( *MR_RSQ_1BOX_DY );
  bkg_2d_1mu->Add(MR_RSQ_1BOX_W);
  bkg_2d_1mu->Add(MR_RSQ_1BOX_Z);
  bkg_2d_1mu->Add(TT_2D_1mu);
  
  TH2F* bkg_2d_2mu = new TH2F( *MR_RSQ_2BOX_DY );
  bkg_2d_2mu->Add(MR_RSQ_2BOX_W);
  bkg_2d_2mu->Add(MR_RSQ_2BOX_Z);
  bkg_2d_2mu->Add(TT_2D_2mu);
  
  
  std::vector<TH1F*> dyplots = dy->DoubleMuBoxPlots();

  //TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  //C1->cd();
  /*
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
  */

  const char* data_file = "/media/data/cmorgoth/Data/DMData/FullHTMHTRereco/HTMHT_ABCD_FullLumi20128TeV.root";

  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  
  std::vector<TH1F*> data_h = data->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (data_h[i])->Integral() << std::endl;
  
  TH1F* MR_22_data = new TH1F( *data_h[4] );
  std::cout << "Data MR 2BOX: " << MR_22_data->Integral() << std::endl;
  MR_22_data->Sumw2();
  
  TH1F* MR_11_data = new TH1F( *data_h[2] );
  std::cout << "Data MR 1BOX: " << MR_11_data->Integral() << std::endl;
  MR_11_data->Sumw2();
  
  TH1F* MR_00_data = new TH1F( *data_h[0] );
  std::cout << "Data MR 0BOX: " << MR_00_data->Integral() << std::endl;
  MR_00_data->Sumw2();

  TH1F* RSQ_22_data = new TH1F( *data_h[5] );
  std::cout << "Data RSQ 2BOX: " << RSQ_22_data->Integral() << std::endl;
  RSQ_22_data->Sumw2();
  
  TH1F* RSQ_11_data = new TH1F( *data_h[3] );
  std::cout << "Data RSQ 1BOX: " << RSQ_11_data->Integral() << std::endl;
  RSQ_11_data->Sumw2();
  
  TH1F* RSQ_00_data = new TH1F( *data_h[1] );
  std::cout << "Data RSQ 0BOX: " << RSQ_00_data->Integral() << std::endl;
  RSQ_00_data->Sumw2();
  
  std::vector<TH2F*> data_h2D = data->Plot_2DRazor();
  for(int i = 0; i < 3; i++)std::cout << i << " " << (data_h2D[i])->Integral() << std::endl;
  
  TH2F* data_2d_0mu = new TH2F( *data_h2D[0] );
  TH2F* data_2d_1mu = new TH2F( *data_h2D[1] );
  TH2F* data_2d_2mu = new TH2F( *data_h2D[2] );

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
  
  double tt_mu[3], dy_mu[3], w_mu[3], z_mu[3];
  double tt_mu_E[3], dy_mu_E[3], w_mu_E[3], z_mu_E[3];
  
  tt_mu[0] = TT_2D_0mu->IntegralAndError(1,4, tt_mu_E[0]);
  tt_mu[1] = TT_2D_1mu->IntegralAndError(1,4, tt_mu_E[1]);
  tt_mu[2] = TT_2D_2mu->IntegralAndError(1,4, tt_mu_E[2]);

  dy_mu[0] = MR_RSQ_0BOX_DY->IntegralAndError(1,4, dy_mu_E[0]);
  dy_mu[1] = MR_RSQ_1BOX_DY->IntegralAndError(1,4, dy_mu_E[1]);
  dy_mu[2] = MR_RSQ_2BOX_DY->IntegralAndError(1,4, dy_mu_E[2]);

  w_mu[0] = MR_RSQ_0BOX_W->IntegralAndError(1,4, w_mu_E[0]);
  w_mu[1] = MR_RSQ_1BOX_W->IntegralAndError(1,4, w_mu_E[1]);
  w_mu[2] = MR_RSQ_2BOX_W->IntegralAndError(1,4, w_mu_E[2]);
  
  z_mu[0] = MR_RSQ_0BOX_Z->IntegralAndError(1,4, z_mu_E[0]);
  z_mu[1] = MR_RSQ_1BOX_Z->IntegralAndError(1,4, z_mu_E[1]);
  z_mu[2] = MR_RSQ_2BOX_Z->IntegralAndError(1,4, z_mu_E[2]);
  
  double tot_0mu = tt_mu[0]+dy_mu[0]+w_mu[0]+z_mu[0];
  double tot_1mu = tt_mu[1]+dy_mu[1]+w_mu[1]+z_mu[1];
  double tot_2mu = tt_mu[2]+dy_mu[2]+w_mu[2]+z_mu[2];

  double tot_0mu_E = sqrt(tt_mu_E[0]*tt_mu_E[0] + dy_mu_E[0]*dy_mu_E[0] + w_mu_E[0]*w_mu_E[0] + z_mu_E[0]*z_mu_E[0]);
  double tot_1mu_E = sqrt(tt_mu_E[1]*tt_mu_E[1] + dy_mu_E[1]*dy_mu_E[1] + w_mu_E[1]*w_mu_E[1] + z_mu_E[1]*z_mu_E[1]);
  double tot_2mu_E = sqrt(tt_mu_E[2]*tt_mu_E[2] + dy_mu_E[2]*dy_mu_E[2] + w_mu_E[2]*w_mu_E[2] + z_mu_E[2]*z_mu_E[2]);
  
  ofs << "\\begin{table}[htdp]\n\\caption{default}\n\\begin{center}\n\\begin{tabular}{|c|c|c|c|}\n\\hline\n";
  
  ofs << "\t&\t$0-\\mu BOX$\t&\t$1-\\mu BOX$\t&\t$2-\\mu BOX$\\\\\n\\hline";
  
  ofs << "$t\\bar{t}$ + Jets\t&\t" << tt_mu[0] << "$\\pm$" << tt_mu_E[0] << "\t&\t" << tt_mu[1] << "$\\pm$" << tt_mu_E[1] << "\t&\t" << tt_mu[2] << "$\\pm$" << tt_mu_E[2] <<"\\\\\n";
  
  ofs << "\\hline\n$Z(ll)$ + Jets\t&\t" << dy_mu[0] << "$\\pm$" << dy_mu_E[0] << "\t&\t" << dy_mu[1] << "$\\pm$" << dy_mu_E[1] << "\t&\t" << dy_mu[2] << "$\\pm$" << dy_mu_E[2] <<"\\\\\n";
  
  ofs << "\\hline\n$Z(\\nu\\bar{\\nu})$ + Jets\t&\t" << z_mu[0] << "$\\pm$" << z_mu_E[0] << "\t&\t" << z_mu[1] << "$\\pm$" << z_mu_E[1] << "\t&\t" << z_mu[2] << "$\\pm$" << z_mu_E[2] << "\\\\\n";
  
  ofs << "\\hline\n$W(l\\nu)$ + Jets\t&\t" <<  w_mu[0] << "$\\pm$" << w_mu_E[0] << "\t&\t" << w_mu[1] << "$\\pm$" << w_mu_E[1] << "\t&\t" << w_mu[2] << "$\\pm$" << w_mu_E[2] <<"\\\\\n";
  
  ofs << "\\hline\n\\hline\nTotal MC\t&\t" << tot_0mu << "$\\pm$" << tot_0mu_E << "\t&\t" << tot_1mu << "$\\pm$" << tot_1mu_E << "\t&\t" << tot_2mu << "$\\pm$" << tot_2mu_E << "\\\\\n";
  ofs << "\\hline\nData\t&\t" << data_2d_0mu->Integral() << "\t&\t" << data_2d_1mu->Integral() << "\t&\t" << data_2d_2mu->Integral() << "\\\\\n\\hline";
  ofs << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
  ofs.close();

  TFile* f1 = new TFile("TwoBtag_LOXsec_NewBinning_New.root","RECREATE");

  //Total Bkg
  R2_bkg_0mu->Write("R2_bkg_0mu");
  R2_bkg_1mu->Write("R2_bkg_1mu");
  R2_bkg_2mu->Write("R2_bkg_2mu");
  
  MR_bkg_0mu->Write("MR_bkg_0mu");
  MR_bkg_1mu->Write("MR_bkg_1mu");
  MR_bkg_2mu->Write("MR_bkg_2mu");
  
  bkg_2d_0mu->Write("bkg_2d_0mu");
  bkg_2d_1mu->Write("bkg_2d_1mu");
  bkg_2d_2mu->Write("bkg_2d_2mu");
  
  //Individual Histos
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
  MR_RSQ_1BOX_W->Write("W_2d_1mu");
  MR_RSQ_0BOX_W->Write("W_2d_0mu");
  
  RSQ_22_data->Write("data_R2_2mu");
  RSQ_11_data->Write("data_R2_1mu");
  RSQ_00_data->Write("data_R2_0mu");
  
  data_2d_0mu->Write("data_2d_0mu");
  data_2d_1mu->Write("data_2d_1mu");
  data_2d_2mu->Write("data_2d_2mu");

  f1->Close();
  
  return 0;
  
}  






