#include <iostream>
#include <sstream>
#include <fstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "DM_WJetsHTBinsV3.hh"
#include "DM_ZJetsNuNuV3.hh"
#include "DM_DY_HTBinsV3.hh"
#include "DM_TT_LSLHV3.hh"
#include "THStack.h"
#include "TString.h"

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
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};        
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.}; 


int main(){
  
  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;

  std::ofstream ofs("Yields_LO_TwoBtag_NewBinning_New.tex", std::ofstream::out);
  
 
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
  
  ///////////////////////////////////
  ///////////DY/////////////////////
  /////////////////////////////////
  
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
  
  //std::vector<TH1F*> TTplots = TT->DoubleMuBoxPlots();
  
  
  return 0;
  
}  






