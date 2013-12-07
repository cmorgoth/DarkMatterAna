#include <iostream>
#include <sstream>
#include <fstream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"
#include "TString.h"
#include "DM_DataV2.hh"

using namespace std;

//5x5 v2
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.50};         
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600., 900., 3500.}; 

int main(){
  
  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;

  std::ofstream ofs("Yields_LO_TwoBtag_NewBinning_New.tex", std::ofstream::out);
 
  //const char* data_file = "/media/data/cmorgoth/Data/DMData/FullHTMHTRereco/HTMHT_ABCD_FullLumi20128TeV.root";
  const char* data_file = "/media/data/cmorgoth/Data/DMData/ZJets/BtagCorrMC/ZJetsToNuNu_400_HT_inf_pu_mu_LooseBtag.root";

  Data* data = new Data(data_file, 2);
  data->PrintEvents();
  std::cout << "debu-1" << std::endl;
  data->SetBtagCut(bL,bM,bT);
  
  std::cout << "debu0" << std::endl;
  std::vector<TH1F*> data_h = data->Plot_1DRazor();
  std::cout << "debu1" << std::endl;
  for(int i = 0; i < 6; i++)std::cout << i << " " << (data_h[i])->Integral() << std::endl;
  std::cout << "debu2" << std::endl;
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
  
 
  return 0;
  
}  






