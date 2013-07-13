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
  
  CreateStackPlots();
  //CreateRatioPlots();
  //Create2DPlots();
  //CreateMetPlots();
  
  return 0;
  
}  






