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


//const float WJets::RSQ_BinArr[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.05, 1.2, 1.5};
//const float WJets::MR_BinArr[] = {150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 1100., 1500.};
//const double BaseDM::RSQ_BinArr[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, \
//				   0.85, 0.90, 0.95, 1.0, 1.05, 1.2, 1.5};
//const double BaseDM::MR_BinArr[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 1100., 1500.};

//const float ZJetsNuNu::RSQ_BinArr[] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.80, \
//				       0.85, 0.90, 0.95, 1.0, 1.05, 1.2, 1.5};
//const float ZJetsNuNu::MR_BinArr[] = {200., 250., 300., 350., 400., 450., 500., 550., 600., 700., 800., 1100., 1500.};
					   
//const double BaseDM::RSQ_BinArr[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};

const double BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const double BaseDM::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 2500.};

//const float ZJetsNuNu::RSQ_BinArr[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float ZJetsNuNu::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float ZJetsNuNu::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 2500.};

//const float WJetsHTBins::RSQ_BinArr[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float WJetsHTBins::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float WJetsHTBins::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 2500.};

const float DY::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float DY::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 2500.};

const float TTJets::RSQ_BinArr[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5};
const float TTJets::MR_BinArr[] = {200., 300., 400., 500., 600., 900., 2500.};

int main(){
  
  CreateStackPlots();
  CreateRatioPlots();
  Create2DPlots();
  CreateMetPlots();
  
  return 0;
  
}  






