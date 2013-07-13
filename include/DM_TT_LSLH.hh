#ifndef DM_TT_LSLH_HH
#define DM_TT_LSLH_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <vector>
#include <math.h>


class TTJets{
  
public:
  
  
  static const int MR_Bins = 6;
  static const int RSQ_Bins = 6;
  
  static const double sigma0 = 136.3*1000.0*0.105;// (LO PREP)
  static const double sigma1 = 136.3*1000.0*0.438;//(LO PREP)
  static const double sigma2 = 136.3*1000.0*0.457;//fb (LO PREP)
  static const float Lumi = 5.;//fb^{-1}
  
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];

  static const int btagIndex = 3;
  
  TTJets();
  TTJets(int );
  TTJets(const char* );
  
  ~TTJets();
  
  TH1F PlotMR_2Box();
  TH1F PlotMR_1Box();
  TH1F PlotMR_0Box();
    
  TH1F PlotRSQ_2Box();
  TH1F PlotRSQ_1Box();
  TH1F PlotRSQ_0Box();
  
  TH2F PlotRSQ_vs_MR_0Box();
  TH2F PlotRSQ_vs_MR_1Box();
  TH2F PlotRSQ_vs_MR_2Box();

  bool pfJetPassCSVM(double );
  int pfJetPassCSVM(double*, int);
  
  std::vector<TH2F*> Plot_2DRazor();
  std::vector<TH1F*> Plot_1DRazor();
  std::vector<TH1F*> PlotMETmag();
  
  bool PrintEvents();
  bool SetStatus();
  bool SetMetStatus();
  bool SetStatus1();
  bool SetMetStatus1();
  bool SetStatus2();
  bool SetMetStatus2();

  double HLTscale(double, double);
  double HLTscaleEle(double, double);
  
  bool SetBtagCut(int a, int b, int c){nBtagCut[0] = a; nBtagCut[1] = b; nBtagCut[2] = c;};

private:
  
  TTree* T;//TTJets_Leptonic
  TTree* T1;//TTJets_SemiLeptonic
  TTree* T2;//TTJets_Hadronic
  
  TFile* F;
  TFile* F1;
  TFile* F2;
  
  TEfficiency* eff;
  TEfficiency* eff_ele;

  int metIndex;
  float MRMin;
  float RSQMin;
    
  float weight0;
  float weight1;
  float weight2;
  
  bool fBtag[5];
  TString BtagBranch[2];

  int nBtagCut[3];
    
};

#endif
