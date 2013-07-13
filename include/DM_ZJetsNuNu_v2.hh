#ifndef DM_ZJETSNUNU_HH
#define DM_ZJETSNUNU_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <vector>
#include <math.h>


class ZJetsNuNu{
  
public:

  
  static const int MR_Bins = 6;
  static const int RSQ_Bins = 6;
  
  
  static const double sigma0 = 381.2*1000.;// (LO PREP)
  static const double sigma1 = 160.3*1000.;//(LO PREP)
  static const double sigma2 = 41.49*1000.;//fb (LO PREP)
  static const double sigma3 = 5.274*1000.;//fb (LO PREP)
  static const float Lumi = 5.;//fb^{-1}
    
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];
  
  static const int btagIndex = 0;//0->Veto Btag(Loose), 1-> Btag(Loose) >=1, 2-> BtagTight >=1 
  
  ZJetsNuNu();
  ZJetsNuNu(int );
  ZJetsNuNu(const char* );
  
  ~ZJetsNuNu();
  
  TH1F PlotMR_1Box();
  TH1F PlotMR_0Box();
  
  TH1F PlotRSQ_1Box();
  TH1F PlotRSQ_0Box();
  
  TH2F PlotRSQ_vs_MR_0Box();
  TH2F PlotRSQ_vs_MR_1Box();

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
  bool SetStatus3();
  bool SetMetStatus3();

  double HLTscale(double, double);
  double HLTscaleEle(double, double);
  
private:
  
  TTree* T;//HT_50_100
  TTree* T1;//HT_100_200
  TTree* T2;//HT_200_400
  TTree* T3;//HT_400_inf
  
  TFile* F;
  TFile* F1;
  TFile* F2;
  TFile* F3;
  
  TH2F* hlt;
  TH2F* hlt_ele;
  TEfficiency* eff;
  TEfficiency* eff_ele;

  int metIndex;
  float MRMin;
  float RSQMin;

    
  float weight0;
  float weight1;
  float weight2;
  float weight3;

  bool fBtag[3];
  TString BtagBranch;
    
};

#endif
