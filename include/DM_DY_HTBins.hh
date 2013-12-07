#ifndef DM_DY_HH_HT
#define DM_DY_HH_HT 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <vector>
#include <math.h>
#include "TLorentzVector.h"

class DY{
  
public:

  
  //static const int MR_Bins = 4;
  static const int MR_Bins = 5;
  //static const int RSQ_Bins = 4;
  static const int RSQ_Bins = 5;

  
  static const double sigma0 = 19.73*1000.;// (LO PREP)
  static const double sigma1 = 2.826*1000.;//(LO PREP)
  
  /*
  static const double sigma0 = 1.188*19.73*1000.;// (NNLO)
  static const double sigma1 = 1.188*2.826*1000.;//(NNLO)
  */
  
  //static const float Lumi = 19.6;//fb^{-1}
  static const float Lumi = 19.364;
  
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];
  
  static const int btagIndex = 0;//0->Veto Btag(Loose), 1-> Btag(Loose) >=1, 2-> BtagTight >=1 
  
  DY();
  DY(int );
  DY(const char* );
  
  ~DY();
  
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
  std::vector<TH1F*> DoubleMuBoxPlots();
  
  bool PrintEvents();
  bool SetStatus();
  bool SetMetStatus();
  bool SetStatus1();
  bool SetMetStatus1();

  double HLTscale(double, double);
  double HLTscaleEle(double, double);
  
  bool SetBtagCut(int a, int b, int c){nBtagCut[0]=a; nBtagCut[1]=b; nBtagCut[2]=c;};
  
private:
  
  TTree* T;//HT_200_140
  TTree* T1;//HT_400_inf
    
  TFile* F;
  TFile* F1;
    
  TH2F* hlt;
  TH2F* hlt_ele;
  TEfficiency* eff;
  TEfficiency* eff_ele;

  int metIndex;
  float MRMin;
  float RSQMin;
    
  float weight0;
  float weight1;
  
  bool fBtag[5];
  int nBtagCut[3];
  TString BtagBranch;
    
};

#endif
