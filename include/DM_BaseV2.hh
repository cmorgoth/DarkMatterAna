#ifndef DM_BaseV2_HH
#define DM_BaseV2_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <vector>
#include <math.h>


class BaseDM{
  
public:
  
  TChain* T;
  TChain* effT;

  static const int MR_Bins = 5;
  static const int RSQ_Bins = 5;
  
  //static const float Lumi = 19.6;//fb^{-1}
  static const float Lumi = 19.364;
  
  static const int btagIndex = 0;//0->Veto Btag(Loose), 1-> Btag(Loose) >=1, 2-> BtagTight >=1
  
  static const float RSQ_BinArr[RSQ_Bins+1];
  static const float MR_BinArr[MR_Bins+1];
  
  
  BaseDM();
  BaseDM(int , TString);
  BaseDM(const char* );
  
  ~BaseDM();
  
  virtual TH1F PlotMR_1Box();
  virtual TH1F PlotMR_0Box();
  
  virtual TH1F PlotRSQ_1Box();
  virtual TH1F PlotRSQ_0Box();
  
  virtual TH2F PlotRSQ_vs_MR_0Box();
  virtual TH2F PlotRSQ_vs_MR_1Box();

  virtual bool pfJetPassCSVM(double );
  virtual int pfJetPassCSVM(double*, int);

  virtual std::vector<TH2F*> Plot_2DRazor();
  virtual std::vector<TH1F*> Plot_1DRazor();
  virtual std::vector<TH1F*> PlotMETmag();

  virtual bool PrintEvents();
  virtual bool SetStatus();
  virtual bool SetMetStatus();
  
  virtual bool SetBtagCut(int a, int b, int c){nBtagCut[0]=a; nBtagCut[1]=b; nBtagCut[2]=c;};
  virtual double HLTscale(double, double);
  virtual double HLTscaleEle(double, double);

private:
  
  //TChain* T;
  //TChain* effT;
  
  TH2F* hlt;
  TH2F* hlt_ele;
  TEfficiency* eff;
  TEfficiency* eff_ele;
  
  //int metIndex;
  //float MRMin;
  //float RSQMin;
  //TString pName;
    
protected:

  int metIndex;
  float MRMin;
  float RSQMin;
  TString pName;
  
  bool fBtag[5];
  TString BtagBranch;
  int nBtagCut[3];
  
};

#endif
