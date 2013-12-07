#ifndef DM_DATAV2_HH
#define DM_DATAV2_HH 1

#include "TH1F.h"
#include "DM_BaseV2.hh"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <math.h>

class Data: public BaseDM{
  
public:
   
  Data();
  Data(const char* FileName, int MetIndex);
  //Data(char const*, int );
  
  ~Data();

  
  TH1F PlotMR_1Box();
  TH1F PlotMR_0Box();
  
  TH1F PlotRSQ_1Box();
  TH1F PlotRSQ_0Box();
  
  TH2F PlotRSQ_vs_MR_0Box();
  TH2F PlotRSQ_vs_MR_1Box();
  
  bool pfJetPassCSVM(double );
  int pfJetPassCSVM(double*, int);
  
  std::vector<TH2F*> Plot_2DRazor();
  std::vector<TH1F*> Plot_1DRazor();
  std::vector<TH1F*> PlotMETmag();
  
  bool PrintEvents();
  bool SetStatus();
  bool SetMetStatus();
    
  //private:
  
  //TChain* T;
    
};


#endif
