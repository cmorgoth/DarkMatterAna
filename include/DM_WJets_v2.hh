#ifndef DM_WJETS_HH
#define DM_WJETS_HH 1

#include "DM_Base.hh"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <stdio.h>
#include <math.h>
//using namespace std;

class WJets : public BaseDM{
  
public:
  
  //static const double sigma = 36257.2*1000.;//fb (NNLO)
  static const double sigma = 30400*1000.;//fb (LO)  
  
  WJets();
  WJets(const char*, int  );
  
  ~WJets();
 
private:
  
  TTree* T;
  TFile* F;

  //TString pname;
  
};



#endif
