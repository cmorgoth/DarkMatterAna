#ifndef DM_DYJETS_LL
#define DM_DYJETS_LL 1

#include "DM_Base.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <stdio.h>
#include <math.h>

class DYJets: public BaseDM{
  
public:
  
  //static const double sigma = 3503.71*1000.;//fb (NNLO)
  static const double sigma = 2950*1000.;//fb (LO)
  
  
  DYJets();
  DYJets(const char*, int );
  
  ~DYJets();
  
private:
  
  TTree* T;
  TFile* F;
  
};


#endif
