#ifndef DM_TTJETS_HH
#define DM_TTJETS_HH 1

#include "DM_Base.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"

//using namespace std;

class TTJets: public BaseDM{
  
public:
  
  //static const double sigma = 234.0*1000.;//fb (NNLO)
  //static const double sigma = 211.0*1000.;//fb (NLO PREP)
  static const double sigma = 136.3*1000.;//fb (LO) from Monojet People 
  
  
  TTJets();
  TTJets(const char*, int );
  
  ~TTJets();
   
private:
  
  TTree* T;
  TFile* F;

  //float MRMin;
  //float RSQMin;
  //float weight;
  
};



#endif
