#ifndef DM_DATA
#define DM_DATA 1

#include "TH1F.h"
#include "DM_Base.hh"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include <math.h>

class Data: public BaseDM{
  
public:
   
  Data();
  Data(const char*, int );
  
  ~Data();
  
private:
  
  TTree* T;
  TFile* F;
  
};


#endif
