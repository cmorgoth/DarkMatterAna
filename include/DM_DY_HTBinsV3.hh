#ifndef DM_DYV3_HH
#define DM_DYV3_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "DM_BaseV2.hh"
#include "TEfficiency.h"
#include <vector>
#include <math.h>

class DY: public BaseDM{
  
public:
  
  DY();
  DY(int );
  
  ~DY();
  
};


#endif
