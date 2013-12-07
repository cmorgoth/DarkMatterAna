#ifndef DM_ZJETSNUNUV3_HH
#define DM_ZJETSNUNUV3_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "DM_BaseV2.hh"
#include "TEfficiency.h"
#include <vector>
#include <math.h>

class ZJetsNuNu: public BaseDM{
  
public:
  
  ZJetsNuNu();
  ZJetsNuNu(int );
  
  ~ZJetsNuNu();
  
};


#endif
