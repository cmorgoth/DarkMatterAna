#ifndef DM_WJetsHTBinsV3_HH
#define DM_WJetsHTBinsV3_HH 1

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "DM_BaseV2.hh"
#include "TEfficiency.h"
#include <vector>
#include <math.h>

class WJetsHTBins: public BaseDM{
  
public:
  
  WJetsHTBins();
  WJetsHTBins(int );
  
  ~WJetsHTBins();
  
};


#endif
