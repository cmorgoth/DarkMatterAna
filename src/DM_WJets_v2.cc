#include "DM_WJets_v2.hh"
#include <iostream>
#include <iomanip>

WJets::WJets(){
  //Do nothing
};

WJets::WJets(const char* FileName, int MetIndex ): BaseDM( FileName, "WJets", sigma, MetIndex ){
  
  //F = new TFile(FileName);
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  //CalcWeight();
};

WJets::~WJets(){
  delete T;
};


