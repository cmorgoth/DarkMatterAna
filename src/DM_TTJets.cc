#include "DM_TTJets.hh"
#include <iostream>
#include <iomanip>

TTJets::TTJets(){
  
};
TTJets::TTJets(const char* FileName, int MetIndex ): BaseDM( FileName, "TTJets", sigma, MetIndex ){
  
  //F = new TFile(FileName);
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  //CalcWeight();
  
  //MRMin = 200.0;//nominal 150.
  //RSQMin = 0.5;//nomimal 0.5
};
TTJets::~TTJets(){
  delete T;
};

