#include "DM_T2CC.hh"
#include <iostream>
#include <iomanip>

T2cc::T2cc(){
  
};
T2cc::T2cc(const char* FileName, int MetIndex ): BaseDM( FileName, "T2cc", sigma, MetIndex ){
  
  //F = new TFile(FileName);
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  //CalcWeight();
  
  //MRMin = 200.0;//nominal 150.
  //RSQMin = 0.5;//nomimal 0.5
};
T2cc::~T2cc(){
  delete T;
};

