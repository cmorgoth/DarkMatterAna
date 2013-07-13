#include "DM_Data.hh"
#include <iostream>
#include <iomanip>

Data::Data(){
  
};

Data::Data(const char* FileName, int MetIndex ): BaseDM( FileName, "Data", MetIndex ){
  
  //F = new TFile(FileName);
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
  //CalcWeight();
};

Data::~Data(){
  delete T;
};

