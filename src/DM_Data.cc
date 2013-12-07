#include "DM_Data.hh"
#include <iostream>
#include <iomanip>

Data::Data(){
  
};

Data::Data(const char* FileName, int MetIndex ): BaseDM( FileName, "Data", MetIndex ){
  
  F = TFile::Open(FileName);
  T = (TTree*)F->Get("outTree");
};

Data::~Data(){
  delete T;
};

