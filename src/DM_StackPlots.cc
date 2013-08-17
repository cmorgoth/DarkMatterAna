#include "DM_StackPlots.hh"
#include "DM_1DRatio.hh"

void CreateStackPlots(){
 
  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;
  WJetsHTBins* W = new WJetsHTBins( 2 );
  ///////////////////////////////////////////
  ///////////////WJETS 1D HISTOS/////////////
  ///////////////////////////////////////////
  //W->SetBtagIndex(0);
  W->SetBtagCut(bL,bM,bT);
  W->PrintEvents();
  
  std::vector<TH1F*> Wjets = W->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (Wjets[i])->Integral() << std::endl;

  
  TH1F* RSQ_0 = new TH1F( *Wjets[1] );                     
  std::cout << "WJets R2 0BOX: " << RSQ_0->Integral() << std::endl;                                                     
                                                                                                                 
  TH1F* RSQ_1 = new TH1F( *Wjets[3] );                                                                           
  std::cout << "WJets R2 1BOX: " << RSQ_1->Integral() << std::endl;                                                      
                                                                                                                         
  TH1F* MR_0 = new TH1F( *Wjets[0] );                                                                              
  std::cout << "WJets MR 0BOX: " << MR_0->Integral() << std::endl;                                                        
  TH1F* MR_1 = new TH1F( *Wjets[2] );                                                                              
  std::cout << "WJets MR 1BOX: " << MR_1->Integral() << std::endl; 

  
  TH1F* W_MR0 = new TH1F( W->PlotMR_0Box() );
  std::cout << "debug -1" << std::endl;
  TH1F* W_R20 =new TH1F( W->PlotRSQ_0Box() );
  TH1F* W_MR1 =new TH1F( W->PlotMR_1Box() );
  TH1F* W_R21 =new TH1F( W->PlotRSQ_1Box() );

  TLegend* leg;
  std::cout << "debug 0" << std::endl;
  /////////////////////////////////////////////
  ///////////////ZJetsNuNu ANA/////////////////
  /////////////////////////////////////////////
  
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  Z->PrintEvents();

  std::vector<TH1F*> Zjets = Z->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (Zjets[i])->Integral() << std::endl;
  
  TH1F* MR_11 = new TH1F( *Zjets[2] );
  std::cout << "ZJets MR 1BOX: " << MR_11->Integral() << std::endl;

  TH1F* MR_00 = new TH1F( *Zjets[0] );
  std::cout << "ZJets MR 0BOX: " << MR_00->Integral() << std::endl;

  TH1F* RSQ_11 = new TH1F( *Zjets[3] );
  std::cout << "ZJets RSQ 1BOX: " << RSQ_11->Integral() << std::endl;

  TH1F* RSQ_00 = new TH1F( *Zjets[1] );
  std::cout << "ZJets RSQ 0BOX: " << RSQ_00->Integral() << std::endl;
  
  TH1F* Z_MR0 =new TH1F( Z->PlotMR_0Box() );
  TH1F* Z_R20 =new TH1F( Z->PlotRSQ_0Box() );
  TH1F* Z_MR1 =new TH1F( Z->PlotMR_1Box() );
  TH1F* Z_R21 =new TH1F( Z->PlotRSQ_1Box() );

  //////////////////////////////////////////////////////
  ///////////////// Drell-Yan///////////////////////////
  /////////////////////////////////////////////////////

  DY* dy = new DY( 2 );
  dy->SetBtagCut(bL,bM,bT);
  dy->PrintEvents();
  
  std::vector<TH1F*> dy_jets = dy->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (dy_jets[i])->Integral() << std::endl;
  
  TH1F* MR_dy_22 = new TH1F( *dy_jets[4] );
  std::cout << "dy Jets MR 2BOX: " << MR_dy_22->Integral() << std::endl;

  TH1F* MR_dy_11 = new TH1F( *dy_jets[2] );
  std::cout << "dy Jets MR 1BOX: " << MR_dy_11->Integral() << std::endl;

  TH1F* MR_dy_00 = new TH1F( *dy_jets[0] );
  std::cout << "dy Jets MR 0BOX: " << MR_dy_00->Integral() << std::endl;

  TH1F* RSQ_dy_22 = new TH1F( *dy_jets[5] );
  std::cout << "dy Jets RSQ 2BOX: " << RSQ_dy_22->Integral() << std::endl;

  TH1F* RSQ_dy_11 = new TH1F( *dy_jets[3] );
  std::cout << "dy Jets RSQ 1BOX: " << RSQ_dy_11->Integral() << std::endl;

  TH1F* RSQ_dy_00 = new TH1F( *dy_jets[1] );
  std::cout << "dy Jets RSQ 0BOX: " << RSQ_dy_00->Integral() << std::endl;

  TH1F* dy_MR0 =new TH1F( dy->PlotMR_0Box() );
  TH1F* dy_R20 =new TH1F( dy->PlotRSQ_0Box() );
  TH1F* dy_MR1 =new TH1F( dy->PlotMR_1Box() );
  TH1F* dy_R21 =new TH1F( dy->PlotRSQ_1Box() );
  TH1F* dy_MR2 =new TH1F( dy->PlotMR_2Box() );
  TH1F* dy_R22 =new TH1F( dy->PlotRSQ_2Box() );

  //////////////////////////////////////////////////
  /////////////////////////////////////////////////
  /////////////////TTbar + Jets////////////////////
  ////////////////////////////////////////////////
  ///////////////////////////////////////////////
  
  
  TTJets* TT = new TTJets(2);
  TT->SetBtagCut(bL,bM,bT);
  TT->PrintEvents();

  std::vector<TH1F*> TTjets = TT->Plot_1DRazor();
  for(int i = 0; i < 6; i++)std::cout << i << " " << (TTjets[i])->Integral() << std::endl;  
  
  TH1F* MR_22_TT = new TH1F( *TTjets[4] );
  std::cout << "TTJets MR 2BOX: " << MR_22_TT->Integral() << std::endl;

  TH1F* MR_11_TT = new TH1F( *TTjets[2] );
  std::cout << "TTJets MR 1BOX: " << MR_11_TT->Integral() << std::endl;
  
  TH1F* MR_00_TT = new TH1F( *TTjets[0] );
  std::cout << "TTJets MR 0BOX: " << MR_00_TT->Integral() << std::endl;
  
  TH1F* RSQ_22_TT = new TH1F( *TTjets[5] );
  std::cout << "TTJets RSQ 2BOX: " << RSQ_22_TT->Integral() << std::endl;
  
  TH1F* RSQ_11_TT = new TH1F( *TTjets[3] );
  std::cout << "TTJets RSQ 1BOX: " << RSQ_11_TT->Integral() << std::endl;
  
  TH1F* RSQ_00_TT = new TH1F( *TTjets[1] );
  std::cout << "TTJets RSQ 0BOX: " << RSQ_00_TT->Integral() << std::endl;

  
  TH1F* TT_MR0 =new TH1F( TT->PlotMR_0Box() );
  std::cout << "blah 0" << std::endl;
  TH1F* TT_R20 =new TH1F( TT->PlotRSQ_0Box() );
  std::cout << "blah 1" << std::endl;
  TH1F* TT_MR1 =new TH1F( TT->PlotMR_1Box() );
  std::cout << "blah 2" << std::endl;
  TH1F* TT_R21 =new TH1F( TT->PlotRSQ_1Box() );
  TH1F* TT_MR2 =new TH1F( TT->PlotMR_2Box() );
  TH1F* TT_R22 =new TH1F( TT->PlotRSQ_2Box() );

  ////////////////////////////////////////////////////////////////
  //////////////////// DATA//////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  //const char* data_file = "/afs/cern.ch/work/c/cpena/DarkMatter/CMSSW_5_2_3/src/VecbosApp/53X/HTMHT_Run2012A_ILV/out/HTMHT_ILV_Run2012AB.root";
  
  const char* data_file = "/media/data/cmorgoth/Data/DMData/HTMHT_ILV_Run2012AB.root";
  
  //const char* data_file = "/afs/cern.ch/work/c/cpena/DarkMatter/CMSSW_5_2_3/src/VecbosApp/afs/DoubleMu-Run2012AB/out/DoubleMu_ILV_AB.root";
  
  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  data->PrintEvents();
  
  TH1F* MR_22_data = new TH1F( data->PlotMR_2Box() );
  std::cout << "Data MR 2BOX: " << MR_22_data->Integral() << std::endl;
  MR_22_data->Sumw2();
    
  TH1F* MR_11_data = new TH1F( data->PlotMR_1Box() );
  std::cout << "Data MR 1BOX: " << MR_11_data->Integral() << std::endl;
  MR_11_data->Sumw2();
    
  TH1F* MR_00_data = new TH1F( data->PlotMR_0Box() );
  std::cout << "Data MR 0BOX: " << MR_00_data->Integral() << std::endl;
  MR_00_data->Sumw2();
    
  TH1F* RSQ_22_data = new TH1F( data->PlotRSQ_2Box() );
  std::cout << "Data RSQ 2BOX: " << RSQ_22_data->Integral() << std::endl;
  RSQ_22_data->Sumw2();
    
  TH1F* RSQ_11_data = new TH1F( data->PlotRSQ_1Box() );
  std::cout << "Data RSQ 1BOX: " << RSQ_11_data->Integral() << std::endl;
  RSQ_11_data->Sumw2();
  
  TH1F* RSQ_00_data = new TH1F( data->PlotRSQ_0Box() );
  std::cout << "Data RSQ 0BOX: " << RSQ_00_data->Integral() << std::endl;
  RSQ_00_data->Sumw2();
  

  TH1F* sum0BoxMR = new TH1F( *Wjets[0] );
  sum0BoxMR->Add(dy_jets[0]);
  sum0BoxMR->Add(TTjets[0]);
  sum0BoxMR->Add(Zjets[0]);

  float sumerr = .0;
  float sum = .0;
  for(int i = 1; i < sum0BoxMR->GetSize()-1; i++){
    std::cout << "error bin: " << i << " is: " << sum0BoxMR->GetBinError(i) <<  " " << sum0BoxMR->GetBinContent(i) << std::endl;
    sumerr += sum0BoxMR->GetBinError(i)*sum0BoxMR->GetBinError(i);
    sum += sum0BoxMR->GetBinContent(i);
  } 
  
  double ERROR = 0.0;
  double Int_Gral = sum0BoxMR->IntegralAndError(1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl; 
  std::cout << "errr: " << sqrt(sumerr) << std::endl;
  std::cout << "tot: " << sum << std::endl;

  TH1F* sum0BoxR2 = new TH1F( *Wjets[1] );
  sum0BoxR2->Add(dy_jets[1]);
  sum0BoxR2->Add(TTjets[1]);
  sum0BoxR2->Add(Zjets[1]);

  sumerr = 0;
  sum = .0;
  for(int i = 1; i < sum0BoxR2->GetSize()-1; i++){
    std::cout << "error bin: " << i << " is: " << sum0BoxR2->GetBinError(i) <<  " " << sum0BoxR2->GetBinContent(i) << std::endl;
    sumerr += sum0BoxR2->GetBinError(i)*sum0BoxR2->GetBinError(i);
    sum += sum0BoxR2->GetBinContent(i);
  }
  ERROR = 0.0;
  Int_Gral = sum0BoxR2->IntegralAndError(1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "errr: " << sqrt(sumerr) << std::endl;
  std::cout << "tot: " << sum << std::endl;

  TH1F* sum1BoxMR = new TH1F( *Wjets[2] );
  sum1BoxMR->Add(dy_jets[2]);
  sum1BoxMR->Add(TTjets[2]);
  sum1BoxMR->Add(Zjets[2]);

  sumerr = .0;
  sum = .0;
  for(int i = 1; i < sum1BoxMR->GetSize()-1; i++){
    std::cout << "error bin: " << i << " is: " << sum1BoxMR->GetBinError(i) <<  " " << sum1BoxMR->GetBinContent(i) << std::endl;
    sumerr += sum1BoxMR->GetBinError(i)*sum1BoxMR->GetBinError(i);
    sum += sum1BoxMR->GetBinContent(i);
  }
  ERROR = 0.0;
  Int_Gral = sum1BoxMR->IntegralAndError(1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "errr: " << sqrt(sumerr) << std::endl;
  std::cout << "tot: " << sum << std::endl;

  TH1F* sum2BoxMR = new TH1F( *Wjets[4] );
  sum2BoxMR->Add(dy_jets[4]);
  sum2BoxMR->Add(TTjets[4]);
  sum2BoxMR->Add(Zjets[4]);

  sumerr = .0;
  sum = .0;
  for(int i = 1; i < sum2BoxMR->GetSize()-1; i++){
    std::cout << "error bin: " << i << " is: " << sum2BoxMR->GetBinError(i) <<  " " << sum2BoxMR->GetBinContent(i) << std::endl;
    sumerr += sum2BoxMR->GetBinError(i)*sum2BoxMR->GetBinError(i);
    sum += sum2BoxMR->GetBinContent(i);
  }
  ERROR = 0.0;
  Int_Gral = sum2BoxMR->IntegralAndError(1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "errr: " << sqrt(sumerr) << std::endl;
  std::cout << "tot: " << sum << std::endl;


  ///////////////////////////////////////////
  /////////// individual error /////////////
  //////////////////////////////////////////
  float sumw_W[6], sumw_Z[6], sumw_TT[6], sumw_DY[6];
  for(int i = 0; i < 6; i++){
    sumw_W[i] = 0.0;
    sumw_TT[i] = 0.0;
    sumw_Z[i] = 0.0;
    sumw_DY[i] = 0.0;
  }
  
  for(int i = 0; i < 6; i++){
    for(int j = 1; j < Wjets[i]->GetSize()-1; j++){
      std::cout << "Wjets "  << i << " " << Wjets[i]->GetBinError(j) << " " << Wjets[i]->GetBinContent(j) << std::endl;
      std::cout << "Zjets "  << i << " " << Zjets[i]->GetBinError(j) << " " << Zjets[i]->GetBinContent(j) << std::endl;
      std::cout << "TTjets "  << i << " " << TTjets[i]->GetBinError(j) << " " << TTjets[i]->GetBinContent(j) << std::endl;
      std::cout << "DYjets "  << i << " " << dy_jets[i]->GetBinError(j) << " " << dy_jets[i]->GetBinContent(j) << std::endl;
      sumw_W[i] += Wjets[i]->GetBinError(j)*Wjets[i]->GetBinError(j);
      sumw_TT[i] += TTjets[i]->GetBinError(j)*TTjets[i]->GetBinError(j);
      sumw_Z[i] += Zjets[i]->GetBinError(j)*Zjets[i]->GetBinError(j);
      sumw_DY[i] += dy_jets[i]->GetBinError(j)*dy_jets[i]->GetBinError(j);
    }
    std::cout << "W: " << sqrt(sumw_W[i]) << std::endl;
    std::cout << "Z: " << sqrt(sumw_Z[i]) << std::endl;
    std::cout << "TT: " << sqrt(sumw_TT[i]) << std::endl;
    std::cout << "DY: " << sqrt(sumw_DY[i]) << std::endl;
  }
  
  
  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  ////////////////////////////////////////////////////////
  /////////////////RSQ O muon Box Stack plots////////////
  //////////////////////////////////////////////////////
  
  THStack* stack1 = new THStack("stack1", "");
  
  float dataMC_ratio = MR_00_data->GetBinContent(1)/( MR_00->GetBinContent(1) + MR_0->GetBinContent(1) + MR_00_TT->GetBinContent(1) + MR_dy_00->GetBinContent(1) );
  
  std::cout <<  "  data/mc first bin: " <<  dataMC_ratio << std::endl;
  
  //RSQ_00->Scale( dataMC_ratio );
  //RSQ_0->Scale( dataMC_ratio );
  //RSQ_00_TT->Scale( dataMC_ratio );
  //RSQ_dy_00->Scale( dataMC_ratio );
  
  RSQ_00_data->SetMarkerStyle(20);
  RSQ_00_data->SetLineColor(1);
  RSQ_00_data->SetMarkerSize(1.5);
  
  TT_R20->SetFillColor(kPink+9);
  dy_R20->SetFillColor(kViolet+9);
  Z_R20->SetFillColor(kYellow-4);
  W_R20->SetFillColor(kSpring+4);

  leg = new TLegend(0.7,0.7,0.9,0.92);
  
  leg->AddEntry(W_R20,"W + jets","f");
  leg->AddEntry(Z_R20,"Z(#nu#bar{#nu}) + jets","f");
  leg->AddEntry(TT_R20,"t #bar{t} + jets","f");
  leg->AddEntry(dy_R20,"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(RSQ_00_data,"Data","lep");
  
  //leg->SetHeader("R^{2} Signal Region");
  
  W_R20->SetTitle("");
  W_R20->SetStats(0);
  TT_R20->SetTitle("");
  TT_R20->SetStats(0);
  dy_R20->SetTitle("");
  dy_R20->SetStats(0);
  Z_R20->SetTitle("");
  Z_R20->SetStats(0);
  RSQ_00_data->SetTitle("");
  RSQ_00_data->SetStats(0);
  
  stack1->Add(dy_R20);//DY
  stack1->Add(TT_R20);//TTbar
  stack1->Add(Z_R20);//Wjets
  stack1->Add(W_R20);//ZJets
  
  stack1->Draw();
  //( (TAxis*)( stack1->GetXaxis() ) )->SetLimits(0.5, 2.5);
  ( (TAxis*)( stack1->GetXaxis() ) )->SetTitle("R^{2}");
  stack1->Draw();
  RSQ_00_data->Draw("same");
  
    
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/RSQ_0BOX_Stack.pdf");
  C1->SaveAs("StackPlots/RSQ_0BOX_Stack.png");
  //delete leg;
  
  TH1F* aux2 = new TH1F( *RSQ_dy_00 );
  aux2->Sumw2();
  aux2->Add( RSQ_00_TT, 1);
  aux2->Add( RSQ_0, 1);
  aux2->Add( RSQ_00, 1);
  RatioPlotsV2(stack1, aux2, RSQ_00_data, "MC 0 #mu BOX", "Data 0 #mu BOX", "StackPlots/Data_MC_RSQ_0BOX", "RSQ", leg);
  delete leg;
  delete aux2;
  
  /////////////////////////////////////////////////////////
  /////////////////MR O muon Box Stacl plots//////////////
  ///////////////////////////////////////////////////////
  THStack* stackMR_0Box = new THStack("stackMR_0Box", "");
  
  std::cout << MR_00->GetBinContent(0)  << "  " << RSQ_00->GetBinContent(1) << std::endl;
  std::cout << MR_00->GetBinContent(2)  << "  " << RSQ_00->GetBinContent(2) << std::endl;
  
  //float dataMC_ratio = MR_00_data->GetBinContent(1)/( MR_00->GetBinContent(1) + MR_0->GetBinContent(1) + MR_00_TT->GetBinContent(1) + MR_dy_00->GetBinContent(1) );
  
  //std::cout <<  "  data/mc first bin: " <<  dataMC_ratio << std::endl;
  
  //MR_00->Scale( dataMC_ratio );
  //MR_0->Scale( dataMC_ratio );
  //MR_00_TT->Scale( dataMC_ratio );
  //MR_dy_00->Scale( dataMC_ratio );
  
  TT_MR0->SetTitle("");
  TT_MR0->SetStats(0);
  Z_MR0->SetTitle("");
  Z_MR0->SetStats(0);
  W_MR0->SetTitle("");
  W_MR0->SetStats(0);
  dy_MR0->SetTitle("");
  dy_MR0->SetStats(0);
  MR_00_data->SetTitle("");
  MR_00_data->SetStats(0);
  
  MR_00_data->SetMarkerStyle(20);
  MR_00_data->SetLineColor(1);
  MR_00_data->SetMarkerSize(1.5);
  
  TT_MR0->SetFillColor(kPink+9);
  dy_MR0->SetFillColor(kViolet+9);
  W_MR0->SetFillColor(kYellow-4);
  Z_MR0->SetFillColor(kSpring+4);
  
  
  leg = new TLegend(0.70,0.7,0.9,0.90);
  
  leg->AddEntry(W_MR0,"W + jets","f");
  leg->AddEntry(Z_MR0,"Z(#nu#nu) + jets","f");
  leg->AddEntry(TT_MR0,"t #bar{t} + jets","f");
  leg->AddEntry(dy_MR0,"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(MR_00_data,"Data","lep");
  
  leg->SetHeader("M_{R} Signal Region");
  
  stackMR_0Box->Add(dy_MR0);//DY
  stackMR_0Box->Add(TT_MR0);//TTbar
  stackMR_0Box->Add(W_MR0);//Wjets
  stackMR_0Box->Add(Z_MR0);//ZJets
  
  stackMR_0Box->Draw();
  ( (TAxis*)( stackMR_0Box->GetXaxis() ) )->SetTitle("M_{R} (GeV)");
  stackMR_0Box->Draw();
  MR_00_data->Draw("same");
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/MR_0BOX_Stack.pdf");
  C1->SaveAs("StackPlots/MR_0BOX_Stack.png");

  aux2 = new TH1F( *MR_dy_00 );
  aux2->Sumw2();
  aux2->Add( MR_00_TT, 1);
  aux2->Add( MR_0, 1);
  aux2->Add( MR_00, 1);
  RatioPlotsV2( stackMR_0Box, aux2, MR_00_data, "MC 0 #mu BOX", "Data 0 #mu BOX", "StackPlots/Data_MC_MR_0BOX", "MR", leg);
  delete leg;
  delete aux2; 
  
  ////////////////////////////////////////////////////////
  /////////////////RSQ 1 muon Box Stacl plots////////////
  //////////////////////////////////////////////////////
  THStack* stackRSQ_1Box = new THStack("stackRSQ_1Box", "");
  
  // dataMC_ratio = MR_11_data->GetBinContent(1)/( MR_11->GetBinContent(1) + MR_1->GetBinContent(1) + MR_11_TT->GetBinContent(1) + MR_dy_11->GetBinContent(1) );

  //RSQ_11->Scale( dataMC_ratio );
  //RSQ_1->Scale( dataMC_ratio );
  //RSQ_11_TT->Scale( dataMC_ratio );
  //RSQ_dy_11->Scale( dataMC_ratio );

  TT_R21->SetTitle("");
  TT_R21->SetStats(0);
  dy_R21->SetTitle("");
  dy_R21->SetStats(0);
  W_R21->SetTitle("");
  W_R21->SetStats(0);
  Z_R21->SetTitle("");
  Z_R21->SetStats(0);
  RSQ_11_data->SetTitle("");
  RSQ_11_data->SetStats(0);
  
  RSQ_11_data->SetMarkerStyle(20);
  RSQ_11_data->SetLineColor(1);
  RSQ_11_data->SetMarkerSize(1.5);
  
  TT_R21->SetFillColor(kPink+9);
  dy_R21->SetFillColor(kViolet+9);
  W_R21->SetFillColor(kYellow-4);
  Z_R21->SetFillColor(kSpring+4);

  leg = new TLegend(0.7,0.7,0.9,0.92);
  
  leg->AddEntry(W_R21,"W + jets","f");
  leg->AddEntry(dy_R21,"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(TT_R21,"t #bar{t} + jets","f");
  leg->AddEntry(RSQ_11_data,"Data","lep");
  leg->SetHeader("R^{2} 1 #mu Control Region");
  
  stackRSQ_1Box->Add(dy_R21);//DY
  stackRSQ_1Box->Add(TT_R21);//TTbar
  stackRSQ_1Box->Add(W_R21);//Wjets
    
  stackRSQ_1Box->Draw();
  //( (TAxis*)( stackRSQ_1Box->GetXaxis() ) )->SetLimits(0.5, 2.5);
  ( (TAxis*)( stackRSQ_1Box->GetXaxis() ) )->SetTitle("R^{2} ");
  stackRSQ_1Box->Draw();
  RSQ_11_data->Draw("same");
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/RSQ_1BOX_Stack.pdf");
  C1->SaveAs("StackPlots/RSQ_1BOX_Stack.png");
      
  aux2 = new TH1F( *RSQ_dy_11 );
  aux2->Sumw2();
  aux2->Add( RSQ_11_TT, 1);
  aux2->Add( RSQ_1, 1);
  aux2->Add( RSQ_11, 1);
  RatioPlotsV2( stackRSQ_1Box, aux2, RSQ_11_data, "MC 1 #mu BOX", "Data 1 #mu BOX", "StackPlots/Data_MC_RSQ_1BOX", "RSQ", leg);
  delete leg;
  delete aux2;

  /////////////////////////////////////////////////////////
  /////////////////MR 1 muon Box Stacl plots//////////////
  ///////////////////////////////////////////////////////
  THStack* stackMR_1Box = new THStack("stackMR_1Box", "");
  
  //MR_11->Scale( dataMC_ratio );
  //MR_1->Scale( dataMC_ratio );
  //MR_11_TT->Scale( dataMC_ratio );
  //MR_dy_11->Scale( dataMC_ratio );

  W_MR1->SetTitle("");
  W_MR1->SetStats(0);
  Z_MR1->SetTitle("");
  Z_MR1->SetStats(0);
  dy_MR1->SetTitle("");
  dy_MR1->SetStats(0);
  TT_MR1->SetTitle("");
  TT_MR1->SetStats(0);
  MR_11_data->SetTitle("");
  MR_11_data->SetStats(0);
  
  MR_11_data->SetMarkerStyle(20);
  MR_11_data->SetLineColor(1);
  MR_11_data->SetMarkerSize(1.5);
  
  TT_MR1->SetFillColor(kPink+9);
  dy_MR1->SetFillColor(kViolet+9);
  W_MR1->SetFillColor(kYellow-4);
  Z_MR1->SetFillColor(kSpring+4);

  leg = new TLegend(0.7,0.7,0.9,0.90);
  
  leg->AddEntry(W_MR1,"W + jets","f");
  leg->AddEntry(dy_MR1,"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(TT_MR1,"t #bar{t} + jets","f");
  leg->AddEntry(Z_MR1,"Z(#nu#bar{#nu}) + jets","f");
  leg->AddEntry(MR_11_data,"Data","lep");
  
  //leg->SetHeader("M_{R} 1 #mu Control Region");
  
  stackMR_1Box->Add(dy_MR1);//DY
  stackMR_1Box->Add(TT_MR1);//TTbar
  stackMR_1Box->Add(W_MR1);//Wjets
  stackMR_1Box->Add(Z_MR1);//ZJets
  
  stackMR_1Box->Draw();
  ( (TAxis*)( stackMR_1Box->GetXaxis() ) )->SetTitle("M_{R} (GeV)");
  stackMR_1Box->Draw();
  MR_11_data->Draw("same");
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/MR_1BOX_Stack.pdf");
  C1->SaveAs("StackPlots/MR_1BOX_Stack.png");
    
  aux2 = new TH1F( *MR_dy_11 );
  aux2->Sumw2();
  aux2->Add( MR_11_TT, 1);
  aux2->Add( MR_1, 1);
  aux2->Add( MR_11, 1);
  RatioPlotsV2( stackMR_1Box, aux2, MR_11_data, "MC 1 #mu BOX", "Data 1 #mu BOX", "StackPlots/Data_MC_MR_1BOX", "MR", leg);
  delete leg;
  delete aux2;

  ////////////////////////////////////////////////////////
  /////////////////RSQ 2 muon Box Stack plots////////////
  //////////////////////////////////////////////////////
  THStack* stackRSQ_2Box = new THStack("stackRSQ_2Box", "");
  
  //dataMC_ratio = MR_22_data->GetBinContent(1)/(  MR_22_TT->GetBinContent(1) + MR_dy_22->GetBinContent(1) );
  
  //RSQ_22->Scale( dataMC_ratio );
  //RSQ_2->Scale( dataMC_ratio );
  //RSQ_22_TT->Scale( dataMC_ratio );
  //RSQ_dy_22->Scale( dataMC_ratio );
  
  //RSQ_2->SetTitle("");
  //RSQ_2->SetStats(0);
  //RSQ_22->SetTitle("");
  //RSQ_22->SetStats(0);
  TT_R22->SetTitle("");
  TT_R22->SetStats(0);
  dy_R22->SetTitle("");
  dy_R22->SetStats(0);
  RSQ_22_data->SetTitle("");
  RSQ_22_data->SetStats(0);
  
  RSQ_22_data->SetMarkerStyle(20);
  RSQ_22_data->SetLineColor(1);
  RSQ_22_data->SetMarkerSize(1.5);
  
  TT_R22->SetFillColor(kPink+9);
  dy_R22->SetFillColor(kViolet+9);
  
  leg = new TLegend(0.7,0.7,0.9,0.92);
  
  leg->AddEntry(TT_R22,"t #bar{t} + jets","f");
  leg->AddEntry(dy_R22,"Z/#gamma^{*}(ll) + jets","f");
  //leg->AddEntry(RSQ_22,"Z(#nu#bar{#nu}) + jets","f");
  //leg->AddEntry(RSQ_2,"W + jets","f");
  leg->AddEntry(RSQ_22_data,"Data","lep");
  
  //leg->SetHeader("R^{2} 2 #mu Control Region");
  
  stackRSQ_2Box->Add(TT_R22);//TTbar
  stackRSQ_2Box->Add(dy_R22);//DY
  //stackRSQ_2Box->Add(RSQ_11);//ZJets
  //stackRSQ_2Box->Add(RSQ_1);//Wjets
  
  stackRSQ_2Box->Draw();
  //( (TAxis*)( stackRSQ_2Box->GetXaxis() ) )->SetLimits(0.5, 2.5);
  ( (TAxis*)( stackRSQ_2Box->GetXaxis() ) )->SetTitle("R^{2} ");
  stackRSQ_2Box->Draw();
  RSQ_22_data->Draw("same");
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/RSQ_2BOX_Stack_NoClean.pdf");
  C1->SaveAs("StackPlots/RSQ_2BOX_Stack_NoClean.png");
  

  aux2 = new TH1F( *RSQ_dy_22 );
  aux2->Sumw2();
  aux2->Add( RSQ_22_TT, 1);
  //aux2->Add( RSQ_2, 1);
  //aux2->Add( RSQ_22, 1);                                                                                                        
  RatioPlotsV2( stackRSQ_2Box, aux2, RSQ_22_data, "MC 2 #mu BOX", "Data 2 #mu BOX", "StackPlots/Data_MC_RSQ_2BOX", "RSQ", leg);
  delete leg;
  delete aux2;

  /////////////////////////////////////////////////////////
  /////////////////MR 2 muon Box Stack plots//////////////
  ///////////////////////////////////////////////////////
  THStack* stackMR_2Box = new THStack("stackMR_2Box", "");
  
  //MR_22->Scale( dataMC_ratio );
  //MR_2->Scale( dataMC_ratio );
  //MR_22_TT->Scale( dataMC_ratio );
  //MR_dy_22->Scale( dataMC_ratio );

  //MR_2->SetTitle("");
  //MR_2->SetStats(0);
  //MR_22->SetTitle("");
  //MR_22->SetStats(0);
  TT_MR2->SetTitle("");
  TT_MR2->SetStats(0);
  dy_MR2->SetTitle("");
  dy_MR2->SetStats(0);
  MR_22_data->SetTitle("");
  MR_22_data->SetStats(0);
  
  MR_22_data->SetMarkerStyle(20);
  MR_22_data->SetLineColor(1);
  MR_22_data->SetMarkerSize(1.5);
  
  TT_MR2->SetFillColor(kPink+9);
  dy_MR2->SetFillColor(kViolet+9);
  //MR_1->SetFillColor(6);
  //MR_11->SetFillColor(kRed);

  leg = new TLegend(0.7,0.7,0.9,0.9);
  
  leg->AddEntry(TT_MR2,"t #bar{t} + jets","f");
  leg->AddEntry(dy_MR2,"Z/#gamma^{*}(ll) + jets","f");
  //leg->AddEntry(MR_22,"Z(#nu#bar{#nu}) + jets","f");
  //leg->AddEntry(MR_2,"W + jets","f");
  leg->AddEntry(MR_22_data,"Data","lep");
  
  //leg->SetHeader("M_{R} 2 #mu Control Region");
  
  stackMR_2Box->Add(TT_MR2);//TTbar
  stackMR_2Box->Add(dy_MR2);//DY
  //stackMR_1Box->Add(MR_22);//ZJets
  //stackMR_1Box->Add(MR_2);//Wjets
  
  stackMR_2Box->Draw();
  ( (TAxis*)( stackMR_2Box->GetXaxis() ) )->SetTitle("R^{2} ");
  
  stackMR_2Box->Draw();
  MR_22_data->Draw("same");
  leg->Draw();
  C1->SetLogy();
  C1->SaveAs("StackPlots/MR_2BOX_Stack.pdf");
  C1->SaveAs("StackPlots/MR_2BOX_Stack.png");
  
  aux2 = new TH1F( *MR_dy_22 );
  aux2->Sumw2();
  aux2->Add( MR_22_TT, 1);
  //aux2->Add( MR_22, 1);
  //aux2->Add( MR_2, 1);
  RatioPlotsV2( stackMR_2Box, aux2, MR_22_data, "MC 2 #mu BOX", "Data 2 #mu BOX", "StackPlots/Data_MC_MR_2BOX", "MR", leg);
  delete leg;
  delete aux2;

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ////////////////////Saving Histos in File////////////////////////////
  ///////////////////////  Finishing  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  TFile* f2 = new TFile("stack1_TightBtag.root","RECREATE");
  
  stack1->Write();
  stackMR_0Box->Write();

  stackRSQ_1Box->Write();
  stackMR_1Box->Write();
  stackRSQ_2Box->Write();
  stackMR_2Box->Write();
  for(int i = 0; i < 6; i++){
    (Wjets[i])->Write();
    (dy_jets[i])->Write();
    (Zjets[i])->Write();
    (TTjets[i])->Write();
    
  }

  MR_22_data->Write("MR_data_2box");
  MR_11_data->Write("MR_data_1box");
  MR_00_data->Write("MR_data_0box");
  RSQ_22_data->Write("R2_data_2box");
  RSQ_11_data->Write("R2_data_1box");
  RSQ_00_data->Write("R2_data_0box");


  
  f2->Close();
  delete C1;
  delete W;
  delete TT;
  delete Z;
  delete dy;
  delete data;



}
