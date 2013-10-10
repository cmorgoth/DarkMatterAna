#include "DM_METPlots.hh"
#include "DM_1DRatio.hh"

void CreateMetPlots(){

  TLegend* lege;
  
  TCanvas* C0 = new TCanvas("C0", "C0", 1024, 1024);
  C0->cd();

  int bL, bM, bT;
  bL = bM = 0;
  bT = 0;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
  //////////////////////////////////////////// W + Jets //////////////////////////////////////////////////////////////   
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   

  WJetsHTBins* W = new WJetsHTBins( 2 );
  W->SetBtagCut(bL,bM,bT);
  std::vector<TH1F*> MET_W = W->PlotMETmag();
  
  int ctr = 0;
  TString metpicname;
  for( std::vector<TH1F*>::iterator itW = MET_W.begin(); itW != MET_W.end(); ++itW ){
    (*itW)->Draw();
    metpicname = TString(Form("METplots/W_metXbox%d.pdf",ctr));
    C0->SaveAs(metpicname);
    metpicname = TString(Form("METplots/W_metXbox%d.png",ctr));
    C0->SaveAs(metpicname);
    ctr++;
  }
  
  /////////////////////////////////////////////
  ///////////////ZJetsNuNu ANA/////////////////
  /////////////////////////////////////////////
  
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  std::vector<TH1F*> MET_Z = Z->PlotMETmag();
  
  ctr = 0;
  for( std::vector<TH1F*>::iterator itZ = MET_Z.begin(); itZ != MET_Z.end(); ++itZ ){
    (*itZ)->Draw();
    metpicname = TString(Form("METplots/Z_metXbox%d.pdf",ctr));
    C0->SaveAs(metpicname);
    metpicname = TString(Form("METplots/Z_metXbox%d.png",ctr));
    C0->SaveAs(metpicname);
    ctr++;
  }

  
  //////////////////////////////////////////////////////
  ///////////////// Drell-Yan///////////////////////////
  /////////////////////////////////////////////////////
  
  DY* dy = new DY( 2 );
  dy->SetBtagCut(bL,bM,bT);
  std::vector<TH1F*> MET_DY = dy->PlotMETmag();
  
  
  ctr = 0;
  for( std::vector<TH1F*>::iterator itDY = MET_DY.begin(); itDY != MET_DY.end(); ++itDY ){
    (*itDY)->Draw();
    metpicname = TString(Form("METplots/DY_metXbox%d.pdf",ctr));
    C0->SaveAs(metpicname);
    metpicname = TString(Form("METplots/DY_metXbox%d.png",ctr));
    C0->SaveAs(metpicname);
    ctr++;
  }
  
  //C0->SaveAs("MET_DY.pdf");
  
  //////////////////////////////////////////////////
  /////////////////////////////////////////////////
  /////////////////TTbar + Jets////////////////////
  ////////////////////////////////////////////////
  ///////////////////////////////////////////////
  
  //const char* ttjets_file = "root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/PFJets80/TTJets80_Tot.root";
  //const char* ttjets_file = "root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/NoTauVeto/TTJets_NoTauVeto_Total.root";
  
  TTJets* TT = new TTJets( 2 );
  TT->SetBtagCut(bL,bM,bT);
  std::vector<TH1F*> MET_TT = TT->PlotMETmag();
  
  ctr = 0;
  for( std::vector<TH1F*>::iterator itTT = MET_TT.begin(); itTT != MET_TT.end(); ++itTT ){
    (*itTT)->Draw();
    metpicname = TString(Form("METplots/TT_metXbox%d.pdf",ctr));
    C0->SaveAs(metpicname);
    metpicname = TString(Form("METplots/TT_metXbox%d.png",ctr));
    C0->SaveAs(metpicname);
    ctr++;
  }
  
  ////////////////////////////////////////////////////////////////
  //////////////////// DATA//////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  //const char* data_file = "/media/data/cmorgoth/Data/DMData/HTMHT_ILV_Run2012AB.root";
  const char* data_file = "/media/data/cmorgoth/Data/DMData/FullHTMHTRereco/HTMHT_ABCD_FullLumi20128TeV.root";
    
  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  std::vector<TH1F*> MET_Data = data->PlotMETmag();
  
  ctr = 0;
  for( std::vector<TH1F*>::iterator itData = MET_Data.begin(); itData != MET_Data.end(); ++itData ){
    (*itData)->Draw();
    metpicname = TString(Form("METplots/Data_metXbox%d.pdf",ctr));
    C0->SaveAs(metpicname);
    metpicname = TString(Form("METplots/Data_metXbox%d.png",ctr));
    C0->SaveAs(metpicname);
    ctr++;
  }

  
  THStack* stack1 = new THStack("stack1", "");
  
  MET_TT[1]->SetFillColor(kPink+9);
  MET_TT[1]->SetStats(0);
  MET_TT[1]->SetYTitle("Events/20 GeV");
  MET_DY[1]->SetFillColor(kViolet+9);
  MET_DY[1]->SetStats(0);
  MET_DY[1]->SetYTitle("Events/20 GeV");
  MET_W[1]->SetFillColor(kYellow-4);
  MET_W[1]->SetStats(0);
  MET_W[1]->SetYTitle("Events/20 GeV");
  MET_Z[1]->SetFillColor(kSpring+4);
  MET_Z[1]->SetStats(0);
  MET_Z[1]->SetYTitle("Events/20 GeV");
  
  stack1->Add(MET_DY[1]);
  stack1->Add(MET_TT[1]);                                                                              
  stack1->Add(MET_W[1]);     
  stack1->Add(MET_Z[1]);
    
  TH1F* aux2 = new TH1F(*MET_DY[1]);
  aux2->Sumw2();
  aux2->Add( MET_TT[1], 1);
  aux2->Add( MET_W[1], 1);
  aux2->Add( MET_Z[1], 1);
  
  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[1]->Sumw2();
  lege->AddEntry(MET_W[1],"W + jets","f");
  lege->AddEntry(MET_Z[1],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[1],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[1],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[1],"Data","lep");


  stack1->Draw();
  stack1->GetXaxis()->SetTitle("Events/20 GeV");
  gPad->Modified();
  
  MET_Data[1]->Sumw2();
  MET_Data[1]->SetMarkerStyle(20);
  MET_Data[1]->SetLineColor(1);
  MET_Data[1]->SetMarkerSize(1.5);
  MET_Data[1]->SetStats(0);
  //MET_Data[1]->SetYTitle("Events/20 GeV");
  MET_Data[1]->Draw("same");

  std::cout << "deb 0" << std::endl;
  RatioPlotsV2( stack1, aux2, MET_Data[1], "MET MC 0 #mu BOX", "MET Data 0 #mu BOX", "METplots/Data_MC_MET_01BOX","MET", lege);
    
  std::cout << "deb 1" << std::endl;
  lege->Draw();
  C0->SetLogy();
  C0->SaveAs("METplots/Met1Stack0Mu.pdf");
  C0->SaveAs("METplots/Met1Stack0Mu.png");

  delete lege;
  delete aux2;

  THStack* stack2 = new THStack("stack2", "");

  MET_TT[0]->SetFillColor(kPink+9);
  MET_TT[0]->SetStats(0);
  MET_DY[0]->SetFillColor(kViolet+9);
  MET_DY[0]->SetStats(0);
  MET_W[0]->SetFillColor(kYellow-4);
  MET_W[0]->SetStats(0);
  MET_Z[0]->SetFillColor(kSpring+4);
  MET_Z[0]->SetStats(0);

  stack2->Add(MET_DY[0]);
  stack2->Add(MET_TT[0]);
  stack2->Add(MET_W[0]);
  stack2->Add(MET_Z[0]);
  
  aux2 = new TH1F(*MET_DY[0]);
  aux2->Sumw2();
  aux2->Add( MET_TT[0], 1);
  aux2->Add( MET_W[0], 1);
  aux2->Add( MET_Z[0], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[0]->Sumw2();
  lege->AddEntry(MET_W[0],"W + jets","f");
  lege->AddEntry(MET_Z[0],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[0],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[0],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[0],"Data","lep");
  
  stack2->Draw();
  MET_Data[0]->Sumw2();
  MET_Data[0]->SetMarkerStyle(20);
  MET_Data[0]->SetLineColor(1);
  MET_Data[0]->SetMarkerSize(1.5);
  MET_Data[0]->SetStats(0);
  MET_Data[0]->Draw("same");
  
  RatioPlotsV2( stack2, aux2, MET_Data[0], "MET MC 0 #mu BOX", "MET Data 0 #mu BOX", "METplots/Data_MC_MET_00BOX","MET", lege);

  C0->SaveAs("METplots/Met0Stack0Mu.pdf");
  C0->SaveAs("METplots/Met0Stack0Mu.png");

  delete aux2;
  delete lege;

  THStack* stack3 = new THStack("stack3", "");

  MET_TT[3]->SetFillColor(kPink+9);
  MET_TT[3]->SetStats(0);
  MET_DY[3]->SetFillColor(kViolet+9);
  MET_DY[3]->SetStats(0);
  MET_W[3]->SetFillColor(kYellow-4);
  MET_W[3]->SetStats(0);
  MET_Z[3]->SetFillColor(kSpring+4);
  MET_Z[3]->SetStats(0);

  stack3->Add(MET_DY[3]);
  stack3->Add(MET_TT[3]);
  stack3->Add(MET_W[3]);
  stack3->Add(MET_Z[3]);

  aux2 = new TH1F(*MET_DY[3]);
  aux2->Sumw2();
  aux2->Add( MET_TT[3], 1);
  aux2->Add( MET_W[3], 1);
  aux2->Add( MET_Z[3], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[3]->Sumw2();
  lege->AddEntry(MET_W[3],"W + jets","f");
  lege->AddEntry(MET_Z[3],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[3],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[3],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[3],"Data","lep");

  stack3->Draw();
  MET_Data[3]->Sumw2();
  MET_Data[3]->SetMarkerStyle(20);
  MET_Data[3]->SetLineColor(1);
  MET_Data[3]->SetMarkerSize(1.5);
  MET_Data[3]->Draw("same");

  RatioPlotsV2( stack3, aux2, MET_Data[3], "MET MC 1 #mu BOX", "MET Data 1 #mu BOX", "METplots/Data_MC_MET_10BOX","MET", lege);

  C0->SaveAs("METplots/Met0Stack1Mu.pdf");
  C0->SaveAs("METplots/Met0Stack1Mu.png");
  
  delete aux2;
  delete lege;
  
  THStack* stack4 = new THStack("stack4", "");

  MET_TT[4]->SetFillColor(kPink+9);
  MET_TT[4]->SetStats(0);
  MET_DY[4]->SetFillColor(kViolet+9);
  MET_DY[4]->SetStats(0);
  MET_W[4]->SetFillColor(kYellow-4);
  MET_W[4]->SetStats(0);
  MET_Z[4]->SetFillColor(kSpring+4);
  MET_Z[4]->SetStats(0);

  stack4->Add(MET_DY[4]);
  stack4->Add(MET_TT[4]);
  stack4->Add(MET_W[4]);
  stack4->Add(MET_Z[4]);
  
  aux2 = new TH1F(*MET_DY[4]);
  aux2->Sumw2();
  aux2->Add( MET_TT[4], 1);
  aux2->Add( MET_W[4], 1);
  aux2->Add( MET_Z[4], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[4]->Sumw2();
  lege->AddEntry(MET_W[4],"W + jets","f");
  lege->AddEntry(MET_Z[4],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[4],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[4],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[4],"Data","lep");
  
  stack4->Draw();
  MET_Data[4]->Sumw2();
  MET_Data[4]->SetMarkerStyle(20);
  MET_Data[4]->SetLineColor(1);
  MET_Data[4]->SetMarkerSize(1.5);
  MET_Data[4]->Draw("same");

  RatioPlotsV2( stack4, aux2, MET_Data[4], "MET MC 1 #mu BOX", "MET Data 1 #mu BOX", "METplots/Data_MC_MET_11BOX","MET", lege);

  C0->SaveAs("METplots/Met1Stack1Mu.pdf");
  C0->SaveAs("METplots/Met1Stack1Mu.png");

  THStack* stack6 = new THStack("stack6", "");

  MET_TT[6]->SetFillColor(kPink+9);
  MET_TT[6]->SetStats(0);
  MET_DY[6]->SetFillColor(kViolet+9);
  MET_DY[6]->SetStats(0);
  MET_W[6]->SetFillColor(kYellow-4);
  MET_W[6]->SetStats(0);
  MET_Z[6]->SetFillColor(kSpring+4);
  MET_Z[6]->SetStats(0);

  stack6->Add(MET_DY[6]);
  stack6->Add(MET_TT[6]);
  stack6->Add(MET_W[6]);
  stack6->Add(MET_Z[6]);

  aux2 = new TH1F(*MET_DY[6]);
  aux2->Sumw2();
  aux2->Add( MET_TT[6], 1);
  aux2->Add( MET_W[6], 1);
  aux2->Add( MET_Z[6], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[6]->Sumw2();
  lege->AddEntry(MET_W[6],"W + jets","f");
  lege->AddEntry(MET_Z[6],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[6],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[6],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[6],"Data","lep");

  stack6->Draw();
  MET_Data[6]->Sumw2();
  MET_Data[6]->SetMarkerStyle(20);
  MET_Data[6]->SetLineColor(1);
  MET_Data[6]->SetMarkerSize(1.5);
  MET_Data[6]->Draw("same");

  RatioPlotsV2( stack6, aux2, MET_Data[6], "MET MC 2 #mu BOX", "MET Data 2 #mu BOX", "METplots/Data_MC_MET_20BOX","MET", lege);

  C0->SaveAs("METplots/Met0Stack2Mu.pdf");
  C0->SaveAs("METplots/Met0Stack2Mu.png");

  delete aux2;
  delete lege;

  THStack* stack7 = new THStack("stack7", "");

  MET_TT[7]->SetFillColor(kPink+9);
  MET_TT[7]->SetStats(0);
  MET_DY[7]->SetFillColor(kViolet+9);
  MET_DY[7]->SetStats(0);
  MET_W[7]->SetFillColor(kYellow-4);
  MET_W[7]->SetStats(0);
  MET_Z[7]->SetFillColor(kSpring+4);
  MET_Z[7]->SetStats(0);

  stack7->Add(MET_DY[7]);
  stack7->Add(MET_TT[7]);
  stack7->Add(MET_W[7]);
  stack7->Add(MET_Z[7]);
  
  aux2 = new TH1F(*MET_DY[7]);
  aux2->Sumw2();
  aux2->Add( MET_TT[7], 1);
  aux2->Add( MET_W[7], 1);
  aux2->Add( MET_Z[7], 1);
  
  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[7]->Sumw2();
  lege->AddEntry(MET_W[7],"W + jets","f");
  lege->AddEntry(MET_Z[7],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[7],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[7],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[7],"Data","lep");

  stack7->Draw();
  MET_Data[7]->Sumw2();
  MET_Data[7]->SetMarkerStyle(20);
  MET_Data[7]->SetLineColor(1);
  MET_Data[7]->SetMarkerSize(1.5);
  MET_Data[7]->Draw("same");

  RatioPlotsV2( stack7, aux2, MET_Data[7], "MET MC 2 #mu BOX", "MET Data 2 #mu BOX", "METplots/Data_MC_MET_21BOX","MET", lege);

  C0->SaveAs("METplots/Met1Stack2Mu.pdf");
  C0->SaveAs("METplots/Met1Stack2Mu.png");
  
  delete aux2;
  delete lege;

  THStack* stack9 = new THStack("stack9", "");

  MET_TT[9]->SetFillColor(kPink+9);
  MET_TT[9]->SetStats(0);
  MET_DY[9]->SetFillColor(kViolet+9);
  MET_DY[9]->SetStats(0);
  MET_W[9]->SetFillColor(kYellow-4);
  MET_W[9]->SetStats(0);
  MET_Z[9]->SetFillColor(kSpring+4);
  MET_Z[9]->SetStats(0);

  stack9->Add(MET_DY[9]);
  stack9->Add(MET_TT[9]);
  stack9->Add(MET_W[9]);
  stack9->Add(MET_Z[9]);

  aux2 = new TH1F(*MET_DY[9]);
  aux2->Sumw2();
  aux2->Add( MET_TT[9], 1);
  aux2->Add( MET_W[9], 1);
  aux2->Add( MET_Z[9], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[9]->Sumw2();
  lege->AddEntry(MET_W[9],"W + jets","f");
  lege->AddEntry(MET_Z[9],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[9],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[9],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[9],"Data","lep");

  stack9->Draw();
  MET_Data[9]->Sumw2();
  MET_Data[9]->SetMarkerStyle(20);
  MET_Data[9]->SetLineColor(1);
  MET_Data[9]->SetMarkerSize(1.5);
  MET_Data[9]->Draw("same");
  
  RatioPlotsV2( stack9, aux2, MET_Data[9], "NJETS MC 0 #mu BOX", "NJETS Data 0 #mu BOX", "METplots/Data_MC_NJETS_0BOX","NJETS", lege);

  delete aux2;
  delete lege;
  
  THStack* stack10 = new THStack("stack10", "");

  MET_TT[10]->SetFillColor(kPink+9);
  MET_TT[10]->SetStats(0);
  MET_DY[10]->SetFillColor(kViolet+9);
  MET_DY[10]->SetStats(0);
  MET_W[10]->SetFillColor(kYellow-4);
  MET_W[10]->SetStats(0);
  MET_Z[10]->SetFillColor(kSpring+4);
  MET_Z[10]->SetStats(0);

  stack10->Add(MET_DY[10]);
  stack10->Add(MET_TT[10]);
  stack10->Add(MET_W[10]);
  stack10->Add(MET_Z[10]);

  aux2 = new TH1F(*MET_DY[10]);
  aux2->Sumw2();
  aux2->Add( MET_TT[10], 1);
  aux2->Add( MET_W[10], 1);
  aux2->Add( MET_Z[10], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[10]->Sumw2();
  lege->AddEntry(MET_W[10],"W + jets","f");
  lege->AddEntry(MET_Z[10],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[10],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[10],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[10],"Data","lep");

  stack10->Draw();
  MET_Data[10]->Sumw2();
  MET_Data[10]->SetMarkerStyle(20);
  MET_Data[10]->SetLineColor(1);
  MET_Data[10]->SetMarkerSize(1.5);
  MET_Data[10]->Draw("same");

  RatioPlotsV2( stack10, aux2, MET_Data[10], "NJETS MC 1 #mu BOX", "NJETS Data 1 #mu BOX", "METplots/Data_MC_NJETS_1BOX","NJETS", lege);

  delete aux2;
  delete lege;

  THStack* stack11 = new THStack("stack11", "");

  MET_TT[11]->SetFillColor(kPink+9);
  MET_TT[11]->SetStats(0);
  MET_DY[11]->SetFillColor(kViolet+9);
  MET_DY[11]->SetStats(0);
  MET_W[11]->SetFillColor(kYellow-4);
  MET_W[11]->SetStats(0);
  MET_Z[11]->SetFillColor(kSpring+4);
  MET_Z[11]->SetStats(0);

  stack11->Add(MET_DY[11]);
  stack11->Add(MET_TT[11]);
  stack11->Add(MET_W[11]);
  stack11->Add(MET_Z[11]);
  

  aux2 = new TH1F(*MET_DY[11]);
  aux2->Sumw2();
  aux2->Add( MET_TT[11], 1);
  aux2->Add( MET_W[11], 1);
  aux2->Add( MET_Z[11], 1);

  lege = new TLegend(0.68,0.68,0.88,0.92);
  MET_Data[11]->Sumw2();
  lege->AddEntry(MET_W[11],"W + jets","f");
  lege->AddEntry(MET_Z[11],"Z(#nu#nu) + jets","f");
  lege->AddEntry(MET_TT[11],"t #bar{t} + jets","f");
  lege->AddEntry(MET_DY[11],"Z/#gamma^{*}(ll) + jets","f");
  lege->AddEntry(MET_Data[11],"Data","lep");

  stack11->Draw();
  MET_Data[11]->Sumw2();
  MET_Data[11]->SetMarkerStyle(20);
  MET_Data[11]->SetLineColor(1);
  MET_Data[11]->SetMarkerSize(1.5);
  MET_Data[11]->Draw("same");

  RatioPlotsV2( stack11, aux2, MET_Data[11], "NJETS MC  2 #mu BOX", "NJETS Data 2 #mu BOX", "METplots/Data_MC_NJETS_2BOX","NJETS", lege);

  delete aux2;
  delete lege;


  TCanvas* C1 = new TCanvas("C1", "C1", 1024, 1024);
  C1->cd();
  
  
  delete C0;
  delete C1;
  delete W;
  delete TT;
  delete Z;
  delete dy;
  delete data;
  delete stack1;
  delete stack2;
  delete stack3;
  delete stack4;
  delete stack6;
  delete stack7;
  
}
