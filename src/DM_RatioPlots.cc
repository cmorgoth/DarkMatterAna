#include "DM_RatioPlots.hh"
#include "DM_1DRatio.hh"
#include "DM_2DRatio.hh"

void CreateRatioPlots(){

  TCanvas* C0 = new TCanvas("C0", "C0", 1024, 1024);
  C0->cd();  
  TLegend* leg;
  
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
  std::vector<TH1F*> Wjets = W->Plot_1DRazor();
  //////////////////////////////////////////////////////
  ///////////////WJETS 1D HISTOS////////////////////////
  //////////////////////////////////////////////////////
  
  std::cout << "debug 0" << std::endl;
  TH1F* RSQ_W00 = new TH1F( *Wjets[1] );
  
  RSQ_W00->Scale(1./RSQ_W00->Integral());
  std::cout << "debug 1" << std::endl;
  TH1F* RSQ_W11 = new TH1F( *Wjets[3] );
  
  RSQ_W11->Scale(1./RSQ_W11->Integral());
  std::cout << "debug 2" << std::endl;
  TH1F* MR_W00 = new TH1F( *Wjets[0] );
  
  MR_W00->Scale(1./MR_W00->Integral());
  std::cout << "debug 3" << std::endl;
  TH1F* MR_W11 = new TH1F( *Wjets[2] );
  
  MR_W11->Scale(1./MR_W11->Integral());

  ////////
  ////////
  std::cout << "debug 4" << std::endl;
  RatioPlots( MR_W00, MR_W11, "W + Jets (l#nu) 0 #mu BOX", "W + Jets (l#nu) 1 #mu BOX", "RatioPlots/WJetsMR_TauVeto", "MR");
  RatioPlots( RSQ_W00, RSQ_W11, "W + Jets (l#nu) 0 #mu BOX", "W + Jets (l#nu) 1 #mu BOX", "RatioPlots/WJetsRSQ_TauVeto", "RSQ");
  
  std::cout << "debug 5" << std::endl;  
  ///////////////////////////////////////////////////////
  /////////////////WJETS 2D HISTOS///////////////////////
  ///////////////////////////////////////////////////////
  TH2F*  MR_R2_1BOX = new TH2F( W->PlotRSQ_vs_MR_1Box() );
  
  std::cout << "debug 6" << std::endl;
  TH2F*  MR_RSQ_0BOX = new TH2F( W->PlotRSQ_vs_MR_0Box() );
  
  std::cout << "debug 7" << std::endl;
  TH2F* wjRatio = (TH2F*)RatioPlotsTwoD( MR_RSQ_0BOX, MR_R2_1BOX, "RatioPlots/Wjets2D", "Wjets2dRatio");

  std::cout << "debug 8" << std::endl;	
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Z(NuNu) + Jets ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

  

  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  Z->SetBtagCut(bL,bM,bT);
  ///////////////////////////////////////////
  ///////////////ZJETS(NuNu) 1D HISTOS///////
  ///////////////////////////////////////////
  std::vector<TH1F*> Zjets = Z->Plot_1DRazor();

  TH1F* MR_Z11 = new TH1F( *Zjets[2] );                                                                            
  MR_Z11->Scale( 1./MR_Z11->Integral() );                                                                                 
                                                                                                                          
  TH1F* MR_Z00 = new TH1F( *Zjets[0] );                                                                            
  MR_Z00->Scale( 1./MR_Z00->Integral() );                                                                                
                                                                                                                          
  TH1F* RSQ_Z11 = new TH1F( *Zjets[1] );                                                                          
  RSQ_Z11->Scale( 1./RSQ_Z11->Integral() );                                                                     
                                                                                                                          
  TH1F* RSQ_Z00 = new TH1F( *Zjets[3] );                                                                          
  RSQ_Z00->Scale( 1./RSQ_Z00->Integral() ); 
  
  ///////////////////////////////////////////
  //////////////////2D///////////////////////
  //////////////////////////////////////////
  
  TH2F*  MR_RSQ_0BOX_Z = new TH2F( Z->PlotRSQ_vs_MR_0Box() );
  MR_RSQ_0BOX_Z->Sumw2();
  
  TH2F*  MR_RSQ_1BOX_Z = new TH2F( Z->PlotRSQ_vs_MR_1Box() );
  MR_RSQ_1BOX_Z->Sumw2();
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// ZGamma(ll) + Jets /////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  DY* dy = new DY(2);
  dy->SetBtagCut(bL,bM,bT);
  ////////////////////////////////////////////////
  /////////////Zgamma*(ll)////////////////////////
  ////////////////////////////////////////////////
  
  std::vector<TH1F*> dy_jets = dy->Plot_1DRazor();
  
  TH1F* MR_DY22 = new TH1F( *dy_jets[4] );                                                                          
  MR_DY22->Scale( 1./MR_DY22->Integral() );                                                                              
                                                                                                                          
  TH1F* MR_DY11 = new TH1F( *dy_jets[2] );                                                                          
  MR_DY11->Scale( 1./MR_DY11->Integral() );                                                                               
                                                                                                                          
  TH1F* MR_DY00 = new TH1F( *dy_jets[0] );                                                                          
  MR_DY00->Scale( 1./MR_DY00->Integral() );                                                                               
                                                                                                                          
  TH1F* RSQ_DY22 = new TH1F( *dy_jets[5] );                                                                        
  RSQ_DY22->Scale( 1./RSQ_DY22->Integral() );                                                                             
                                                                                                                          
  TH1F* RSQ_DY11 = new TH1F( *dy_jets[3] );                                                                        
  RSQ_DY11->Scale( 1./RSQ_DY11->Integral() );                                                                             
                                                                                                                          
  TH1F* RSQ_DY00 = new TH1F( *dy_jets[1] );                                                                        
  RSQ_DY00->Scale( 1./RSQ_DY00->Integral() );
  
  ///////////////////////////////////////////
  //////////////////2D///////////////////////
  //////////////////////////////////////////
  
  TH2F*  MR_RSQ_0BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_0Box() );
  MR_RSQ_0BOX_DY->Sumw2();
  TH2F*  MR_RSQ_1BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_1Box() );
  MR_RSQ_1BOX_DY->Sumw2();
  TH2F*  MR_RSQ_2BOX_DY = new TH2F( dy->PlotRSQ_vs_MR_2Box() );
  MR_RSQ_2BOX_DY->Sumw2();
  
  std::cout << "DY0mu: " << MR_RSQ_0BOX_DY->Integral() << std::endl;
  std::cout << "DY1mu: " << MR_RSQ_1BOX_DY->Integral() << std::endl;
  std::cout << "DY2mu: " << MR_RSQ_2BOX_DY->Integral() << std::endl;
  
  //////////////////////////////////////
  ///////////DY 1mu/2mu BOX////////////
  /////////////////////////////////////
  
  RatioPlots( MR_DY11, MR_DY22, "Z/#gamma^{*}(ll) + Jets 1 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/DYJets_MR_1to2box", "MR");
  RatioPlots( RSQ_DY11, RSQ_DY22, "Z/#gamma^{*}(ll) + Jets 1 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/DYJets_RSQ_1to2box", "RSQ");
  
  TH2F* DY_1to2 = (TH2F*)RatioPlotsTwoD( MR_RSQ_1BOX_DY, MR_RSQ_2BOX_DY, "RatioPlots/DY2dRatio_1to2box", "DY2dRatio_1to2box");
  
  //////////////////////////////////////
  ///////////DY 0mu/2mu BOX////////////
  /////////////////////////////////////
  
  RatioPlots( MR_DY00, MR_DY22, "Z/#gamma^{*}(ll) + Jets 0 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/DYJets_MR_0to2box", "MR");
  RatioPlots( RSQ_DY00, RSQ_DY22, "Z/#gamma^{*}(ll) + Jets 0 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/DYJets_RSQ_0to2box", "RSQ");
  
  TH2F* DY_0to2 = (TH2F*)RatioPlotsTwoD( MR_RSQ_0BOX_DY, MR_RSQ_2BOX_DY, "RatioPlots/DY2dRatio_0to2box", "DY2dRatio_0to2box");
  
  //////////////////////////////////////
  ///////////DY-Znunu 0mu/2mu BOX///////
  /////////////////////////////////////
  
  RatioPlots( MR_Z00, MR_DY22, "Z(#nu #bar{#nu}) + Jets 0 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/Znunu_DYJets_MR_0to2box", "MR");
  RatioPlots( RSQ_Z00, RSQ_DY22, "Z(#nu #bar{#nu}) + Jets 0 #mu BOX", "Z/#gamma^{*}(ll) + Jets 2 #mu BOX", "RatioPlots/Znunu_DYJets_RSQ_0to2box", "RSQ");
  
  TH2F* Z_DY_0to2 = (TH2F*)RatioPlotsTwoD( MR_RSQ_0BOX_Z, MR_RSQ_2BOX_DY, "RatioPlots/Znunu_DYJets_2dRatio_1to2box", "Znunu_DYJets_2dRatio_0to2box");
  
  std::cout << "DY 0mu/2mu --> " << DY_0to2->Integral() << std::endl;
  std::cout << "Z DY 0mu/2mu --> " << Z_DY_0to2->Integral() << std::endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// TTbar + Jets //////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTJets* TT = new TTJets(2);
  TT->SetBtagCut(bL,bM,bT);
  /////////////////////////////////////////////////
  /////////////////TTbar + Jets 1D/////////////////
  ////////////////////////////////////////////////

  std::vector<TH1F*> TTjets = TT->Plot_1DRazor();
  
  TH1F* MR_22_TT = new TH1F( *TTjets[4] );                                                                         
  MR_22_TT->Scale( 1./MR_22_TT->Integral() );                                                                             
                                                                                                                          
  TH1F* MR_11_TT = new TH1F( *TTjets[2] );                                                                         
  MR_11_TT->Scale( 1./MR_11_TT->Integral() );                                                                             
                                                                                                                          
  TH1F* MR_00_TT = new TH1F( *TTjets[0] );                                                                         
  MR_00_TT->Scale( 1./MR_00_TT->Integral() );                                                                             
                                                                                                                          
  TH1F* RSQ_22_TT = new TH1F( *TTjets[5] );                                                                       
  RSQ_22_TT->Scale( 1./RSQ_22_TT->Integral() );                                                                           
                                                                                                                          
  TH1F* RSQ_11_TT = new TH1F( *TTjets[3] );                                                                       
  RSQ_11_TT->Scale( 1./RSQ_11_TT->Integral() );                                                                           
                                                                                                                          
  TH1F* RSQ_00_TT = new TH1F( *TTjets[1] );                                                                       
  RSQ_00_TT->Scale( 1./RSQ_00_TT->Integral() );
  
  
  ////////////////////////////////////////
  //////////////////2D TT+jets////////////
  ////////////////////////////////////////
  
  TH2F*  MR_RSQ_1BOX_TT = new TH2F( TT->PlotRSQ_vs_MR_1Box( ) );
  MR_RSQ_1BOX_TT->Sumw2();
  TH2F*  MR_RSQ_2BOX_TT = new TH2F( TT->PlotRSQ_vs_MR_2Box( ) );
  MR_RSQ_2BOX_TT->Sumw2();
  TH2F*  MR_RSQ_0BOX_TT = new TH2F( TT->PlotRSQ_vs_MR_0Box( ) ); 
  MR_RSQ_0BOX_TT->Sumw2();
  
  ////////////////////////////////////////////////////////////////
  //////////////////// DATA//////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  std::cout << "========= debug 0 =======" << std::endl;

  const char* data_file = "/media/data/cmorgoth/Data/DMData/HTMHT_ILV_Run2012AB.root";

  Data* data = new Data(data_file, 2);
  data->SetBtagCut(bL,bM,bT);
  
  TH1F* MR_22_data = new TH1F( data->PlotMR_2Box() );
  MR_22_data->Sumw2();
  MR_22_data->Scale ( 1./MR_22_data->Integral() );
  
  TH1F* MR_11_data = new TH1F( data->PlotMR_1Box() );
  MR_11_data->Sumw2();
  MR_11_data->Scale ( 1./MR_11_data->Integral() );
  
  TH1F* MR_00_data = new TH1F( data->PlotMR_0Box() );
  MR_00_data->Sumw2();
  MR_00_data->Scale ( 1./MR_00_data->Integral() );
  
  TH1F* RSQ_22_data = new TH1F( data->PlotRSQ_2Box() );
  RSQ_22_data->Sumw2();
  RSQ_22_data->Scale ( 1./RSQ_22_data->Integral() );  
  
  TH1F* RSQ_11_data = new TH1F( data->PlotRSQ_1Box() );
  RSQ_11_data->Sumw2();
  RSQ_11_data->Scale ( 1./RSQ_11_data->Integral() );  
  
  TH1F* RSQ_00_data = new TH1F( data->PlotRSQ_0Box() );
  RSQ_00_data->Sumw2();
  RSQ_00_data->Scale ( 1./RSQ_00_data->Integral() );
  
  ///////////////////////////////////////////
  //////////////////2D///////////////////////
  //////////////////////////////////////////
  
  TH2F*  data0 = new TH2F( data->PlotRSQ_vs_MR_0Box( ) );
    
  TH2F*  data1 = new TH2F( data->PlotRSQ_vs_MR_1Box( ) );
    
  TH2F*  data2 = new TH2F( data->PlotRSQ_vs_MR_2Box( ) );
  
  std::cout << "========= debug 1 =======" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////// Predicting BKG///////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////


  TH2F* data_0mu = new TH2F(*data0);
  TH2F* data_1mu = new TH2F(*data1);
  TH2F* data_2mu = new TH2F(*data2);
  
  TH2F* tt_0mu = new TH2F(*MR_RSQ_0BOX_TT);
  TH2F* tt_1mu = new TH2F(*MR_RSQ_1BOX_TT);
  TH2F* tt_2mu = new TH2F(*MR_RSQ_2BOX_TT);

  TH2F* dy_0mu = new TH2F(*MR_RSQ_0BOX_DY);
  TH2F* dy_1mu = new TH2F(*MR_RSQ_1BOX_DY);
  TH2F* dy_2mu = new TH2F(*MR_RSQ_2BOX_DY);

  TH2F* z_0mu = new TH2F(*MR_RSQ_0BOX_Z);
  TH2F* z_1mu = new TH2F(*MR_RSQ_1BOX_Z);

  TH2F* w_0mu = new TH2F(*MR_RSQ_0BOX);
  TH2F* w_1mu = new TH2F(*MR_RSQ_0BOX);

  std::cout << "========= debug 2 =======" << std::endl;
  
  TH2F* Nhat_DY_1mu = new TH2F( *data2 );//This will be the prediction for the Drell-Yan in the 1-mu box
  
  Nhat_DY_1mu->Add(MR_RSQ_2BOX_TT, -1.0);//Subtracting ttbar background from MC

  std::cout << "======= Prediction Drell-Yan 2-mu:: " << Nhat_DY_1mu->Integral() << std::endl;
  Nhat_DY_1mu->Multiply( DY_1to2 );
  std::cout << "======= Prediction Drell-Yan 1-mu: " << Nhat_DY_1mu->Integral() << std::endl;
      
  double sum_dy = .0;
  double sum_dy_err = .0;
  for(int i = 1; i <= Nhat_DY_1mu->GetNbinsX(); i++){
    for(int j = 1; j <= Nhat_DY_1mu->GetNbinsY(); j++){
      sum_dy_err += Nhat_DY_1mu->GetBinError(i,j)*Nhat_DY_1mu->GetBinError(i,j);
      sum_dy += Nhat_DY_1mu->GetBinContent(i,j);
    }
  }
  
  double ERROR = 0.0;
  double Int_Gral = Nhat_DY_1mu->IntegralAndError(1, 6, 1, 6, ERROR);
  
  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;
  
  TH2F* data_wj1mu = new TH2F( *data1 );//This will be the estimate for Wjets in the 1-mu box 
  std::cout << "======= Total Data 1 muon Box: " << data_wj1mu->Integral() << std::endl;
  
  data_wj1mu->Add(Nhat_DY_1mu, -1.0);//Subtracting the Drell-Yan estimation from the 1-mu
  data_wj1mu->Add(MR_RSQ_1BOX_TT, -1.0);//Subtracting the ttbar from MC
  
  std::cout << "====== Bkg prediction for W+jets in the 1-mu box: " << data_wj1mu->Integral() << std::endl;

  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= data_wj1mu->GetNbinsX(); i++){
    for(int j = 1; j <= data_wj1mu->GetNbinsY(); j++){
      sum_dy_err += data_wj1mu->GetBinError(i,j)*data_wj1mu->GetBinError(i,j);
      sum_dy += data_wj1mu->GetBinContent(i,j);
    }
  }
  
  ERROR = 0.0;
  Int_Gral = data_wj1mu->IntegralAndError(1, 6, 1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;


  std::cout << "========= TOTAL Bkg prediction 1 muon box: " << data_wj1mu->Integral() + Nhat_DY_1mu->Integral() + MR_RSQ_1BOX_TT->Integral()  << std::endl;
  
  TH2F* bkg_wj_0mu = new TH2F( *data_wj1mu );//Background Estimation for Wjets in the 1-mu box 
  bkg_wj_0mu->Multiply( wjRatio );//wjRatio is W + jest 0mu/1mu
 
  std::cout << "========= Bkg Prediction for W+jets in the 0 muon box: " << bkg_wj_0mu->Integral() << std::endl;

  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= bkg_wj_0mu->GetNbinsX(); i++){
    for(int j = 1; j <= bkg_wj_0mu->GetNbinsY(); j++){
      sum_dy_err += bkg_wj_0mu->GetBinError(i,j)*bkg_wj_0mu->GetBinError(i,j);
      sum_dy += bkg_wj_0mu->GetBinContent(i,j);
    }
  }
  
  ERROR = 0.0;
  Int_Gral = bkg_wj_0mu->IntegralAndError(1, 6, 1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;
  
  TH2F* bkg_ZDY = new TH2F( *data2 );//This will be Background Estimation for Drell-Yan in the 0-mu box
  bkg_ZDY->Add(MR_RSQ_2BOX_TT, -1.0);//Subtracting ttbar from MC

  TH2F* zero_2muon = new TH2F( *Z_DY_0to2 );//Drell-Yan Estimation 2-mu box
  bkg_ZDY->Multiply( Z_DY_0to2 );//Z(nunu) Estimation 0-mu box
  std::cout << "Z(nunu) Prediction 0Muon box : " << bkg_ZDY->Integral() << std::endl;

  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= bkg_ZDY->GetNbinsX(); i++){
    for(int j = 1; j <= bkg_ZDY->GetNbinsY(); j++){
      sum_dy_err += bkg_ZDY->GetBinError(i,j)*bkg_ZDY->GetBinError(i,j);
      sum_dy += bkg_ZDY->GetBinContent(i,j);
    }
  }

  ERROR = 0.0;
  Int_Gral = bkg_ZDY->IntegralAndError(1, 6, 1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;


  zero_2muon->Add( DY_0to2, 1.0);//Combined Ratio (Z(nunu)+Z(ll))0mu/Z(ll)2mu
  TH2F* bkg_ZDY_comb = new TH2F( *data2 );
  bkg_ZDY_comb->Add(MR_RSQ_2BOX_TT, -1.0);//Subtracting ttbar from MC
  DY_0to2->Multiply( bkg_ZDY_comb );
  std::cout << "====== Drell-Yan Prediction 0-mu box " << DY_0to2->Integral() << std::endl;
  
  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= DY_0to2->GetNbinsX(); i++){
    for(int j = 1; j <= DY_0to2->GetNbinsY(); j++){
      sum_dy_err += DY_0to2->GetBinError(i,j)*DY_0to2->GetBinError(i,j);
      sum_dy += DY_0to2->GetBinContent(i,j);
    }
  }
  ERROR = 0.0;
  Int_Gral = DY_0to2->IntegralAndError(1, 6, 1, 6, ERROR);
  
  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;

  
  zero_2muon->Multiply( bkg_ZDY_comb );
  
  std::cout << "====== Total(DY+Z(nunu)) Bkg contribution 0-mu box: " << zero_2muon->Integral() << std::endl;
  
  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= zero_2muon->GetNbinsX(); i++){
    for(int j = 1; j <= zero_2muon->GetNbinsY(); j++){
      sum_dy_err += zero_2muon->GetBinError(i,j)*zero_2muon->GetBinError(i,j);
      sum_dy += zero_2muon->GetBinContent(i,j);
    }
  }
  ERROR = 0.0;
  Int_Gral = zero_2muon->IntegralAndError(1, 6, 1, 6, ERROR);

  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl;
  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;


  TH2F* tot_bkg_0mu = new TH2F( *zero_2muon );
  TH2F* tt_0box = new TH2F( *MR_RSQ_0BOX_TT );
  tot_bkg_0mu->Add(bkg_wj_0mu,1.0);
  tot_bkg_0mu->Add(tt_0box,1.0);
  
  std::cout << "====== Total Bkg contribution 0 muon box: " << tot_bkg_0mu->Integral() << std::endl;
  ERROR = 0.0;
  Int_Gral = tot_bkg_0mu->IntegralAndError(1, 6, 1, 6, ERROR); 
  
  std::cout << "Integral: " << Int_Gral << " Error: " << ERROR << std::endl; 

  sum_dy = .0;
  sum_dy_err = .0;
  for(int i = 1; i <= tot_bkg_0mu->GetNbinsX(); i++){
    for(int j = 1; j <= tot_bkg_0mu->GetNbinsY(); j++){
      sum_dy_err += tot_bkg_0mu->GetBinError(i,j)*tot_bkg_0mu->GetBinError(i,j);
      sum_dy += tot_bkg_0mu->GetBinContent(i,j);
    }
  }

  std::cout << "--- Error: " << sqrt(sum_dy_err) << std::endl;
  std::cout << "--- sum from bins: " << sum_dy << std::endl;

  TFile* f2 = new TFile("2DHistos_TightBtag.root", "recreate");
  data0->Write("data");
  TH2F* MC = new TH2F( *MR_RSQ_0BOX );
  MC->Add(MR_RSQ_0BOX_Z);
  MC->Add(MR_RSQ_0BOX_DY);
  MC->Add(MR_RSQ_0BOX_TT);

  std::cout << "Data: " << data0->Integral() << " MC: " << MC->Integral() << std::endl;

  TH1F* bkg_projection_MR =  (TH1F*)tot_bkg_0mu->ProjectionX("M_{R}", 0, -1, "eo"); 

  TH1F* bkg_projection_R2 =  (TH1F*)tot_bkg_0mu->ProjectionY("R^{2}", 0, -1, "eo");
  
  //RatioPlots( , bkg_projection_MR, "0 #mu BOX", "0 #mu BOX Pred", "RatioPlots/ZeroMu_Pre", "MR");
  //RatioPlots( , bkg_projection_MR, "0 #mu BOX", "0 #mu BOX Pred", "RatioPlots/ZeroMu_Pre", "RSQ");

  MC->Write("bkg");
  bkg_projection_MR->Write("bkg_MR");
  bkg_projection_R2->Write("bkg_R2");

  TH1F* data_proj_MR = (TH1F*)data_0mu->ProjectionX("MR_data_0mu", 0, -1, "eo");
  TH1F* data_proj_R2 = (TH1F*)data_0mu->ProjectionY("R2_data_0mu", 0, -1, "eo");
  
  RatioPlots( data_proj_MR, bkg_projection_MR, "0 #mu BOX", "0 #mu BOX Pred", "RatioPlots/ZeroMu_Pre_MR", "MR");
  RatioPlots( data_proj_R2, bkg_projection_R2, "0 #mu BOX", "0 #mu BOX Pred", "RatioPlots/ZeroMu_Pre_R2", "RSQ");

  data_proj_MR->Write("data_MR");
  data_proj_R2->Write("data_R2");
  
  data_0mu->Write("data_0mu");
  data_1mu->Write("data_1mu");
  data_2mu->Write("data_2mu");
  
  tt_0mu->Write("tt_0mu");
  tt_1mu->Write("tt_1mu");
  tt_2mu->Write("tt_2mu");

  dy_0mu->Write("dy_0mu");
  dy_1mu->Write("dy_1mu");
  dy_2mu->Write("dy_2mu");

  w_0mu->Write("w_0mu");
  w_1mu->Write("w_1mu");
  
  z_0mu->Write("z_0mu");
  z_1mu->Write("z_1mu");
  
  f2->Close();
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// FINALIZING //////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  delete C0;
  delete W;
  delete TT;
  delete Z;
  delete dy;
  
  
  
  
};
