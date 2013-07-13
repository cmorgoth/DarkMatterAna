#include "DM_RatioPlots.hh"
#include "DM_1DRatio.hh"
#include "DM_2DRatio.hh"
#include "DM_2D_MR_RSQ_Dist.hh"

void Create2DPlots(){

  TCanvas* C0 = new TCanvas("C0", "C0", 1024, 1024);
  C0->cd();  
  TLegend* leg;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// W + Jets //////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  WJetsHTBins* W = new WJetsHTBins( 2 );
  
      
  ///////////////////////////////////////////////////////
  /////////////////WJETS 2D HISTOS///////////////////////
  ///////////////////////////////////////////////////////
  std::vector<TH2F*> Razor2D_W = W->Plot_2DRazor();
  
  int ctr = 0;
  TString Razor2Dpicname;
  TString box;
  TLatex *t;
  for( std::vector<TH2F*>::iterator itW = Razor2D_W.begin(); itW != Razor2D_W.end(); ++itW ){
    C0->SetRightMargin(0.125);
    (*itW)->SetStats(0);
    (*itW)->SetTitle("");
    (*itW)->SetXTitle("M_{R}");
    (*itW)->SetYTitle("R^{2}");
    (*itW)->Draw("colz");
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.73,0.87,"CMS Simulation");
    t->DrawLatex(0.73,0.83,"#sqrt{s} = 8 TeV");
    t->DrawLatex(0.73,0.79,"W(l#nu) + Jets");
    box = TString(Form("%d #mu BOX",ctr));
    t->DrawLatex(0.73,0.75, box);
    //(*itW)->Draw("colz");
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/W_Razor2D_box%d.pdf",ctr));
    C0->SaveAs(Razor2Dpicname);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/W_Razor2D_box%d.png",ctr));
    C0->SaveAs(Razor2Dpicname);
    ctr++;
    delete t;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Z(NuNu) + Jets ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ZJetsNuNu* Z = new ZJetsNuNu( 2 );
  
  std::vector<TH2F*> Razor2D_Z = Z->Plot_2DRazor();

  ctr = 0;
  for( std::vector<TH2F*>::iterator itZ = Razor2D_Z.begin(); itZ != Razor2D_Z.end(); ++itZ ){
    C0->SetRightMargin(0.125);
    (*itZ)->SetStats(0);
    (*itZ)->SetTitle("");
    (*itZ)->SetXTitle("M_{R}");
    (*itZ)->SetYTitle("R^{2}");
    (*itZ)->Draw("colz");
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.73,0.87,"CMS Simulation");
    t->DrawLatex(0.73,0.83,"#sqrt{s} = 8 TeV");
    t->DrawLatex(0.73,0.79,"Z(#nu#bar{#nu}) + Jets");
    box = TString(Form("%d #mu BOX",ctr));
    t->DrawLatex(0.73,0.75, box);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/Z_Razor2D_box%d.pdf",ctr));
    C0->SaveAs(Razor2Dpicname);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/Z_Razor2D_box%d.png",ctr));
    C0->SaveAs(Razor2Dpicname);
    ctr++;
    delete t;
  }

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// ZGamma(ll) + Jets /////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  DY* dy = new DY(2);
  
  std::vector<TH2F*> Razor2D_DY = dy->Plot_2DRazor();
  
  ctr = 0;
  for( std::vector<TH2F*>::iterator itDY = Razor2D_DY.begin(); itDY != Razor2D_DY.end(); ++itDY ){
    C0->SetRightMargin(0.125);
    (*itDY)->SetStats(0);
    (*itDY)->SetTitle("");
    (*itDY)->SetXTitle("M_{R}");
    (*itDY)->SetYTitle("R^{2}");
    (*itDY)->Draw("colz");
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.73,0.87,"CMS Simulation");
    t->DrawLatex(0.73,0.83,"#sqrt{s} = 8 TeV");
    t->DrawLatex(0.73,0.79,"Z/#gamma^{*}(ll) + Jets");
    box = TString(Form("%d #mu BOX",ctr));
    t->DrawLatex(0.73,0.75, box);

    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/DY_Razor2D_box%d.pdf",ctr));
    C0->SaveAs(Razor2Dpicname);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/DY_Razor2D_box%d.png",ctr));
    C0->SaveAs(Razor2Dpicname);
    ctr++;
  }
  
  
  ///////////////////////////////////////////////
  //////////////////////////////////////////////
  /////////////SMS//////////////////////////////
  /////////////////////////////////////////////
  /////////////////////////////////////////////

  /*
  T2cc* sms = new T2cc("/afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src/VecbosApp/V07_5_3_X/T2ccNoFilter/out/T2ccNoFilter_Total.root", 2);
  
  ///////////////////////////////////////////                                                                    
  //////////////////2D///////////////////////                                         
  //////////////////////////////////////////                                                                                    

  TH2F*  MR_R2_0BOX_SMS = new TH2F( );
  sms->PlotRSQ_vs_MR_0Box( MR_R2_0BOX_SMS );
  MR_R2_0BOX_SMS->Draw("colz");
  C0->SaveAs("MR_RSQ_PLOTS/MR_RSQ_SMS_0box.pdf");

  TH2F*  MR_R2_1BOX_SMS = new TH2F( );
  sms->PlotRSQ_vs_MR_1Box( MR_R2_1BOX_SMS );
  MR_R2_1BOX_SMS->Draw("colz");
  C0->SaveAs("MR_RSQ_PLOTS/MR_RSQ_SMS_1box.pdf");
  
  TH2F*  MR_R2_2BOX_SMS = new TH2F( );
  sms->PlotRSQ_vs_MR_2Box( MR_R2_2BOX_SMS );
  MR_R2_2BOX_SMS->Draw("colz");
  C0->SaveAs("MR_RSQ_PLOTS/MR_RSQ_SMS_2box.pdf");
  */
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// TTbar + Jets //////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  const char* ttjets_file = "root://eoscms//eos/cms/store/group/phys_susy/razor/Cristian_DM/PFJets80/TTJets80_Tot.root";

  TTJets* TT = new TTJets(2);
  
  std::vector<TH2F*> Razor2D_TT = TT->Plot_2DRazor();

  ctr = 0;
  for( std::vector<TH2F*>::iterator itTT = Razor2D_TT.begin(); itTT != Razor2D_TT.end(); ++itTT ){
    C0->SetRightMargin(0.125);
    (*itTT)->SetStats(0);
    (*itTT)->SetTitle("");
    (*itTT)->SetXTitle("M_{R}");
    (*itTT)->SetYTitle("R^{2}");
    (*itTT)->Draw("colz");
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.73,0.87,"CMS Simulation");
    t->DrawLatex(0.73,0.83,"#sqrt{s} = 8 TeV");
    t->DrawLatex(0.73,0.79,"t#bar{t} + Jets");
    box = TString(Form("%d #mu BOX",ctr));
    t->DrawLatex(0.73,0.75, box);
    
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/TT_Razor2D_box%d.pdf",ctr));
    C0->SaveAs(Razor2Dpicname);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/TT_Razor2D_box%d.png",ctr));
    C0->SaveAs(Razor2Dpicname);
    ctr++;
  }
  
  ////////////////////////////////////////////////////////////////
  //////////////////// DATA//////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  const char* data_file = "/afs/cern.ch/work/c/cpena/DarkMatter/CMSSW_5_2_3/src/VecbosApp/53X/HTMHT_Run2012A_ILV/out/HTMHT_ILV_Run2012AB.root";

  Data* data = new Data(data_file, 2);
        
  std::vector<TH2F*> Razor2D_data = data->Plot_2DRazor();

  ctr = 0;
  for( std::vector<TH2F*>::iterator itdata = Razor2D_data.begin(); itdata != Razor2D_data.end(); ++itdata ){
    C0->SetRightMargin(0.125);
    (*itdata)->SetStats(0);
    (*itdata)->SetTitle("");
    (*itdata)->SetXTitle("M_{R}");
    (*itdata)->SetYTitle("R^{2}");
    (*itdata)->Draw("colz");
    t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.03);
    t->DrawLatex(0.73,0.87,"CMS Simulation");
    t->DrawLatex(0.73,0.83,"#sqrt{s} = 8 TeV");
    t->DrawLatex(0.73,0.79,"Data + Jets");
    box = TString(Form("%d #mu BOX",ctr));
    t->DrawLatex(0.73,0.75, box);


    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/data_Razor2D_box%d.pdf",ctr));
    C0->SaveAs(Razor2Dpicname);
    Razor2Dpicname = TString(Form("MR_RSQ_PLOTS/data_Razor2D_box%d.png",ctr));
    C0->SaveAs(Razor2Dpicname);
    ctr++;
  }

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// FINALIZING //////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  delete C0;
  delete W;
  delete TT;
  delete Z;
  delete dy;
  delete data;
  
  
  
  
};
