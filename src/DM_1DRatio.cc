#include "DM_1DRatio.hh"

int RatioPlots(TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType"){
	
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  TLegend* leg;
  
  TH1F*  RATIO;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200,3500);
    h2->GetXaxis()->SetRangeUser(200,3500);
    RATIO->GetXaxis()->SetRangeUser(200,3500);

  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 2.5);
    h2->GetXaxis()->SetRangeUser(0.5, 2.5);
    RATIO->GetXaxis()->SetRangeUser(0.5, 2.5);

  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  std::cout << "=====================Dividing Histograms=====================" << std::endl;
  RATIO->Divide(h1, h2, 1, 1, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 2.05);

  h1->SetLineColor(2);
  h1->SetMarkerSize(1.);
  h1->SetMarkerColor(2);
  h1->SetMarkerStyle(20);
  h1->SetFillColor(2);
  
  h2->SetLineColor(4);
  h2->SetMarkerSize(1.);
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(20);
  h2->SetFillColor(4);
  
  h1->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");
  h2->SetStats(0);
  h2->SetXTitle( type );
  h1->SetXTitle( type );
  
  //RATIO->GetXaxis()->SetLabelFont(63); //font in pixels
  //RATIO->GetXaxis()->SetLabelSize(16); //in pixels
  //RATIO->GetYaxis()->SetLabelFont(63); //font in pixels
  //RATIO->GetYaxis()->SetLabelSize(16); //in pixels

  //RATIO->GetYaxis()->SetTitleOffset(0.24);
  //RATIO->GetXaxis()->SetTitleOffset(0.29);

  //RATIO->GetYaxis()->SetTitleSize(0.14);
  //RATIO->GetXaxis()->SetTitleSize(0.12);

  std::cout << "GET X Title Size:  " << RATIO->GetYaxis()->GetTitleSize() << std::endl;
   
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  h1->DrawCopy();
  h2->Draw("same");
  C->cd();
  
  leg = new TLegend(0.55, 0.7, 0.89, 0.9);//(xmin, ymin, xmax, ymax)
  leg->AddEntry(h1, label + " " + h1Name ,"f");
  leg->AddEntry(h2, label + " " + h2Name ,"f");
  //leg->SetHeader("");
  leg->SetTextSize(.022);
  leg->SetFillColor(0);
  leg->Draw();
  pad1->SetLogy();
  C->Update();
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  //RATIO->GetXaxis()->SetLabelSize(0.05);
  //RATIO->GetYaxis()->SetLabelSize(0.05);
  //RATIO->GetYaxis()->SetTitleOffset(0.35);
  //RATIO->GetXaxis()->SetTitleOffset(0.4);
  //RATIO->GetYaxis()->SetTitleSize(0.11);
  //RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetLabelSize(0.1);
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);

  /*
  pad2->SetLogy();
  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->SetLabelSize(0.2);
  RATIO->GetXaxis()->SetLabelSize(0.046);
  RATIO->GetXaxis()->SetTitleOffset(.9);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Ratio");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(1.3);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);
  */
 
  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  delete leg;
  delete C;
  delete RATIO;
  
  return 0;
  
};


int RatioPlotsV2(THStack* s, TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", TLegend* le = NULL){

  TCanvas* C = new TCanvas("C", "C      ", 400, 500);
  C->cd();

  TH1F*  RATIO;
  TString label;
  if(type == "MR"){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::MR_Bins, BaseDM::MR_BinArr);
    label = "M_{R}";
    h1->GetXaxis()->SetRangeUser(200.,3500.);
    h2->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetXaxis()->SetRangeUser(200.,3500.);
    RATIO->GetYaxis()->SetRangeUser(.0, 2.0);
    s->SetMaximum(100000.);
  }else if(type == "RSQ" ){
    RATIO = new TH1F("RATIO", fname + "_" + type , BaseDM::RSQ_Bins, BaseDM::RSQ_BinArr);
    label = "R^{2}";
    h1->GetXaxis()->SetRangeUser(0.5, 2.50);
    h2->GetXaxis()->SetRangeUser(0.5, 2.50);
    RATIO->GetXaxis()->SetRangeUser(0.5, 2.50);
    RATIO->GetYaxis()->SetRangeUser(.0, 2.0);
    s->SetMaximum(100000.);
  }else if(type == "MET"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 50, 0, 1000);
    label = "#slash{E}_{T}  GeV";
    s->SetMaximum(10000.);
  }else if(type == "NJETS"){
    RATIO = new TH1F("RATIO", fname + "_" + type , 9, 1, 10);
    label = "Jet Multiplicity";
    s->SetMaximum(100000.);
  }else{
    delete RATIO;
    delete C;
    std::cout << "Unknown Type, please use: MR or RSQ" << std::endl;
    return -1;
  }
  
  std::cout << "=====================Dividing Histograms=====================" << std::endl;
  RATIO->Divide(h2, h1, 1, 1, "");
  RATIO->GetYaxis()->SetRangeUser(.0, 3.05);
  h2->SetMarkerSize(.7);
  h2->SetStats(0);
  
  s->SetMinimum(1.);
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  s->Draw();
  if(type == "MET"){
    s->GetYaxis()->SetTitle("Events/20 GeV");
  }else{
    s->GetYaxis()->SetTitle("Events");
  }
  s->GetYaxis()->SetTitleOffset(1.25);
  gPad->Modified();
  h2->SetStats(0);
  h2->Draw("same");
  C->cd();
  
  le->SetFillColor(0);
  le->SetBorderSize(0);
  le->SetTextSize(0.02);
  le->Draw();
  pad1->SetLogy();
  C->Update();
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.25,0.87,"CMS Preliminary");
  t->DrawLatex(0.25,0.83,"#sqrt{s} = 8 TeV");
  t->DrawLatex(0.25,0.77,"#int L dt = 5 fb^{-1}");
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
  pad2->SetTopMargin(0.008);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  RATIO->SetLineColor(4);
  RATIO->SetStats(0);
  RATIO->SetTitle("");
  RATIO->GetXaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetLabelSize(0.1); 
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetXaxis()->SetTitleOffset(0.88);
  RATIO->GetYaxis()->SetTitleSize(0.11);
  RATIO->GetXaxis()->SetTitleSize(0.11);
  RATIO->SetXTitle( label );
  RATIO->SetYTitle("Data/MC");
  RATIO->SetLineColor(4);
  RATIO->SetMarkerSize(.7);
  RATIO->SetMarkerColor(4);
  RATIO->SetMarkerStyle(20);
  RATIO->SetFillColor(4);
  RATIO->Draw();
  C->cd();
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");

  delete C;
  delete RATIO;

  return 0;

};
