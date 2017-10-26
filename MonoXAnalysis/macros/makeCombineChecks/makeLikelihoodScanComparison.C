#include "../CMS_lumi.h"

void makeLikelihoodScanComparison(string inputFileName_1, string inputFileName_2, string nuisanceName, string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFile_1 = TFile::Open(inputFileName_1.c_str(),"READ");
  TFile* inputFile_2 = TFile::Open(inputFileName_2.c_str(),"READ");
  
  TTree* limit_1 = (TTree*) inputFile_1->Get("limit");
  TTree* limit_2 = (TTree*) inputFile_2->Get("limit");

  TGraph* graph_1 = new TGraph();
  TGraph* graph_2 = new TGraph();

  double mh;
  float r, deltaNLL;
  limit_1->SetBranchAddress("mh",&mh);
  limit_1->SetBranchAddress(nuisanceName.c_str(),&r);
  limit_1->SetBranchAddress("deltaNLL",&deltaNLL);

  for(int entries = 0; entries < limit_1->GetEntries(); entries++){
    limit_1->GetEntry(entries);
    graph_1->SetPoint(entries,r,2*deltaNLL);
  }

  limit_2->SetBranchAddress("mh",&mh);
  limit_2->SetBranchAddress(nuisanceName.c_str(),&r);
  limit_2->SetBranchAddress("deltaNLL",&deltaNLL);

  for(int entries = 0; entries < limit_2->GetEntries(); entries++){
    limit_2->GetEntry(entries);
    graph_2->SetPoint(entries,r,2*deltaNLL);
  }

  graph_1->Sort();
  graph_2->Sort();

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();
  graph_1->GetXaxis()->SetTitle(nuisanceName.c_str());
  graph_1->GetXaxis()->SetTitleOffset(1.1);
  graph_1->GetXaxis()->SetTitleSize(graph_1->GetXaxis()->GetTitleSize()*0.9);
  graph_1->GetYaxis()->SetTitleSize(graph_1->GetYaxis()->GetTitleSize()*0.9);
  graph_1->GetYaxis()->SetTitle("-2*#Delta Log(L)");  

  graph_1->SetLineColor(kBlue);
  graph_2->SetLineColor(kRed);
  graph_1->SetMarkerColor(kBlue);
  graph_2->SetMarkerColor(kRed);
  graph_1->SetLineWidth(2);
  graph_2->SetLineWidth(2);
  graph_2->SetLineStyle(7);
  graph_1->SetMarkerStyle(20);
  graph_2->SetMarkerStyle(24);
  graph_1->SetMarkerSize(1);
  graph_2->SetMarkerSize(1);

  graph_1->Draw("APL");
  graph_1->GetYaxis()->SetRangeUser(0,TMath::MaxElement(graph_1->GetN(),graph_1->GetY())*1.25);
  graph_2->Draw("PLsame");

  CMS_lumi(canvas,"35.9");

  TLegend leg (0.3,0.7,0.7,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(graph_1,"Split V-QCD / V-EW","PL");
  leg.AddEntry(graph_2,"Merge V-QCD / V-EW","PL");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/scan_"+nuisanceName+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/scan_"+nuisanceName+".pdf").c_str(),"pdf");
  
}
