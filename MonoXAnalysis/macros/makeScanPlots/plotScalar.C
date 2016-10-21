#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include <iostream>

int mmed(double mh, int code){
    if (code == 800) return ((int)(mh-80000000000))/10000; 
    if (code == 801) return ((int)(mh-80100000000))/10000; 
    if (code == 805) return ((int)(mh-80500000000))/10000; 
    if (code == 806) return ((int)(mh-80600000000))/10000; 
    return -1;
}

int mdm(double mh, int code){
    if (code == 800) return (mh-80000000000)  - ( ((Int_t)(mh-80000000000))/10000 )*10000;
    if (code == 801) return (mh-80100000000)  - ( ((Int_t)(mh-80100000000))/10000 )*10000;
    if (code == 805) return (mh-80500000000)  - ( ((Int_t)(mh-80500000000))/10000 )*10000;
    if (code == 806) return (mh-80600000000)  - ( ((Int_t)(mh-80600000000))/10000 )*10000;
    return -1;
}

int code(double mh){
    return (int)(mh/100000000);
}

void plotScalar(string inputDIR, string outputDIR, string coupling = "025", string energy = "13") {
  
  
  // Set the color palette
  bool useNicksPalette = true;
  int ncontours = 999;
  if (useNicksPalette) {
    
    TColor::InitializeColors();
    Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
    Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
    Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
    TColor::CreateGradientColorTable(9, stops, red, green, blue, ncontours);
  }
  else gStyle->SetPalette(70);
  
  gStyle->SetNumberContours(ncontours);
  
  // This is where all the plots are made
  TFile *file = new TFile((inputDIR+"/limits_805_"+string(coupling)+"_v1_"+string(energy)+"TeV.root").c_str());    
  TTree *tree = (TTree*)file->Get("limit");
  
  TFile* file2 = new TFile((inputDIR+"/scalar_out.root").c_str());
  
  TGraph2D* grexp = new TGraph2D();
  TGraph2D* grexp_up = new TGraph2D();
  TGraph2D* grexp_down = new TGraph2D();
  TGraph2D* grobs = new TGraph2D();
  TGraph2D* grobu = new TGraph2D();
  TGraph2D* grobd = new TGraph2D();
  
  double mh;
  double limit;
  float quantile;
  
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quantile);
  
  int expcounter = 0;
  int exp_up_counter = 0;
  int exp_down_counter = 0;
  int obscounter = 0;

  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    
    int c = code(mh);
    int medmass = mmed(mh, c);
    int dmmass = mdm(mh, c);
    
    if (medmass < 2* dmmass) continue;
    
    if (quantile == 0.5) {
      expcounter++;
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
    }
    
    
    if (quantile < 0.17 && quantile > 0.14 ) {
      exp_up_counter++;
      grexp_up->SetPoint(exp_up_counter, double(medmass), double(dmmass), limit);      
    }

    if (quantile < 0.85 && quantile > 0.83 ) {
      exp_down_counter++;
      grexp_down->SetPoint(exp_down_counter, double(medmass), double(dmmass), limit);      
    }
    
    if (quantile == -1) {
      
      obscounter++;
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      grobu->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobd->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
      
    }
  }
  tree->ResetBranchAddresses();
  
  TH2D* hexp = new TH2D("hexp", "", 250, 0, 5000, 75, 0, 1500);
  TH2D* hexp_up = new TH2D("hexp_up", "", 250, 0, 5000, 75, 0, 1500);
  TH2D* hexp_down = new TH2D("hexp_down", "", 250, 0, 5000, 75, 0, 1500);
  TH2D* hobs = new TH2D("hobs", "", 250, 0, 5000, 75, 0, 1500);
  TH2D* hobu = new TH2D("hobu", "", 250, 0, 5000, 75, 0, 1500);
  TH2D* hobd = new TH2D("hobd", "", 250, 0, 5000, 75, 0, 1500);
  
  for (int i = 1; i <= 250; i++) {
    for (int j = 1; j <= 75; j++) {
      hexp_up->SetBinContent(i, j, grexp_up->Interpolate(double(i)*20., double(j)*20.));
      hexp_down->SetBinContent(i, j, grexp_down->Interpolate(double(i)*20., double(j)*20.));
      hexp->SetBinContent(i, j, grexp->Interpolate(double(i)*20., double(j)*20.));
      hobs->SetBinContent(i, j, grobs->Interpolate(double(i)*20., double(j)*20.));
      hobu->SetBinContent(i, j, grobu->Interpolate(double(i)*20., double(j)*20.));
      hobd->SetBinContent(i, j, grobd->Interpolate(double(i)*20., double(j)*20.));
    }
  }

  TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
  TH2* hexp2_up = (TH2*)hexp_up->Clone("hexp2_up");
  TH2* hexp2_down = (TH2*)hexp_down->Clone("hexp2_down");
  TH2* hobs2 = (TH2*)hobs->Clone("hobs2");
  TH2* hobu2 = (TH2*)hobu->Clone("hobu2");
  TH2* hobd2 = (TH2*)hobd->Clone("hobd2");
  
  for (int i = 1; i <= hexp2_up->GetNbinsX(); i++) {
    for (int j = 1; j <= hexp2_up->GetNbinsY(); j++) {
      if (hexp2_up->GetBinContent(i, j) <= 0) hexp2_up->SetBinContent(i, j, 100.);
    }
  }
  
  for (int i = 1; i <= hexp2_down->GetNbinsX(); i++) {
        for (int j = 1; j <= hexp2_down->GetNbinsY(); j++) {
	  if (hexp2_down->GetBinContent(i, j) <= 0) hexp2_down->SetBinContent(i, j, 100.);
        }
    }
  
  for (int i = 1; i <= hexp2->GetNbinsX(); i++) {
    for (int j = 1; j <= hexp2->GetNbinsY(); j++) {
      if (hexp2->GetBinContent(i, j) <= 0) hexp2->SetBinContent(i, j, 100.);
    }
  }
  
  for (int i = 1; i <= hobs2->GetNbinsX(); i++) {
        for (int j = 1; j <= hobs2->GetNbinsY(); j++) {
	  if (hobs2->GetBinContent(i, j) <= 0) hobs2->SetBinContent(i, j, 100.);
        }
  }
  for (int i = 1; i <= hobu2->GetNbinsX(); i++) {
    for (int j = 1; j <= hobu2->GetNbinsY(); j++) {
      if (hobu2->GetBinContent(i, j) <= 0) hobu2->SetBinContent(i, j, 100.);
    }
    }
  for (int i = 1; i <= hobd2->GetNbinsX(); i++) {
    for (int j = 1; j <= hobd2->GetNbinsY(); j++) {
      if (hobd2->GetBinContent(i, j) <= 0) hobd2->SetBinContent(i, j, 100.);
    }
  }
  
  hexp2->SetContour(2);
  hexp2->SetContourLevel(1, 1);    
  hexp2_up->SetContour(2);
  hexp2_up->SetContourLevel(1, 1);
  hexp2_down->SetContour(2);
  hexp2_down->SetContourLevel(1, 1);
  hobs2->SetContour(2);
  hobs2->SetContourLevel(1, 1);
  hobu2->SetContour(2);
  hobu2->SetContourLevel(1, 1);
  hobd2->SetContour(2);
  hobd2->SetContourLevel(1, 1);
    
  // All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 800);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(20., 10., 200., 100.0, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetTitle("M_{med} [GeV]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.15);
  
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetLabelSize(0.035,"Z");
  
  frame->Draw();
  
  hexp2->SetLineStyle(2);
  hexp2->SetLineWidth(3);
  
  hexp2_up->SetLineStyle(2);
  hexp2_up->SetLineWidth(1);
  hexp2_down->SetLineStyle(2);
  hexp2_down->SetLineWidth(1);
  
  
  hobs2->SetLineWidth(3);
  hobu2->SetLineWidth(1);
  hobd2->SetLineWidth(1);
  hobu2->SetLineColor(kOrange);
  hobd2->SetLineColor(kOrange);
  
  hexp->SetMinimum(0.1);
  hexp->SetMaximum(5.);
  hobs->SetMinimum(0.1);
  hobs->SetMaximum(5.);
  hobu->SetMinimum(0.1);
  hobu->SetMaximum(5.);
  hobd->SetMinimum(0.1);
  hobd->SetMaximum(5.);
    
  hobs2->SetLineColor(2);
  hexp2->SetLineColor(kBlue);
  hobu2->SetLineColor(2);
  hobd2->SetLineColor(2);
  
  hobs->Draw("COLZ SAME");
  hexp2->Draw("CONT3 SAME");
  hexp2_up->SetLineColor(kBlue);
  hexp2_up->Draw("CONT3 SAME");
  hexp2_down->SetLineColor(kBlue);
  hexp2_down->Draw("CONT3 SAME");
  hexp2->Draw("CONT3 SAME");
  hobs2->SetLineColor(kRed);
  hobs2->Draw("CONT3 SAME");

  TLegend *leg = new TLegend(0.36,0.62,0.80,0.88,NULL,"brNDC");

  leg->AddEntry(hexp2,"Median Expected ("+TString(energy)+" TeV) 95% CL","L");
  leg->AddEntry(hexp2_up,"Expected ("+TString(energy)+" TeV) #pm 1#sigma_{experiment} ","L");
  leg->SetFillColor(0);
  leg->Draw("SAME");
  
  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (energy == "8"){
    tex->DrawLatex(0.68,0.96,"19.6 fb^{-1} (8 TeV)");
  } 
  if (energy == "13"){
    tex->DrawLatex(0.68,0.96,"7.6 fb^{-1} (13 TeV)");
  }
  if (energy == "8+13"){
    tex->DrawLatex(0.68,0.96,"7.6 fb^{-1} (13 TeV) + 19.6 fb^{-1} (8 TeV)");
  }
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(270);
  tex2->DrawLatex(0.965,0.93,"Observed    #sigma_{95% CL}/#sigma_{th}");
  
  TLatex* texCMS = new TLatex(0.22,0.96,"#bf{CMS} #it{Preliminary}");
  texCMS->SetNDC();
  texCMS->SetTextFont(42);
  texCMS->SetLineWidth(2);
  texCMS->SetTextSize(0.042); texCMS->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.36,0.90,"#bf{Scalar-vector, Dirac, g_{q} = 1, g_{DM} = 1}");
  }
  else
    tex->DrawLatex(0.36,0.90,"#bf{Scalar-vector, Dirac, g_{q} = 0.25, g_{DM} = 1}");
  
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();

  canvas->SaveAs((outputDIR+"/scan_scalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_scalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.png").c_str());

  TFile* outputFile = new TFile((outputDIR+"/fullLikelihood_scan_scalar.root").c_str(),"RECREATE");
  outputFile->cd();
  hexp->Write("scan_expected");
  hobs->Write("scan_observed");
  hexp2->Write("contour_expected");
  hobs2->Write("contour_observed");
  outputFile->Write();
}
