#include "CMS_lumi.h"

void makeHiggsInvisiblePlot(){

  setTDRStyle();

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.06);
  canvas->SetRightMargin(0.05);

  TH1F* limitValueExpected = new TH1F("limitValueExpected","",3,0,3);
  limitValueExpected->SetBinContent(1,1.13);
  limitValueExpected->SetBinContent(2,1.42);
  limitValueExpected->SetBinContent(3,0.86);

  TH1F* limitValueObserved = new TH1F("limitValueObserved","",3,0,3);
  limitValueObserved->SetBinContent(1,1.62);
  limitValueObserved->SetBinContent(2,1.02);
  limitValueObserved->SetBinContent(3,0.89);

				    
  TGraphAsymmErrors* limitErr_1s = new TGraphAsymmErrors();
  limitErr_1s->SetPoint(1,0.5,1.13);
  limitErr_1s->SetPoint(2,1.5,1.42);
  limitErr_1s->SetPoint(3,2.5,0.86);
  limitErr_1s->SetPointError(1,0.5,0.5,fabs(1.13-0.78),fabs(1.68-1.13));
  limitErr_1s->SetPointError(2,0.5,0.5,fabs(1.42-1.00),fabs(1.42-2.06));
  limitErr_1s->SetPointError(3,0.5,0.5,fabs(0.61-0.89),fabs(0.89-1.25));

  TGraphAsymmErrors* limitErr_2s = new TGraphAsymmErrors();
  limitErr_2s->SetPoint(1,0.5,1.13);
  limitErr_2s->SetPoint(2,1.5,1.42);
  limitErr_2s->SetPoint(3,2.5,0.86);
  limitErr_2s->SetPointError(1,0.5,0.5,fabs(1.13-0.57),fabs(2.42-1.13));
  limitErr_2s->SetPointError(2,0.5,0.5,fabs(1.42-0.74),fabs(1.42-2.85));
  limitErr_2s->SetPointError(3,0.5,0.5,fabs(0.45-0.89),fabs(0.89-1.75));
  
  limitValueExpected->GetXaxis()->SetBinLabel(1,"monojet");
  limitValueExpected->GetXaxis()->SetBinLabel(2,"mono-V");
  limitValueExpected->GetXaxis()->SetBinLabel(3,"combined");

  limitValueExpected->SetLineColor(kBlack);
  limitValueExpected->SetLineWidth(2);
  limitValueExpected->SetLineStyle(2);
  limitValueExpected->GetYaxis()->SetRangeUser(0.,4.5);
  limitValueExpected->GetYaxis()->SetLabelSize(0.035);
  limitValueExpected->GetXaxis()->SetLabelSize(0.05);
  limitValueExpected->GetYaxis()->SetTitle("95% C.L. upper limit on #mu = #sigma x BR/#sigma_{SM}");
  limitValueExpected->GetYaxis()->SetTitleSize(0.038);
  limitValueExpected->GetYaxis()->SetTitleOffset(1.25);

  limitValueExpected->Draw("hist");
  limitErr_1s->SetFillColor(kGreen);
  limitErr_1s->SetLineWidth(0);
  limitErr_2s->SetFillColor(kYellow);
  limitErr_2s->Draw("2same");
  limitErr_1s->Draw("2same");
  limitValueExpected->Draw("hist same");
  limitValueExpected->SetMarkerSize(1.2);
  limitValueExpected->SetMarkerStyle(24);
  limitValueExpected->SetMarkerColor(kBlack);
  limitValueExpected->Draw("P same");

  limitValueObserved->SetLineColor(kBlack);
  limitValueObserved->SetMarkerColor(kBlack);
  limitValueObserved->SetMarkerSize(1.2);
  limitValueObserved->SetMarkerStyle(20);
  limitValueObserved->SetLineWidth(2);
  limitValueObserved->Draw("hist same");
  limitValueObserved->Draw("P same");

  TLegend* leg = new TLegend(0.5,0.66,0.9,0.9);
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(limitValueObserved,"Oberved CL_{s}","PL");
  leg->AddEntry(limitValueExpected,"Median Expected CL_{s}","PL");
  leg->AddEntry(limitErr_1s,"Expected CL_{s} #pm 1#sigma","F");
  leg->AddEntry(limitErr_2s,"Expected CL_{s} #pm 2#sigma","F");
  leg->Draw("same");
  CMS_lumi(canvas,"2.30",false,2);

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs("higgsInv_brazilian.png","png");
  canvas->SaveAs("higgsInv_brazilian.pdf","pdf");
   
  

}
