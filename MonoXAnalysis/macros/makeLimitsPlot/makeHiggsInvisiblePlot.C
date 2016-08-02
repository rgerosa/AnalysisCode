#include "../CMS_lumi.h"

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
  limitValueExpected->SetBinContent(1,0.8477);
  limitValueExpected->SetBinContent(2,0.7207);
  limitValueExpected->SetBinContent(3,0.5645);

  TH1F* limitValueObserved = new TH1F("limitValueObserved","",3,0,3);
  limitValueObserved->SetBinContent(1,0.4826);
  limitValueObserved->SetBinContent(2,1.1692);
  limitValueObserved->SetBinContent(3,0.4371);
				    
  TGraphAsymmErrors* limitErr_1s = new TGraphAsymmErrors();
  limitErr_1s->SetPoint(1,0.5,0.8477);
  limitErr_1s->SetPoint(2,1.5,0.7207);
  limitErr_1s->SetPoint(3,2.5,0.5645);
  limitErr_1s->SetPointError(1,0.5,0.5,fabs(0.8477-0.5842),fabs(0.8477-1.2666));
  limitErr_1s->SetPointError(2,0.5,0.5,fabs(0.7207-0.5144),fabs(0.7207-1.0223));
  limitErr_1s->SetPointError(3,0.5,0.5,fabs(0.5645-0.3994),fabs(0.5645-0.8052));

  TGraphAsymmErrors* limitErr_2s = new TGraphAsymmErrors();
  limitErr_2s->SetPoint(1,0.5,0.8477);
  limitErr_2s->SetPoint(2,1.5,0.7207);
  limitErr_2s->SetPoint(3,2.5,0.5645);
  limitErr_2s->SetPointError(1,0.5,0.5,fabs(0.8477-0.4188),fabs(0.8477-1.8555));
  limitErr_2s->SetPointError(2,0.5,0.5,fabs(0.7207-0.3843),fabs(0.7207-1.3923));
  limitErr_2s->SetPointError(3,0.5,0.5,fabs(0.5645-0.2988),fabs(0.5645-1.1078));
  
  limitValueExpected->GetXaxis()->SetBinLabel(1,"monojet");
  limitValueExpected->GetXaxis()->SetBinLabel(2,"mono-V");
  limitValueExpected->GetXaxis()->SetBinLabel(3,"combined");

  limitValueExpected->SetLineColor(kBlack);
  limitValueExpected->SetLineWidth(2);
  limitValueExpected->SetLineStyle(2);
  limitValueExpected->GetYaxis()->SetRangeUser(0.,2.5);
  limitValueExpected->GetYaxis()->SetLabelSize(0.035);
  limitValueExpected->GetXaxis()->SetLabelSize(0.05);
  limitValueExpected->GetYaxis()->SetTitle("95% C.L. upper limit on #mu = #sigma x BR/#sigma_{SM}");
  limitValueExpected->GetYaxis()->SetTitleSize(0.038);
  limitValueExpected->GetYaxis()->SetTitleOffset(1.25);

  limitValueExpected->Draw("hist");
  limitErr_1s->SetFillColor(kGreen);
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
  CMS_lumi(canvas,"12.9",false);

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs("higgsInv_brazilian.png","png");
  canvas->SaveAs("higgsInv_brazilian.pdf","pdf");
   
  

}
