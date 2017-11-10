#include "../CMS_lumi.h"

void makeHiggsInvisiblePlotNoVBF(bool blind = false, bool addPreliminary = false){


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

  TGraphErrors* limitValueExpected = new TGraphErrors();
  limitValueExpected->SetPoint(0,0.5,0.576);
  limitValueExpected->SetPoint(1,1.5,0.451);  
  limitValueExpected->SetPoint(2,2.5,0.404);
  limitValueExpected->SetPointError(0,0.5,0.);  
  limitValueExpected->SetPointError(1,0.5,0.);
  limitValueExpected->SetPointError(2,0.5,0.);

  TGraphErrors* limitValueObserved = new TGraphErrors();
  limitValueObserved->SetPoint(0,0.5,0.741);
  limitValueObserved->SetPoint(1,1.5,0.495);
  limitValueObserved->SetPoint(2,2.5,0.526);
  limitValueObserved->SetPointError(0,0.5,0.);  
  limitValueObserved->SetPointError(1,0.5,0.);
  limitValueObserved->SetPointError(2,0.5,0.);

  TGraphAsymmErrors* limitErr_1s = new TGraphAsymmErrors();
  limitErr_1s->SetPoint(0,0.5,0.576);
  limitErr_1s->SetPoint(1,1.5,0.4511);
  limitErr_1s->SetPoint(2,2.5,0.4042968);
  limitErr_1s->SetPointError(0,0.5,0.5,fabs(0.576-0.399),fabs(0.576-0.858));
  limitErr_1s->SetPointError(1,0.5,0.5,fabs(0.451-0.322),fabs(0.451-0.640));
  limitErr_1s->SetPointError(2,0.5,0.5,fabs(0.404-0.288),fabs(0.404-0.576));

  TGraphAsymmErrors* limitErr_2s = new TGraphAsymmErrors();
  limitErr_2s->SetPoint(0,0.5,0.576);
  limitErr_2s->SetPoint(1,1.5,0.4511);
  limitErr_2s->SetPoint(2,2.5,0.4042968);
  limitErr_2s->SetPointError(0,0.5,0.5,fabs(0.576-0.294),fabs(0.576-1.246));
  limitErr_2s->SetPointError(1,0.5,0.5,fabs(0.451-0.241),fabs(0.451-0.877));
  limitErr_2s->SetPointError(2,0.5,0.5,fabs(0.404-0.214),fabs(0.404-0.788));

  TH1* frame = (TH1*) canvas->DrawFrame(0.,0.,3,1.5);
  frame->SetBins(3,0,3);
  frame->GetXaxis()->SetBinLabel(1,"Monojet");
  frame->GetXaxis()->SetBinLabel(2,"Mono-V");
  frame->GetXaxis()->SetBinLabel(3,"Combined");
  frame->GetYaxis()->SetRangeUser(0.,1.8);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitle("95% CL upper limit on #sigma x B(H #rightarrow inv.)/#sigma_{SM}");
  frame->GetYaxis()->SetTitleSize(0.038);
  frame->GetYaxis()->SetTitleOffset(1.30);
  frame->Draw();

  limitValueExpected->SetLineColor(kBlack);
  limitValueExpected->SetLineWidth(2);
  limitValueExpected->SetLineStyle(2);
  limitValueExpected->Draw("PE0same");
  limitErr_1s->SetFillColor(kGreen+1);
  limitErr_1s->SetLineColor(kGreen+1);
  limitErr_2s->SetFillColor(kOrange);
  limitErr_2s->SetLineColor(kOrange);
  limitErr_2s->Draw("2same");
  limitErr_1s->Draw("2same");
  limitValueExpected->Draw("PE0 same");
  limitValueExpected->SetMarkerSize(1.2);
  limitValueExpected->SetMarkerStyle(24);
  limitValueExpected->SetMarkerColor(kBlack);
  limitValueExpected->Draw("PE0 same");

  if(not blind){
    limitValueObserved->SetLineColor(kBlack);
    limitValueObserved->SetMarkerColor(kBlack);
    limitValueObserved->SetMarkerSize(1.2);
    limitValueObserved->SetMarkerStyle(20);
    limitValueObserved->SetLineWidth(2);
    limitValueObserved->Draw("PE0 same");
  }

  TLegend* leg = new TLegend(0.58,0.63,0.9,0.9);
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  if(not blind)
    leg->AddEntry(limitValueObserved,"Observed","PL");
  leg->AddEntry(limitValueExpected,"Median expected","PL");
  leg->AddEntry(limitErr_1s,"68% expected","F");
  leg->AddEntry(limitErr_2s,"95% expected","F");
  leg->Draw("same");

  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",false,true,false);
  else
    CMS_lumi(canvas,"35.9",false,false,false);

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs("higgsInv_brazilian.png","png");
  canvas->SaveAs("higgsInv_brazilian.pdf","pdf");
   
  

}
