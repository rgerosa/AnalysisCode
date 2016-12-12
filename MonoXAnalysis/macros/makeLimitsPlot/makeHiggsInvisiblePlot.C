#include "../CMS_lumi.h"

void makeHiggsInvisiblePlot(bool blind = false){

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
  limitValueExpected->SetPoint(0,0.5,0.5254);
  limitValueExpected->SetPoint(1,1.5,0.4824);  
  limitValueExpected->SetPoint(2,2.5,0.1753);
  limitValueExpected->SetPoint(3,3.5,0.1528);
  limitValueExpected->SetPointError(0,0.5,0.);  
  limitValueExpected->SetPointError(1,0.5,0.);
  limitValueExpected->SetPointError(2,0.5,0.);
  limitValueExpected->SetPointError(3,0.5,0.);

  TGraphErrors* limitValueObserved = new TGraphErrors();
  limitValueObserved->SetPoint(0,0.5,0.4826);
  limitValueObserved->SetPoint(1,1.5,1.1692);
  limitValueObserved->SetPoint(2,2.5,0.1528);
  limitValueObserved->SetPoint(3,3.5,0.1528);
  limitValueObserved->SetPointError(0,0.5,0.);  
  limitValueObserved->SetPointError(1,0.5,0.);
  limitValueObserved->SetPointError(2,0.5,0.);
  limitValueObserved->SetPointError(3,0.5,0.);

				    
  TGraphAsymmErrors* limitErr_1s = new TGraphAsymmErrors();
  limitErr_1s->SetPoint(0,0.5,0.5254);
  limitErr_1s->SetPoint(1,1.5,0.4824);
  limitErr_1s->SetPoint(2,2.5,0.1753);
  limitErr_1s->SetPoint(3,3.5,0.1528);
  limitErr_1s->SetPointError(0,0.5,0.5,fabs(0.5254-0.3762),fabs(0.5254-0.7411));
  limitErr_1s->SetPointError(1,0.5,0.5,fabs(0.4824-0.3337),fabs(0.4824-0.7228));
  limitErr_1s->SetPointError(2,0.5,0.5,fabs(0.1753-0.1261),fabs(0.1753-0.2459));
  limitErr_1s->SetPointError(3,0.5,0.5,fabs(0.1528-0.1094),fabs(0.1528-0.2144));

  TGraphAsymmErrors* limitErr_2s = new TGraphAsymmErrors();
  limitErr_2s->SetPoint(0,0.5,0.5254);
  limitErr_2s->SetPoint(1,1.5,0.4824);
  limitErr_2s->SetPoint(2,2.5,0.1753);
  limitErr_2s->SetPoint(3,3.5,0.1528);
  limitErr_2s->SetPointError(0,0.5,0.5,fabs(0.5254-0.2822),fabs(0.5254-1.0121));
  limitErr_2s->SetPointError(1,0.5,0.5,fabs(0.4824-0.2459),fabs(0.4824-1.0071));
  limitErr_2s->SetPointError(2,0.5,0.5,fabs(0.1753-0.0955),fabs(0.1753-0.3300));
  limitErr_2s->SetPointError(3,0.5,0.5,fabs(0.1528-0.0821),fabs(0.1528-0.2897));

  TH1* frame = (TH1*) canvas->DrawFrame(0.,0.,4,1.5);
  frame->SetBins(4,0,4);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetBinLabel(1,"mono-V");
  frame->GetXaxis()->SetBinLabel(2,"monojet");
  frame->GetXaxis()->SetBinLabel(3,"VBF");
  frame->GetXaxis()->SetBinLabel(4,"Combination");
  frame->GetYaxis()->SetRangeUser(0.,1.5);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitle("95% C.L. upper limit on #mu = #sigma x BR/#sigma_{SM}");
  frame->GetYaxis()->SetTitleSize(0.038);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->Draw();

  limitValueExpected->SetLineColor(kBlack);
  limitValueExpected->SetLineWidth(2);
  limitValueExpected->SetLineStyle(2);
  limitValueExpected->Draw("PE0same");
  limitErr_1s->SetFillColor(kGreen);
  limitErr_2s->SetFillColor(kYellow);
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

  TLegend* leg = new TLegend(0.5,0.66,0.9,0.9);
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  if(not blind)
    leg->AddEntry(limitValueObserved,"Oberved CL_{s}","PL");
  leg->AddEntry(limitValueExpected,"Median Expected CL_{s}","PL");
  leg->AddEntry(limitErr_1s,"Expected CL_{s} #pm 1#sigma","F");
  leg->AddEntry(limitErr_2s,"Expected CL_{s} #pm 2#sigma","F");
  leg->Draw("same");
  CMS_lumi(canvas,"35.9",false);

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs("higgsInv_brazilian.png","png");
  canvas->SaveAs("higgsInv_brazilian.pdf","pdf");
   
  

}
