#include "../CMS_lumi.h"

void makeHiggsInvisibleComparison(){

  setTDRStyle();

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 675);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);

  TH1F* limitValueExpected = new TH1F("limitValueExpected","",4,0,4);
  limitValueExpected->SetBinContent(1,0.8477);
  limitValueExpected->SetBinContent(2,1.0820);
  limitValueExpected->SetBinContent(3,1.1758);
  limitValueExpected->SetBinContent(4,1.0664);
  //  limitValueExpected->SetBinContent(5,0.8867);
  //  limitValueExpected->SetBinContent(2,0.9610);
  //  limitValueExpected->SetBinContent(3,1.0273);
  //  limitValueExpected->SetBinContent(4,0.9570);

  TH1F* limitValueObserved = new TH1F("limitValueObserved","",4,0,4);
  limitValueObserved->SetBinContent(1,0.4826);
  limitValueObserved->SetBinContent(2,0.5686);
  limitValueObserved->SetBinContent(3,0.7528);
  limitValueObserved->SetBinContent(4,0.7752);
  //  limitValueObserved->SetBinContent(5,0.4561);
  //  limitValueObserved->SetBinContent(2,0.7628);
  //  limitValueObserved->SetBinContent(3,1.2895);
  //  limitValueObserved->SetBinContent(4,0.7882);
				    
  TGraphAsymmErrors* limitErr_1s = new TGraphAsymmErrors();
  limitErr_1s->SetPoint(1,0.5,0.8477);
  limitErr_1s->SetPoint(2,1.5,1.0820);
  limitErr_1s->SetPoint(3,2.5,1.1758);
  limitErr_1s->SetPoint(4,3.5,1.0664);
  //  limitErr_1s->SetPoint(5,4.5,0.8867);
  //  limitErr_1s->SetPoint(2,1.5,0.9610);
  //  limitErr_1s->SetPoint(3,2.5,1.0273);
  //  limitErr_1s->SetPoint(4,3.5,0.9570);

  limitErr_1s->SetPointError(1,0.5,0.5,fabs(0.8477-0.5842),fabs(0.8477-1.2666));
  limitErr_1s->SetPointError(2,0.5,0.5,fabs(1.0820-0.7484),fabs(1.0820-1.6125));
  limitErr_1s->SetPointError(3,0.5,0.5,fabs(1.1758-0.8133),fabs(1.1758-1.7522));
  limitErr_1s->SetPointError(4,0.5,0.5,fabs(1.0664-0.7350),fabs(1.0664-1.5977));
  //  limitErr_1s->SetPointError(5,0.5,0.5,fabs(0.8867-0.6133),fabs(0.8867-1.3215));
  //  limitErr_1s->SetPointError(2,0.5,0.5,fabs(0.9610-0.6680),fabs(0.9610-1.4262));
  //  limitErr_1s->SetPointError(3,0.5,0.5,fabs(1.0273-0.7081),fabs(1.0273-1.5351));
  //  limitErr_1s->SetPointError(4,0.5,0.5,fabs(0.9570-0.6643),fabs(0.9570-1.4262));

  TGraphAsymmErrors* limitErr_2s = new TGraphAsymmErrors();
  limitErr_2s->SetPoint(1,0.5,0.8477);
  limitErr_2s->SetPoint(2,1.5,1.0820);
  limitErr_2s->SetPoint(3,2.5,1.1758);
  limitErr_2s->SetPoint(4,3.5,1.0664);
  //  limitErr_2s->SetPoint(5,4.5,0.8867);
  //  limitErr_2s->SetPoint(2,1.5,0.9610);
  //  limitErr_2s->SetPoint(3,2.5,1.0273);
  //  limitErr_2s->SetPoint(4,3.5,0.9570);

  limitErr_2s->SetPointError(1,0.5,0.5,fabs(0.8477-0.4188),fabs(0.8477-1.8555));
  limitErr_2s->SetPointError(2,0.5,0.5,fabs(1.0820-0.5516),fabs(1.0820-2.3662));
  limitErr_2s->SetPointError(3,0.5,0.5,fabs(1.1758-0.5994),fabs(1.1758-2.5427));
  limitErr_2s->SetPointError(4,0.5,0.5,fabs(1.0664-0.5395),fabs(1.0664-2.3367));
  //  limitErr_2s->SetPointError(5,0.5,0.5,fabs(0.8867-0.4520),fabs(0.8867-1.9284));
  //  limitErr_2s->SetPointError(2,1.5,0.5,fabs(0.9610-0.4916),fabs(0.9610-2.0813));
  //  limitErr_2s->SetPointError(3,2.5,0.5,fabs(1.0273-0.5197),fabs(1.0273-2.2489));
  //  limitErr_2s->SetPointError(4,3.5,0.5,fabs(0.9570-0.4916),fabs(0.9570-2.0697));
  
  limitValueExpected->GetXaxis()->SetBinLabel(1,"nominal");
  limitValueExpected->GetXaxis()->SetBinLabel(2,"Z#mu#mu+W#mu#nu only");
  limitValueExpected->GetXaxis()->SetBinLabel(3,"Zee+We#nu only");
  limitValueExpected->GetXaxis()->SetBinLabel(4,"ele+#gamma");
  //  limitValueExpected->GetXaxis()->SetBinLabel(5,"#mu+ele only");
  //  limitValueExpected->GetXaxis()->SetBinLabel(2,"MET<800");
  //  limitValueExpected->GetXaxis()->SetBinLabel(3,"MET>350");
  //  limitValueExpected->GetXaxis()->SetBinLabel(4,"merge MET>800");

  limitValueExpected->SetLineColor(kBlack);
  limitValueExpected->SetLineWidth(2);
  limitValueExpected->SetLineStyle(2);
  limitValueExpected->GetYaxis()->SetRangeUser(0.,3.5);
  limitValueExpected->GetYaxis()->SetLabelSize(0.035);
  limitValueExpected->GetXaxis()->SetLabelSize(0.04);
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
