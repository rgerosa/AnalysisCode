#include "../CMS_lumi.h"

enum class Sample {znn, wjet, gam};

void makeKFactotComparison(string outputDirectory, Sample sample){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());

  TFile* kFactorFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis//data/kFactors/uncertainties_EWK_24bins.root");

  TH1F* QCD_NLO = NULL;
  TH1F* QCD_LO  = NULL;
  TH1F* QCD_EWK_NLO  = NULL;

  if(sample == Sample::znn){
    QCD_NLO = (TH1F*) kFactorFile->Get("ZJets_012j_NLO/nominal");
    QCD_LO  = (TH1F*) kFactorFile->Get("ZJets_LO/inv_pt");
    QCD_EWK_NLO = (TH1F*) kFactorFile->Get("EWKcorr/Z");
  }
  else if(sample == Sample::wjet){
    QCD_NLO = (TH1F*) kFactorFile->Get("WJets_012j_NLO/nominal");
    QCD_LO  = (TH1F*) kFactorFile->Get("WJets_LO/inv_pt");
    QCD_EWK_NLO = (TH1F*) kFactorFile->Get("EWKcorr/W");    
  }
  else if(sample == Sample::gam){
    QCD_NLO = (TH1F*) kFactorFile->Get("GJets_01j_NLO/nominal");
    QCD_LO  = (TH1F*) kFactorFile->Get("GJets_LO/inv_pt");
    QCD_EWK_NLO = (TH1F*) kFactorFile->Get("EWKcorr/photon");    
  }

  QCD_NLO->Rebin(2);
  QCD_LO->Rebin(2);
  QCD_EWK_NLO->Rebin(2);
  
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetLogy();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.27);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->SetBottomMargin(0.35);

  TH1F* frame =  canvas->DrawFrame(QCD_NLO->GetXaxis()->GetXmin(),0.00001,1090,100, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV]");
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame ->Draw();

  QCD_LO->SetLineColor(kBlack);
  QCD_LO->SetMarkerColor(kBlack);
  QCD_LO->SetLineWidth(2);
  QCD_NLO->SetLineColor(kRed);
  QCD_NLO->SetMarkerColor(kRed);
  QCD_NLO->SetLineWidth(2);
  QCD_EWK_NLO->SetLineColor(kBlue);
  QCD_EWK_NLO->SetMarkerColor(kBlue);
  QCD_EWK_NLO->SetLineWidth(2);

  QCD_LO->Draw("hist same");
  QCD_NLO->Draw("hist same");
  QCD_EWK_NLO->Draw("hist same");

  CMS_lumi(canvas,"");

  TLegend* leg = new TLegend(0.48, 0.66, 0.90, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  if(sample == Sample::znn)
    leg->AddEntry((TObject*)(0),"Z+jets","");
  else if(sample == Sample::wjet)
    leg->AddEntry((TObject*)(0),"W+jets","");
  else if(sample == Sample::gam)
    leg->AddEntry((TObject*)(0),"#gamma+jets","");

  leg->AddEntry(QCD_EWK_NLO,"NLO QCD + EWK","FL");
  leg->AddEntry(QCD_NLO,"NLO-QCD","FL");
  leg->AddEntry(QCD_LO,"LO","FL");  
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");
    
  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1F* frame2 =  pad2->DrawFrame(QCD_NLO->GetXaxis()->GetXmin(),0.5,1090,2.0, "");
  
  frame2->GetXaxis()->SetLabelSize(0.12);
  frame2->GetXaxis()->SetTitleSize(0.15);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.12);
  frame2->GetYaxis()->SetTitleSize(0.15);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetXaxis()->SetTitle("boson p_{T} [GeV]");
  frame2->GetYaxis()->SetTitle("#sigma_{NLO}/#sigma_{LO}");
  frame2->GetYaxis()->SetNdivisions(505);
  frame2->Draw();

  TH1* ratio = (TH1*) QCD_NLO->Clone("ratio");
  ratio->Divide(QCD_LO);

  TH1* ratio2 = (TH1*) QCD_EWK_NLO->Clone("ratio2");
  ratio2->Divide(QCD_LO);

  TF1* line = new TF1("line","1",QCD_NLO->GetXaxis()->GetXmin(),QCD_NLO->GetXaxis()->GetXmax());
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  
  ratio->SetLineWidth(1);
  ratio->SetMarkerSize(1);
  ratio->SetMarkerStyle(20);

  ratio2->SetLineWidth(1);
  ratio2->SetMarkerSize(1);
  ratio2->SetMarkerStyle(20);

  ratio->Draw("Psame");
  ratio2->Draw("Psame");
  line->Draw("same");
  pad2->RedrawAxis("sameaxis");
  
  if(sample == Sample::znn){
    canvas->SaveAs((outputDirectory+"/Zjets_NLO_QCD_EWK.pdf").c_str(),"pdf");
    canvas->SaveAs((outputDirectory+"/Zjets_NLO_QCD_EWK.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/Zjets_NLO_QCD_EWK.root").c_str(),"root");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDirectory+"/Wjets_NLO_QCD_EWK.pdf").c_str(),"pdf");
    canvas->SaveAs((outputDirectory+"/Wjets_NLO_QCD_EWK.png").c_str(),"png");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDirectory+"/Gjets_NLO_QCD_EWK.pdf").c_str(),"pdf");
    canvas->SaveAs((outputDirectory+"/Gjets_NLO_QCD_EWK.png").c_str(),"png");
  }
}
