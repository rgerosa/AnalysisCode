void makeKFactotComparison(string outputDirectory){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  system(("mkdir -p "+outputDirectory).c_str());

  TFile* kFactorFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis//data/kFactors/uncertainties_EWK_24bins.root");
  
  TH1F* ZJets_QCD_NLO = (TH1F*) kFactorFile->Get("ZJets_012j_NLO/nominal");
  TH1F* ZJets_QCD_LO  = (TH1F*) kFactorFile->Get("ZJets_LO/inv_pt");
  TH1F* ZJets_EWKCorr = (TH1F*) kFactorFile->Get("EWKcorr/Z");

  TH1F* WJets_QCD_NLO = (TH1F*) kFactorFile->Get("WJets_012j_NLO/nominal");
  TH1F* WJets_QCD_LO  = (TH1F*) kFactorFile->Get("WJets_LO/inv_pt");
  TH1F* WJets_EWKCorr = (TH1F*) kFactorFile->Get("EWKcorr/W");

  ZJets_EWKCorr->Divide(ZJets_QCD_NLO);
  WJets_EWKCorr->Divide(WJets_QCD_NLO);

  ZJets_QCD_NLO->Divide(ZJets_QCD_LO);
  WJets_QCD_NLO->Divide(WJets_QCD_LO);

  TCanvas* canvas = new TCanvas("canvas","",500,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.27);
  pad2->SetTickx();
  pad2->SetTicky();

  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();

  TH1F* frame =  pad1->DrawFrame(ZJets_QCD_NLO->GetXaxis()->GetXmin(),0.4,ZJets_QCD_NLO->GetXaxis()->GetXmax(),2.0, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("K-factors");
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame ->Draw();
  
  ZJets_QCD_NLO->SetLineColor(kRed);
  ZJets_QCD_NLO->SetLineWidth(2);

  TH1F* ZJets_QCD_NLO_band = (TH1F*) ZJets_QCD_NLO->Clone("ZJets_QCD_NLO_band");
  ZJets_QCD_NLO_band->SetFillColor(kRed);
  ZJets_QCD_NLO_band->SetFillStyle(3001);

  WJets_QCD_NLO->SetLineColor(kBlue);
  WJets_QCD_NLO->SetLineWidth(2);

  TH1F* WJets_QCD_NLO_band = (TH1F*) WJets_QCD_NLO->Clone("WJets_QCD_NLO_band");
  WJets_QCD_NLO_band->SetFillColor(kBlue);
  WJets_QCD_NLO_band->SetFillStyle(3001);

  ZJets_QCD_NLO_band->Draw("E2 same");
  WJets_QCD_NLO_band->Draw("E2 same");
  ZJets_QCD_NLO->Draw("hist same");
  WJets_QCD_NLO->Draw("hist same");

  TLegend* leg = new TLegend(0.48, 0.66, 0.90, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(ZJets_QCD_NLO,"ZJets: NLO/LO QCD","FL");
  leg->AddEntry(WJets_QCD_NLO,"WJets: NLO/LO QCD","FL");
  
  leg->Draw("same");
  pad1->RedrawAxis("sameaxis");
    
  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1F* frame2 =  pad2->DrawFrame(ZJets_QCD_NLO->GetXaxis()->GetXmin(),0.5,ZJets_QCD_NLO->GetXaxis()->GetXmax(),1.75, "");
  
  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleSize(0.13);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleSize(0.10);
  frame2->GetXaxis()->SetTitle("boson p_{T} [GeV]");
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->GetYaxis()->SetTitle("ZJets/WJets");
  frame2->GetYaxis()->SetTitleOffset(0.5);
  frame2->Draw();

  TH1* ratio = (TH1*) ZJets_QCD_NLO->Clone("ratio");
  ratio->Divide(WJets_QCD_NLO);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineWidth(1);
  ratio->SetMarkerSize(0.8);
  ratio->SetMarkerStyle(20);
  ratio->Draw("Psame");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDirectory+"/ZoverW_NLO_QCD.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/ZoverW_NLO_QCD.png").c_str(),"png");

  //make QCD*EWK
  TH1F* ZJets_EWK_QCD = (TH1F*) ZJets_QCD_NLO->Clone("ZJets_EWK_QCD");
  for(int iBin = 0; iBin < ZJets_EWK_QCD->GetNbinsX(); iBin++)
    ZJets_EWK_QCD->SetBinContent(iBin+1,ZJets_QCD_NLO->GetBinContent(iBin+1)*ZJets_EWKCorr->GetBinContent(iBin+1));

  TH1F* WJets_EWK_QCD = (TH1F*) WJets_QCD_NLO->Clone("WJets_EWK_QCD");
  for(int iBin = 0; iBin < WJets_EWK_QCD->GetNbinsX(); iBin++)
    WJets_EWK_QCD->SetBinContent(iBin+1,WJets_QCD_NLO->GetBinContent(iBin+1)*WJets_EWKCorr->GetBinContent(iBin+1));
 
  canvas->cd();
  pad1->cd();
  frame->Draw();

  TH1F* ZJets_EWK_QCD_band = (TH1F*) ZJets_EWK_QCD->Clone("ZJets_EWK_QCD_band");
  ZJets_EWK_QCD_band->SetFillColor(kRed);
  ZJets_EWK_QCD_band->SetFillStyle(3001);

  TH1F* WJets_EWK_QCD_band = (TH1F*) WJets_EWK_QCD->Clone("WJets_EWK_QCD_band");
  WJets_EWK_QCD_band->SetFillColor(kBlue);
  WJets_EWK_QCD_band->SetFillStyle(3001);
  
  ZJets_EWK_QCD_band->Draw("E2 same");
  WJets_EWK_QCD_band->Draw("E2 same");

  ZJets_EWK_QCD->Draw("hist same");
  WJets_EWK_QCD->Draw("hist same");

  leg->Clear();
  leg->AddEntry(ZJets_EWK_QCD,"ZJets: NLO/LO QCD+EWK","FL");
  leg->AddEntry(WJets_EWK_QCD,"WJets: NLO/LO QCD+EWK","FL");
  leg->Draw("same");

  canvas->cd();
  pad2->cd();
  frame2->Draw();

  ratio = (TH1*) ZJets_EWK_QCD->Clone("ratio");
  ratio->Divide(WJets_EWK_QCD);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineWidth(1);
  ratio->SetMarkerSize(0.8);
  ratio->SetMarkerStyle(20);
  ratio->Draw("Psame");

  canvas->SaveAs((outputDirectory+"/ZoverW_NLO_EWK_QCD.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/ZoverW_NLO_EWK_QCD.png").c_str(),"png");


}
