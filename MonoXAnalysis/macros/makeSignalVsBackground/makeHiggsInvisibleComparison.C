#include "../CMS_lumi.h"

void plotHist(TCanvas* canvas, TH1* vbf, TH1* ggh, TH1* zvv, TH1* wjet, TH1* ewkz, TH1* ewkw, string outputDIR, string postfix, string postfix_latex, bool isLog = true){

  zvv->Add(wjet);
  ewkz->Add(ewkw);

  vbf->Scale(1./vbf->Integral());
  ggh->Scale(1./ggh->Integral());
  zvv->Scale(1./zvv->Integral());
  ewkz->Scale(1./ewkz->Integral());

  vbf->SetLineColor(kBlack);
  vbf->SetLineWidth(3);

  ggh->SetLineColor(kGray+2);
  ggh->SetLineWidth(3);

  zvv->SetLineColor(kRed);
  zvv->SetFillColor(kRed);
  zvv->SetFillStyle(3001);
  zvv->SetLineWidth(2);
  ewkz->SetLineColor(kBlue);
  ewkz->SetFillColor(kBlue);
  ewkz->SetFillStyle(3001);
  ewkz->SetLineWidth(2);
  

  vbf->GetXaxis()->SetTitle(postfix_latex.c_str());
  vbf->GetYaxis()->SetTitle("a.u.");

  if(isLog){
    canvas->SetLogy();
    vbf->GetYaxis()->SetRangeUser(0.001,10);
  }
  else
    vbf->GetYaxis()->SetRangeUser(0,0.5);

  vbf->Draw("hist");
  zvv->Draw("hist same");
  ewkz->Draw("hist same");
  ggh->Draw("hist same");
  vbf->Draw("hist same");

  TLegend leg (0.6,0.5,0.85,0.85);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(vbf,"VBF H_{125}","L");
  leg.AddEntry(ggh,"ggH H_{125}","L");
  leg.AddEntry(zvv,"V+jets QCD","L");
  leg.AddEntry(ewkz,"V+jets EW","L");
  leg.Draw("same");  

  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  canvas->SetLogy(0);

}

void makeHiggsInvisibleComparison(string inputFileName_sig, string inputFileName_bkg, string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TFile* inputFile_sig = TFile::Open(inputFileName_sig.c_str(),"READ");
  TFile* inputFile_bkg = TFile::Open(inputFileName_bkg.c_str(),"READ");
  inputFile_sig->cd();
  inputFile_bkg->cd();
  
  TH1F* vbf_met = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_met");
  TH1F* vbf_mjj = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_mjj");
  TH1F* vbf_jetpt = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_jetpt");
  TH1F* vbf_jetpt2 = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_jetpt2");
  TH1F* vbf_detajj = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_detajj");
  TH1F* vbf_dphijj = (TH1F*) inputFile_sig->Get("vbfH/vbfHhist_125_dphiJJ");

  TH1F* ggh_met = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_met");
  TH1F* ggh_mjj = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_mjj");
  TH1F* ggh_jetpt = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_jetpt");
  TH1F* ggh_jetpt2 = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_jetpt2");
  TH1F* ggh_detajj = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_detajj");
  TH1F* ggh_dphijj = (TH1F*) inputFile_sig->Get("ggH/ggHhist_125_dphiJJ");

  TH1F* zvv_met = (TH1F*) inputFile_bkg->Get("SR/zinvhist_met");
  TH1F* zvv_mjj = (TH1F*) inputFile_bkg->Get("SR/zinvhist_mjj");
  TH1F* zvv_jetpt = (TH1F*) inputFile_bkg->Get("SR/zinvhist_jetpt");
  TH1F* zvv_jetpt2 = (TH1F*) inputFile_bkg->Get("SR/zinvhist_jetpt2");
  TH1F* zvv_detajj = (TH1F*) inputFile_bkg->Get("SR/zinvhist_detajj");
  TH1F* zvv_dphijj = (TH1F*) inputFile_bkg->Get("SR/zinvhist_dphiJJ");

  TH1F* wjet_met = (TH1F*) inputFile_bkg->Get("SR/wjethist_met");
  TH1F* wjet_mjj = (TH1F*) inputFile_bkg->Get("SR/wjethist_mjj");
  TH1F* wjet_jetpt = (TH1F*) inputFile_bkg->Get("SR/wjethist_jetpt");
  TH1F* wjet_jetpt2 = (TH1F*) inputFile_bkg->Get("SR/wjethist_jetpt2");
  TH1F* wjet_detajj = (TH1F*) inputFile_bkg->Get("SR/wjethist_detajj");
  TH1F* wjet_dphijj = (TH1F*) inputFile_bkg->Get("SR/wjethist_dphiJJ");

  TH1F* ewkz_met = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_met");
  TH1F* ewkz_mjj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_mjj");
  TH1F* ewkz_jetpt = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_jetpt");
  TH1F* ewkz_jetpt2 = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_jetpt2");
  TH1F* ewkz_detajj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_detajj");
  TH1F* ewkz_dphijj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgzhist_dphiJJ");

  TH1F* ewkw_met = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_met");
  TH1F* ewkw_mjj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_mjj");
  TH1F* ewkw_jetpt = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_jetpt");
  TH1F* ewkw_jetpt2 = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_jetpt2");
  TH1F* ewkw_detajj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_detajj");
  TH1F* ewkw_dphijj = (TH1F*) inputFile_bkg->Get("SR/ewkbkgwhist_dphiJJ");

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();

  plotHist(canvas,vbf_met,ggh_met,zvv_met,wjet_met,ewkz_met,ewkw_met,outputDIR,"met","E_{T}^{miss} [GeV]");
  plotHist(canvas,vbf_mjj,ggh_mjj,zvv_mjj,wjet_mjj,ewkz_mjj,ewkw_mjj,outputDIR,"mjj","M_{jj} [GeV]");
  plotHist(canvas,vbf_jetpt,ggh_jetpt,zvv_jetpt,wjet_jetpt,ewkz_jetpt,ewkw_jetpt,outputDIR,"jetpt","p_{T}^{j1} [GeV]");
  plotHist(canvas,vbf_jetpt2,ggh_jetpt2,zvv_jetpt2,wjet_jetpt2,ewkz_jetpt2,ewkw_jetpt2,outputDIR,"jetpt2","p_{T}^{j2} [GeV]");
  plotHist(canvas,vbf_detajj,ggh_detajj,zvv_detajj,wjet_detajj,ewkz_detajj,ewkw_detajj,outputDIR,"detajj","#Delta#eta_{jj}",false);
  plotHist(canvas,vbf_dphijj,ggh_dphijj,zvv_dphijj,wjet_dphijj,ewkz_dphijj,ewkw_dphijj,outputDIR,"dphijj","#Delta#phi_{jj}",false);

}
