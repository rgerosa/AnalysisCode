#include "../CMS_lumi.h"

enum class Sample {zw,zgamma};

static float minX = 200;
static float maxX = 1400;

///////////////////
void makeVJetsNLOComparisonWithTheory(string inputFileNameNum, string inputFileNameDen, string outputDirectory, Sample sample){

  system(("mkdir -p "+outputDirectory).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);


  TFile* inputFileNum = TFile::Open(inputFileNameNum.c_str(),"READ");
  TFile* inputFileDen = TFile::Open(inputFileNameDen.c_str(),"READ");

  ////// ----------
  TH1F* num_cms_lhe  = (TH1F*) inputFileNum->Get("bosonpt_zvv_cms_lhe");
  TH1F* num_cms_ps   = (TH1F*) inputFileNum->Get("bosonpt_zvv_cms_ps");
  TH1F* num_theo     = (TH1F*) inputFileNum->Get("bosonpt_zvv_theo");

  TH1F* den_cms_lhe  = (TH1F*) inputFileDen->Get("bosonpt_wjet_cms_lhe");
  TH1F* den_cms_ps   = (TH1F*) inputFileDen->Get("bosonpt_wjet_cms_lhe");
  TH1F* den_theo     = (TH1F*) inputFileDen->Get("bosonpt_wjet_theo");

  if(den_cms_lhe == 0)
    den_cms_lhe  = (TH1F*) inputFileDen->Get("bosonpt_gam_cms_lhe");
  if(den_cms_ps == 0)
    den_cms_ps   = (TH1F*) inputFileDen->Get("bosonpt_gam_cms_ps");
  if(den_theo == 0)
    den_theo     = (TH1F*) inputFileDen->Get("bosonpt_gam_theo");

  ///// ----------
  TH1F* ratio_cms_lhe = (TH1F*) num_cms_lhe->Clone("ratio_cms_lhe");
  ratio_cms_lhe->Divide(den_cms_lhe);

  TH1F* ratio_cms_ps = (TH1F*) num_cms_ps->Clone("ratio_cms_ps");
  ratio_cms_ps->Divide(den_cms_ps);

  TH1F* ratio_theo = (TH1F*) num_theo->Clone("ratio_theo");
  ratio_theo->Divide(den_theo);

  ///// ----------
  for(int iBin = 0; iBin < ratio_cms_lhe->GetNbinsX(); iBin++)
    ratio_cms_lhe->SetBinError(iBin+1,0.);
  for(int iBin = 0; iBin < ratio_cms_ps->GetNbinsX(); iBin++)
    ratio_cms_ps->SetBinError(iBin+1,0.);
  for(int iBin = 0; iBin < ratio_theo->GetNbinsX(); iBin++)
    ratio_theo->SetBinError(iBin+1,0.);


  // sety few style stuff
  ratio_cms_lhe->SetLineWidth(2);
  ratio_cms_lhe->SetLineColor(kBlack);

  ratio_cms_ps->SetLineWidth(2);
  ratio_cms_ps->SetLineColor(kRed);

  ratio_theo->SetLineWidth(2);
  ratio_theo->SetLineColor(kBlue);

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();
  
  ratio_cms_lhe->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  if(sample == Sample::zw)
    ratio_cms_lhe->GetYaxis()->SetTitle("Z+jet/W+jets");
  else if(sample == Sample::zgamma)
    ratio_cms_lhe->GetYaxis()->SetTitle("Z+jet/#gamam+jets");

  ratio_cms_lhe->GetXaxis()->SetRangeUser(minX,maxX);
  ratio_cms_lhe->Draw("hist");
  ratio_cms_ps->Draw("hist same");
  ratio_theo->Draw("hist same");

  if(sample == Sample::zw){
    canvas->SaveAs((outputDirectory+"/zw_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/zw_ratio.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::zgamma){
    canvas->SaveAs((outputDirectory+"/zgamma_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/zgamma_ratio.pdf").c_str(),"pdf");
  }
}
