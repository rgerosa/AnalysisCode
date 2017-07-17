#include "../CMS_lumi.h"

enum class Sample {zjet,wjet};

void makeControlRegionComparison(string inputFileName, string outputDIR, string observable, string observableLatex, Sample sample, bool isEWK){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TH1* signalregion  = NULL;
  TH1* controlregion1 = NULL;
  TH1* controlregion2 = NULL;

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();

  if(sample == Sample::zjet){
    if(not isEWK){
      signalregion   = (TH1*) inputFile->FindObjectAny(("zinvhist_"+observable).c_str());
      controlregion1 = (TH1*) inputFile->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
      controlregion2 = (TH1*) inputFile->FindObjectAny(("vllbkghistzee_"+observable).c_str());
    }
    else{
      signalregion   = (TH1*) inputFile->FindObjectAny(("ewkbkgzhist_"+observable).c_str());
      controlregion1 = (TH1*) inputFile->FindObjectAny(("ewkzbkghistzmm_"+observable).c_str());
      controlregion2 = (TH1*) inputFile->FindObjectAny(("ewkzbkghistzee_"+observable).c_str());
    }
  }
  else if(sample == Sample::wjet){
    if(not isEWK){
      signalregion   = (TH1*) inputFile->FindObjectAny(("wjethist_"+observable).c_str());
      controlregion1 = (TH1*) inputFile->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
      controlregion2 = (TH1*) inputFile->FindObjectAny(("vlbkghistwen_"+observable).c_str());
    }
    else{
      signalregion   = (TH1*) inputFile->FindObjectAny(("ewkbkgwhist_"+observable).c_str());
      controlregion1 = (TH1*) inputFile->FindObjectAny(("ewkwbkghistwmn_"+observable).c_str());
      controlregion2 = (TH1*) inputFile->FindObjectAny(("ewkwbkghistwen_"+observable).c_str());
    }
  }

  signalregion->Scale(1./signalregion->Integral());
  controlregion1->Scale(1./controlregion1->Integral());
  controlregion2->Scale(1./controlregion2->Integral());

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  canvas->SetBottomMargin(0.3);
  canvas->cd();
  
  signalregion->GetXaxis()->SetTitle(observableLatex.c_str());
  signalregion->GetYaxis()->SetTitle("a.u.");
  signalregion->GetXaxis()->SetTitleSize(0);
  signalregion->GetXaxis()->SetLabelSize(0);

  signalregion->SetLineColor(kBlack);
  signalregion->SetLineWidth(2);

  controlregion1->SetLineColor(kRed);
  controlregion1->SetLineWidth(2);

  controlregion2->SetLineColor(kBlue);
  controlregion2->SetLineWidth(2);

  signalregion->Draw("hist");
  controlregion1->Draw("hist same");
  controlregion2->Draw("hist same");

  TLegend leg (0.7,0.7,0.92,0.92);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(sample == Sample::zjet){
    leg.AddEntry(signalregion,"Z #rightarrow #nu#nu SR","L");
    leg.AddEntry(controlregion1,"Z #rightarrow #mu#mu CR","L");
    leg.AddEntry(controlregion2,"Z #rightarrow ee CR","L");
  }
  else{
    leg.AddEntry(signalregion,"W+jets SR","L");
    leg.AddEntry(controlregion1,"W #rightarrow #mu#nu CR","L");
    leg.AddEntry(controlregion2,"W #rightarrow e#nu CR","L");
  }

  leg.Draw("same");
  CMS_lumi(canvas,"35.9");
  
  signalregion->GetYaxis()->SetRangeUser(min(signalregion->GetMinimum(),min(controlregion1->GetMinimum(),controlregion2->GetMinimum()))*0.1,
					 max(signalregion->GetMaximum(),max(controlregion1->GetMaximum(),controlregion2->GetMaximum()))*10);

  
  canvas->SetLogy();

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* ratio1 = (TH1*) controlregion1->Clone("ratio1");
  TH1* ratio2 = (TH1*) controlregion2->Clone("ratio2");
  
  ratio1->Divide(signalregion);
  ratio2->Divide(signalregion);
    
  ratio1->GetYaxis()->SetNdivisions(5);
  ratio1->GetYaxis()->CenterTitle();
  ratio1->GetYaxis()->SetTitleOffset(1.5);
  ratio1->GetYaxis()->SetLabelSize(0.04);
  ratio1->GetYaxis()->SetTitleSize(0.04);
  ratio1->GetXaxis()->SetLabelSize(0.04);
  ratio1->GetXaxis()->SetTitleSize(0.05);
  ratio1->GetXaxis()->SetTitleOffset(1.1);
  ratio1->GetXaxis()->SetTitle(observableLatex.c_str());
  ratio1->GetYaxis()->SetTitle("Ratio");
  ratio1->Draw("hist");
  ratio2->Draw("hist same");
  if(sample == Sample::zjet)
    ratio1->GetYaxis()->SetRangeUser(0.8,1.2);
  else if(sample == Sample::wjet)
    ratio1->GetYaxis()->SetRangeUser(0.8,1.2);
  
  TH1* uncertainty = (TH1*) controlregion1->Clone("uncertainty");
  uncertainty->Divide(signalregion);
  uncertainty->SetMarkerSize(0);
  for(int iBin = 0; iBin < uncertainty->GetNbinsX(); iBin++)
    uncertainty->SetBinContent(iBin+1,1);

  uncertainty->SetFillColor(kGray+1);
  uncertainty->SetFillStyle(1001);
  uncertainty->Draw("E2same");
  ratio1->Draw("hist same");
  ratio2->Draw("hist same");

  if(not isEWK){
    if(sample == Sample::zjet){
      canvas->SaveAs((outputDIR+"/zjets_comparison_"+observable+".png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/zjets_comparison_"+observable+".pdf").c_str(),"pdf");
    }
    else if(sample == Sample::wjet){
      canvas->SaveAs((outputDIR+"/wjets_comparison_"+observable+".png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/wjets_comparison_"+observable+".pdf").c_str(),"pdf");
    }
  }
  else{
    if(sample == Sample::zjet){
      canvas->SaveAs((outputDIR+"/zjets_comparison_"+observable+"_ewk.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/zjets_comparison_"+observable+"_ewk.pdf").c_str(),"pdf");
    }
    else if(sample == Sample::wjet){
      canvas->SaveAs((outputDIR+"/wjets_comparison_"+observable+"_ewk.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/wjets_comparison_"+observable+"_ewk.pdf").c_str(),"pdf");
    }
  }
}
