#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void plotHistograms(TCanvas* canvas, TH1* histo1, TH1* histo2, string outputDIR, string postfix){

  histo1->SetMarkerColor(kBlack);
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2);
  histo1->SetMarkerSize(1);
  histo1->SetMarkerStyle(20);

  histo2->SetMarkerColor(kRed);
  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2);

  histo1->GetXaxis()->SetTitle("Recoil [GeV]");
  histo1->GetYaxis()->SetTitle("Uncertainty");

  histo1->Draw("P");
  histo2->Draw("hist same");

  histo1->GetYaxis()->SetRangeUser(0.95,1.05);
  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".C").c_str(),"C");

  

}

void compareTheoryShape(string fileName1, string fileName2, bool isWZ, bool isZG, string outputDIR, Category category){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* file1 = TFile::Open(fileName1.c_str());
  TFile* file2 = TFile::Open(fileName2.c_str());

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();

  if(isWZ){

    TH1* ZW_QCDScale_met = (TH1*) file1->Get("ZW_QCDScale_met");
    TH1* ZW_QCDShape_met = (TH1*) file1->Get("ZW_QCDShape_met");
    TH1* ZW_QCDProcess_met = (TH1*) file1->Get("ZW_QCDProcess_met");
    TH1* ZW_NNLOEWK_met = (TH1*) file1->Get("ZW_NNLOEWK_met");
    TH1* ZW_Sudakov1_met = (TH1*) file1->Get("ZW_Sudakov1_met");
    TH1* ZW_Sudakov2_met = (TH1*) file1->Get("ZW_Sudakov2_met");
    TH1* ZW_NNLOMiss1_met = (TH1*) file1->Get("ZW_NNLOMiss1_met");
    TH1* ZW_NNLOMiss2_met = (TH1*) file1->Get("ZW_NNLOMiss2_met");
    TH1* ZW_MIX_met = (TH1*) file1->Get("ZW_MIX_met");
    TH1* ZW_PDF_met = (TH1*) file1->Get("ZW_PDF_met");

    if(category == Category::monoV){
      ZW_QCDScale_met = (TH1*) file1->Get("ZW_QCDScale_met_monov");
      ZW_QCDShape_met = (TH1*) file1->Get("ZW_QCDShape_met_monov");
      ZW_QCDProcess_met = (TH1*) file1->Get("ZW_QCDProcess_met_monov");
      ZW_NNLOEWK_met = (TH1*) file1->Get("ZW_NNLOEWK_met_monov");
      ZW_Sudakov1_met = (TH1*) file1->Get("ZW_Sudakov1_met_monov");
      ZW_Sudakov2_met = (TH1*) file1->Get("ZW_Sudakov2_met_monov");
      ZW_NNLOMiss1_met = (TH1*) file1->Get("ZW_NNLOMiss1_met_monov");
      ZW_NNLOMiss2_met = (TH1*) file1->Get("ZW_NNLOMiss2_met_monov");
      ZW_MIX_met = (TH1*) file1->Get("ZW_MIX_met_monov");
      ZW_PDF_met = (TH1*) file1->Get("ZW_PDF_met_monov");      
    }

    TH1* ZW_QCDScale_met2 = (TH1*) file2->FindObjectAny("ZW_QCDScale_met");
    TH1* ZW_QCDShape_met2 = (TH1*) file2->FindObjectAny("ZW_QCDShape_met");
    TH1* ZW_QCDProcess_met2 = (TH1*) file2->FindObjectAny("ZW_QCDProcess_met");
    TH1* ZW_NNLOEWK_met2 = (TH1*) file2->FindObjectAny("ZW_NNLOEWK_met");
    TH1* ZW_Sudakov1_met2 = (TH1*) file2->FindObjectAny("ZW_Sudakov1_met");
    TH1* ZW_Sudakov2_met2 = (TH1*) file2->FindObjectAny("ZW_Sudakov2_met");
    TH1* ZW_NNLOMiss1_met2 = (TH1*) file2->FindObjectAny("ZW_NNLOMiss1_met");
    TH1* ZW_NNLOMiss2_met2 = (TH1*) file2->FindObjectAny("ZW_NNLOMiss2_met");
    TH1* ZW_MIX_met2 = (TH1*) file2->FindObjectAny("ZW_MIX_met");
    TH1* ZW_PDF_met2 = (TH1*) file2->FindObjectAny("ZW_PDF_met");

    for(int iBin = 1; iBin <= ZW_QCDScale_met2->GetNbinsX(); iBin++){
      ZW_QCDScale_met2->AddBinContent(iBin);
      ZW_QCDShape_met2->AddBinContent(iBin);
      ZW_QCDProcess_met2->AddBinContent(iBin);
      ZW_NNLOEWK_met2->AddBinContent(iBin);
      ZW_Sudakov1_met2->AddBinContent(iBin);
      ZW_Sudakov2_met2->AddBinContent(iBin);
      ZW_NNLOMiss1_met2->AddBinContent(iBin);
      ZW_NNLOMiss2_met2->AddBinContent(iBin);
      ZW_MIX_met2->AddBinContent(iBin);
      ZW_PDF_met2->AddBinContent(iBin);
    }


    plotHistograms(canvas, ZW_QCDScale_met, ZW_QCDScale_met2, "outputTheoryComparison", "ZW_QCDScale");
    plotHistograms(canvas, ZW_QCDShape_met, ZW_QCDShape_met2, "outputTheoryComparison", "ZW_QCDShape");
    plotHistograms(canvas, ZW_QCDProcess_met, ZW_QCDProcess_met2, "outputTheoryComparison", "ZW_QCDProcess");
    plotHistograms(canvas, ZW_NNLOEWK_met, ZW_NNLOEWK_met2, "outputTheoryComparison", "ZW_NNLOEWK");
    plotHistograms(canvas, ZW_Sudakov1_met, ZW_Sudakov1_met2, "outputTheoryComparison", "ZW_Sudakov1");
    plotHistograms(canvas, ZW_Sudakov2_met, ZW_Sudakov2_met2, "outputTheoryComparison", "ZW_Sudakov2");
    plotHistograms(canvas, ZW_NNLOMiss1_met, ZW_NNLOMiss1_met2, "outputTheoryComparison", "ZW_NNLOMiss1");
    plotHistograms(canvas, ZW_NNLOMiss2_met, ZW_NNLOMiss2_met2, "outputTheoryComparison", "ZW_NNLOMiss2");
    plotHistograms(canvas, ZW_MIX_met, ZW_MIX_met2, "outputTheoryComparison", "ZW_MIX");
    plotHistograms(canvas, ZW_PDF_met, ZW_PDF_met2, "outputTheoryComparison", "ZW_PDF");
    
  }

  if(isZG){

    TH1* ZG_QCDScale_met = (TH1*) file1->Get("ZG_QCDScale_met");
    TH1* ZG_QCDShape_met = (TH1*) file1->Get("ZG_QCDShape_met");
    TH1* ZG_QCDProcess_met = (TH1*) file1->Get("ZG_QCDProcess_met");
    TH1* ZG_NNLOEWK_met = (TH1*) file1->Get("ZG_NNLOEWK_met");
    TH1* ZG_Sudakov1_met = (TH1*) file1->Get("ZG_Sudakov1_met");
    TH1* ZG_Sudakov2_met = (TH1*) file1->Get("ZG_Sudakov2_met");
    TH1* ZG_NNLOMiss1_met = (TH1*) file1->Get("ZG_NNLOMiss1_met");
    TH1* ZG_NNLOMiss2_met = (TH1*) file1->Get("ZG_NNLOMiss2_met");
    TH1* ZG_MIX_met = (TH1*) file1->Get("ZG_MIX_met");
    TH1* ZG_PDF_met = (TH1*) file1->Get("ZG_PDF_met");

    if(category == Category::monoV){
      ZG_QCDScale_met = (TH1*) file1->Get("ZG_QCDScale_met_monov");
      ZG_QCDShape_met = (TH1*) file1->Get("ZG_QCDShape_met_monov");
      ZG_QCDProcess_met = (TH1*) file1->Get("ZG_QCDProcess_met_monov");
      ZG_NNLOEWK_met = (TH1*) file1->Get("ZG_NNLOEWK_met_monov");
      ZG_Sudakov1_met = (TH1*) file1->Get("ZG_Sudakov1_met_monov");
      ZG_Sudakov2_met = (TH1*) file1->Get("ZG_Sudakov2_met_monov");
      ZG_NNLOMiss1_met = (TH1*) file1->Get("ZG_NNLOMiss1_met_monov");
      ZG_NNLOMiss2_met = (TH1*) file1->Get("ZG_NNLOMiss2_met_monov");
      ZG_MIX_met = (TH1*) file1->Get("ZG_MIX_met_monov");
      ZG_PDF_met = (TH1*) file1->Get("ZG_PDF_met_monov");
    }

    TH1* ZG_QCDScale_met2 = (TH1*) file2->FindObjectAny("ZG_QCDScale_met");
    TH1* ZG_QCDShape_met2 = (TH1*) file2->FindObjectAny("ZG_QCDShape_met");
    TH1* ZG_QCDProcess_met2 = (TH1*) file2->FindObjectAny("ZG_QCDProcess_met");
    TH1* ZG_NNLOEWK_met2 = (TH1*) file2->FindObjectAny("ZG_NNLOEWK_met");
    TH1* ZG_Sudakov1_met2 = (TH1*) file2->FindObjectAny("ZG_Sudakov1_met");
    TH1* ZG_Sudakov2_met2 = (TH1*) file2->FindObjectAny("ZG_Sudakov2_met");
    TH1* ZG_NNLOMiss1_met2 = (TH1*) file2->FindObjectAny("ZG_NNLOMiss1_met");
    TH1* ZG_NNLOMiss2_met2 = (TH1*) file2->FindObjectAny("ZG_NNLOMiss2_met");
    TH1* ZG_MIX_met2 = (TH1*) file2->FindObjectAny("ZG_MIX_met");
    TH1* ZG_PDF_met2 = (TH1*) file2->FindObjectAny("ZG_PDF_met");

    for(int iBin = 1; iBin <= ZG_QCDScale_met2->GetNbinsX(); iBin++){
      ZG_QCDScale_met2->AddBinContent(iBin);
      ZG_QCDShape_met2->AddBinContent(iBin);
      ZG_QCDProcess_met2->AddBinContent(iBin);
      ZG_NNLOEWK_met2->AddBinContent(iBin);
      ZG_Sudakov1_met2->AddBinContent(iBin);
      ZG_Sudakov2_met2->AddBinContent(iBin);
      ZG_NNLOMiss1_met2->AddBinContent(iBin);
      ZG_NNLOMiss2_met2->AddBinContent(iBin);
      ZG_MIX_met2->AddBinContent(iBin);
      ZG_PDF_met2->AddBinContent(iBin);
    }

    plotHistograms(canvas, ZG_QCDScale_met, ZG_QCDScale_met2, "outputTheoryComparison", "ZG_QCDScale");
    plotHistograms(canvas, ZG_QCDShape_met, ZG_QCDShape_met2, "outputTheoryComparison", "ZG_QCDShape");
    plotHistograms(canvas, ZG_QCDProcess_met, ZG_QCDProcess_met2, "outputTheoryComparison", "ZG_QCDProcess");
    plotHistograms(canvas, ZG_NNLOEWK_met, ZG_NNLOEWK_met2, "outputTheoryComparison", "ZG_NNLOEWK");
    plotHistograms(canvas, ZG_Sudakov1_met, ZG_Sudakov1_met2, "outputTheoryComparison", "ZG_Sudakov1");
    plotHistograms(canvas, ZG_Sudakov2_met, ZG_Sudakov2_met2, "outputTheoryComparison", "ZG_Sudakov2");
    plotHistograms(canvas, ZG_NNLOMiss1_met, ZG_NNLOMiss1_met2, "outputTheoryComparison", "ZG_NNLOMiss1");
    plotHistograms(canvas, ZG_NNLOMiss2_met, ZG_NNLOMiss2_met2, "outputTheoryComparison", "ZG_NNLOMiss2");
    plotHistograms(canvas, ZG_MIX_met, ZG_MIX_met2, "outputTheoryComparison", "ZG_MIX");
    plotHistograms(canvas, ZG_PDF_met, ZG_PDF_met2, "outputTheoryComparison", "ZG_PDF");
  }
}
