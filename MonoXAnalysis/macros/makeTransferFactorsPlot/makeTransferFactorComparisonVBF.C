#include "../CMS_lumi.h"

void plotRatioComparison(TCanvas* canvas,TH1F* qcd_1, TH1F* qcd_2, TH1F* ewk_1, TH1F* ewk_2, string postfix, string label, string outputDIR){

  canvas->cd();

  qcd_1->GetXaxis()->SetTitle("M_{jj} [GeV]");
  qcd_1->GetYaxis()->SetTitle(label.c_str());
  qcd_1->GetXaxis()->SetTitleOffset(1.2);
  qcd_1->GetYaxis()->SetTitleOffset(1.2);

  qcd_1->SetLineColor(kBlack);
  qcd_1->SetLineWidth(2);
  qcd_1->SetMarkerColor(kBlack);
  qcd_1->SetMarkerStyle(20);
  qcd_1->SetMarkerSize(1.1);
  qcd_1->Draw("hist");

  qcd_2->SetLineColor(kBlack);
  qcd_2->SetLineWidth(2);
  qcd_2->SetMarkerColor(kBlack);
  qcd_2->SetMarkerStyle(20);
  qcd_2->SetMarkerSize(1.1);
  qcd_2->Draw("EPsame");

  ewk_1->SetLineColor(kRed);
  ewk_1->SetLineWidth(2);
  ewk_1->SetMarkerColor(kRed);
  ewk_1->SetMarkerStyle(20);
  ewk_1->SetMarkerSize(1.1);
  ewk_1->Draw("hist same");

  ewk_2->SetLineColor(kRed);
  ewk_2->SetLineWidth(2);
  ewk_2->SetMarkerColor(kRed);
  ewk_2->SetMarkerStyle(20);
  ewk_2->SetMarkerSize(1.1);
  ewk_2->Draw("EPsame");

  qcd_1->GetYaxis()->SetRangeUser(min(qcd_1->GetMinimum(),min(qcd_2->GetMinimum(),min(ewk_1->GetMinimum(),ewk_2->GetMinimum())))*0.75,
				  max(qcd_1->GetMaximum(),max(qcd_2->GetMaximum(),max(ewk_1->GetMaximum(),ewk_2->GetMaximum())))*1.5);


  TLegend leg (0.6,0.65,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.AddEntry(qcd_1,"QCD, w. Pre-firing","L");
  leg.AddEntry(ewk_1,"EW, w. Pre-firing","L");
  leg.AddEntry(qcd_2,"QCD, w/o Pre-firing","PEL");
  leg.AddEntry(ewk_2,"EW, w/o Pre-firing","PEL");
  leg.Draw("same");

  CMS_lumi(canvas,"35.9");
  
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}

void makeTransferFactorComparisonVBF(string file1, string file2, string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  
  TFile* inputFile_1 = TFile::Open(file1.c_str());
  TFile* inputFile_2 = TFile::Open(file2.c_str());

  TH1F* qcd_zvv_1 = (TH1F*) inputFile_1->FindObjectAny("zinvhist_mjj");
  TH1F* qcd_zvv_2 = (TH1F*) inputFile_2->FindObjectAny("zinvhist_mjj");
  TH1F* ewk_zvv_1 = (TH1F*) inputFile_1->FindObjectAny("ewkbkgzhist_mjj");
  TH1F* ewk_zvv_2 = (TH1F*) inputFile_2->FindObjectAny("ewkbkgzhist_mjj");

  TH1F* qcd_wjet_1 = (TH1F*) inputFile_1->FindObjectAny("wjethist_mjj");
  TH1F* qcd_wjet_2 = (TH1F*) inputFile_2->FindObjectAny("wjethist_mjj");
  TH1F* ewk_wjet_1 = (TH1F*) inputFile_1->FindObjectAny("ewkbkgwhist_mjj");
  TH1F* ewk_wjet_2 = (TH1F*) inputFile_2->FindObjectAny("ewkbkgwhist_mjj");

  TH1F* qcd_zmm_1 = (TH1F*) inputFile_1->FindObjectAny("vllbkghistzmm_mjj");
  TH1F* qcd_zmm_2 = (TH1F*) inputFile_2->FindObjectAny("vllbkghistzmm_mjj");
  TH1F* ewk_zmm_1 = (TH1F*) inputFile_1->FindObjectAny("ewkzbkghistzmm_mjj");
  TH1F* ewk_zmm_2 = (TH1F*) inputFile_2->FindObjectAny("ewkzbkghistzmm_mjj");

  TH1F* qcd_zee_1 = (TH1F*) inputFile_1->FindObjectAny("vllbkghistzee_mjj");
  TH1F* qcd_zee_2 = (TH1F*) inputFile_2->FindObjectAny("vllbkghistzee_mjj");
  TH1F* ewk_zee_1 = (TH1F*) inputFile_1->FindObjectAny("ewkzbkghistzee_mjj");
  TH1F* ewk_zee_2 = (TH1F*) inputFile_2->FindObjectAny("ewkzbkghistzee_mjj");

  TH1F* qcd_wmn_1 = (TH1F*) inputFile_1->FindObjectAny("vlbkghistwmn_mjj");
  TH1F* qcd_wmn_2 = (TH1F*) inputFile_2->FindObjectAny("vlbkghistwmn_mjj");
  TH1F* ewk_wmn_1 = (TH1F*) inputFile_1->FindObjectAny("ewkwbkghistwmn_mjj");
  TH1F* ewk_wmn_2 = (TH1F*) inputFile_2->FindObjectAny("ewkwbkghistwmn_mjj");

  TH1F* qcd_wen_1 = (TH1F*) inputFile_1->FindObjectAny("vlbkghistwen_mjj");
  TH1F* qcd_wen_2 = (TH1F*) inputFile_2->FindObjectAny("vlbkghistwen_mjj");
  TH1F* ewk_wen_1 = (TH1F*) inputFile_1->FindObjectAny("ewkwbkghistwen_mjj");
  TH1F* ewk_wen_2 = (TH1F*) inputFile_2->FindObjectAny("ewkwbkghistwen_mjj");

  TH1F* qcd_zvv_1_temp = (TH1F*) qcd_zvv_1->Clone("qcd_zvv_1_temp");
  TH1F* qcd_zvv_2_temp = (TH1F*) qcd_zvv_2->Clone("qcd_zvv_2_temp");
  TH1F* ewk_zvv_1_temp = (TH1F*) ewk_zvv_1->Clone("ewk_zvv_1_temp");
  TH1F* ewk_zvv_2_temp = (TH1F*) ewk_zvv_2->Clone("ewk_zvv_2_temp");

  qcd_zvv_1_temp->Divide(qcd_wjet_1);
  qcd_zvv_2_temp->Divide(qcd_wjet_2);
  ewk_zvv_1_temp->Divide(ewk_wjet_1);
  ewk_zvv_2_temp->Divide(ewk_wjet_2);

  TH1F* qcd_zll_1 = (TH1F*) qcd_zmm_1->Clone("qcd_zll_1");
  TH1F* qcd_zll_2 = (TH1F*) qcd_zmm_2->Clone("qcd_zll_2");
  qcd_zll_1->Add(qcd_zee_1);
  qcd_zll_2->Add(qcd_zee_2);
  qcd_zvv_1->Divide(qcd_zll_1);
  qcd_zvv_2->Divide(qcd_zll_2);

  TH1F* ewk_zll_1 = (TH1F*) ewk_zmm_1->Clone("ewk_zll_1");
  TH1F* ewk_zll_2 = (TH1F*) ewk_zmm_2->Clone("ewk_zll_2");
  ewk_zll_1->Add(ewk_zee_1);
  ewk_zll_2->Add(ewk_zee_2);
  ewk_zvv_1->Divide(ewk_zll_1);
  ewk_zvv_2->Divide(ewk_zll_2);

  TH1F* qcd_wln_1 = (TH1F*) qcd_wmn_1->Clone("qcd_wln_1");
  TH1F* qcd_wln_2 = (TH1F*) qcd_wmn_2->Clone("qcd_wln_2");
  qcd_wln_1->Add(qcd_wen_1);
  qcd_wln_2->Add(qcd_wen_2);
  qcd_wjet_1->Divide(qcd_wmn_1);
  qcd_wjet_2->Divide(qcd_wmn_2);

  TH1F* ewk_wln_1 = (TH1F*) ewk_wmn_1->Clone("ewk_wln_1");
  TH1F* ewk_wln_2 = (TH1F*) ewk_wmn_2->Clone("ewk_wln_2");
  ewk_wln_1->Add(ewk_wen_1);
  ewk_wln_2->Add(ewk_wen_2);
  ewk_wjet_1->Divide(ewk_wmn_1);
  ewk_wjet_2->Divide(ewk_wmn_2);

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();

  plotRatioComparison(canvas,qcd_zvv_1_temp,qcd_zvv_2_temp,ewk_zvv_1_temp,ewk_zvv_2_temp,"zvv_wjet_ratio","Z(#nu#nu)+jets / W+jets",outputDIR);
  plotRatioComparison(canvas,qcd_zvv_1,qcd_zvv_2,ewk_zvv_1,ewk_zvv_2,"zll_znn_ratio", "Z(#nu#nu)+jets / Z(ll)+jets",outputDIR);
  plotRatioComparison(canvas,qcd_wjet_1,qcd_wjet_2,ewk_wjet_1,ewk_wjet_2,"wln_wjet_ratio"," W+jets / W(l#nu)+jets",outputDIR);

  
}
