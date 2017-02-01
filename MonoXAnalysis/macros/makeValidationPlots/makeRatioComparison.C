#include "../makeTemplates/histoUtils.h"
#include "../makeTemplates/histoUtils2D.h"
#include "../CMS_lumi.h"

void makePlot(TH1* histoData, TH1* histoMC, const string & observable, const Category & category, const string & observableLatex, const string & postfix, const int & rebinFactor){

  // final plot                                                                                                                                                                                     
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,650);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  TPad *pad1 = new TPad(("pad1"+postfix).c_str(),"",0,0.30,1,1);
  pad1->SetTickx();
  pad1->SetTicky();
  pad1->SetBottomMargin(0.02);

  TPad *pad2 = new TPad(("pad2"+postfix).c_str(),"",0,0.,1,0.30);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->SetTopMargin(0.08);
  pad2->SetBottomMargin(0.3);

  // Draw Pad1                                                                                                                                                                                       
  pad1->Draw();
  canvas->cd();
  pad2->Draw();
  canvas->cd();

  if(rebinFactor > 1){
    histoData->Rebin(rebinFactor);
    histoMC->Rebin(rebinFactor);
  }
  pad1->cd();

  TH1* frame  = pad1->DrawFrame(histoData->GetXaxis()->GetBinLowEdge(1),0.,histoData->GetXaxis()->GetBinLowEdge(histoData->GetNbinsX()+1),0.2, "");
  frame->GetXaxis()->SetTitle(observableLatex.c_str());
  if(TString(postfix).Contains("ZG"))
    frame->GetYaxis()->SetTitle("Ratio Z/#gamma");
  else if(TString(postfix).Contains("ZW"))
    frame->GetYaxis()->SetTitle("Ratio Z/W");
  else if(TString(postfix).Contains("WG"))
    frame->GetYaxis()->SetTitle("Ratio W/#gamma");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->Draw();
  CMS_lumi(pad1,"36.4",true);

  float maxdata  = -1;
  float mindata  = 1000000;
  for(int iBin = 0; iBin < histoData->GetNbinsX(); iBin++){
    if(histoData->GetBinContent(iBin+1)+histoData->GetBinError(iBin+1) >= maxdata)
      maxdata = histoData->GetBinContent(iBin+1)+histoData->GetBinError(iBin+1);
    if(histoData->GetBinContent(iBin+1)-histoData->GetBinError(iBin+1) <= mindata)
      mindata = histoData->GetBinContent(iBin+1)-histoData->GetBinError(iBin+1);
  }
  float maxmc  = -1;
  float minmc  = 1000000;
  for(int iBin = 0; iBin < histoMC->GetNbinsX(); iBin++){
    if(histoMC->GetBinContent(iBin+1)+histoMC->GetBinError(iBin+1) >= maxmc)
      maxmc = histoMC->GetBinContent(iBin+1)+histoMC->GetBinError(iBin+1);
    if(histoMC->GetBinContent(iBin+1)-histoMC->GetBinError(iBin+1) <= maxmc)
      minmc = histoMC->GetBinContent(iBin+1)-histoMC->GetBinError(iBin+1);
  }

  frame->GetYaxis()->SetRangeUser(min(mindata,minmc)*0.5,max(maxdata,maxmc)*1.5);

  // histo style                                                                                                                                                                                      
  histoData->SetLineColor(kBlack);
  histoData->SetLineWidth(2);
  histoData->SetMarkerColor(kBlack);
  histoData->SetMarkerStyle(20);
  histoData->SetMarkerSize(1.);

  histoMC->SetLineColor(kRed);
  histoMC->SetLineWidth(2);
  histoMC->SetMarkerSize(0);

  TLegend leg (0.5,0.5,0.7,0.7);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(histoData,"Moriond analysis","L");
  leg.AddEntry(histoMC,"ICHEP analysis","L");
  leg.Draw("same");

  histoMC->Draw("HIST same");
  histoData->Draw("HIST same");  
  pad1->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->cd();

  TH1* frame2 = NULL;
  if(category == Category::monojet)
    frame2 = pad2->DrawFrame(histoData->GetXaxis()->GetBinLowEdge(1), 0.75, histoData->GetXaxis()->GetBinLowEdge(histoData->GetNbinsX()+1), 1.25, "");
  else if(category == Category::monoV)
    frame2 = pad2->DrawFrame(histoData->GetXaxis()->GetBinLowEdge(1), 0.75, histoData->GetXaxis()->GetBinLowEdge(histoData->GetNbinsX()+1), 1.25, "");
  else if(category == Category::VBF)
    frame2 = pad2->DrawFrame(histoData->GetXaxis()->GetBinLowEdge(1), 0.75, histoData->GetXaxis()->GetBinLowEdge(histoData->GetNbinsX()+1), 1.25, "");

  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleSize(0.13);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleSize(0.10);
  frame2->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->SetTitleOffset(0.5);
  frame2->Draw();

  TH1* ratio         = (TH1*) histoData->Clone(("ratio_"+postfix).c_str());
  TH1* ratiod        = (TH1*) histoMC->Clone(("ratiod_"+postfix).c_str());
  
  for (int i = 1; i <= ratiod->GetNbinsX(); i++) ratiod->SetBinError(i, 0);

  ratio->Divide(ratiod);
  ratio->Draw("HIST same");

  canvas->SaveAs((postfix+".png").c_str(),"png");
  canvas->SaveAs((postfix+".pdf").c_str(),"pdf");
}




void makeRatioComparison(string file_1, string file_2, string observable, Category category, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile_1 = TFile::Open(file_1.c_str(),"OPEN");
  TFile* inputFile_2 = TFile::Open(file_2.c_str(),"OPEN");

  TH1* data_zmm_1 = (TH1*) inputFile_1->FindObjectAny(("datahistzmm_"+observable).c_str());
  TH1* data_zee_1 = (TH1*) inputFile_1->FindObjectAny(("datahistzee_"+observable).c_str());
  TH1* data_wmn_1 = (TH1*) inputFile_1->FindObjectAny(("datahistwmn_"+observable).c_str());
  TH1* data_wen_1 = (TH1*) inputFile_1->FindObjectAny(("datahistwen_"+observable).c_str());
  TH1* data_gam_1 = (TH1*) inputFile_1->FindObjectAny(("datahistgam_"+observable).c_str());

  TH1* mc_zmm_1 = (TH1*) inputFile_1->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
  TH1* mc_zee_1 = (TH1*) inputFile_1->FindObjectAny(("vllbkghistzee_"+observable).c_str());
  TH1* mc_wmn_1 = (TH1*) inputFile_1->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
  TH1* mc_wen_1 = (TH1*) inputFile_1->FindObjectAny(("vlbkghistwen_"+observable).c_str());
  TH1* mc_gam_1 = (TH1*) inputFile_1->FindObjectAny(("gbkghistgam_"+observable).c_str());

  mc_zmm_1->Add((TH1*) inputFile_1->FindObjectAny(("vlbkghistzmm_"+observable).c_str()));
  mc_zmm_1->Add((TH1*) inputFile_1->FindObjectAny(("dbkghistzmm_"+observable).c_str()));
  mc_zmm_1->Add((TH1*) inputFile_1->FindObjectAny(("tbkghistzmm_"+observable).c_str()));
  mc_zmm_1->Add((TH1*) inputFile_1->FindObjectAny(("qbkghistzmm_"+observable).c_str()));
  mc_zmm_1->Add((TH1*) inputFile_1->FindObjectAny(("gbkghistzmm_"+observable).c_str()));

  mc_zee_1->Add((TH1*) inputFile_1->FindObjectAny(("vlbkghistzee_"+observable).c_str()));
  mc_zee_1->Add((TH1*) inputFile_1->FindObjectAny(("dbkghistzee_"+observable).c_str()));
  mc_zee_1->Add((TH1*) inputFile_1->FindObjectAny(("tbkghistzee_"+observable).c_str()));
  mc_zee_1->Add((TH1*) inputFile_1->FindObjectAny(("qbkghistzee_"+observable).c_str()));
  mc_zee_1->Add((TH1*) inputFile_1->FindObjectAny(("gbkghistzee_"+observable).c_str()));

  mc_wmn_1->Add((TH1*) inputFile_1->FindObjectAny(("vllbkghistwmn_"+observable).c_str()));
  mc_wmn_1->Add((TH1*) inputFile_1->FindObjectAny(("dbkghistwmn_"+observable).c_str()));
  mc_wmn_1->Add((TH1*) inputFile_1->FindObjectAny(("tbkghistwmn_"+observable).c_str()));
  mc_wmn_1->Add((TH1*) inputFile_1->FindObjectAny(("qbkghistwmn_"+observable).c_str()));
  mc_wmn_1->Add((TH1*) inputFile_1->FindObjectAny(("gbkghistwmn_"+observable).c_str()));

  mc_wen_1->Add((TH1*) inputFile_1->FindObjectAny(("vllbkghistwen_"+observable).c_str()));
  mc_wen_1->Add((TH1*) inputFile_1->FindObjectAny(("dbkghistwen_"+observable).c_str()));
  mc_wen_1->Add((TH1*) inputFile_1->FindObjectAny(("tbkghistwen_"+observable).c_str()));
  mc_wen_1->Add((TH1*) inputFile_1->FindObjectAny(("qbkghistwen_"+observable).c_str()));
  mc_wen_1->Add((TH1*) inputFile_1->FindObjectAny(("gbkghistwen_"+observable).c_str()));
  
  mc_gam_1->Add((TH1*) inputFile_1->FindObjectAny(("qbkghistgam_"+observable).c_str()));
  mc_gam_1->Add((TH1*) inputFile_1->FindObjectAny(("vgbkghistgam_"+observable).c_str()));
  mc_gam_1->Add((TH1*) inputFile_1->FindObjectAny(("vlbkghistgam_"+observable).c_str()));

  // take ratios
  TH1* zgamma_data_zmm_1 = (TH1*) data_zmm_1->Clone("zgamma_data_zmm_1");
  zgamma_data_zmm_1->Divide(data_gam_1);
  TH1* zgamma_data_zee_1 = (TH1*) data_zee_1->Clone("zgamma_data_zee_1");
  zgamma_data_zee_1->Divide(data_gam_1);
  TH1* zgamma_mc_zmm_1 = (TH1*) mc_zmm_1->Clone("zgamma_mc_zmm_1");
  zgamma_mc_zmm_1->Divide(mc_gam_1);
  TH1* zgamma_mc_zee_1 = (TH1*) mc_zee_1->Clone("zgamma_mc_zee_1");
  zgamma_mc_zee_1->Divide(mc_gam_1);

  TH1* wgamma_data_wmn_1 = (TH1*) data_wmn_1->Clone("wgamma_data_wmn_1");
  wgamma_data_wmn_1->Divide(data_gam_1);
  TH1* wgamma_data_wen_1 = (TH1*) data_wen_1->Clone("wgamma_data_wen_1");
  wgamma_data_wen_1->Divide(data_gam_1);
  TH1* wgamma_mc_wmn_1 = (TH1*) mc_wmn_1->Clone("wgamma_mc_wmn_1");
  wgamma_mc_wmn_1->Divide(mc_gam_1);
  TH1* wgamma_mc_wen_1 = (TH1*) mc_wen_1->Clone("wgamma_mc_wen_1");
  wgamma_mc_wen_1->Divide(mc_gam_1);


  TH1* data_zmm_2 = (TH1*) inputFile_2->FindObjectAny(("datahistzmm_"+observable).c_str());
  TH1* data_zee_2 = (TH1*) inputFile_2->FindObjectAny(("datahistzee_"+observable).c_str());
  TH1* data_wmn_2 = (TH1*) inputFile_2->FindObjectAny(("datahistwmn_"+observable).c_str());
  TH1* data_wen_2 = (TH1*) inputFile_2->FindObjectAny(("datahistwen_"+observable).c_str());
  TH1* data_gam_2 = (TH1*) inputFile_2->FindObjectAny(("datahistgam_"+observable).c_str());

  TH1* mc_zmm_2 = (TH1*) inputFile_2->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
  TH1* mc_zee_2 = (TH1*) inputFile_2->FindObjectAny(("vllbkghistzee_"+observable).c_str());
  TH1* mc_wmn_2 = (TH1*) inputFile_2->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
  TH1* mc_wen_2 = (TH1*) inputFile_2->FindObjectAny(("vlbkghistwen_"+observable).c_str());
  TH1* mc_gam_2 = (TH1*) inputFile_2->FindObjectAny(("gbkghistgam_"+observable).c_str());

  mc_zmm_2->Add((TH1*) inputFile_2->FindObjectAny(("vlbkghistzmm_"+observable).c_str()));
  mc_zmm_2->Add((TH1*) inputFile_2->FindObjectAny(("dbkghistzmm_"+observable).c_str()));
  mc_zmm_2->Add((TH1*) inputFile_2->FindObjectAny(("tbkghistzmm_"+observable).c_str()));
  mc_zmm_2->Add((TH1*) inputFile_2->FindObjectAny(("qbkghistzmm_"+observable).c_str()));
  mc_zmm_2->Add((TH1*) inputFile_2->FindObjectAny(("gbkghistzmm_"+observable).c_str()));

  mc_zee_2->Add((TH1*) inputFile_2->FindObjectAny(("vlbkghistzee_"+observable).c_str()));
  mc_zee_2->Add((TH1*) inputFile_2->FindObjectAny(("dbkghistzee_"+observable).c_str()));
  mc_zee_2->Add((TH1*) inputFile_2->FindObjectAny(("tbkghistzee_"+observable).c_str()));
  mc_zee_2->Add((TH1*) inputFile_2->FindObjectAny(("qbkghistzee_"+observable).c_str()));
  mc_zee_2->Add((TH1*) inputFile_2->FindObjectAny(("gbkghistzee_"+observable).c_str()));

  mc_wmn_2->Add((TH1*) inputFile_2->FindObjectAny(("vllbkghistwmn_"+observable).c_str()));
  mc_wmn_2->Add((TH1*) inputFile_2->FindObjectAny(("dbkghistwmn_"+observable).c_str()));
  mc_wmn_2->Add((TH1*) inputFile_2->FindObjectAny(("tbkghistwmn_"+observable).c_str()));
  mc_wmn_2->Add((TH1*) inputFile_2->FindObjectAny(("qbkghistwmn_"+observable).c_str()));
  mc_wmn_2->Add((TH1*) inputFile_2->FindObjectAny(("gbkghistwmn_"+observable).c_str()));

  mc_wen_2->Add((TH1*) inputFile_2->FindObjectAny(("vllbkghistwen_"+observable).c_str()));
  mc_wen_2->Add((TH1*) inputFile_2->FindObjectAny(("dbkghistwen_"+observable).c_str()));
  mc_wen_2->Add((TH1*) inputFile_2->FindObjectAny(("tbkghistwen_"+observable).c_str()));
  mc_wen_2->Add((TH1*) inputFile_2->FindObjectAny(("qbkghistwen_"+observable).c_str()));
  mc_wen_2->Add((TH1*) inputFile_2->FindObjectAny(("gbkghistwen_"+observable).c_str()));

  mc_gam_2->Add((TH1*) inputFile_2->FindObjectAny(("qbkghistgam_"+observable).c_str()));
  mc_gam_2->Add((TH1*) inputFile_2->FindObjectAny(("vgbkghistgam_"+observable).c_str()));
  mc_gam_2->Add((TH1*) inputFile_2->FindObjectAny(("vlbkghistgam_"+observable).c_str()));

  // take ratios
  TH1* zgamma_data_zmm_2 = (TH1*) data_zmm_2->Clone("zgamma_data_zmm_2");
  zgamma_data_zmm_2->Divide(data_gam_2);
  TH1* zgamma_mc_zmm_2 = (TH1*) mc_zmm_2->Clone("zgamma_mc_zmm_2");
  zgamma_mc_zmm_2->Divide(mc_gam_2);

  makePlot(zgamma_data_zmm_1,zgamma_data_zmm_2,observable,category,"Recoil [GeV]","ZG_zmm_data",2);
  makePlot(zgamma_mc_zmm_1,zgamma_mc_zmm_2,observable,category,"Recoil [GeV]","ZG_zmm_mc",2);

  // take ratios
  TH1* zgamma_data_zee_2 = (TH1*) data_zee_2->Clone("zgamma_data_zee_2");
  zgamma_data_zee_2->Divide(data_gam_2);
  TH1* zgamma_mc_zee_2 = (TH1*) mc_zee_2->Clone("zgamma_mc_zee_2");
  zgamma_mc_zee_2->Divide(mc_gam_2);

  makePlot(zgamma_data_zee_1,zgamma_data_zee_2,observable,category,"Recoil [GeV]","ZG_zee_data",2);
  makePlot(zgamma_mc_zee_1,zgamma_mc_zee_2,observable,category,"Recoil [GeV]","ZG_zee_mc",2);
  

  // take ratios
  TH1* wgamma_data_wmn_2 = (TH1*) data_wmn_2->Clone("wgamma_data_wmn_2");
  wgamma_data_wmn_2->Divide(data_gam_2);
  TH1* wgamma_mc_wmn_2 = (TH1*) mc_wmn_2->Clone("wgamma_mc_wmn_2");
  wgamma_mc_wmn_2->Divide(mc_gam_2);

  makePlot(wgamma_data_wmn_1,wgamma_data_wmn_2,observable,category,"Recoil [GeV]","WG_wmn_data",2);
  makePlot(wgamma_mc_wmn_1,wgamma_mc_wmn_2,observable,category,"Recoil [GeV]","WG_wmn_mc",2);

  // take ratios
  TH1* wgamma_data_wen_2 = (TH1*) data_wen_2->Clone("wgamma_data_wen_2");
  wgamma_data_wen_2->Divide(data_gam_2);
  TH1* wgamma_mc_wen_2 = (TH1*) mc_wen_2->Clone("wgamma_mc_wen_2");
  wgamma_mc_wen_2->Divide(mc_gam_2);

  makePlot(wgamma_data_wen_1,wgamma_data_wen_2,observable,category,"Recoil [GeV]","WG_wen_data",2);
  makePlot(wgamma_mc_wen_1,wgamma_mc_wen_2,observable,category,"Recoil [GeV]","WG_wen_mc",2);


}
