#include "../CMS_lumi.h"

void plotHist(TH1F* mc1, TH1F* mc2,TCanvas* canvas, int ratioType, string observable){

  canvas->cd();
  TH1F* ratio_mc = (TH1F*) mc1->Clone("ratio_mc");
  ratio_mc->Divide(mc2);
  ratio_mc->SetLineColor(kRed);
  ratio_mc->SetLineWidth(2);

  TH1* histoMCband = (TH1*) ratio_mc->Clone("histoMCband");
  histoMCband->SetFillColor(kGray);
  histoMCband->GetXaxis()->SetTitle(observable.c_str());
  histoMCband->GetYaxis()->SetTitle("Ratio");
  histoMCband->GetYaxis()->SetRangeUser(0.85,1.15);
  histoMCband->Draw("E2");
  ratio_mc->Draw("hist same");

  CMS_lumi(canvas,"2.6",true);
  
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(ratioType == 1){
    leg->AddEntry(ratio_mc,"Z/#gamma MC","FL");
  }
  else if(ratioType == 2){
    leg->AddEntry(ratio_mc,"W/#gamma MC","FL");
  }
  else if(ratioType == 3){
    leg->AddEntry(ratio_mc,"Z/W MC","FL");
  }
  leg->Draw("same");
  
  if(ratioType == 1){
    canvas->SaveAs(("Zgamma_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("Zgamma_"+observable+".pdf").c_str(),"pdf");
  }
  else if(ratioType == 2){
    canvas->SaveAs(("Wgamma_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("Wgamma_"+observable+".pdf").c_str(),"pdf");
  }
  else if(ratioType == 3){
    canvas->SaveAs(("ZW_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("ZW_"+observable+".pdf").c_str(),"pdf");
  }

}


void plotHist( TH1F* data1, TH1F* data2, TH1F* mc1, TH1F* mc2,TCanvas* canvas, int ratioType, string observable){

  canvas->cd();
  TH1F* ratio_data = (TH1F*) data1->Clone("ratio_data");
  ratio_data->Divide(data2);
  TH1F* ratio_mc = (TH1F*) mc1->Clone("ratio_mc");
  ratio_mc->Divide(mc2);

  ratio_data->GetXaxis()->SetTitle(observable.c_str());
  ratio_data->GetYaxis()->SetTitle("Ratio");
  ratio_data->SetMarkerColor(kBlack);
  ratio_data->SetLineColor(kBlack);
  ratio_data->SetMarkerStyle(20);
  ratio_data->SetMarkerSize(1);
  ratio_data->GetYaxis()->SetRangeUser(0.65,1.35);
  ratio_data->Draw("EP");
  ratio_mc->SetLineColor(kRed);
  ratio_mc->SetLineWidth(2);

  TH1* histoMCband = (TH1*) ratio_mc->Clone("histoMCband");
  histoMCband->SetFillColor(kGray);
  histoMCband->Draw("E2 same");
  ratio_mc->Draw("hist same");
  ratio_data->Draw("EP same");

  CMS_lumi(canvas,"2.6",true);
  
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(ratioType == 1){
    leg->AddEntry(ratio_data,"Z/#gamma Data","PE");
    leg->AddEntry(ratio_mc,"Z/#gamma MC","FL");
  }
  else if(ratioType == 2){
    leg->AddEntry(ratio_data,"W/#gamma Data","PE");
    leg->AddEntry(ratio_mc,"W/#gamma MC","FL");
  }
  else if(ratioType == 3){
    leg->AddEntry(ratio_data,"Z/W Data","PE");
    leg->AddEntry(ratio_mc,"Z/W MC","FL");
  }
  leg->Draw("same");
  
  if(ratioType == 1){
    canvas->SaveAs(("Zgamma_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("Zgamma_"+observable+".pdf").c_str(),"pdf");
  }
  else if(ratioType == 2){
    canvas->SaveAs(("Wgamma_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("Wgamma_"+observable+".pdf").c_str(),"pdf");
  }
  else if(ratioType == 3){
    canvas->SaveAs(("ZW_"+observable+".png").c_str(),"png");
    canvas->SaveAs(("ZW_"+observable+".pdf").c_str(),"pdf");
  }

}


void makeQGLComparison(string templateFile, string observable, bool useData){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TFile* inputTemplate = TFile::Open(templateFile.c_str());

  // ZMM control region
  TH1F* vllzmm_mc  = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
  TH1F* vlzmm_mc   = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzmm_"+observable).c_str());
  TH1F* topzmm_mc  = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzmm_"+observable).c_str());
  TH1F* dbzmm_mc   = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzmm_"+observable).c_str());
  TH1F* datazmm    = (TH1F*) inputTemplate->FindObjectAny(("datahistzmm_"+observable).c_str());
  datazmm->Add(vlzmm_mc,-1);
  datazmm->Add(topzmm_mc,-1);
  datazmm->Add(dbzmm_mc,-1);

  // ZEE control region
  TH1F* vllzee_mc  = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistzee_"+observable).c_str());
  TH1F* vlzee_mc   = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzee_"+observable).c_str());
  TH1F* topzee_mc  = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzee_"+observable).c_str());
  TH1F* dbzee_mc   = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzee_"+observable).c_str());
  TH1F* datazee    = (TH1F*) inputTemplate->FindObjectAny(("datahistzee_"+observable).c_str());
  datazee->Add(vlzee_mc,-1);
  datazee->Add(topzee_mc,-1);
  datazee->Add(dbzee_mc,-1);

  // Add Zmm + Zee
  TH1F* datazll = (TH1F*) datazmm->Clone("datazll");
  datazll->Add(datazee);
  TH1F* vllzll_mc = (TH1F*) vllzmm_mc->Clone("vllzll_mc");
  vllzll_mc->Add(vllzee_mc);

  // WMN control region
  TH1F* vllwmn_mc  = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwmn_"+observable).c_str());
  TH1F* vlwmn_mc   = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
  TH1F* topwmn_mc  = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwmn_"+observable).c_str());
  TH1F* dbwmn_mc   = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwmn_"+observable).c_str());
  TH1F* qcdwmn_mc  = (TH1F*) inputTemplate->FindObjectAny(("qbkghistwmn_"+observable).c_str());
  TH1F* datawmn    = (TH1F*) inputTemplate->FindObjectAny(("datahistwmn_"+observable).c_str());
  datawmn->Add(vllwmn_mc,-1);
  datawmn->Add(topwmn_mc,-1);
  datawmn->Add(dbwmn_mc,-1);
  datawmn->Add(qcdwmn_mc,-1);

  // WEN control region
  TH1F* vllwen_mc  = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwen_"+observable).c_str());
  TH1F* vlwen_mc   = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistwen_"+observable).c_str());
  TH1F* topwen_mc  = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwen_"+observable).c_str());
  TH1F* dbwen_mc   = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwen_"+observable).c_str());
  TH1F* qcdwen_mc  = (TH1F*) inputTemplate->FindObjectAny(("qbkghistwen_"+observable).c_str());
  TH1F* datawen    = (TH1F*) inputTemplate->FindObjectAny(("datahistwen_"+observable).c_str());
  datawen->Add(vllwen_mc,-1);
  datawen->Add(topwen_mc,-1);
  datawen->Add(dbwen_mc,-1);
  datawen->Add(qcdwen_mc,-1);

  // Add Wmn + Wen
  TH1F* datawln = (TH1F*) datawmn->Clone("datawln");
  datawln->Add(datawen);
  TH1F* vlwln_mc = (TH1F*) vlwmn_mc->Clone("vlwln_mc");
  vlwln_mc->Add(vlwln_mc);
  
  // Gam control region
  TH1F* gjetsgam_mc = (TH1F*) inputTemplate->FindObjectAny(("gbkghistgam_"+observable).c_str());
  TH1F* qcdgam_mc   = (TH1F*) inputTemplate->FindObjectAny(("qbkghistgam_"+observable).c_str());
  TH1F* datagam     = (TH1F*) inputTemplate->FindObjectAny(("datahistgam_"+observable).c_str());
  datagam->Add(qcdgam_mc,-1);
    
  TCanvas* canvas = new TCanvas("canvas","",600,625);

  datazll->Scale(1./datazll->Integral());
  datagam->Scale(1./datagam->Integral());
  datawln->Scale(1./datawln->Integral());
  vllzll_mc->Scale(1./vllzll_mc->Integral());
  gjetsgam_mc->Scale(1./gjetsgam_mc->Integral());
  vlwln_mc->Scale(1./vlwln_mc->Integral());

  if(useData){
    //Plot Z/gamma QGL 
    plotHist(datazll,datagam,vllzll_mc,gjetsgam_mc,canvas,1,observable);
    //Plot W/gamma QGL 
    plotHist(datawln,datagam,vlwln_mc,gjetsgam_mc,canvas,2,observable);
    //Plot Z/W QGL
    plotHist(datazll,datawln,vllzll_mc,vlwln_mc,canvas,3,observable);
  }
  else{
    //Plot Z/gamma QGL 
    plotHist(vllzll_mc,gjetsgam_mc,canvas,1,observable);
    //Plot W/gamma QGL 
    plotHist(vlwln_mc,gjetsgam_mc,canvas,2,observable);
    //Plot Z/W QGL
    plotHist(vllzll_mc,vlwln_mc,canvas,3,observable);
  }
}
