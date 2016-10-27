#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static string lumi_ = "12.9";

void plotTransferFactor(TCanvas* canvas,TH1* histo1, TH1* histo2, const string & xAxisName, const string & yAxisName, const string & outputDIR, const string & postfix){

  histo1->GetXaxis()->SetTitle(xAxisName.c_str());
  histo1->GetYaxis()->SetTitle(yAxisName.c_str());
  histo1->GetYaxis()->SetTitleOffset(1.2);
  histo1->GetXaxis()->SetTitleOffset(1.1);
  histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetTitleSize()*1.2);
  histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetTitleSize()*1.2);
  histo1->GetYaxis()->CenterTitle();

  histo1->SetLineColor(kRed);  
  histo1->SetLineWidth(2);
  if(histo2 != NULL){
    histo2->SetLineColor(kBlue);  
    histo2->SetLineWidth(2);
  }
  TH1F* temp = (TH1F*) histo1->Clone(Form("%s_temp",histo1->GetName()));  
  temp->SetFillColor(kRed);
  temp->SetFillStyle(3001);
  temp->SetMarkerSize(0);
  
  TH1F* temp2 = NULL;
  if(histo2 != NULL){
    temp2 = (TH1F*) histo2->Clone(Form("%s_temp",histo2->GetName()));  
    temp2->SetFillColor(kBlue);
    temp2->SetFillStyle(3004);
    temp2->SetMarkerSize(0);
  }
  if(histo2 != NULL){
    histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),histo2->GetMinimum())*0.7,max(histo1->GetMaximum(),histo2->GetMaximum())*1.3);
  }
  else
    histo1->GetYaxis()->SetRangeUser(histo1->GetMinimum()*0.5, histo1->GetMaximum()*2);

  histo1->Draw("hist");
  temp->Draw("E2 same");
  histo1->Draw("hist same");
  if(histo2 != NULL){
    temp2->Draw("E2 same");
    histo1->Draw("hist same");
    histo2->Draw("hist same");
  }
  CMS_lumi(canvas,lumi_,true);
  
  canvas->SaveAs((outputDIR+"/transfer_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/transfer_"+postfix+".pdf").c_str(),"pdf");
}


void makeTransferFactorComparison(string inputDIR, string outputDIR, string observable, string postfix = "", string axisName = ""){

  gROOT->SetBatch(kTRUE);
  initializeBinning();
  setTDRStyle();
  
  system(("mkdir -p "+outputDIR).c_str());

  TFile* zmm_file_qcd = TFile::Open((inputDIR+"/zmmcor"+postfix+".root").c_str());
  TFile* zee_file_qcd = TFile::Open((inputDIR+"/zeecor"+postfix+".root").c_str());
  TFile* wmn_file_qcd = TFile::Open((inputDIR+"/wmncor"+postfix+".root").c_str());
  TFile* wen_file_qcd = TFile::Open((inputDIR+"/wencor"+postfix+".root").c_str());
  TFile* zgm_file_qcd = TFile::Open((inputDIR+"/gamcorewk"+postfix+".root").c_str());
  TFile* zwj_file_qcd = TFile::Open((inputDIR+"/zwjcorewk"+postfix+".root").c_str());

  TFile* zmm_file_ewk = TFile::Open((inputDIR+"/zewkmmcor"+postfix+".root").c_str());
  TFile* zee_file_ewk = TFile::Open((inputDIR+"/zewkeecor"+postfix+".root").c_str());
  TFile* wmn_file_ewk = TFile::Open((inputDIR+"/wewkmncor"+postfix+".root").c_str());
  TFile* wen_file_ewk = TFile::Open((inputDIR+"/wewkencor"+postfix+".root").c_str());
  TFile* zwj_file_ewk = TFile::Open((inputDIR+"/zwjewkcor"+postfix+".root").c_str());

  //Make plots
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  
  // take numerator and denominator zmm for qcd and ewk components
  TH1F* rzmm_qcd = (TH1F*) zmm_file_qcd->Get(("zmmcorhist_"+observable).c_str());
  TH1F* rzmm_ewk = (TH1F*) zmm_file_ewk->Get(("zewkmmcorhist_"+observable).c_str());
  plotTransferFactor(canvas,rzmm_qcd,rzmm_ewk,axisName,"R(Z#rightarrow #nu#nu/Z#rightarrow #mu#mu)",outputDIR,"zmm_"+observable);

  // take numerator and denominator zee for qcd and ewk components
  TH1F* rzee_qcd = (TH1F*) zee_file_qcd->Get(("zeecorhist_"+observable).c_str());
  TH1F* rzee_ewk = (TH1F*) zee_file_ewk->Get(("zewkeecorhist_"+observable).c_str());
  plotTransferFactor(canvas,rzee_qcd,rzee_ewk,axisName,"R(Z#rightarrow #nu#nu/Z#rightarrow ee)",outputDIR,"zee_"+observable);

  // take numerator and denominator wmn for qcd and ewk components
  TH1F* rwmn_qcd = (TH1F*) wmn_file_qcd->Get(("wmncorhist_"+observable).c_str());
  TH1F* rwmn_ewk = (TH1F*) wmn_file_ewk->Get(("wewkmncorhist_"+observable).c_str());
  plotTransferFactor(canvas,rwmn_qcd,rwmn_ewk,axisName,"R(W+jets/W#rightarrow #mu#nu)",outputDIR,"wmn_"+observable);

  // take numerator and denominator wln for qcd and ewk components
  TH1F* rwen_qcd = (TH1F*) wen_file_qcd->Get(("wencorhist_"+observable).c_str());
  TH1F* rwen_ewk = (TH1F*) wen_file_ewk->Get(("wewkencorhist_"+observable).c_str());
  plotTransferFactor(canvas,rwen_qcd,rwen_ewk,axisName,"R(W+jets/W#rightarrow e#nu)",outputDIR,"wen_"+observable);

  // take numerator and denominator wln for qcd and ewk components
  TH1F* rgam_qcd = (TH1F*) zgm_file_qcd->Get(("gamcorewkhist_"+observable).c_str());
  plotTransferFactor(canvas,rgam_qcd,NULL,axisName,"R(Z#rightarrow #nu#nu/#gamma)",outputDIR,"gam_"+observable);

  // take numerator and denominator wln for qcd and ewk components
  TH1F* rzwj_qcd = (TH1F*) zwj_file_qcd->Get(("zwjcorewkhist_"+observable).c_str());
  TH1F* rzwj_ewk = (TH1F*) zwj_file_ewk->Get(("zwjewkcorhist_"+observable).c_str());
  plotTransferFactor(canvas,rzwj_qcd,rzwj_ewk,axisName,"R(Z#rightarrow #nu#nu/W+jets)",outputDIR,"zwj_"+observable);

  // Make Zmm / Wmn (control regions)
  TH1F* zmm_qcd = (TH1F*) zmm_file_qcd->Get(("dhist_zmm_"+observable).c_str());
  TH1F* zmm_ewk = (TH1F*) zmm_file_ewk->Get(("dhist_ewk_zmm_"+observable).c_str());
  TH1F* wmn_qcd = (TH1F*) wmn_file_qcd->Get(("dhist_wmn_"+observable).c_str());
  TH1F* wmn_ewk = (TH1F*) wmn_file_ewk->Get(("dhist_ewk_wmn_"+observable).c_str());
  TH1F* zll_qcd = (TH1F*) zmm_qcd->Clone("zll_qcd");
  TH1F* zll_ewk = (TH1F*) zmm_ewk->Clone("zll_ewk");
  TH1F* wln_qcd = (TH1F*) wmn_qcd->Clone("wln_qcd");
  TH1F* wln_ewk = (TH1F*) wmn_ewk->Clone("wln_ewk");
  zmm_qcd->Divide(wmn_qcd);
  zmm_ewk->Divide(wmn_ewk);
  plotTransferFactor(canvas,zmm_qcd,zmm_ewk,axisName,"R(Z#rightarrow #mu#mu/W#rightarrow #mu#nu)",outputDIR,"zmmwmn_"+observable);

  // Make Zee / Wen (control regions)
  TH1F* zee_qcd = (TH1F*) zee_file_qcd->Get(("dhist_zee_"+observable).c_str());
  TH1F* zee_ewk = (TH1F*) zee_file_ewk->Get(("dhist_ewk_zee_"+observable).c_str());
  TH1F* wen_qcd = (TH1F*) wen_file_qcd->Get(("dhist_wen_"+observable).c_str());
  TH1F* wen_ewk = (TH1F*) wen_file_ewk->Get(("dhist_ewk_wen_"+observable).c_str());
  zll_qcd->Add(zee_qcd);
  zll_ewk->Add(zee_ewk);
  wln_qcd->Add(wen_qcd);
  wln_ewk->Add(wen_ewk);
  zee_qcd->Divide(wen_qcd);
  zee_ewk->Divide(wen_ewk);
  plotTransferFactor(canvas,zee_qcd,zee_ewk,axisName,"R(Z#rightarrow ee/W#rightarrow e#nu)",outputDIR,"zeewen_"+observable);

  // Make Zll / Wln (control regions)
  zll_qcd->Divide(wln_qcd);
  zll_ewk->Divide(wln_ewk);
  plotTransferFactor(canvas,zll_qcd,zll_ewk,axisName,"R(Z#rightarrow ll/W#rightarrow l#nu)",outputDIR,"zllwln_"+observable);

  // Transfer factor Zvv QCD / Zvv EWK SR 
  TH1F* zvv_qcd = (TH1F*) zmm_file_qcd->Get(("nhist_zmm_"+observable).c_str());
  TH1F* zvv_ewk = (TH1F*) zmm_file_ewk->Get(("nhist_ewk_zmm_"+observable).c_str());
  TH1F* wjet_qcd = (TH1F*) wmn_file_qcd->Get(("nhist_wmn_"+observable).c_str());
  TH1F* wjet_ewk = (TH1F*) wmn_file_ewk->Get(("nhist_ewk_wmn_"+observable).c_str());
  
  TH1F* rzvv  = (TH1F*) zvv_ewk->Clone("rzvv");
  TH1F* rwjet = (TH1F*) wjet_ewk->Clone("rwjet");
  rzvv->Divide(zvv_qcd);
  rwjet->Divide(wjet_qcd);
  
  plotTransferFactor(canvas,rzvv,rwjet,axisName,"R(V_{EWK}/V_{QCD})",outputDIR,"zoverz_"+observable);

  // Z/W LO vs NLO for QCD
  TFile* zwj_file_qcd_nlo = TFile::Open((inputDIR+"/zwjcorqcd"+postfix+".root").c_str());
  TFile* zwj_file_qcd_lo = TFile::Open((inputDIR+"/zwjcor"+postfix+".root").c_str());

  TH1F* rzwj_qcd_nlo = (TH1F*) zwj_file_qcd_nlo->Get(("zwjcorqcdhist_"+observable).c_str());
  TH1F* rzwj_qcd_lo = (TH1F*) zwj_file_qcd_lo->Get(("zwjcorhist_"+observable).c_str());

  canvas->cd();
  rzwj_qcd_nlo->SetFillColor(0);
  rzwj_qcd_nlo->SetFillStyle(0);
  rzwj_qcd_nlo->SetLineColor(kBlack);
  rzwj_qcd_nlo->SetLineWidth(2);
  rzwj_qcd_lo->SetLineColor(kBlue);
  rzwj_qcd_lo->SetLineWidth(2);
  rzwj_qcd_nlo->GetYaxis()->SetTitle("Z/W Ratio");
  rzwj_qcd->Draw("hist");
  rzwj_qcd_nlo->Draw("hist same");
  rzwj_qcd_lo->Draw("hist same");

  TLegend leg (0.2,0.6,0.4,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(rzwj_qcd,"Z/W NLO QCD x EWK","L");
  leg.AddEntry(rzwj_qcd_nlo,"Z/W NLO QCD","L");
  leg.AddEntry(rzwj_qcd_lo,"Z/W LO QCD","L");
  leg.Draw("same");

  CMS_lumi(canvas,lumi_,"true");

  canvas->SaveAs((outputDIR+"/zwjcomparison.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/zwjcomparison.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/zwjcomparison.root").c_str(),"root");
  
}
