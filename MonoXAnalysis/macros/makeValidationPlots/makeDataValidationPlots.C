#include "../makeTemplates/histoUtils.h"
#include "../makeTemplates/histoUtils2D.h"
#include "../CMS_lumi.h"

float musf = 0.02;
float elsf = 0.02;
float phsf = 0.02;
float mutrack = 0.01;
float eltrack = 0.01;
float mettrig = 0.01;
float eltrig = 0.02;
float phtrig = 0.02;
float lepveto = 0.03;
float inflateWZ_ewk = sqrt(2);


void makePlot(TH1* histoData, TH1* histoMC,const string & observable, const Category & category, const string & observableLatex, const string & postfix, const int & rebinFactor){

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

  if(rebinFactor > 1){
    histoData->Rebin(rebinFactor);
    histoMC->Rebin(rebinFactor);
  }

  pad1->cd();
  vector<double> bins = selectBinning(observable,category);

  TH1* frame  = pad1->DrawFrame(bins.front(),0.,bins.back(),0.2, "");
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

  TH1* histoMCband = (TH1*) histoMC->Clone("histoMCband");
  histoMCband->SetFillColor(kGray);
  histoMCband->Draw("E2same");
  histoMC->Draw("HIST same");
  histoData->Draw("PESAME");

  TLegend* leg = new TLegend(0.18, 0.66, 0.45, 0.92);
  if(TString(postfix).Contains("ZG_mm")){
    leg->AddEntry(histoData,"Z(#mu#mu)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"Z(#mu#mu)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("ZG_ee")){
    leg->AddEntry(histoData,"Z(ee)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"Z(ee)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("ZG_ll")){
    leg->AddEntry(histoData,"Z(ll)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"Z(ll)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("ZW_mm")){
    leg->AddEntry(histoData,"Z(#mu#mu)/W(#mu#nu) Data","PL");    
    leg->AddEntry(histoMCband,"Z(#mu#mu)/W(#mu#nu) MC","FL");
  }
  else if(TString(postfix).Contains("ZW_ee")){
    leg->AddEntry(histoMCband,"Z(ee)/W(e#nu) MC","FL");
    leg->AddEntry(histoData,"Z(ee)/W(e#nu) Data","PL");    
    leg->AddEntry(histoMCband,"Z(ee)/W(e#nu) MC","FL");
  }
  else if(TString(postfix).Contains("ZW_ll")){
    leg->AddEntry(histoData,"Z(ll)/W(l#nu) Data","PL");    
    leg->AddEntry(histoMCband,"Z(ll)/W(l#nu) MC","FL");
  }
  else if(TString(postfix).Contains("WG_m")){
    leg->AddEntry(histoData,"W(#mu#nu)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"W(#mu#nu)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("WG_e")){
    leg->AddEntry(histoData,"W(e#nu)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"W(e#nu)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("WG_l")){
    leg->AddEntry(histoData,"W(l#nu)/#gamma Data","PL");    
    leg->AddEntry(histoMCband,"W(l#nu)/#gamma MC","FL");
  }


  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->cd();

  TH1* frame2 = NULL;
  if(category == Category::monojet)
    frame2 = pad2->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");
  else if(category == Category::monoV)
    frame2 = pad2->DrawFrame(bins.front(), 0.0, bins.back(), 2.0, "");
  else if(category == Category::VBF)
    frame2 = pad2->DrawFrame(bins.front(), 0.0, bins.back(), 2.0, "");
  

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
  TH1* ratioD        = (TH1*) histoMC->Clone(("ratioD_"+postfix).c_str());
  TH1* ratioD_band   = (TH1*) histoMCband->Clone(("ratioD_band_"+postfix).c_str());

  for (int i = 1; i <= ratiod->GetNbinsX(); i++) ratiod->SetBinError(i, 0);
  
  ratio->Divide(ratiod);
  ratioD->Divide(ratiod);
  ratioD_band->Divide(ratiod);

  ratioD_band->Draw("E2 same");
  ratioD->Draw("HIST same");
  ratio->Draw("PE same");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((postfix+".png").c_str(),"png");
  canvas->SaveAs((postfix+".pdf").c_str(),"pdf");

}


void makeDataValidationPlots(string inputFileName, Category category, string observable, string observableLatex, bool addWgamma = false, int rebinFactor = 1){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  initializeBinning();

  // open the input file with all the templates
  TFile* inputFile = TFile::Open(inputFileName.c_str());

  TH1* data_zmm = (TH1*) inputFile->FindObjectAny(("datahistzmm_"+observable).c_str());
  TH1* data_zee = (TH1*) inputFile->FindObjectAny(("datahistzee_"+observable).c_str());
  TH1* data_wen = (TH1*) inputFile->FindObjectAny(("datahistwen_"+observable).c_str());
  TH1* data_wmn = (TH1*) inputFile->FindObjectAny(("datahistwmn_"+observable).c_str());
  TH1* data_gam = (TH1*) inputFile->FindObjectAny(("datahistgam_"+observable).c_str());
    
  // ZMM control region  
  TH1* vllbkg_zmm = (TH1*) inputFile->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
  TH1* vlbkg_zmm = (TH1*) inputFile->FindObjectAny(("vlbkghistzmm_"+observable).c_str());
  TH1* dbbkg_zmm = (TH1*) inputFile->FindObjectAny(("dbkghistzmm_"+observable).c_str());
  TH1* ttbkg_zmm = (TH1*) inputFile->FindObjectAny(("tbkghistzmm_"+observable).c_str());
  TH1* ewkwbkg_zmm = (TH1*) inputFile->FindObjectAny(("ewkwbkghistzmm_"+observable).c_str());
  TH1* ewkzbkg_zmm = (TH1*) inputFile->FindObjectAny(("ewkzbkghistzmm_"+observable).c_str());
  TH1* gambkg_zmm = (TH1*) inputFile->FindObjectAny(("gbkghistzmm_"+observable).c_str());
  TH1* qcdbkg_zmm = (TH1*) inputFile->FindObjectAny(("qbkghistzmm_"+observable).c_str());
  vllbkg_zmm->Add(vlbkg_zmm);
  vllbkg_zmm->Add(dbbkg_zmm);
  vllbkg_zmm->Add(ttbkg_zmm);
  if(category == Category::VBF){
    vllbkg_zmm->Add(ewkwbkg_zmm);
    vllbkg_zmm->Add(ewkzbkg_zmm);
  }
  vllbkg_zmm->Add(gambkg_zmm);
  vllbkg_zmm->Add(qcdbkg_zmm);

  // ZEE  control region
  TH1* vllbkg_zee = (TH1*) inputFile->FindObjectAny(("vllbkghistzee_"+observable).c_str());
  TH1* vlbkg_zee = (TH1*) inputFile->FindObjectAny(("vlbkghistzee_"+observable).c_str());
  TH1* dbbkg_zee = (TH1*) inputFile->FindObjectAny(("dbkghistzee_"+observable).c_str());
  TH1* ttbkg_zee = (TH1*) inputFile->FindObjectAny(("tbkghistzee_"+observable).c_str());
  TH1* ewkwbkg_zee = (TH1*) inputFile->FindObjectAny(("ewkwbkghistzee_"+observable).c_str());
  TH1* ewkzbkg_zee = (TH1*) inputFile->FindObjectAny(("ewkzbkghistzee_"+observable).c_str());
  TH1* gambkg_zee = (TH1*) inputFile->FindObjectAny(("gbkghistzee_"+observable).c_str());
  TH1* qcdbkg_zee = (TH1*) inputFile->FindObjectAny(("qbkghistzee_"+observable).c_str());
  vllbkg_zee->Add(vlbkg_zee);
  vllbkg_zee->Add(dbbkg_zee);
  vllbkg_zee->Add(ttbkg_zee);
  if(category == Category::VBF){
    vllbkg_zee->Add(ewkwbkg_zee);
    vllbkg_zee->Add(ewkzbkg_zee);
  }
  vllbkg_zee->Add(gambkg_zee);
  vllbkg_zee->Add(qcdbkg_zee);

  // WEN  control region
  TH1* vlbkg_wen   = (TH1*) inputFile->FindObjectAny(("vlbkghistwen_"+observable).c_str());
  TH1* vllbkg_wen  = (TH1*) inputFile->FindObjectAny(("vllbkghistwen_"+observable).c_str());
  TH1* dbbkg_wen   = (TH1*) inputFile->FindObjectAny(("dbkghistwen_"+observable).c_str());
  TH1* ttbkg_wen   = (TH1*) inputFile->FindObjectAny(("tbkghistwen_"+observable).c_str());
  TH1* qbkg_wen    = (TH1*) inputFile->FindObjectAny(("qbkghistwen_"+observable).c_str());
  TH1* gambkg_wen  = (TH1*) inputFile->FindObjectAny(("gbkghistwen_"+observable).c_str());
  TH1* ewkwbkg_wen = (TH1*) inputFile->FindObjectAny(("ewkwbkghistwen_"+observable).c_str());
  TH1* ewkzbkg_wen = (TH1*) inputFile->FindObjectAny(("ewkzbkghistwen_"+observable).c_str());
  vlbkg_wen->Add(vllbkg_wen);
  vlbkg_wen->Add(dbbkg_wen);
  vlbkg_wen->Add(ttbkg_wen);
  vlbkg_wen->Add(qbkg_wen);
  vlbkg_wen->Add(gambkg_wen);
  if(category == Category::VBF){
    vlbkg_wen->Add(ewkwbkg_wen);
    vlbkg_wen->Add(ewkzbkg_wen);
  }
  // WMN  control region
  TH1* vlbkg_wmn   = (TH1*) inputFile->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
  TH1* vllbkg_wmn  = (TH1*) inputFile->FindObjectAny(("vllbkghistwmn_"+observable).c_str());
  TH1* dbbkg_wmn   = (TH1*) inputFile->FindObjectAny(("dbkghistwmn_"+observable).c_str());
  TH1* ttbkg_wmn   = (TH1*) inputFile->FindObjectAny(("tbkghistwmn_"+observable).c_str());
  TH1* qbkg_wmn    = (TH1*) inputFile->FindObjectAny(("qbkghistwmn_"+observable).c_str());
  TH1* gambkg_wmn  = (TH1*) inputFile->FindObjectAny(("gbkghistwmn_"+observable).c_str());
  TH1* ewkwbkg_wmn = (TH1*) inputFile->FindObjectAny(("ewkwbkghistwmn_"+observable).c_str());
  TH1* ewkzbkg_wmn = (TH1*) inputFile->FindObjectAny(("ewkzbkghistwmn_"+observable).c_str());
  vlbkg_wmn->Add(vllbkg_wmn);
  vlbkg_wmn->Add(dbbkg_wmn);
  vlbkg_wmn->Add(ttbkg_wmn);
  vlbkg_wmn->Add(qbkg_wmn);
  vlbkg_wmn->Add(gambkg_wmn);
  if(category == Category::VBF){
    vlbkg_wmn->Add(ewkwbkg_wmn);
    vlbkg_wmn->Add(ewkzbkg_wmn);
  }

  // GAM  control region
  TH1* gbkg_gam   = NULL;
  TH1* qbkg_gam   = NULL;
  TH1* vgbkg_gam   = NULL;
  TH1* vlbkg_gam   = NULL;
  if(category != Category::VBF){
    gbkg_gam = (TH1*) inputFile->FindObjectAny(("gbkghistgam_"+observable).c_str());
    qbkg_gam = (TH1*) inputFile->FindObjectAny(("qbkghistgam_"+observable).c_str());
    vgbkg_gam = (TH1*) inputFile->FindObjectAny(("vgbkghistgam_"+observable).c_str());
    vlbkg_gam = (TH1*) inputFile->FindObjectAny(("vlbkghistgam_"+observable).c_str());
    gbkg_gam->Add(qbkg_gam);
    gbkg_gam->Add(vgbkg_gam);
    gbkg_gam->Add(vlbkg_gam);
  }

  //SYS Unc on ratios
  TH1*  ZG_ewk = NULL;
  TH1*  ZG_re1 = NULL;
  TH1*  ZG_re2 = NULL;
  TH1*  ZG_fa1 = NULL;
  TH1*  ZG_fa2 = NULL;
  TH1*  ZG_pdf = NULL;
  TH1*  ZG_fp  = NULL;
  if(category != Category::VBF){
    ZG_ewk = (TH1*)inputFile->FindObjectAny(("ZG_EWK_"+observable).c_str());
    ZG_re1 = (TH1*)inputFile->FindObjectAny(("ZG_RenScale1_"+observable).c_str());
    ZG_re2 = (TH1*)inputFile->FindObjectAny(("ZG_RenScale2_"+observable).c_str());
    ZG_fa1 = (TH1*)inputFile->FindObjectAny(("ZG_FactScale1_"+observable).c_str());
    ZG_fa2 = (TH1*)inputFile->FindObjectAny(("ZG_FactScale2_"+observable).c_str());
    ZG_pdf = (TH1*)inputFile->FindObjectAny(("ZG_PDF_"+observable).c_str());
    ZG_fp  = (TH1*)inputFile->FindObjectAny(("ZG_Footprint_"+observable).c_str());
  }

  TH1*  ZW_ewk = (TH1*)inputFile->FindObjectAny(("ZW_EWK_"+observable).c_str());
  TH1*  ZW_re1 = (TH1*)inputFile->FindObjectAny(("ZW_RenScale1_"+observable).c_str());
  TH1*  ZW_re2 = (TH1*)inputFile->FindObjectAny(("ZW_RenScale2_"+observable).c_str());
  TH1*  ZW_fa1 = (TH1*)inputFile->FindObjectAny(("ZW_FactScale1_"+observable).c_str());
  TH1*  ZW_fa2 = (TH1*)inputFile->FindObjectAny(("ZW_FactScale2_"+observable).c_str());
  TH1*  ZW_pdf = (TH1*)inputFile->FindObjectAny(("ZW_PDF_"+observable).c_str());

  TH1*  WG_ewk = NULL;
  TH1*  WG_re1 = NULL;
  TH1*  WG_re2 = NULL;
  TH1*  WG_fa1 = NULL;
  TH1*  WG_fa2 = NULL;
  TH1*  WG_pdf = NULL;
  TH1*  WG_fp  = NULL;
  if(category != Category::VBF){
    WG_ewk = (TH1*)inputFile->FindObjectAny(("WG_EWK_"+observable).c_str());
    WG_re1 = (TH1*)inputFile->FindObjectAny(("WG_RenScale1_"+observable).c_str());
    WG_re2 = (TH1*)inputFile->FindObjectAny(("WG_RenScale2_"+observable).c_str());
    WG_fa1 = (TH1*)inputFile->FindObjectAny(("WG_FactScale1_"+observable).c_str());
    WG_fa2 = (TH1*)inputFile->FindObjectAny(("WG_FactScale2_"+observable).c_str());
    WG_pdf = (TH1*)inputFile->FindObjectAny(("WG_PDF_"+observable).c_str());
    WG_fp = (TH1*)inputFile->FindObjectAny(("WG_Footprint_"+observable).c_str());
  }

  TH1* ZWData_ee = (TH1*) data_zee->Clone("ZWData_e");
  ZWData_ee->Divide(data_wen);
  TH1* ZWData_mm = (TH1*) data_zmm->Clone("ZWData_m");
  ZWData_mm->Divide(data_wmn);
  TH1* ZWData_ll = (TH1*) data_zmm->Clone("ZWData_l");
  ZWData_ll->Add(data_zee);
  TH1* temp = (TH1*) data_wmn->Clone("temp"); 
  temp->Add(data_wen);
  ZWData_ll->Divide(temp);
  
  TH1* ZWMC_ee = (TH1*) vllbkg_zee->Clone("ZWMC_e");
  ZWMC_ee->Divide(vlbkg_wen);
  TH1* ZWMC_mm = (TH1*) vllbkg_zmm->Clone("ZWMC_m");
  ZWMC_mm->Divide(vlbkg_wmn);
  TH1* ZWMC_ll = (TH1*) vllbkg_zmm->Clone("ZWMC_l");
  ZWMC_ll->Add(vllbkg_zee);
  temp = (TH1*) vlbkg_wmn->Clone("temp"); 
  temp->Add(vlbkg_wen);
  ZWMC_ll->Divide(temp);

  for(int iBin = 0; iBin < ZWMC_mm->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_mm->GetBinError(iBin+1)*ZWMC_mm->GetBinError(iBin+1);
    err += pow(ZWMC_mm->GetBinContent(iBin+1)*musf,2);
    err += pow(ZWMC_mm->GetBinContent(iBin+1)*mutrack,2);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    if(category == Category::VBF)
      ZWMC_mm->SetBinError(iBin+1,inflateWZ_ewk*sqrt(err));
    else
      ZWMC_mm->SetBinError(iBin+1,sqrt(err));

  }
  
  for(int iBin = 0; iBin < ZWMC_ee->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ee->GetBinError(iBin+1)*ZWMC_ee->GetBinError(iBin+1);
    err += pow(ZWMC_ee->GetBinContent(iBin+1)*elsf,2);
    err += pow(ZWMC_ee->GetBinContent(iBin+1)*eltrack,2);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    if(category == Category::VBF)
      ZWMC_ee->SetBinError(iBin+1,inflateWZ_ewk*sqrt(err));
    else
      ZWMC_ee->SetBinError(iBin+1,sqrt(err));
  }

  for(int iBin = 0; iBin < ZWMC_ll->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ll->GetBinError(iBin+1)*ZWMC_ll->GetBinError(iBin+1);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*elsf,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*eltrack,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*musf,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*mutrack,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*lepveto,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*mettrig,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*eltrig,2);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    if(category == Category::VBF)
      ZWMC_ll->SetBinError(iBin+1,inflateWZ_ewk*sqrt(err));
    else
      ZWMC_ll->SetBinError(iBin+1,sqrt(err));
  }

  TH1* ZGData_mm = NULL;
  TH1* ZGData_ee = NULL;
  TH1* ZGData_ll = NULL;
  TH1* WGData_m  = NULL;
  TH1* WGData_e  = NULL;
  TH1* WGData_l  = NULL;

  TH1* ZGMC_mm = NULL;
  TH1* ZGMC_ee = NULL;
  TH1* ZGMC_ll = NULL;
  TH1* WGMC_m = NULL;
  TH1* WGMC_e = NULL;
  TH1* WGMC_l = NULL;

  if(category != Category::VBF){
    //Ratios Data
    ZGData_mm = (TH1*) data_zmm->Clone("ZGData_mm");
    ZGData_mm->Divide(data_gam);
    ZGData_ee = (TH1*) data_zee->Clone("ZGData_ee");
    ZGData_ee->Divide(data_gam);
    ZGData_ll = (TH1*) data_zmm->Clone("ZGData_ll");
    ZGData_ll->Add(data_zee);
    ZGData_ll->Divide(data_gam);
    
    WGData_m = (TH1*) data_wmn->Clone("WGData_m");
    WGData_m->Divide(data_gam);
    WGData_e = (TH1*) data_wen->Clone("WGData_e");
    WGData_e->Divide(data_gam);
    WGData_l = (TH1*) data_wmn->Clone("WGData_l");
    WGData_l->Add(data_wen);
    WGData_l->Divide(data_gam);
    
    //Ratios MC
    ZGMC_mm = (TH1*) vllbkg_zmm->Clone("ZGMC_mm");
    ZGMC_mm->Divide(gbkg_gam);
    ZGMC_ee = (TH1*) vllbkg_zee->Clone("ZGMC_ee");
    ZGMC_ee->Divide(gbkg_gam);
    ZGMC_ll = (TH1*) vllbkg_zmm->Clone("ZGMC_ll");
    ZGMC_ll->Add(vllbkg_zee);
    ZGMC_ll->Divide(gbkg_gam);
    
    WGMC_m = (TH1*) vlbkg_wmn->Clone("WGMC_m");
    WGMC_m->Divide(gbkg_gam);
    WGMC_e = (TH1*) vlbkg_wen->Clone("WGMC_e");
    WGMC_e->Divide(gbkg_gam);
    WGMC_l = (TH1*) vlbkg_wmn->Clone("WGMC_l");
    WGMC_l->Add(vlbkg_wen);
    WGMC_l->Divide(gbkg_gam);
    
    //Add systematic uncertainties
    for(int iBin = 0; iBin < ZGMC_mm->GetNbinsX(); iBin++){
      double err = 0.;
      err += ZGMC_mm->GetBinError(iBin+1)*ZGMC_mm->GetBinError(iBin+1);
      err += pow(ZGMC_mm->GetBinContent(iBin+1)*musf*2,2);
      err += pow(ZGMC_mm->GetBinContent(iBin+1)*mutrack*2,2);
      err += pow(ZGMC_mm->GetBinContent(iBin+1)*mettrig,2);
      err += pow(ZGMC_mm->GetBinContent(iBin+1)*phtrig,2);
      err += pow(ZGMC_mm->GetBinContent(iBin+1)*phsf,2);
      err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
      ZGMC_mm->SetBinError(iBin+1,sqrt(err));
    }
    
    for(int iBin = 0; iBin < ZGMC_ee->GetNbinsX(); iBin++){
      double err = 0.;
      err += ZGMC_ee->GetBinError(iBin+1)*ZGMC_ee->GetBinError(iBin+1);
      err += pow(ZGMC_ee->GetBinContent(iBin+1)*elsf*2,2);
      err += pow(ZGMC_ee->GetBinContent(iBin+1)*eltrack*2,2);
      err += pow(ZGMC_ee->GetBinContent(iBin+1)*phtrig,2);
      err += pow(ZGMC_ee->GetBinContent(iBin+1)*phsf,2);
      err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
      ZGMC_ee->SetBinError(iBin+1,sqrt(err));
    }
    
    for(int iBin = 0; iBin < ZGMC_ll->GetNbinsX(); iBin++){
      double err = 0.;
      err += ZGMC_ll->GetBinError(iBin+1)*ZGMC_ll->GetBinError(iBin+1);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*musf*2,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*mutrack*2,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*mettrig,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*elsf*2,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*eltrack*2,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*phtrig,2);
      err += pow(ZGMC_ll->GetBinContent(iBin+1)*phsf,2);
      err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
      ZGMC_ll->SetBinError(iBin+1,sqrt(err));
    }
    
    
    if(addWgamma){
      
      for(int iBin = 0; iBin < WGMC_m->GetNbinsX(); iBin++){
	double err = 0.;
	err += WGMC_m->GetBinError(iBin+1)*WGMC_m->GetBinError(iBin+1);
	err += pow(WGMC_m->GetBinContent(iBin+1)*musf,2);
	err += pow(WGMC_m->GetBinContent(iBin+1)*mutrack,2);
	err += pow(WGMC_m->GetBinContent(iBin+1)*phtrig,2);
	err += pow(WGMC_m->GetBinContent(iBin+1)*phsf,2);
	err += pow(WG_ewk->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_re1->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_re2->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_fa1->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_fa2->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_pdf->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	err += pow(WG_fp->GetBinContent(iBin+1)*WGMC_m->GetBinContent(iBin+1), 2);
	WGMC_m->SetBinError(iBin+1,0,sqrt(err));
      }
      
      for(int iBin = 0; iBin < WGMC_e->GetNbinsX(); iBin++){
	double err = 0.;
	err += WGMC_e->GetBinError(iBin+1)*WGMC_e->GetBinError(iBin+1);
	err += pow(WGMC_e->GetBinContent(iBin+1)*elsf,2);
	err += pow(WGMC_e->GetBinContent(iBin+1)*eltrack,2);
	err += pow(WGMC_e->GetBinContent(iBin+1)*phtrig,2);
	err += pow(WGMC_e->GetBinContent(iBin+1)*phsf,2);
	err += pow(WG_ewk->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_re1->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_re2->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_fa1->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_fa2->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_pdf->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	err += pow(WG_fp->GetBinContent(iBin+1)*WGMC_e->GetBinContent(iBin+1), 2);
	WGMC_e->SetBinError(iBin+1,sqrt(err));
      }
      
      for(int iBin = 0; iBin < WGMC_l->GetNbinsX(); iBin++){
	double err = 0.;
	err += WGMC_l->GetBinError(iBin+1)*WGMC_l->GetBinError(iBin+1);
	err += pow(WGMC_l->GetBinContent(iBin+1)*elsf,2);
	err += pow(WGMC_l->GetBinContent(iBin+1)*eltrack,2);
	err += pow(WGMC_l->GetBinContent(iBin+1)*musf,2);
	err += pow(WGMC_l->GetBinContent(iBin+1)*mutrack,2);
	err += pow(WGMC_l->GetBinContent(iBin+1)*phtrig,2);
	err += pow(WGMC_l->GetBinContent(iBin+1)*phsf,2);
	err += pow(WG_ewk->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_re1->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_re2->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_fa1->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_fa2->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_pdf->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	err += pow(WG_fp->GetBinContent(iBin+1)*WGMC_l->GetBinContent(iBin+1), 2);
	WGMC_l->SetBinError(iBin+1,sqrt(err));
      }
    }
  }
  // make plots
  
  makePlot(ZWData_mm,ZWMC_mm,observable,category,observableLatex,"ZW_mm",rebinFactor);
  makePlot(ZWData_ee,ZWMC_ee,observable,category,observableLatex,"ZW_ee",rebinFactor);
  makePlot(ZWData_ll,ZWMC_ll,observable,category,observableLatex,"ZW_ll",rebinFactor);  

  if(category != Category::VBF){
    makePlot(ZGData_mm,ZGMC_mm,observable,category,observableLatex,"ZG_mm",rebinFactor);  
    makePlot(ZGData_ee,ZGMC_ee,observable,category,observableLatex,"ZG_ee",rebinFactor);
    makePlot(ZGData_ll,ZGMC_ll,observable,category,observableLatex,"ZG_ll",rebinFactor);
    
    if(addWgamma){
      makePlot(WGData_m,WGMC_m,observable,category,observableLatex,"WG_m",rebinFactor);
      makePlot(WGData_e,WGMC_e,observable,category,observableLatex,"WG_e",rebinFactor);
      makePlot(WGData_l,WGMC_l,observable,category,observableLatex,"WG_l",rebinFactor);    
    }
  }
}
