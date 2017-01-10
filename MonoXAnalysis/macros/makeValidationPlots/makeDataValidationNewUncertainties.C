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

void makeUncertaintyPlot (TH1* histo1, TH1* histo2, TH1* histo3, TH1* histo4, const string & observableLatex, const string & postfix){

  TH1* histo1_temp = (TH1*) histo1->Clone("histo1_temp");
  TH1* histo2_temp = (TH1*) histo2->Clone("histo2_temp");
  TH1* histo3_temp = (TH1*) histo3->Clone("histo3_temp");
  TH1* histo4_temp = (TH1*) histo4->Clone("histo4_temp");

  for(int iBin = 0; iBin < histo1_temp->GetNbinsX()+1; iBin++)
    histo1_temp->SetBinContent(iBin,1-fabs(1-histo1_temp->GetBinContent(iBin)));
  for(int iBin = 0; iBin < histo2_temp->GetNbinsX()+1; iBin++)
    histo2_temp->SetBinContent(iBin,1-fabs(1-histo2_temp->GetBinContent(iBin)));
  for(int iBin = 0; iBin < histo3_temp->GetNbinsX()+1; iBin++)
    histo3_temp->SetBinContent(iBin,1-fabs(1-histo3_temp->GetBinContent(iBin)));
  for(int iBin = 0; iBin < histo4_temp->GetNbinsX()+1; iBin++)
    histo4_temp->SetBinContent(iBin,1-fabs(1-histo4_temp->GetBinContent(iBin)));

  

  // final plot
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,625);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  histo1->GetXaxis()->SetRangeUser(200,2000);
  histo1->GetYaxis()->SetRangeUser(0.9,1.1);
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2);
  histo1_temp->SetLineColor(kBlack);
  histo1_temp->SetLineWidth(2);
  histo1->GetXaxis()->SetTitle(observableLatex.c_str());
  histo1->GetYaxis()->SetTitle("1 #pm #sigma_{TH}");
  histo1->GetXaxis()->SetTitleOffset(1.1);
  histo1->GetYaxis()->SetTitleOffset(1.1);
  
  histo1->Draw("hist");
  histo1_temp->Draw("hist same");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2);
  histo2_temp->SetLineColor(kRed);
  histo2_temp->SetLineWidth(2);

  histo2->Draw("hist same");
  histo2_temp->Draw("hist same");

  histo3->SetLineColor(kBlue);
  histo3->SetLineWidth(2);
  histo3_temp->SetLineColor(kBlue);
  histo3_temp->SetLineWidth(2);

  histo3->Draw("hist same");
  histo3_temp->Draw("hist same");

  histo4->SetLineColor(kGreen+1);
  histo4->SetLineWidth(2);
  histo4_temp->SetLineColor(kGreen+1);
  histo4_temp->SetLineWidth(2);

  histo4->Draw("hist same");
  histo4_temp->Draw("hist same");

  TLegend leg (0.2,0.7,0.5,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(histo1,"QCD scale","L");
  leg.AddEntry(histo2,"NLO EWK","L");
  leg.AddEntry(histo3,"NNLO Sudakov","L");
  leg.AddEntry(histo4,"QCD-EWK Mixing","L");
  leg.Draw("same");
  

  canvas->SaveAs(("histo_"+postfix+".png").c_str(),"png");
  canvas->SaveAs(("histo_"+postfix+".pdf").c_str(),"pdf");

}

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
    frame2 = pad2->DrawFrame(bins.front(), 0.7, bins.back(), 1.3, "");
  else if(category == Category::monoV)
    frame2 = pad2->DrawFrame(bins.front(), 0.7, bins.back(), 1.3, "");
  else if(category == Category::VBF)
    frame2 = pad2->DrawFrame(bins.front(), 0.7, bins.back(), 1.3, "");
  

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


void makeDataValidationNewUncertainties(string inputFileName, string observable, string observableLatex, int rebinFactor = 1){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  initializeBinning();

  // open the input file with all the templates
  TFile* inputFile = TFile::Open(inputFileName.c_str());
  TFile* kfactorFile_zll = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist/eej.root");
  TFile* kfactorFile_wln = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist/evj.root");

  TH1* data_zmm = (TH1*) inputFile->FindObjectAny(("datahistzmm_"+observable).c_str());
  TH1* data_zee = (TH1*) inputFile->FindObjectAny(("datahistzee_"+observable).c_str());
  TH1* data_wen = (TH1*) inputFile->FindObjectAny(("datahistwen_"+observable).c_str());
  TH1* data_wmn = (TH1*) inputFile->FindObjectAny(("datahistwmn_"+observable).c_str());
    
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

  //// Make Ratios 
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

  // Uncertainties for Z and W --> taken from files
  TH1*  Zll_qcd    = (TH1*) kfactorFile_zll->FindObjectAny("eej_pTV_d1K_NLO");
  TH1*  Zll_nloewk = (TH1*) kfactorFile_zll->FindObjectAny("eej_pTV_kappa_NLO_EW");
  TH1*  Zll_sudewk = (TH1*) kfactorFile_zll->FindObjectAny("eej_pTV_kappa_NNLO_Sud");
  TH1*  Zll_nlokfact = (TH1*) kfactorFile_zll->FindObjectAny("eej_pTV_K_NLO");
  TH1*  Zll_lokfact  = (TH1*) kfactorFile_zll->FindObjectAny("eej_pTV_K_LO");

  TH1*  Wln_qcd    = (TH1*) kfactorFile_wln->FindObjectAny("evj_pTV_d1K_NLO");
  TH1*  Wln_nloewk = (TH1*) kfactorFile_wln->FindObjectAny("evj_pTV_kappa_NLO_EW");
  TH1*  Wln_sudewk = (TH1*) kfactorFile_wln->FindObjectAny("evj_pTV_kappa_NNLO_Sud");
  TH1*  Wln_nlokfact = (TH1*) kfactorFile_wln->FindObjectAny("evj_pTV_K_NLO");
  TH1*  Wln_lokfact  = (TH1*) kfactorFile_wln->FindObjectAny("evj_pTV_K_LO");

  // Uncertainty in the multiplicative approach
  TH1* ZW_qcd_up = (TH1*) Zll_qcd->Clone("ZW_qcd_up");
  ZW_qcd_up->Reset("ICES");
  TH1* ZW_qcd_dw = (TH1*) Zll_qcd->Clone("ZW_qcd_dw");
  ZW_qcd_dw->Reset("ICES");
  for(int iBin = 0; iBin < ZW_qcd_up->GetNbinsX()+1; iBin++){
    ZW_qcd_dw->SetBinContent(iBin,(Zll_nlokfact->GetBinContent(iBin)-Zll_qcd->GetBinContent(iBin))*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)));
    ZW_qcd_up->SetBinContent(iBin,(Zll_nlokfact->GetBinContent(iBin)+Zll_qcd->GetBinContent(iBin))*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)));
  }

  for(int iBin = 0; iBin < ZW_qcd_up->GetNbinsX()+1; iBin++){
    ZW_qcd_dw->SetBinContent(iBin,ZW_qcd_dw->GetBinContent(iBin)/((Wln_nlokfact->GetBinContent(iBin)-Wln_qcd->GetBinContent(iBin))*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin))));
    ZW_qcd_up->SetBinContent(iBin,ZW_qcd_up->GetBinContent(iBin)/((Wln_nlokfact->GetBinContent(iBin)+Wln_qcd->GetBinContent(iBin))*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin))));
  }

  TH1* ZW_qcd = (TH1*) ZW_qcd_up->Clone("ZW_qcd");
  ZW_qcd->Reset("ICES");
  for(int iBin=0; iBin < ZW_qcd->GetNbinsX()+1; iBin++)
    ZW_qcd->SetBinContent(iBin,1+fabs((ZW_qcd_up->GetBinContent(iBin)-ZW_qcd_dw->GetBinContent(iBin))));
			  
  // NLO EWK uncertainty --> correlation 100%
  TH1* ZW_nloewk_up = (TH1*) Zll_nloewk->Clone("ZW_nloewk");
  ZW_nloewk_up->Reset("ICES");
  TH1* ZW_nloewk_dw = (TH1*) Zll_nloewk->Clone("ZW_nloewk");
  ZW_nloewk_dw->Reset("ICES");

  for(int iBin = 0; iBin < ZW_nloewk_up->GetNbinsX()+1; iBin++){
    ZW_nloewk_up->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)+0.1*Zll_nloewk->GetBinContent(iBin)));
    ZW_nloewk_dw->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)-0.1*Zll_nloewk->GetBinContent(iBin)));
  }
  for(int iBin = 0; iBin < ZW_nloewk_up->GetNbinsX()+1; iBin++){
    ZW_nloewk_up->SetBinContent(iBin,ZW_nloewk_up->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)+0.1*Wln_nloewk->GetBinContent(iBin))));
    ZW_nloewk_dw->SetBinContent(iBin,ZW_nloewk_dw->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)-0.1*Wln_nloewk->GetBinContent(iBin))));
  }
  
  TH1* ZW_nloewk = (TH1*) ZW_nloewk_up->Clone("ZW_nloewk");
  ZW_nloewk->Reset("ICES");
  for(int iBin=0; iBin < ZW_nloewk->GetNbinsX()+1; iBin++)
    ZW_nloewk->SetBinContent(iBin,1+fabs((ZW_nloewk_up->GetBinContent(iBin)-ZW_nloewk_dw->GetBinContent(iBin)))/2);
  
  // NNLO Sudakov uncertainty --> correlation 100%
  TH1* ZW_sudewk_up = (TH1*) Zll_nloewk->Clone("ZW_sudewk");
  ZW_sudewk_up->Reset("ICES");
  TH1* ZW_sudewk_dw = (TH1*) Zll_nloewk->Clone("ZW_sudewk");
  ZW_sudewk_dw->Reset("ICES");

  for(int iBin = 0; iBin < ZW_sudewk_up->GetNbinsX()+1; iBin++){
    ZW_sudewk_up->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)+2/3*Zll_sudewk->GetBinContent(iBin)*Zll_nloewk->GetBinContent(iBin)));
    ZW_sudewk_dw->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin)-2/3*Zll_sudewk->GetBinContent(iBin)*Zll_nloewk->GetBinContent(iBin)));
  }
  for(int iBin = 0; iBin < ZW_sudewk_up->GetNbinsX()+1; iBin++){
    ZW_sudewk_up->SetBinContent(iBin,ZW_sudewk_up->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)+2.3*Wln_sudewk->GetBinContent(iBin)*Wln_nloewk->GetBinContent(iBin))));
    ZW_sudewk_dw->SetBinContent(iBin,ZW_sudewk_dw->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)-2/3*Wln_sudewk->GetBinContent(iBin)*Wln_nloewk->GetBinContent(iBin))));
  }
  
  TH1* ZW_sudewk = (TH1*) ZW_sudewk_up->Clone("ZW_sudewk");
  ZW_sudewk->Reset("ICES");
  for(int iBin=0; iBin < ZW_sudewk->GetNbinsX()+1; iBin++)
    ZW_sudewk->SetBinContent(iBin,1+fabs((ZW_sudewk_up->GetBinContent(iBin)-ZW_sudewk_dw->GetBinContent(iBin)))/2);

  // Mixed QCD EKW
  TH1* ZW_ewkqcd_up = (TH1*) Zll_qcd->Clone("Zll_ewkqcd_up");
  ZW_ewkqcd_up->Reset("ICES");
  TH1* ZW_ewkqcd_dw = (TH1*) Zll_qcd->Clone("Zll_ewkqcd_dw");
  ZW_ewkqcd_dw->Reset("ICES");

  for(int iBin=0; iBin < ZW_ewkqcd_up->GetNbinsX()+1; iBin++){
    ZW_ewkqcd_up->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_lokfact->GetBinContent(iBin)/(Zll_lokfact->GetBinContent(iBin)+0.5*(Zll_nlokfact->GetBinContent(iBin)-Zll_lokfact->GetBinContent(iBin)))*(Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin))));
    ZW_ewkqcd_dw->SetBinContent(iBin,Zll_nlokfact->GetBinContent(iBin)*(1+Zll_lokfact->GetBinContent(iBin)/(Zll_lokfact->GetBinContent(iBin)-0.5*(Zll_nlokfact->GetBinContent(iBin)-Zll_lokfact->GetBinContent(iBin)))*(Zll_nloewk->GetBinContent(iBin)+Zll_sudewk->GetBinContent(iBin))));
  }

  for(int iBin=0; iBin < ZW_ewkqcd_up->GetNbinsX()+1; iBin++){
    ZW_ewkqcd_up->SetBinContent(iBin,ZW_ewkqcd_up->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_lokfact->GetBinContent(iBin)/(Wln_lokfact->GetBinContent(iBin)+0.5*(Wln_nlokfact->GetBinContent(iBin)-Wln_lokfact->GetBinContent(iBin)))*(Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)))));
    ZW_ewkqcd_dw->SetBinContent(iBin,ZW_ewkqcd_dw->GetBinContent(iBin)/(Wln_nlokfact->GetBinContent(iBin)*(1+Wln_lokfact->GetBinContent(iBin)/(Wln_lokfact->GetBinContent(iBin)-0.5*(Wln_nlokfact->GetBinContent(iBin)-Wln_lokfact->GetBinContent(iBin)))*(Wln_nloewk->GetBinContent(iBin)+Wln_sudewk->GetBinContent(iBin)))));
  }
  
  TH1* ZW_ewkqcd = (TH1*) ZW_sudewk_up->Clone("ZW_ewkqcd");
  ZW_ewkqcd->Reset("ICES");
  for(int iBin=0; iBin < ZW_ewkqcd->GetNbinsX()+1; iBin++)
    ZW_ewkqcd->SetBinContent(iBin,1+fabs((ZW_ewkqcd_up->GetBinContent(iBin)-ZW_ewkqcd_dw->GetBinContent(iBin)))/2);
  
  // Plot theory uncertainty
  makeUncertaintyPlot(ZW_qcd,ZW_nloewk,ZW_sudewk,ZW_ewkqcd,"boson p_{T} [GeV]","ZW_unc");
  
  for(int iBin = 0; iBin < ZWMC_mm->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_mm->GetBinError(iBin+1)*ZWMC_mm->GetBinError(iBin+1);
    err += pow(ZWMC_mm->GetBinContent(iBin+1)*musf,2);
    err += pow(ZWMC_mm->GetBinContent(iBin+1)*mutrack,2);
    err += pow(fabs(1-ZW_qcd->GetBinContent(ZW_qcd->FindBin(ZWMC_mm->GetBinCenter(iBin+1))))*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_nloewk->GetBinContent(ZW_nloewk->FindBin(ZWMC_mm->GetBinCenter(iBin+1))))*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_sudewk->GetBinContent(ZW_sudewk->FindBin(ZWMC_mm->GetBinCenter(iBin+1))))*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_ewkqcd->GetBinContent(ZW_ewkqcd->FindBin(ZWMC_mm->GetBinCenter(iBin+1))))*ZWMC_mm->GetBinContent(iBin+1), 2);
    ZWMC_mm->SetBinError(iBin+1,sqrt(err));
  }
  
  for(int iBin = 0; iBin < ZWMC_ee->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ee->GetBinError(iBin+1)*ZWMC_ee->GetBinError(iBin+1);
    err += pow(ZWMC_ee->GetBinContent(iBin+1)*elsf,2);
    err += pow(ZWMC_ee->GetBinContent(iBin+1)*eltrack,2);
    err += pow(fabs(1-ZW_qcd->GetBinContent(ZW_qcd->FindBin(ZWMC_ee->GetBinCenter(iBin+1))))*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_nloewk->GetBinContent(ZW_nloewk->FindBin(ZWMC_ee->GetBinCenter(iBin+1))))*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_sudewk->GetBinContent(ZW_sudewk->FindBin(ZWMC_ee->GetBinCenter(iBin+1))))*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_ewkqcd->GetBinContent(ZW_ewkqcd->FindBin(ZWMC_ee->GetBinCenter(iBin+1))))*ZWMC_ee->GetBinContent(iBin+1), 2);
    ZWMC_ee->SetBinError(iBin+1,sqrt(err));
  }

  for(int iBin = 0; iBin < ZWMC_ll->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ll->GetBinError(iBin+1)*ZWMC_ll->GetBinError(iBin+1);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*elsf,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*eltrack,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*musf,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*mutrack,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*mettrig,2);
    err += pow(ZWMC_ll->GetBinContent(iBin+1)*eltrig,2);
    err += pow(fabs(1-ZW_qcd->GetBinContent(ZW_qcd->FindBin(ZWMC_ll->GetBinCenter(iBin+1))))*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_nloewk->GetBinContent(ZW_nloewk->FindBin(ZWMC_ll->GetBinCenter(iBin+1))))*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_sudewk->GetBinContent(ZW_sudewk->FindBin(ZWMC_ll->GetBinCenter(iBin+1))))*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(fabs(1-ZW_ewkqcd->GetBinContent(ZW_ewkqcd->FindBin(ZWMC_ll->GetBinCenter(iBin+1))))*ZWMC_ll->GetBinContent(iBin+1), 2);
    ZWMC_ll->SetBinError(iBin+1,sqrt(err));
  }

  makePlot(ZWData_mm,ZWMC_mm,observable,Category::monojet,observableLatex,"ZW_mm",rebinFactor);
  makePlot(ZWData_ee,ZWMC_ee,observable,Category::monojet,observableLatex,"ZW_ee",rebinFactor);
  makePlot(ZWData_ll,ZWMC_ll,observable,Category::monojet,observableLatex,"ZW_ll",rebinFactor);  
}
