#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void prepostSig_COMB_withSignals(string fitFilename, 
				 string templateFile,
				 string observable, 
				 Category category) {

  string dir = "ch1_ch1";
  if(category == Category::monoV)
    dir = "ch2_ch1";
  
      
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->SetLineColor(0);
  pad2->SetGridy();

  TColor *color; // for color definition with alpha                                                                                                                             
  TFile* pfile   = new TFile(fitFilename.c_str());
  TFile* monoj_v, *monow_v, *monoz_v, *monoj_av, *monow_av, *monoz_av, *higgs;
  if(category == Category::monoV){
    monoj_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonov_13TeV_v1.root","READ");
    monow_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonov_13TeV_v1.root","READ");
    monoz_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonov_13TeV_v1.root","READ");
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonov_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonov_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonov_13TeV_v1.root","READ");
    higgs = new TFile("~/work/MONOJET_ANALYSIS/CMSSW_7_4_16/src/AnalysisCode/MonoXAnalysis/macros/monoV_hinv_forCombination/templates_met_v2.root","READ");
  }
  else{
    monoj_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonojet_13TeV_v1.root","READ");
    monow_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonojet_13TeV_v1.root","READ");
    monoz_v = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonojet_13TeV_v1.root","READ");
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");
    higgs = new TFile("~/work/MONOJET_ANALYSIS/CMSSW_7_4_16/src/AnalysisCode/MonoXAnalysis/macros/monoj_hinv_forCombination/templates_met_v2.root","READ");
  }

  // in case of b-only fit just dispaly three possible signal on the stack
  TH1* mjhist_v = NULL;
  TH1* mwhist_v = NULL;
  TH1* mzhist_v = NULL;
  TH1* mjhist_av = NULL;
  TH1* mwhist_av = NULL;
  TH1* mzhist_av = NULL;

  TH1* ggHhist = NULL;
  TH1* vbfhist = NULL;
  TH1* wHhist = NULL;
  TH1* zHhist = NULL;

  mjhist_v = (TH1*) monoj_v->FindObjectAny("signal_signal_80018000001");
  mwhist_v = (TH1*) monow_v->FindObjectAny("signal_signal_80018000001");
  mzhist_v = (TH1*) monoz_v->FindObjectAny("signal_signal_80018000001");
  mjhist_v->Scale(1.,"width");
  mwhist_v->Scale(1.,"width");
  mzhist_v->Scale(1.,"width");

  mjhist_v->Add(mwhist_v);
  mjhist_v->Add(mzhist_v);

  mjhist_av = (TH1*) monoj_av->FindObjectAny("signal_signal_80116000001");
  mwhist_av = (TH1*) monow_av->FindObjectAny("signal_signal_80116000001");
  mzhist_av = (TH1*) monoz_av->FindObjectAny("signal_signal_80116000001");
  mjhist_av->Scale(1.,"width");
  mwhist_av->Scale(1.,"width");
  mzhist_av->Scale(1.,"width");
  mjhist_av->Add(mwhist_av);
  mjhist_av->Add(mzhist_av);
  
  ggHhist = (TH1*) higgs->FindObjectAny("ggHhist_125_met");
  vbfhist = (TH1*) higgs->FindObjectAny("vbfHhist_125_met");
  wHhist = (TH1*) higgs->FindObjectAny("wHhist_125_met");
  zHhist = (TH1*) higgs->FindObjectAny("zHhist_125_met");
  ggzHhist = (TH1*) higgs->FindObjectAny("ggzHhist_125_met");
  ggHhist->Scale(1.,"width");
  vbfhist->Scale(1.,"width");
  wHhist->Scale(1.,"width");
  zHhist->Scale(1.,"width");
  ggzHhist->Scale(1.,"width");
  if(category == Category::monoV){
    mjhist_v->Scale(1.68);
    mjhist_av->Scale(1.68);
    ggHhist->Scale(5.55);
    vbfhist->Scale(5.05);
    wHhist->Scale(5.05);
    zHhist->Scale(5.05);
    ggzHhist->Scale(5.05);
    ggHhist->Add(vbfhist);
    ggHhist->Add(wHhist);
    ggHhist->Add(zHhist);
    ggHhist->Add(ggzHhist);
  }
  else{
    mjhist_v->Scale(1.68);
    mjhist_av->Scale(1.68);
    ggHhist->Scale(5.97);
    vbfhist->Scale(5.38);
    wHhist->Scale(5.38);
    zHhist->Scale(5.38);
    ggzHhist->Scale(5.38);
    ggHhist->Add(vbfhist);
    ggHhist->Add(wHhist);
    ggHhist->Add(zHhist);
    ggHhist->Add(ggzHhist);
  }
  TH1* znhist = NULL;
  TH1* zlhist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* qchist = NULL;
  TH1* gmhist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* tohist = NULL;
  TH1* tphist = NULL;
  TH1* sighist = NULL;

  znhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Znunu").c_str());    
  zlhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/ZJets").c_str());    
  wlhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/WJets").c_str());    
  tthist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Top").c_str());    
  dihist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Dibosons").c_str());    
  ewkwhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/EWKW").c_str());    
  ewkzhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/EWKZ").c_str());    
  qchist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/QCD").c_str());    
  gmhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/GJets").c_str());    
  tohist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/total_background").c_str());    
  tphist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());    

  TH1* dthist = NULL;
  TFile* dfile = new TFile(templateFile.c_str(),"READ");
  dthist = (TH1*)dfile->FindObjectAny(("datahist_"+observable).c_str());
  dthist->Scale(1.0,"width");
  
  //signal style  
  if(mjhist_v){
    mjhist_v->SetFillColor(0);
    mjhist_v->SetFillStyle(0);
    mjhist_v->SetLineColor(kBlue);
    mjhist_v->SetLineWidth(3);
    mjhist_v->SetMarkerSize(0);
  }
  
  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kBlack);
    ggHhist->SetLineWidth(3);
    ggHhist->SetMarkerSize(0);
  }
    
  if(mjhist_av){
    mjhist_av->SetFillColor(0);
    mjhist_av->SetFillStyle(0);
    mjhist_av->SetLineColor(kBlue);
    mjhist_av->SetLineWidth(3);
    mjhist_av->SetMarkerSize(0);
  }

  qchist->SetFillColor(TColor::GetColor("#F1F1F2"));
  qchist->SetLineColor(kBlack);

  gmhist->SetFillColor(TColor::GetColor("#9A9EAB"));
  gmhist->SetLineColor(TColor::GetColor("#9A9EAB"));
  zlhist->SetFillColor(TColor::GetColor("#9A9EAB"));  
  zlhist->SetLineColor(kBlack);
  zlhist->Add(gmhist);

  znhist->SetFillColor(TColor::GetColor("#3A8C4C"));
  znhist->SetLineColor(kBlack);

  wlhist->SetFillColor(TColor::GetColor("#FAAF08"));
  wlhist->SetLineColor(kBlack);

  dihist->SetFillColor(TColor::GetColor("#4897D8"));
  dihist->SetLineColor(kBlack);

  tthist->SetFillColor(TColor::GetColor("#CF3721"));
  tthist->SetLineColor(kBlack);
  
  // make the stack
  THStack* stack = new THStack("stack", "stack");
  stack->Add(qchist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(dihist);
  stack->Add(wlhist);
  stack->Add(znhist);


  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*250);
  else
    frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*500);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);
  frame ->Draw();

  CMS_lumi(canvas,"36.2");

  TLatex* categoryLabel = new TLatex();
  categoryLabel->SetNDC();
  categoryLabel->SetTextSize(0.5*canvas->GetTopMargin());
  categoryLabel->SetTextFont(42);
  categoryLabel->SetTextAlign(11);
  if(category == Category::monojet)
    categoryLabel ->DrawLatex(0.175,0.80,"monojet");
  else if(category == Category::monoV)
    categoryLabel ->DrawLatex(0.175,0.80,"mono-V");
  else if(category == Category::VBF)
    categoryLabel ->DrawLatex(0.175,0.80,"VBF");
  categoryLabel->Draw("same");

  stack ->Draw("HIST SAME");
  //mjhist_v->Draw("hist same");
  mjhist_av->Draw("hist same");
  ggHhist->Draw("hist same");

  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetLineWidth(1);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("PE SAME");
  
  TLegend* leg = new TLegend(0.50, 0.55, 0.92, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);  

  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  // ratio plot
  pad2->Draw();
  pad2->cd();

  TH1* frame2 = (TH1*) dthist->Clone("frame2");
  frame2->Reset();
  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(510);
  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");

  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->Draw("AXIS");
  frame2->Draw("AXIG same");

  // for post-fit pre-fit data/mc
  TH1* dphist = (TH1*)dthist->Clone("dahist");
  TH1* dahist = (TH1*)dthist->Clone("dahist");

  dphist->SetLineColor(kRed);
  dphist->SetMarkerColor(kRed);

  dahist->SetLineColor(kBlue);
  dahist->SetMarkerColor(kBlue);

  dphist->SetMarkerSize(1);
  dphist->SetMarkerStyle(24);
  dahist->SetMarkerSize(1);
  dahist->SetMarkerStyle(20);

  TH1* mphist = (TH1*) tphist->Clone("mphist");
  TH1* mchist = (TH1*) tphist->Clone("mchist");
  TH1* unhist = (TH1*) tphist->Clone("unhist");
  mchist->Reset();
  unhist->Reset();
  mchist->Add(qchist);
  mchist->Add(wlhist);
  mchist->Add(zlhist);
  mchist->Add(tthist);
  mchist->Add(dihist);
  mchist->Add(znhist);


  for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0);
  for (int i = 1; i <= mphist->GetNbinsX(); i++) mphist->SetBinError(i, 0);

  // ratio data/post-fit
  dahist->Divide(mchist);
  // ratio data/pre-fit
  dphist->Divide(mphist);
  
  tohist->Divide(mchist);
  tohist->SetLineColor(0);
  tohist->SetMarkerColor(0);
  tohist->SetMarkerSize(0);
  tohist->SetFillColor(kGray);
  
  dahist->SetStats(kFALSE);
  dphist->SetStats(kFALSE);
  
  // line at 1
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetLineWidth(2);
  unhist->SetFillColor(0);
  
  dahist->GetXaxis()->SetLabelOffset(999999);
  dahist->GetXaxis()->SetLabelSize(0);
  dahist->GetXaxis()->SetTitleOffset(999999);
  dahist->GetXaxis()->SetTitleSize(0);
    
  tohist->Draw("E2 SAME");
  unhist->Draw("SAME");
  dphist->Draw("P0E1 SAME");
  dahist->Draw("P0E1 SAME");

  TLegend* leg2 = new TLegend(0.14,0.24,0.40,0.28,NULL,"brNDC");
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(dphist,"pre-fit","PLE");
  leg2->AddEntry(dahist,"post-fit","PLE");
  leg2->Draw("same");


  canvas->cd();
  leg->AddEntry(dthist, "Data", "PEL");

  leg->AddEntry(znhist, "Z #rightarrow #nu#nu", "F");
  leg->AddEntry(wlhist, "W #rightarrow l#nu", "F");
  leg->AddEntry(dihist, "WW/WZ/ZZ", "F");
  leg->AddEntry(tthist, "Top Quark", "F");
  leg->AddEntry(zlhist, "Z/#gamma #rightarrow ll, #gamma+jets", "F");
  leg->AddEntry(qchist, "QCD", "F");
  leg->AddEntry(mjhist_av, "Axial-Vector, M_{med} = 1.6 TeV","L");
  leg->AddEntry(ggHhist,   "Higgs invisible, m_{H} = 125 GeV","L");
  
  leg->Draw("SAME");  
  pad2->RedrawAxis("sameaxis");
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();

  if(category == Category::monojet){
    canvas->SaveAs("postfit_sig_MJ.pdf");
    canvas->SaveAs("postfit_sig_MJ.png");
    canvas->SaveAs("postfit_sig_MJ.root");
  }
  else{
    canvas->SaveAs("postfit_sig_MV.pdf");
    canvas->SaveAs("postfit_sig_MV.png");
  }
}

