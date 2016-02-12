#include "CMS_lumi.h"
#include "makehist.h"

float minTau2Tau1 = 0.1;

void makeControlPlots(string templateFileName, 
		      int category, 
		      string observable, string observableLatex, 
		      string controlRegion, 
		      bool blind, bool isLog,
		      bool plotResonant = false,
		      string mediatorMass = "1000",
		      string DMMass = "50",
		      int signalScale = 100) {

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.28);
  pad2->SetTickx();
  pad2->SetTicky();

  TFile* inputFile = new TFile(templateFileName.c_str());

  TH1* datahist = NULL;
  TH1* qcdhist  = NULL;
  TH1* vllhist  = NULL;
  TH1* vnnhist  = NULL;
  TH1* vlhist   = NULL;
  TH1* dbhist   = NULL;
  TH1* tophist  = NULL;
  TH1* tophist_matched   = NULL;
  TH1* tophist_unmatched = NULL;
  TH1* gamhist  = NULL;

  TH1* monoJhist  = NULL;
  TH1* monoWhist  = NULL;
  TH1* monoZhist  = NULL;

  // take the templates
  if(controlRegion == "gam"){
    datahist = (TH1*)inputFile->Get(("datahistgam_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghistgam_"+observable).c_str());
    gamhist  = (TH1*)inputFile->Get(("gbkghistgam_"+observable).c_str());
  }
  else if(controlRegion == "zmm"){
    datahist = (TH1*)inputFile->Get(("datahistzmm_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghistzmm_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghistzmm_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghistzmm_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghistzmm_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghistzmm_"+observable).c_str());  
  }
  else if(controlRegion == "zee"){
    datahist = (TH1*)inputFile->Get(("datahistzee_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghistzee_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghistzee_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghistzee_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghistzee_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghistzee_"+observable).c_str());  
  }
  else if(controlRegion == "wmn"){
    datahist = (TH1*)inputFile->Get(("datahistwmn_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghistwmn_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghistwmn_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghistwmn_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghistwmn_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghistwmn_"+observable).c_str());  
  }
  else if(controlRegion == "wen"){
    datahist = (TH1*)inputFile->Get(("datahistwen_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghistwen_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghistwen_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghistwen_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghistwen_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghistwen_"+observable).c_str());  
  }
  else if(controlRegion == "topmu" and plotResonant){
    datahist = (TH1*)inputFile->Get(("datahisttopmu_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghisttopmu_"+observable).c_str());
    tophist_matched   = (TH1*)inputFile->Get(("tbkghist_matchedtopmu_"+observable).c_str());
    tophist_unmatched = (TH1*)inputFile->Get(("tbkghist_unmatchedtopmu_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghisttopmu_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghisttopmu_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghisttopmu_"+observable).c_str());  
  }
  else if(controlRegion == "topmu" and not plotResonant){
    datahist = (TH1*)inputFile->Get(("datahisttopmu_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghisttopmu_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghisttopmu_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghisttopmu_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghisttopmu_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghisttopmu_"+observable).c_str());  
  }

  else if(controlRegion == "topel" and not plotResonant){
    datahist = (TH1*)inputFile->Get(("datahisttopel_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghisttopel_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghisttopel_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghisttopel_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghisttopel_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghisttopel_"+observable).c_str());  
  }

  else if(controlRegion == "topel" and plotResonant){
    datahist = (TH1*)inputFile->Get(("datahisttopel_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghisttopel_"+observable).c_str());
    tophist_matched   = (TH1*)inputFile->Get(("tbkghist_matchedtopel_"+observable).c_str());
    tophist_unmatched = (TH1*)inputFile->Get(("tbkghist_unmatchedtopel_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("vlbkghisttopel_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("vllbkghisttopel_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghisttopel_"+observable).c_str());  
  }

  else if(controlRegion == "SR"){
    datahist = (TH1*)inputFile->Get(("datahist_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->Get(("qbkghist_"+observable).c_str());
    tophist  = (TH1*)inputFile->Get(("tbkghist_"+observable).c_str());
    vlhist   = (TH1*)inputFile->Get(("wjethist_"+observable).c_str());
    vllhist  = (TH1*)inputFile->Get(("zjethist_"+observable).c_str());
    vnnhist  = (TH1*)inputFile->Get(("zinvhist_"+observable).c_str());
    dbhist   = (TH1*)inputFile->Get(("dbkghist_"+observable).c_str());  
    
    monoJhist = (TH1*)inputFile->Get(("monoJhist_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    monoWhist = (TH1*)inputFile->Get(("monoWhist_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    monoZhist = (TH1*)inputFile->Get(("monoZhist_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());    
  }
  
  // BLIND OPTION
  if(blind and controlRegion == "SR"){

    for (int i = 0; i <= datahist->GetNbinsX(); i++) {
      double yield = 0.0;
      yield += vllhist->GetBinContent(i);
      yield += vlhist->GetBinContent(i);
      yield += tophist->GetBinContent(i);
      yield += dbhist->GetBinContent(i);
      yield += qcdhist->GetBinContent(i);
      yield += vnnhist->GetBinContent(i);
      datahist->SetBinContent(i, yield);
      datahist->SetBinError(i, 0.);
    }
  }
  
  //SCALE BIN WIDTH
  if(datahist)
    datahist->Scale(1.0,"width");
  if(qcdhist)
    qcdhist->Scale(1.0,"width");
  if(tophist)
    tophist->Scale(1.0,"width");
  if(tophist_matched)
    tophist_matched->Scale(1.0,"width");
  if(tophist_unmatched)
    tophist_unmatched->Scale(1.0,"width");
  if(vlhist)
    vlhist->Scale(1.0,"width");
  if(vllhist)
    vllhist->Scale(1.0,"width");
  if(vnnhist)
    vnnhist->Scale(1.0,"width");
  if(dbhist)
    dbhist->Scale(1.0,"width");
  if(gamhist)
    gamhist->Scale(1.0,"width");

  if(controlRegion == "SR")
    qcdhist->Scale(2.);

  if(monoJhist){
    monoJhist->Scale(1.0,"width");
  }
  if(monoWhist){
    monoWhist->Scale(1.0,"width");
    monoWhist->Scale(signalScale);
  }
  if(monoZhist){
    monoZhist->Scale(1.0,"width");
    monoZhist->Scale(signalScale);
  }

  // set colors
  if(datahist){
    datahist->SetLineColor(kBlack);
    datahist->SetLineWidth(2);
    datahist->SetMarkerColor(kBlack);
    datahist->SetMarkerStyle(20);
    datahist->SetMarkerSize(1.);
  }

  if(vnnhist)  {
    vnnhist->SetFillColor(kGreen+1);
    vnnhist->SetLineColor(kBlack);
  }
  if(vllhist){
    vllhist->SetFillColor(kCyan);
    vllhist->SetLineColor(kBlack);
  }
  if(vlhist){
    vlhist->SetFillColor(kRed);
    vlhist->SetLineColor(kBlack);
  }
  if(tophist){
    tophist->SetFillColor(kBlue);
    tophist->SetLineColor(kBlack);
  }
  if(tophist_matched){
    tophist_matched->SetFillColor(kBlue);
    tophist_matched->SetLineColor(kBlack);
  }
  if(tophist_unmatched){
    tophist_unmatched->SetFillColor(kOrange);
    tophist_unmatched->SetLineColor(kBlack);
  }
  if(dbhist){
    dbhist->SetFillColor(kViolet);
    dbhist->SetLineColor(kBlack);
  }
  if(qcdhist) {
    qcdhist->SetFillColor(kGray+1);
    qcdhist->SetLineColor(kBlack);
  }
  if(gamhist){
    gamhist->SetFillColor(kOrange);
    gamhist->SetLineColor(kBlack);
  }

  if(monoJhist){
    monoJhist->SetFillColor(0);
    monoJhist->SetFillStyle(0);
    monoJhist->SetLineColor(kBlack);
    monoJhist->SetLineWidth(2);
  }

  if(monoWhist){
    monoWhist->SetFillColor(0);
    monoWhist->SetFillStyle(0);
    monoWhist->SetLineColor(kBlack);
    monoWhist->SetLineWidth(2);
    monoWhist->SetLineStyle(7);
  }

  if(monoZhist){
    monoZhist->SetFillColor(0);
    monoZhist->SetFillStyle(0);
    monoZhist->SetLineColor(kBlack);
    monoZhist->SetLineWidth(2);
    monoZhist->SetLineStyle(4);
  }

  THStack* stack = new THStack("stack", "stack");
  if(controlRegion == "gam"){
    stack->Add(qcdhist);
    stack->Add(gamhist);
  }
  else if(controlRegion == "zmm" or controlRegion == "zee"){
    stack->Add(qcdhist);
    stack->Add(vlhist);
    stack->Add(tophist);
    stack->Add(dbhist);
    stack->Add(vllhist);
  }
  else if(controlRegion == "wmn" or controlRegion == "wen"){
    stack->Add(qcdhist);
    stack->Add(vllhist);
    stack->Add(tophist);
    stack->Add(dbhist);
    stack->Add(vlhist);
  }
  else if((controlRegion == "topmu" or controlRegion == "topel") and not plotResonant){
    stack->Add(qcdhist);
    stack->Add(vllhist);
    stack->Add(dbhist);
    stack->Add(vlhist);
    stack->Add(tophist);    
  }
  else if((controlRegion == "topmu" or controlRegion == "topel") and plotResonant){
    stack->Add(qcdhist);
    stack->Add(vllhist);
    stack->Add(dbhist);
    stack->Add(vlhist);
    stack->Add(tophist_unmatched);    
    stack->Add(tophist_matched);    
  }
  else if(controlRegion == "SR"){
    stack->Add(qcdhist);
    stack->Add(dbhist);
    stack->Add(tophist);
    stack->Add(vllhist);
    stack->Add(vlhist);
    stack->Add(vnnhist);
  }

  TH1* frame = NULL;
  vector<float> bins = selectBinning(observable,category);

  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);  
  pad1->Draw();
  pad1->cd();

  float xMin = bins.front();
  if(observable == "tau2tau1")
    xMin = minTau2Tau1;
  float xMax = bins.back();

  if(category <= 1 and isLog)
    frame =  pad1->DrawFrame(xMin, 1.5e-4, xMax, datahist->GetMaximum()*1000000, "");
  else if(category <= 1 and not isLog)
    frame =  pad1->DrawFrame(xMin, 1.5e-4, xMax, datahist->GetMaximum()*1.5, "");
  else if(category > 1 and isLog)
    frame =  pad1->DrawFrame(xMin, 1.5e-4, xMax, datahist->GetMaximum()*1000000, "");
  else
    frame =  pad1->DrawFrame(xMin, 1.5e-4, xMax, datahist->GetMaximum()*2.5, "");
    
  frame->GetXaxis()->SetTitle(observableLatex.c_str());
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);

  frame ->Draw();
  CMS_lumi(pad1, 4, 0, true);
  stack ->Draw("HIST SAME");
  datahist->Draw("P SAME");

  if(controlRegion == "SR"){
    monoJhist->Draw("hist same");
    monoWhist->Draw("hist same");
    monoZhist->Draw("hist same");
  }

 TLegend* leg = NULL;
  if(controlRegion == "gam")
    leg = new TLegend(0.58, 0.66, 0.85, 0.92);
  else if(controlRegion == "SR" and isLog){
    leg = new TLegend(0.32, 0.42, 0.85, 0.92);
    if(observable != "njet")
      frame->GetYaxis()->SetRangeUser(1.5e-4,datahist->GetMaximum()*1000000);
    else
      frame->GetYaxis()->SetRangeUser(1.5e-4,datahist->GetMaximum()*1000000000);
  }
  else if(controlRegion == "SR" and not isLog){
    leg = new TLegend(0.32, 0.42, 0.85, 0.92);
    frame->GetYaxis()->SetRangeUser(1.5e-4,datahist->GetMaximum()*2.5);
  }
  else
    leg = new TLegend(0.58, 0.42, 0.85, 0.92);

  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  if(controlRegion == "gam"){
    leg->AddEntry(datahist, "Data","PL");
    leg->AddEntry(gamhist, "#gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "zmm"){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(vllhist, "Z#rightarrow #mu#mu","F");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "zee"){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(vllhist, "Z #rightarrow ee","F");
    leg->AddEntry(vlhist, "W #rightarrow e #nu","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "wmn"){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow #mu#mu","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "wen"){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(vlhist, "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow ee","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "topmu" and plotResonant){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(tophist_matched, "Top Resonant","F");
    leg->AddEntry(tophist_unmatched, "Top non Resonant","F");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow #mu#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "topmu" and not plotResonant){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow #mu#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "topel" and not plotResonant){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(vlhist, "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow e#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "topel" and plotResonant){
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(tophist_matched, "Top Resonant","F");
    leg->AddEntry(tophist_unmatched, "Top non Resonant","F");
    leg->AddEntry(vlhist, "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow e#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "SR"){
    leg->SetNColumns(2);
    leg->AddEntry(datahist, "Data");
    leg->AddEntry(vnnhist, "Z(#nu#nu)","F");
    leg->AddEntry(vlhist, "W(l#nu)", "F");
    leg->AddEntry(vllhist, "Z(ll)", "F");
    leg->AddEntry(tophist, "Top", "F");
    leg->AddEntry(dbhist, "Dibosons", "F");
    leg->AddEntry(qcdhist, "QCD", "F");
    TString mass = TString::Format("%.1f TeV",stof(mediatorMass)/1000); 
    leg->AddEntry(monoJhist, ("Mono-J M_{Med} = "+string(mass)).c_str(), "L");
    leg->AddEntry(monoWhist, ("Mono-W M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
    leg->AddEntry(monoZhist, ("Mono-Z M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");    
  }  

  leg->Draw("SAME");
  
  pad1->RedrawAxis("sameaxis");
  if(isLog)
    pad1->SetLogy();

  // make data/MC ratio plot
  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 = NULL;
  if(category <= 1)
    frame2 =  pad2->DrawFrame(xMin, 0.8, xMax, 1.2, "");
  else if(category > 1)
    frame2 =  pad2->DrawFrame(xMin, 0.25, xMax, 1.75, "");


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
  

  TH1* nhist = (TH1*) datahist->Clone("datahist_tot");
  TH1* unhist = (TH1*) datahist->Clone("unhist");
  TH1* dhist = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot");
  TH1* dhist_p = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot_p");

  nhist->SetStats(kFALSE);
  nhist->SetLineColor(kBlack);
  nhist->SetMarkerColor(kBlack);
  nhist->SetMarkerSize(1.0);

  // set to zero for plotting reasons of error bar and error band
  for (int i = 1; i <= dhist->GetNbinsX(); i++) dhist->SetBinError(i, 0);

  nhist->Divide(dhist);
  dhist_p->Divide(dhist);

  dhist_p->SetLineColor(0);
  dhist_p->SetMarkerColor(0);
  dhist_p->SetMarkerSize(0);
  dhist_p->SetFillColor(kGray);
  
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetFillColor(0);
  
  nhist->Draw("PE SAME");
  dhist_p->Draw("E2 SAME");
  unhist->Draw("SAME");
  nhist->Draw("PE SAME");
  
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((observable+"_"+controlRegion+".png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+".pdf").c_str());
  //  canvas->SaveAs((observable+"_"+controlRegion+".root").c_str());
}

