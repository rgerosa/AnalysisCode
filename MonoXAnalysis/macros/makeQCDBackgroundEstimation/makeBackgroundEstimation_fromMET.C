#include "../CMS_lumi.h"

void makeRatioFit(TH1* ratio, const string & outputDIR, const string & xAxisTitle, const string & yAxisTitle){

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  canvas->SetRightMargin(0.06);
  ratio->GetXaxis()->SetTitle(xAxisTitle.c_str());
  ratio->GetYaxis()->SetTitle(yAxisTitle.c_str());
  ratio->GetXaxis()->SetNdivisions(505);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerSize(1.0);
  ratio->SetMarkerStyle(20);
  ratio->Draw("EP");

  canvas->SetLogy(0);

  CMS_lumi(canvas,"35.9");

  ratio->GetYaxis()->SetRangeUser(0,ratio->GetMaximum()*2.0);
  ratio->GetYaxis()->SetTitleOffset(1.10);
  ratio->GetXaxis()->SetTitleOffset(1.10);
  ratio->GetYaxis()->SetLabelSize(0.035);

  TF1* fitExp = new TF1("fitExp","[0]*TMath::Exp([1]*x)",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
  TFitResultPtr result = ratio->Fit(fitExp,"RMS");
  TF1* func = ratio->GetFunction("fitExp");

  func->SetLineColor(kBlue);
  func->SetLineWidth(2);
  ratio->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+Form("%s",ratio->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+Form("%s",ratio->GetName())+".pdf").c_str(),"pdf");

  if(canvas) delete canvas;

}

void plotComparison(TH1* histo1, TH1* histo2, const string & outputDIR, const string & postfix){

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->cd();

  TH1F* ratio =  (TH1F*) histo2->Clone("ratio");
  ratio->Divide(histo1);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerSize(1);

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2);
  histo2->GetXaxis()->SetTitle("");
  histo2->GetXaxis()->SetNdivisions(505);
  histo2->GetXaxis()->SetLabelSize(0);
  histo2->GetYaxis()->SetTitle("Events");
  histo2->GetYaxis()->SetTitleOffset(1.20);
  histo2->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),histo2->GetMinimum())*0.1,
                                   max(histo1->GetMaximum(),histo2->GetMaximum())*100);

  if(histo1->GetMinimum() == 0 or histo2->GetMinimum() == 0)
    histo2->GetYaxis()->SetRangeUser(0.0001,max(histo1->GetMaximum(),histo2->GetMaximum())*100);


  histo2->Draw("hist");

  histo1->SetMarkerColor(kBlack);
  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1);
  histo1->Draw("EPsame");

  TLegend* leg = new TLegend(0.55,0.67,0.9,0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo1,"QCD-MC Prediction","EP");
  leg->AddEntry(histo2,"Data Driven Prediction","L");
  leg->Draw("same");

  CMS_lumi(canvas,"35.9");

  canvas->SetLogy();
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  pad2->Draw();
  pad2->cd();

  ratio->Draw("EP");
  ratio->GetYaxis()->SetRangeUser(min(0.,ratio->GetMinimum())*0.8,ratio->GetMaximum()*1.3);
  ratio->GetXaxis()->SetTitle("M_{jj} [GeV]");
  ratio->GetYaxis()->SetTitle("DD/MC");
  ratio->GetYaxis()->SetTitleSize(0.04);
  ratio->GetXaxis()->SetTitleSize(0.04);
  ratio->GetYaxis()->SetTitleOffset(1.20);
  ratio->GetXaxis()->SetTitleOffset(1.20);
  ratio->GetYaxis()->SetLabelSize(0.035);
  ratio->GetXaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetNdivisions(505);

  TF1* line = new TF1("line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->Draw("Lsame");
  canvas->SaveAs((outputDIR+"/distribution_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/distribution_"+postfix+".pdf").c_str(),"pdf");

  if(canvas) delete canvas;
  if(line) delete line;
}


void makeBackgroundEstimation_fromMET(string inputTemplateFile, string jetHTfile, string outputDIR){

  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);

  ///////// take input from control region
  TFile* input_QCDCR = TFile::Open(inputTemplateFile.c_str(),"READ");

  TH1F* data_QCDCR   = (TH1F*) input_QCDCR->Get("QCD/datahistqcd_mjj");
  TH1F* qcdbkg_QCDCR = (TH1F*) input_QCDCR->Get("QCD/qbkghistqcd_mjj");
  TH1F* zvvbkg_QCDCR = (TH1F*) input_QCDCR->Get("QCD/vnnbkghistqcd_mjj");
  TH1F* wjetbkg_QCDCR = (TH1F*) input_QCDCR->Get("QCD/vlbkghistqcd_mjj");

  zvvbkg_QCDCR->Add((TH1F*) input_QCDCR->Get("QCD/ewkbkgzhistqcd_mjj"));
  wjetbkg_QCDCR->Add((TH1F*) input_QCDCR->Get("QCD/ewkbkgwhistqcd_mjj"));

  TH1F* minorbkg_QCDCR = (TH1F*) input_QCDCR->Get("QCD/vllbkghistqcd_mjj");

  minorbkg_QCDCR->Add((TH1F*) input_QCDCR->Get("QCD/tbkghistqcd_mjj"));
  minorbkg_QCDCR->Add((TH1F*) input_QCDCR->Get("QCD/dbkghistqcd_mjj"));
  
  TH1F* qcdbkg_QCDCR_den = (TH1F*) input_QCDCR->Get("QCD/qbkghistqcd_mjj_v4");

  TFile* inputJetHT      = TFile::Open(jetHTfile.c_str(),"READ");
  TH1F* qcdbkg_QCDCR_num = (TH1F*) inputJetHT->Get("QCD_MC_v2/histo_QCD_v2_mjj_regionD");
  TH1F* qcdbkg_QCDCR_SR  = (TH1F*) inputJetHT->Get("QCD_MC/histo_QCD_mjj_regionD");

  /////////
  TH1F* transferFactor = (TH1F*) qcdbkg_QCDCR_num->Clone("transferFactor");
  transferFactor->Divide(qcdbkg_QCDCR_den);
  TGraph* transferFactor_graph = new TGraph(transferFactor);  

  vector<TH1F*> transferFactor_binUp;
  vector<TH1F*> transferFactor_binDw;
  vector<TGraph*> transferFactor_graph_binUp;
  vector<TGraph*> transferFactor_graph_binDw;

  for(int iBin = 0; iBin < transferFactor->GetNbinsX(); iBin++){
    transferFactor_binUp.push_back((TH1F*) transferFactor->Clone(Form("transferFactor_bin%d_Up",iBin+1)));
    transferFactor_binDw.push_back((TH1F*) transferFactor->Clone(Form("transferFactor_bin%d_Dw",iBin+1)));
    transferFactor_binUp.back()->SetBinContent(iBin+1,transferFactor->GetBinContent(iBin+1)+transferFactor->GetBinError(iBin+1));
    transferFactor_binDw.back()->SetBinContent(iBin+1,transferFactor->GetBinContent(iBin+1)-transferFactor->GetBinError(iBin+1));
    transferFactor_graph_binUp.push_back(new TGraph(transferFactor_binUp.back()));
    transferFactor_graph_binDw.push_back(new TGraph(transferFactor_binDw.back()));
  }

  makeRatioFit(transferFactor,outputDIR,"M_{jj} [GeV]","Transfer factor");
  
  /////// estimation central value
  TH1F* estimationInSR = (TH1F*) qcdbkg_QCDCR_SR->Clone("estimationInSR");
  estimationInSR->Reset();
  for(int iBin = 0; iBin < estimationInSR->GetNbinsX(); iBin++){
    float binCenter = estimationInSR->GetBinCenter(iBin+1);
    estimationInSR->SetBinContent(iBin+1,transferFactor_graph->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-wjetbkg_QCDCR->GetBinContent(iBin+1)-zvvbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)));
  }

  ///////// systematic uncertainties
  float scale_zjet = 0.2;
  float scale_wjet = 0.2;

  TH1F* estimationInSR_zjet_up = (TH1F*) qcdbkg_QCDCR_SR->Clone("estimationInSR_zjet_up");
  TH1F* estimationInSR_zjet_dw = (TH1F*) qcdbkg_QCDCR_SR->Clone("estimationInSR_zjet_dw");
  TH1F* estimationInSR_wjet_up = (TH1F*) qcdbkg_QCDCR_SR->Clone("estimationInSR_wjet_up");
  TH1F* estimationInSR_wjet_dw = (TH1F*) qcdbkg_QCDCR_SR->Clone("estimationInSR_wjet_dw");

  for(int iBin = 0; iBin < estimationInSR->GetNbinsX(); iBin++){
    float binCenter = estimationInSR->GetBinCenter(iBin+1);
    estimationInSR_zjet_up->SetBinContent(iBin+1,transferFactor_graph->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-wjetbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-
											zvvbkg_QCDCR->GetBinContent(iBin+1)*(1+scale_zjet)));
    estimationInSR_zjet_dw->SetBinContent(iBin+1,transferFactor_graph->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-wjetbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-
											zvvbkg_QCDCR->GetBinContent(iBin+1)*(1-scale_zjet)));
    estimationInSR_wjet_up->SetBinContent(iBin+1,transferFactor_graph->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-zvvbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-
											wjetbkg_QCDCR->GetBinContent(iBin+1)*(1+scale_wjet)));
    estimationInSR_wjet_dw->SetBinContent(iBin+1,transferFactor_graph->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-zvvbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-
											wjetbkg_QCDCR->GetBinContent(iBin+1)*(1-scale_wjet)));
  }

  estimationInSR_zjet_up->SetLineColor(kBlue);
  estimationInSR_zjet_dw->SetLineColor(kBlue);
  estimationInSR_wjet_up->SetLineColor(kCyan+1);
  estimationInSR_wjet_dw->SetLineColor(kCyan+1);
  estimationInSR_wjet_up->SetLineStyle(7);
  estimationInSR_wjet_dw->SetLineStyle(7);
  estimationInSR_zjet_up->SetLineWidth(2);
  estimationInSR_zjet_dw->SetLineWidth(2);
  estimationInSR_wjet_up->SetLineWidth(2);
  estimationInSR_wjet_dw->SetLineWidth(2);

 ///////
  vector<TH1F*> estimationInSR_binByBinUp;
  vector<TH1F*> estimationInSR_binByBinDw;
  for(int iBin = 0; iBin < estimationInSR->GetNbinsX(); iBin++){
    estimationInSR_binByBinUp.push_back((TH1F*) estimationInSR->Clone(Form("estimationInSR_bin%d_Up",iBin)));
    estimationInSR_binByBinDw.push_back((TH1F*) estimationInSR->Clone(Form("estimationInSR_bin%d_Dw",iBin)));
    float binCenter = estimationInSR->GetBinCenter(iBin+1);
    int bin         = qcdbkg_QCDCR_num->FindBin(binCenter);
    estimationInSR_binByBinUp.back()->SetBinContent(iBin+1,transferFactor_graph_binUp.at(bin-1)->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-wjetbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-zvvbkg_QCDCR->GetBinContent(iBin+1)));
    estimationInSR_binByBinDw.back()->SetBinContent(iBin+1,transferFactor_graph_binDw.at(bin-1)->Eval(binCenter)*(data_QCDCR->GetBinContent(iBin+1)-wjetbkg_QCDCR->GetBinContent(iBin+1)-minorbkg_QCDCR->GetBinContent(iBin+1)-zvvbkg_QCDCR->GetBinContent(iBin+1)));
  }

  //////
  TH1F* estimationInSR_bin_by_bin_up =  (TH1F*) estimationInSR->Clone("estimationInSR_bin_by_bin_up");
  TH1F* estimationInSR_bin_by_bin_dw =  (TH1F*) estimationInSR->Clone("estimationInSR_bin_by_bin_dw");
  for(int iBin = 0; iBin < estimationInSR->GetNbinsX(); iBin++){
    estimationInSR_bin_by_bin_up->SetBinContent(iBin+1,estimationInSR_binByBinUp.at(iBin)->GetBinContent(iBin+1));
    estimationInSR_bin_by_bin_dw->SetBinContent(iBin+1,estimationInSR_binByBinDw.at(iBin)->GetBinContent(iBin+1));
  }
  
  estimationInSR->SetLineColor(kBlack);
  estimationInSR->SetMarkerColor(kBlack);
  estimationInSR->SetMarkerStyle(20);
  estimationInSR->SetMarkerSize(1);
  
  estimationInSR_bin_by_bin_up->SetLineColor(kBlack);
  estimationInSR_bin_by_bin_dw->SetLineColor(kBlack);
  estimationInSR_bin_by_bin_up->SetLineWidth(2);
  estimationInSR_bin_by_bin_dw->SetLineWidth(2);						    

  // make closure test uncertainty band                                                                                                                                                               
  TH1F* estimationInSR_closureTest_up = (TH1F*) estimationInSR->Clone("estimationInSR_closureTest_up");
  TH1F* estimationInSR_closureTest_dw = (TH1F*) estimationInSR->Clone("estimationInSR_closureTest_dw");

  estimationInSR_closureTest_up->SetBinContent(1,estimationInSR_closureTest_up->GetBinContent(1)*1.50);
  estimationInSR_closureTest_up->SetBinContent(2,estimationInSR_closureTest_up->GetBinContent(2)*1.50);
  estimationInSR_closureTest_up->SetBinContent(3,estimationInSR_closureTest_up->GetBinContent(3)*1.50);
  estimationInSR_closureTest_up->SetBinContent(4,estimationInSR_closureTest_up->GetBinContent(4)*1.50);
  estimationInSR_closureTest_up->SetBinContent(5,estimationInSR_closureTest_up->GetBinContent(5)*1.50);
  estimationInSR_closureTest_up->SetBinContent(6,estimationInSR_closureTest_up->GetBinContent(6)*1.50);
  estimationInSR_closureTest_up->SetBinContent(7,estimationInSR_closureTest_up->GetBinContent(7)*1.50);
  estimationInSR_closureTest_up->SetBinContent(8,estimationInSR_closureTest_up->GetBinContent(8)*1.50);
  estimationInSR_closureTest_up->SetBinContent(9,estimationInSR_closureTest_up->GetBinContent(9)*1.50);
  estimationInSR_closureTest_dw->SetBinContent(1,estimationInSR_closureTest_dw->GetBinContent(1)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(2,estimationInSR_closureTest_dw->GetBinContent(2)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(3,estimationInSR_closureTest_dw->GetBinContent(3)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(4,estimationInSR_closureTest_dw->GetBinContent(4)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(5,estimationInSR_closureTest_dw->GetBinContent(5)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(6,estimationInSR_closureTest_dw->GetBinContent(6)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(7,estimationInSR_closureTest_dw->GetBinContent(7)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(8,estimationInSR_closureTest_dw->GetBinContent(8)*0.50);
  estimationInSR_closureTest_dw->SetBinContent(9,estimationInSR_closureTest_dw->GetBinContent(9)*0.50);

  estimationInSR_closureTest_up->SetLineColor(kRed);
  estimationInSR_closureTest_dw->SetLineColor(kRed);
  estimationInSR_closureTest_up->SetLineWidth(2);
  estimationInSR_closureTest_dw->SetLineWidth(2);

  // Plot
  plotComparison(qcdbkg_QCDCR_SR,estimationInSR,outputDIR,"estimationSR");
  
  //// Plot
  // make final plot with systematics                                                                                                                                                                 
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->cd();

  estimationInSR->GetXaxis()->SetTitle("");
  estimationInSR->GetXaxis()->SetNdivisions(505);
  estimationInSR->GetXaxis()->SetLabelSize(0);
  estimationInSR->GetYaxis()->SetLabelSize(0.035);
  estimationInSR->GetYaxis()->SetTitle("Events");
  estimationInSR->GetYaxis()->SetTitleOffset(1.20);
  estimationInSR->GetYaxis()->SetRangeUser(0.001,estimationInSR_bin_by_bin_up->GetMaximum()*1000);
  canvas->SetLogy();
  estimationInSR->Draw("P");

  estimationInSR_closureTest_up->Draw("hist same");
  estimationInSR_closureTest_dw->Draw("hist same");

  estimationInSR_bin_by_bin_up->Draw("hist same");
  estimationInSR_bin_by_bin_dw->Draw("hist same");

  estimationInSR_zjet_up->Draw("hist same");
  estimationInSR_zjet_dw->Draw("hist same");
  estimationInSR_wjet_up->Draw("hist same");
  estimationInSR_wjet_dw->Draw("hist same");
  estimationInSR->Draw("Psame");


  TLegend* leg = new TLegend(0.6,0.55,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(estimationInSR,"Central Prediction","P");
  leg->AddEntry(estimationInSR_bin_by_bin_up,"Bin-By-Bin statistical","L");
  leg->AddEntry(estimationInSR_closureTest_up,"Closure test","L");
  leg->AddEntry(estimationInSR_zjet_up,"Z+jets Up/Dw","L");
  leg->AddEntry(estimationInSR_wjet_dw,"W+jets Up/Dw","L");
  leg->Draw("same");
  CMS_lumi(canvas,"35.9");

  canvas->SetLogy();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  pad2->Draw();
  pad2->cd();

  TH1F* ratio_zjet_up = (TH1F*) estimationInSR_zjet_up->Clone("ratio_zjet_up");
  ratio_zjet_up->Divide(estimationInSR);
  TH1F* ratio_zjet_dw = (TH1F*) estimationInSR_zjet_dw->Clone("ratio_zjet_dw");
  ratio_zjet_dw->Divide(estimationInSR);

  TH1F* ratio_wjet_up = (TH1F*) estimationInSR_wjet_up->Clone("ratio_wjet_up");
  ratio_wjet_up->Divide(estimationInSR);
  TH1F* ratio_wjet_dw = (TH1F*) estimationInSR_wjet_dw->Clone("ratio_wjet_dw");
  ratio_wjet_dw->Divide(estimationInSR);

  TH1F* ratio_closure_up = (TH1F*) estimationInSR_closureTest_up->Clone("ratio_closure_up");
  ratio_closure_up->Divide(estimationInSR);
  TH1F* ratio_closure_dw = (TH1F*) estimationInSR_closureTest_dw->Clone("ratio_closure_dw");
  ratio_closure_dw->Divide(estimationInSR);

  TH1F* ratio_bin_by_bin_up =  (TH1F*) estimationInSR_bin_by_bin_up->Clone("ratio_bin_by_bin_up");
  TH1F* ratio_bin_by_bin_dw =  (TH1F*) estimationInSR_bin_by_bin_dw->Clone("ratio_bin_by_bin_dw");
  ratio_bin_by_bin_up->Divide(estimationInSR);
  ratio_bin_by_bin_dw->Divide(estimationInSR);
  

  ratio_zjet_up->GetXaxis()->SetTitle("M_{jj} [GeV]");
  ratio_zjet_up->GetYaxis()->SetTitle("DD/MC");
  ratio_zjet_up->GetYaxis()->SetTitleSize(0.04);
  ratio_zjet_up->GetYaxis()->SetLabelSize(0.035);
  ratio_zjet_up->GetXaxis()->SetTitleSize(0.04);
  ratio_zjet_up->GetYaxis()->SetTitleOffset(1.20);
  ratio_zjet_up->GetXaxis()->SetTitleOffset(1.20);
  ratio_zjet_up->GetYaxis()->SetLabelSize(0.035);
  ratio_zjet_up->GetXaxis()->SetNdivisions(505);
  ratio_zjet_up->GetYaxis()->SetNdivisions(505);
  ratio_zjet_up->GetYaxis()->SetRangeUser(0,3.0);
  ratio_zjet_up->Draw("hist");
  ratio_zjet_dw->Draw("hist same");
  ratio_wjet_up->Draw("hist same");
  ratio_wjet_dw->Draw("hist same");
  ratio_closure_up->Draw("hist same");
  ratio_closure_dw->Draw("hist same");
  ratio_bin_by_bin_up->Draw("hist same");
  ratio_bin_by_bin_dw->Draw("hist same");

  canvas->SaveAs((outputDIR+"/distribution_estimationSR_systematics.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/distribution_estimationSR_systematics.pdf").c_str(),"pdf");

  if(canvas) delete canvas;

  // make a reduced output file for the QCD in SR                                                                                                                                                      
  TFile* outputFile = new TFile((outputDIR+"/templates_QCD_DD.root").c_str(),"RECREATE");
  outputFile->cd();

  estimationInSR->Write("template_QCD_SR_fromDD");
  outputFile->mkdir("BinByBin");
  outputFile->cd("BinByBin");

  int iBin = 0;
  for(auto hist : estimationInSR_binByBinUp){
    hist->Write(Form("template_QCD_SR_fromDD_bin_%d_statUp",iBin));
    iBin++;
  }
  iBin = 0;
  for(auto hist : estimationInSR_binByBinDw){
    hist->Write(Form("template_QCD_SR_fromDD_bin_%d_statDown",iBin));
    iBin++;
  }

  outputFile->cd();
  outputFile->cd();
  outputFile->mkdir("ZJets");
  outputFile->cd("ZJets");
  estimationInSR_zjet_up->Write("template_QCD_SR_ZJetsUp");
  estimationInSR_zjet_dw->Write("template_QCD_SR_ZJetsDown");
  outputFile->cd();
  outputFile->mkdir("WJets");
  outputFile->cd("WJets");
  estimationInSR_wjet_up->Write("template_QCD_SR_WJetsUp");
  estimationInSR_wjet_dw->Write("template_QCD_SR_WJetsDown");
  outputFile->cd();
  outputFile->Close();
}
