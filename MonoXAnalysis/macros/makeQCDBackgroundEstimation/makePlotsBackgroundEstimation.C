#include "../CMS_lumi.h"

////////////----------
void plotDistributions(TH1F* histo1, TH1F* histo2, TH1F* histo3, const string & outputDIR, const string & postfix){

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();
  
  histo1->SetLineColor(kBlack);
  histo1->SetMarkerStyle(20);
  histo1->SetMarkerSize(1.0);
  histo1->GetXaxis()->SetTitle("M_{jj} [GeV]");
  histo1->GetYaxis()->SetTitle("Events");
  histo1->GetXaxis()->SetTitleOffset(1.10);
  histo1->GetYaxis()->SetTitleOffset(1.20);
  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.1,
				   max(histo1->GetMaximum(),max(histo2->GetMaximum(),histo3->GetMaximum()))*100);
  
  histo1->Draw("EP");
  TString name (histo2->GetName());
  if(name.Contains("WJets"))
    histo2->SetFillColor(TColor::GetColor("#FAAF08"));
  else if(name.Contains("ZJets"))
    histo2->SetFillColor(TColor::GetColor("#4D975D"));
  
  TString name2 (histo3->GetName());
  if(name2.Contains("WJets"))
    histo3->SetFillColor(TColor::GetColor("#FAAF08"));
  else if(name2.Contains("ZJets"))
    histo3->SetFillColor(TColor::GetColor("#4D975D"));
  
  histo2->SetLineColor(kBlack);
  histo3->SetLineColor(kBlack);
  
  histo3->Add(histo2);
  histo3->Draw("hist same");
  histo2->Draw("hist same");
  histo1->Draw("EP same");

  TLegend* leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo1,"Data","EP");
  if(name.Contains("WJets"))
    leg->AddEntry(histo2,"W+jets","F");
  else
    leg->AddEntry(histo2,"Z+jets","F");
  if(name2.Contains("ZJets"))
    leg->AddEntry(histo3,"Z+jets","F");
  else
    leg->AddEntry(histo3,"W+jets","FL");
  leg->Draw("same");
  canvas->SetLogy();
  CMS_lumi(canvas,"35.9");
  canvas->RedrawAxis("same");
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".pdf").c_str(),"pdf");
  delete canvas;
}

void plotTransfer(TH1F* histo, TGraph* graph, const string & outputDIR, const string & postfix){

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();
  
  histo->SetLineColor(kBlack);
  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.0);
  histo->GetXaxis()->SetTitle("M_{jj} [GeV]");
  histo->GetYaxis()->SetTitle("Transfer factor");
  histo->GetXaxis()->SetTitleOffset(1.10);
  histo->GetYaxis()->SetTitleOffset(1.20);
  
  histo->Draw("EP");
  graph->SetLineColor(kRed);
  graph->SetLineWidth(2);
  graph->Draw("Lsame");

  canvas->SetLogy();
  CMS_lumi(canvas,"35.9");
  canvas->RedrawAxis("same");
  canvas->SaveAs((outputDIR+"/transfer_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/transfer_"+postfix+".pdf").c_str(),"pdf");
  delete canvas;
  

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

  TLegend* leg = new TLegend(0.6,0.65,0.9,0.9);
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


////////////----------
void makePlotsBackgroundEstimation(string inputFileName, string outputDIR){

  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");

  TH1F* dataRegion_A = (TH1F*) inputFile->Get("Data/histo_Data_mjj_regionA");
  TH1F* dataRegion_B = (TH1F*) inputFile->Get("Data/histo_Data_mjj_regionB");

  TH1F* wjetsRegion_A = (TH1F*) inputFile->Get("WJets_MC/histo_WJets_mjj_regionA");
  TH1F* wjetsRegion_B = (TH1F*) inputFile->Get("WJets_MC/histo_WJets_mjj_regionB");

  TH1F* zjetsRegion_A = (TH1F*) inputFile->Get("ZJets_MC/histo_ZJets_mjj_regionA");
  TH1F* zjetsRegion_B = (TH1F*) inputFile->Get("ZJets_MC/histo_ZJets_mjj_regionB");

  // to understand data subtraction
  plotDistributions(dataRegion_A,wjetsRegion_A,zjetsRegion_A,outputDIR,"regionA");
  plotDistributions(dataRegion_B,wjetsRegion_B,zjetsRegion_B,outputDIR,"regionB");
  
  TH1F* transferQCD_DB = (TH1F*) inputFile->Get("Transfer_QCD/transferQCD_DB_0");
  TH1F* transferQCD_CA = (TH1F*) inputFile->Get("Transfer_QCD/transferQCD_CA_0");
  TGraph* transferQCD_DB_graph = (TGraph*) inputFile->Get("Transfer_QCD/Graph_from_transferQCD_DB_0_graph"); 
  TGraph* transferQCD_CA_graph = (TGraph*) inputFile->Get("Transfer_QCD/Graph_from_transferQCD_CA_0_graph"); 
  
  plotTransfer(transferQCD_DB,transferQCD_DB_graph,outputDIR,"regionBD");
  plotTransfer(transferQCD_CA,transferQCD_CA_graph,outputDIR,"regionAC");

  TH1F* qcdMCRegion_C = (TH1F*) inputFile->Get("QCD_MC/histo_QCD_mjj_regionC");
  TH1F* qcdMCRegion_D = (TH1F*) inputFile->Get("QCD_MC/histo_QCD_mjj_regionD");
  TH1F* estimationRegion_C = (TH1F*) inputFile->Get("Estimation/histoQCD_estimatedInC_0_graph");
  TH1F* estimationRegion_D = (TH1F*) inputFile->Get("Estimation/histoQCD_estimatedInD_0_graph");

  plotComparison(qcdMCRegion_C,estimationRegion_C,outputDIR,"closureC");
  plotComparison(qcdMCRegion_D,estimationRegion_D,outputDIR,"estimationSR");

  // systematic uncertainties
  vector<TH1F*> estimationRegion_D_binByBinUp;
  vector<TH1F*> estimationRegion_D_binByBinDw;

  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_0_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_1_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_2_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_3_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_4_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_5_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_6_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_7_Up"));
  estimationRegion_D_binByBinUp.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_8_Up"));

  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_0_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_1_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_2_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_3_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_4_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_5_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_6_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_7_Dw"));
  estimationRegion_D_binByBinDw.push_back((TH1F*) inputFile->Get("BinByBin/histoQCD_estimatedInD_0_Bin_8_Dw"));

  TH1F* estimationRegion_D_bin_by_bin_up =  (TH1F*) estimationRegion_D->Clone("estimationRegion_D_bin_by_bin_up");
  TH1F* estimationRegion_D_bin_by_bin_dw =  (TH1F*) estimationRegion_D->Clone("estimationRegion_D_bin_by_bin_dw");
  for(int iBin = 0; iBin < estimationRegion_D->GetNbinsX(); iBin++){
    estimationRegion_D_bin_by_bin_up->SetBinContent(iBin+1,estimationRegion_D_binByBinUp.at(iBin)->GetBinContent(iBin+1));
    estimationRegion_D_bin_by_bin_dw->SetBinContent(iBin+1,estimationRegion_D_binByBinDw.at(iBin)->GetBinContent(iBin+1));
  }

  estimationRegion_D->SetLineColor(kBlack);
  estimationRegion_D->SetMarkerColor(kBlack);
  estimationRegion_D->SetMarkerStyle(20);
  estimationRegion_D->SetMarkerSize(1);

  estimationRegion_D_bin_by_bin_up->SetLineColor(kBlack);
  estimationRegion_D_bin_by_bin_dw->SetLineColor(kBlack);
  estimationRegion_D_bin_by_bin_up->SetLineWidth(2);
  estimationRegion_D_bin_by_bin_dw->SetLineWidth(2);

  // make closure test uncertainty band 
  TH1F* estimationRegion_D_closureTest_up = (TH1F*) estimationRegion_D->Clone("estimationRegion_D_closureTest_up");
  TH1F* estimationRegion_D_closureTest_dw = (TH1F*) estimationRegion_D->Clone("estimationRegion_D_closureTest_dw");

  estimationRegion_D_closureTest_up->SetBinContent(1,estimationRegion_D_closureTest_up->GetBinContent(1)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(2,estimationRegion_D_closureTest_up->GetBinContent(2)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(3,estimationRegion_D_closureTest_up->GetBinContent(3)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(4,estimationRegion_D_closureTest_up->GetBinContent(4)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(5,estimationRegion_D_closureTest_up->GetBinContent(5)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(6,estimationRegion_D_closureTest_up->GetBinContent(6)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(7,estimationRegion_D_closureTest_up->GetBinContent(7)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(8,estimationRegion_D_closureTest_up->GetBinContent(8)*1.50);
  estimationRegion_D_closureTest_up->SetBinContent(9,estimationRegion_D_closureTest_up->GetBinContent(9)*1.50);
  estimationRegion_D_closureTest_dw->SetBinContent(1,estimationRegion_D_closureTest_dw->GetBinContent(1)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(2,estimationRegion_D_closureTest_dw->GetBinContent(2)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(3,estimationRegion_D_closureTest_dw->GetBinContent(3)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(4,estimationRegion_D_closureTest_dw->GetBinContent(4)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(5,estimationRegion_D_closureTest_dw->GetBinContent(5)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(6,estimationRegion_D_closureTest_dw->GetBinContent(6)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(7,estimationRegion_D_closureTest_dw->GetBinContent(7)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(8,estimationRegion_D_closureTest_dw->GetBinContent(8)*0.50);
  estimationRegion_D_closureTest_dw->SetBinContent(9,estimationRegion_D_closureTest_dw->GetBinContent(9)*0.50);

  estimationRegion_D_closureTest_up->SetLineColor(kRed);
  estimationRegion_D_closureTest_dw->SetLineColor(kRed);
  estimationRegion_D_closureTest_up->SetLineWidth(2);
  estimationRegion_D_closureTest_dw->SetLineWidth(2);

  TH1F* estimationRegion_D_zjet_up = (TH1F*) inputFile->Get("ZJets/histoQCD_estimatedInD_zjet_up_0");
  TH1F* estimationRegion_D_zjet_dw = (TH1F*) inputFile->Get("ZJets/histoQCD_estimatedInD_zjet_dw_0");
  TH1F* estimationRegion_D_wjet_up = (TH1F*) inputFile->Get("WJets/histoQCD_estimatedInD_wjet_up_0");
  TH1F* estimationRegion_D_wjet_dw = (TH1F*) inputFile->Get("WJets/histoQCD_estimatedInD_wjet_dw_0");

  estimationRegion_D_zjet_up->SetLineColor(kBlue);
  estimationRegion_D_zjet_dw->SetLineColor(kBlue);
  estimationRegion_D_wjet_up->SetLineColor(kCyan+1);
  estimationRegion_D_wjet_dw->SetLineColor(kCyan+1);
  estimationRegion_D_wjet_up->SetLineStyle(7);
  estimationRegion_D_wjet_dw->SetLineStyle(7);
  estimationRegion_D_zjet_up->SetLineWidth(2);
  estimationRegion_D_zjet_dw->SetLineWidth(2);
  estimationRegion_D_wjet_up->SetLineWidth(2);
  estimationRegion_D_wjet_dw->SetLineWidth(2);
  
  // make final plot with systematics
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);  
  canvas->cd();
 
  estimationRegion_D->GetXaxis()->SetTitle("");
  estimationRegion_D->GetXaxis()->SetNdivisions(505);
  estimationRegion_D->GetXaxis()->SetLabelSize(0);
  estimationRegion_D->GetYaxis()->SetTitle("Events");
  estimationRegion_D->GetYaxis()->SetTitleOffset(1.20);
  estimationRegion_D->GetYaxis()->SetRangeUser(0.001,estimationRegion_D_zjet_dw->GetMaximum()*100);
  canvas->SetLogy();
  estimationRegion_D->Draw("P");
  estimationRegion_D_bin_by_bin_up->Draw("hist same");
  estimationRegion_D_bin_by_bin_dw->Draw("hist same");
  
  estimationRegion_D_zjet_up->Draw("hist same");
  estimationRegion_D_zjet_dw->Draw("hist same");  
  estimationRegion_D_wjet_up->Draw("hist same");
  estimationRegion_D_wjet_dw->Draw("hist same");

  estimationRegion_D_closureTest_up->Draw("hist same");
  estimationRegion_D_closureTest_dw->Draw("hist same");
  estimationRegion_D->Draw("Psame");


  TLegend* leg = new TLegend(0.6,0.55,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(estimationRegion_D,"Central Prediction","P");
  leg->AddEntry(estimationRegion_D_bin_by_bin_up,"Bin-By-Bin statistical","L");
  leg->AddEntry(estimationRegion_D_closureTest_up,"Closure test","L");
  leg->AddEntry(estimationRegion_D_zjet_up,"Z+jets Up/Dw","L");
  leg->AddEntry(estimationRegion_D_wjet_dw,"W+jets Up/Dw","L");
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

  TH1F* ratio_zjet_up = (TH1F*) estimationRegion_D_zjet_up->Clone("ratio_zjet_up");
  ratio_zjet_up->Divide(estimationRegion_D);
  TH1F* ratio_zjet_dw = (TH1F*) estimationRegion_D_zjet_dw->Clone("ratio_zjet_dw");
  ratio_zjet_dw->Divide(estimationRegion_D);

  TH1F* ratio_wjet_up = (TH1F*) estimationRegion_D_wjet_up->Clone("ratio_wjet_up");
  ratio_wjet_up->Divide(estimationRegion_D);
  TH1F* ratio_wjet_dw = (TH1F*) estimationRegion_D_wjet_dw->Clone("ratio_wjet_dw");
  ratio_wjet_dw->Divide(estimationRegion_D);

  TH1F* ratio_closure_up = (TH1F*) estimationRegion_D_closureTest_up->Clone("ratio_closure_up");
  ratio_closure_up->Divide(estimationRegion_D);
  TH1F* ratio_closure_dw = (TH1F*) estimationRegion_D_closureTest_dw->Clone("ratio_closure_dw");
  ratio_closure_dw->Divide(estimationRegion_D);

  TH1F* ratio_bin_by_bin_up =  (TH1F*) estimationRegion_D->Clone("ratio_bin_by_bin_up");
  TH1F* ratio_bin_by_bin_dw =  (TH1F*) estimationRegion_D->Clone("ratio_bin_by_bin_dw");
  for(int iBin = 0; iBin < estimationRegion_D->GetNbinsX(); iBin++){
    ratio_bin_by_bin_up->SetBinContent(iBin+1,estimationRegion_D_binByBinUp.at(iBin)->GetBinContent(iBin+1));
    ratio_bin_by_bin_dw->SetBinContent(iBin+1,estimationRegion_D_binByBinDw.at(iBin)->GetBinContent(iBin+1));
  }
  ratio_bin_by_bin_up->Divide(estimationRegion_D);
  ratio_bin_by_bin_dw->Divide(estimationRegion_D);

  ratio_zjet_up->GetXaxis()->SetTitle("M_{jj} [GeV]");
  ratio_zjet_up->GetYaxis()->SetTitle("DD/MC");
  ratio_zjet_up->GetYaxis()->SetTitleSize(0.04);
  ratio_zjet_up->GetXaxis()->SetTitleSize(0.04);
  ratio_zjet_up->GetYaxis()->SetTitleOffset(1.20);
  ratio_zjet_up->GetXaxis()->SetTitleOffset(1.20);
  ratio_zjet_up->GetYaxis()->SetLabelSize(0.035);
  ratio_zjet_up->GetXaxis()->SetNdivisions(505);
  ratio_zjet_up->GetYaxis()->SetNdivisions(505);
  ratio_zjet_up->GetYaxis()->SetRangeUser(0,2.5);
  ratio_zjet_up->Draw("hist");
  ratio_zjet_dw->Draw("hist same");
  ratio_wjet_up->Draw("hist same");
  ratio_wjet_dw->Draw("hist same");
  ratio_closure_up->Draw("hist same");
  ratio_closure_dw->Draw("hist same");
  ratio_bin_by_bin_up->Draw("hist same");
  ratio_bin_by_bin_dw->Draw("hist same");
  
  canvas->SaveAs((outputDIR+"/distribution_estimationSR_systematics.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/distribution_estimationSR_systeamtics.pdf").c_str(),"pdf");

  if(canvas) delete canvas;
    
}
