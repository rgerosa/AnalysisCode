#include "../CMS_lumi.h"

enum Sample {zmm, wmn, sig};

void plotHistogram(TCanvas* canvas, TH1* histoNominal, TH1* histoUp, TH1* histoDw, string outputDIR, string postfix, string observableLatex){

  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  histoNominal->SetLineColor(kBlack);
  histoNominal->SetLineWidth(2);

  histoUp->SetLineColor(kRed);
  histoUp->SetLineWidth(2);

  histoDw->SetLineColor(kBlue);
  histoDw->SetLineWidth(2);

  histoNominal->GetYaxis()->SetTitle("Events");
  histoNominal->GetXaxis()->SetTitleSize(0);
  histoNominal->GetXaxis()->SetLabelSize(0);
  histoNominal->GetYaxis()->SetTitleOffset(1.45);

  histoNominal->GetYaxis()->SetRangeUser(min(histoNominal->GetMinimum(),min(histoDw->GetMinimum(),histoUp->GetMinimum()))*0.9,
					 max(histoNominal->GetMaximum(),max(histoDw->GetMaximum(),histoUp->GetMaximum()))*1.1);
  
  histoNominal->Draw("hist");
  histoUp->Draw("hist same");
  histoDw->Draw("hist same");

  CMS_lumi(canvas,"35.9");

  TLegend leg (0.70,0.70,0.92,0.92);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(histoNominal,"Nominal","L");
  leg.AddEntry(histoUp,"JES UP","L");
  leg.AddEntry(histoDw,"JES DW","L");
  leg.Draw("same");

  TPad* pad2 = new TPad(("pad2"+postfix).c_str(),"pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* ratio_up = (TH1*) histoUp->Clone("ratio_up");
  ratio_up->SetLineColor(kRed);
  ratio_up->SetLineWidth(2);
  ratio_up->Divide(histoNominal);

  ratio_up->GetYaxis()->SetTitle("Ratio");
  ratio_up->GetYaxis()->SetTitleOffset(1.45);
  ratio_up->GetXaxis()->SetTitle(observableLatex.c_str());   
  ratio_up->GetXaxis()->SetNdivisions(510);
  ratio_up->GetYaxis()->SetNdivisions(504);
  ratio_up->Draw("hist");

  TH1* ratio_dw = (TH1*) histoDw->Clone("ratio_dw");
  ratio_dw->SetLineColor(kBlue);
  ratio_dw->SetLineWidth(2);
  ratio_dw->Divide(histoNominal);

  ratio_dw->GetYaxis()->SetTitle("Ratio");
  ratio_dw->GetXaxis()->SetTitle(observableLatex.c_str());
  ratio_dw->Draw("hist same");
  
  ///
  ratio_up->GetYaxis()->SetRangeUser(min(ratio_dw->GetMinimum(),ratio_up->GetMinimum())*0.9,max(ratio_up->GetMaximum(),ratio_dw->GetMaximum())*1.1);

  TF1 line ("line","1",ratio_up->GetBinLowEdge(1),ratio_up->GetBinLowEdge(ratio_up->GetNbinsX()+1));
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw("hist same");

  canvas->RedrawAxis("sameaxis");
  
  canvas->cd();

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  histoNominal->GetYaxis()->SetRangeUser(min(histoNominal->GetMinimum(),min(histoDw->GetMinimum(),histoUp->GetMinimum()))*0.01,
					 max(histoNominal->GetMaximum(),max(histoDw->GetMaximum(),histoUp->GetMaximum()))*100);
  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/"+postfix+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+"_log.pdf").c_str(),"pdf");
  
  canvas->SetLogy(0);

}


////
void plotRatios(TCanvas* canvas, TH1* ratioNominal, TH1* ratioUp, TH1* ratioDw, string outputDIR, string postfix, string observableLatex){

  canvas->cd();
  canvas->SetLogy(0);

  ratioNominal->SetLineColor(kBlack);
  ratioNominal->SetLineWidth(2);
  ratioNominal->SetMarkerSize(1);
  ratioNominal->SetMarkerStyle(20);

  ratioNominal->GetXaxis()->SetTitleSize(0);
  ratioNominal->GetXaxis()->SetLabelSize(0);
  ratioNominal->GetYaxis()->SetTitle("Transfer Factor");

  ratioNominal->GetYaxis()->SetRangeUser(min(ratioNominal->GetMinimum(),min(ratioDw->GetMinimum(),ratioUp->GetMinimum()))*0.75,
  					 max(ratioNominal->GetMaximum(),max(ratioDw->GetMaximum(),ratioUp->GetMaximum()))*1.5);

  
  ratioNominal->Draw("EP");

  TGraphAsymmErrors* graphError = new TGraphAsymmErrors();  

  for(int iBin = 0; iBin < ratioNominal->GetNbinsX(); iBin++){
    graphError->SetPoint(iBin,ratioNominal->GetBinCenter(iBin+1),ratioNominal->GetBinContent(iBin+1));
    if(ratioDw->GetBinContent(iBin+1) < ratioNominal->GetBinContent(iBin+1) and ratioUp->GetBinContent(iBin+1) > ratioNominal->GetBinContent(iBin+1))
      graphError->SetPointError(iBin,ratioNominal->GetBinWidth(iBin+1)/2,ratioNominal->GetBinWidth(iBin+1)/2,
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioDw->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)),
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioUp->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)));    
    else if(ratioDw->GetBinContent(iBin+1) > ratioNominal->GetBinContent(iBin+1) and ratioUp->GetBinContent(iBin+1) < ratioNominal->GetBinContent(iBin+1))
      graphError->SetPointError(iBin,ratioNominal->GetBinWidth(iBin+1)/2,ratioNominal->GetBinWidth(iBin+1)/2,
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioUp->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)),
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioDw->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)));    
    else if(ratioDw->GetBinContent(iBin+1) > ratioNominal->GetBinContent(iBin+1) and ratioUp->GetBinContent(iBin+1) > ratioNominal->GetBinContent(iBin+1))
      graphError->SetPointError(iBin,ratioNominal->GetBinWidth(iBin+1)/2,ratioNominal->GetBinWidth(iBin+1)/2,0,
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioUp->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)+
				     pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioDw->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)));    
    else if(ratioDw->GetBinContent(iBin+1) < ratioNominal->GetBinContent(iBin+1) and ratioUp->GetBinContent(iBin+1) < ratioNominal->GetBinContent(iBin+1))
      graphError->SetPointError(iBin,ratioNominal->GetBinWidth(iBin+1)/2,ratioNominal->GetBinWidth(iBin+1)/2,
				sqrt(pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioUp->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)+
				     pow(ratioNominal->GetBinError(iBin+1),2)+pow(ratioDw->GetBinContent(iBin+1)-ratioNominal->GetBinContent(iBin+1),2)),0);    
    
  }
  

  graphError->SetFillColor(33);
  graphError->SetFillStyle(1001);
  graphError->Draw("2same");
  ratioNominal->Draw("EPsame");
  CMS_lumi(canvas,"35.9");
  
  
  TLegend leg (0.65,0.65,0.92,0.92);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(ratioNominal,"Ratio Nominal","L");
  leg.AddEntry(graphError,"JES + Stat Unc.","F");
  leg.Draw("same");


  TPad* pad2 = new TPad(("pad2"+postfix).c_str(),"pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  canvas->cd();
  pad2->Draw();
  pad2->cd();


  TH1* ratio_up = (TH1*) ratioUp->Clone("ratio_up");
  ratio_up->SetLineColor(kRed);
  ratio_up->SetLineWidth(2);
  ratio_up->Divide(ratioNominal);

  ratio_up->GetYaxis()->SetTitle("JES SYS");
  ratio_up->GetYaxis()->SetTitleOffset(1.45);
  ratio_up->GetXaxis()->SetTitle(observableLatex.c_str());   
  ratio_up->GetXaxis()->SetNdivisions(510);
  ratio_up->GetYaxis()->SetNdivisions(504);
  ratio_up->Draw("hist");

  TH1* ratio_dw = (TH1*) ratioDw->Clone("ratio_dw");
  ratio_dw->SetLineColor(kBlue);
  ratio_dw->SetLineWidth(2);
  ratio_dw->Divide(ratioNominal);

  ratio_dw->GetYaxis()->SetTitle("Ratio");
  ratio_dw->GetXaxis()->SetTitle(observableLatex.c_str());
  ratio_dw->Draw("hist same");
  
  ///
  ratio_up->GetYaxis()->SetRangeUser(min(ratio_dw->GetMinimum(),ratio_up->GetMinimum())*0.9,max(ratio_up->GetMaximum(),ratio_dw->GetMaximum())*1.1);

  TF1 line ("line","1",ratio_up->GetBinLowEdge(1),ratio_up->GetBinLowEdge(ratio_up->GetNbinsX()+1));
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw("hist same");

  canvas->RedrawAxis("sameaxis");
  
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  canvas->SetLogy();

  ratioNominal->GetYaxis()->SetRangeUser(min(ratioNominal->GetMinimum(),min(ratioDw->GetMinimum(),ratioUp->GetMinimum()))*0.01,
  					 max(ratioNominal->GetMaximum(),max(ratioDw->GetMaximum(),ratioUp->GetMaximum()))*100);

  canvas->SaveAs((outputDIR+"/"+postfix+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+"_log.pdf").c_str(),"pdf");

  canvas->SetLogy(0);

}

////////////////
void makeTransferFactorJesEffect(string fileNominal, string fileUp, string fileDw, string outputDIR, string observable, string observableLatex, Sample sample){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  
  TFile* inputFileNominal = TFile::Open(fileNominal.c_str());
  TFile* inputFileUp = TFile::Open(fileUp.c_str());
  TFile* inputFileDw = TFile::Open(fileDw.c_str());

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);  

  string controlregion;
  if (sample == Sample::zmm)
    controlregion = "zmm";
  else if (sample == Sample::wmn)
    controlregion = "wmn";
  else if (sample == Sample::sig)
    controlregion = "zwj";

  // numerator
  TH1* numeratorNominal = (TH1*) inputFileNominal->Get(("nhist_"+controlregion+"_"+observable).c_str());
  TH1* numeratorUp = (TH1*) inputFileUp->Get(("nhist_"+controlregion+"_jesUp_"+observable).c_str());
  TH1* numeratorDw = (TH1*) inputFileDw->Get(("nhist_"+controlregion+"_jesDw_"+observable).c_str());

  // denominator
  TH1* denominatorNominal = (TH1*) inputFileNominal->Get(("dhist_"+controlregion+"_"+observable).c_str());
  TH1* denominatorUp = (TH1*) inputFileUp->Get(("dhist_"+controlregion+"_jesUp_"+observable).c_str());
  TH1* denominatorDw = (TH1*) inputFileDw->Get(("dhist_"+controlregion+"_jesDw_"+observable).c_str());

  TH1* ratioNominal = (TH1*) numeratorNominal->Clone("ratioNominal");
  ratioNominal->Divide(denominatorNominal);  
  TH1* ratioUp = (TH1*) numeratorUp->Clone("ratioUp");
  ratioUp->Divide(denominatorUp);
  TH1* ratioDw = (TH1*) numeratorDw->Clone("ratioDw");
  ratioDw->Divide(denominatorDw);

  ////
  plotHistogram(canvas,numeratorNominal,numeratorUp,numeratorDw,outputDIR,"num_"+controlregion+"_"+observable,observableLatex);
  plotHistogram(canvas,denominatorNominal,denominatorUp,denominatorDw,outputDIR,"den_"+controlregion+"_"+observable,observableLatex);
  plotRatios(canvas,ratioNominal,ratioUp,ratioDw,outputDIR,"ratio_"+controlregion+"_"+observable,observableLatex);

  /// printout
  double numeratorErr = 0;
  double numeratorInt = numeratorNominal->IntegralAndError(1,numeratorNominal->GetNbinsX(),numeratorErr);
  double numeratorJesUpErr   = 0;
  double numeratorJesUpInt   = numeratorUp->IntegralAndError(1,numeratorUp->GetNbinsX(),numeratorJesUpErr);
  double numeratorJesDwErr   = 0;
  double numeratorJesDwInt   = numeratorDw->IntegralAndError(1,numeratorDw->GetNbinsX(),numeratorJesDwErr);
  
  cout<<"Numerator: nominal "<<numeratorInt<<" \\pm "<<numeratorErr<<" up variation "<<numeratorJesUpInt<<" \\pm "<<numeratorJesUpErr<<" dw variation "<<numeratorJesDwInt<<" \\pm "<<numeratorJesDwErr<<endl;

  double denominatorErr = 0;
  double denominatorInt = denominatorNominal->IntegralAndError(1,denominatorNominal->GetNbinsX(),denominatorErr);
  double denominatorJesUpErr   = 0;
  double denominatorJesUpInt   = denominatorUp->IntegralAndError(1,denominatorUp->GetNbinsX(),denominatorJesUpErr);
  double denominatorJesDwErr   = 0;
  double denominatorJesDwInt   = denominatorDw->IntegralAndError(1,denominatorDw->GetNbinsX(),denominatorJesDwErr);
  
  cout<<"Denominator: nominal "<<denominatorInt<<" \\pm "<<denominatorErr<<" up variation "<<denominatorJesUpInt<<" \\pm "<<denominatorJesUpErr<<" dw variation "<<denominatorJesDwInt<<" \\pm "<<denominatorJesDwErr<<endl;

  double ratioErr = 0;
  double ratioInt = ratioNominal->IntegralAndError(1,ratioNominal->GetNbinsX(),ratioErr);
  double ratioJesUpErr   = 0;
  double ratioJesUpInt   = ratioUp->IntegralAndError(1,ratioUp->GetNbinsX(),ratioJesUpErr);
  double ratioJesDwErr   = 0;
  double ratioJesDwInt   = ratioDw->IntegralAndError(1,ratioDw->GetNbinsX(),ratioJesDwErr);
  
  cout<<"Ratio: nominal "<<numeratorInt/denominatorInt<<" \\pm "<<sqrt(numeratorErr*numeratorErr/(denominatorInt*denominatorInt)+denominatorErr*denominatorErr*(numeratorInt*numeratorInt)/(pow(denominatorInt,4)))
      <<" up variation "<<numeratorJesUpInt/denominatorJesUpInt<<" \\pm "<<sqrt(numeratorJesUpErr*numeratorJesUpErr/(denominatorJesUpInt*denominatorJesUpInt)+denominatorJesUpErr*denominatorJesUpErr*(numeratorJesUpInt*numeratorJesUpInt)/(pow(denominatorJesUpInt,4)))
      <<" dw variation "<<numeratorJesDwInt/denominatorJesDwInt<<" \\pm "<<sqrt(numeratorJesDwErr*numeratorJesDwErr/(denominatorJesDwInt*denominatorJesDwInt)+denominatorJesDwErr*denominatorJesDwErr*(numeratorJesDwInt*numeratorJesDwInt)/(pow(denominatorJesDwInt,4)))<<endl;

  
  
}
