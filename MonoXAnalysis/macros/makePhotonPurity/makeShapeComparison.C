//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../CMS_lumi.h"   // with this header I can't compile macro with .L makeShapeComparison.C.C++

using namespace std;

static vector<int> ptBins = {175,200,225,250,280,320,370,420,1000};  // must be consistent with input file content

//==============================================

void checkHistogramInFile(TH1F* h, const string& histoName, const string& fileName) {

  if (!h || h == NULL) {
    cout << "Error: histogram '" << histoName << "' not found in file '" << fileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

}


//==============================================

void drawPlotBkg(TH1F* h1, TH1F* h2, 
	      const string& xAxisName = "", 
	      const string& canvasTitle = "", 
	      const string& outputDIR = "./", 
	      const float& lumi = 36.2,
	      const int& ptmin = 0,
	      const int& ptmax = 1000
	      ) {

  gROOT->SetBatch(kTRUE);

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  TCanvas *canvas = new TCanvas("canvas","",700,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.06);

  // for ratio plot
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetLeftMargin(0.15);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  // this for ratio plot as well
  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);
  
  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);

  h1->SetLineWidth(2);
  h2->SetLineWidth(2);

  h1->GetXaxis()->SetTitle(xAxisName.c_str());
  h1->GetXaxis()->SetLabelSize(0);

  h1->GetYaxis()->SetTitle("a.u");
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(1.03);

  h1->SetTitle("");
  h1->SetStats(0);

  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
 
  Double_t maximumYaxisValue = -100000.0;
  Double_t minimumYaxisValue = 100000.0;
  maximumYaxisValue = TMath::Max(h1->GetBinContent(h1->GetMaximumBin()),h2->GetBinContent(h2->GetMaximumBin()));
  minimumYaxisValue = TMath::Min(h1->GetBinContent(h1->GetMinimumBin()),h2->GetBinContent(h2->GetMinimumBin()));

  h1->SetMaximum(maximumYaxisValue * 1.2);
  h1->SetMinimum(0);

  h1->Draw("EP");
  h2->Draw("HIST SAME");  // error in ratio plot

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.90);
  leg->AddEntry((TObject*)0,Form("%d < p_{T} < %d",ptmin,ptmax),"");
  leg->AddEntry(h1,"background data","PLE");
  leg->AddEntry(h2,"background QCD","L");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend
  leg->SetFillColor(0);

  // set x axis range to have histogram cover (almost) whole length on the right side
  // can be useful when the canvas's x axis range is decided a priori (e.g. 0 - 1000) but the plot starts being empty at some point (e.g. at 500)
  // Int_t lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  // h1->GetXaxis()->SetRange(0,lastVisibleBin);

  canvas->RedrawAxis("sameaxis");

  CMS_lumi(canvas,Form("%.1f",lumi));

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.0,2.0);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle("Ratio");
  frame->GetXaxis()->SetTitle("Photon Isolation [GeV]");
  
  TH1F* ratio = (TH1F*) h1->Clone("ratio");
  TH1F* denFit_noerr = (TH1F*) h2->Clone("denFit_noerr");
  TH1F* denFit = (TH1F*) h2->Clone("denFit");
  for(int iBin = 1; iBin < denFit->GetNbinsX()+1; iBin++)
    denFit_noerr->SetBinError(iBin,0.);
  
  ratio->Divide(denFit_noerr);
  denFit->Divide(denFit_noerr);
  denFit->SetFillColor(kGray);
  frame->Draw();
  ratio->Draw("EPsame");
  denFit->Draw("E2same");

  TF1* line = new TF1("line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  // canvas->SaveAs( (outputDIR + canvasTitle + ".pdf").c_str() );
  // canvas->SaveAs( (outputDIR + canvasTitle + ".png").c_str() );

  canvas->cd();

  // now log scale plot
  h1->SetMaximum( maximumYaxisValue * 100 );
  h1->SetMinimum( 0.8 * TMath::Max( 0.001, h1->GetMinimum() ) );  

  // lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  // //if (lastVisibleBin < h3->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h3->FindLastBinAbove(h1->GetMinimum());
  // h1->GetXaxis()->SetRange(0,lastVisibleBin);

  canvas->Update();
  canvas->SetLogy();
  canvas->SaveAs( (outputDIR + canvasTitle + "_logy.pdf").c_str() );
  canvas->SaveAs( (outputDIR + canvasTitle + "_logy.png").c_str() );
  canvas->SetLogy(0);

  delete canvas;
  delete leg;  

}


//==============================================

void drawPlotSig(TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4,
	      const string& xAxisName = "", 
	      const string& canvasTitle = "", 
	      const string& outputDIR = "./", 
	      const float& lumi = 36.2,
	      const int& ptmin = 0,
	      const int& ptmax = 1000
	      ) {

  gROOT->SetBatch(kTRUE);

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  TCanvas *canvas = new TCanvas("canvas","",700,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.06);
  
  // for ratio plot
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetLeftMargin(0.15);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);

  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kBlue);
  h4->SetLineColor(kGreen+2);


  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  h4->SetLineWidth(2);

  h1->GetXaxis()->SetTitle(xAxisName.c_str());
  h1->GetXaxis()->SetLabelSize(0);

  h1->GetYaxis()->SetTitle("a.u");
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(1.03);

  h1->SetTitle("");
  h1->SetStats(0);

  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1);

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
  h4->Scale(1./h4->Integral());
  
  Double_t maximumYaxisValue = -100000.0;
  Double_t minimumYaxisValue = 100000.0;
  maximumYaxisValue = TMath::Max(h1->GetBinContent(h1->GetMaximumBin()),h2->GetBinContent(h2->GetMaximumBin()));
  maximumYaxisValue = TMath::Max(maximumYaxisValue,h3->GetBinContent(h3->GetMaximumBin()));
  maximumYaxisValue = TMath::Max(maximumYaxisValue,h4->GetBinContent(h4->GetMaximumBin()));
  minimumYaxisValue = TMath::Min(h1->GetBinContent(h1->GetMinimumBin()),h2->GetBinContent(h2->GetMinimumBin()));
  minimumYaxisValue = TMath::Min(minimumYaxisValue,h3->GetBinContent(h3->GetMinimumBin()));
  minimumYaxisValue = TMath::Min(minimumYaxisValue,h4->GetBinContent(h4->GetMinimumBin()));

  h1->SetMaximum(maximumYaxisValue * 1.2);
  h1->SetMinimum(0);

  h1->Draw("EP");
  h2->Draw("EP SAME");  
  h3->Draw("HIST SAME");  
  h4->Draw("HIST SAME");  

  TLegend *leg = new TLegend(0.4,0.6,0.9,0.90);
  leg->AddEntry((TObject*)0,Form("%d < p_{T} < %d",ptmin,ptmax),"");
  leg->AddEntry(h1,"signal data RND04","PLE");
  leg->AddEntry(h2,"signal data RND08","PLE");
  leg->AddEntry(h3,"signal #gamma+jets matched","L");
  leg->AddEntry(h4,"signal #gamma+jets RND04","L");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend
  leg->SetFillColor(0);

  // set x axis range to have histogram cover (almost) whole length on the right
  // can be useful when the canvas's x axis range is decided a priori (e.g. 0 - 1000) but the plot starts being empty at some point (e.g. at 500)
  // Int_t lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h3->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h3->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h4->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h4->FindLastBinAbove(h1->GetMinimum());
  // h1->GetXaxis()->SetRange(0,lastVisibleBin);

  canvas->RedrawAxis("sameaxis");

  CMS_lumi(canvas,Form("%.1f",lumi));

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.0,2.0);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle("data/MC");  // title based on the ratios below
  frame->GetXaxis()->SetTitle("Photon Isolation [GeV]");
  
  // do ratio of data / MC matched using data RND04 and RND08
  TH1F* ratio1 = (TH1F*) h1->Clone("ratio1");
  TH1F* ratio2 = (TH1F*) h2->Clone("ratio2");
  TH1F* denFit_noerr = (TH1F*) h3->Clone("denFit_noerr");
  TH1F* denFit = (TH1F*) h3->Clone("denFit");
  for(int iBin = 1; iBin < denFit->GetNbinsX()+1; iBin++)
    denFit_noerr->SetBinError(iBin,0.);
  
  ratio1->Divide(denFit_noerr);
  ratio2->Divide(denFit_noerr);
  denFit->Divide(denFit_noerr);
  denFit->SetFillColor(kGray);
  frame->Draw();
  ratio1->Draw("EPsame");
  ratio2->Draw("EPsame");
  denFit->Draw("E2same");

  TF1* line = new TF1("line","1",ratio1->GetXaxis()->GetBinLowEdge(1),ratio1->GetXaxis()->GetBinLowEdge(ratio1->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio1->Draw("EPsame");
  ratio2->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");


  // canvas->SaveAs( (outputDIR + canvasTitle + ".pdf").c_str() );
  // canvas->SaveAs( (outputDIR + canvasTitle + ".png").c_str() );

  canvas->cd();

  // now log scale plot
  h1->SetMaximum( maximumYaxisValue * 100 );
  h1->SetMinimum( 0.8 * TMath::Max( 0.000001, h1->GetMinimum() ) );  

  // lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h3->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h3->FindLastBinAbove(h1->GetMinimum());
  // if (lastVisibleBin < h4->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h4->FindLastBinAbove(h1->GetMinimum());
  // h1->GetXaxis()->SetRange(0,lastVisibleBin);

  canvas->RedrawAxis("sameaxis");

  canvas->SetLogy();
  canvas->SaveAs( (outputDIR + canvasTitle + "_logy.pdf").c_str() );
  canvas->SaveAs( (outputDIR + canvasTitle + "_logy.png").c_str() );
  canvas->SetLogy(0);

  delete canvas;
  delete leg;  

}


//=====================

void makeShapeComparison(const string& inputDIR = "./", const float& lumi = 36.2) {

  // assume the input file is called PhotonPurityFitResult.root, but can choose its location with inputDIR (default is current directory)

  string outputDIR = inputDIR + "shapeComparison/";
  if (outputDIR != "./") system(("mkdir -p "+outputDIR).c_str());
  string bkgDIR = "/background/";
  string sigDIR = "/signal/";
  system(("mkdir -p "+outputDIR+bkgDIR).c_str());
  system(("mkdir -p "+outputDIR+sigDIR).c_str());

  // style                                                                                                                                                                  
  gROOT->SetBatch(kTRUE);
  //setTDRStyle();

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()

  string inputfilename = inputDIR + "PhotonPurityFitResult.root"; 
  string infileDir = "Templates/";

  // open file
  TFile* inputfile = TFile::Open(inputfilename.c_str());
  if (!inputfile || !inputfile->IsOpen()) {
    cout << "******************************* " << endl;
    cout << "Error opening file \"" << inputfile << "\".\nApplication will be terminated." << endl;
    cout << "******************************* "<< endl;
    exit(EXIT_FAILURE);
  }


  TH1F* hist = NULL; // dummy pointer to get histogram in file, the we clone it 
  int nptBins = ( (int)ptBins.size() ) - 1; 

  // loop on pt bins
  for (int ipt = 0; ipt < nptBins; ipt++) {

    cout << "bin " << ipt << " : [" << ptBins[ipt] << ", " << ptBins[ipt+1]  << "]" << endl;

    TH1F *hBkg_data = NULL;
    TH1F *hBkg_mcQCD = NULL;
    TH1F *hSig_data_RND04 = NULL;
    TH1F *hSig_data_RND08 = NULL;
    TH1F *hSig_mcGJETS_matched = NULL;
    TH1F *hSig_mcGJETS_RND04 = NULL;
   
    string hname_Bkg_data = string( (char*)Form("backgroundTemplate_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string hname_Bkg_mcQCD = string( (char*)Form("backgroundTemplate_qcd_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string hname_Sig_data_RND04 = string( (char*)Form("signalTemplateRND04_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string hname_Sig_data_RND08 = string( (char*)Form("signalTemplateRND08_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string hname_Sig_mcGJETS_matched = string( (char*)Form("signalTemplate_gjets_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string hname_Sig_mcGJETS_RND04 = string( (char*)Form("signalTemplateRND04_gjets_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    
    hist = (TH1F*) inputfile->Get((infileDir+hname_Bkg_data).c_str());
    checkHistogramInFile(hist,hname_Bkg_data,inputfilename);
    hBkg_data = (TH1F*) hist->Clone();
    //cout << hBkg_data->GetName() << endl;

    hist = (TH1F*) inputfile->Get((infileDir+hname_Bkg_mcQCD).c_str());
    checkHistogramInFile(hist,hname_Bkg_mcQCD,inputfilename);
    hBkg_mcQCD = (TH1F*) hist->Clone();
    //cout << hBkg_mcQCD->GetName() << endl;

    hist = (TH1F*) inputfile->Get((infileDir+hname_Sig_data_RND04).c_str());
    checkHistogramInFile(hist,hname_Sig_data_RND04,inputfilename);
    hSig_data_RND04 = (TH1F*) hist->Clone();
    //cout << hSig_data_RND04->GetName() << endl;

    hist = (TH1F*) inputfile->Get((infileDir+hname_Sig_data_RND08).c_str());
    checkHistogramInFile(hist,hname_Sig_data_RND08,inputfilename);
    hSig_data_RND08 = (TH1F*) hist->Clone();
    //cout << hSig_data_RND08->GetName() << endl;

    hist = (TH1F*) inputfile->Get((infileDir+hname_Sig_mcGJETS_matched).c_str());
    checkHistogramInFile(hist,hname_Sig_mcGJETS_matched,inputfilename);
    hSig_mcGJETS_matched = (TH1F*) hist->Clone();
    //cout << hSig_mcGJETS_matched->GetName() << endl;

    hist = (TH1F*) inputfile->Get((infileDir+hname_Sig_mcGJETS_RND04).c_str());
    checkHistogramInFile(hist,hname_Sig_mcGJETS_RND04,inputfilename);
    hSig_mcGJETS_RND04 = (TH1F*) hist->Clone();
    //cout << hSig_mcGJETS_RND04->GetName() << endl;


    //drawPlotBkg(TH1F* h1, TH1F* h2, const string& xAxisName = "", const string& canvasTitle = "", const string& outputDIR = "./")
    drawPlotBkg(hBkg_data, hBkg_mcQCD, 
	     hBkg_data->GetXaxis()->GetTitle(), 
	     Form("comp_bkg_data_mcQCD_pt%dTo%d",ptBins[ipt],ptBins[ipt+1]), 
	     outputDIR+bkgDIR,
	     lumi,
	     ptBins[ipt],ptBins[ipt+1]);

    drawPlotSig(hSig_data_RND04,hSig_data_RND08,hSig_mcGJETS_matched,hSig_mcGJETS_RND04, 
	     hSig_data_RND04->GetXaxis()->GetTitle(), 
	     Form("comp_sig_data_mcGJETS_pt%dTo%d",ptBins[ipt],ptBins[ipt+1]), 
	     outputDIR+sigDIR,
	     lumi,
	     ptBins[ipt],ptBins[ipt+1]);


    cout << endl;
    

  }
  // end loop on pt bins

} 
