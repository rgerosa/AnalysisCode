#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include <RooStats/RooStatsUtils.h>
#include <RooStats/HLFactory.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>
#include <RooRandom.h>

#include "makePhotonPurityPDFs.h"
#include "../CMS_lumi.h"

using namespace std;

static vector<int> ptBins = {175,200,225,250,280,320,370,420,480,1000};

static bool debug = false;

//===========================================================0


void checkWorkspaceInFile(RooWorkspace * ws, const string& wsName, const string& fileName) {

  if (!ws || ws == NULL) {
    cout << "Error: workspace '" << wsName << "' not found in file '" << fileName << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

}

//===========================================================0

void fillTGraphFromWorkspace(TGraphAsymmErrors* gr, 
			     double& meanPurityValue, double& meanPurityUncertainty,
			     TFile* inputfile = NULL, const string& inputfilename = "", 
			     const string& wsName = "", 
			     const int& ipt = -1, const int& ptMin = 0, const int& ptMax = 0) 
{

    RooWorkspace* ws = (RooWorkspace*) inputfile->Get( (wsName).c_str() ); 
    checkWorkspaceInFile(ws, wsName, inputfilename);
    RooRealVar* purity = (RooRealVar*) ws->var("photonPurity");
    purity->Print();

    gr->SetPoint(ipt,double((ptMax+ptMin)/2.), purity->getVal());
    gr->SetPointError(ipt,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());

    // get error as mean value between up and down uncertainty
    double unc = 0.5 * ( fabs(purity->getErrorLo()) + purity->getErrorHi() );
    if (unc == 0.) {
      cout << "WARNING: uncertainty for pt [" << ptMin << "," << ptMax << "] and workspace " << wsName << " is 0. Setting it to 0.0015" << endl;
      unc = 0.0015;
    }

    // compute weight of mean average and add term in the numerator of weighted average formula
    double weight = 1. / (unc * unc);

    meanPurityValue += purity->getVal() * weight;
    meanPurityUncertainty += weight;

    if (debug) {
      cout << "purity->getVal() = " << purity->getVal() << "   weight = " << weight << endl;
      cout << "meanPurityValue = " << meanPurityValue << "   meanPurityUncertainty = " << meanPurityUncertainty << endl;
    }

}

//===========================================================0


void meanPhotonPurity(const string& inputDIR = "./", const string& outputDIR = "./", const double& lumi = 36.2, const bool& makeFitBasedOnlyOnTemplates = false) {

  if (outputDIR!= "./") system(("mkdir -p "+outputDIR).c_str());

  gROOT->ProcessLine(".L makePhotonPurityPDFs.cc+");
  gSystem->Load("makePhotonPurityPDFs_cc.so");

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()                                           

  string inputfilename = inputDIR + "PhotonPurityFitResult.root";

  // open file                                                                                                                                  
  TFile* inputfile = TFile::Open(inputfilename.c_str());
  if (!inputfile || !inputfile->IsOpen()) {
    cout << "******************************* " << endl;
    cout << "Error opening file \"" << inputfile << "\".\nApplication will be terminated." << endl;
    cout << "******************************* "<< endl;
    exit(EXIT_FAILURE);
  }

  TGraphAsymmErrors* rnd04 = new TGraphAsymmErrors();  // this is the nominal value
  TGraphAsymmErrors* rnd08 = new TGraphAsymmErrors();
  TGraphAsymmErrors* gjets = new TGraphAsymmErrors();
  TGraphAsymmErrors* qcd = new TGraphAsymmErrors();
  TGraphAsymmErrors* rnd04_altSig = new TGraphAsymmErrors();
  TGraphAsymmErrors* rnd04_altBkg = new TGraphAsymmErrors();

  TGraphAsymmErrors* finalPurity = new TGraphAsymmErrors();

  vector<TGraphAsymmErrors*> vec;
  vec.push_back(rnd04);
  vec.push_back(rnd08);
  vec.push_back(gjets);
  vec.push_back(qcd);
  vec.push_back(rnd04_altSig);
  vec.push_back(rnd04_altBkg);

  int nptBins = ( (int) ptBins.size() ) - 1;

  // loop on pt bins                                                                                                             
  for (int ipt = 0; ipt < nptBins; ipt++) {

    int ptMin = ptBins[ipt];
    int ptMax = ptBins[ipt+1];

    double meanPurity = 0;
    double uncertaintyMeanPurity = 0;

    cout << endl; 
    cout << endl; 
    cout << "pt : [min,max] = " << ptMin << " - " << ptMax << endl;

    string wsName_rnd04 = string( (char*)Form("wsRND04_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string wsName_rnd08 = string( (char*)Form("wsRND08_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string wsName_gjets = string( (char*)Form("ws_gjets_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string wsName_qcd = string( (char*)Form("ws_qcd_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string wsName_altSig = string( (char*)Form("wsRND04_altSig_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;
    string wsName_altBkg = string( (char*)Form("wsRND04_altBkg_data_pt_%d_%d",ptBins[ipt],ptBins[ipt+1]) ) ;


    fillTGraphFromWorkspace(rnd04, meanPurity, uncertaintyMeanPurity, inputfile, inputfilename, wsName_rnd04, ipt, ptMin, ptMax);
    fillTGraphFromWorkspace(rnd08, meanPurity, uncertaintyMeanPurity, inputfile, inputfilename, wsName_rnd08, ipt, ptMin, ptMax);
    fillTGraphFromWorkspace(gjets, meanPurity, uncertaintyMeanPurity, inputfile, inputfilename, wsName_gjets, ipt, ptMin, ptMax);
    fillTGraphFromWorkspace(qcd, meanPurity, uncertaintyMeanPurity,   inputfile, inputfilename, wsName_qcd,   ipt, ptMin, ptMax);
    fillTGraphFromWorkspace(rnd04_altSig, meanPurity, uncertaintyMeanPurity, inputfile, inputfilename, wsName_altSig, ipt, ptMin, ptMax);
    fillTGraphFromWorkspace(rnd04_altBkg, meanPurity, uncertaintyMeanPurity, inputfile, inputfilename, wsName_altBkg, ipt, ptMin, ptMax);

    // at this point uncertaintyMeanPurity is the sum of weights in the weighted average formula (see inside fillTGraphFromWorkspace(...) )
    // while meanPurity is the numerator in the same formula, which must be divided by the sum of weights

    uncertaintyMeanPurity = 1./uncertaintyMeanPurity;  // do inversion here
    meanPurity *= uncertaintyMeanPurity;               // compute weighted average completing the formula
    uncertaintyMeanPurity = sqrt(uncertaintyMeanPurity);  // uncertainty on the average is the square root of the inverse of the sum of weights

    cout << "pt : [min,max] = " << ptMin << " - " << ptMax << "    purity --> " << meanPurity << " +/- " << uncertaintyMeanPurity << endl;

    finalPurity->SetPoint(ipt, double((ptMax+ptMin)/2.), meanPurity);
    finalPurity->SetPointError(ipt, double((ptMax-ptMin)/2.), double((ptMax-ptMin)/2.), uncertaintyMeanPurity, uncertaintyMeanPurity);

  }

  cout << endl;
  cout << endl;

  // plot on canvas

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetRightMargin(0.06);

  finalPurity->SetLineColor(kBlack);
  finalPurity->SetMarkerColor(kBlack);
  finalPurity->SetMarkerStyle(20);
  finalPurity->SetMarkerSize(1);
  finalPurity->GetYaxis()->SetRangeUser(0.5,1.2);
  finalPurity->GetXaxis()->SetTitle("photon p_{T} [GeV]");
  finalPurity->GetYaxis()->SetTitle("photon purity");
  finalPurity->SetFillColor(kOrange+1);
  finalPurity->Draw("AP");
  finalPurity->Draw("E2same");
  finalPurity->Draw("Psame");

  CMS_lumi(canvas,Form("%.1f",lumi));    

  canvas->SaveAs((outputDIR+"/photonPurityFinal_weightedAverage.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/photonPurityFinal_weightedAverage.pdf").c_str(),"pdf");

  TGraphAsymmErrors* gr = (TGraphAsymmErrors*) inputfile->Get( "purity" );
  if (!gr || gr == NULL) {
    cout << "Error: TGraphAsymmErrors 'purity' not found in file '" << inputfilename << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  // to get stat error we use rnd04 graph,; to get syst use "purity"
  // simply we take the purity mean value in each pt bin as the weighted average
  TGraphAsymmErrors* finalPurity_totUnc = new TGraphAsymmErrors(*finalPurity);
  TGraphAsymmErrors* finalPurity_statOnly_rnd04 = new TGraphAsymmErrors(*finalPurity);

  for (int i = 0; i < finalPurity_totUnc->GetN(); i++) {

    finalPurity_totUnc->SetPointEYlow(i, gr->GetErrorYlow(i));    
    finalPurity_totUnc->SetPointEYhigh(i, gr->GetErrorYhigh(i));    
    finalPurity_statOnly_rnd04->SetPointEYlow(i, rnd04->GetErrorYlow(i));    
    finalPurity_statOnly_rnd04->SetPointEYhigh(i, rnd04->GetErrorYhigh(i));    

  } 


  TCanvas* canvas2 = new TCanvas("canvas2","",600,650);
  canvas2->cd();
  canvas2->SetTickx(1);
  canvas2->SetTicky(1);
  canvas2->cd();
  canvas2->SetRightMargin(0.06);

  finalPurity_statOnly_rnd04->SetLineColor(kBlack);
  finalPurity_statOnly_rnd04->SetMarkerColor(kBlack);
  finalPurity_statOnly_rnd04->SetMarkerStyle(20);
  finalPurity_statOnly_rnd04->SetMarkerSize(1);

  finalPurity_totUnc->SetLineColor(kBlack);
  finalPurity_totUnc->SetMarkerColor(kBlack);
  finalPurity_totUnc->SetMarkerStyle(20);
  finalPurity_totUnc->SetMarkerSize(1);
  finalPurity_statOnly_rnd04->GetYaxis()->SetRangeUser(0.5,1.2);
  finalPurity_statOnly_rnd04->GetXaxis()->SetTitle("photon p_{T} [GeV]");
  finalPurity_statOnly_rnd04->GetYaxis()->SetTitle("photon purity");
  finalPurity_totUnc->SetFillColor(kOrange+1);
  finalPurity_statOnly_rnd04->Draw("AP");
  finalPurity_totUnc->Draw("E2same");
  finalPurity_statOnly_rnd04->Draw("Psame");

  CMS_lumi(canvas2,Form("%.1f",lumi));    

  canvas2->SaveAs((outputDIR+"/photonPurityFinal_weightedAverage_err.png").c_str(),"png");
  canvas2->SaveAs((outputDIR+"/photonPurityFinal_weightedAverage_err.pdf").c_str(),"pdf");

  double* purity = finalPurity_totUnc->GetY();
  double* stat_high = finalPurity_statOnly_rnd04->GetEYhigh();
  double* stat_low= finalPurity_statOnly_rnd04->GetEYlow();
  double* tot_high = finalPurity_totUnc->GetEYhigh();
  double* tot_low= finalPurity_totUnc->GetEYlow();
  double syst_high = 0;
  double syst_low = 0;

  cout << endl;
  cout << endl;
  cout << "PURITY TABLE WITH UNCERTAINTIES" << endl;
  cout << "-------------------------------------------------" << endl;
  for (int i = 0; i < finalPurity_totUnc->GetN(); i++) {

    syst_high = sqrt(tot_high[i] * tot_high[i] - stat_high[i] * stat_high[i]);    
    syst_low = sqrt(tot_low[i] * tot_low[i] - stat_low[i] * stat_low[i]);    
    cout << " purity --> " << purity[i] << " + " << stat_high[i] << " - " << stat_low[i] << " (stat.)  + " << syst_high << " - " << syst_low << " (syst.)" << endl;  

  } 
  cout << "-------------------------------------------------" << endl;

  // save two new TGraphAsymmErrors in another file
  string outputfilename = outputDIR + "finalPurity_weightedAverage.root";
  TFile* outputFile = new TFile((outputfilename).c_str(),"RECREATE");
  if (!outputFile || !outputFile->IsOpen()) {
    cout << "******************************* " << endl;
    cout << "Error opening file \"" << outputfilename << "\".\nApplication will be terminated." << endl;
    cout << "******************************* "<< endl;
    exit(EXIT_FAILURE);
  }
  outputFile->cd();
  finalPurity_totUnc->Write("purity_weightedAverage_totUnc");
  finalPurity_statOnly_rnd04->Write("purity_weightedAverage_statOnly_rnd04");
  outputFile->Close();
  delete outputFile;

}
