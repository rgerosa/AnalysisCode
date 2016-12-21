#ifndef FitPurityNewId_h
#define FitPurityNewId_h

//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TChain.h>
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
#include <TStyle.h>
#include <TString.h>
#include <TVirtualFitter.h>

//gSystem->AddIncludePath("-I/wherever/RooFit/is/include");
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
#include <iomanip> 
#include <algorithm>  // to use the "reverse" function to reverse the order in the array
#include <Rtypes.h> // to use kColor

/* using namespace RooFit; */
/* using namespace RooStats; */

void AddData(RooWorkspace*, TTree*);
void runfits(RooWorkspace* w,std::string region, std::string phid, std::string ptcut);
void PlotBkgTemplates(RooWorkspace* w);
void AddSignalMCRND(RooWorkspace* w);
void AddSignalMC(RooWorkspace* w);
void AddBkgMC(RooWorkspace* w);
void ReduceSignalMCRND(RooWorkspace* w, TString ptrange, TString suffix);
void ReduceBkgMC(RooWorkspace* w, TString ptrange, TString suffix);
void PlotSignalTemplateMCs(RooWorkspace* w);
void PlotBkgTemplateMCs(RooWorkspace* w);
void PlotSignalTemplates(RooWorkspace* w);
void PlotBkgTemplates(RooWorkspace* w);
RooRealVar* ModelFit(RooWorkspace* w);
void PlotAllSystUnc(RooWorkspace* w,std::string region, std::string id);
void ComputeAllSystUnc(RooWorkspace* w, std::string region, std::string id);
void ComputeSystUnc(RooWorkspace*w, std::string region, std::string id, RooHistPdf* genSig, RooHistPdf* genBkg, RooHistPdf* fitSig, RooHistPdf* fitBkg, TString Tsig, TString Tbkg);
//RooFitResult*  ModelFit(RooWorkspace*);
void AddSignalTemplateRND04(RooWorkspace* w, TTree* tree);
void AddSignalTemplateRND08(RooWorkspace* w, TTree* tree);
void AddBkgTemplate(RooWorkspace* w, TTree* tree);
void AddSignalMCRNDRC08(RooWorkspace* w);
void runAllfits(std::string phid);


#endif 
