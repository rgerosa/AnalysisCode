#include "makePhotonPurityUtils.h"
#include <fstream>

// to decide which data to run
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};
// photon pt bins
static vector<float> ptBins     = {175,200,225,250,280,320,370,420,480,1000};
// photon isolation info
static vector<int>   nBinPhotonIso = {30, 30, 30, 30, 30, 25, 25, 20, 20};  // size of these vectors is ptBins.size()-1
static vector<float> photonIsoMax  = {20, 20, 20, 20, 20, 20, 20, 20, 20};
static vector<float> photonIsoMin  = { 0,  0,  0,  0,  0,  0,  0,  0,  0};

///////////////////////////////////////////////////////////////
// photon isolation info with variable bin width
static bool useTestVector = false;
static vector<float> testVector = {0.0, 0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 8.0, 8.8, 9.6, 10.4, 11.2, 12.0, 12.8, 13.6, 14.4, 15.2, 16.0, 16.8, 17.6, 18.4, 19.2, 20.0};
//--------------------------------------------
static vector<float> photonIsoVariableBinWidth_lowPt = {0.0, 0.5, 0.8, 1.1, 1.4, 
							1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 
							4.8, 5.4, 6.0, 6.6, 7.2, 7.8, 8.4, 9.0, 9.6, 10.2, 10.8, 11.4, 12.0,
							12.6, 13.2, 13.8, 14.4, 15.0, 15.6, 16.2, 16.8, 17.4, 18.0, 18.6, 19.2, 19.8, 20.4, 21};//, 21.6, 
//22.4, 23.2, 24.0};
//--------------------------------------------
static float mediumPt = 279.0;
static vector<float> photonIsoVariableBinWidth_mediumPt = {0.0, 0.5, 0.8, 1.1, 1.4, 
							   1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 
							   4.8, 5.4, 6.0, 6.6, 7.2, 7.8, 8.4, 9.0, 9.6, 10.2, 10.8, 11.4, 12.0,
							   12.8, 13.6, 14.4, 15.2, 16.0, 16.8, 17.6, 18.4, 19.2, 20.0, 
							   21.0};//, 22.0, 23.0, 24.0};
//--------------------------------------------
static float highPt = 419.0;
static vector<float> photonIsoVariableBinWidth_highPt = {0.0, 0.6, 1.0, 1.4, 1.8, 2.2,
							 2.8, 3.4, 4.0, 4.6, 5.2, 5.8, 
							 6.6, 7.4, 8.2, 9.0, 9.8, 10.6, 11.4, 12.2, 13.0,
							 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0};//,
							 //6.6, 7.4, 8.2, 9.0, 
							 //10.0, 11.0, 12.5, 14.0,  
							 //15.5, 17.0, 21.0};
static float veryHighPt = 9999.0;
static vector<float> photonIsoVariableBinWidth_veryHighPt = {0.0, 0.6, 1.0, 1.4, 1.8, 2.2,
							     2.8, 3.4, 4.0, 4.6, 5.2, 5.8, 
							     6.6, 7.4, 8.2, 9.0, 9.8, 10.6, 11.4}; 

//------------------------------------------------

// useful information from previous vector
// static int nBinsPhIsoVarWidth_lowPt = (int) photonIsoVariableBinWidth_lowPt.size() -1;
// static float phIsoVarWidth_lowPt_min = photonIsoVariableBinWidth_lowPt.front();
// static float phIsoVarWidth_lowPt_max = photonIsoVariableBinWidth_lowPt.back();
// //--------------------------------------------
// static int nBinsPhIsoVarWidth_mediumPt = (int) photonIsoVariableBinWidth_mediumPt.size() -1;
// static float phIsoVarWidth_mediumPt_min = photonIsoVariableBinWidth_mediumPt.front();
// static float phIsoVarWidth_mediumPt_max = photonIsoVariableBinWidth_mediumPt.back();
// //--------------------------------------------
// static int nBinsPhIsoVarWidth_highPt = (int) photonIsoVariableBinWidth_highPt.size() -1;
// static float phIsoVarWidth_highPt_min = photonIsoVariableBinWidth_highPt.front();
// static float phIsoVarWidth_highPt_max = photonIsoVariableBinWidth_highPt.back();
//////////////////////////////////////////////////////////////

int getNbinsPhIsoVarWidth(const float pt) {

  if (useTestVector)      return ((int) testVector.size()) -1; 
  if (pt > veryHighPt)    return ((int) photonIsoVariableBinWidth_veryHighPt.size()) -1;
  else if (pt > highPt)   return ((int) photonIsoVariableBinWidth_highPt.size()) -1;
  else if (pt > mediumPt) return ((int) photonIsoVariableBinWidth_mediumPt.size()) -1;
  else                    return ((int) photonIsoVariableBinWidth_lowPt.size()) -1;

}

//========================================================

float* getPhIsoBinArray(const float pt) {

  if (useTestVector)      return testVector.data();
  if (pt > veryHighPt)    return photonIsoVariableBinWidth_veryHighPt.data();
  else if (pt > highPt)   return photonIsoVariableBinWidth_highPt.data();
  else if (pt > mediumPt) return photonIsoVariableBinWidth_mediumPt.data();
  else                    return photonIsoVariableBinWidth_lowPt.data();

}

//-------------------------------------------

vector<float> getPhIsoBinVector(const float pt) {

  if (useTestVector)      return testVector;
  if (pt > veryHighPt)    return photonIsoVariableBinWidth_veryHighPt;
  else if (pt > highPt)   return photonIsoVariableBinWidth_highPt;
  else if (pt > mediumPt) return photonIsoVariableBinWidth_mediumPt;
  else                    return photonIsoVariableBinWidth_lowPt;

}

//--------------------------------------------

//////////////////////////////////
// dividing by bin width (necessary if using variable bin width)
// careful, when taking integral must then multiply by bin width
void divideBinContentByBinWidth(TH1* h) {

  h->Sumw2();

  for (int i = 1; i < h->GetNbinsX()+1; i++) {

    h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
    h->SetBinError(i, h->GetBinError(i)/h->GetBinWidth(i));

  }

}
//////////////////////////////////


// debug mode
//static bool debug = false;  // implemented as an option
static bool saveHistograms = true;
// k-facotr file
static string kfactorFileName = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";

void makePhotonPurityFit(string inputDirectory, // directory with dataFiles (file names must match with RunEra: all matched files, and only them, are used)
			 float  lumi, // luminosity
			 string outputDIR,
			 bool   addSystematics = false,
			 string inputDirectorySignalMC = "",   // here there must be only MC signal root files
			 string inputDirectoryBackgroundMC = "",  // here there must be only MC background root files
			 bool   makeFitBasedOnlyOnTemplates = false,
			 bool   uniformIsoBinning = true,     // decide to use or not a uniform or variable width binning (bin width depends on photon pt)
			 bool   debug = false
			 ){

  system(("mkdir -p "+outputDIR).c_str());

  // style
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()

  //from twiki https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
  photonID mediumID (0.0396,0.01022,0.01400,0.441,2.725,0.0148,0.000017,2.571,0.0047); // set wp for medium id

  // bins for purity and histograms
  vector<fitPurity> dataHisto;
  vector<fitPurity> signalTemplateRND04_data;
  vector<fitPurity> signalTemplateRND08_data;
  vector<fitPurity> backgroundTemplate_data;
  vector<fitPurity> signalTemplate_gjets;
  vector<fitPurity> signalTemplateRND04_gjets;
  vector<fitPurity> backgroundTemplate_qcd;

  // variable bin width histograms, used to estimate an additional systematic uncertainty
  // vector<fitPurity> dataHisto_varBinWidth;
  // vector<fitPurity> signalTemplateRND04_varBinWidth_data;
  // vector<fitPurity> backgroundTemplate_varBinWidth_data;

  for(size_t ibin = 0; ibin < ptBins.size()-1; ibin++){

    ///// create histograms
    if (uniformIsoBinning) {

      dataHisto.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
				    new TH1F(Form("dataHisto_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
					     nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
      signalTemplateRND04_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("signalTemplateRND04_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							    nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
      signalTemplateRND08_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("signalTemplateRND08_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							    nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
      backgroundTemplate_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						  new TH1F(Form("backgroundTemplate_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							   nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  

    } else {

      cout << "pt " << ptBins.at(ibin) << "   nIsoBins " << getNbinsPhIsoVarWidth(ptBins.at(ibin)) << endl;
      cout << "iso binning --> ";
      vector<float> vecptr = getPhIsoBinVector(ptBins.at(ibin));
      for (int i = 0; i <= getNbinsPhIsoVarWidth(ptBins.at(ibin)); i++) {
	cout << vecptr.at(i) << " ";
      }
      cout << endl;      

      dataHisto.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
				    new TH1F(Form("dataHisto_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
					     getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));  
      signalTemplateRND04_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("signalTemplateRND04_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							    getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));  
      signalTemplateRND08_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("signalTemplateRND08_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							    getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));  
      backgroundTemplate_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						  new TH1F(Form("backgroundTemplate_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							   getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));  

    } 
    /////
    dataHisto.back().phHisto->Sumw2();
    signalTemplateRND04_data.back().phHisto->Sumw2();
    signalTemplateRND08_data.back().phHisto->Sumw2();
    backgroundTemplate_data.back().phHisto->Sumw2();

    
    if(addSystematics){ // create alternative templates for signal and background

      if (uniformIsoBinning) {

	signalTemplate_gjets.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						 new TH1F(Form("signalTemplate_gjets_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));
	signalTemplateRND04_gjets.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						      new TH1F(Form("signalTemplateRND04_gjets_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));
	backgroundTemplate_qcd.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("backgroundTemplate_qcd_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));

      } else {

	signalTemplate_gjets.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						 new TH1F(Form("signalTemplate_gjets_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							  getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));
	signalTemplateRND04_gjets.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						      new TH1F(Form("signalTemplateRND04_gjets_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							       getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));
	backgroundTemplate_qcd.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						   new TH1F(Form("backgroundTemplate_qcd_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",
							    getNbinsPhIsoVarWidth(ptBins.at(ibin)),getPhIsoBinArray(ptBins.at(ibin)))));

      }

  
      //////
      signalTemplate_gjets.back().phHisto->Sumw2();
      signalTemplateRND04_gjets.back().phHisto->Sumw2();
      backgroundTemplate_qcd.back().phHisto->Sumw2();
    }

  }



  // add specific files to the chain --> data
  TChain* chain_data = new TChain("tree/tree");
  system(("ls "+inputDirectory+" | grep root > file.list").c_str());
  ifstream infile("file.list");
  string line;
  if(infile.is_open()){
    while(getline(infile,line)){
      if(line == "") continue;
      if(not TString(line).Contains("root")) continue;
      bool found = false;
      for(auto era : RunEra){
	if(TString(line).Contains(era.c_str()))
	  found = true;
      }
      if(found){
	cout<<"Add data file in chain "<<inputDirectory+"/"+line<<endl;
	chain_data->Add((inputDirectory+"/"+line).c_str());
      }
    }
  }
  system("rm file.list");

  
  //files for MC background and signal
  TChain* chain_gjets = new TChain("tree/tree");
  TChain* chain_qcd   = new TChain("tree/tree");
  TChain* genchain_gjets = new TChain("gentree/gentree");
  TChain* genchain_qcd   = new TChain("gentree/gentree");
  TFile*  kfactorFile = NULL;
  vector<TH1*> khists;

  if(addSystematics){

    cout<<"Add gamma+jets file in chain "<<inputDirectorySignalMC+"/*root"<<endl;
    chain_gjets->Add((inputDirectorySignalMC+"/*root").c_str());
    genchain_gjets->Add((inputDirectorySignalMC+"/*root").c_str());

    cout<<"Add qcd multijets file in chain "<<inputDirectoryBackgroundMC+"/*root"<<endl;
    chain_qcd->Add((inputDirectoryBackgroundMC+"/*root").c_str());    
    genchain_qcd->Add((inputDirectoryBackgroundMC+"/*root").c_str());    

    // k-factor for gjets
    kfactorFile = TFile::Open(kfactorFileName.c_str());
    TH1*  alohist   = (TH1*) kfactorFile->Get("GJets_LO/inv_pt_G");
    TH1*  anlohist  = (TH1*) kfactorFile->Get("GJets_1j_NLO/nominal_G");
    TH1*  aewkhist  = (TH1*) kfactorFile->Get("EWKcorr/photon");
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);

    khists.push_back(aewkhist);
    khists.push_back(anlohist);    
  }
  
  // output file
  TFile* outputFile = new TFile((outputDIR+"/PhotonPurityFitResult.root").c_str(),"RECREATE");
  outputFile->cd();

  // fillHistograms for data  
  fillDataHistograms(chain_data,Sample::data,dataHisto,signalTemplateRND04_data,signalTemplateRND08_data,backgroundTemplate_data,mediumID);  

  // to calculate mean pt
  float mean = 0;
  for(int ibin = 0; ibin < dataHisto.size(); ibin++){
    mean = dataHisto.at(ibin).ptMean/dataHisto.at(ibin).phHisto->Integral();
    dataHisto.at(ibin).ptMean = mean;
  }
  for(int ibin = 0; ibin < signalTemplateRND04_data.size(); ibin++){
    mean = signalTemplateRND04_data.at(ibin).ptMean/signalTemplateRND04_data.at(ibin).phHisto->Integral();
    signalTemplateRND04_data.at(ibin).ptMean = mean;
  }
  for(int ibin = 0; ibin < signalTemplateRND08_data.size(); ibin++){
    mean = signalTemplateRND08_data.at(ibin).ptMean/signalTemplateRND08_data.at(ibin).phHisto->Integral();
    signalTemplateRND08_data.at(ibin).ptMean = mean;
  }
  for(int ibin = 0; ibin < backgroundTemplate_data.size(); ibin++){
    mean = backgroundTemplate_data.at(ibin).ptMean/backgroundTemplate_data.at(ibin).phHisto->Integral();
    backgroundTemplate_data.at(ibin).ptMean = mean;
  }

  if (debug) {  // print binning only once

    cout << "===== makePhotonPurityFit.C: data iso binning         --> ";
    for (int i = 1; i < (dataHisto.at(0).phHisto->GetNbinsX()+2); i++) {
      cout << dataHisto.at(0).phHisto->GetXaxis()->GetBinLowEdge(i) << " ";
    }
    cout << endl;
    cout << "===== makePhotonPurityFit.C: sig template iso binning --> ";
    for (int i = 1; i < (signalTemplateRND04_data.at(0).phHisto->GetNbinsX()+2); i++) {
      cout << signalTemplateRND04_data.at(0).phHisto->GetXaxis()->GetBinLowEdge(i) << " ";
    }
    cout << endl;
    cout << "===== makePhotonPurityFit.C: bkg template iso binning --> ";
    for (int i = 1; i < (backgroundTemplate_data.at(0).phHisto->GetNbinsX()+2); i++) {
      cout << backgroundTemplate_data.at(0).phHisto->GetXaxis()->GetBinLowEdge(i) << " ";
    }
    cout << endl;

  }


  // Build Model for fit
  vector<RooWorkspace*> workspaceRND04;
  vector<RooWorkspace*> workspaceRND08;
  vector<RooWorkspace*> workspaceRND04_altSig;
  vector<RooWorkspace*> workspaceRND04_altBkg;

  TGraphAsymmErrors* photonPurityRND04 = new TGraphAsymmErrors();
  TGraphAsymmErrors* photonPurityRND08 = new TGraphAsymmErrors();

  TGraphAsymmErrors* photonPurityRND04_altSig = new TGraphAsymmErrors();
  TGraphAsymmErrors* photonPurityRND04_altBkg = new TGraphAsymmErrors();

  TCanvas* canvas = new TCanvas("canvas","canvas",600,700);
  canvas->cd();
  
  cout<<"##### Start fitting data ... "<<endl;
  outputFile->cd();
  // fit of real data with data-driven templates
  for(size_t isize = 0; isize < dataHisto.size(); isize++){

    int ptMin = int(dataHisto.at(isize).ptMin);
    int ptMax = int(dataHisto.at(isize).ptMax); 

    cout<<"####### Fit bin RND04 : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
    //create workspace and fit for RND = 0.4
    workspaceRND04.push_back(new RooWorkspace(Form("wsRND04_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_pt_%d_%d",ptMin,ptMax)));
    makePurityFit(workspaceRND04.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
    RooRealVar* purity = workspaceRND04.back()->var("photonPurity");
    photonPurityRND04->SetPoint(isize,double((ptMax+ptMin)/2.),purity->getVal());   
    photonPurityRND04->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());
    // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
    //   photonPurityRND04->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
    // save in the output file
    plotFitResult(canvas,dataHisto.at(isize).phHisto,workspaceRND04.back(),outputDIR,ptMin,ptMax,"RND04",lumi,uniformIsoBinning);
    
    //create workspace and fit for RND = 0.8
    cout<<"####### Fit bin RND08 : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
    workspaceRND08.push_back(new RooWorkspace(Form("wsRND08_data_pt_%d_%d",ptMin,ptMax),Form("wsRND08_pt_%d_%d",ptMin,ptMax)));    
    makePurityFit(workspaceRND08.back(),dataHisto.at(isize),signalTemplateRND08_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
    purity = workspaceRND08.back()->var("photonPurity");
    photonPurityRND08->SetPoint(isize,double(ptMax+ptMin)/2.,purity->getVal());
    photonPurityRND08->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());
    // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
    //   photonPurityRND08->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
    plotFitResult(canvas,dataHisto.at(isize).phHisto,workspaceRND08.back(),outputDIR,ptMin,ptMax,"RND08",lumi,uniformIsoBinning);

    if(addSystematics and makeFitBasedOnlyOnTemplates == false){

      cout<<"####### Fit bin alt sig : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
      workspaceRND04_altSig.push_back(new RooWorkspace(Form("wsRND04_altSig_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_altSig_pt_%d_%d",ptMin,ptMax)));    
      makePurityFit(workspaceRND04_altSig.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates,true,false);
      RooRealVar* purity = workspaceRND04_altSig.back()->var("photonPurity");
      photonPurityRND04_altSig->SetPoint(isize,double((ptMax+ptMin)/2),purity->getVal());
      photonPurityRND04_altSig->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());
      // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
      // 	photonPurityRND04_altSig->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
      plotFitResult(canvas,dataHisto.at(isize).phHisto,workspaceRND04_altSig.back(),outputDIR,ptMin,ptMax,"RND04_altSig",lumi,uniformIsoBinning);

      cout<<"####### Fit bin alt bkg : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
      workspaceRND04_altBkg.push_back(new RooWorkspace(Form("wsRND04_altBkg_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_altBkg_pt_%d_%d",ptMin,ptMax)));    
      makePurityFit(workspaceRND04_altBkg.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates,false,true);
      purity = workspaceRND04_altBkg.back()->var("photonPurity");
      photonPurityRND04_altBkg->SetPoint(isize,double((ptMax+ptMin)/2),purity->getVal());
      photonPurityRND04_altBkg->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());
      // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
      // 	photonPurityRND04_altBkg->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
      plotFitResult(canvas,dataHisto.at(isize).phHisto,workspaceRND04_altBkg.back(),outputDIR,ptMin,ptMax,"RND04_altBkg",lumi,uniformIsoBinning);
    }
  }

  // save workspaces
  for(auto ws: workspaceRND04)
    ws->Write();
  for(auto ws: workspaceRND08)
    ws->Write();
  for(auto ws: workspaceRND04_altSig)
    ws->Write();
  for(auto ws: workspaceRND04_altBkg)
    ws->Write();

  // Build Model for fit
  vector<RooWorkspace*> workspace_gjets; // alternative signal
  vector<RooWorkspace*> workspace_qcd;   // alternative background

  TGraphAsymmErrors* photonPurity_gjets = new TGraphAsymmErrors();
  TGraphAsymmErrors* photonPurity_qcd = new TGraphAsymmErrors();

  if(addSystematics){ // implement sys uncertainties using MC based templates to fit data

    // fillHistograms for gamma+jets --> temp values                                                                                                          
                                   
    fillMCHistograms(chain_gjets,Sample::gjets,signalTemplate_gjets,mediumID,khists,lumi,genchain_gjets);                                   
    fillMCHistograms(chain_gjets,Sample::gjets,signalTemplateRND04_gjets,mediumID,khists,lumi,genchain_gjets,true);
    // to calculate mean pt                                                                                                                                                 

    for(size_t ibin = 0; ibin < signalTemplate_gjets.size(); ibin++) {
      signalTemplate_gjets.at(ibin).ptMean = signalTemplate_gjets.at(ibin).ptMean/signalTemplate_gjets.at(ibin).phHisto->Integral();
    }
    // to calculate mean pt                                                                                                                                                
                      
    for(size_t ibin = 0; ibin < signalTemplateRND04_gjets.size(); ibin++) {
      signalTemplateRND04_gjets.at(ibin).ptMean = signalTemplateRND04_gjets.at(ibin).ptMean/signalTemplateRND04_gjets.at(ibin).phHisto->Integral();
    }
    // fillHistograms for qcd                                                                                                                        
                                  
    fillMCHistograms(chain_qcd,Sample::qcd,backgroundTemplate_qcd,mediumID,khists,lumi,genchain_qcd);
    // to calculate mean pt                                                                                                                                  
                                     
    for(size_t ibin = 0; ibin <  backgroundTemplate_qcd.size(); ibin++) {
      backgroundTemplate_qcd.at(ibin).ptMean = backgroundTemplate_qcd.at(ibin).ptMean/backgroundTemplate_qcd.at(ibin).phHisto->Integral();
    }
    // make alternative fits
    outputFile->cd();
    for(size_t isize = 0; isize < dataHisto.size(); isize++){
      
      int ptMin = int(dataHisto.at(isize).ptMin);
      int ptMax = int(dataHisto.at(isize).ptMax);

      cout<<"####### Fit bin Gamma+jets : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;      
      //create workspace and fit (this is not for RND = 0.4) using gamma+jets MC                                                                                          
                                     
      workspace_gjets.push_back(new RooWorkspace(Form("ws_gjets_pt_%d_%d",ptMin,ptMax),Form("ws_gjets_%d_%d",ptMin,ptMax)));
      makePurityFit(workspace_gjets.back(),dataHisto.at(isize),signalTemplate_gjets.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
      RooRealVar* purity = workspace_gjets.back()->var("photonPurity");
      photonPurity_gjets->SetPoint(isize,double((ptMax+ptMin)/2),purity->getVal());
      photonPurity_gjets->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());      
      // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
      // 	photonPurity_gjets->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
      // plot fit result:                                                                                                                                                   
      plotFitResult(canvas,dataHisto.at(isize).phHisto,workspace_gjets.back(),outputDIR,ptMin,ptMax,"gjets",lumi,uniformIsoBinning);

      //create workspace and fit for QCD MC /sigma_ietaieta sideband)                                                                                                       
                        
      workspace_qcd.push_back(new RooWorkspace(Form("ws_qcd_pt_%d_%d",ptMin,ptMax),Form("ws_qcd_%d_%d",ptMin,ptMax)));
      makePurityFit(workspace_qcd.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_qcd.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
      purity = workspace_qcd.back()->var("photonPurity");
      photonPurity_qcd->SetPoint(isize,double((ptMax+ptMin)/2),purity->getVal());
      photonPurity_qcd->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(purity->getErrorLo()),purity->getErrorHi());      
      // if(purity->getErrorLo() > 0.01 or purity->getErrorHi() > 0.01)
      // 	photonPurity_qcd->SetPointError(isize,double(ptMax-ptMin)/2.,double(ptMax-ptMin)/2.,fabs(0.01),0.01);
      // plot fit result:                                                                                                                                                   
      plotFitResult(canvas,dataHisto.at(isize).phHisto,workspace_qcd.back(),outputDIR,ptMin,ptMax,"qcd",lumi,uniformIsoBinning);
    }
  }

  for(auto ws : workspace_gjets)
    ws->Write();
  for(auto ws : workspace_qcd)
    ws->Write();

    
  // plot purity result 
  TCanvas* canvas2 = new TCanvas("canvas2","",600,650);
  canvas2->cd();
  canvas2->SetTickx(1);
  canvas2->SetTicky(1);
  canvas2->cd();
  canvas2->SetRightMargin(0.06);

  photonPurityRND04->SetLineColor(kBlack);
  photonPurityRND04->SetMarkerColor(kBlack);
  photonPurityRND04->SetMarkerStyle(20);
  photonPurityRND04->SetMarkerSize(1);
  photonPurityRND04->GetXaxis()->SetTitle("photon p_{T} [GeV]");
  photonPurityRND04->GetYaxis()->SetTitle("photon purity");
  photonPurityRND04->GetYaxis()->SetRangeUser(0.5,1.2);
  photonPurityRND04->Draw("AP");

  photonPurityRND08->SetLineColor(kRed);
  photonPurityRND08->SetMarkerColor(kRed);
  photonPurityRND08->SetMarkerStyle(24);
  photonPurityRND08->SetMarkerSize(1);
  photonPurityRND08->Draw("Psame");

  if(addSystematics){
    photonPurity_gjets->SetLineColor(kBlue);
    photonPurity_gjets->SetMarkerColor(kBlue);
    photonPurity_gjets->SetMarkerStyle(24);
    photonPurity_gjets->SetMarkerSize(1);
    photonPurity_gjets->Draw("Psame");

    photonPurity_qcd->SetLineColor(kCyan+1);
    photonPurity_qcd->SetMarkerColor(kCyan+1);
    photonPurity_qcd->SetMarkerStyle(24);
    photonPurity_qcd->SetMarkerSize(1);
    photonPurity_qcd->Draw("Psame");
  }

  CMS_lumi(canvas2,Form("%.1f",lumi));
  
  TLegend leg (0.3,0.3,0.6,0.6);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(photonPurityRND04,"Purity with #DeltaR = 0.4","PLE");
  leg.AddEntry(photonPurityRND08,"Purity with #DeltaR = 0.8","PLE");
  if(addSystematics){
    leg.AddEntry(photonPurity_gjets,"Purity from #gamma+jets","PLE");
    leg.AddEntry(photonPurity_qcd,"Purity from QCD","PLE");
  }
  leg.Draw("same");

  canvas2->SaveAs((outputDIR+"/photonPurityComparison.png").c_str(),"png");
  canvas2->SaveAs((outputDIR+"/photonPurityComparison.pdf").c_str(),"pdf");

  if(addSystematics and makeFitBasedOnlyOnTemplates == false){
    canvas2->cd();
    photonPurityRND04->Draw("AP");
    photonPurityRND04_altSig->SetLineColor(kBlue);
    photonPurityRND04_altSig->SetMarkerColor(kBlue);
    photonPurityRND04_altSig->SetMarkerStyle(24);
    photonPurityRND04_altSig->SetMarkerSize(1);
    photonPurityRND04_altSig->Draw("Psame");

    photonPurityRND04_altBkg->SetLineColor(kRed);
    photonPurityRND04_altBkg->SetMarkerColor(kRed);
    photonPurityRND04_altBkg->SetMarkerStyle(24);
    photonPurityRND04_altBkg->SetMarkerSize(1);
    photonPurityRND04_altBkg->Draw("Psame");
    
    CMS_lumi(canvas2,Form("%.1f",lumi));
    
    TLegend leg (0.3,0.3,0.6,0.6);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(photonPurityRND04,"Purity with #DeltaR = 0.4","PLE");
    leg.AddEntry(photonPurityRND04_altSig,"Purity with CB","PLE");
    leg.AddEntry(photonPurityRND04_altBkg,"Purity with Power-Law","PLE");
    leg.Draw("same");

    canvas2->SaveAs((outputDIR+"/photonPurityComparison_2.png").c_str(),"png");
    canvas2->SaveAs((outputDIR+"/photonPurityComparison_2.pdf").c_str(),"pdf");
  }
  
  //make final plot with statistical and sys errors
  TGraphAsymmErrors* finalPurity = new TGraphAsymmErrors();

  if(addSystematics){

    cout<<"###### Purity Result with systematics "<<endl;
    for(int ipoint = 0; ipoint < photonPurityRND04->GetN(); ipoint++){
      
      // nominal value and uncertainty
      double x, y, err_x_low, err_x_high, err_y_low, err_y_high;
      photonPurityRND04->GetPoint(ipoint,x,y);
      err_x_low  = photonPurityRND04->GetErrorXlow(ipoint);
      err_x_high = photonPurityRND04->GetErrorXhigh(ipoint);
      err_y_low  = photonPurityRND04->GetErrorYlow(ipoint);
      err_y_high = photonPurityRND04->GetErrorYhigh(ipoint);
      finalPurity->SetPoint(ipoint,x,y);

      double err_y_low_stat = err_y_low;
      double err_y_high_stat = err_y_high;

      //sys variations
      double x_rnd08, y_rnd08;
      photonPurityRND08->GetPoint(ipoint,x_rnd08,y_rnd08);

      double x_gjets, y_gjets;
      photonPurity_gjets->GetPoint(ipoint,x_gjets,y_gjets);
      
      double x_qcd, y_qcd;
      photonPurity_qcd->GetPoint(ipoint,x_qcd,y_qcd);
      

      if(makeFitBasedOnlyOnTemplates == true){
	double err_y_low_sys  = sqrt(fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y));	
	double err_y_high_sys = sqrt(fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y));	

	// inflate a bit if it is too low
	// if(err_y_low_sys < 0.015)  err_y_low_sys = err_y_low_sys+0.008;
	// if(err_y_high_sys < 0.015) err_y_high_sys = err_y_high_sys+0.008;

	finalPurity->SetPointError(ipoint,err_x_low,err_x_high,sqrt(err_y_low_stat*err_y_low_stat+err_y_low_sys*err_y_low_sys),sqrt(err_y_high_stat*err_y_high_stat+err_y_high_sys*err_y_high_sys));
	
	cout<<"#### Central value: x = "<<x<<" purity = "<<y<<" - "<<err_y_low_stat<<" + "<<err_y_high_stat<<" (stat) "<<" - "<<err_y_low_sys<<" + "<<err_y_high_sys<<endl;      

      }
      else{

	//sys variations
	double x_altSig, y_altSig;
	photonPurityRND04_altSig->GetPoint(ipoint,x_altSig,y_altSig);       
	
	//sys variations
	double x_altBkg, y_altBkg;
	photonPurityRND04_altBkg->GetPoint(ipoint,x_altBkg,y_altBkg);	

	double err_y_low_sys  = sqrt(fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y)+
				     fabs(y_altSig-y)*fabs(y_altSig-y)+fabs(y_altBkg-y)*fabs(y_altBkg-y));
	
	double err_y_high_sys = sqrt(fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y)+
				     fabs(y_altSig-y)*fabs(y_altSig-y)+fabs(y_altBkg-y)*fabs(y_altBkg-y));

	// inflate a bit if it is too low
	// if(err_y_low_sys < 0.015) err_y_low_sys = err_y_low_sys+0.008;
	// if(err_y_high_sys < 0.015) err_y_high_sys = err_y_high_sys+0.008;
	
	finalPurity->SetPointError(ipoint,err_x_low,err_x_high,sqrt(err_y_low_stat*err_y_low_stat+err_y_low_sys*err_y_low_sys),sqrt(err_y_high_stat*err_y_high_stat+err_y_high_sys*err_y_high_sys));

	cout<<"#### Central value: x = "<<x<<" purity = "<<y<<" - "<<err_y_low_stat<<" + "<<err_y_high_stat<<" (stat) "<<" - "<<err_y_low_sys<<" + "<<err_y_high_sys<<endl;      
	
      }
    }
    
    canvas2->cd();
    finalPurity->SetLineColor(kBlack);
    finalPurity->SetMarkerColor(kBlack);
    finalPurity->SetMarkerStyle(20);
    finalPurity->SetMarkerSize(1);
    finalPurity->GetYaxis()->SetRangeUser(0.5,1.2);
    finalPurity->SetFillColor(kOrange+1);
    photonPurityRND04->Draw("AP");
    finalPurity->Draw("E2same");
    photonPurityRND04->Draw("Psame");

    CMS_lumi(canvas2,Form("%.1f",lumi));

    canvas2->SaveAs((outputDIR+"/photonPurityFinal.png").c_str(),"png");
    canvas2->SaveAs((outputDIR+"/photonPurityFinal.pdf").c_str(),"pdf");    
    finalPurity->Write("purity");
  }
  
  if(saveHistograms){    
    for(size_t isize = 0; isize < dataHisto.size(); isize++){      
      
      if(not outputFile->GetDirectory("Templates"))
	outputFile->mkdir("Templates");
      outputFile->cd("Templates");

      dataHisto.at(isize).phHisto->Write();
      signalTemplateRND04_data.at(isize).phHisto->Write();
      signalTemplateRND08_data.at(isize).phHisto->Write();
      backgroundTemplate_data.at(isize).phHisto->Write();

      if(addSystematics){
	signalTemplate_gjets.at(isize).phHisto->Write();
	signalTemplateRND04_gjets.at(isize).phHisto->Write();
	backgroundTemplate_qcd.at(isize).phHisto->Write();
      }      
    }
    outputFile->cd();
  }
  outputFile->Close();
}
