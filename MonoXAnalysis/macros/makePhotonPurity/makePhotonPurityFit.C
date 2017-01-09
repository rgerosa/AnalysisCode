#include "makePhotonPurityUtils.h"
#include <fstream>

// to decide which data to run
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};
// photon pt bins
static vector<float> ptBins = {175,200,225,250,280,320,375,425,1000};
// photon isolation info
static vector<int>   nBinPhotonIso = {30,30,30,30,30,25,25,25,25};
static vector<float> photonIsoMax  = {20,20,20,20,20,20,20,20,20};
static vector<float> photonIsoMin  = {0,0,0,0,0,0,0,0,0,0};
// debug mode
static bool debug = false;
static bool saveHistograms = true;
// k-facotr file
static string kfactorFileName = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";


void makePhotonPurityFit(string inputDirectory, // directory with dataFiles
			 float  lumi, // luminosity
			 string outputDIR,
			 bool   addSystematics = false,
			 string inputDirectorySignalMC = "",
			 string inputDirectoryBackgroundMC = "",
			 bool   makeFitBasedOnlyOnTemplates = false
			 ){

  system(("mkdir -p "+outputDIR).c_str());

  // style
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  //from twiki https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
  photonID mediumID (0.0396,0.01022,0.01400,0.441,2.725,0.0148,0.000017,2.571,0.0047); // set wp for medium id

  // bins for purity and histograms
  vector<fitPurity> dataHisto;
  vector<fitPurity> signalTemplateRND04_data;
  vector<fitPurity> signalTemplateRND08_data;
  vector<fitPurity> backgroundTemplate_data;
  vector<fitPurity> signalTemplate_gjets;
  vector<fitPurity> backgroundTemplate_qcd;

  for(size_t ibin = 0; ibin < ptBins.size()-1; ibin++){

    ///// create histograms
    dataHisto.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
				  new TH1F(Form("dataHisto_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
    signalTemplateRND04_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						 new TH1F(Form("signalTemplateRND04_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
    signalTemplateRND08_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						 new TH1F(Form("signalTemplateRND08_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
    backgroundTemplate_data.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						new TH1F(Form("backgroundTemplate_data_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
    /////
    dataHisto.back().phHisto->Sumw2();
    signalTemplateRND04_data.back().phHisto->Sumw2();
    signalTemplateRND08_data.back().phHisto->Sumw2();
    backgroundTemplate_data.back().phHisto->Sumw2();

    if(addSystematics){ // create alternative templates for signal and background

      signalTemplate_gjets.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
					       new TH1F(Form("signalTemplate_gjets_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));
      backgroundTemplate_qcd.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
						 new TH1F(Form("backgroundTemplate_qcd_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso.at(ibin),photonIsoMin.at(ibin),photonIsoMax.at(ibin))));  
      //////
      signalTemplate_gjets.back().phHisto->Sumw2();
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

    // k-factor gor gjets
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

  // Build Model for fit
  vector<RooWorkspace*> worksapceRND04;
  vector<RooWorkspace*> worksapceRND08;
  vector<RooWorkspace*> worksapceRND04_altSig;
  vector<RooWorkspace*> worksapceRND04_altBkg;

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
    worksapceRND04.push_back(new RooWorkspace(Form("wsRND04_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_pt_%d_%d",ptMin,ptMax)));

    if(saveHistograms){
      dataHisto.at(isize).phHisto->Write();
      signalTemplateRND04_data.at(isize).phHisto->Write();
      signalTemplateRND08_data.at(isize).phHisto->Write();
      backgroundTemplate_data.at(isize).phHisto->Write();
    }

    makePurityFit(worksapceRND04.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
    RooRealVar* purity = worksapceRND04.back()->var("photonPurity");
    photonPurityRND04->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
    photonPurityRND04->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());
    // save in the output file
    plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapceRND04.back(),outputDIR,ptMin,ptMax,"RND04",lumi);
    
    //create workspace and fit for RND = 0.8
    cout<<"####### Fit bin RND08 : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
    worksapceRND08.push_back(new RooWorkspace(Form("wsRND08_data_pt_%d_%d",ptMin,ptMax),Form("wsRND08_pt_%d_%d",ptMin,ptMax)));    
    makePurityFit(worksapceRND08.back(),dataHisto.at(isize),signalTemplateRND08_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates);
    purity = worksapceRND08.back()->var("photonPurity");
    photonPurityRND08->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
    photonPurityRND08->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());
    plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapceRND08.back(),outputDIR,ptMin,ptMax,"RND08",lumi);

    if(addSystematics and makeFitBasedOnlyOnTemplates == false){
      cout<<"####### Fit bin alt sig : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
      worksapceRND04_altSig.push_back(new RooWorkspace(Form("wsRND04_altSig_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_altSig_pt_%d_%d",ptMin,ptMax)));    
      makePurityFit(worksapceRND04_altSig.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates,true,false);
      purity = worksapceRND04_altSig.back()->var("photonPurity");
      photonPurityRND04_altSig->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
      photonPurityRND04_altSig->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());
      plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapceRND04_altSig.back(),outputDIR,ptMin,ptMax,"RND04_altSig",lumi);
      
      cout<<"####### Fit bin alt bkg : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
      worksapceRND04_altBkg.push_back(new RooWorkspace(Form("wsRND04_altBkg_data_pt_%d_%d",ptMin,ptMax),Form("wsRND04_altBkg_pt_%d_%d",ptMin,ptMax)));    
      makePurityFit(worksapceRND04_altBkg.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,makeFitBasedOnlyOnTemplates,false,true);
      purity = worksapceRND04_altBkg.back()->var("photonPurity");
      photonPurityRND04_altBkg->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
      photonPurityRND04_altBkg->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());
      plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapceRND04_altBkg.back(),outputDIR,ptMin,ptMax,"RND04_altBkg",lumi);
    }
  }

  // save workspaces
  for(auto ws: worksapceRND04)
    ws->Write();
  for(auto ws: worksapceRND08)
    ws->Write();
  for(auto ws: worksapceRND04_altSig)
    ws->Write();
  for(auto ws: worksapceRND04_altBkg)
    ws->Write();

  // Build Model for fit
  vector<RooWorkspace*> worksapce_gjets; // alternative signal
  vector<RooWorkspace*> worksapce_qcd;   // alternative background

  TGraphAsymmErrors* photonPurity_gjets = new TGraphAsymmErrors();
  TGraphAsymmErrors* photonPurity_qcd = new TGraphAsymmErrors();

  if(addSystematics){ // implement sys uncertainties using MC based templates to fit data

    // fillHistograms for gamma+jets --> temp values                                                                                                                                             
    fillMCHistograms(chain_gjets,Sample::gjets,signalTemplate_gjets,mediumID,khists,lumi,genchain_gjets);
    // to calculate mean pt                                                                                                                                                                       
    for(size_t ibin = 0; ibin < signalTemplate_gjets.size(); ibin++)
      signalTemplate_gjets.at(ibin).ptMean = signalTemplate_gjets.at(ibin).ptMean/signalTemplate_gjets.at(ibin).phHisto->Integral();
    
    // fillHistograms for qcd                                                                                                                                                          
    fillMCHistograms(chain_qcd,Sample::qcd,backgroundTemplate_qcd,mediumID,khists,lumi,genchain_qcd);
    // to calculate mean pt                                                                                                                                                                       
    for(size_t ibin = 0; ibin <  backgroundTemplate_qcd.size(); ibin++)
      backgroundTemplate_qcd.at(ibin).ptMean = backgroundTemplate_qcd.at(ibin).ptMean/backgroundTemplate_qcd.at(ibin).phHisto->Integral();

    // make alternative fits
    outputFile->cd();
    for(size_t isize = 0; isize < dataHisto.size(); isize++){
      
      int ptMin = int(dataHisto.at(isize).ptMin);
      int ptMax = int(dataHisto.at(isize).ptMax);

      if(saveHistograms){
	signalTemplate_gjets.at(isize).phHisto->Write();
	backgroundTemplate_qcd.at(isize).phHisto->Write();
      }

      cout<<"####### Fit bin Gamma+jets : ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;      
      //create workspace and fit for RND = 0.4 using gamma+jets MC                                                                                                                               
      worksapce_gjets.push_back(new RooWorkspace(Form("ws_gjets_pt_%d_%d",ptMin,ptMax),Form("ws_gjets_%d_%d",ptMin,ptMax)));

      makePurityFit(worksapce_gjets.back(),dataHisto.at(isize),signalTemplate_gjets.at(isize),backgroundTemplate_data.at(isize),mediumID,debug,true);
      RooRealVar* purity = worksapce_gjets.back()->var("photonPurity");
      photonPurity_gjets->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
      photonPurity_gjets->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());      
      // plot fit result:                                                                                                                                                                         
      plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapce_gjets.back(),outputDIR,ptMin,ptMax,"gjets",lumi);

      //create workspace and fit for RND = 0.4 using gamma+jets MC                                                                                                                               
      worksapce_qcd.push_back(new RooWorkspace(Form("ws_qcd_pt_%d_%d",ptMin,ptMax),Form("ws_qcd_%d_%d",ptMin,ptMax)));
      makePurityFit(worksapce_qcd.back(),dataHisto.at(isize),signalTemplateRND04_data.at(isize),backgroundTemplate_qcd.at(isize),mediumID,debug,true);
      purity = worksapce_qcd.back()->var("photonPurity");
      photonPurity_qcd->SetPoint(isize,(ptMax+ptMin)/2,purity->getVal());
      photonPurity_qcd->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,fabs(purity->getErrorLo()),purity->getErrorHi());      
      // plot fit result:                                                                                                                                                                         
      plotFitResult(canvas,dataHisto.at(isize).phHisto,worksapce_qcd.back(),outputDIR,ptMin,ptMax,"qcd",lumi);
    }
  }
  
  for(auto ws : worksapce_gjets)
    ws->Write();
  for(auto ws : worksapce_qcd)
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

    photonPurityRND04_altBkg->SetLineColor(kCyan+1);
    photonPurityRND04_altBkg->SetMarkerColor(kCyan+1);
    photonPurityRND04_altBkg->SetMarkerStyle(24);
    photonPurityRND04_altBkg->SetMarkerSize(1);
    photonPurityRND04_altBkg->Draw("Psame");
    
    CMS_lumi(canvas2,Form("%.1f",lumi));
    
    TLegend leg (0.3,0.3,0.6,0.6);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(photonPurityRND04,"Purity with #DeltaR = 0.4","PLE");
    leg.AddEntry(photonPurityRND04_altSig,"Purity with Landau","PLE");
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
      cout<<"#### Central value: x = "<<x<<" purity = "<<y<<" err_low "<<err_y_low<<" err_high "<<err_y_high<<endl;
      
      finalPurity->SetPoint(ipoint,x,y);

      //sys variations
      double x_rnd08, y_rnd08;
      photonPurityRND08->GetPoint(ipoint,x_rnd08,y_rnd08);
      cout<<"#### Variation 1 with RND 0.8: (purity - purity_rnd08)/purity "<<fabs(y-y_rnd08)/y<<endl;

      double x_gjets, y_gjets;
      photonPurity_gjets->GetPoint(ipoint,x_gjets,y_gjets);
      cout<<"#### Variation 2 with gamma+jets: (purity - purity_gjets)/purity "<<fabs(y-y_gjets)/y<<endl;
      
      double x_qcd, y_qcd;
      photonPurity_qcd->GetPoint(ipoint,x_qcd,y_qcd);
      cout<<"#### Variation 3 with QCD EM enriched: (purity - purity_qcd)/purity "<<fabs(y-y_qcd)/y<<endl;      
      

      if(makeFitBasedOnlyOnTemplates == true){
	err_y_low  = sqrt(err_y_low*err_y_low+fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y));
	
	err_y_high = sqrt(err_y_high*err_y_high+fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y));	

	finalPurity->SetPointError(ipoint,err_x_low,err_x_high,err_y_low,err_y_high);
	cout<<"#### Total uncertainty : low "<<err_y_low<<" high "<<err_y_high<<endl;      
      }
      else{

	//sys variations
	double x_altSig, y_altSig;
	photonPurityRND04_altSig->GetPoint(ipoint,x_altSig,y_altSig);       
	cout<<"#### Variation 4 with Landau: (purity - purity_altSig)/purity "<<fabs(y-y_altSig)/y<<endl;
	
	//sys variations
	double x_altBkg, y_altBkg;
	photonPurityRND04_altBkg->GetPoint(ipoint,x_altBkg,y_altBkg);	
	cout<<"#### Variation 5 with Power Law: (purity - purity_altBkg)/purity "<<fabs(y-y_altBkg)/y<<endl;


	err_y_low  = sqrt(err_y_low*err_y_low+fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y)+
			  fabs(y_altSig-y)*fabs(y_altSig-y)+fabs(y_altBkg-y)*fabs(y_altBkg-y));
	
	err_y_high = sqrt(err_y_high*err_y_high+fabs(y_rnd08-y)*fabs(y_rnd08-y)+fabs(y_gjets-y)*fabs(y_gjets-y)+fabs(y_qcd-y)*fabs(y_qcd-y)+
			  fabs(y_altSig-y)*fabs(y_altSig-y)+fabs(y_altBkg-y)*fabs(y_altBkg-y));
	
	finalPurity->SetPointError(ipoint,err_x_low,err_x_high,err_y_low,err_y_high);
	cout<<"#### Total uncertainty : low "<<err_y_low<<" high "<<err_y_high<<endl;      
	
      }
    }
    
    canvas2->cd();
    finalPurity->SetLineColor(kBlack);
    finalPurity->SetMarkerColor(kBlack);
    finalPurity->SetMarkerStyle(20);
    finalPurity->SetMarkerSize(1);
    finalPurity->GetYaxis()->SetRangeUser(0.5,1.2);
    finalPurity->SetFillColor(kRed);
    photonPurityRND04->Draw("AP");
    finalPurity->Draw("E2same");
    photonPurityRND04->Draw("Psame");

    CMS_lumi(canvas2,Form("%.1f",lumi));

    canvas2->SaveAs((outputDIR+"/photonPurityFinal.png").c_str(),"png");
    canvas2->SaveAs((outputDIR+"/photonPurityFinal.pdf").c_str(),"pdf");
    
    finalPurity->Write("purity");
    outputFile->Close();
  }
}
