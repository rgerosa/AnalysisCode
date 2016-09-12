///////// Basic options
static bool  debug = false; // for debugging printout
static float muMin = 0.001; // min mu for limit scam
static float muMax = 0.001; // max mu for limit scan
static float muPoint = 300; // number of points in the limit scan
static bool  calculateExpSigma = true; // if calculate also 1s and 2s band on limit values
static bool  doLikelihoodScan  = true; // to perform a likelihood scan around minimum [-muMax,muMax] in muPoint
  
// multipy by bin width
void multipleByBinWidth (TH1F* hist){

  for(int iBin = 0; iBin < hist->GetNbinsX(); iBin++){
    hist->SetBinContent(iBin+1,hist->GetBinContent(iBin+1)*hist->GetBinWidth(iBin+1));
    hist->SetBinError(iBin+1,hist->GetBinError(iBin+1)*hist->GetBinWidth(iBin+1));
  } 
}

// multiply by bin width  a covariance from ML fit in combine
void multipleByBinWidth (TH2F* hist, TH1F* histo1){

  for(int iBin = 0; iBin < hist->GetNbinsX(); iBin++){
    for(int jBin = 0; jBin < hist->GetNbinsY(); jBin++){
      hist->SetBinContent(iBin+1,jBin+1,hist->GetBinContent(iBin+1,jBin+1)*histo1->GetXaxis()->GetBinWidth(iBin+1)*histo1->GetXaxis()->GetBinWidth(jBin+1));
      hist->SetBinError(iBin+1,jBin+1,hist->GetBinError(iBin+1,jBin+1)*histo1->GetXaxis()->GetBinWidth(iBin+1)*histo1->GetXaxis()->GetBinWidth(jBin+1));
    } 
  }
}

// reset a list of parameters to the original values
void resetRooArgList(RooArgList & list1, const RooArgList & list2){

  for(int ibin = 0; ibin < list1.getSize(); ibin++){
    float originalVal = ((RooRealVar*) list1.at(ibin))->getVal();
    RooRealVar* var1 = (RooRealVar*) list1.at(ibin);
    RooRealVar* var2 = (RooRealVar*) list2.at(ibin);
    var1->setVal(var2->getVal());
    var1->setRange(var2->getMin(),var2->getMax());
    if(debug)
      cout<<"resetRooArgList --> bin "<<ibin<<" original "<<originalVal<<" reset to "<<((RooRealVar*) list1.at(ibin))->getVal()<<endl;
  }
  
}

// calculate CLS using asymptotic formula for both observed and expected case, giving a quantile. Store also -Log(L) from S+B best fit
float calculateCLsLimit(RooAbsReal *nllSB, RooAbsReal *nllAsimov, RooRealVar *mu, double muVal, bool isObserved, double & nllData, double & nllA, const double & quantile = 0){

  mu->setConstant(false); 

  if(debug){
    cout<<"#################################"<<endl;
    cout<<"#################################"<<endl;
    if(isObserved)
      cout<<"##### Calculate Observed CLs ####"<<endl;
    else
      cout<<"##### Calculate Expected CLs ####"<<endl;
    cout<<"#################################"<<endl;
    cout<<"#################################"<<endl;
  }

  //make S+B fit of asimov toy --> used as background only hypothesis
  RooMinimizer mAsimovSB(*nllAsimov);
  if(not debug){
    mAsimovSB.setVerbose(kFALSE);
    mAsimovSB.setPrintLevel(-1);
  }
  mAsimovSB.minimize("Minuit2","minimize"); 
  if(debug)
    cout<<"##### -Log(L) Asimov Toy S+B fit Minimize Step: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
  mAsimovSB.minimize("Minuit2","migrad"); 
  if(debug)
    cout<<"##### -Log(L) Asimov Toy S+B fit Migrad Step: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
  mAsimovSB.minimize("Minuit2","migradimproved"); 
  if(debug)
    cout<<"##### -Log(L) Asimov Toy S+B fit Migrad Improved: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
  mAsimovSB.minos(RooArgSet(*mu));
  double minNLLAsimov_SB = nllAsimov->getVal();
  double muAsimov_SB     = mu->getVal();
  double muAsimov_SB_errUp     = mu->getErrorHi();
  double muAsimov_SB_errDw     = mu->getErrorLo();
  nllA = muAsimov_SB;
  if(debug)
    cout<<"##### -Log(L) Asimov Toy S+B fit Minos Step: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
  
  //make S+B fit of asimov toy fixed mu
  mu->setConstant(true);  mu->setVal(muVal);
  RooMinimizer mAsimovSfix(*nllAsimov);
  if(not debug){
    mAsimovSfix.setVerbose(kFALSE);
    mAsimovSfix.setPrintLevel(-1);
  }
  mAsimovSfix.minimize("Minuit2","minimize"); 
  if(debug)
    cout<<"##### -Log(L) Asimov Toy Sfix+B fit Minimize Step: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<endl;
  mAsimovSfix.minimize("Minuit2","migrad"); 
  if(debug)
    cout<<"##### -Log(L) Asimov Toy Sfix+B fit Migrad Step: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<endl;
  mAsimovSfix.minimize("Minuit2","migradimproved"); 
  double minNLLAsimov_Sfix = nllAsimov->getVal();
  double muAsimov_Sfix     = mu->getVal();
  if(debug)
    cout<<"##### -Log(L) Asimov Toy Sfix+B fit Migrad Improved: "<<nllAsimov->getVal()<<" mu "<<mu->getVal()<<endl;

  
  //make S+B again on data
  mu->setConstant(false); 
  RooMinimizer mSB(*nllSB);
  if(not debug){
    mSB.setVerbose(kFALSE);
    mSB.setPrintLevel(-1);
  }
  double minNLL_SB = 0 ;
  double mu_SB     = 0 ;

  if(isObserved){
    mSB.minimize("Minuit2","minimize");
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Minimize Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
    mSB.minimize("Minuit2","migrad");
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Migrad Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
    mSB.minimize("Minuit2","migradimproved");
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Migrad improved Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
    mSB.minos(RooArgSet(*mu));
    minNLL_SB = nllSB->getVal();
    mu_SB = mu->getVal();
    nllData   = minNLL_SB;
    if(debug)
      cout<<"##### S+B fit on data : -Log(L) "<<nllSB->getVal()<<" mu "<<mu->getVal()<<" err up "<<mu->getErrorHi()<<" err dw "<<mu->getErrorLo()<<endl;
  }
  

  //make S+B again on data
  mu->setConstant(true); mu->setVal(muVal);
  RooMinimizer mSfix(*nllSB);
  if(not debug){
    mSfix.setVerbose(kFALSE);
    mSfix.setPrintLevel(-1);
  }
  double minNLL_Sfix = 0;
  if(isObserved){
    mSfix.minimize("Minuit2","minimize");
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Minimize Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
    mSfix.minimize("Minuit2","migrad");
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Migrad Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
    mSfix.minimize("Minuit2","migradimproved");
    minNLL_Sfix = nllSB->getVal();
    if(debug)
      cout<<"##### -Log(L) Data S+B fit Migrad improved Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
  }

  
  //Calculate LHC test statistics
  double qmu = 2*(minNLL_Sfix - minNLL_SB); 
  if (qmu < 0){
    cerr<<" Problem -->negative q_mu --> something bad going on in the fit "<<endl;
    qmu = 0;
  }
  if (muVal < mu_SB){
    qmu = 0;
    cerr<<" Mu reference is lower than the one extracted from the fit "<<endl;
  }
  
  if(debug and isObserved)
    cout<<"LHC test statistics: Data wrt to mu = "<<muVal<<" qmu "<<qmu<<" -Log(L) @ mu = 1 "<<minNLL_Sfix<<" -Log(L) @ mu = mu_best "<<minNLL_SB<<endl;

  if(mu_SB < 0){ // re-do the fit fixing the strenght to zero according to LHC test stat
    mu->setVal(0);
    mu->setConstant(kTRUE);
    if(debug)
      cout<<"### Constraint mu to be zero --> re-do fits "<<endl;
    mSB.minimize("Minuit2","minimizer");
    if(debug)
      cout<<"##### -Log(L) Data S = 0 +B fit Minimize Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
    mSB.minimize("Minuit2","migrad");
    if(debug)
      cout<<"##### -Log(L) Data S = 0 +B fit Migrad Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
    mSB.minimize("Minuit2","migrad improved");
    qmu = 2*(minNLL_Sfix - nllSB->getVal());
    if(debug)
      cout<<"##### -Log(L) Data S = 0 +B fit Migrad Improved Step: "<<nllSB->getVal()<<" mu "<<mu->getVal()<<endl;
    if(debug)
      cout<<"LHC test statistics: Data wrt to mu = "<<muVal<<" qmu "<<qmu<<" -Log(L) @ mu = 1 "<<minNLL_Sfix<<" -Log(L) @ mu = mu_best "<<nllSB->getVal()<<endl;
  }
  
  
  double qA  = 2*(minNLLAsimov_Sfix - minNLLAsimov_SB); 
  if (qA < 0){
    qA = 0; // shouldn't this always be 0?
    if(debug)
      cerr<<"Test statistics for Asimov fit is less than zero --> check what's going on "<<endl;
  }
  
  if(debug)
    cout<<"LHC test statistics: Asimov wrt to mu = "<<muVal<<" qA "<<qA<<" -Log(L) @ mu = 1 "<<minNLLAsimov_Sfix<<" -Log(L) @ mu best "<<minNLLAsimov_SB<<endl;
  mu->setConstant(false);

  if(isObserved){

    double CLsb = ROOT::Math::normal_cdf_c(TMath::Sqrt(qmu));
    double CLb  = ROOT::Math::normal_cdf(TMath::Sqrt(qA)-TMath::Sqrt(qmu));
    
    if (qmu > qA) {
      // In this region, things are tricky
      double mos = TMath::Sqrt(qA); // mu/sigma
      CLsb = ROOT::Math::normal_cdf_c( (qmu + qA)/(2*mos) );
      CLb  = ROOT::Math::normal_cdf_c( (qmu - qA)/(2*mos) );
    }
    
    double CLs  = (CLb == 0 ? 0 : CLsb/CLb);
    if(debug)
      cout<<"Observed CLsb "<<CLsb<<" CLb "<<CLb<<" CLs "<<CLs<<endl;
    return CLs;
  }
  else{

    double N    = ROOT::Math::normal_quantile(quantile,1.0);
    double CLb  = quantile;
    double CLsb = ROOT::Math::normal_cdf_c( TMath::Sqrt(qA) - N, 1.);
    double CLs  = (CLb != 0 ? CLsb/CLb : 0); 
    if(debug)
      cout<<"Expected  CLsb "<<CLsb<<" CLb "<<CLb<<" CLs "<<CLs<<endl;
    return CLs;
  }

  return -1;
  
}

/////
double makeFitFixedStrenght(RooAbsReal *nll, RooRealVar *mu, double muVal){

  mu->setConstant(true);
  mu->setVal(muVal);
  RooMinimizer minim (*nll);
  if(not debug){
    minim.setVerbose(kFALSE);
    minim.setPrintLevel(-1);
  }
  
  minim.minimize("Minuit2","minimize");
  if(debug)
    cout<<"##### -Log(L) at a fixed signal strenght "<<muVal<<" : "<<nll->getVal()<<endl;
  minim.minimize("Minuit2","migrad");
  if(debug)
    cout<<"##### -Log(L) at a fixed signal strenght "<<muVal<<" : "<<nll->getVal()<<endl;

  return nll->getVal();
  
}


// take the data from imput workspace
TH1F* importDataHistogram(TFile* dataFile, const int & category){
    
  RooWorkspace* dataWS = NULL;
  RooRealVar* met_obs  = NULL;

  ///////
  if(category == 1){
    dataWS = (RooWorkspace*) dataFile->Get("SR_MJ");
    if(dataWS == 0 or dataWS == NULL){
      dataWS = (RooWorkspace*) dataFile->Get("combinedws");
    }  
  }
  else if(category == 2){    
    if(dataWS == 0 or dataWS == NULL){
      dataWS = (RooWorkspace*) dataFile->Get("SR_MV");
    }
    if(dataWS == 0 or dataWS == NULL){
      dataWS = (RooWorkspace*) dataFile->Get("combinedws");
    }  
  }

  ///////
  if(category == 1){
    met_obs = (RooRealVar*) dataWS->obj("met_MJ");
    if(met_obs == 0 or dataWS == NULL)
      met_obs = (RooRealVar*) dataWS->obj("met_monojet");
  }
  else if(category == 2){
      met_obs = (RooRealVar*) dataWS->obj("met_MV");
      if(met_obs == 0 or dataWS == NULL)
	met_obs = (RooRealVar*) dataWS->obj("met_monov");
  }
    
  ///////
  RooDataHist* data = NULL;
  if(category == 1){
    data = (RooDataHist*) dataWS->obj("data_obs_SR_MJ");
    if(data == 0 or data == NULL)
      data = (RooDataHist*) dataWS->obj("monojet_signal_data");
  }
  else if(category == 2){
    data = (RooDataHist*) dataWS->obj("data_obs_SR_MV");
    if(data == 0 or data == NULL)
      data = (RooDataHist*) dataWS->obj("monov_signal_data");
  }

  ///////
  return (TH1F*) data->createHistogram("data_obs",*met_obs); 
}

// import signal template --> total signl
TH1F* importSignalHistogram( const vector<string> & signalROOTFile, const string & signalID){

  TH1F*  signal  = NULL;
  for(auto name :  signalROOTFile){
    TFile* file = TFile::Open(name.c_str(),"READ");
    RooWorkspace* ws = (RooWorkspace*) file->Get("combinedws");
    RooDataHist*  histo = NULL;
    RooRealVar*   met = NULL;
    met = (RooRealVar*) ws->obj("met_monojet");
    if(met == NULL or met == 0)
      met = (RooRealVar*) ws->obj("met_monov");

    histo = (RooDataHist*) ws->obj(("monojet_signal_signal_"+signalID).c_str());
    if(histo == NULL or histo == 0)
      histo = (RooDataHist*) ws->obj(("monov_signal_signal_"+signalID).c_str());
    
    if(debug){
      cout<<"file name "<<name<<" histo "<<histo<<" met "<<met<<endl;
    }
    
    if(signal == NULL or signal == 0){
      signal = (TH1F*) histo->createHistogram(histo->GetName(),*met);
      if(debug)
	cout<<"signal integral "<<signal->Integral()<<endl;
    }
    else{
      TH1F* temp = (TH1F*) histo->createHistogram(histo->GetName(),*met);
      signal->Add(temp);
      if(debug)
	cout<<"signal integral "<<signal->Integral()<<endl;
    }
  }


  return signal;
}

/// main function to run the analysis
void limitSimplifiedLikelihood(string dataWorkspace, string combineMLFitRootFile, vector<string> signalROOTFile, string signalID, int category, string outputDirectory){

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc493/libHiggsAnalysisCombinedLimit.so");

  // workspace with data distribution
  TFile* dataFile = TFile::Open(dataWorkspace.c_str(),"READ");
  TH1F* data_obs  = NULL;
  data_obs  = importDataHistogram(dataFile,category); 
  
  if(debug){   
    for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
      cout<<"Data OBS: iBin "<<iBin+1<<" low edge "<<data_obs->GetXaxis()->GetBinLowEdge(iBin+1)<<" upper edge "<<data_obs->GetXaxis()->GetBinLowEdge(iBin+2)<<" content "<<data_obs->GetBinContent(iBin+1)<<endl;
    }
  }

  // Create a RooDataSet for the observed events
  RooCategory* templateBin = new RooCategory("bin","");
  RooRealVar*  observation =  new RooRealVar("observed","",1);
  
  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    templateBin->defineType(Form("bin_%d",iBin+1),iBin+1);
    templateBin->setIndex(iBin+1);
  }

  if(debug){
    templateBin->Print();
  }
 
  RooArgSet*  obsargset = new RooArgSet(*observation,*templateBin);
  RooDataSet* obsdata   = new RooDataSet("combinedData","Data in all Bins",*obsargset);

  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    observation->setVal(data_obs->GetBinContent(iBin+1));
    templateBin->setIndex(iBin+1);
    obsdata->add(RooArgSet(*observation,*templateBin));
  }

  //signal template
  TH1F*  signal  = NULL;
  signal         = importSignalHistogram(signalROOTFile,signalID);

  // Signal
  RooRealVar* mu = new RooRealVar("mu","",0.,-1e4,1e4);
  RooArgList obsSignal;
  for(int iBin = 0; iBin < signal->GetNbinsX(); iBin++){
    RooFormulaVar *mean = new RooFormulaVar(Form("signal_exp_%d",iBin+1),Form("%f*@0",signal->GetBinContent(iBin+1)),RooArgList(*mu));
    if(debug)
      cout<<"Signal exp: iBin "<<iBin<<" low edge "<<signal->GetXaxis()->GetBinLowEdge(iBin+1)<<" upper edge "<<signal->GetXaxis()->GetBinLowEdge(iBin+2)<<" content "<<signal->GetBinContent(iBin+1)<<" value "<<mean->getVal()<<endl;
    obsSignal.add(*mean);
  }

  if(debug)
    obsSignal.Print();
 
  
  //take histograms from combine ML fit
  TFile* combineMLFitFile = TFile::Open(combineMLFitRootFile.c_str(),"READ");  
  TH1F* background = NULL;
  if(category == 1)
    background = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_background"); // take post fit background
  else if(category == 2)
    background = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch2_ch1/total_background"); // take post fit background

  if(background == NULL or background == 0){
    if(category == 1)
      background = (TH1F*) combineMLFitFile->Get("shapes_fit_b/monojet_signal/total_background");
    else if(category == 2)
      background = (TH1F*) combineMLFitFile->Get("shapes_fit_b/monov_signal/total_background");
  }  
  multipleByBinWidth(background);

  TH2F* correlation = NULL;
  if(category == 1)
    correlation = (TH2F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_covar");
  else if(category == 2)
    correlation = (TH2F*) combineMLFitFile->Get("shapes_fit_b/ch2_ch1/total_covar");

  if(correlation == NULL or correlation == 0){
    if(category == 1)
      correlation = (TH2F*) combineMLFitFile->Get("shapes_fit_b/monojet_signal/total_covar");
    else
      correlation = (TH2F*) combineMLFitFile->Get("shapes_fit_b/monov_signal/total_covar");
  }
  
  multipleByBinWidth(correlation,background);
  
  //Create a vector with: central background prediction fixed --> one RooRealVar per bin, one nuisance per bin, one uncertainty per bin
  RooArgList centralBkg;
  RooArgList centralBkgAsimov;
  RooArgList obsBkg;
  TMatrixDSym covariance(background->GetNbinsX());

  for(int iBin = 0; iBin < background->GetNbinsX(); iBin++){

    RooRealVar *mean = new RooRealVar(Form("exp_bin_%d_In",iBin+1),"",background->GetBinContent(iBin+1));
    mean->setConstant(kTRUE);
    centralBkg.add(*mean);

    RooRealVar *meanAsimov = new RooRealVar(Form("exp_bin_%d_In_Asimov",iBin+1),"",background->GetBinContent(iBin+1));
    meanAsimov->setConstant(kTRUE);
    centralBkgAsimov.add(*meanAsimov);

    RooRealVar* nuis = new RooRealVar(Form("exp_bin_%d",iBin+1),"",background->GetBinContent(iBin+1),background->GetBinContent(iBin+1)/10,background->GetBinContent(iBin+1)*10);
    obsBkg.add(*nuis);

    for(int jBin = 0; jBin < background->GetNbinsX(); jBin++){
      covariance[iBin][jBin] = correlation->GetBinContent(iBin+1,jBin+1);
    }

    if(debug)
      cout<<"Background pre-fit: iBin "<<iBin<<" low edge "<<background->GetXaxis()->GetBinLowEdge(iBin+1)<<" upper edge "<<background->GetXaxis()->GetBinLowEdge(iBin+2)<<" content "<<mean->getVal()<<endl;
    
  }

  if(debug){
    centralBkg.Print();
    obsBkg.Print();
  }
  
  //Create the likelihood --> one poisson for each bin
  RooArgList binSum;
  RooArgList binPoission;
  templateBin->setIndex(1);
  RooSimultaneous combined_pdf("combined_pdf","",*templateBin); // one category for each bin and one pdf for each bin
 
  for(int iBin = 0; iBin < signal->GetNbinsX(); iBin++){
    RooAddition *sum  = new RooAddition(Form("splusb_bin_%d",iBin+1),"",RooArgList(*((RooRealVar*)(obsSignal.at(iBin))),*((RooRealVar*)(obsBkg.at(iBin)))));
    if(debug)
      cout<<"iBin "<<iBin+1<<" signal "<<((RooRealVar*)(obsSignal.at(iBin)))->getVal()<<" bkg "<<((RooRealVar*)(obsBkg.at(iBin)))->getVal()<<" sum "<<sum->getVal()<<" S/sqrt(B) "<<((RooRealVar*)(obsSignal.at(iBin)))->getVal()/sqrt(((RooRealVar*)(obsBkg.at(iBin)))->getVal())<<endl;
    RooPoisson  *pois = new RooPoisson(Form("pdf_bin_%d",iBin+1),"",*observation,*sum); 
    combined_pdf.addPdf(*pois,Form("bin_%d",iBin+1));
    binSum.add(*sum);
    binPoission.add(*pois);
  }
  
  if(debug){
    binSum.Print();
    binPoission.Print();
    combined_pdf.Print();
  }
  
  RooMultiVarGaussian constraint_pdf("constraint_pdf","",obsBkg,centralBkg,covariance);
  RooAbsReal *nllSB = combined_pdf.createNLL(*obsdata,RooFit::ExternalConstraints(RooArgList(constraint_pdf)),RooFit::Verbose(-1));
  if(debug){
    constraint_pdf.Print();
    nllSB->Print();
  }

  // make S+B fit
  RooMinimizer mSB(*nllSB);
  if(not debug){
    mSB.setVerbose(kFALSE);
    mSB.setPrintLevel(-1);
  }
  mSB.minimize("Minuit2","minimizer");
  if(debug)
    cout<<"##### S+B fit Minimizer step : -Log(L) "<<nllSB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<<endl;
  mSB.minimize("Minuit2","migrad");
  if(debug)
    cout<<"##### S+B fit Migrad step : -Log(L) "<<nllSB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<<endl;
  mSB.minimize("Minuit2","migradimproved");
  if(debug)
    cout<<"##### S+B fit Migrad improved step : -Log(L) "<<nllSB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<< endl;
  mSB.minos(RooArgSet(*mu));
  cout<<"##### S+B fit Minos step : -Log(L) "<<nllSB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<<endl;
  double nllFit_SB = nllSB->getVal();
  double muFit_SB  = mu->getVal();
  double muFit_SB_errUp  = mu->getErrorHi();
  double muFit_SB_errDw  = mu->getErrorLo();

  // fix muMin and muMax dynamically
  if(muFit_SB < 0) muMin = 0.001;
  else if( muFit_SB + 10*muFit_SB_errDw < 0) muMin = 0.001;
  else muMin = muFit_SB +  10*muFit_SB_errDw;
  muMax = fabs(muFit_SB) + 10*muFit_SB_errUp;

  if(debug){
    cout<<"########### After S+B fit ############"<<endl;
    combined_pdf.Print();
    constraint_pdf.Print();
    cout<<"nllFit "<<nllFit_SB<<" muFit "<<muFit_SB<<" + "<<muFit_SB_errUp<<" - "<<muFit_SB_errDw<<endl;
    cout<<"######################################"<<endl;
  }

  //Make B-only fit and prepare asimov toy
  mu->setVal(0);
  mu->setConstant(true);  
  resetRooArgList(obsBkg,centralBkg);

  RooAbsReal *nllB = combined_pdf.createNLL(*obsdata,RooFit::ExternalConstraints(RooArgList(constraint_pdf)),RooFit::Verbose(-1));  
  RooMinimizer mB(*nllB);
  if(not debug){
    mB.setVerbose(kFALSE);
    mB.setPrintLevel(-1);
  }
  mB.minimize("Minuit2","minimizer");
  if(debug)
    cout<<"##### B fit Minimizer step : -Log(L) "<<nllB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<<endl;
  mB.minimize("Minuit2","migrad");
  if(debug)
    cout<<"##### B fit Migrad step : -Log(L) "<<nllB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<<endl;
  mB.minimize("Minuit2","migradimproved");
  cout<<"##### B fit Migrad improved step : -Log(L) "<<nllB->getVal()<<" mu : "<<mu->getVal()<<" err dw "<<mu->getErrorLo()<<" err up "<<mu->getErrorHi()<< endl;
  double nllFit_B = nllB->getVal();
  
  if(debug){
    cout<<"########### After B-only fit ############"<<endl;
    combined_pdf.Print();
    constraint_pdf.Print();
    cout<<"nllFit "<<nllFit_B<<endl;
    cout<<"#########################################"<<endl;
  }

  //Build asimov data
  RooDataSet asimovdata("AsimovData","",*obsargset);
  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    templateBin->setIndex(iBin+1);
    // make integer to mimic data
    double expectation = double(int(((RooRealVar*) binSum.at(iBin))->getVal()+0.5));
    if(debug)
      cout<<"Asimov daaset: ibin "<<iBin<<" sum S+B, with S=0 "<<((RooRealVar*) binSum.at(iBin))->getVal()<<" expectation "<<expectation<<endl;
    ((RooRealVar*) centralBkgAsimov.at(iBin))->setVal(expectation);
    observation->setVal(expectation);
    asimovdata.add(RooArgSet(*observation,*templateBin));
  }
    
  RooMultiVarGaussian asimov_constraint_pdf("asimov_constraint_pdf","",obsBkg,centralBkgAsimov,covariance);
  RooAbsReal *nllAsimov = combined_pdf.createNLL(asimovdata,RooFit::ExternalConstraints(RooArgList(asimov_constraint_pdf)),RooFit::Verbose(-1)); // no need to refit, we have done the b-only fit estimating parameters, they toy is just fixing the central value without randomizing and re-doing the fit

  if(debug){
    cout<<"########### Asimov dataset ############"<<endl;
    asimovdata.Print();
    asimov_constraint_pdf.Print();
    combined_pdf.Print();
    cout<<"#######################################"<<endl;
  }
    
  system(("mkdir -p "+outputDirectory).c_str());

  string cat;
  if (category == 1)
    cat = "monojet";
  else if(category == 2)
    cat = "monov";

  TFile* outputFile = new TFile((outputDirectory+"/simplifiedLikelihood_"+cat+"_"+signalID+".root").c_str(),"RECREATE");
  outputFile->cd();

  double mh      = 0;
  float mu_      = 0;
  float muErrUp_ = 0;
  float muErrDw_ = 0;
  float NLL_SB   = 0;
  float NLL_B    = 0;
  float CLsObs   = 0;
  float CLsExp   = 0;
  float limitObs = 0;
  float limitExp = 0;
  float limitExpUp1s = 0;
  float limitExpUp2s = 0;
  float limitExpDw1s = 0;
  float limitExpDw2s = 0;
  vector<float> nllObserved = {};
  vector<float> nllExpected = {};
  vector<float> muVal = {};

  TTree* tree = new TTree("limit","limit");
  tree->Branch("mh",&mh,"mh/D");
  tree->Branch("mu",&mu_,"mu/F");
  tree->Branch("muErrUp",&muErrUp_,"muErrUp/F");
  tree->Branch("muErrDw",&muErrDw_,"muErrDw/F");
  tree->Branch("NLL_SB",&NLL_SB,"NLL_SB/F");
  tree->Branch("NLL_B",&NLL_B,"NLL_B/F");
  tree->Branch("CLsObs",&CLsObs,"CLsObs/F");
  tree->Branch("CLsExp",&CLsExp,"CLsExp/F");
  tree->Branch("limitObs",&limitObs,"limitObs/F");
  tree->Branch("limitExp",&limitExp,"limitExp/F");
  tree->Branch("limitExpUp1s",&limitExpUp1s,"limitExpUp1s/F");
  tree->Branch("limitExpUp2s",&limitExpUp2s,"limitExpUp2s/F");
  tree->Branch("limitExpDw1s",&limitExpDw1s,"limitExpDw1s/F");
  tree->Branch("limitExpDw2s",&limitExpDw2s,"limitExpDw2s/F");
  tree->Branch("muVal","vector<float>",&muVal);
  tree->Branch("nllObserved","vector<float>",&nllObserved);
  tree->Branch("nllExpected","vector<float>",&nllExpected);
    
  mh     = atof(signalID.c_str());
  NLL_SB = nllFit_SB;
  NLL_B  = nllFit_B;
  mu_    = muFit_SB;
  muErrUp_ = muFit_SB_errUp;
  muErrDw_ = muFit_SB_errDw;

  double nllData   = 0;
  double nllA      = 0;

  CLsObs = calculateCLsLimit(nllSB,nllAsimov,mu,1,true,nllData,nllA);
  cout<<"#### -Log(L) data : "<<nllData<<" -Log(L) Asimov : "<<nllA+NLL_B<<endl;
  CLsExp = calculateCLsLimit(nllSB,nllAsimov,mu,1,false,nllData,nllA,0.5);
  cout<<"#### Get observed and expected CLs values : observed "<<CLsObs<<" expected "<<CLsExp<<endl;
  
  // Calculate upper limit on signal strenght --> find CLs value at 95% of the range -> associated mu --> store also the nLL for Asimov and Data --> Likelihood scan
  double nllDataScan = 0;
  double nllAsimovScan = 0;
  double clsObs = 0;
  double clsExp = 0;
  double clsExp1sUp = 0;
  double clsExp1sDw = 0;
  double clsExp2sUp = 0;
  double clsExp2sDw = 0;

  TGraph* graphObserved = new TGraph();
  TGraph* graphExpected = new TGraph();
  TGraph* graphExpected1sUp = new TGraph();
  TGraph* graphExpected2sUp = new TGraph();
  TGraph* graphExpected1sDw = new TGraph();
  TGraph* graphExpected2sDw = new TGraph();

  cout<<"#### Start Scan for limit computation : muMin "<<muMin<<" muMax "<<muMax<<endl;  
  for(int istep = 0; istep <= muPoint; istep++){
    clsObs = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,true,nllDataScan,nllAsimovScan);    
    clsExp = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.5);
    if(debug)
      cout<<"#### Step: "<<istep<<" muMin "<<muMin<<" muMax "<<muMax<<" muVal "<<muMin+istep*(muMax-muMin)/muPoint<<" ---> nllObs "<<nllDataScan-nllData<<" nllExp "<<nllAsimovScan-nllA-NLL_B<<" clsObs "<<clsObs<<" clsExp "<<clsExp<<endl;
    graphObserved->SetPoint(istep,clsObs,muMin+istep*(muMax-muMin)/muPoint);
    graphExpected->SetPoint(istep,clsExp,muMin+istep*(muMax-muMin)/muPoint);
    if(calculateExpSigma){
      clsExp1sUp = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.84);
      graphExpected1sUp->SetPoint(istep,clsExp1sUp,muMin+istep*(muMax-muMin)/muPoint);
      clsExp2sUp = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.975);
      graphExpected2sUp->SetPoint(istep,clsExp2sUp,muMin+istep*(muMax-muMin)/muPoint);
      clsExp1sDw = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.16);
      graphExpected1sDw->SetPoint(istep,clsExp1sDw,muMin+istep*(muMax-muMin)/muPoint);
      clsExp2sDw = calculateCLsLimit(nllSB,nllAsimov,mu,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.025);
      graphExpected2sDw->SetPoint(istep,clsExp2sDw,muMin+istep*(muMax-muMin)/muPoint);
    }
  }

  limitObs = graphObserved->Eval(1.-0.95);
  limitExp = graphExpected->Eval(1.-0.95);

  cout<<"##### limitObs mu < "<<limitObs<<endl;
  cout<<"##### limitExp mu < "<<limitExp<<endl;

  if(calculateExpSigma){
    limitExpDw2s = graphExpected2sDw->Eval(1-0.95);
    cout<<"##### limitExp 2sigma dw mu < "<<limitExpDw2s<<endl;
    limitExpDw1s = graphExpected1sDw->Eval(1-0.95);
    cout<<"##### limitExp 1sigma dw mu < "<<limitExpDw1s<<endl;
    limitExpUp1s = graphExpected1sUp->Eval(1-0.95);
    cout<<"##### limitExp 1sigma up mu < "<<limitExpUp1s<<endl;
    limitExpUp2s = graphExpected2sUp->Eval(1-0.95);
    cout<<"##### limitExp 2sigma up mu < "<<limitExpUp2s<<endl;

  }

  // make likelihood scan
  if(doLikelihoodScan){
    cout<<"Make Likelihood Scan  ---> "<<endl;
    nllA = makeFitFixedStrenght(nllAsimov,mu,0);
    for(int istep = 0; istep <= muPoint; istep++){
      muVal.push_back(-muMax+istep*(muMax+muMax)/muPoint);
      nllDataScan = makeFitFixedStrenght(nllSB,mu,muVal.back());
      nllAsimovScan = makeFitFixedStrenght(nllAsimov,mu,muVal.back());
      nllObserved.push_back(nllDataScan-nllData);
      nllExpected.push_back(nllAsimovScan-nllA);
    }
  }

  tree->Fill();
  tree->Write();
  outputFile->Close();
  
}
