//////// Basic options
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
float calculateCLsLimit(RooAbsReal *nllSB, RooAbsReal *nllAsimov, RooWorkspace* ws_sb, RooWorkspace* ws_asimov, double muVal, bool isObserved, double & nllData, double & nllA, const double & quantile = 0){

  RooRealVar* mu = ws_asimov->var("mu");
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
  mu = ws_sb->var("mu");
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
  mu = ws_asimov->var("mu");
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
double makeFitFixedStrenght(RooAbsReal *nll, RooWorkspace *ws, double muVal){

  RooRealVar* mu = ws->var("mu");
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
TH1F* importSignalHistogram( const vector<string> & signalROOTFile, const string & signalID, const int & category){

  TH1F*  signal  = NULL;
  for(auto name :  signalROOTFile){
    TFile* file = TFile::Open(name.c_str(),"READ");
    RooWorkspace* ws = (RooWorkspace*) file->Get("combinedws");
    RooDataHist*  histo = NULL;
    RooRealVar*   met = NULL;
    if(category == 1)
      met = (RooRealVar*) ws->obj("met_monojet");
    else if(category == 2)
      met = (RooRealVar*) ws->obj("met_monov");

    if(category == 1)
      histo = (RooDataHist*) ws->obj(("monojet_signal_signal_"+signalID).c_str());
    else if(category == 2)
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


// make a likelihood for a given event category
void makeLikelihood (RooWorkspace* model_sb, RooWorkspace* model_b, TH1* data_obs, TH1* signal, TH1* background, TH2* correlation, const string & postfix = ""){

  // Create a RooDataSet for the observed events
  RooCategory* templateBin = new RooCategory(Form("bin_%s",postfix.c_str()),"");
  RooRealVar*  observation = new RooRealVar(Form("observed_%s",postfix.c_str()),"",1);
  
  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    templateBin->defineType(Form("bin_%s_%d",postfix.c_str(),iBin+1),iBin+1);
    templateBin->setIndex(iBin+1);
  }
  
  if(debug){
    templateBin->Print();
  }
  
  RooArgSet*  obsargset = new RooArgSet(*observation,*templateBin);
  obsargset->setName(Form("obs_%s",postfix.c_str()));
  
  RooDataSet* obsdata   = new RooDataSet(Form("combinedData_%s",postfix.c_str()),"Data in all Bins",*obsargset);  
  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    observation->setVal(data_obs->GetBinContent(iBin+1));
    templateBin->setIndex(iBin+1);
    obsdata->add(RooArgSet(*observation,*templateBin));
  }

  // imports
  model_sb->import(*obsdata);
  model_b->import(*obsdata);
  
  RooArgList obsSignal;
  obsSignal.setName(Form("obsSignal_%s",postfix.c_str()));
  RooRealVar* mu = NULL;
  if(model_sb->var("mu"))
    mu = model_sb->var("mu");
  else
    mu = new RooRealVar("mu","",0.,-1e4,1e4);

  for(int iBin = 0; iBin < signal->GetNbinsX(); iBin++){
    RooFormulaVar *mean = new RooFormulaVar(Form("signal_exp_%s_%d",postfix.c_str(),iBin+1),Form("%f*@0",signal->GetBinContent(iBin+1)),RooArgList(*mu));
    if(debug)
      cout<<"Signal exp: iBin "<<iBin<<" low edge "<<signal->GetXaxis()->GetBinLowEdge(iBin+1)<<" upper edge "<<signal->GetXaxis()->GetBinLowEdge(iBin+2)<<" content "<<signal->GetBinContent(iBin+1)<<" value "<<mean->getVal()<<endl;
    obsSignal.add(*mean);
  }

  //Create a vector with: central background prediction fixed --> one RooRealVar per bin, one nuisance per bin, one uncertainty per bin
  RooArgList centralBkg;
  centralBkg.setName(Form("centralBkg_%s",postfix.c_str()));
  RooArgList centralBkgAsimov;
  centralBkgAsimov.setName(Form("centralBkgAsimov_%s",postfix.c_str()));
  RooArgList obsBkg;
  obsBkg.setName(Form("obsBkg_%s",postfix.c_str()));
  TMatrixDSym covariance(background->GetNbinsX());
  
  for(int iBin = 0; iBin < background->GetNbinsX(); iBin++){
    
    RooRealVar *mean = new RooRealVar(Form("exp_bin_%s_%d_In",postfix.c_str(),iBin+1),"",background->GetBinContent(iBin+1));
    mean->setConstant(kTRUE);
    centralBkg.add(*mean);
    
    RooRealVar *meanAsimov = new RooRealVar(Form("exp_bin_%s_%d_In_Asimov",postfix.c_str(),iBin+1),"",background->GetBinContent(iBin+1));
    meanAsimov->setConstant(kTRUE);
    centralBkgAsimov.add(*meanAsimov);

    RooRealVar* nuis = new RooRealVar(Form("exp_bin_%s_%d",postfix.c_str(),iBin+1),"",background->GetBinContent(iBin+1),background->GetBinContent(iBin+1)/10,background->GetBinContent(iBin+1)*10);
    obsBkg.add(*nuis);

    model_sb->import(*meanAsimov);
    model_b->import(*meanAsimov);

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
  binSum.setName(Form("binSum_%s",postfix.c_str()));
  RooArgList binPoission;
  templateBin->setIndex(1);
  RooSimultaneous* combined_pdf = new RooSimultaneous(Form("combined_pdf_%s",postfix.c_str()),"",*templateBin); // one category for each bin and one pdf for each bin
  

  for(int iBin = 0; iBin < signal->GetNbinsX(); iBin++){
    RooAddition *sum  = new RooAddition(Form("splusb_bin_%s_%d",postfix.c_str(),iBin+1),"",RooArgList(*((RooRealVar*)(obsSignal.at(iBin))),*((RooRealVar*)(obsBkg.at(iBin)))));
    if(debug)
      cout<<"iBin "<<iBin+1<<" signal "<<((RooRealVar*)(obsSignal.at(iBin)))->getVal()<<" bkg "<<((RooRealVar*)(obsBkg.at(iBin)))->getVal()<<" sum "<<sum->getVal()<<" S/sqrt(B) "<<((RooRealVar*)(obsSignal.at(iBin)))->getVal()/sqrt(((RooRealVar*)(obsBkg.at(iBin)))->getVal())<<endl;
    observation->setVal(data_obs->GetBinContent(iBin+1));
    RooPoisson  *pois = new RooPoisson(Form("pdf_bin_%s_%d",postfix.c_str(),iBin+1),"",*observation,*sum); 
    combined_pdf->addPdf(*pois,Form("bin_%s_%d",postfix.c_str(),iBin+1));
    binSum.add(*sum);
    binPoission.add(*pois);  
  }
  
  if(debug){
    binSum.Print();
    binPoission.Print();
    combined_pdf->Print();
  }

  model_sb->import(*combined_pdf);
  model_b->import(*combined_pdf);

  RooMultiVarGaussian* constraint_pdf = new RooMultiVarGaussian(Form("constraint_%s_pdf",postfix.c_str()),"",obsBkg,centralBkg,covariance);
  // S+B likelihood
  RooAbsReal* nllSB = combined_pdf->createNLL(*obsdata,RooFit::ExternalConstraints(RooArgList(*constraint_pdf)),RooFit::Verbose(-1));
  nllSB->SetName(Form("nllSB_%s",postfix.c_str()));

  //Make B-only fit and prepare asimov toy
  resetRooArgList(obsBkg,centralBkg);
  RooAbsReal* nllB = combined_pdf->createNLL(*obsdata,RooFit::ExternalConstraints(RooArgList(*constraint_pdf)),RooFit::Verbose(-1));  
  nllB->SetName(Form("nllB_%s",postfix.c_str()));
  model_sb->import(*nllSB);
  model_b->import(*nllB);

  // define sets
  model_sb->defineSet(Form("obs_%s",postfix.c_str()),*obsargset);
  model_b->defineSet(Form("obs_%s",postfix.c_str()),*obsargset);
  model_sb->defineSet(Form("obsSignal_%s",postfix.c_str()),obsSignal);
  model_b->defineSet(Form("obsSignal_%s",postfix.c_str()),obsSignal);
  model_sb->defineSet(Form("obsBkg_%s",postfix.c_str()),obsBkg);
  model_b->defineSet(Form("obsBkg_%s",postfix.c_str()),obsBkg);
  model_sb->defineSet(Form("centralBkg_%s",postfix.c_str()),centralBkg);
  model_b->defineSet(Form("centralBkg_%s",postfix.c_str()),centralBkg);
  model_sb->defineSet(Form("centralBkgAsimov_%s",postfix.c_str()),centralBkgAsimov);
  model_b->defineSet(Form("centralBkgAsimov_%s",postfix.c_str()),centralBkgAsimov);
  model_sb->defineSet(Form("binSum_%s",postfix.c_str()),binSum);
  model_b->defineSet(Form("binSum_%s",postfix.c_str()),binSum);  

  if(debug){
    constraint_pdf->Print();
    nllSB->Print();
    nllB->Print();
    model_sb->Print();
    model_b->Print();
  }
  
}

// in order to make the asimov likelihood
void  makeAsimovLikelihood(RooWorkspace* ws_asimov, TH1* data_obs, RooWorkspace* ws, TMatrixDSym* covariance, const string & postfix = ""){

  RooCategory* templateBin = (RooCategory*) ws->obj(Form("bin_%s",postfix.c_str()));
  RooRealVar*  observation = ws->var(Form("observed_%s",postfix.c_str()));
  RooArgList*  binSum      = (RooArgList*)  ws->set(Form("binSum_%s",postfix.c_str()));
  RooArgList*  centralBkgAsimov = (RooArgList*) ws->set(Form("centralBkgAsimov_%s",postfix.c_str())); 
  RooArgList*  obsBkg      = (RooArgList*) ws->set(Form("obsBkg_%s",postfix.c_str())); 
  RooDataSet*  asimovdata  = new RooDataSet(Form("AsimovData_%s",postfix.c_str()),"",*ws->set(Form("obs_%s",postfix.c_str())));

  for(int iBin = 0; iBin < data_obs->GetNbinsX(); iBin++){
    templateBin->setIndex(iBin+1);
    // make integer to mimic data
    double expectation = double(int(((RooRealVar*) binSum->at(iBin))->getVal()+0.5));
    if(debug)
      cout<<"Asimov daaset: ibin "<<iBin<<" sum S+B, with S=0 "<<((RooRealVar*) binSum->at(iBin))->getVal()<<" expectation "<<expectation<<endl;
    ((RooRealVar*) centralBkgAsimov->at(iBin))->setVal(expectation);
    observation->setVal(expectation);
    asimovdata->add(RooArgSet(*observation,*templateBin));
  }
  ws_asimov->import(*asimovdata);
    
  RooMultiVarGaussian asimov_constraint_pdf(Form("asimov_constraint_%s_pdf",postfix.c_str()),"",*obsBkg,*centralBkgAsimov,*covariance);
  RooSimultaneous* combined_pdf = (RooSimultaneous*) ws->pdf(Form("combined_pdf_%s",postfix.c_str()));
  ws_asimov->import(*combined_pdf);
  RooAbsReal* nll = combined_pdf->createNLL(*asimovdata,RooFit::ExternalConstraints(RooArgList(asimov_constraint_pdf)),RooFit::Verbose(-1));
  nll->SetName(Form("nllAsimov_%s",postfix.c_str()));
  ws_asimov->import(*nll);

}


/// main function to run the analysis
void limitSimplifiedLikelihood(string dataWorkspace, string combineMLFitRootFile, vector<string> signalROOTFileMonoJet, vector<string> signalROOTFileMonoV, string signalID, int category, string outputDirectory){

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc493/libHiggsAnalysisCombinedLimit.so");

  // workspace with data distribution
  TFile* dataFile         = TFile::Open(dataWorkspace.c_str(),"READ");
  TH1F* data_obs_monojet  = NULL;
  TH1F* data_obs_monov    = NULL;
  if(category == 0){ // take both mono-jet and mono-v data
    data_obs_monojet  = importDataHistogram(dataFile,1); 
    data_obs_monov    = importDataHistogram(dataFile,2); 
  }
  else if(category == 1)
    data_obs_monojet  = importDataHistogram(dataFile,category); 
  else if(category == 2)
    data_obs_monov    = importDataHistogram(dataFile,category); 
 
  //signal template
  TH1F*  signal_monojet  = NULL;
  TH1F*  signal_monov    = NULL;
  if(category == 0){// take both mono-jet and mono-v
    signal_monojet    = importSignalHistogram(signalROOTFileMonoJet,signalID,1);
    signal_monov      = importSignalHistogram(signalROOTFileMonoV,signalID,2);
  }
  else if(category == 1)
    signal_monojet  = importSignalHistogram(signalROOTFileMonoJet,signalID,category);
  else if(category == 2)
    signal_monov    = importSignalHistogram(signalROOTFileMonoV,signalID,category);
     
  //take histograms from combine ML fit
  TFile* combineMLFitFile = TFile::Open(combineMLFitRootFile.c_str(),"READ");  
  TH1F* background_monojet = NULL;
  TH1F* background_monov   = NULL;
  if(category == 0){
    background_monojet = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_background"); // take post fit background
    background_monov   = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_background"); // take post fit background
  }
  else if(category == 1)
    background_monojet = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_background"); // take post fit background
  else if(category == 2)
    background_monov   = (TH1F*) combineMLFitFile->Get("shapes_fit_b/ch2_ch1/total_background"); // take post fit background

  if(background_monojet == NULL or background_monojet == 0){
    if(category == 0 or category == 1)
      background_monojet = (TH1F*) combineMLFitFile->Get("shapes_fit_b/monojet_signal/total_background");
  }  
  
  if(background_monov == NULL or background_monov == 0){
    if(category == 0 or category == 2)
      background_monov = (TH1F*) combineMLFitFile->Get("shapes_fit_b/monov_signal/total_background");
  }  
  
  if(background_monojet)
    multipleByBinWidth(background_monojet);
  if(background_monov)
    multipleByBinWidth(background_monov);

  TH2F* correlation_monojet = NULL;
  TH2F* correlation_monov   = NULL;
  if(category == 0 or category == 1)
    correlation_monojet = (TH2F*) combineMLFitFile->Get("shapes_fit_b/ch1_ch1/total_covar");
  if(category == 0 or category == 2)
    correlation_monov   = (TH2F*) combineMLFitFile->Get("shapes_fit_b/ch2_ch1/total_covar");

  if(correlation_monojet == NULL or correlation_monojet == 0){
    if(category == 0 or category == 1)
      correlation_monojet = (TH2F*) combineMLFitFile->Get("shapes_fit_b/monojet_signal/total_covar");
  }
  if(correlation_monov == NULL or correlation_monov == 0){
    if(category == 0 or category == 2)
      correlation_monov = (TH2F*) combineMLFitFile->Get("shapes_fit_b/monov_signal/total_covar");
  }
  if(background_monojet)
    multipleByBinWidth(correlation_monojet,background_monojet);
  if(background_monov)
    multipleByBinWidth(correlation_monov,background_monov);

  // Build the individual likelihoods for each category
  RooAbsReal* nllSB_monojet = NULL;
  RooAbsReal* nllB_monojet  = NULL;
  RooAbsReal* nllSB_monov   = NULL;
  RooAbsReal* nllB_monov    = NULL;
  RooWorkspace* ws_sb  = NULL;
  RooWorkspace* ws_b   = NULL;
  if(category == 0 or category == 1){
    ws_sb = new RooWorkspace("ws_sb","ws_sb");
    ws_b  = new RooWorkspace("ws_b","ws_b");
    makeLikelihood(ws_sb,ws_b,data_obs_monojet,signal_monojet,background_monojet,correlation_monojet,"monojet");
    nllSB_monojet = (RooAbsReal*) ws_sb->obj("nllSB_monojet");
    nllB_monojet = (RooAbsReal*) ws_b->obj("nllB_monojet");
  }
  
  if(category == 0 or category == 2){
    if(ws_sb == NULL)
      ws_sb = new RooWorkspace("ws_sb","ws_sb");
    if(ws_b == NULL)
    ws_b  = new RooWorkspace("ws_b","ws_b");
    makeLikelihood(ws_sb,ws_b,data_obs_monov,signal_monov,background_monov,correlation_monov,"monov");
    nllSB_monov = (RooAbsReal*) ws_sb->obj("nllSB_monov");
    nllB_monov  = (RooAbsReal*) ws_b->obj("nllB_monov");
  }

  RooAddition* nllSB = NULL;
  RooAddition* nllB  = NULL;

  if(category == 1){
    nllSB = dynamic_cast<RooAddition*>(nllSB_monojet);
    nllB  = dynamic_cast<RooAddition*>(nllB_monojet);
  }
  else if(category == 2){
    nllSB = dynamic_cast<RooAddition*>(nllSB_monov);
    nllB  = dynamic_cast<RooAddition*>(nllB_monov);
  }
  else if(category == 0){
    nllSB = new RooAddition("nllSBsum","",RooArgList(*nllSB_monojet,*nllSB_monov));
    nllB = new RooAddition("nllBsum","",RooArgList(*nllB_monojet,*nllB_monov));
  }

  ///////////////// make S+B fit
  RooRealVar* mu = ws_sb->var("mu");
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
    cout<<"nllFit "<<nllFit_SB<<" muFit "<<muFit_SB<<" + "<<muFit_SB_errUp<<" - "<<muFit_SB_errDw<<endl;
    cout<<"######################################"<<endl;
  }

  
  /////////// make b-only fit
  mu = ws_b->var("mu");
  mu->setConstant(true);
  mu->setVal(0);
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
    cout<<"nllFit "<<nllFit_B<<endl;
    cout<<"#########################################"<<endl;
  }
  
  //make covariance
  TMatrixDSym* covariance_monojet = NULL;
  if(correlation_monojet != NULL){
    covariance_monojet =  new TMatrixDSym(background_monojet->GetNbinsX());
    for(int ibin = 0; ibin < correlation_monojet->GetNbinsX(); ibin++){
      for(int jbin = 0; jbin < correlation_monojet->GetNbinsY(); jbin++){
	(*covariance_monojet)[ibin][jbin] = correlation_monojet->GetBinContent(ibin+1,jbin+1);
      }
    }
  }
  TMatrixDSym* covariance_monov = NULL;
  if(correlation_monov != NULL){
    covariance_monov =  new TMatrixDSym(background_monov->GetNbinsX());
    for(int ibin = 0; ibin < correlation_monov->GetNbinsX(); ibin++){
      for(int jbin = 0; jbin < correlation_monov->GetNbinsY(); jbin++){
	(*covariance_monov)[ibin][jbin] = correlation_monov->GetBinContent(ibin+1,jbin+1);
      }
    }
  }

  //Build asimov data
  RooWorkspace* ws_asimov = NULL;
  RooAbsReal* nllAsimov_monojet = NULL;
  RooAbsReal* nllAsimov_monov   = NULL;

  if(category == 0 or category == 1){
    ws_asimov = new RooWorkspace("ws_asimov","ws_asimov");
    makeAsimovLikelihood(ws_asimov,data_obs_monojet,ws_b,covariance_monojet,"monojet");
    nllAsimov_monojet = (RooAbsReal*) ws_asimov->obj("nllAsimov_monojet");
  }
  if(category == 0 or category == 2){
    if(ws_asimov == NULL)
      ws_asimov = new RooWorkspace("ws_asimov","ws_asimov");
    makeAsimovLikelihood(ws_asimov,data_obs_monov,ws_b,covariance_monov,"monov");
    nllAsimov_monov = (RooAbsReal*) ws_asimov->obj("nllAsimov_monov");
  }

  RooAddition* nllAsimov = NULL;
  if(category == 1)
    nllAsimov = dynamic_cast<RooAddition*>(nllAsimov_monojet);
  else if(category == 2)
    nllAsimov = dynamic_cast<RooAddition*>(nllAsimov_monov);
  else if(category == 0)
   nllAsimov = new RooAddition("nllAsimovSum","",RooArgList(*nllAsimov_monojet,*nllAsimov_monov));
    
  system(("mkdir -p "+outputDirectory).c_str());

  string cat;
  if (category == 1)
    cat = "monojet";
  else if(category == 2)
    cat = "monov";
  else if(category == 0)
    cat =  "combined";

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
  
  CLsObs = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,1,true,nllData,nllA);
  cout<<"#### -Log(L) data : "<<nllData<<" -Log(L) Asimov : "<<nllA+NLL_B<<endl;
  CLsExp = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,1,false,nllData,nllA,0.5);
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
    clsObs = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,true,nllDataScan,nllAsimovScan);    
    clsExp = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.5);    
    if(debug)
      cout<<"#### Step: "<<istep<<" muMin "<<muMin<<" muMax "<<muMax<<" muVal "<<muMin+istep*(muMax-muMin)/muPoint<<" ---> nllObs "<<nllDataScan-nllData<<" nllExp "<<nllAsimovScan-nllA-NLL_B<<" clsObs "<<clsObs<<" clsExp "<<clsExp<<endl;
    graphObserved->SetPoint(istep,clsObs,muMin+istep*(muMax-muMin)/muPoint);
    graphExpected->SetPoint(istep,clsExp,muMin+istep*(muMax-muMin)/muPoint);
    
    if(calculateExpSigma){
      clsExp1sUp = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.84);
      graphExpected1sUp->SetPoint(istep,clsExp1sUp,muMin+istep*(muMax-muMin)/muPoint);
      clsExp2sUp = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.975);
      graphExpected2sUp->SetPoint(istep,clsExp2sUp,muMin+istep*(muMax-muMin)/muPoint);
      clsExp1sDw = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.16);
      graphExpected1sDw->SetPoint(istep,clsExp1sDw,muMin+istep*(muMax-muMin)/muPoint);
      clsExp2sDw = calculateCLsLimit(nllSB,nllAsimov,ws_sb,ws_asimov,muMin+istep*(muMax-muMin)/muPoint,false,nllDataScan,nllAsimovScan,0.025);
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
    nllA = makeFitFixedStrenght(nllAsimov,ws_asimov,0);
    for(int istep = 0; istep <= muPoint; istep++){
      muVal.push_back(-muMax+istep*(muMax+muMax)/muPoint);
      nllDataScan   = makeFitFixedStrenght(nllSB,ws_sb,muVal.back());
      nllAsimovScan = makeFitFixedStrenght(nllAsimov,ws_asimov,muVal.back());
      nllObserved.push_back(nllDataScan-nllData);
      nllExpected.push_back(nllAsimovScan-nllA);
    }
  }
  tree->Fill();
  tree->Write();
  outputFile->Close();
}
