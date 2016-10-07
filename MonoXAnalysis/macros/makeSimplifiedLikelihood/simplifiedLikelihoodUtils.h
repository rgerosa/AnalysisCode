// multipy by bin width                                                                                                                                                                             
void multipleByBinWidth (TH1F* hist){

  for(int iBin = 0; iBin < hist->GetNbinsX(); iBin++){
    hist->SetBinContent(iBin+1,hist->GetBinContent(iBin+1)*hist->GetBinWidth(iBin+1));
    hist->SetBinError(iBin+1,hist->GetBinError(iBin+1)*hist->GetBinWidth(iBin+1));
  }
}

void multipleByBinWidth (TH1F* hist, TH1F* hist_monojet, TH1F* hist_monoV){
  int iBin = 0;
  for(; iBin < hist_monojet->GetNbinsX(); iBin++){
    hist->SetBinContent(iBin+1,hist->GetBinContent(iBin+1)*hist_monojet->GetBinWidth(iBin+1));
    hist->SetBinError(iBin+1,hist->GetBinError(iBin+1)*hist_monojet->GetBinWidth(iBin+1));
    if(debug)
      cout<<"Bin monojet "<<iBin<<" content "<<hist->GetBinContent(iBin+1)<<" width "<<hist_monojet->GetBinWidth(iBin+1)<<endl;
  }
  for(int iBinX = 0; iBinX < hist_monoV->GetNbinsX(); iBinX++, iBin++){
    hist->SetBinContent(iBin+1,hist->GetBinContent(iBin+1)*hist_monoV->GetBinWidth(iBinX+1));
    hist->SetBinError(iBin+1,hist->GetBinError(iBin+1)*hist_monoV->GetBinWidth(iBinX+1));
    if(debug)
      cout<<"Bin monoV "<<iBin<<" content "<<hist->GetBinContent(iBin+1)<<" width "<<hist_monoV->GetBinWidth(iBinX+1)<<endl;
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

// multiply by bin width  a covariance from ML fit in combine                                                                                                                                      
void multipleByBinWidth (TH2F* hist, TH1F* histo_monojet, TH1F* histo_monov){

  for(int iBin = 0; iBin < hist->GetNbinsX(); iBin++){
    for(int jBin = 0; jBin < hist->GetNbinsY(); jBin++){
      if(iBin < histo_monojet->GetNbinsX() and jBin < histo_monojet->GetNbinsY()){
	hist->SetBinContent(iBin+1,jBin+1,hist->GetBinContent(iBin+1,jBin+1)*histo_monojet->GetXaxis()->GetBinWidth(iBin+1)*histo_monojet->GetXaxis()->GetBinWidth(jBin+1));
	hist->SetBinError(iBin+1,jBin+1,hist->GetBinError(iBin+1,jBin+1)*histo_monojet->GetXaxis()->GetBinWidth(iBin+1)*histo_monojet->GetXaxis()->GetBinWidth(jBin+1));
      }
      else if(iBin < histo_monojet->GetNbinsX() and jBin >= histo_monojet->GetNbinsY()){
	hist->SetBinContent(iBin+1,jBin+1,hist->GetBinContent(iBin+1,jBin+1)*histo_monojet->GetXaxis()->GetBinWidth(iBin+1)*histo_monov->GetXaxis()->GetBinWidth(jBin+1-histo_monojet->GetNbinsY()));
	hist->SetBinError(iBin+1,jBin+1,hist->GetBinError(iBin+1,jBin+1)*histo_monojet->GetXaxis()->GetBinWidth(iBin+1)*histo_monov->GetXaxis()->GetBinWidth(jBin+1-histo_monojet->GetNbinsY()));
      }
      else if(iBin >= histo_monojet->GetNbinsX() and jBin < histo_monojet->GetNbinsY()){
	hist->SetBinContent(iBin+1,jBin+1,hist->GetBinContent(iBin+1,jBin+1)*histo_monov->GetXaxis()->GetBinWidth(iBin+1-histo_monojet->GetNbinsX())*histo_monojet->GetXaxis()->GetBinWidth(jBin+1));
	hist->SetBinError(iBin+1,jBin+1,hist->GetBinError(iBin+1,jBin+1)*histo_monov->GetXaxis()->GetBinWidth(iBin+1-histo_monojet->GetNbinsX())*histo_monojet->GetXaxis()->GetBinWidth(jBin+1));
      }
      else if(iBin >= histo_monojet->GetNbinsX() and jBin >= histo_monojet->GetNbinsY()){
	hist->SetBinContent(iBin+1,jBin+1,hist->GetBinContent(iBin+1,jBin+1)*histo_monov->GetXaxis()->GetBinWidth(iBin+1-histo_monojet->GetNbinsX())*histo_monov->GetXaxis()->GetBinWidth(jBin+1-histo_monojet->GetNbinsY()));
	hist->SetBinError(iBin+1,jBin+1,hist->GetBinError(iBin+1,jBin+1)*histo_monov->GetXaxis()->GetBinWidth(iBin+1-histo_monojet->GetNbinsX())*histo_monov->GetXaxis()->GetBinWidth(jBin+1-histo_monojet->GetNbinsY()));
      }
      else{
	cerr<<"Problem with binning in total covariance matrix : x-axis bin "<<iBin<<" y-axis bin "<<jBin<<endl;
      }
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

////                                                                                                                                                                                                 
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

  RooWorkspace* dataWS  = NULL;
  RooRealVar*   met_obs = NULL;

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

/////////
TH1F* importMergedDataHistogram(TFile* dataFile){

  ///////                                                                                                                                                                                           
  TH1F* dataMJ = importDataHistogram(dataFile,1);
  dataMJ->SetName("data_obs_monoj");
  TH1F* dataMV = importDataHistogram(dataFile,2);
  dataMJ->SetName("data_obs_monoV");

  TH1F* data = new TH1F("data_obs","",dataMJ->GetNbinsX()+dataMV->GetNbinsX(),0,dataMJ->GetNbinsX()+dataMV->GetNbinsX()+1);
  data->Sumw2();
  int iBin = 0;
  for( ; iBin < dataMJ->GetNbinsX(); iBin++){
    data->SetBinContent(iBin+1,dataMJ->GetBinContent(iBin+1));
    data->SetBinError(iBin+1,dataMJ->GetBinError(iBin+1));
  }	
  for(int iBinX = 0; iBinX < dataMV->GetNbinsX()+1; iBinX++, iBin++){
    data->SetBinContent(iBin+1,dataMV->GetBinContent(iBinX+1));
    data->SetBinError(iBin+1,dataMV->GetBinError(iBinX+1));
  }	

  return data;

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


// import signal template --> total signl                                                                                                                                                        
TH1F* importMergedSignalHistogram( const vector<string> & signalMonojetROOTFile, const vector<string> & signalMonoVROOTFile,  const string & signalID){

  TH1F*  signal_monojet  = importSignalHistogram(signalMonojetROOTFile,signalID,1);
  TH1F*  signal_monoV    = importSignalHistogram(signalMonoVROOTFile,signalID,2);

  TH1F* signal = new TH1F("signal_signal","",signal_monojet->GetNbinsX()+signal_monoV->GetNbinsX(),0,signal_monojet->GetNbinsX()+signal_monoV->GetNbinsX()+1);
  signal->Sumw2();
  int iBin = 0;
  for( ; iBin < signal_monojet->GetNbinsX(); iBin++){
    signal->SetBinContent(iBin+1,signal_monojet->GetBinContent(iBin+1));
    signal->SetBinError(iBin+1,signal_monojet->GetBinError(iBin+1));
  }	
  for(int iBinX = 0; iBinX < signal_monoV->GetNbinsX()+1; iBinX++, iBin++){
    signal->SetBinContent(iBin+1,signal_monoV->GetBinContent(iBinX+1));
    signal->SetBinError(iBin+1,signal_monoV->GetBinError(iBinX+1));
  }	
  
  if(debug){
    cout<<"monojet signal "<<endl;
    for(int iBinX = 0; iBinX < signal_monojet->GetNbinsX(); iBinX++){
      cout<<"Bin content "<<iBinX+1<<" "<<signal_monojet->GetBinContent(iBinX+1)<<endl;
    }
    cout<<"monoV signal "<<endl;
    for(int iBinX = 0; iBinX < signal_monoV->GetNbinsX(); iBinX++){
      cout<<"Bin content "<<iBinX+1<<" "<<signal_monoV->GetBinContent(iBinX+1)<<endl;
    }
    cout<<"total signal "<<endl;
    for(int iBinX = 0; iBinX < signal->GetNbinsX(); iBinX++){
      cout<<"Bin content "<<iBinX+1<<" "<<signal->GetBinContent(iBinX+1)<<endl;
    }
  }
  return signal;
}

// reduce the big covariance matrix
TH2F* importCorrelationMatrix (TH2F* inputCovariance, const string & monojetLabel, const string & monovLabel){

  int nBinsMonojet = 0;
  int nBinsMonoV   = 0;

  for(int iBinX = 0; iBinX < inputCovariance->GetNbinsX(); iBinX++){
    if(TString(inputCovariance->GetXaxis()->GetBinLabel(iBinX+1)).Contains(monojetLabel.c_str()))
      nBinsMonojet++;
    if(TString(inputCovariance->GetXaxis()->GetBinLabel(iBinX+1)).Contains(monovLabel.c_str()))
      nBinsMonoV++;    
  }

  if(nBinsMonojet == 0 or nBinsMonoV == 0)
    return NULL;

  TH2F* outputCovariance = new TH2F("reducedCovariance","",nBinsMonojet+nBinsMonoV,0,nBinsMonojet+nBinsMonoV+1,nBinsMonojet+nBinsMonoV,0,nBinsMonojet+nBinsMonoV+1);
  outputCovariance->Sumw2();
  int jBinX = 0;
  for(int iBinX = 0; iBinX < inputCovariance->GetNbinsX(); iBinX++){
    int jBinY = 0;
    if(TString(inputCovariance->GetXaxis()->GetBinLabel(iBinX+1)).Contains(monojetLabel.c_str()) or TString(inputCovariance->GetXaxis()->GetBinLabel(iBinX+1)).Contains(monovLabel.c_str())){
      for(int iBinY = 0; iBinY < inputCovariance->GetNbinsY(); iBinY++){
	if(TString(inputCovariance->GetYaxis()->GetBinLabel(iBinY+1)).Contains(monojetLabel.c_str()) or TString(inputCovariance->GetYaxis()->GetBinLabel(iBinY+1)).Contains(monovLabel.c_str())){
	  if(debug)
	    cout<<" X-axis label "<<inputCovariance->GetXaxis()->GetBinLabel(iBinX+1)<<" Y-axis label "<<inputCovariance->GetYaxis()->GetBinLabel(iBinY+1)<<" binx "<<jBinX+1<<" biny "<<jBinY<<" value "<<inputCovariance->GetBinContent(iBinX+1,iBinY+1)<<endl;
	  outputCovariance->SetBinContent(jBinX+1,jBinY+1,inputCovariance->GetBinContent(iBinX+1,iBinY+1));
	  jBinY++;
	}
      }
      jBinX++;
    }
  }

  return outputCovariance;

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
    mu = new RooRealVar("mu","",1.,-1e4,1e4);

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
      if(not skipCorrelations)
	covariance[iBin][jBin] = correlation->GetBinContent(iBin+1,jBin+1);
      else{
	if(jBin == iBin)
	  covariance[iBin][jBin] = correlation->GetBinContent(iBin+1,jBin+1);
	else
	  covariance[iBin][jBin] = 0;
      }
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
    double expectation = double(int(((RooRealVar*) binSum->at(iBin))->getVal()));
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
