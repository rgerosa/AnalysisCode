void makeYieldTableFromMLFit(string workspaceFileName, string inputFileName, string category, bool isZeynep){

  if(category != "monojet" and category != "monoV"){
    cerr<<"Problem in category definition check --> exit "<<endl;
    return;
  }

  TFile* inputFile = TFile::Open(inputFileName.c_str());
  TFile* workspaceFile = TFile::Open(workspaceFileName.c_str());
  RooWorkspace* ws = NULL;

  if(not isZeynep and category == "monojet")
    ws = (RooWorkspace*) workspaceFile->Get("SR_MJ");
  else if(not isZeynep and category == "monoV")
    ws = (RooWorkspace*) workspaceFile->Get("SR_MV");
  else if(isZeynep and category == "monojet")
    ws = (RooWorkspace*) workspaceFile->Get("combinedws");
  else if(isZeynep and category == "monoV")
    ws = (RooWorkspace*) workspaceFile->Get("combinedws");

  TH1* datahist = NULL;
  if(not isZeynep and category == "monojet"){
    RooDataHist* datatemp = (RooDataHist*) ws->data("data_obs_SR_MJ");
    RooRealVar* met_obs  = (RooRealVar*) ws->var("met_MJ");
    datahist = datatemp->createHistogram("data_obs",*met_obs);
  }
  else if(not isZeynep and category == "monoV"){
    RooDataHist* datatemp = (RooDataHist*) ws->data("data_obs_SR_MV");
    RooRealVar* met_obs  = (RooRealVar*) ws->var("met_MV");
    datahist = datatemp->createHistogram("data_obs",*met_obs);
  }
  else if(isZeynep and category == "monojet"){
    RooDataHist* datatemp = (RooDataHist*) ws->data("monojet_signal_data");
    RooRealVar* met_obs  = (RooRealVar*) ws->var("met_monojet");
    datahist = datatemp->createHistogram("data_obs",*met_obs);
  }
  else if(not isZeynep and category == "monoV"){
    RooDataHist* datatemp = (RooDataHist*) ws->data("monov_signal_data");
    RooRealVar* met_obs  = (RooRealVar*) ws->var("met_monov");
    datahist = datatemp->createHistogram("data_obs",*met_obs);
  }

  TH1* zvvhist = NULL;
  TH1* wjethist = NULL;
  TH1* tophist = NULL;
  TH1* dibosonhist = NULL;
  TH1* qcdhist = NULL;
  TH1* dyhist = NULL;
  TH1* gammahist = NULL;
  TH1* totalhist = NULL;

  if(not isZeynep and category == "monojet"){
    zvvhist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/Znunu");
    wjethist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/WJets");
    dyhist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/ZJets");
    tophist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/Top");
    dibosonhist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/Dibosons");
    qcdhist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/QCD");
    gammahist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/GJets");
    totalhist = (TH1*)inputFile->Get("shapes_fit_b/ch1_ch1/total_background");
  }
  else if(not isZeynep and category == "monoV"){
    zvvhist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/Znunu");
    wjethist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/WJets");
    dyhist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/ZJets");
    tophist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/Top");
    dibosonhist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/Dibosons");
    qcdhist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/QCD");
    gammahist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/GJets");
    totalhist = (TH1*)inputFile->Get("shapes_fit_b/ch2_ch1/total_background");
  }
  else if(isZeynep and category == "monojet"){
    zvvhist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/zjets");
    wjethist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/wjets");
    dyhist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/zll");
    tophist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/top");
    dibosonhist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/diboson");
    qcdhist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/qcd");
    gammahist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/gjets");
    totalhist = (TH1*)inputFile->Get("shapes_fit_b/monojet_signal/total_background");
  }
  else if(isZeynep and category == "monoV"){
    zvvhist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/zjets");
    wjethist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/wjets");
    dyhist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/zll");
    tophist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/top");
    dibosonhist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/diboson");
    qcdhist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/qcd");
    gammahist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/gjets");
    totalhist = (TH1*)inputFile->Get("shapes_fit_b/monov_signal/total_background");
  }

  TH1* other = (TH1*) qcdhist->Clone("other");
  other->Add(gammahist);
  other->Add(dyhist);

  // make the table
  ofstream outputfile;
  outputfile.open(Form("yield_%s.txt",category.c_str()));
  outputfile<<"$E_{T}^{miss}$ (GeV) & Observed & $Z \\rightarrow \\nu\\nu$+jets & $W \\rightarrow \\ell\\nu$+jets & Top & Dibosons & Other & Total Bkg. \\\\"<<endl; 
  for(int ibin = 0; ibin < totalhist->GetNbinsX(); ibin++){
    outputfile<<Form("%d-%d",int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)),int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2)))<<" & ";
    outputfile<<Form("%d",int(datahist->GetBinContent(ibin+1)))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",zvvhist->GetBinContent(ibin+1)*zvvhist->GetBinWidth(ibin+1),zvvhist->GetBinError(ibin+1)*zvvhist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",wjethist->GetBinContent(ibin+1)*wjethist->GetBinWidth(ibin+1),wjethist->GetBinError(ibin+1)*wjethist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",tophist->GetBinContent(ibin+1)*tophist->GetBinWidth(ibin+1),tophist->GetBinError(ibin+1)*tophist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",dibosonhist->GetBinContent(ibin+1)*dibosonhist->GetBinWidth(ibin+1),dibosonhist->GetBinError(ibin+1)*dibosonhist->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",other->GetBinContent(ibin+1)*other->GetBinWidth(ibin+1),other->GetBinError(ibin+1)*other->GetBinWidth(ibin+1))<<" & ";
    outputfile<<Form("%f $\\pm$ %f",totalhist->GetBinContent(ibin+1)*totalhist->GetBinWidth(ibin+1),totalhist->GetBinError(ibin+1)*totalhist->GetBinWidth(ibin+1))<<" \\\\ ";
    outputfile<<"\n";
  }
  outputfile.close();
}

