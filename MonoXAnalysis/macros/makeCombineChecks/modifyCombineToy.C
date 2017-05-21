#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float shift = 0.97;
static bool  applyFirstBinOnly = true;

void modifyCombineToy(string inputCombineFileName, string outputFileName, bool isHiggsInvisible, string channelToModify){

  TFile* inputFile = TFile::Open(inputCombineFileName.c_str());
  RooDataSet* data_obs = (RooDataSet*) inputFile->Get("toys/toy_asimov");
  RooDataHist* data_hist = data_obs->binnedClone();

  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  outputFile->mkdir("toys");
  outputFile->cd("toys");

  RooArgList* observables = new RooArgList(*data_hist->get());
  RooRealVar* obs_monojet = NULL;
  RooRealVar* obs_monov   = NULL;
  RooCategory* channel    = NULL;

  if(not isHiggsInvisible){
    obs_monojet = (RooRealVar*) observables->find("met_monojet");
    obs_monov = (RooRealVar*) observables->find("met_monov");
    channel = (RooCategory*) observables->find("CMS_channel");
  }
  else{
    obs_monojet = (RooRealVar*) observables->find("met_MJ");
    obs_monov = (RooRealVar*) observables->find("met_MV");
    channel = (RooCategory*) observables->find("CMS_channel");
  }

  RooDataHist* data_ch1_ch1 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 0");
  RooDataHist* data_ch1_ch2 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 1");
  RooDataHist* data_ch1_ch3 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 2");
  RooDataHist* data_ch1_ch4 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 3");
  RooDataHist* data_ch1_ch5 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 4");
  RooDataHist* data_ch1_ch6 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monojet),"CMS_channel == 5");

  RooDataHist* data_ch2_ch1 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 6");
  RooDataHist* data_ch2_ch2 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 7");
  RooDataHist* data_ch2_ch3 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 8");
  RooDataHist* data_ch2_ch4 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 9");
  RooDataHist* data_ch2_ch5 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 10");
  RooDataHist* data_ch2_ch6 = (RooDataHist*) data_hist->reduce(RooArgSet(*obs_monov),"CMS_channel == 11");
  
  ///////////////////
  if(channelToModify == "ZCR"){

    TH1F* data_ch1_ch2_hist = (TH1F*) data_ch1_ch2->createHistogram("data_ch1_ch2_hist",*obs_monojet);
    TH1F* data_ch1_ch5_hist = (TH1F*) data_ch1_ch5->createHistogram("data_ch1_ch5_hist",*obs_monojet);

    if(applyFirstBinOnly){
      data_ch1_ch2_hist->SetBinContent(1,data_ch1_ch2_hist->GetBinContent(1)*shift);
      data_ch1_ch5_hist->SetBinContent(1,data_ch1_ch5_hist->GetBinContent(1)*shift);
    }
    else{
      data_ch1_ch2_hist->Scale(shift);
      data_ch1_ch5_hist->Scale(shift);
    }

    RooDataHist* data_ch1_ch2_modified = new RooDataHist("data_ch1_ch2_modified","",*obs_monojet,data_ch1_ch2_hist);
    RooDataHist* data_ch1_ch5_modified = new RooDataHist("data_ch1_ch5_modified","",*obs_monojet,data_ch1_ch5_hist);

    map<string,RooDataHist*> map;
    map["ch1_ch1"] = data_ch1_ch1;
    map["ch1_ch2"] = data_ch1_ch2_modified;
    map["ch1_ch3"] = data_ch1_ch3;
    map["ch1_ch4"] = data_ch1_ch4;
    map["ch1_ch5"] = data_ch1_ch5_modified;
    map["ch1_ch6"] = data_ch1_ch6;
    map["ch2_ch1"] = data_ch2_ch1;
    map["ch2_ch2"] = data_ch2_ch2;
    map["ch2_ch3"] = data_ch2_ch3;
    map["ch2_ch4"] = data_ch2_ch4;
    map["ch2_ch5"] = data_ch2_ch5;
    map["ch2_ch6"] = data_ch2_ch6;    

    RooDataHist* data_obs_modified = new RooDataHist("data_obs","",RooArgList(*obs_monojet,*obs_monov),*channel,map);
    data_obs_modified->Write("toy_asimov");

  }
  ///////////////
  if(channelToModify == "WCR"){

    TH1F* data_ch1_ch3_hist = (TH1F*) data_ch1_ch3->createHistogram("data_ch1_ch3_hist",*obs_monojet);
    TH1F* data_ch1_ch6_hist = (TH1F*) data_ch1_ch6->createHistogram("data_ch1_ch6_hist",*obs_monojet);

    if(applyFirstBinOnly){
      data_ch1_ch3_hist->SetBinContent(1,data_ch1_ch3_hist->GetBinContent(1)*shift);
      data_ch1_ch6_hist->SetBinContent(1,data_ch1_ch6_hist->GetBinContent(1)*shift);
    }
    else{
      data_ch1_ch3_hist->Scale(shift);
      data_ch1_ch6_hist->Scale(shift);
    }

    RooDataHist* data_ch1_ch3_modified = new RooDataHist("data_ch1_ch3_modified","",*obs_monojet,data_ch1_ch3_hist);
    RooDataHist* data_ch1_ch6_modified = new RooDataHist("data_ch1_ch6_modified","",*obs_monojet,data_ch1_ch6_hist);

    map<string,RooDataHist*> map;
    map["ch1_ch1"] = data_ch1_ch1;
    map["ch1_ch2"] = data_ch1_ch2;
    map["ch1_ch3"] = data_ch1_ch3_modified;
    map["ch1_ch4"] = data_ch1_ch4;
    map["ch1_ch5"] = data_ch1_ch5;
    map["ch1_ch6"] = data_ch1_ch6_modified;
    map["ch2_ch1"] = data_ch2_ch1;
    map["ch2_ch2"] = data_ch2_ch2;
    map["ch2_ch3"] = data_ch2_ch3;
    map["ch2_ch4"] = data_ch2_ch4;
    map["ch2_ch5"] = data_ch2_ch5;
    map["ch2_ch6"] = data_ch2_ch6;    

    RooDataHist* data_obs_modified = new RooDataHist("data_obs","",RooArgList(*obs_monojet,*obs_monov),*channel,map);
    data_obs_modified->Write("toy_asimov");

  }

  ///////
  if(channelToModify == "GCR"){

    TH1F* data_ch1_ch4_hist = (TH1F*) data_ch1_ch4->createHistogram("data_ch1_ch4_hist",*obs_monojet);
    if(applyFirstBinOnly){
      data_ch1_ch4_hist->SetBinContent(1,data_ch1_ch4_hist->GetBinContent(1)*shift);
    }
    else{
      data_ch1_ch4_hist->Scale(shift);
    }

    RooDataHist* data_ch1_ch4_modified = new RooDataHist("data_ch1_ch4_modified","",*obs_monojet,data_ch1_ch4_hist);

    map<string,RooDataHist*> map;
    map["ch1_ch1"] = data_ch1_ch1;
    map["ch1_ch2"] = data_ch1_ch2;
    map["ch1_ch3"] = data_ch1_ch3;
    map["ch1_ch4"] = data_ch1_ch4_modified;
    map["ch1_ch5"] = data_ch1_ch5;
    map["ch1_ch6"] = data_ch1_ch6;
    map["ch2_ch1"] = data_ch2_ch1;
    map["ch2_ch2"] = data_ch2_ch2;
    map["ch2_ch3"] = data_ch2_ch3;
    map["ch2_ch4"] = data_ch2_ch4;
    map["ch2_ch5"] = data_ch2_ch5;
    map["ch2_ch6"] = data_ch2_ch6;    

    RooDataHist* data_obs_modified = new RooDataHist("data_obs","",RooArgList(*obs_monojet,*obs_monov),*channel,map);
    data_obs_modified->Write("toy_asimov");

  }
  //////
  if(channelToModify == "SR"){
    TH1F* data_ch1_ch1_hist = (TH1F*) data_ch1_ch1->createHistogram("data_ch1_ch1_hist",*obs_monojet);
    if(applyFirstBinOnly){
      data_ch1_ch1_hist->SetBinContent(1,data_ch1_ch1_hist->GetBinContent(1)*shift);
    }
    else{
      data_ch1_ch1_hist->Scale(shift);
    }

    RooDataHist* data_ch1_ch1_modified = new RooDataHist("data_ch1_ch1_modified","",*obs_monojet,data_ch1_ch1_hist);

    
    map<string,RooDataHist*> map;
    map["ch1_ch1"] = data_ch1_ch1_modified;
    map["ch1_ch2"] = data_ch1_ch2;
    map["ch1_ch3"] = data_ch1_ch3;
    map["ch1_ch4"] = data_ch1_ch4;
    map["ch1_ch5"] = data_ch1_ch5;
    map["ch1_ch6"] = data_ch1_ch6;
    map["ch2_ch1"] = data_ch2_ch1;
    map["ch2_ch2"] = data_ch2_ch2;
    map["ch2_ch3"] = data_ch2_ch3;
    map["ch2_ch4"] = data_ch2_ch4;
    map["ch2_ch5"] = data_ch2_ch5;
    map["ch2_ch6"] = data_ch2_ch6;
    RooDataHist* data_obs_modified = new RooDataHist("data_obs","",RooArgList(*obs_monojet,*obs_monov),*channel,map);
    data_obs_modified->Write("toy_asimov");
    
  }

  outputFile->Close();
  
}
