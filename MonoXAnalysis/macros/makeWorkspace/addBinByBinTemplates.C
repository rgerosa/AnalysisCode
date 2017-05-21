#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) < 0)
      histo->SetBinContent(iBin+1,0);
  }
}

/// function to create a RooDataHist from TH1F and import it in a workspace                                                                                                                          
void addTemplate(string procname,
                 const RooArgList& varlist,
                 RooWorkspace& ws,
                 TH1F* hist,
                 const bool & isCoutAndCount = false) {

  if(hist == 0 || hist == NULL){
    cerr<<"addTemplate: found an null pointer --> check "<<endl;
    return;
  }

  if(hist->Integral() == 0)
    addDummyBinContent(hist); // avoind empty histograms in the workspace                                                                                                                             

  checkNegativeBin(hist);

  if(not isCoutAndCount){
    RooDataHist rhist((procname).c_str(), "", varlist, hist);
    ws.import(rhist);
  }
  else{
    RooRealVar *var = dynamic_cast<RooRealVar*>(varlist.at(0));
    vector<double> bin;
    bin.push_back(var->getMin());
    bin.push_back(var->getMax());
    TH1F* hist_temp = new TH1F(Form("%s_cutAndCount",hist->GetName()),"",bin.size()-1,&bin[0]);
    hist_temp->SetBinContent(var->getBins(),hist->Integral(hist->FindBin(var->getMin()),hist->FindBin(var->getMax())));
    RooDataHist rhist((procname).c_str(), "", varlist,hist_temp);
    ws.import(rhist);
  }
}

void generateStatTemplate(string procname,
			  string postfix,
                          const RooArgList& varlist,
                          RooWorkspace& ws,
                          TH1* histo,
                          float scaleUncertainty = 1.,
                          const bool & isCutAndCount = false){
  vector<TH1F*> histStatUp;
  vector<TH1F*> histStatDw;

  if(histo == 0 || histo == NULL){
    cerr<<"generateStatTemplate: found an null pointer --> check "<<endl;
    return;
  }
  if(histo->Integral() == 0)
    addDummyBinContent(histo);

  if(not isCutAndCount){
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statUp_tmp",procname.c_str(),postfix.c_str(),iBin+1)));
      histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statDown_tmp",procname.c_str(),postfix.c_str(),iBin+1)));
    }

    for( size_t iHisto =0; iHisto < histStatUp.size(); iHisto++){
      histStatUp.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)+histo->GetBinError(iHisto+1)*scaleUncertainty);
      if(histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty >= 0.)
	histStatDw.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty);
      else
        histStatDw.at(iHisto)->SetBinContent(iHisto+1,0.);
    }
  }
  else{ // for cut and count stat uncertainty                                                                                                                                                        
    histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_statUp_tmp",procname.c_str(),postfix.c_str())));
    histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_statDown_tmp",procname.c_str(),postfix.c_str())));
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histStatUp.at(0)->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)+histo->GetBinError(iBin+1)*scaleUncertainty);
      histStatDw.at(0)->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)-histo->GetBinError(iBin+1)*scaleUncertainty);
    }
  }

  for(size_t iHisto =0; iHisto < histStatUp.size(); iHisto++) {
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping                                                                                              
    if(not isCutAndCount)
      name << procname << "_"<< postfix << "_CMS_bin" << iHisto+1 << "_statUp";
    else
      name << procname << "_"<< postfix << "_CMS_bin_statUp";

    addTemplate(name.str(),varlist,ws,histStatUp[iHisto],isCutAndCount);
  }
  for(size_t iHisto =0; iHisto < histStatDw.size(); iHisto++){
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping                                                                                               
    if(not isCutAndCount)
      name << procname << "_" << postfix << "_CMS_bin" << iHisto+1 << "_statDown";
    else
      name << procname << "_" << postfix << "_CMS_bin_statDown";

    addTemplate(name.str(),varlist,ws,histStatDw[iHisto],isCutAndCount);
  }
}


/// main function
void addBinByBinTemplates(string inputFileName, Category category, string signalName){

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  RooWorkspace* w = (RooWorkspace*) inputFile->Get("combinedws");  

  TString outName (inputFileName.c_str());
  outName.ReplaceAll(".root","");  
  TFile* outputFile = new TFile(outName+"_updated.root","RECREATE");
  RooWorkspace* outw = new RooWorkspace("combinedws");  
  
  // take whole list of histograms in the file
  std::list<RooAbsData*> dataset = w->allData();

  RooRealVar* var = NULL;
  if(category == Category::monojet)
    var = (RooRealVar*) w->var("met_monojet");
  else if(category == Category::monoV)
    var = (RooRealVar*) w->var("met_monov");

  RooArgList varList (*var);

  // name as postfix
  string postfix = signalName;
  if(category == Category::monojet)
    postfix += "_SR_MJ";
  else if(category == Category::monoV)
    postfix += "_SR_MV";
    

  for (std::list<RooAbsData *>::iterator data=dataset.begin(); data != dataset.end(); ++data){
    cout<<"Generate bin-by-bin variations for "<<(*data)->GetName()<<endl;
    TH1F* histo = (TH1F*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
    // fix on the fly problems in the templates
    checkNegativeBin(histo);
    if(histo->Integral() == 0) addDummyBinContent(histo);
    smoothEmptyBins(histo,2); // do some smoothing
    addTemplate(string((*data)->GetName()),varList,*outw,histo,false);
    // produce the bin-by-bin stat variations
    generateStatTemplate(string((*data)->GetName()),postfix,varList,*outw,histo,1,false);
  }
  
  outw->Write();
  outputFile->Close();
   

}
