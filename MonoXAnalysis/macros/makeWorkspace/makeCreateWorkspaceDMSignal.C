#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>
#include <sstream>

#include "TSystem.h"
#include "workspaceUtils.h"

using namespace std;

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


static bool addBinByBinMCUncertainty = true; // add bin-by-bin uncertainties

// function to create workspace, to be run from a release which has the combine package
void makeCreateWorkspaceDMSignal(string   inputName,                        // input template file
				 Category category,                         // analysis category
				 string   observable  = "met",            // observable 1D or 2D
				 string   interaction = "Vector",         // Interaction type
				 string   signalType  = "MonoJ",
				 bool     isLongLived = false
				 ){

  if(category != Category::monojet and category != Category::monoV){
    cout<<"Code only valid for monojet and monoV category --> please check"<<endl;
    return;
  }

  if(signalType != "MonoJ" and signalType != "MonoW" and signalType != "MonoZ"){
    cout<<"Invalid signal type --> please check"<<endl;
    return;
  }

  if(interaction != "Vector" and interaction != "Axial" and interaction != "Scalar" and interaction != "PseudoScalar"){
    cout<<"Interaction type don't recognized --> please check "<<endl;
    return;
  }

  // identifier code for histogram name
  int code = 0;
  if(interaction == "Vector") code = 1;
  else if(interaction == "Axial") code = 2;
  else if(interaction == "Scalar") code = 3;
  else if(interaction == "PseudoScalar") code = 4;

  string cat_string;
  if(category == Category::monojet) cat_string =  "catmonojet";
  else if(category == Category::monoV) cat_string =  "catmonov";

  // basic loads
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooMsgService::instance().setSilentMode(kTRUE); 
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  // to load all the variables information
  initializeBinning();

  // for templates and sys naming
  string suffix;
  if(category == Category::monojet)    suffix = "monojet";
  else if(category == Category::monoV) suffix = "monov";
    
  // create the output workspac
  cout<<"Create output file ..."<<endl;

  string outputName;
  if(not isLongLived)
    outputName = signalType+"_"+to_string(code)+"_"+cat_string+"_13TeV.root"; // output workspace name    
  else
    outputName = signalType+"_"+to_string(code)+"_"+cat_string+"_LL_13TeV.root"; // output workspace name    

  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  outfile->cd();
  

  // Select observable and binning
  cout<<"Load binning and observable ..."<<endl;
  vector<double> bins = selectBinning(observable,category);  
  RooBinning *binning = NULL;
  double xMin = 0, xMax = 0;

  if(not bins.empty()){ // non empty
    xMin = bins.at(0);
    xMax = bins.back();
    binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str()); 
    }
  else{
    
    bins.clear();
    bin2D bin = selectBinning2D(observable,category);    
    if(not bin.binX.empty() and not bin.binY.empty()){ // in case of 2D analysis --> unrolled histo
      xMin = 0.;
      xMax = double((bin.binX.size()-1)*(bin.binY.size()-1));
      for(size_t iBin = xMin; iBin <= xMax ; iBin++)
	bins.push_back(double(iBin));
      binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str()); 
    }      
    else
      cout<<"Binning not implemented for the observable "<<observable<<" --> please define it "<<endl;
  }

  /// Make observable
  RooRealVar* met = new RooRealVar((observable+"_"+suffix).c_str(),"",xMin,xMax);;    
  met->setBinning(*binning);
  RooArgList vars(*met);

  // Templates
  cout<<"Open inputFils ..."<<endl;
  TFile* templatesfile = TFile::Open(inputName.c_str());

  ///////////////////////////////////////
  // -------- SIGNAL REGION  -------- ///
  ///////////////////////////////////////
  cout<<"Make Signal templates ..."<<endl;
  // create a workspace for the signal region
  RooWorkspace wspace_SR(("SR_"+suffix).c_str(),(suffix+"_SR").c_str());

  string postfix;
  if(category == Category::monojet) postfix = signalType+"_SR_MJ";
  else if(category == Category::monoV) postfix = signalType+"_SR_MV";

  // Loop on the input file -->
  TIter next(templatesfile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if(not TString(key->GetClassName()).Contains("TH1")) continue;

    // split the name
    stringstream histoName(key->GetName());
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(histoName, segment,'_')){
      seglist.push_back(segment);
    }

    string signal = seglist.at(0); // this is MonoJ, MonoW or MonoZ indipendently from DM s-channel or LL samples
    string medMass;
    string dmMass;
    string ctau;

    if(signal != signalType){
      cout<<"Different signal type template found in the file --> skip "<<endl;
      continue;
    }

    if(isLongLived)
      medMass = seglist.at(2);
    else
      medMass = seglist.at(2);

    if(isLongLived)
      dmMass = seglist.at(4);
    else
      dmMass = seglist.at(5);

    if(isLongLived)
      ctau = seglist.at(6);

    TH1F* histo = (TH1F*) templatesfile->Get(key->GetName()); // take the object
    
    // format the right digits
    std::stringstream s;
    s << std::setfill('0') << std::setw(4) << atoi(medMass.c_str());
    s << std::setfill('0') << std::setw(4) << atoi(dmMass.c_str());
    if(isLongLived)
      s << std::setfill('0') << std::setw(6) << atoi(ctau.c_str());

    addTemplate(signalType+"_"+to_string(code)+s.str(),vars,wspace_SR,histo,false);                                                                                                      
    if(addBinByBinMCUncertainty)                                                                                                                                                                     
      generateStatTemplate(signalType+"_"+to_string(code)+s.str(),postfix,vars,wspace_SR,histo,1,false);                                                                                          
                               
  }

  outfile->cd();
  wspace_SR.Write();    
  outfile->Close();
  return;  
}
