#include "../../../../macros/makeTemplates/histoUtils.h"

void modify_workspace(string workspaceFile, Category category, float lumiScale = 3.32){

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  float efficiencyMV = 0.9;
  float efficiencyMJ = 0.9;
  
  TFile* inputFile = TFile::Open(workspaceFile.c_str(),"READ");
  RooWorkspace* w = NULL;
  if(category == Category::monojet)
    w = (RooWorkspace*) inputFile->Get("SR_MJ");
  else if(category == Category::monoV)
    w = (RooWorkspace*) inputFile->Get("SR_MV");
  std::list<RooAbsData*> dataset = w->allData();

  RooRealVar* var = NULL;
  if(category == Category::monojet)
    var = (RooRealVar*) w->var("met_MJ");
  else if(category == Category::monoV)
    var = (RooRealVar*) w->var("met_MV");

  for (std::list<RooAbsData *>::iterator data=dataset.begin(); data != dataset.end(); ++data){
    if(category == Category::monoV){
	cout<<(*data)->GetName()<<endl;
      if(TString((*data)->GetName()).Contains("ggH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	histo->Scale(efficiencyMV*1.1*lumiScale);
	RooDataHist hist(Form("%s_scaled",(*data)->GetName()),"",RooArgList(*var), histo);
	w->import(hist);
      }
        
      else if(TString((*data)->GetName()).Contains("WH") or TString((*data)->GetName()).Contains("qqH") or TString((*data)->GetName()).Contains("ZH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	histo->Scale(efficiencyMV*lumiScale);
	RooDataHist hist(Form("%s_scaled",(*data)->GetName()),"",RooArgList(*var), histo);
	w->import(hist);
      }
    }
    else if(category == Category::monojet){
      if(TString((*data)->GetName()).Contains("ggH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	histo->Scale(efficiencyMJ*1.1*lumiScale);
	RooDataHist hist(Form("%s_scaled",(*data)->GetName()),"",RooArgList(*var), histo);
	w->import(hist);
      }
      else if(TString((*data)->GetName()).Contains("WH") or TString((*data)->GetName()).Contains("qqH") or TString((*data)->GetName()).Contains("ZH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	histo->Scale(efficiencyMJ*lumiScale);
	RooDataHist hist(Form("%s_scaled",(*data)->GetName()),"",RooArgList(*var), histo);
	w->import(hist);
      }    
    }
  }  

  TFile* output  = NULL;
  if(category == Category::monojet)
    output = new TFile("output_MJ.root","RECREATE");
  else
    output = new TFile("output_MV.root","RECREATE");
  w->Write();
  output->Close();

}
