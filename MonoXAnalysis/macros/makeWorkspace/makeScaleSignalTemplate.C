#include "../makeTemplates/histoUtils.h"

void makeScaleSignalTemplate(string workspaceFile, Category category, string outputName, float lumiScale = 13.2){
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* inputFile = TFile::Open(workspaceFile.c_str(),"READ");
  RooWorkspace* w = NULL;
  RooWorkspace* outw = NULL;

  TFile* output  = new TFile(outputName.c_str(),"RECREATE");
  
  if(category == Category::monojet){
    w = (RooWorkspace*) inputFile->Get("SR_MJ");
    outw = new RooWorkspace("SR_MJ");
  }
  else if(category == Category::monoV){
    w = (RooWorkspace*) inputFile->Get("SR_MV");
    outw = new RooWorkspace("SR_MV");
  }
  std::list<RooAbsData*> dataset = w->allData();
  
  RooRealVar* var = NULL;
  if(category == Category::monojet)
    var = (RooRealVar*) w->var("met_MJ");
  else if(category == Category::monoV)
    var = (RooRealVar*) w->var("met_MV");

  for (std::list<RooAbsData *>::iterator data=dataset.begin(); data != dataset.end(); ++data){
    if(category == Category::monoV or category == Category::monojet){
      cout<<(*data)->GetName()<<endl;
      if(TString((*data)->GetName()).Contains("ggH") or TString((*data)->GetName()).Contains("WH") or TString((*data)->GetName()).Contains("qqH") or TString((*data)->GetName()).Contains("ZH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	histo->Scale(lumiScale);
	RooDataHist hist(Form("%s",(*data)->GetName()),"",RooArgList(*var), histo);
	outw->import(hist);
      }
      else{
	outw->import(**data);
      }    
    }
  }

  outw->Write();
  output->Close();

}
