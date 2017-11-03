#include "../makeTemplates/histoUtils.h"

void makeScaleSignalTemplate(string workspaceFile, Category category, string outputName, float lumiScale = 13.2, bool applyCorrectionFromFile = false){
  
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

  // to apply correction to signal yields as ratio/variation on Zvv backgroud --> for correcting Higgs templates after removing VBF evets  
  TFile* inputCorrection_1 = NULL;
  TFile* inputCorrection_2 = NULL;
  TH1F* correction = NULL;

  if(applyCorrectionFromFile){

    if(category == Category::monojet){
      inputCorrection_1 = TFile::Open("templates_monojet/templates_bkg_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monojet/templates_bkg.root","READ");
    }
    else if(category == Category::monoV){
      inputCorrection_1 = TFile::Open("templates_monoV/templates_bkg_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monoV/templates_bkg.root","READ");
    }

    TH1F* zvv_1 = (TH1F*) inputCorrection_1->Get("SR/zinvhist_met");
    TH1F* zvv_2 = (TH1F*) inputCorrection_2->Get("SR/zinvhist_met");
    correction = (TH1F*) zvv_1->Clone("correction");
    correction->Divide(zvv_2);
  }

  output->cd();
  for (std::list<RooAbsData *>::iterator data=dataset.begin(); data != dataset.end(); ++data){
    if(category == Category::monoV or category == Category::monojet){
      cout<<(*data)->GetName()<<endl;
      if(TString((*data)->GetName()).Contains("ggH") or TString((*data)->GetName()).Contains("WH") or TString((*data)->GetName()).Contains("qqH") or TString((*data)->GetName()).Contains("ZH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	if(not applyCorrectionFromFile)
	  histo->Scale(lumiScale);
	else{
	  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
	    histo->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)*correction->GetBinContent(iBin+1));
	    histo->SetBinError(iBin+1,histo->GetBinError(iBin+1)*correction->GetBinContent(iBin+1));
	  }
	}
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
