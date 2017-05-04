#include "../makeTemplates/histoUtils.h"

void applyGGHiggsNNLO(string workspaceFile, Category category, string outputName){
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* inputFile = TFile::Open(workspaceFile.c_str(),"READ");
  RooWorkspace* w = NULL;
  RooWorkspace* outw = NULL;

  ///// --------- 
  TFile* inputGGHPt = TFile::Open("../../data/HiggsPT/ggH_NNLO_phill.root","READ");
  TH1F* MG_NNLO_FT  = (TH1F*) inputGGHPt->Get("MG_NNLO_FT");
  TH1F* Powheg      = (TH1F*) inputGGHPt->Get("Powheg");
  MG_NNLO_FT->Divide(Powheg);

  TFile* output  = new TFile(outputName.c_str(),"RECREATE");
  
  ///// ---------
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
      if(TString((*data)->GetName()).Contains("ggH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
	  histo->SetBinContent(iBin,histo->GetBinContent(iBin)*MG_NNLO_FT->GetBinContent(MG_NNLO_FT->FindBin(histo->GetBinCenter(iBin+1))));	  
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
