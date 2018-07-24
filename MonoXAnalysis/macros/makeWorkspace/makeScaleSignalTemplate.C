#include "../makeTemplates/histoUtils.h"

void makeScaleSignalTemplate(string workspaceFile, Category category, string outputName, float lumiScale = 1, bool applyCorrectionFromFile = false, bool isSignal = false){
  
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
    if(category == Category::monojet and not isSignal){
      inputCorrection_1 = TFile::Open("templates_monojet_noVBF/templates_bkg_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monojet_noVBF/templates_bkg.root","READ");
    }
    else if(category == Category::monoV and not isSignal){
      inputCorrection_1 = TFile::Open("templates_monoV_noVBF/templates_bkg_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monoV_noVBF/templates_bkg.root","READ");
    }
    else if(category == Category::monojet and isSignal){
      inputCorrection_1 = TFile::Open("templates_monojet_noVBF/templates_signal_higgs_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monojet_noVBF/templates_signal_higgs.root","READ");
    }
    else if(category == Category::monoV and isSignal){
      inputCorrection_1 = TFile::Open("templates_monoV_noVBF/templates_signal_higgs_noVBF.root","READ");
      inputCorrection_2 = TFile::Open("templates_monoV_noVBF/templates_signal_higgs.root","READ");
    }
  }
  output->cd();
  for (std::list<RooAbsData *>::iterator data=dataset.begin(); data != dataset.end(); ++data){
    if(category == Category::monoV or category == Category::monojet){
      if(TString((*data)->GetName()).Contains("ggH") or TString((*data)->GetName()).Contains("WH") or TString((*data)->GetName()).Contains("qqH") or TString((*data)->GetName()).Contains("ZH")){
	TH1D* histo = (TH1D*) (*data)->createHistogram(Form("%s_temp",(*data)->GetName()),*var);
	if(applyCorrectionFromFile){
	  TH1F* histo_1 = NULL;
	  TH1F* histo_2 = NULL;
	  TH1F* correction = NULL;
	  if(not isSignal){
	    histo_1 = (TH1F*) inputCorrection_1->Get("SR/zinvhist_met");
	    histo_2 = (TH1F*) inputCorrection_2->Get("SR/zinvhist_met");	    
	    correction = (TH1F*) histo_1->Clone("correction_bkg");
	    correction->Divide(histo_2);
	  }
	  else{
	    if(TString((*data)->GetName()).Contains("ggH")){
	      histo_1 = (TH1F*) inputCorrection_1->Get("ggH/ggHhist_125_met");
	      histo_2 = (TH1F*) inputCorrection_2->Get("ggH/ggHhist_125_met");
	      correction = (TH1F*) histo_1->Clone("correction_ggH");
	      correction->Divide(histo_2);   
	    }
	    else if(TString((*data)->GetName()).Contains("qqH")){
	      histo_1 = (TH1F*) inputCorrection_1->Get("vbfH/vbfHhist_125_met");
	      histo_2 = (TH1F*) inputCorrection_2->Get("vbfH/vbfHhist_125_met");
	      correction = (TH1F*) histo_1->Clone("correction_qqH");
	      correction->Divide(histo_2);   
	    }
	    else if(TString((*data)->GetName()).Contains("WH")){
	      histo_1 = (TH1F*) inputCorrection_1->Get("wH/wHhist_125_met");
	      histo_2 = (TH1F*) inputCorrection_2->Get("wH/wHhist_125_met");
	      correction = (TH1F*) histo_1->Clone("correction_WH");
	      correction->Divide(histo_2);   
	    }
	    else if(TString((*data)->GetName()).Contains("ZH") and not TString((*data)->GetName()).Contains("ggZH")){
	      histo_1 = (TH1F*) inputCorrection_1->Get("zH/zHhist_125_met");
	      histo_2 = (TH1F*) inputCorrection_2->Get("zH/zHhist_125_met");
	      correction = (TH1F*) histo_1->Clone("correction_ZH");
	      correction->Divide(histo_2);   
	    }
	    else if(TString((*data)->GetName()).Contains("ggZH")){
	      histo_1 = (TH1F*) inputCorrection_1->Get("ggZH/ggZHhist_125_met");
	      histo_2 = (TH1F*) inputCorrection_2->Get("ggZH/ggZHhist_125_met");
	      correction = (TH1F*) histo_1->Clone("correction_ggZH");
	      correction->Divide(histo_2);   	      
	    }
	  }	  	  
	  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
	    if(correction){
		histo->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)*correction->GetBinContent(iBin+1));
		histo->SetBinError(iBin+1,histo->GetBinError(iBin+1)*correction->GetBinContent(iBin+1));
	    }
	  }
	}

	if(TString((*data)->GetName()).Contains("qqH"))
	   histo->Scale(lumiScale);
	
	   
	if(TString((*data)->GetName()).Contains("hptUp") or TString((*data)->GetName()).Contains("hptDown")){
	  (*data)->SetName(TString((*data)->GetName()).ReplaceAll("hptUp","ggH_pTUp"));
	  (*data)->SetName(TString((*data)->GetName()).ReplaceAll("hptDown","ggH_pTDown"));
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
