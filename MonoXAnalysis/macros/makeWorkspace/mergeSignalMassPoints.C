#include "../CMS_lumi.h"

//inputFileName1 : file to be filled with elements from 2
//inputFileName2 : file to be used to full 1
void mergeSignalMassPoints(string inputFileName1, string inputFileName2, string outputFileName, float scaleHistograms = 1){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  //  RooMsgService::instance().setSilentMode(kTRUE);
  //  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  TFile* inputFile1 = TFile::Open(inputFileName1.c_str(),"READ");
  TFile* inputFile2 = TFile::Open(inputFileName2.c_str(),"READ");

  RooWorkspace* ws1 = (RooWorkspace*) inputFile1->Get("combinedws");
  RooWorkspace* ws2 = (RooWorkspace*) inputFile2->Get("combinedws");

  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  RooWorkspace* ws_out = new RooWorkspace("combinedws","combinedws");

  // take the variable
  RooRealVar* var1 = (RooRealVar*) ws1->var("met_monojet");
  if(var1 == 0 or var1 == NULL)
    var1 = (RooRealVar*) ws1->var("met_monov");

  RooRealVar* var2 = (RooRealVar*) ws2->var("met_monojet");
  if(var2 == 0 or var2 == NULL)
    var2 = (RooRealVar*) ws2->var("met_monov");

  // take list of histograms
  outputFile->cd();
  list<RooAbsData*> list2 = ws2->allData();
  list<RooAbsData*> list1 = ws1->allData();

  cout<<"input file 1 size: "<<list1.size()<<endl;
  cout<<"input file 2 size: "<<list2.size()<<endl;

  int originalGoodMassPoints = 0;
  for(auto dataset1 : list1){
    // skip empty templates
    if(double(dataset1->sumEntries()) <= double(0.0001)){
      continue;
    }
    originalGoodMassPoints++;
    ws_out->import(*dataset1);
  }
  cout<<"originalGoodMassPoints "<<originalGoodMassPoints<<endl;

  //////
  int massPointsFixedNormalization = 0;
  for(auto dataset1 : list1){
    if(double(dataset1->sumEntries()) <= double(0.0001)){
      for(auto dataset2 : list2){
	if(string(dataset2->GetName()) == string(dataset1->GetName())){
	  // make histogram
	  TH1D* histo = (TH1D*) (*dataset2).createHistogram(Form("%s_temp",(*dataset2).GetName()),*var2);
	  histo->Scale(scaleHistograms);
	  RooDataHist hist(Form("%s",(*dataset2).GetName()),"",RooArgList(*var1), histo);
	  massPointsFixedNormalization++;
	  cout<<"dataset2 "<<dataset2->sumEntries()<<" histo "<<histo->Integral()<<" hist "<<hist.sumEntries()<<endl;
	  ws_out->import(hist);
	}
      }
    }
  }

  cout<<"massPointsFixedNormalization "<<massPointsFixedNormalization<<endl;

  int newMassPointsAdded = 0;
  for(auto dataset2 : list2){    
    // search if this dataset is there in 1
    bool isfound = false;
    for(auto dataset1 : list1){      
      if(string(dataset2->GetName()) == string(dataset1->GetName())){
	isfound = true;
	break;
      }
    }
    
    if(isfound == true) continue;
    // make histogram
    TH1D* histo = (TH1D*) (*dataset2).createHistogram(Form("%s_temp",(*dataset2).GetName()),*var1);
    histo->Scale(scaleHistograms);
    RooDataHist hist(Form("%s",(*dataset2).GetName()),"",RooArgList(*var1), histo);
    newMassPointsAdded++;
    ws_out->import(hist);    
    // search if the dataset is there in the other file
  }

  cout<<"newMassPointsAdded "<<newMassPointsAdded<<endl;

  cout<<"Output file: number of templates "<<ws_out->allData().size()<<endl;
  
  outputFile->cd();
  ws_out->Write();
  outputFile->Close();
  

} 
