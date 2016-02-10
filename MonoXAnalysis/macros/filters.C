#include <vector>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"

void filters(){}

// function tht takes in input a tree and counts the total of negative and positive signs
double sumwgt(TTree* tree) {

  if(tree == NULL or tree == 0)
    return 1.;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sumwgt function --> loop on gentree"<<std::endl;

  TTreeReader treeReader(tree);
  TTreeReaderValue<double> wgtsign (treeReader,"wgtsign");

  double weightsum = 0.;
  float nEvents = 0;
  while(treeReader.Next()){    
    if(int(nEvents) %100000 == 0){
      std::cout<<"Events "<<nEvents/tree->GetEntries()*100<<"%"<<std::endl;
    }
    weightsum += (*wgtsign);
    nEvents++;
  }

  std::cout<<"sumwgt function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;
  
  return weightsum;
}

//function to apply pileup re-weight, i.e. return an histogram given by the ratio between pileup in data and mc
TH1D* pileupwgt(TTree* tree, std::string scenario = ""){

  if(tree == NULL or tree == 0)
    return 0;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"pileupwgt function --> start"<<std::endl;
  
  // PU data  
  TH1D* histoPUData;
  TFile* fileInputData;
  
  if(scenario == "silver"){
    fileInputData = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/JSON_PILEUP/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.root");
    histoPUData   = (TH1D*) fileInputData->Get("pileup");
  }
  else{
    fileInputData = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/JSON_PILEUP/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.root");
    histoPUData   = (TH1D*) fileInputData->Get("pileup");
  }
  histoPUData->SetName("histoPUData");
  // scale to unit
  histoPUData->Scale(1./histoPUData->Integral());
  
  // PU MC
  TH1D* histoPUMC = (TH1D*) histoPUData->Clone("histoPUMC");
  // fill an histogram with the same binning and properties of data one
  tree->Draw("putrue >> histoPUMC","","goff");
  // scale to unit
  histoPUMC->Scale(1./histoPUMC->Integral());
 
  // PU ratio
  TH1D* puRatio = (TH1D*) histoPUData->Clone("histoRatio");
  puRatio->Divide(histoPUMC);

  std::cout<<"pileupwgt function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  return puRatio;

}

//function to return the pileup weights
void btagWeights(TTree* tree, TH2F* eff_b, TH2F* eff_c, TH2F* eff_ucsdg){

  std::cout<<"###################################"<<std::endl;
  std::cout<<"btagWeights function --> start"<<std::endl;

  // smooth input histograms to deal with low statistics
  eff_b->Smooth();
  eff_c->Smooth();
  eff_ucsdg->Smooth();


  // make graph 2D for interpolation for central values
  TGraph2D* graph_b = new TGraph2D(eff_b);
  TGraph2D* graph_c = new TGraph2D(eff_c);
  TGraph2D* graph_ucsdg = new TGraph2D(eff_ucsdg);

  double jetPtMin  = 20;
  double jetEtaMax = 2.4;
  double btagWP    = 0.89; // medium working point
  
  // decleare reader
  TTreeReader myReader(tree);
  TTreeReaderValue<std::vector<double> > centraljetpt         (myReader,"centraljetpt");
  TTreeReaderValue<std::vector<double> > centraljeteta        (myReader,"centraljeteta");
  TTreeReaderValue<std::vector<double> > centraljetHFlav      (myReader,"centraljetHFlav");
  TTreeReaderValue<std::vector<double> > centraljetbtag       (myReader,"centraljetbtag");
  TTreeReaderValue<std::vector<double> > centraljetBtagSF     (myReader,"centraljetBtagSF");
  TTreeReaderValue<std::vector<double> > centraljetBtagSFUp   (myReader,"centraljetBtagSFUp");
  TTreeReaderValue<std::vector<double> > centraljetBtagSFDown (myReader,"centraljetBtagSFDown");


  // add a b-tag weight branch for medium wp
  double wgtbtag, wgtbtagUp, wgtbtagDown;
  TBranch* bwgtbtag     = tree->Branch("wgtbtag", &wgtbtag, "wgtbtag/D");  
  TBranch* bwgtbtagUp   = tree->Branch("wgtbtagUp", &wgtbtagUp, "wgtbtagUp/D");  
  TBranch* bwgtbtagDown = tree->Branch("wgtbtagDown", &wgtbtagDown, "wgtbtagDown/D");  

  int nEvents = 0;

  while(myReader.Next()){

    // recipe for b-tag weight on https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
    vector<double> PMC;
    vector<double> PDATA ;
    vector<double> PDATAErrUp ;
    vector<double> PDATAErrDw ;
    double efficiency = 1.;

    if(nEvents %100000 == 0){
      std::cout<<"Events "<<float(nEvents)/tree->GetEntries()*100<<"%"<<std::endl;
    }      

    // take only events in the acceptance region for b-tagging   
    for(size_t iJet = 0; iJet < centraljetpt->size(); iJet++){            
      if(centraljetpt->at(iJet) > jetPtMin and fabs(centraljeteta->at(iJet)) < jetEtaMax){
	// get efficiency
	if(centraljetHFlav->at(iJet) == 5)
	  efficiency     = graph_b->Interpolate(centraljetpt->at(iJet) ,fabs(centraljeteta->at(iJet)));
	else if(centraljetHFlav->at(iJet) == 4)
	  efficiency     = graph_c->Interpolate(centraljetpt->at(iJet),fabs(centraljeteta->at(iJet)));
	else
	  efficiency     = graph_ucsdg->Interpolate(centraljetpt->at(iJet),fabs(centraljeteta->at(iJet)));

	//checks
	if(efficiency > 1) efficiency = 1;
	if(efficiency < 0) efficiency = 0;
	
	// check b-tag value
	if(centraljetbtag->at(iJet) > btagWP){
	  PMC.push_back(efficiency);
	  PDATA.push_back(efficiency*(centraljetBtagSF->at(iJet)));		  
	  PDATAErrUp.push_back(efficiency*(centraljetBtagSFUp->at(iJet)));
	  PDATAErrDw.push_back(efficiency*(centraljetBtagSFDown->at(iJet)));
	}
	else{
	  PMC.push_back((1-efficiency));
	  PDATA.push_back((1-efficiency*(centraljetBtagSF->at(iJet))));
	  PDATAErrUp.push_back(1-efficiency*(centraljetBtagSFUp->at(iJet)));
	  PDATAErrDw.push_back(1-efficiency*(centraljetBtagSFDown->at(iJet)));
	}	
      }
    }

    // Fill branch
    double Num = 1.;
    double Den = 1.;

    for(size_t iJet = 0; iJet < PMC.size(); iJet++){
      Num *= PDATA.at(iJet);
      Den *= PMC.at(iJet);
    }

    if(Den != 0)
      wgtbtag = Num/Den;
    else
      wgtbtag = 1.;

    vector<double> wbtagErr;
    double NumMax = 1.;
    double NumMin = 1.;
    for(size_t iErr = 0; iErr < PDATAErrUp.size(); iErr++){
      if(PDATAErrUp.at(iErr) > PDATAErrDw.at(iErr)){
	NumMax *= PDATAErrUp.at(iErr);
	NumMin *= PDATAErrDw.at(iErr);
      }
      else{
	NumMax *= PDATAErrDw.at(iErr);
	NumMin *= PDATAErrUp.at(iErr);	
      }
    }

    if(Den != 0){
      wgtbtagUp   = NumMax/Den;
      wgtbtagDown = NumMin/Den;
    }
    else{
      wgtbtagUp   = wgtbtag;
      wgtbtagDown = wgtbtag;
    }

    bwgtbtag->Fill();
    bwgtbtagUp->Fill();
    bwgtbtagDown->Fill();
    
    nEvents++;
  }

  std::cout<<"btagWeights function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;

}


// function that take as input a ROOT file
void sigfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sigfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "sigtree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree = NULL;
  TH1D*  puRatio = NULL;
  double wgtsum;
  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    // calculate sum of weights
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  // accept signal region event if loose muons, electrons, taus and phototns == 0 and triggered by hltmet90
  const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && hltmet90 > 0";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree"); 
  treedir->cd();
  // copy the tree applying a selection
  std::cout<<"sigfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"sigfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  
  // add a weight sum branch  
  TBranch* bwgtsum, *bwgtpileup;
  double wgtpileup = 1;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");      
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");
    
    std::cout<<"sigfilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) 
	std::cout<<"Events "<<float(nEvents)/outtree->GetEntries()*100<<"%"<<std::endl;
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();    
      nEvents++;
    }
  }

  
  if(applyBTagWeights and isMC){

    std::cout<<"sigfilter --> apply btag-weight"<<std::endl; 
    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);

    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }
  
  // write tree in the file
  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;
  
}

// function to apply Zmumu selections
void zmmfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"zmmfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zmmtree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree =  NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }
    
  // Selections 2 loose muons, 0 loose ele, taus and photons, m(mumu) = m(z) [60,120], mupt > 20, one of the two tight ... no trigger requirement
  const char* cut = "nmuons == 2 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && zmass > 60 && zmass < 120 && mu1pt > 20 && (mu1id >= 1 || mu2id >= 1)";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"zmmfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"zmmfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"zmmfilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){    
      if(nEvents %100000 == 0) 
	std::cout<<"Events "<<float(nEvents)/outtree->GetEntries()*100<<"%"<<std::endl;
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();    
      nEvents++;
    }
  }

  if(applyBTagWeights and isMC){
    
    std::cout<<"zmmfilter --> apply btag-weight"<<std::endl;
 
    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;

}

// function to apply Zee selections
void zeefilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"zeefilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zeetree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");

  TTree* intree =  NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }
  
  const char* cut = "nmuons == 0 && nelectrons == 2 && ntaus == 0 && nphotons == 0 && zeemass > 60 && zeemass < 120 && el1pt > 40 && (el1id >= 1 || el2id >= 1)";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"zeefilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"zeefilter --> outtree events "<<outtree->GetEntries()<<std::endl;


  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"zeefilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) 
	std::cout<<"Events "<<float(nEvents)/outtree->GetEntries()*100<<"%"<<std::endl;
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();    
      nEvents++;
    }
  }

  if(applyBTagWeights and isMC){
    
    std::cout<<"zeefilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;

}

// function to apply Wmunu selections
void wmnfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"wmnfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wmntree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree =  NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  const char* cut = "nmuons == 1 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && mu1pt > 20 && mu1id >= 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"wmnfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"wmnfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);    
    TTreeReaderValue<int> putrue(myReader,"putrue");
    std::cout<<"wmnfilter --> apply sumwgt and puweight"<<std::endl;    
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) 
	std::cout<<"Events "<<float(nEvents)/outtree->GetEntries()*100<<"%"<<std::endl;
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();          
      nEvents++;      
    }    
  }

  if(applyBTagWeights and isMC){

    std::cout<<"wmnfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);

  }

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;

}

// function to apply Wenu selections
void wenfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"wenfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wentree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree = NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  const char* cut = "nmuons == 0 && nelectrons == 1 && ntaus == 0 && nphotons == 0 && el1pt > 40 && el1id >= 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"wenfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"wenfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"wenfilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) std::cout<<"Event: "<<float(nEvents)/outtree->GetEntries()<<std::flush;
	bwgtsum->Fill();
	wgtpileup = puRatio->GetBinContent(*putrue);
	bwgtpileup->Fill();    
	nEvents++;
    }
  }

  if(applyBTagWeights and isMC){
    
    std::cout<<"wenfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;
}

// function to apply photon+jets selections
void gamfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"gamfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "gamtree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree = NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  // medium id + pt + veto
  const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 1 && phpt > 175 && phidm == 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"gamfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"gamfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"gamfilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) std::cout<<"Event: "<<float(nEvents)/outtree->GetEntries()<<std::flush;
	bwgtsum->Fill();
	wgtpileup = puRatio->GetBinContent(*putrue);
	bwgtpileup->Fill();    
	nEvents++;
    }
  }

  if(applyBTagWeights and isMC){
    
    std::cout<<"gamfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }  

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;

}



// function to apply photon+jets selections
void topfilter( std::string inputFileName,  std::string outputFileName, bool isMC, bool applyBTagWeights, bool storeGenTree = false) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"topfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "toptree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  TTree* intree = NULL;
  double wgtsum;
  TH1D*  puRatio = NULL;
  double wgtpileup = 1;

  if(isMC){
    intree = (TTree*)infile->Get("gentree/gentree");
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  // medium id + pt + veto
  const char* cut = "(nmuons > 0 || nelectrons > 0) && ntaus == 0 && nphotons == 0 && nbjets > 0 && (el1id >= 1 || mu1id >=1)";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"topfilter --> apply signal region preselection"<<std::endl;
  TTree* outtree = frtree->CopyTree(cut);
  std::cout<<"topfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum;
  TBranch* bwgtpileup ;
  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"topfilter --> apply sumwgt and puweight"<<std::endl;
    int nEvents=0;
    while(myReader.Next()){
      if(nEvents %100000 == 0) std::cout<<"Event: "<<float(nEvents)/outtree->GetEntries()<<std::flush;
	bwgtsum->Fill();
	wgtpileup = puRatio->GetBinContent(*putrue);
	bwgtpileup->Fill();    
	nEvents++;
    }
  }

  if(applyBTagWeights and isMC){
    
    std::cout<<"topfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b, *eff_Num_c, *eff_Num_ucsdg;
    TH2F*  eff_Denom_b, *eff_Denom_c, *eff_Denom_ucsdg;
    eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
    eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
    eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
    // take denominators
    eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
    eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
    eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);
    
    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);
  }  

  outfile->cd();
  if(storeGenTree)
    intree->CloneTree()->Write();
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;

}

