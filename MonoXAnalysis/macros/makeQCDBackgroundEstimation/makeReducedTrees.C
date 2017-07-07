#include <vector>
#include <string>
#include <fstream>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"

void makeReducedTrees(){}

// function taking a path and a flag for EOS cern and returns a list of files as vector of strings
vector<string> fileListForChain (const std::string path, bool isEOS){

  vector<string> outputFileList;
  if(isEOS)
    system(("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+path+" -name \"*.root\" | grep -v failed > file.temp").c_str());
  else
    system(("find "+path+" -name \"*.root\" | grep -v failed > file.temp").c_str());

  
  ifstream infile;
  string line;
  infile.open("file.temp",ifstream::in);   
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      if(isEOS)
	outputFileList.push_back("root://eoscms.cern.ch//"+line);
      else
	outputFileList.push_back(line);
    }
  }
  infile.close();
  system("rm file.temp");

  return outputFileList;
}

// sum xsec values inside gen-tree branch                                                                                                                                                              
double sumxsec(TChain* tree){

  if(tree == NULL or tree == 0)
    return 1.;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sumxsec function --> loop on gentree: nEvents = "<<tree->GetEntries()<<std::endl;
  std::cout<<"###################################"<<std::endl;

  TTreeReader treeReader(tree);
  TTreeReaderValue<float> wgtsign (treeReader,"lheXSEC");

  double xsecsum = 0.;
  double nEvents   = 0.;

  while(treeReader.Next()){
    if(int(nEvents) %100000 == 0){
      std::cout.flush();
      std::cout<<"\r"<<"Events "<<nEvents/tree->GetEntries()*100<<"%";
    }
    xsecsum += (*wgtsign);
    nEvents++;
  }
  std::cout<<std::endl;
  std::cout<<"###################################"<<std::endl;
  std::cout<<"sumxsec function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  return xsecsum;

}

// function tht takes in input a tree and counts the total of negative and positive signs
double sumwgt(TChain* tree) {

  if(tree == NULL or tree == 0)
    return 1.;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sumwgt function --> loop on gentree: nEntries "<<tree->GetEntries()<<std::endl;
  std::cout<<"###################################"<<std::endl;

  TTreeReader treeReader(tree);
  TTreeReaderValue<float> wgt (treeReader,"wgt");
  double weightsum = 0.;
  double nEvents   = 0.;

  while(treeReader.Next()){    
    if(int(nEvents) %100000 == 0){
      std::cout.flush();
      std::cout<<"\r "<<"Events "<<nEvents/tree->GetEntries()*100<<"%";
    }
    weightsum += (*wgt);
    nEvents++;
  }

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sumwgt function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;
 
  return weightsum;
  
}

//function to apply pileup re-weight, i.e. return an histogram given by the ratio between pileup in data and mc
TH1D* pileupwgt(TChain* tree, std::string scenario = ""){

  if(tree == NULL or tree == 0)
    return 0;

  std::cout<<"###################################"<<std::endl;
  std::cout<<"pileupwgt function --> start"<<std::endl;
  std::cout<<"###################################"<<std::endl;
  
  // PU data  
  TH1D*  histoPUData;
  TFile* fileInputData;
  
  fileInputData = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/JSON_PILEUP/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root");
  histoPUData   = (TH1D*) fileInputData->Get("pileup");

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

  std::cout<<"###################################"<<std::endl;
  std::cout<<"pileupwgt function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  return puRatio;

}

//function to return the pileup weights
void btagWeights(TTree* tree, TH2F* eff_b, TH2F* eff_c, TH2F* eff_ucsdg){

  std::cout<<"###################################"<<std::endl;
  std::cout<<"btagWeights function --> start"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  // make graph 2D for interpolation for central values
  TGraph2D* graph_b = new TGraph2D(eff_b);
  TGraph2D* graph_c = new TGraph2D(eff_c);
  TGraph2D* graph_ucsdg = new TGraph2D(eff_ucsdg);

  double jetPtMin  = 20;
  double jetEtaMax = 2.4;
  double btagCSVMediumWP = 0.8484; // medium working point
  double btagMVAMediumWP = 0.4432; // medium working point
  
  // decleare reader
  TTreeReader myReader(tree);
  TTreeReaderValue<std::vector<float> > combinejetpt         (myReader,"combinejetpt");
  TTreeReaderValue<std::vector<float> > combinejeteta        (myReader,"combinejeteta");
  TTreeReaderValue<std::vector<float> > combinejetHFlav      (myReader,"combinejetHFlav");
  TTreeReaderValue<std::vector<float> > combinejetbtag       (myReader,"combinejetbtag");
  TTreeReaderValue<std::vector<float> > combinejetBtagSF     (myReader,"combinejetBtagSF");
  TTreeReaderValue<std::vector<float> > combinejetBtagSFUp   (myReader,"combinejetBtagSFUp");
  TTreeReaderValue<std::vector<float> > combinejetBtagSFDown (myReader,"combinejetBtagSFDown");

  // add a b-tag weight branch for medium wp
  float wgtbtag, wgtbtagUp, wgtbtagDown;
  TBranch* bwgtbtag     = tree->Branch("wgtbtag", &wgtbtag, "wgtbtag/F");  
  TBranch* bwgtbtagUp   = tree->Branch("wgtbtagUp", &wgtbtagUp, "wgtbtagUp/F");  
  TBranch* bwgtbtagDown = tree->Branch("wgtbtagDown", &wgtbtagDown, "wgtbtagDown/F");  

  long int nEvents = 0;
  while(myReader.Next()){

    // recipe for b-tag weight on https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
    vector<float> PMC;
    vector<float> PDATA ;
    vector<float> PDATAErrUp ;
    vector<float> PDATAErrDw ;
    float efficiency = 1.;

    if(nEvents %10000 == 0){
      std::cout<<"Events "<<float(nEvents)/tree->GetEntries()*100<<"%"<<std::endl;
    }      

    // take only events in the acceptance region for b-tagging   
    for(size_t iJet = 0; iJet < combinejetpt->size(); iJet++){            
      if(combinejetpt->at(iJet) < jetPtMin) continue;
      if(fabs(combinejeteta->at(iJet)) > jetEtaMax) continue;

      // get efficiency --> get a value for the efficiency interpolating the histogram
      if(combinejetHFlav->at(iJet) == 5)
	efficiency     = graph_b->Interpolate(combinejetpt->at(iJet) ,fabs(combinejeteta->at(iJet)));
      else if(combinejetHFlav->at(iJet) == 4)
	efficiency     = graph_c->Interpolate(combinejetpt->at(iJet),fabs(combinejeteta->at(iJet)));
      else
	efficiency     = graph_ucsdg->Interpolate(combinejetpt->at(iJet),fabs(combinejeteta->at(iJet)));

      //checks
      if(efficiency > 1){
	cerr<<"Problem b-tagging efficiency larger than 1 for this bin: pt "<<combinejetpt->at(iJet)<<" eta "<<fabs(combinejeteta->at(iJet)) <<endl;
	efficiency = 1;
      }
      if(efficiency < 0){
	cerr<<"Problem b-tagging efficiency larger smaller than 0 for this bin: pt "<<combinejetpt->at(iJet)<<" eta "<<fabs(combinejeteta->at(iJet)) <<endl;
	efficiency = 0;
      }
	
      // check b-tag value
      if(combinejetbtag->at(iJet) > btagCSVMediumWP){
	PMC.push_back(efficiency);
	PDATA.push_back(efficiency*combinejetBtagSF->at(iJet));		  
	PDATAErrUp.push_back(efficiency*combinejetBtagSFUp->at(iJet));
	PDATAErrDw.push_back(efficiency*combinejetBtagSFDown->at(iJet));
      }
      else{
	PMC.push_back(1-efficiency);
	PDATA.push_back(1-efficiency*combinejetBtagSF->at(iJet));
	PDATAErrUp.push_back(1-efficiency*combinejetBtagSFUp->at(iJet));
	PDATAErrDw.push_back(1-efficiency*combinejetBtagSFDown->at(iJet));
	}	
    }
  
  
    // Fill branch
    float Num = 1.;
    float Den = 1.;

    for(size_t iJet = 0; iJet < PMC.size(); iJet++){
      Num *= PDATA.at(iJet); // probability of data
      Den *= PMC.at(iJet);   // probability of MC
    }
    
    if(Den != 0)
      wgtbtag = Num/Den;
    else{
      cerr<<"Problem : null probability for denominator PMC --> fix weight to 1 "<<endl;
      wgtbtag = 1.;
    }
    
    // uncertainty 
    vector<float> wbtagErr;
    float NumMax = 1.;
    float NumMin = 1.;
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

  std::cout<<"###################################"<<std::endl;
  std::cout<<"btagWeights function --> end"<<std::endl;
  std::cout<<"###################################"<<std::endl;

}


// function that take as input a ROOT file
void qcdfilter( std::string inputFileName,  // name of a single file or directory path
		std::string outputFileName,  // output file name --> single file
		bool isMC, // is data or MC
		bool applyBTagWeights, // store b-tag weights
		bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory
		bool isEOS, // if the directory is on eos or not
		int  xsType = 0, 		
		bool storeGenTree = false, // store gentree in the output
		bool dropHLTFilter = false
		) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"qcdfilter --> start function"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "qcdtree.root";

  TFile*  infile = NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;

  // single file as input --> load the tree
  if(not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  // fill the chain with the files
  else{    
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList){
      frtree->Add(name.c_str());    
    }
  }
  
  TChain* intree  = new TChain("gentree/gentree");

  TH1D*   puRatio = NULL;
  double  wgtsum;

  // only for MC take the gentree to evaluate sum of weights and true pileup one
  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      for(auto name : fileList)
	intree->Add(name.c_str());    
    }    
    // calculate sum of weights
    wgtsum = sumwgt(intree);
    // caluclate puweight
    puRatio = pileupwgt(intree);
  }

  string cut = "nmuons == 0 && nelectrons == 0 && nphotons == 0 && nbjetslowpt == 0 && ntausold == 0 && (flaghbhenoise == 1 && flaghbheiso == 1 && flageebadsc == 1 && flagecaltp == 1 && flaggoodvertices == 1 && flagglobaltighthalo == 1 && flagbadchpf == 1 && flagbadpfmu == 1) && combinejetpt[0] > 80 && abs(combinejeteta[0]) < 4.7";
  if (not dropHLTFilter)
    cut += " && (hltPFHT350 > 0 || hltPFHT400 > 0 || hltPFHT475 > 0 || hltPFHT600 > 0 || hltPFHT650 > 0 || hltPFHT800 > 0 || hltPFHT900 > 0)";
  

  // output file
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree"); 
  treedir->cd();

  // copy the tree applying a selection
  std::cout<<"qcdfilter --> apply signal region preselection"<<std::endl;
  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  frtree->SetBranchStatus("*SubJet*",0);
  frtree->SetBranchStatus("*softDrop*",0);
  frtree->SetBranchStatus("mu*",0);
  frtree->SetBranchStatus("el*",0);

  // copy output tree from input apply selections and discarding useless branches
  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"qcdfilter --> outtree events "<<outtree->GetEntries()<<std::endl;
  
  // add a weight sum branch  
  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;

  float wgtpileup = 1;
  double xsec;

  // add wgtsum, pileupwgt and xsec if needed
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;
  
  if(isMC){
    bwgtsum    = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");      

    if(xsType > 0 and isMC)
      bxsec = outtree->Branch("xsec", &xsec, "xsec/F");
    
    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");
    
    std::cout<<"qcdfilter --> apply sumwgt and puweight"<<std::endl;
    long int nEvents=0;
    while(myReader.Next()){
      std::cout.flush();
      if(nEvents %100000 == 0) 
	std::cout<<"\r"<<"Event: "<<100*float(nEvents)/outtree->GetEntries()<<" %";
      bwgtsum->Fill();
      wgtpileup = puRatio->GetBinContent(puRatio->FindBin(*putrue));
      bwgtpileup->Fill();    
      if(xsType > 0 and isMC)
	bxsec->Fill();
      nEvents++;
    }
  }
  
  if(applyBTagWeights and isMC){

    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

    // in case it's not a list of files
    if(not isInputDirectory){
      eff_Num_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"); 
      eff_Num_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"); 
      eff_Num_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"); 
      // take denominators
      eff_Denom_b = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"); 
      eff_Denom_c = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"); 
      eff_Denom_ucsdg = (TH2F*) infile->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"); 
    }
    else{

      // loop on the file --> open each at a time     
      TFile* file_temp = 0;
      for(auto name : fileList){
	file_temp = TFile::Open(name.c_str());
	if(not eff_Num_b){
	  eff_Num_b = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b")->Clone("eff_Num_b");
	  eff_Num_b->SetDirectory(0);
	}
	else
	  eff_Num_b->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_b"));

	if(not eff_Num_c){
          eff_Num_c = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c")->Clone("eff_Num_c");
	  eff_Num_c->SetDirectory(0);
	}
	else
          eff_Num_c->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_c"));

	if(not eff_Num_ucsdg){
          eff_Num_ucsdg = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg")->Clone("eff_Num_ucsdg");
	  eff_Num_ucsdg->SetDirectory(0);
	}
	else
          eff_Num_ucsdg->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Num_ucsdg"));

	if(not eff_Denom_b){
	  eff_Denom_b = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b")->Clone("eff_Denom_b");
	  eff_Denom_b->SetDirectory(0);
	}
	else
	  eff_Denom_b->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_b"));

	if(not eff_Denom_c){
          eff_Denom_c = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c")->Clone("eff_Denom_c");
	  eff_Denom_c->SetDirectory(0);
	}
	else
          eff_Denom_c->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_c"));

	if(not eff_Denom_ucsdg){
          eff_Denom_ucsdg = (TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg")->Clone("eff_Denom_ucdsg");
	  eff_Denom_ucsdg->SetDirectory(0);
	}
	else
          eff_Denom_ucsdg->Add((TH2F*) file_temp->Get("btageff/eff_pfCombinedInclusiveSecondaryVertexV2BJetTags_Medium_Denom_ucsdg"));	

	// close the file
	file_temp->Close();
      }
    }

    // compute efficiency
    TH2F* eff_b = (TH2F*) eff_Num_b->Clone("eff_b");
    TH2F* eff_c = (TH2F*) eff_Num_c->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);

    eff_b->Smooth();
    eff_c->Smooth();
    eff_ucsdg->Smooth();

    btagWeights(outtree,eff_b,eff_c,eff_ucsdg);


  }
  
  // write tree in the file  
  outfile->cd();
  if(storeGenTree){
    TDirectoryFile* gentreedir = new TDirectoryFile("gentree", "gentree"); 
    gentreedir->cd();
    intree->CloneTree()->Write();
  }
  treedir->cd();
  outtree->Write();
  outfile->Close();
  std::cout<<"###################################"<<std::endl;
}
