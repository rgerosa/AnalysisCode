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

void filters(){}

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
  
  fileInputData = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/JSON_PILEUP/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.root");
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
void sigfilter( std::string inputFileName,  // name of a single file or directory path
		std::string outputFileName,  // output file name --> single file
		bool isMC, // is data or MC
		bool applyBTagWeights, // store b-tag weights
		bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory
		bool isEOS, // if the directory is on eos or not
		int  xsType = 0, 		
		bool storeGenTree = false, // store gentree in the output
		bool dropPuppiBranches = true,
		bool dropPuppiBoostedJets = true,
		bool dropSubJetsBranches = true,
		bool dropHLTFilter = false,
		string metCut = "175"
		) {

  // if xsType = 0: means keep the value inside of xsec branch in the standard tree
  // if xsType = 1: evaluate the xs as sumwgt/Nevents (useful for powheg+minlo mono-jet samples)
  // if xsType = 2: take the cross section from gentree lheXSEC read from the lhe file, useful for madgraph samples

  std::cout<<"###################################"<<std::endl;
  std::cout<<"sigfilter --> start function"<<std::endl;
  std::cout<<"###################################"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "sigtree.root";

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

  // accept signal region event if loose muons, electrons, taus and phototns == 0 and triggered by hltmet90
  string cut = "nmuons == 0 && nelectrons == 0 && nphotons == 0 && t1pfmet > "+metCut;
  if (not dropHLTFilter)
    cut += " && (hltmet90 > 0 || hltmet100 > 0 || hltmet110 > 0 || hltmet120 > 0 || hltmetwithmu90 > 0 ||  hltmetwithmu100 > 0 || hltmetwithmu110 || hltmetwithmu120 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltjetmet > 0 )";
      
  // output file
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree"); 
  treedir->cd();

  // copy the tree applying a selection
  std::cout<<"sigfilter --> apply signal region preselection"<<std::endl;
  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  // copy output tree from input apply selections and discarding useless branches
  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"sigfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  // add a weight sum branch  
  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;

  double wgtpileup = 1;
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
    
    std::cout<<"sigfilter --> apply sumwgt and puweight"<<std::endl;
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
    TH2F* eff_c = (TH2F*) eff_Num_b->Clone("eff_c");
    TH2F* eff_ucsdg = (TH2F*) eff_Num_ucsdg->Clone("eff_ucsdg");
    eff_b->Divide(eff_Denom_b);
    eff_c->Divide(eff_Denom_c);
    eff_ucsdg->Divide(eff_Denom_ucsdg);

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

// function to apply Zmumu selections
void zmmfilter(std::string inputFileName,  // name of a single file or directory path                                                                                           
	       std::string outputFileName,  // output file name --> single file                                                                                                 
	       bool isMC, // is data or MC                                                                                                                                      
	       bool applyBTagWeights, // store b-tag weights                                                                                                                    
	       bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                       
	       bool isEOS, // if the directory is on eos or not                                                                                                                 
	       int  xsType = 0,
	       bool storeGenTree = false, // store gentree in the output                                                                                                      
	       bool dropPuppiBranches = true,
	       bool dropPuppiBoostedJets = true,
	       bool dropSubJetsBranches = true,
	       bool dropHLTFilter = false,
	       string metCut = "175"
	       ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"zmmfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zmmtree.root";

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }
  

  TChain* intree  = new TChain("gentree/gentree");
  TH1D*  puRatio = NULL;
  double wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                                 
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                       
    puRatio = pileupwgt(intree);
  }

    
  // Selections 2 loose muons, 0 loose ele, taus and photons, m(mumu) = m(z) [60,120], mupt > 20, one of the two tight
  string cut = "nmuons == 2 && nelectrons == 0 && nphotons == 0 && zmass > 60 && zmass < 120 && ((mu1pt > 20 && mu1id >= 1) || (mu2pt > 20 && mu2id >= 1)) && t1mumet > "+metCut+" && (mu1pid != mu2pid)";
  if(not dropHLTFilter)
    cut += " && (hltmet90 > 0 || hltmet100 > 0 || hltmet110 > 0 || hltmet120 > 0 || hltmetwithmu90 || hltmetwithmu100 || hltmetwithmu110 ||  hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0 || hltsinglemu > 0 || hltjetmet > 0)";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"zmmfilter --> apply signal region preselection"<<std::endl;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);


  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"zmmfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  // add a weight sum branch                                                                                                                                                    
  TBranch* bwgtsum;
  TBranch *bwgtpileup; 
  TBranch *bxsec = NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){

    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"zmmfilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"zmmfilter --> apply btag-weight"<<std::endl;
 
    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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
        file_temp->Close();
      }
    }
    
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

// function to apply Zee selections
void zeefilter(std::string inputFileName,  // name of a single file or directory path                                                                                           
	       std::string outputFileName,  // output file name --> single file                                                                                                 
	       bool isMC, // is data or MC                                                                                                                                      
	       bool applyBTagWeights, // store b-tag weights                                                                                                                    
	       bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                       
	       bool isEOS, // if the directory is on eos or not                                                                                                                 
	       int  xsType = 0,
	       bool storeGenTree = false, // store gentree in the output                                                                                                        
	       bool isSinglePhoton = false, // to use also single photon trigger when running on data: singleEle trigger on singleEle dataset, photon trigger on photon dataset
	       bool isJetHT  = false, // to use also single photon trigger when running on data: singleEle trigger on singleEle dataset, photon trigger on photon dataset
	       bool isDoubleEG = false,
	       bool dropPuppiBranches = true,
	       bool dropPuppiBoostedJets = true,
	       bool dropSubJetsBranches = true,
	       bool dropHLTFilter = false,
	       string metCut = "175"
	       ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"zeefilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zeetree.root";

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList ;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }

  TChain* intree  = new TChain("gentree/gentree");
  TH1D*  puRatio = NULL;
  double wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                                 
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                       
    puRatio = pileupwgt(intree);
  }

  string cut = "nmuons == 0 && nelectrons == 2  && nphotons == 0 && zeemass > 60 && zeemass < 120 && ((el1pt > 40 && el1id >= 1) || (el2pt > 40 && el2id >= 1)) && t1elmet > "+metCut+" && (el1pid != el2pid)";
  
  // trigger 
  if(not isMC and not isSinglePhoton and not isJetHT and not isDoubleEG and not dropHLTFilter)
    cut += " && (hltsingleel || hltelnoiso)";
  else if(not isMC and isSinglePhoton and not isJetHT and not isDoubleEG and not dropHLTFilter)
    cut += " && ( hltphoton165 || hltphoton175 ) && (hltsingleel == 0 && hltelnoiso == 0)";      
  else if(not isMC and isJetHT and not isDoubleEG and not isSinglePhoton and not dropHLTFilter)
    cut += " && (hltPFHT800 || hltEcalHT800) && (hltsingleel == 0 && hltelnoiso == 0)";    
  else if(not isMC and not isJetHT and isDoubleEG and not isSinglePhoton and not dropHLTFilter)
    cut += " && (hltPFHT800 || hltEcalHT800) && (hltsingleel == 0 && hltelnoiso == 0)";    
  else if(isMC and not dropHLTFilter)
    cut += " && (hltsingleel > 0 || hltphoton175 > 0 || hltphoton165 > 0 || hltelnoiso > 0 || hltPFHT800 > 0 || hltEcalHT800 > 0)";
  

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"zeefilter --> apply signal region preselection"<<std::endl;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"zeefilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;

  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"zeefilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"zeefilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }

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

// function to apply Wmunu selections
void wmnfilter(std::string inputFileName,  // name of a single file or directory path                                                                                           
	       std::string outputFileName,  // output file name --> single file                                                                                                 
	       bool isMC, // is data or MC                                                                                                                                      
	       bool applyBTagWeights, // store b-tag weights                                                                                                                    
	       bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                       
	       bool isEOS, // if the directory is on eos or not                                                                                                                 
	       int  xsType = 0,
	       bool storeGenTree = false, // store gentree in the output                                                                                                        
	       bool dropPuppiBranches = true,
	       bool dropPuppiBoostedJets = true,
	       bool dropSubJetsBranches = true,
	       bool dropHLTFilter = false,
	       string metCut = "175"
	       ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"wmnfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wmntree.root";

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }
  
  TChain* intree  = new TChain("gentree/gentree");
  TH1D*   puRatio = NULL;
  double  wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                                 
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                       
    puRatio = pileupwgt(intree);
  }


  string cut = "nmuons == 1 && nelectrons == 0  && nphotons == 0 && mu1pt > 20 && mu1id >= 1 && t1mumet > "+metCut;
  if(not dropHLTFilter)
    cut += " && (hltmet90 > 0 || hltmet100 > 0 || hltmet110 > 0 || hltmet120 > 0 || hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0 || hltmetwithmu100 > 0 || hltmetwithmu110 > 0 || hltsinglemu > 0 || hltjetmet > 0)";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"wmnfilter --> apply signal region preselection"<<std::endl;
  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);

  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"wmnfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  if(isMC){

    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);    
    TTreeReaderValue<int> putrue(myReader,"putrue");
    std::cout<<"wmnfilter --> apply sumwgt and puweight"<<std::endl;    
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

    std::cout<<"wmnfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }
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

// function to apply Wenu selections
void wenfilter(std::string inputFileName,  // name of a single file or directory path                                                                                           
	       std::string outputFileName,  // output file name --> single file                                                                                                 
	       bool isMC, // is data or MC                                                                                                                                      
	       bool applyBTagWeights, // store b-tag weights                                                                                                                    
	       bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                       
	       bool isEOS, // if the directory is on eos or not                                                                                                                 
	       int  xsType = 0,
	       bool storeGenTree = false, // store gentree in the output                                                                                                        
	       bool isSinglePhoton = false,
	       bool isJetHT = false,
	       bool isDoubleEG = false,
	       bool dropPuppiBranches = true,
	       bool dropPuppiBoostedJets = true,
	       bool dropSubJetsBranches = true,
	       bool dropHLTFilter = false,
	       string metCut = "175") {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"wenfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wentree.root";

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }


  TChain* intree  = new TChain("gentree/gentree");
  TH1D*   puRatio = NULL;
  double  wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                               
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                       
    puRatio = pileupwgt(intree);
  }


  string cut = "nmuons == 0 && nelectrons == 1  && nphotons == 0 && el1pt > 40 && el1id >= 1 && t1elmet > "+metCut;

  if(not isMC and not isSinglePhoton and not isJetHT and not isDoubleEG and not dropHLTFilter)
      cut += " && (hltsingleel >0 || hltelnoiso > 0)";  
  else if(not isMC and isSinglePhoton and not isJetHT and not isDoubleEG and not dropHLTFilter)
    cut += " && (hltsingleel == 0 && hltelnoiso == 0) && (hltphoton165 > 0 || hltphoton175 > 0)";
  else if(not isMC and isJetHT and not isSinglePhoton and not isDoubleEG and not dropHLTFilter)
    cut += " && (hltsingleel == 0 && hltelnoiso == 0) && (hltEcalHT800 > 0 || hltPFHT800 > 0)";
  else if(not isMC and not isJetHT and not isSinglePhoton and isDoubleEG and not dropHLTFilter)
    cut += " && (hltsingleel == 0 && hltelnoiso == 0) && (hltEcalHT800 > 0 || hltPFHT800 > 0)";
  else if(isMC and not dropHLTFilter)
    cut += " && (hltsingleel > 0 || hltelnoiso > 0 || hltphoton165 > 0 || hltphoton175 > 0 || hltEcalHT800 > 0 || hltPFHT800 > 0)";

  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"wenfilter --> apply signal region preselection"<<std::endl;


  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"wenfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec =  NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){
    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"wenfilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"wenfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }

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

// function to apply photon+jets selections
void gamfilter(std::string inputFileName,  // name of a single file or directory path                                                                                           
	       std::string outputFileName,  // output file name --> single file                                                                                                 
	       bool isMC, // is data or MC                                                                                                                                      
	       bool applyBTagWeights, // store b-tag weights                                                                                                                    
	       bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                       
	       bool isEOS, // if the directory is on eos or not                                                                                                                
	       int  xsType = 0,
	       bool storeGenTree = false, // store gentree in the output                                                                                                   
	       bool isJetHT = false,
	       bool isDoubleEG = false,
	       bool dropPuppiBranches = true,
	       bool dropPuppiBoostedJets = true,
	       bool dropSubJetsBranches = true,
	       bool dropHLTFilter = false,
	       string metCut = "175"
	       ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"gamfilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "gamtree.root";
  

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }

  TChain* intree  = new TChain("gentree/gentree");
  TH1D*   puRatio = NULL;
  double  wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                             
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                       
    puRatio = pileupwgt(intree);
  }


  // medium id + pt + veto
  string cut = "nmuons == 0 && nelectrons == 0  && nphotons == 1 && phpt > 120 && phidm == 1 && t1phmet > "+metCut;
  if(not isMC and not dropHLTFilter and not isJetHT and not isDoubleEG)
    cut += " && (hltphoton165 || hltphoton175 || hltphoton120)";
  else if(not isMC and not dropHLTFilter and isJetHT and not isDoubleEG)
    cut += " && (hltphoton165 == 0 && hltphoton175 == 0) && (hltEcalHT800 || hltPFHT800 > 0)";
  else if(not isMC and not dropHLTFilter and not isJetHT and isDoubleEG)
    cut += " && (hltphoton165 == 0 && hltphoton175 == 0) && (hltEcalHT800 || hltPFHT800 > 0)";
  else if(isMC and not dropHLTFilter)
    cut += " && (hltphoton165 || hltphoton175 || hltphoton120 || hltEcalHT800 || hltPFHT800 > 0)";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"gamfilter --> apply signal region preselection"<<std::endl;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"gamfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){

    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"gamfilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"gamfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }

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



// function to apply photon+jets selections
void topmufilter(std::string inputFileName,  // name of a single file or directory path                                                                                        
		 std::string outputFileName,  // output file name --> single file                                                                                            
		 bool isMC, // is data or MC                                                                                                                                  
		 bool applyBTagWeights, // store b-tag weights                                                                                                                
		 bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                   
		 bool isEOS, // if the directory is on eos or not                                                                                                             
		 int  xsType = 0,
		 bool storeGenTree = false, // store gentree in the output                                                                                                     
		 bool dropPuppiBranches = true,
		 bool dropPuppiBoostedJets = true,
		 bool dropSubJetsBranches = true,
		 bool dropHLTFilter = false,
		 string metCut = "175"
	       ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"topmufilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "topmutree.root";

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if( not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }
  

  TChain* intree  = new TChain("gentree/gentree");
  TH1D*   puRatio = NULL;
  double  wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }
    // calculate sum of weights                                                                                                                                             
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                      
    puRatio = pileupwgt(intree);
  }


  // one tight muon + b-jet --> semi-leptonic ttbar events
  string cut = "nmuons == 1 && nelectrons == 0  && nphotons == 0 && nbjets > 0 && mu1id >=1 && mu1pt > 20 && t1mumet > "+metCut;
  if(not dropHLTFilter)
    cut += " && (hltmet90 > 0 || hltmet100 > 0 || hltmet110 > 0 || hltmet120 > 0 || hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0 || hltmetwithmu100 > 0 || hltmetwithmu110 > 0 || hltsinglemu > 0 || hltjetmet > 0)";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"topfilter --> apply signal region preselection"<<std::endl;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);

  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"topfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){

    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"topfilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"topfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights
    // take numerators for b-tag
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }
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



// function to apply photon+jets selections
void topelfilter(std::string inputFileName,  // name of a single file or directory path                                                                                        
		 std::string outputFileName,  // output file name --> single file                                                                                            
		 bool isMC, // is data or MC                                                                                                                                  
		 bool applyBTagWeights, // store b-tag weights                                                                                                                
		 bool isInputDirectory, // to tell wether the inputFileName is a single file or a directory                                                                   
		 bool isEOS, // if the directory is on eos or not                                                                                                             
		 int  xsType = 0,
		 bool storeGenTree = false, // store gentree in the output                                                                                                     
		 bool isSinglePhoton = false, 
		 bool isJetHT  = false, 
		 bool isDoubleEG  = false, 
		 bool dropPuppiBranches = true,
		 bool dropPuppiBoostedJets = true,
		 bool dropSubJetsBranches = true,
		 bool dropHLTFilter = false,
		 string metCut = "175"
		 ) {

  std::cout<<"###################################"<<std::endl;
  std::cout<<"topmufilter --> start function"<<std::endl;

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "topeltree.root";
  

  TFile* infile  =  NULL;
  TChain* frtree = new TChain("tree/tree");
  vector<string> fileList;
  if(not isInputDirectory){
    infile = TFile::Open(inputFileName.c_str());
    frtree = (TChain*)infile->Get("tree/tree");
  }
  else{
    frtree = new TChain("tree/tree");
    fileList = fileListForChain(inputFileName,isEOS);
    for(auto name : fileList)
      frtree->Add(name.c_str());
  }

  TChain* intree  = new TChain("gentree/gentree");
  TH1D*   puRatio = NULL;
  double  wgtsum;

  if(isMC){
    if(not isInputDirectory)
      intree = (TChain*)infile->Get("gentree/gentree");
    else{
      intree = new TChain("gentree/gentree");
      for(auto name : fileList)
        intree->Add(name.c_str());
    }

    // calculate sum of weights                                                                                                                                                
    wgtsum = sumwgt(intree);
    // caluclate puweight                                                                                                                                                     
    puRatio = pileupwgt(intree);
  }


  // one tight muon + b-jet --> semi-leptonic ttbar events
  string cut = "nmuons == 0 && nelectrons == 1  && nphotons == 0 && nbjets > 0 && el1id >=1 && el1pt > 40 && t1elmet > "+metCut;
  // trigger requirements
  if(not isMC and not dropHLTFilter and not isSinglePhoton and not isJetHT and not isDoubleEG)
    cut +=  " && (hltsingleel || hltelnoiso)";
  else if(not isMC and not dropHLTFilter and isSinglePhoton and not isJetHT and not isDoubleEG)
    cut +=  " && (hltsingleel == 0 && hltelnoiso == 0) && (hltphoton165 || hltphoton175)";
  else if(not isMC and not dropHLTFilter and not isSinglePhoton and isJetHT and not isDoubleEG)
    cut +=  " && (hltsingleel == 0 && hltelnoiso == 0) && (hltPFHT800 || hltEcalHT800)";    
  else if(not isMC and not dropHLTFilter and not isSinglePhoton and not isJetHT and isDoubleEG)
    cut +=  " && (hltsingleel == 0 && hltelnoiso == 0) && (hltPFHT800 || hltEcalHT800)";    
  else if(isMC and not dropHLTFilter)
    cut +=  " && (hltsingleel || hltelnoiso || hltphoton165 || hltphoton175 || hltPFHT800 || hltEcalHT800)";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  std::cout<<"topfilter --> apply signal region preselection"<<std::endl;

  if(xsType > 0 and isMC)
    frtree->SetBranchStatus("xsec",0);
  if(dropPuppiBranches){
    frtree->SetBranchStatus("*puppi*",0);
    frtree->SetBranchStatus("combinePuppi*",0);
    frtree->SetBranchStatus("Puppi*",0);
    frtree->SetBranchStatus("incPuppi*",0);
  }
  if(dropPuppiBoostedJets){
    frtree->SetBranchStatus("boostedPuppi*",0);
    frtree->SetBranchStatus("prunedPuppi*",0);
    frtree->SetBranchStatus("softDropPuppi*",0);
  }
  if(dropSubJetsBranches){
    frtree->SetBranchStatus("*SubJet*",0);
  }
  frtree->SetBranchStatus("emu*",0);
  frtree->SetBranchStatus("taumu*",0);
  frtree->SetBranchStatus("taue*",0);


  TTree* outtree = frtree->CopyTree(cut.c_str());
  std::cout<<"topfilter --> outtree events "<<outtree->GetEntries()<<std::endl;

  TBranch* bwgtsum = NULL;
  TBranch* bwgtpileup = NULL;
  TBranch* bxsec = NULL;
  double wgtpileup = 1;
  double xsec;
  if(xsType == 1 and isMC)
    xsec = (wgtsum/intree->GetEntries())*1000;
  else if(xsType == 2 and isMC)
    xsec = (sumxsec(intree)/intree->GetEntries())*1000;

  if(isMC){

    bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
    bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/F");  
    if(xsType > 0 and isMC)
      bxsec      = outtree->Branch("xsec", &xsec, "xsec/F");

    TTreeReader myReader(outtree);
    TTreeReaderValue<int> putrue(myReader,"putrue");

    std::cout<<"topfilter --> apply sumwgt and puweight"<<std::endl;
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
    
    std::cout<<"topfilter --> apply btag-weight"<<std::endl;

    // applying b-tag weights --> taking numerators and denominators for efficiency 
    TH2F*  eff_Num_b = 0;
    TH2F*  eff_Num_c = 0;
    TH2F*  eff_Num_ucsdg = 0;
    TH2F*  eff_Denom_b = 0;
    TH2F*  eff_Denom_c = 0;
    TH2F*  eff_Denom_ucsdg = 0;

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

        file_temp->Close();
      }
    }
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

