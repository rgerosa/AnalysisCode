// function that takes in input a tree and counts the total of negative and positive signs
double sumwgt(TTree* tree) {

  if(tree == NULL or tree == 0)
    return 1.;
  
  TTreeReader treeReader(tree);
  TTreeReaderValue<double> wgtsign (treeReader,"wgtsign");

  double weightsum = 0.;
  while(treeReader.Next()){
    weightsum += (*wgtsign);
  }
  return weightsum;
}

//function to apply pileup re-weight, i.e. return an histogram given by the ratio between pileup in data and mc
TH1D* pileupwgt(TTree* tree, std::string scenario){

  if(tree == NULL or tree == 0)
    return 0;
  
  // PU data  
  TH1D* histoPUData;
  TFile* fileInputData;

  if(scenario == "silver"){
    fileInputData = TFile::Open("../data/JSON_PILEUP/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.root");
    histoPUData = (TH1D*) fileInputData->Get("pileup");
  }
  else{
    fileInputData = TFile::Open("../data/JSON_PILEUP/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.root");
    histoPUData = (TH1D*) fileInputData->Get("pileup");
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

  return puRatio;

}

// function that take as input a ROOT file
void sigfilter( std::string inputFileName,  std::string outputFileName ) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "sigtree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");

  // calculate sum of weights
  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
 
 double wgtpileup = 0;

  // accept signal region event if loose muons, electrons, taus and phototns == 0 and triggered by hltmet90
  const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && hltmet90 > 0";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  // copy the tree applying a selection
  TTree* outtree = frtree->CopyTree(cut);
  
  // add a weight sum branch  
  TBranch* bwgtsum    = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  // add a pileup weight branch
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }
  
  // write tree in the file
  outfile->Write();
  
}

// function to apply Zmumu selections
void zmmfilter( std::string inputFileName,  std::string outputFileName) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zmmtree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");

  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
  double wgtpileup = 0;
    
  // Selections 2 loose muons, 0 loose ele, taus and photons, m(mumu) = m(z) [60,120], mupt > 20, one of the two tight ... no trigger requirement
  const char* cut = "nmuons == 2 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && zmass > 60 && zmass < 120 && mu1pt > 20 && (mu1id == 1 || mu2id == 1)";

  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  TTree* outtree = frtree->CopyTree(cut);

  TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }


  outfile->Write();

}

// function to apply Zee selections
void zeefilter( std::string inputFileName,  std::string outputFileName) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "zeetree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  
  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
  double wgtpileup = 0;
  
  const char* cut = "nmuons == 0 && nelectrons == 2 && ntaus == 0 && nphotons == 0 && zeemass > 60 && zeemass < 120 && el1pt > 40 && (el1id == 1 || el2id == 1)";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  TTree* outtree = frtree->CopyTree(cut);
  
  TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }

  outfile->Write();

}

// function to apply Wmunu selections
void wmnfilter( std::string inputFileName,  std::string outputFileName) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wmntree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  
  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
  double wgtpileup = 0;

  const char* cut = "nmuons == 1 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && mu1pt > 20 && mu1id == 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  TTree* outtree = frtree->CopyTree(cut);

  TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }
  outfile->Write();
  
}

// function to apply Wenu selections
void wenfilter( std::string inputFileName,  std::string outputFileName) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "wentree.root";

  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  
  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
  double wgtpileup = 0;
  
  const char* cut = "nmuons == 0 && nelectrons == 1 && ntaus == 0 && nphotons == 0 && el1pt > 40 && el1id == 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  TTree* outtree = frtree->CopyTree(cut);

  TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }

  outfile->Write();
}

// function to apply photon+jets selections
void gamfilter( std::string inputFileName,  std::string outputFileName) {

  if(inputFileName == "")
    inputFileName = "tree.root";
  if(outputFileName == "")
    outputFileName = "gamtree.root";
  
  TFile* infile = TFile::Open(inputFileName.c_str());
  TTree* intree = (TTree*)infile->Get("gentree/gentree");
  TTree* frtree = (TTree*)infile->Get("tree/tree");
  
  double wgtsum = sumwgt(intree);
  // caluclate puweight
  std::string scenario = "silver";
  TH1D*  puRatio = pileupwgt(intree,scenario);
  double wgtpileup = 0;

  // medium id + pt + veto
  const char* cut = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 1 && phpt > 175 && phidm == 1";
  
  TFile* outfile = new TFile(outputFileName.c_str(), "RECREATE");
  outfile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");
  treedir->cd();
  TTree* outtree = frtree->CopyTree(cut);

  TBranch* bwgtsum = outtree->Branch("wgtsum", &wgtsum, "wgtsum/D");
  TBranch* bwgtpileup = outtree->Branch("wgtpileup", &wgtpileup, "wgtpileup/D");  

  TTreeReader myReader(outtree);
  TTreeReaderValue<int> putrue(myReader,"putrue");

  while(myReader.Next()){
    bwgtsum->Fill();
    wgtpileup = puRatio->GetBinContent(*putrue);
    bwgtpileup->Fill();
  }
  
  outfile->Write();

}

