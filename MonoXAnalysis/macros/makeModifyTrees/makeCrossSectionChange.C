
void makeCrossSectionChange(string fileName, string outputName, float newXSEC, string treename = "tree/tree"){

  gROOT->SetBatch(kTRUE);

  TFile* input = TFile::Open(fileName.c_str(),"READ");
  TTree* tree = (TTree*) input->Get(treename.c_str());
  TTree* gentree = (TTree*) input->Get("gentree/gentree");
  
  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  output->cd();

  // clone tree skipping xsec branch
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree"); 
  treedir->cd();
  tree->SetBranchStatus("xsec",0);  
  TTree* outtree = tree->CopyTree("");
  float xsec = 0;
  TBranch* bxsec = outtree->Branch("xsec", &xsec, "xsec/F");

  TTreeReader myReader(outtree);
  while(myReader.Next()){
    xsec = newXSEC*1000;
    bxsec->Fill();
  }

  treedir->cd();
  outtree->Write();
  if(gentree){
    TDirectoryFile* gentreedir = new TDirectoryFile("gentree", "gentree"); 
    gentreedir->cd();
    gentree->CloneTree()->Write();
  }

  output->Close();
  
}
