
void makeMergeExtendedSamples(vector<string> inputFileName, string outputFileName){

  gROOT->SetBatch(kTRUE);

  // take input files
  TChain* inputTree = new TChain("tree/tree");
  TChain* inputGenTree = new TChain("gentree/gentree");

  for(auto file : inputFileName){
    inputTree->Add(file.c_str());
    inputGenTree->Add(file.c_str());
  }

  // create ouytput file
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  TDirectoryFile* treedir = new TDirectoryFile("tree", "tree");   
  treedir->cd();

  // calculate the wgtsum
  TTreeReader myReader(inputTree);
  TTreeReaderValue <double> wgt (myReader,"wgtsum");
  TTreeReaderValue <float> xsec (myReader,"xsec");

  cout<<"Calculate the sumwgt "<<endl;
  string currentFile = "";
  double wgtsum = 0;  
  while(myReader.Next()){
    if(dynamic_cast<TChain*>(myReader.GetTree())->GetFile()->GetName() != currentFile){
      currentFile = dynamic_cast<TChain*>(myReader.GetTree())->GetFile()->GetName();
      wgtsum += *wgt;
    }
    else
      continue;
  }
  
  // copy sipping sumwgt
  cout<<"Clone tree "<<endl;
  inputTree->SetBranchStatus("*wgtsum*",0);
  TTree* outputTree = inputTree->CopyTree("");

  // create the new branch
  TBranch* bwgtsum = outputTree->Branch("wgtsum", &wgtsum, "wgtsum/D");
 
  //
  myReader.SetEntry(0);
  cout<<"Store wgtsum in the tree "<<endl;
  while(myReader.Next()){
    bwgtsum->Fill();
  }

  // 
  cout<<"Store output file "<<endl;
  outputFile->cd();
  treedir->cd();
  outputTree->Write();

  outputFile->cd();
  if(inputGenTree->GetEntries() > 0){
    TDirectoryFile* gentreedir = new TDirectoryFile("gentree", "gentree"); 
    gentreedir->cd();
    inputGenTree->CloneTree()->Write();
  }
  outputFile->Close();

  // remove old files
  for(auto file : inputFileName){
    system(("rm "+file).c_str());
  }
}
