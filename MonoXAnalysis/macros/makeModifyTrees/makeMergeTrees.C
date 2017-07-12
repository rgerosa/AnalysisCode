
void makeMergeTrees(string inputDirectory, bool isEOS, string outputDIR, string outputFileName){
  
  gROOT->SetBatch(kTRUE);

  // take input files
  TChain* inputTree = new TChain("tree/tree");

  vector<string> inputFileName;

  // make a list of files
  if(isEOS)
    system(("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+inputDirectory+" -name \"*.root\" | grep -v failed > file.temp").c_str());
  else
    system(("find "+inputDirectory+" -name \"*.root\" | grep -v failed > file.temp").c_str());

  ifstream infile;
  string line;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      inputFileName.push_back(line);
      if(isEOS)
        inputTree->Add(("root://eoscms.cern.ch//"+line).c_str());
      else
        inputTree->Add(line.c_str());
    }
  }

  system("rm file.temp");
  
  // create ouytput file
  TFile* outputFile = new TFile((outputDIR+"/"+outputFileName).c_str(),"RECREATE");
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
  outputFile->Close();
  
  // remove old files
  if(not isEOS){
    for(auto file : inputFileName){
      system(("rm "+file).c_str());
    }
  }
}
