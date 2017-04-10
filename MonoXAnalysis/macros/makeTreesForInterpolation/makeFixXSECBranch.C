
void makeFixXSECBranch(string inputFileName, string outputFileName, vector<pair<float,float> > massPoint, vector<float> xsecVal){

  if(massPoint.size() != xsecVal.size()){
    cerr<<" mediator mass point , DM mass points and xsec values must have same size --> break "<<endl;
    return;
  }

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("tree");

  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  
  // not copy xsec branch
  tree->SetBranchStatus("xsec",0);
  cout<<"Copying tree "<<endl;
  TTree* treeOut = tree->CopyTree("");

  // set it on again
  tree->SetBranchStatus("xsec",1);
  TTreeReader reader (tree);
  TTreeReaderValue<float>  xsec (reader,"xsec");
  TTreeReaderValue<float>  genMediatorMass (reader,"genMediatorMass");
  TTreeReaderValue<float>  genX1Mass (reader,"genX1Mass");

  // add xsec branch to output tree
  float xsection = 0;
  TBranch* xsec_b = treeOut->Branch("xsec",&xsection,"xsec/F");
  
  cout<<"Looping on events to fix xsec value "<<endl;
  long int nEvents = 0;
  long int nTotal = tree->GetEntries();
  while(reader.Next()){

    if(int(nEvents) %10000 == 0){
      std::cout.flush();
      std::cout<<"\r"<<"Events "<<float(nEvents)/float(nTotal)*100<<"%";
    }
    nEvents++;

    int posFound = -1;
    int pos = 0;
    for(auto medmass : massPoint){
      if(*genMediatorMass == medmass.first and *genX1Mass == medmass.second){
	posFound = pos;
	break;
      }
      pos++;
    }

    if(posFound == -1)
      xsection = *xsec;
    else
      xsection = xsecVal.at(posFound)*1000;

    xsec_b->Fill();    
  }

  treeOut->Write();
  outputFile->Close();

}
