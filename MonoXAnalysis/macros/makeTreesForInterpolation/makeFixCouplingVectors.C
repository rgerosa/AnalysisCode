/// to re-calculate the sumwgt
class massPoint {

public:

  massPoint(){};

  massPoint(float & medmass, float & dmmass){
    medmass_ =  medmass;
    dmmass_  =  dmmass;
  }

  bool operator < (const massPoint & a) const {
    if(medmass_ < a.medmass_) return true;
    else if(medmass_ > a.medmass_) return false;
    else if(medmass_ == a.medmass_){
      if(dmmass_ <= a.dmmass_) return true;
      else if(dmmass_ > a.dmmass_) return false;
    }

    return true;
  }

  bool operator == (const massPoint & a) const{
    if(medmass_ == a.medmass_ and dmmass_ == a.dmmass_) return true;
    else return false;
  }

  bool operator != (const massPoint & a) const{
    if(medmass_ != a.medmass_ or dmmass_ != a.dmmass_) return true;
    else return false;
  }

  float  medmass_;
  float  dmmass_;
  double sumwgt_;
};

//////
void makeFixCouplingVectors(string inputFileName, string outputFileName, bool isSMM){

  // take the original tree
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("tree");
  
  //TTree Reader
  float genMediatorMass, genX1Mass, genWeight;
  vector<float> *gSM = 0;
  vector<float> *gDM = 0;
  vector<float> *couplingwgt = 0;

  tree->SetBranchAddress("genMediatorMass",&genMediatorMass);
  tree->SetBranchAddress("genX1Mass",&genX1Mass);
  tree->SetBranchAddress("genWeight",&genWeight);
  if(not isSMM){
    tree->SetBranchAddress("gSM",&gSM);
    tree->SetBranchAddress("gDM",&gDM);
    tree->SetBranchAddress("couplingwgt",&couplingwgt);
  }
  else{
    tree->SetBranchAddress("theta",&gSM);
    tree->SetBranchAddress("yDM",&gDM);
    tree->SetBranchAddress("couplingwgt",&couplingwgt);
  }
  
  
  vector<massPoint> massList;
  // first loop on the file dumping all mass points                                                                                                                                                    
  cout<<"Start loop on the file to get mass points "<<endl;
  for(long int iEvent = 0; iEvent < tree->GetEntries(); iEvent++){
    tree->GetEntry(iEvent);
    if(massList.size() == 0){
      massList.push_back(massPoint(genMediatorMass,genX1Mass));
      massList.back().sumwgt_ = 0;
      if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
        massList.back().sumwgt_ += genWeight;
    }
    else{
      if(massList.back().medmass_ == genMediatorMass and massList.back().dmmass_ == genX1Mass){
        if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
          massList.back().sumwgt_ += genWeight;
        continue; // assumption is that mass points are in sequence                                                                                                                                    
      }
      else{
        massList.push_back(massPoint(genMediatorMass,genX1Mass));
        massList.back().sumwgt_ = 0;
        if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
          massList.back().sumwgt_ += genWeight;
      }
    }
  }


  cout<<"Fix the coupling vector branch"<<endl;
  // create outputFile
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  
  // not copy xsec branch
  tree->SetBranchStatus("sumwgt",0);
  TTree* treeOut = tree->CloneTree(0);
  // re-add the branch
  double sumwgt;
  treeOut->Branch("sumwgt",&sumwgt,"sumwgt/D");
  long int nEvents = 0;
  long int nTotal = tree->GetEntries();
  for(long int iEvent = 0; iEvent < tree->GetEntries(); iEvent++){

    if(int(nEvents) %10000 == 0){
      std::cout.flush();
      std::cout<<"\r"<<"Events "<<float(nEvents)/float(nTotal)*100<<"%";
    }
    nEvents++;

    tree->GetEntry(iEvent);
    // skip bad mass points                                                                                                                                                                           
    if(gSM->size() != gDM->size() or gSM->size() != couplingwgt->size() or gDM->size() != couplingwgt->size()) continue;
    if(gSM->size() == 0 or gDM->size() == 0 or couplingwgt->size() == 0) continue;

    massPoint a (genMediatorMass,genX1Mass);
    vector<massPoint>::iterator itMass = find(massList.begin(),massList.end(),a);
    sumwgt = (*itMass).sumwgt_;
    treeOut->Fill();    
  }
  cout<<endl;

  treeOut->Write();
  outputFile->Close();
  
}
