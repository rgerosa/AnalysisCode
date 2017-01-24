#include "makeComparisonData.C"

float getNHEA(float eta){
  if(fabs(eta) < 1) return 0.0599;
  if(fabs(eta) > 1 and fabs(eta) < 1.479) return 0.0819;
  else return 1;
}

float getPHEA(float eta){
  if(fabs(eta) < 1) return 0.1271;
  if(fabs(eta) > 1 and fabs(eta) < 1.479) return 0.1101;
  else return 1;
}

void checkPhotonJetsEvents (string inputDir, string eventList, bool isICHEPData){

  gROOT->SetBatch(kTRUE);
  TChain* chain = new TChain("tree/tree","tree/tree");

  if(not isICHEPData)
    chain->Add((inputDir+"/*root").c_str());
  else{
    chain->Add((inputDir+"/*Run2016B*root").c_str());
    chain->Add((inputDir+"/*Run2016C*root").c_str());
    chain->Add((inputDir+"/*Run2016D*root").c_str());
  }

  cout<<"Event in the chain "<<chain->GetEntries()<<endl;
  
  vector<eventID> eventsToSearch;
  string line;
  ifstream inputfile (eventList.c_str());
  if(inputfile.is_open()){
    while(!inputfile.eof()){
      getline(inputfile,line);
      if(line == "") continue;
      stringstream name(line.c_str());
      string segment;
      vector<string> seglist;
      while(getline(name, segment, ':')){
	seglist.push_back(segment);
      }
      if(seglist.size() != 3 ){
	cerr<<"Huston we have a problem with the event specification --> should be run:lumi:event"<<endl;
	continue;
      }
      unsigned int run = stod(seglist.at(0));
      unsigned int lumi = stod(seglist.at(1));
      unsigned int event = stod(seglist.at(2));
      eventsToSearch.push_back(eventID(event,lumi,run));      
    }
  }  
  cout<<"Event to look for are "<<eventsToSearch.size()<<endl;

  TTreeReader reader (chain);

  int foundEvents = 0;
  int passLooseMuonVeto = 0;
  int passLooseElectronVeto = 0;
  int passNPhotons = 0;
  int passKinematicRequirement = 0;
  int passHoverE = 0;
  int passSigmaIetaIeta = 0;
  int passElectronVeto = 0;
  int passPHIsolation = 0;
  int passNHIsolation = 0;
  int passCHIsolation = 0;

  TTreeReaderValue<unsigned int> run         (reader,"run");
  TTreeReaderValue<unsigned int> lumi        (reader,"lumi");
  TTreeReaderValue<unsigned int> event       (reader,"event");
  TTreeReaderValue<float> phpt  (reader,"phPuritypt");
  TTreeReaderValue<float> pheta (reader,"phPurityeta");
  TTreeReaderValue<float> phphi (reader,"phPurityphi");
  TTreeReaderValue<float> phElVeto (reader,"phPurityElectronVeto");
  TTreeReaderValue<float> phPHIso  (reader,"phPHiso");
  TTreeReaderValue<float> phCHIso  (reader,"phCHiso");
  TTreeReaderValue<float> phNHIso  (reader,"phNHiso");
  TTreeReaderValue<float> phHoE    (reader,"phPurityhoe");
  TTreeReaderValue<float> phSieie  (reader,"phPuritysieie");
  TTreeReaderValue<float> phPHIsoRND04 (reader,"phPurityRND04PHiso");
  TTreeReaderValue<float> phPHIsoRND08 (reader,"phPurityRND08PHiso");
  TTreeReaderValue<float> phEAEgamma   (reader,"phPurityEAEGamma");
  TTreeReaderValue<float> rho  (reader,"rho");
  TTreeReaderValue<unsigned int> nphotonsPurity (reader,"nphotonsPurity");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nmuons     (reader,"nmuons");
  TTreeReaderValue<unsigned int> nphotons   (reader,"nphotons");
  TTreeReaderValue<unsigned int> ntausraw   (reader,"ntaus");
  TTreeReaderValue<unsigned int> nbjets     (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> njets  (reader,"njets");
  TTreeReaderValue<UChar_t> hltp165     (reader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175     (reader,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<float> t1met     (reader,"t1pfmet");
  TTreeReaderValue<float> t1metphi  (reader,"t1pfmetphi");

  while(reader.Next()){

    eventID event_tmp(*event,*lumi,*run);
    auto iterator = std::find(eventsToSearch.begin(),eventsToSearch.end(),event_tmp);
    if(iterator != eventsToSearch.end()){
      foundEvents++;
      if(*nmuons != 0) continue;
      passLooseMuonVeto++;
      if(*nelectrons != 0) continue;
      passLooseElectronVeto++;
      if(*phHoE > 0.05) continue;
      passHoverE++;
      if(*phSieie > 0.0102) continue;
      passSigmaIetaIeta++;
      if(*phElVeto != 1) continue;
      passElectronVeto++;
      if(*phCHIso > 1.37) continue;
      passCHIsolation++;
      if(max(0.,double(*phNHIso-*rho*getNHEA(*pheta))) > 1.06*0.014*(*phpt)+0.000019*(*phpt*(*phpt))) continue;
      passNHIsolation++;
      if(max(0.,double(*phPHIso-*rho*getPHEA(*pheta))) > 0.28+0.0053*(*phpt)) continue;
      passPHIsolation++;
      if(*phpt < 175 or fabs(*pheta) > 1.4442) continue;
    }
  }
  cout<<"Found events        = "<<foundEvents<<endl;
  cout<<"Pass muon veto      = "<<passLooseMuonVeto<<endl;
  cout<<"Pass electron veto  = "<<passLooseElectronVeto<<endl;
  cout<<"Pass HoverE         = "<<passHoverE<<endl;
  cout<<"Pass SigIetaIeta    = "<<passSigmaIetaIeta<<endl;
  cout<<"Pass electron veto  = "<<passElectronVeto<<endl;
  cout<<"Pass CH isolation   = "<<passCHIsolation<<endl;
  cout<<"Pass NH isolation   = "<<passNHIsolation<<endl;
  cout<<"Pass PH isolation   = "<<passPHIsolation<<endl;
  cout<<"Pass Nphotons       = "<<passNPhotons<<endl;
}
