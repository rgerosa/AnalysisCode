#include <string>
#include <vector>
#include <sstream>

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLorentzVector.h"

void makeSmallGenTree(string inputDirectory, string interaction, string signalType, string outputDirectory){

  system(("mkdir -p "+outputDirectory).c_str());

  TChain* chain = new TChain("tree/tree");

  string baseDir = "";
  if(interaction == "Vector" and signalType == "MonoJ")
    baseDir = inputDirectory+"DMV_Vector/";
  else if(interaction == "Axial" and signalType == "MonoJ")
    baseDir = inputDirectory+"DMV_Axial/";
  else if(interaction == "Scalar" and signalType == "MonoJ")
    baseDir = inputDirectory+"DMS_Scalar/";
  else if(interaction == "Pseudoscalar" and signalType == "MonoJ")
    baseDir = inputDirectory+"DMS_Pseudoscalar/";
  else if(interaction == "Vector" and signalType == "MonoW")
    baseDir = inputDirectory+"MonoW_Vector/";
  else if(interaction == "Axial" and signalType == "MonoW")
    baseDir = inputDirectory+"MonoW_Axial/";
  else if(interaction == "Scalar" and signalType == "MonoW")
    baseDir = inputDirectory+"MonoW_Scalar/";
  else if(interaction == "Pseudoscalar" and signalType == "MonoW")
    baseDir = inputDirectory+"MonoW_Pseudoscalar/";
  else if(interaction == "Vector" and signalType == "MonoW")
    baseDir = inputDirectory+"MonoZ_Vector/";
  else if(interaction == "Axial" and signalType == "MonoZ")
    baseDir = inputDirectory+"MonoZ_Axial/";
  else if(interaction == "Scalar" and signalType == "MonoZ")
    baseDir = inputDirectory+"MonoZ_Scalar/";
  else if(interaction == "Pseudoscalar" and signalType == "MonoZ")
    baseDir = inputDirectory+"MonoZ_Pseudoscalar/";  
  else if(interaction == "HiggsInv" and signalType == "ggH")
    baseDir = inputDirectory+"ggH/";  
  else if(interaction == "HiggsInv" and signalType == "VBF")
    baseDir = inputDirectory+"qqH/";  

  cout<<"Read all the files from baseDir "<<baseDir<<endl;

  // find all the files inside the base dir
  system(("find "+baseDir+" -name \"*root\" > file.temp").c_str());
  ifstream infile;
  string line;
  infile.open("file.temp");
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      chain->Add(line.c_str());
      cout<<"Add following file into the chain "<<line<<endl;
    }
  }
  infile.close();
  system("rm file.temp");

  // load all the re-weight files
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_12p9fb.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");

  // trigger efficiency for met trigger                                                                                                                                       
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency_12p9.root");
  TEfficiency*  triggermet_eff    = (TEfficiency*) trmfile->Get("trig_eff");
  TGraphAsymmErrors* triggermet = (TGraphAsymmErrors*) triggermet_eff->CreateGraph();

  TFile* outputFile = new TFile((outputDirectory+"/tree_"+interaction+"_"+signalType+".root").c_str(),"RECREATE");
  outputFile->cd();

  TTree* outputTree = new TTree("tree","tree");

  int id;
  double  genVBosonPt, genVBosonEta, genVBosonPhi, genVBosonMass;
  double  genVLepBosonPt, genVLepBosonEta, genVLepBosonPhi, genVLepBosonMass;
  double  genVHadBosonPt, genVHadBosonEta, genVHadBosonPhi, genVHadBosonMass;
  double  genAK8JetPt, genAK8JetEta, genAK8JetPhi, genAK8JetMass, genAK8JetPrunedMass, genAK8JetTau2Tau1;
  double  genAK4JetPt, genAK4JetEta, genAK4JetPhi, genAK4JetMass;
  double  genMediatorPt, genMediatorEta, genMediatorPhi, genMediatorMass, genMediatorRealMass;
  double  genX1Pt, genX1Eta, genX1Phi, genX1Mass, genX1RealMass;
  double  genX2Pt, genX2Eta, genX2Phi, genX2Mass, genX2RealMass;
  double  genMetPt, genMetPhi;
  double  weightPU, weightTurnOn, genWeight;
  double  pfMetPt, pfMetPhi;
  double  recoAK8JetPt,recoAK8JetPrunedMass, recoAK8JetTau2Tau1;

  outputTree->Branch("id",&id,"id/I");  
  outputTree->Branch("genVBosonPt", &genVBosonPt, "genVBosonPt/D");  
  outputTree->Branch("genVBosonEta", &genVBosonEta, "genVBosonEta/D");  
  outputTree->Branch("genVBosonPhi", &genVBosonPhi, "genVBosonPhi/D");  
  outputTree->Branch("genVBosonMass", &genVBosonMass, "genVBosonMass/D");  

  outputTree->Branch("genVLepBosonPt", &genVLepBosonPt, "genVLepBosonPt/D");  
  outputTree->Branch("genVLepBosonEta", &genVLepBosonEta, "genVLepBosonEta/D");  
  outputTree->Branch("genVLepBosonPhi", &genVLepBosonPhi, "genVLepBosonPhi/D");  
  outputTree->Branch("genVLepBosonMass", &genVLepBosonMass, "genVLepBosonMass/D");  

  outputTree->Branch("genVHadBosonPt", &genVHadBosonPt, "genVHadBosonPt/D");  
  outputTree->Branch("genVHadBosonEta", &genVHadBosonEta, "genVHadBosonEta/D");  
  outputTree->Branch("genVHadBosonPhi", &genVHadBosonPhi, "genVHadBosonPhi/D");  
  outputTree->Branch("genVHadBosonMass", &genVHadBosonMass, "genVHadBosonMass/D");  

  outputTree->Branch("recoAK8JetPt", &recoAK8JetPt, "recoAK8JetPt/D");  
  outputTree->Branch("recoAK8JetPrunedMass", &recoAK8JetPrunedMass, "recoAK8JetPrunedMass/D");  
  outputTree->Branch("recoAK8JetTau2Tau1", &recoAK8JetTau2Tau1, "recoAK8JetTau2Tau1/D");  

  outputTree->Branch("genAK8JetPt", &genAK8JetPt, "genAK8JetPt/D");  
  outputTree->Branch("genAK8JetEta", &genAK8JetEta, "genAK8JetEta/D");  
  outputTree->Branch("genAK8JetPhi", &genAK8JetPhi, "genAK8JetPhi/D");  
  outputTree->Branch("genAK8JetMass", &genAK8JetMass, "genAK8JetMass/D");  
  outputTree->Branch("genAK8JetPrunedMass", &genAK8JetPrunedMass, "genAK8JetPrunedMass/D");  
  outputTree->Branch("genAK8JetTau2Tau1", &genAK8JetTau2Tau1, "genAK8JetTau2Tau1/D");  
 
  outputTree->Branch("genAK4JetPt", &genAK4JetPt, "genAK4JetPt/D");  
  outputTree->Branch("genAK4JetEta", &genAK4JetEta, "genAK4JetEta/D");  
  outputTree->Branch("genAK4JetPhi", &genAK4JetPhi, "genAK4JetPhi/D");  
  outputTree->Branch("genAK4JetMass", &genAK4JetMass, "genAK4JetMass/D");  

  outputTree->Branch("genMediatorPt", &genMediatorPt, "genMediatorPt/D");  
  outputTree->Branch("genMediatorEta", &genMediatorEta, "genMediatorEta/D");  
  outputTree->Branch("genMediatorPhi", &genMediatorPhi, "genMediatorPhi/D");  
  outputTree->Branch("genMediatorMass", &genMediatorMass, "genMediatorMass/D");  
  outputTree->Branch("genMediatorRealMass", &genMediatorRealMass, "genMediatorRealMass/D");  

  outputTree->Branch("genX1Pt", &genX1Pt, "genX1Pt/D");  
  outputTree->Branch("genX1Eta", &genX1Eta, "genX1Eta/D");  
  outputTree->Branch("genX1Phi", &genX1Phi, "genX1Phi/D");  
  outputTree->Branch("genX1Mass", &genX1Mass, "genX1Mass/D");  
  outputTree->Branch("genX1RealMass", &genX1RealMass, "genX1RealMass/D");  

  outputTree->Branch("genX2Pt", &genX2Pt, "genX2Pt/D");  
  outputTree->Branch("genX2Eta", &genX2Eta, "genX2Eta/D");  
  outputTree->Branch("genX2Phi", &genX2Phi, "genX2Phi/D");  
  outputTree->Branch("genX2Mass", &genX2Mass, "genX2Mass/D");  
  outputTree->Branch("genX2RealMass", &genX2RealMass, "genX2RealMass/D");  

  outputTree->Branch("genMetPt", &genMetPt, "genMetPt/D");  
  outputTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/D");  
  outputTree->Branch("pfMetPt", &pfMetPt, "pfMetPt/D");  
  outputTree->Branch("pfMetPhi", &pfMetPhi, "pfMetPhi/D");  

  outputTree->Branch("weightPU", &weightPU, "weightPU/D");  
  outputTree->Branch("weightTurnOn", &weightTurnOn, "weightTurnOn/D");  
  outputTree->Branch("genWeight", &genWeight, "genWeight/D");  
  
  // setup read input file                                                                                                                                                  
  TTreeReader myReader(chain);

  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
  TTreeReaderValue<double> wgt   (myReader,"wgt");

  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nphotons  (myReader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons  (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus  (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons  (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets (myReader,"nbjetslowpt");
  TTreeReaderValue<vector<double> > jetpt   (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<double> > jetm    (myReader,"combinejetm");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"combinejetNHfrac");  
  TTreeReaderValue<vector<double> > jetGenpt   (myReader,"combinejetGenpt");
  TTreeReaderValue<vector<double> > jetGeneta  (myReader,"combinejetGeneta");
  TTreeReaderValue<vector<double> > jetGenphi  (myReader,"combinejetGenphi");
  TTreeReaderValue<vector<double> > jetGenm    (myReader,"combinejetGenm");
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > boostedJetGenpt (myReader,"boostedJetGenpt");
  TTreeReaderValue<vector<double> > boostedJetGeneta   (myReader,"boostedJetGeneta");
  TTreeReaderValue<vector<double> > boostedJetGenphi   (myReader,"boostedJetGenphi");
  TTreeReaderValue<vector<double> > boostedJetGenm     (myReader,"boostedJetGenm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > prunedJetGenm      (myReader,"prunedJetGenm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<double> > boostedJetGentau2  (myReader,"boostedJetGentau2");
  TTreeReaderValue<vector<double> > boostedJetGentau1  (myReader,"boostedJetGentau1");
  TTreeReaderValue<double > bosonMass    (myReader,"wzmass_h");
  TTreeReaderValue<double > bosonPt      (myReader,"wzpt_h");
  TTreeReaderValue<double > bosonEta     (myReader,"wzeta_h");
  TTreeReaderValue<double > bosonPhi     (myReader,"wzphi_h");
  TTreeReaderValue<double > bosonMass_lep    (myReader,"wzmass");
  TTreeReaderValue<double > bosonPt_lep      (myReader,"wzpt");
  TTreeReaderValue<double > bosonEta_lep     (myReader,"wzeta");
  TTreeReaderValue<double > bosonPhi_lep     (myReader,"wzphi");
  TTreeReaderValue<double > mediatorMass    (myReader,"dmmass");
  TTreeReaderValue<double > mediatorPt      (myReader,"dmpt");
  TTreeReaderValue<double > mediatorEta     (myReader,"dmeta");
  TTreeReaderValue<double > mediatorPhi     (myReader,"dmphi");
  TTreeReaderValue<double > x1Mass    (myReader,"dmX1mass");
  TTreeReaderValue<double > x1Pt      (myReader,"dmX1pt");
  TTreeReaderValue<double > x1Eta     (myReader,"dmX1eta");
  TTreeReaderValue<double > x1Phi     (myReader,"dmX1phi");
  TTreeReaderValue<double > x2Mass    (myReader,"dmX2mass");
  TTreeReaderValue<double > x2Pt      (myReader,"dmX2pt");
  TTreeReaderValue<double > x2Eta     (myReader,"dmX2eta");
  TTreeReaderValue<double > x2Phi     (myReader,"dmX2phi");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double> genmet        (myReader,"genmet");
  TTreeReaderValue<double> genmetphi     (myReader,"genmetphi");
  TTreeReaderValue<double> jmmdphi (myReader,"incjetmumetdphimin4");

  // loop on the event and apply selections
  cout<<"Events in the chain "<<chain->GetEntries()<<endl;

  TString fileName;
  string dmMass;
  string medMass;
  vector<string> seglist;
  
  while(myReader.Next()){
    
    TString name_tmp(myReader.GetTree()->GetCurrentFile()->GetName());
    name_tmp.ReplaceAll(baseDir.c_str(),"");
    if(fileName != name_tmp){
      seglist.clear();
      fileName = name_tmp;
      stringstream name(fileName.Data());
      string segment;
    
      while(getline(name, segment, '/')){
	seglist.push_back(segment);
      }
      
      TString name_tmp_2 (seglist.at(0).c_str());
      name_tmp_2.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-powheg","");
      name_tmp_2.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-powheg","");
      name_tmp_2.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-madgraph","");
      name_tmp_2.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-JHUGen","");
      name_tmp_2.ReplaceAll("_gSM-0p25_gDM-1p0_13TeV-powheg","");
      name_tmp_2.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-powheg","");
      name_tmp_2.ReplaceAll("_gSM-0p25_gDM-1p0_13TeV-madgraph","");
      name_tmp_2.ReplaceAll("_gSM-0p25_gDM-1p0_13TeV-JHUGen","");
      name_tmp_2.ReplaceAll("_13TeV_powheg_pythia8","");
      vector<string> seglist_2;
      string segment_2;
      stringstream name_2(name_tmp_2.Data());
      if(interaction != "HiggsInv"){
	while(getline(name_2, segment_2, '-')){
	  seglist_2.push_back(segment_2);
	}
	dmMass = seglist_2.back();
	TString tmp (seglist_2.at(seglist_2.size()-2).c_str());
	tmp.ReplaceAll("_Mchi","");
	medMass = tmp;
      }
      else{
	while(getline(name_2, segment_2, '_')){
          seglist_2.push_back(string(TString(segment_2).ReplaceAll("M","").Data()));
        }
        medMass = seglist_2.back();
        dmMass  = "1";
      }
    }
    

    // Set Branches
    id = -1;
    
    genVBosonPt   = 0.;
    genVBosonEta  = 0.;
    genVBosonPhi  = 0.;
    genVBosonMass = 0.;

    genVLepBosonPt   = 0.;
    genVLepBosonEta  = 0.;
    genVLepBosonPhi  = 0.;
    genVLepBosonMass = 0.;

    genVHadBosonPt   = 0.;
    genVHadBosonEta  = 0.;
    genVHadBosonPhi  = 0.;
    genVHadBosonMass = 0.;

    recoAK8JetPt   = 0.;
    genAK8JetPt   = 0.;
    genAK8JetEta  = 0.; 
    genAK8JetPhi  = 0.; 
    genAK8JetMass = 0.; 
    genAK8JetPrunedMass = 0.;
    recoAK8JetPrunedMass = 0.;
    genAK8JetTau2Tau1   = 0.;
    recoAK8JetTau2Tau1   = 0.;
    genAK4JetPt     = 0.; 
    genAK4JetEta    = 0.;
    genAK4JetPhi    = 0.; 
    genAK4JetMass   = 0.;
    genMediatorPt   = 0.;
    genMediatorEta  = 0.; 
    genMediatorPhi  = 0.;
    genMediatorMass = 0.;
    genMediatorRealMass = 0.;
    genX1Pt  = 0.;
    genX1Eta = 0.;
    genX1Phi = 0.;
    genX1Mass = 0.;
    genX1RealMass = 0.;
    genX2Pt  = 0.; 
    genX2Eta = 0.; 
    genX2Phi = 0.; 
    genX2Mass = 0.;
    genX2RealMass = 0.;
    genMetPt = 0.; 
    genMetPhi = 0.;
    pfMetPt = 0.;
    pfMetPhi = 0.;
    weightPU = 1.;
    weightTurnOn = 1.;
    genWeight = 1.;

    // basic monojet
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) id = 0;
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) id = 0;
    if (*njets  < 1) id = 0;
    if (*nbjets > 0) id = 0;    
    if (chfrac->size() == 0 or nhfrac->size() == 0 or jetpt->size() == 0 or jeteta->size() == 0) id = 0;                                                                
    if (chfrac->size() > 0 and chfrac->at(0) < 0.1)   id = 0;
    if (nhfrac->size() > 0 and nhfrac->at(0) > 0.8)   id = 0;
    if (jeteta->size() > 0 and fabs(jeteta->at(0)) > 2.5) id = 0;
    if (jetpt->size()  > 0 and jetpt->at(0)  < 100.)  id = 0;
    if (*jmmdphi < 0.5)  id = 0;
    if (*mmet < 200)     id = 0;    
    if (*nmuons > 0)     id = 0;
    if (*nelectrons > 0) id = 0;
    if (*ntaus > 0)      id = 0;
    if (*nphotons > 0)   id = 0;

    if(id == -1 and boostedJetpt->size()  == 0) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and fabs(boostedJeteta->at(0)) > 2.4) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and boostedJetpt->at(0) < 250.) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and (prunedJetm->at(0) < 65 or prunedJetm->at(0) > 105)) id = 1;

    if(id == -1 and prunedJetm->size() > 0 and prunedJetm->at(0) > 0 and boostedJetpt->at(0) > 250. and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105 and *mmet > 250 and fabs(boostedJeteta->at(0)) < 2.4) id = 2;
    
    if (id == -1 and *mmet > 200 and *mmet < 250) id = 0;

    if(id == -1) cout<<"Huston we have a problem ... "<<*hltm90<<" "<<*fhbhe<<" "<<*fhbiso<<" "<<*feeb<<" "<<*fcsc<<" "<<*njets<<" "<<*nbjets<<" "<<chfrac->at(0)<<" "<<nhfrac->at(0) << " "<<jetpt->at(0)<<" "<<*jmmdphi<<" "<<*mmet<<" "<<boostedJetpt->at(0)<<" "<<boostedJettau2->at(0)/boostedJettau1->at(0)<<" "<<prunedJetm->at(0)<<endl;
      
    // fill branches
    genVBosonPt   = *bosonPt;  
    genVBosonPhi  = *bosonPhi;  
    genVBosonEta  = *bosonEta;  
    genVBosonMass = *bosonMass;  
    
    if(genVBosonPt <= 0){
      genVBosonPt   = *bosonPt_lep;
      genVBosonPhi  = *bosonPhi_lep;
      genVBosonEta  = *bosonEta_lep;
      genVBosonMass = *bosonMass_lep;
    }

    genVLepBosonPt   = *bosonPt_lep;
    genVLepBosonPhi  = *bosonPhi_lep;
    genVLepBosonEta  = *bosonEta_lep;
    genVLepBosonMass = *bosonMass_lep;

    genVHadBosonPt   = *bosonPt;
    genVHadBosonPhi  = *bosonPhi;
    genVHadBosonEta  = *bosonEta;
    genVHadBosonMass = *bosonMass;

    if(boostedJetGenpt->size() > 0){
      genAK8JetPt   = boostedJetGenpt->at(0);
      genAK8JetEta  = boostedJetGeneta->at(0);
      genAK8JetPhi  = boostedJetGenphi->at(0);
      genAK8JetMass = boostedJetGenm->at(0);
      genAK8JetPrunedMass = prunedJetGenm->at(0);
      genAK8JetTau2Tau1   = boostedJetGentau2->at(0)/boostedJetGentau1->at(0);
    }
    
    if(jetGenpt->size() > 0){
      genAK4JetPt   = jetGenpt->at(0);
      genAK4JetEta  = jetGeneta->at(0);
      genAK4JetPhi  = jetGenphi->at(0);
      genAK4JetMass = jetGenm->at(0);
    }
    
    genX1Pt   = *x1Pt;
    genX1Eta  = *x1Eta;
    genX1Phi  = *x1Phi;
    genX1RealMass = *x1Mass;
    genX1Mass = stod(dmMass);
    genX2Pt   = *x2Pt;
    genX2Eta  = *x2Eta;
    genX2Phi  = *x2Phi;
    genX2RealMass = *x2Mass;
    genX2Mass = stod(dmMass);

    TLorentzVector dmX1, dmX2, dmMED;
    dmX1.SetPtEtaPhiM(genX1Pt,genX1Eta,genX1Phi,genX1Mass);
    dmX2.SetPtEtaPhiM(genX2Pt,genX2Eta,genX2Phi,genX2Mass);
    dmMED = dmX1+dmX2;

    genMediatorPt   = *mediatorPt;
    genMediatorEta  = *mediatorEta;
    genMediatorPhi  = *mediatorPhi;
    genMediatorRealMass = *mediatorMass;
    genMediatorMass = stod(medMass);

    weightTurnOn = triggermet->Eval(min(*mmet,triggermet->GetXaxis()->GetXmax()));
    if (*nvtx <= 60)
      weightPU = puhist->GetBinContent(puhist->FindBin(*nvtx));
    
    pfMetPt   = *mmet;
    pfMetPhi  = *mmetphi;
    genMetPt  = *genmet;
    genMetPhi = *genmetphi;

    if(boostedJetpt->size() > 0){
      recoAK8JetPt         = boostedJetpt->at(0);
      recoAK8JetPrunedMass = prunedJetm->at(0);
      recoAK8JetTau2Tau1   = boostedJettau2->at(0)/boostedJettau1->at(0);
    }

    genWeight = *wgt;
    outputTree->Fill();       
  }
  
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}
