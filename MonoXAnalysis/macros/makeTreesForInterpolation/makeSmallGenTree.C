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

//#### selections
static float recoil         = 200;
static float recoilMonoV    = 250;
static float prunedMassMin  = 65.;
static float prunedMassMax  = 105.;
static float tau2tau1       = 0.6;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;
static float detajjVBF      = 3.5;
static float mjjVBF         = 1000;
static float dphijjVBF      = 1.5;
static bool  useMoriondSetup = true;

void makeSmallGenTree(string inputDirectory, string interaction, string signalType, string outputDirectory){

  system(("mkdir -p "+outputDirectory).c_str());

  TChain* chain = new TChain("tree/tree");

  cout<<"Read all the files from inputDirectory "<<inputDirectory<<endl;

  // find all the files inside the base dir
  system(("find "+inputDirectory+" -name \"*root\" > file.temp").c_str());
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

  ////////////
  // load all the weights for the SR selection
  ////////////
  TFile* pufile = NULL;
  if(not useMoriondSetup)
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_12p9fb.root");
  else
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  TH1* puhist = (TH1*) pufile->Get("puhist");

  ////////////
  // trigger efficiency for met trigger                                                                                                                                       
  ////////////
  TFile* triggerfile_MET = NULL;
  vector<TFile*> triggerfile_MET_binned;
  if(useMoriondSetup){
    triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet.root");
    triggerfile_MET_binned.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_0.0_800.0.root"));
    triggerfile_MET_binned.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_800.0_1200.0.root"));
    triggerfile_MET_binned.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1200.0_1700.0.root"));
    triggerfile_MET_binned.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1700.0_3000.0.root"));
  }
  else
    triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/Monojet/metTriggerEfficiency_12p9.root");    
  
  
  TEfficiency*       triggermet       = NULL;
  TGraphAsymmErrors* triggermet_graph = NULL;
  if(triggerfile_MET != NULL){
    triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency_vbf_loose");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency");
    triggermet_graph = triggermet->CreateGraph();
  }

  vector<TF1*> triggermet_func_binned;
  if(triggerfile_MET_binned.size() != 0){
    for(auto ifile : triggerfile_MET_binned)
      triggermet_func_binned.push_back((TF1*) ifile->Get("efficiency_func"));
  }


    
  ////////////
  /// output file
  ////////////

  TFile* outputFile = new TFile((outputDirectory+"/tree_"+interaction+"_"+signalType+".root").c_str(),"RECREATE");
  outputFile->cd();
  TTree* outputTree = new TTree("tree","tree");

  int id;
  float  genVBosonPt, genVBosonEta, genVBosonPhi, genVBosonMass; // vector boson kinematic at generator level if any
  float  genAK8JetPt, genAK8JetEta, genAK8JetPhi, genAK8JetMass, genAK8JetPrunedMass, genAK8JetTau2Tau1; // AK8 jet kinematic at gen level
  float  genLeadAK4JetPt, genLeadAK4JetEta, genLeadAK4JetPhi, genLeadAK4JetMass; // first leading jet gen level
  float  genTrailingAK4JetPt, genTrailingAK4JetEta, genTrailingAK4JetPhi, genTrailingAK4JetMass; // second leading jet at gen level
  float  genMediatorPt, genMediatorEta, genMediatorPhi, genMediatorMass, genMediatorRealMass; // mediator gen level
  float  genMjj,  genDetajj, genDphijj; // generator level info for VBF
  float  genX1Pt, genX1Eta, genX1Phi, genX1Mass, genX1RealMass; // DM candidates
  float  genX2Pt, genX2Eta, genX2Phi, genX2Mass, genX2RealMass; // DM candidates
  float  genMetPt, genMetPhi; // gen-met
  float  weightPU, weightTurnOn, genWeight; // different weights: pileup, trigger-turn on and generator weight
  float  pfMetPt, pfMetPhi; // reco-met and director
  float  recoAK8JetPt,recoAK8JetEta,recoAK8JetPhi,recoAK8JetMass,recoAK8JetPrunedMass, recoAK8JetTau2Tau1; // reco-AK8 info
  float  recoLeadAK4JetPt,recoLeadAK4JetEta,recoLeadAK4JetPhi,recoLeadAK4JetMass; // reco leading jet
  float  recoTrailingAK4JetPt,recoTrailingAK4JetEta,recoTrailingAK4JetPhi,recoTrailingAK4JetMass; // reco leading jet
  float  recoMjj, recoDetajj, recoDphijj;

  outputTree->Branch("id",&id,"id/I");  // id = 0 means event not selections, id = 1 means mono-jet, id = 2 means mono-V, id = 3 means VBF 

  outputTree->Branch("genVBosonPt",   &genVBosonPt,   "genVBosonPt/F");  
  outputTree->Branch("genVBosonEta",  &genVBosonEta,  "genVBosonEta/F");  
  outputTree->Branch("genVBosonPhi",  &genVBosonPhi,  "genVBosonPhi/F");  
  outputTree->Branch("genVBosonMass", &genVBosonMass, "genVBosonMass/F");  

  outputTree->Branch("genAK8JetPt", &genAK8JetPt, "genAK8JetPt/F");  
  outputTree->Branch("genAK8JetEta", &genAK8JetEta, "genAK8JetEta/F");  
  outputTree->Branch("genAK8JetPhi", &genAK8JetPhi, "genAK8JetPhi/F");  
  outputTree->Branch("genAK8JetMass", &genAK8JetMass, "genAK8JetMass/F");  
  outputTree->Branch("genAK8JetPrunedMass", &genAK8JetPrunedMass, "genAK8JetPrunedMass/F");  
  outputTree->Branch("genAK8JetTau2Tau1", &genAK8JetTau2Tau1, "genAK8JetTau2Tau1/F");  

  outputTree->Branch("genLeadAK4JetPt", &genLeadAK4JetPt, "genLeadAK4JetPt/F");  
  outputTree->Branch("genLeadAK4JetEta", &genLeadAK4JetEta, "genLeadAK4JetEta/F");  
  outputTree->Branch("genLeadAK4JetPhi", &genLeadAK4JetPhi, "genLeadAK4JetPhi/F");  
  outputTree->Branch("genLeadAK4JetMass", &genLeadAK4JetMass, "genLeadAK4JetMass/F");  

  outputTree->Branch("genTrailingAK4JetPt", &genTrailingAK4JetPt, "genTrailingAK4JetPt/F");  
  outputTree->Branch("genTrailingAK4JetEta", &genTrailingAK4JetEta, "genTrailingAK4JetEta/F");  
  outputTree->Branch("genTrailingAK4JetPhi", &genTrailingAK4JetPhi, "genTrailingAK4JetPhi/F");  
  outputTree->Branch("genTrailingAK4JetMass", &genTrailingAK4JetMass, "genTrailingAK4JetMass/F");  

  outputTree->Branch("genMjj", &genMjj, "genMjj/F");  
  outputTree->Branch("genDetajj", &genDetajj, "genDetajj/F");  
  outputTree->Branch("genDphijj", &genDphijj, "genDphijj/F");  

  outputTree->Branch("genMediatorPt", &genMediatorPt, "genMediatorPt/F");  
  outputTree->Branch("genMediatorEta", &genMediatorEta, "genMediatorEta/F");  
  outputTree->Branch("genMediatorPhi", &genMediatorPhi, "genMediatorPhi/F");  
  outputTree->Branch("genMediatorMass", &genMediatorMass, "genMediatorMass/F");  
  outputTree->Branch("genMediatorRealMass", &genMediatorRealMass, "genMediatorRealMass/F");  

  outputTree->Branch("genX1Pt", &genX1Pt, "genX1Pt/F");  
  outputTree->Branch("genX1Eta", &genX1Eta, "genX1Eta/F");  
  outputTree->Branch("genX1Phi", &genX1Phi, "genX1Phi/F");  
  outputTree->Branch("genX1Mass", &genX1Mass, "genX1Mass/F");  
  outputTree->Branch("genX1RealMass", &genX1RealMass, "genX1RealMass/F");  

  outputTree->Branch("genX2Pt", &genX2Pt, "genX2Pt/F");  
  outputTree->Branch("genX2Eta", &genX2Eta, "genX2Eta/F");  
  outputTree->Branch("genX2Phi", &genX2Phi, "genX2Phi/F");  
  outputTree->Branch("genX2Mass", &genX2Mass, "genX2Mass/F");  
  outputTree->Branch("genX2RealMass", &genX2RealMass, "genX2RealMass/F");  

  outputTree->Branch("genMetPt", &genMetPt, "genMetPt/F");  
  outputTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");  

  outputTree->Branch("pfMetPt", &pfMetPt, "pfMetPt/F");  
  outputTree->Branch("pfMetPhi", &pfMetPhi, "pfMetPhi/F");  

  outputTree->Branch("recoAK8JetPt", &recoAK8JetPt, "recoAK8JetPt/F");  
  outputTree->Branch("recoAK8JetEta", &recoAK8JetEta, "recoAK8JetEta/F");  
  outputTree->Branch("recoAK8JetPhi", &recoAK8JetPhi, "recoAK8JetPhi/F");  
  outputTree->Branch("recoAK8JetMass", &recoAK8JetMass, "recoAK8JetMass/F");  
  outputTree->Branch("recoAK8JetPrunedMass", &recoAK8JetPrunedMass, "recoAK8JetPrunedMass/F");  
  outputTree->Branch("recoAK8JetTau2Tau1", &recoAK8JetTau2Tau1, "recoAK8JetTau2Tau1/F");  
  
  outputTree->Branch("recoLeadAK4JetPt",&recoLeadAK4JetPt,"recoLeadAK4JetPt/F");
  outputTree->Branch("recoLeadAK4JetEta",&recoLeadAK4JetEta,"recoLeadAK4JetEta/F");
  outputTree->Branch("recoLeadAK4JetPhi",&recoLeadAK4JetPhi,"recoLeadAK4JetPhi/F");
  outputTree->Branch("recoLeadAK4JetMass",&recoLeadAK4JetMass,"recoLeadAK4JetMass/F");

  outputTree->Branch("recoTrailingAK4JetPt",&recoTrailingAK4JetPt,"recoTrailingAK4JetPt/F");
  outputTree->Branch("recoTrailingAK4JetEta",&recoTrailingAK4JetEta,"recoTrailingAK4JetEta/F");
  outputTree->Branch("recoTrailingAK4JetPhi",&recoTrailingAK4JetPhi,"recoTrailingAK4JetPhi/F");
  outputTree->Branch("recoTrailingAK4JetMass",&recoTrailingAK4JetMass,"recoTrailingAK4JetMass/F");

  outputTree->Branch("recoMjj",&recoMjj,"recoMjj/F");
  outputTree->Branch("recoDetajj",&recoDetajj,"recoDetajj/F");
  outputTree->Branch("recoDphijj",&recoDphijj,"recoDphijj/F");

  outputTree->Branch("weightPU", &weightPU, "weightPU/F");  
  outputTree->Branch("weightTurnOn", &weightTurnOn, "weightTurnOn/F");  
  outputTree->Branch("genWeight", &genWeight, "genWeight/F");  
  
  // setup read input file                                                                                                                                                  
  TTreeReader myReader(chain);

  TTreeReaderValue<unsigned int> nvtx (myReader,"nvtx");
  TTreeReaderValue<unsigned int> event (myReader,"event");
  TTreeReaderValue<float>  wgt        (myReader,"wgt");

  TTreeReaderValue<UChar_t> hltm90    (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100   (myReader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110   (myReader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120   (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90  (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm120 (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170 (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300 (myReader,"hltmetwithmu300");

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsct  (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (myReader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (myReader,"flagbadchpf");

  TTreeReaderValue<unsigned int> njets      (myReader,"njets");
  TTreeReaderValue<unsigned int> nphotons   (myReader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (myReader,"ntausrawold");
  TTreeReaderValue<unsigned int> nmuons     (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets     (myReader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nincjets   (myReader,"njetsinc");

  TTreeReaderValue<vector<float> > jetpt   (myReader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<float> > jetm    (myReader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (myReader,"combinejetNHfrac");  

  TTreeReaderValue<vector<float> > jetGenpt   (myReader,"combinejetGenpt");
  TTreeReaderValue<vector<float> > jetGeneta  (myReader,"combinejetGeneta");
  TTreeReaderValue<vector<float> > jetGenphi  (myReader,"combinejetGenphi");
  TTreeReaderValue<vector<float> > jetGenm    (myReader,"combinejetGenm");

  TTreeReaderValue<vector<float> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<float> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<float> > boostedJetGenpt (myReader,"boostedJetGenpt");
  TTreeReaderValue<vector<float> > boostedJetGeneta   (myReader,"boostedJetGeneta");
  TTreeReaderValue<vector<float> > boostedJetGenphi   (myReader,"boostedJetGenphi");
  TTreeReaderValue<vector<float> > boostedJetGenm     (myReader,"boostedJetGenm");
  TTreeReaderValue<vector<float> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<float> > prunedJetGenm   (myReader,"prunedJetGenm");
  TTreeReaderValue<vector<float> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<float> > boostedJetGentau2  (myReader,"boostedJetGentau2");
  TTreeReaderValue<vector<float> > boostedJetGentau1  (myReader,"boostedJetGentau1");

  TTreeReaderValue<float> bosonMass    (myReader,"wzmass_h");
  TTreeReaderValue<float> bosonPt      (myReader,"wzpt_h");
  TTreeReaderValue<float> bosonEta     (myReader,"wzeta_h");
  TTreeReaderValue<float> bosonPhi     (myReader,"wzphi_h");
  TTreeReaderValue<float> mediatorMass    (myReader,"dmmass");
  TTreeReaderValue<float> mediatorPt      (myReader,"dmpt");
  TTreeReaderValue<float> mediatorEta     (myReader,"dmeta");
  TTreeReaderValue<float> mediatorPhi     (myReader,"dmphi");
  TTreeReaderValue<float> x1Mass    (myReader,"dmX1mass");
  TTreeReaderValue<float> x1Pt      (myReader,"dmX1pt");
  TTreeReaderValue<float> x1Eta     (myReader,"dmX1eta");
  TTreeReaderValue<float> x1Phi     (myReader,"dmX1phi");
  TTreeReaderValue<float> x2Mass    (myReader,"dmX2mass");
  TTreeReaderValue<float> x2Pt      (myReader,"dmX2pt");
  TTreeReaderValue<float> x2Eta     (myReader,"dmX2eta");
  TTreeReaderValue<float> x2Phi     (myReader,"dmX2phi");

  TTreeReaderValue<float> mmet       (myReader,"t1mumet");
  TTreeReaderValue<float> mmetphi    (myReader,"t1mumetphi");
  TTreeReaderValue<float> genmet     (myReader,"genmet");
  TTreeReaderValue<float> genmetphi  (myReader,"genmetphi");
  TTreeReaderValue<float> jmmdphi    (myReader,"incjetmumetdphimin4");
 

  // loop on the event and apply selections
  cout<<"Events in the chain "<<chain->GetEntries()<<endl;

  TString fileName;
  string dmMass;
  string medMass;
  vector<string> seglist;
  long int nEvents = 0;
  long int nTotal  = chain->GetEntries();
  while(myReader.Next()){
 
    if(int(nEvents) %10000 == 0){
      std::cout.flush();
      std::cout<<"\r"<<"Events "<<nEvents/nTotal*100<<"%";
    }
    nEvents++;
    // extract real mass value of the sample for the interpolation --> read the string of the tree-name
    TString name_tmp(myReader.GetTree()->GetCurrentFile()->GetName());
    name_tmp.ReplaceAll(inputDirectory.c_str(),"");
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
    
    genVBosonPt   = 0.; genVBosonEta  = 0.; genVBosonPhi  = 0.; genVBosonMass = 0.;
    genAK8JetPt   = 0.; genAK8JetEta  = 0.; genAK8JetPhi  = 0.; genAK8JetMass = 0.; 
    genAK8JetTau2Tau1 = 0.; genAK8JetPrunedMass = 0.;
    genLeadAK4JetPt   = 0.; genLeadAK4JetEta  = 0.; genLeadAK4JetPhi  = 0.; genLeadAK4JetMass = 0.;
    genTrailingAK4JetPt = 0.; genTrailingAK4JetEta = 0.; genTrailingAK4JetPhi = 0.; genTrailingAK4JetMass = 0.;
    genMjj = 0.; genDetajj = 0.; genDphijj = 0.;
    genMediatorPt = 0.; genMediatorEta = 0.; genMediatorPhi  = 0.; genMediatorMass = 0.;
    genMediatorRealMass = 0.;
    genX1Pt  = 0.; genX1Eta = 0.; genX1Phi = 0; genX1Mass = 0.;
    genX2Pt  = 0.; genX2Eta = 0.; genX2Phi = 0.; genX2Mass = 0.;
    genX1RealMass = 0.; genX2RealMass = 0.;
    genMetPt = 0.; genMetPhi = 0.;
    recoAK8JetPt   = 0.; recoAK8JetEta  = 0.; recoAK8JetPhi  = 0.; recoAK8JetMass  = 0.;
    recoAK8JetPrunedMass = 0.; recoAK8JetTau2Tau1   = 0.;
    recoLeadAK4JetPt   = 0.; recoLeadAK4JetEta  = 0.; recoLeadAK4JetPhi  = 0.; recoLeadAK4JetMass  = 0.;
    recoTrailingAK4JetPt   = 0.; recoTrailingAK4JetEta  = 0.; recoTrailingAK4JetPhi  = 0.; recoTrailingAK4JetMass  = 0.;    
    pfMetPt = 0.; pfMetPhi = 0.;
    recoMjj = 0.; recoDetajj = 0.; recoDphijj = 0.;
    weightPU = 1.;
    weightTurnOn = 1.;
    genWeight = 1.;
    
    bool isVBF = false;
    bool isMonoV = false;
    bool isMonoJet = false;
    // basic selections in common between all categories
    // trigger
    if (*hltm90 or *hltm100  or *hltm110  or  *hltm120  or *hltmwm90  or *hltmwm120  or *hltmwm170  or *hltmwm300) {
      if (*fhbhe and *fhbiso  and *feeb  and *fcsct and *fetp  and *fvtx  and *fbadmu  and *fbadch) {
	if (*nbjets  == 0 and *nmuons == 0 and *nelectrons == 0 and *ntaus == 0 and *nphotons == 0){
	  if(jetpt->size() > 0 and *mmet > recoil and *jmmdphi > 0.5){

	    /// check the VBF selection
	    // VBF selections
	    if(jetpt->size()  >= 2 and *nincjets >= 2){
	      if(jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and fabs(jeteta->at(1))< 4.7 and
		 fabs(jeteta->at(0)-jeteta->at(1)) > detajjVBF and jeteta->at(0)*jeteta->at(1) < 0){
		float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
		if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
		if(deltaPhi < dphijjVBF){
		  TLorentzVector jet1, jet2;
		  jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
		  jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
		  if((jet1+jet2).M() > mjjVBF) {
		    if(fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) > 0.1 and nhfrac->at(0) < 0.8)
		      isVBF = true;		    
		    if(fabs(jeteta->at(0)) > 3.0 and fabs(jetpt->at(0)) < 3.2 and nhfrac->at(0) < 0.96)
		      isVBF = true;
		    if(fabs(jeteta->at(0)) > 2.5 and not (fabs(jeteta->at(0)) > 3.0 and fabs(jetpt->at(0)) < 3.2))
		      isVBF = true;
		  }
		}
	      }
	    }
	    
	    // check if mono-V
	    if(not isVBF){      
	      if(jetpt->size()  > 0 and jetpt->at(0)  > 100. and fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) < 0.8 and chfrac->at(0) > 0.1){
		if(boostedJetpt->size()  != 0 and fabs(boostedJeteta->at(0)) < 2.4 and boostedJetpt->at(0) > recoilMonoV and
		   *mmet > recoilMonoV and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1 and
		   prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax){
		  isMonoV = true;
		}
	      }
	    }
	    
	    if(not isVBF and not isMonoV){
	      if(jetpt->size()  > 0 and jetpt->at(0)  > 100. and fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) < 0.8 and chfrac->at(0) > 0.1)
		isMonoJet = true;
	    }

	    // set id
	    if(isVBF)
	      id = 3;
	    else if(isMonoV and not isVBF)
	      id = 2;	    
	    else if(isMonoJet and not isVBF and not isMonoV)
	      id = 1;
	    else if(not isMonoJet and not isVBF and not isMonoV)
	      id = 0;
	    else
	      cerr<<"Problem with nomenclature in the logic of selecting categories --> please check "<<endl;
	  }
	  else id = 0;
	}
	else id = 0;
      }
      else id = 0;
    }else id = 0;
    
    
    if(id == -1) cout<<"Huston we have a problem ... "<<*hltm90<<" "<<*fhbhe<<" "<<*fhbiso<<" "<<*feeb<<" "<<*fcsct<<" "<<*njets<<" "<<*nbjets<<" "<<chfrac->at(0)<<" "<<nhfrac->at(0) << " "<<jetpt->at(0)<<" "<<*jmmdphi<<" "<<*mmet<<" "<<boostedJetpt->at(0)<<" "<<boostedJettau2->at(0)/boostedJettau1->at(0)<<" "<<prunedJetm->at(0)<<endl;
    
    
    // fill branches
    genVBosonPt   = *bosonPt;  
    genVBosonPhi  = *bosonPhi;  
    genVBosonEta  = *bosonEta;  
    genVBosonMass = *bosonMass;  

    if(boostedJetGenpt->size() > 0){
      genAK8JetPt   = boostedJetGenpt->at(0);
      genAK8JetEta  = boostedJetGeneta->at(0);
      genAK8JetPhi  = boostedJetGenphi->at(0);
      genAK8JetMass = boostedJetGenm->at(0);
      genAK8JetPrunedMass = prunedJetGenm->at(0);
      genAK8JetTau2Tau1   = boostedJetGentau2->at(0)/boostedJetGentau1->at(0);
    }
    
    if(jetGenpt->size() > 0){
      genLeadAK4JetPt   = jetGenpt->at(0);
      genLeadAK4JetEta  = jetGeneta->at(0);
      genLeadAK4JetPhi  = jetGenphi->at(0);
      genLeadAK4JetMass = jetGenm->at(0);
    }
    
    if(jetGenpt->size() > 1){
      genTrailingAK4JetPt   = jetGenpt->at(1);
      genTrailingAK4JetEta  = jetGeneta->at(1);
      genTrailingAK4JetPhi  = jetGenphi->at(1);
      genTrailingAK4JetMass = jetGenm->at(1);

      genDetajj = fabs(genLeadAK4JetEta-genTrailingAK4JetEta);
      genDphijj = fabs(genLeadAK4JetPhi-genTrailingAK4JetPhi);
      if(genDphijj > TMath::Pi()) genDphijj = 2*TMath::Pi()-genDphijj;
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(genLeadAK4JetPt,genLeadAK4JetEta,genLeadAK4JetPhi,genLeadAK4JetMass);
      jet2.SetPtEtaPhiM(genTrailingAK4JetPt,genTrailingAK4JetEta,genTrailingAK4JetPhi,genTrailingAK4JetMass);
      genMjj = (jet1+jet2).M();
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

    genMediatorPt   = *mediatorPt;
    genMediatorEta  = *mediatorEta;
    genMediatorPhi  = *mediatorPhi;
    genMediatorRealMass = *mediatorMass;
    genMediatorMass = stod(medMass);

    pfMetPt   = *mmet;
    pfMetPhi  = *mmetphi;
    genMetPt  = *genmet;
    genMetPhi = *genmetphi;

    if(boostedJetpt->size() > 0){
      recoAK8JetPt         = boostedJetpt->at(0);
      recoAK8JetEta        = boostedJeteta->at(0);
      recoAK8JetPhi        = boostedJetphi->at(0);
      recoAK8JetMass       = boostedJetm->at(0);
      recoAK8JetPrunedMass = prunedJetm->at(0);
      recoAK8JetTau2Tau1   = boostedJettau2->at(0)/boostedJettau1->at(0);
    }

    if(jetpt->size() > 0){
      recoLeadAK4JetPt         = jetpt->at(0);
      recoLeadAK4JetEta        = jeteta->at(0);
      recoLeadAK4JetPhi        = jetphi->at(0);
      recoLeadAK4JetMass       = jetm->at(0);
    }

    if(jetpt->size() > 1){
      recoTrailingAK4JetPt         = jetpt->at(1);
      recoTrailingAK4JetEta        = jeteta->at(1);
      recoTrailingAK4JetPhi        = jetphi->at(1);
      recoTrailingAK4JetMass       = jetm->at(1);

      recoDetajj = fabs(recoLeadAK4JetEta-recoTrailingAK4JetEta);
      recoDphijj = fabs(recoLeadAK4JetPhi-recoTrailingAK4JetPhi);
      if(recoDphijj > TMath::Pi()) recoDphijj = 2*TMath::Pi()-recoDphijj;
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(recoLeadAK4JetPt,recoLeadAK4JetEta,recoLeadAK4JetPhi,recoLeadAK4JetMass);
      jet2.SetPtEtaPhiM(recoTrailingAK4JetPt,recoTrailingAK4JetEta,recoTrailingAK4JetPhi,recoTrailingAK4JetMass);
      recoMjj = (jet1+jet2).M();
      
    }

    if(isVBF and useMoriondSetup and jetpt->size() > 1){
      TLorentzVector jet1 ;
      TLorentzVector jet2 ;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      if((jet1+jet2).M() < 800)
	weightTurnOn = triggermet_func_binned.at(0)->Eval(min(double(*mmet),triggermet_func_binned.at(0)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	weightTurnOn = triggermet_func_binned.at(1)->Eval(min(double(*mmet),triggermet_func_binned.at(1)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	weightTurnOn = triggermet_func_binned.at(2)->Eval(min(double(*mmet),triggermet_func_binned.at(2)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1700)
	weightTurnOn = triggermet_func_binned.at(3)->Eval(min(double(*mmet),triggermet_func_binned.at(3)->GetXaxis()->GetXmax()));
    }
    else{
      weightTurnOn = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));
    }

    weightPU     = puhist->GetBinContent(puhist->FindBin(*nvtx));    
    genWeight = *wgt;

    outputTree->Fill();       
  }

  std::cout<<std::endl;

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}
