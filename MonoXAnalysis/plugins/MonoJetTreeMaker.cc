// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

// FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

// includes
#include "AnalysisCode/MonoXAnalysis/interface/GenParticleTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/JetTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/MuonTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/ElectronTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/PhotonTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/TauTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/TriggerTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/VJetTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/RecoilTreeFiller.h"
#include "AnalysisCode/MonoXAnalysis/interface/JetMetDphiTreeFiller.h"

/// other includes
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"


class MonoJetTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {

public:
  explicit MonoJetTreeMaker(const edm::ParameterSet&);
  ~MonoJetTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  void DeclareAndSetBranches(const edm::ParameterSet&);

  // Gen Particles and MC info
  const bool isMC;
  const bool isReMiniAOD;
  const bool useLHEWeights;

  // xsection value
  float xsec;
  
  // Basic information filled in the main analyzers are pileup, MET filters
  const edm::InputTag pileupInfoTag;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
  
  // InputTags for MET filters filters
  const edm::InputTag filterResultsTag;
  const edm::InputTag badChargedCandidateTag;
  const edm::InputTag badPFMuonTag;
  edm::EDGetTokenT<edm::TriggerResults> filterResultsToken;
  edm::EDGetTokenT<bool> badChargedCandidateToken;
  edm::EDGetTokenT<bool> badPFMuonToken;

  // Vertex
  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  // rho
  const edm::InputTag rhoTag;
  edm::EDGetTokenT<double>  rhoToken;
  // puppi stuff
  const bool addPuppiJets;
  const bool addPuppiMET;
  // Substructure stuff
  const bool addSubstructureCHS;
  const bool addSubstructurePuppi;
  const bool useMiniAODSubstructure;

  // inner vectors
  std::vector<std::string>   filterPathsVector;
  std::map<std::string,int>  filterPathsMap;

  // to store all other objects
  /////--
  TriggerTreeFiller* triggerFiller;
  /////--
  MuonTreeFiller* muonFiller;
  /////--
  ElectronTreeFiller* electronFiller;
  /////--
  PhotonTreeFiller* photonFiller;
  /////--
  TauTreeFiller* tauFiller;
  /////--
  RecoilTreeFiller* recoilTreeFiller;
  /////--
  RecoilTreeFiller* recoilPuppiTreeFiller;
  /////--
  GenParticleTreeFiller* genParticleFiller;
  /////--
  VJetTreeFiller* jetAK8Filler;
  VJetTreeFiller* jetPuppiAK8Filler;
  /////--
  JetTreeFiller* jetAK4Filler;
  JetTreeFiller* jetPuppiAK4Filler;
  ////--
  JetMetDphiTreeFiller* jetMetDphiFiller;
  JetMetDphiTreeFiller* puppiJetMetDphiFiller;
  
  // tree
  TTree* tree;
  edm::Service<TFileService> fs;

  // pileup info
  uint32_t puobs,putrue; 
  // event info
  uint32_t event, run, lumi;  
  uint32_t nvtx; 
  //met filters
  uint8_t flagcsctight,flaghbhenoise,flaghbheiso,flageebadsc,flagecaltp,flaggoodvertices, flagglobaltighthalo, flagbadchpf, flagbadpfmu, flagbadmuon, flagduplicatemuon;  
  float rho, puwgt;
  // for fastSIM
  float samplemedM,sampledmM;
  
};


MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
  isMC          (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  isReMiniAOD   (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  useLHEWeights (iConfig.existsAs<bool>("useLHEWeights") ? iConfig.getParameter<bool>("useLHEWeights") : false),
  xsec          (iConfig.existsAs<double>("xsec") ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.),
  pileupInfoTag (iConfig.getParameter<edm::InputTag>("pileup")),
  filterResultsTag       (iConfig.getParameter<edm::InputTag>("filterResults")),
  badChargedCandidateTag (iConfig.getParameter<edm::InputTag>("badChargedCandidate")),
  badPFMuonTag           (iConfig.getParameter<edm::InputTag>("badPFMuon")),  
  verticesTag   (iConfig.getParameter<edm::InputTag>("vertices")),
  rhoTag        (iConfig.getParameter<edm::InputTag>("rho")),
  addPuppiJets(iConfig.existsAs<bool>("addPuppiJets") ? iConfig.getParameter<bool>("addPuppiJets") : false),
  addPuppiMET(iConfig.existsAs<bool>("addPuppiMET") ? iConfig.getParameter<bool>("addPuppiMET") : false),
  addSubstructureCHS(iConfig.existsAs<bool>("addSubstructureCHS") ? iConfig.getParameter<bool>("addSubstructureCHS") : false), 
  addSubstructurePuppi(iConfig.existsAs<bool>("addSubstructurePuppi") ? iConfig.getParameter<bool>("addSubstructurePuppi") : false),
  useMiniAODSubstructure(iConfig.existsAs<bool>("useMiniAODSubstructure") ? iConfig.getParameter<bool>("useMiniAODSubstructure") : false){

  usesResource();
  usesResource("TFileService");

  if(isMC)
    pileupInfoToken = consumes<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileup"));

  filterResultsToken       = consumes<edm::TriggerResults> (filterResultsTag);
  badChargedCandidateToken = consumes<bool>(badChargedCandidateTag);
  badPFMuonToken           = consumes<bool>(badPFMuonTag);

  verticesToken  = consumes<std::vector<reco::Vertex> > (verticesTag);
  rhoToken       = consumes<double>(rhoTag),
  
  // trigger tokens
  filterResultsToken    = consumes<edm::TriggerResults> (filterResultsTag);
  badChargedCandidateToken = consumes<bool>(badChargedCandidateTag);
  badPFMuonToken           = consumes<bool>(badPFMuonTag);
  verticesToken  = consumes<std::vector<reco::Vertex> > (verticesTag);

  // make output tree
  tree = fs->make<TTree>("tree","tree");  
  // Set Basic Branches
  if(tree!=0){
    DeclareAndSetBranches(iConfig);
  }

  // consumes collector
  edm::ConsumesCollector iC = consumesCollector();

  // initialize all pointers
  triggerFiller  = NULL;
  muonFiller     = NULL;
  electronFiller = NULL;
  photonFiller   = NULL;
  tauFiller      = NULL;
  jetAK4Filler   = NULL;
  jetPuppiAK4Filler = NULL;
  recoilTreeFiller = NULL;
  jetMetDphiFiller = NULL;
  recoilPuppiTreeFiller = NULL;
  puppiJetMetDphiFiller = NULL;
  jetAK8Filler = NULL;
  jetPuppiAK8Filler = NULL;
  genParticleFiller = NULL;

  //Create filler objects
  triggerFiller = new TriggerTreeFiller(iConfig,iC,tree);
  triggerFiller->hltPrescaleProvider.reset(new HLTPrescaleProvider(iConfig,iC,*this)); //ND 

  muonFiller     = new MuonTreeFiller(iConfig,iC,tree);
  electronFiller = new ElectronTreeFiller(iConfig,iC,tree);
  photonFiller   = new PhotonTreeFiller(iConfig,iC,tree);
  tauFiller      = new TauTreeFiller(iConfig,iC,tree);
 

  jetAK4Filler   = new JetTreeFiller(iConfig,iC,tree);
  if(addPuppiJets) 
    jetPuppiAK4Filler  = new JetTreeFiller(iConfig,iC,tree,true);

  recoilTreeFiller = new RecoilTreeFiller(iConfig,iC,tree);
  jetMetDphiFiller = new JetMetDphiTreeFiller(iConfig,iC,tree);
  
  if(addPuppiMET)
    recoilPuppiTreeFiller = new RecoilTreeFiller(iConfig,iC,tree,true);   

  if(addPuppiMET and addPuppiJets)
    puppiJetMetDphiFiller = new JetMetDphiTreeFiller(iConfig,iC,tree,true);

  if(addSubstructureCHS or useMiniAODSubstructure)
    jetAK8Filler   = new VJetTreeFiller(iConfig,iC,tree);

  if(addSubstructurePuppi or useMiniAODSubstructure) 
    jetPuppiAK8Filler  = new VJetTreeFiller(iConfig,iC,tree,true);

  if(isMC)
    genParticleFiller = new GenParticleTreeFiller(iConfig,iC,tree);
}


MonoJetTreeMaker::~MonoJetTreeMaker() {}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace reco;
    using namespace std;
    
    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();
    
    // get the rho value
    Handle<double> rhoH;
    iEvent.getByToken(rhoToken, rhoH);
    rho = *rhoH;

    // MET filters
    Handle<TriggerResults> filterResultsH;
    iEvent.getByToken(filterResultsToken, filterResultsH);
    edm::Handle<bool> filterBadChCandH;
    iEvent.getByToken(badChargedCandidateToken,filterBadChCandH);      
    edm::Handle<bool> filterBadPFMuonH;
    iEvent.getByToken(badPFMuonToken,filterBadPFMuonH);
   
    // MET filter info
    flagbadchpf = *filterBadChCandH;
    flagbadpfmu = *filterBadPFMuonH;
    flagcsctight  = 0; flaghbhenoise = 0; flaghbheiso   = 0;
    flageebadsc   = 0; flagecaltp    = 0; flaggoodvertices = 0;
    flagglobaltighthalo = 0;
    flagbadmuon = 0; flagduplicatemuon = 0;
    
    // Which MET filters passed
    if(filterResultsH.isValid()){
      for (size_t i = 0; i < filterPathsVector.size(); i++) {
        if (filterPathsMap[filterPathsVector[i]] == -1) continue;
	if (filterResultsH->accept(filterPathsMap[filterPathsVector[i]]))
	  if (i == 0  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagcsctight  = 1; // CSCTightHaloFilter
        if (i == 1  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbhenoise = 1; // HBHENoiseFilter
        if (i == 2  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flageebadsc   = 1; // eeBadScFilter
	if (i == 3  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbheiso   = 1; // HBHENoiseIsoFilter
	if (i == 4  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagecaltp    = 1; // Dead Ecal TP
	if (i == 5  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaggoodvertices = 1; // Good vertices
	if (i == 6  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagglobaltighthalo = 1; // 2016 global halo filter	
	if (isReMiniAOD){
	  if (i == 7  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagbadmuon = 1;
	  if (i == 8  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagduplicatemuon = 1;
	}
      }
    }


    // PILEUP INFO    
    Handle<vector<PileupSummaryInfo> > pileupInfoH;    
    if(isMC)
      iEvent.getByToken(pileupInfoToken, pileupInfoH);

    puobs  = 0; putrue = 0; puwgt  = 1.;    
    if (pileupInfoH.isValid()) {
      for (auto pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
	if (pileupInfo_iter->getBunchCrossing() == 0) {
	  puobs  = pileupInfo_iter->getPU_NumInteractions();
	  putrue = pileupInfo_iter->getTrueNumInteractions();
	}
      }
    }

    // VERTEX INFO
    Handle<vector<Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);
    
    if(verticesH.isValid()) nvtx = verticesH->size();
    else nvtx = 0;

    // Fill other object fillers
    bool isGoodEvent = true;
    // fill trigger
    if(triggerFiller != NULL) {
      isGoodEvent = triggerFiller->Fill(iEvent,iSetup);    
      if(not isGoodEvent) return;
    }

    // fill muons
    if(muonFiller != NULL) {
      isGoodEvent = muonFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }

    // fill electrons
    if(electronFiller != NULL) {
      isGoodEvent = electronFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }

    // fill photons
    if(photonFiller != NULL) {
      isGoodEvent = photonFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    // fill taus
    if(tauFiller != NULL) {
      isGoodEvent = tauFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    // fill Recoil
    if(recoilTreeFiller != NULL) {
      isGoodEvent = recoilTreeFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }

    
    // fill Puppi recoil
    if(recoilPuppiTreeFiller  != NULL) {
      isGoodEvent = recoilPuppiTreeFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }

    // fill AK4 jets
    if(jetAK4Filler != NULL) {
      isGoodEvent = jetAK4Filler->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    // fill AK4 Puppi jets
    if(jetPuppiAK4Filler != NULL) {
      isGoodEvent = jetPuppiAK4Filler->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }

    if(jetMetDphiFiller != NULL){
      isGoodEvent = jetMetDphiFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    if(puppiJetMetDphiFiller != NULL){
      isGoodEvent = puppiJetMetDphiFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
       
    // fill AK8 jets
    if(jetAK8Filler != NULL){
      isGoodEvent = jetAK8Filler->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    // fill AK8 Puppi jets
    if(jetPuppiAK8Filler != NULL){
      isGoodEvent = jetPuppiAK8Filler->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    if(genParticleFiller != NULL) {
      isGoodEvent = genParticleFiller->Fill(iEvent,iSetup);
      if(not isGoodEvent) return;
    }
    
    // FILL the Event tree    
    tree->Fill();    
    
}    



void MonoJetTreeMaker::beginJob() {}

void MonoJetTreeMaker::DeclareAndSetBranches(const edm::ParameterSet& iConfig){

  // Run, Lumi, Event info
  tree->Branch("event", &event, "event/i");
  tree->Branch("run"  , &run  , "run/i");
  tree->Branch("lumi" , &lumi , "lumi/i");

  bool isTriggerTree         = (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false);
  bool addPhotonIDVariables  = (iConfig.existsAs<bool>("addPhotonIDVariables") ? iConfig.getParameter<bool>("addPhotonIDVariables") : false);
  bool addElectronIDVariables = (iConfig.existsAs<bool>("addElectronIDVariables") ? iConfig.getParameter<bool>("addElectronIDVariables") : false);
  bool isQCDTree             = (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false);
  bool applyDiMuonFilter     = (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false);
  bool applyDiElectronFilter = (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false);
  bool isPhotonPurity        = (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false);
  
  if(not isTriggerTree){
    tree->Branch("puwgt", &puwgt, "puwgt/F");
    tree->Branch("puobs", &puobs, "puobs/I");
    tree->Branch("putrue", &putrue, "putrue/I");
  }
  
  tree->Branch("xsec"                 , &xsec                 , "xsec/F");  
  tree->Branch("nvtx"                 , &nvtx                 , "nvtx/i");  
  
  // MET filters
  tree->Branch("flagcsctight"         , &flagcsctight         , "flagcsctight/b");
  tree->Branch("flaghbhenoise"        , &flaghbhenoise        , "flaghbhenoise/b");
  tree->Branch("flaghbheiso"          , &flaghbheiso          , "flaghbheiso/b");
  tree->Branch("flageebadsc"          , &flageebadsc          , "flageebadsc/b");
  tree->Branch("flagecaltp"           , &flagecaltp           , "flagecaltp/b");
  tree->Branch("flaggoodvertices"     , &flaggoodvertices     , "flaggoodvertices/b");
  tree->Branch("flagglobaltighthalo"  , &flagglobaltighthalo  , "flagglobaltighthalo/b");
  tree->Branch("flagbadchpf"          , &flagbadchpf          , "flagbadchpf/b");
  tree->Branch("flagbadpfmu"          , &flagbadpfmu          , "flagbadpfmu/b");
  if(isReMiniAOD){
    tree->Branch("flagbadmuon"          , &flagbadmuon        , "flagbadmuon/b");
    tree->Branch("flagduplicatemuon"    , &flagduplicatemuon  , "flagduplicatemuon/b");
  }

  if(isPhotonPurity and not isTriggerTree and not isQCDTree)
    tree->Branch("rho",&rho,"rho/F");
  else if(addPhotonIDVariables and not isTriggerTree and not isQCDTree and not applyDiMuonFilter and not applyDiElectronFilter)
    tree->Branch("rho",&rho,"rho/F");
  else if(addElectronIDVariables and not isTriggerTree and not isQCDTree and not applyDiMuonFilter and not applyDiElectronFilter)
    tree->Branch("rho",&rho,"rho/F");
}

void MonoJetTreeMaker::endJob() {}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  triggerFiller->SetTriggerPaths(iRun,iSetup);

  // MET filter Paths
  filterPathsVector.push_back("Flag_CSCTightHalo2015Filter");
  filterPathsVector.push_back("Flag_HBHENoiseFilter");
  filterPathsVector.push_back("Flag_eeBadScFilter");
  filterPathsVector.push_back("Flag_HBHENoiseIsoFilter");
  filterPathsVector.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filterPathsVector.push_back("Flag_goodVertices");
  filterPathsVector.push_back("Flag_globalTightHalo2016Filter");
  filterPathsVector.push_back("Flag_badMuons");
  filterPathsVector.push_back("Flag_duplicateMuons");
  
  HLTConfigProvider fltrConfig;
  bool changedConfig = false;
  fltrConfig.init(iRun, iSetup, filterResultsTag.process(), changedConfig);
  
  for (size_t i = 0; i < filterPathsVector.size(); i++)
    filterPathsMap[filterPathsVector[i]] = -1;

  for(size_t i = 0; i < filterPathsVector.size(); i++){
    TPRegexp pattern(filterPathsVector[i]);
    for(size_t j = 0; j < fltrConfig.triggerNames().size(); j++){
      std::string pathName = fltrConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	filterPathsMap[filterPathsVector[i]] = j;
      }
    }
  }

  // info about MC cross section in case the xsec parsed has a dummy value
  if(isMC and useLHEWeights)
    genParticleFiller->ReadLHERunProduct(iRun,xsec);
}


void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void MonoJetTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {}

DEFINE_FWK_MODULE(MonoJetTreeMaker);
