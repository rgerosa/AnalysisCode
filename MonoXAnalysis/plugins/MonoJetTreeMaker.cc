// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

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

// HLT info
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// Gen Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

// DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Jet corrections
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// b-tagging SF
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"

// Level-1
#include "DataFormats/L1Trigger/interface/EGamma.h" 
#include "DataFormats/L1Trigger/interface/Tau.h"    
#include "DataFormats/L1Trigger/interface/Jet.h"    
#include "DataFormats/L1Trigger/interface/Muon.h"   
#include "DataFormats/L1Trigger/interface/EtSum.h"  
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"

// ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

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
  
  // find photon info
  void findMother(const reco::Candidate*, int &, double &, double &, double &);
  void findFirstNonPhotonMother(const reco::Candidate*, int &, double &, double &, double &);
  // compute isolation
  double computeMuonIso(const reco::Muon &);  
  // to apply jet ID
  bool applyJetID(const pat::Jet &, const std::string &);
  // to apply pileup-jet id
  bool applyPileupJetID(const pat::Jet &, const std::string &, const bool &);
  // fill btag scale factors
  void calculateBtagSF(const pat::Jet &, const std::string &, std::vector<double> &, std::vector<double> &, std::vector<double> &);

  // information for photon purity
  double getGammaEAForPhotonIso(double eta);
  double getGammaNewEAForPhotonIso(double eta);
  double getChargedHadronEAForPhotonIso(double eta);
  double getNeutralHadronEAForPhotonIso(double eta);
  double computeDR(const reco::Candidate *genPart,pat::PhotonRef phot);

  // to dump trigger informations (HLT,L1)
  void fillTriggerObjects(const edm::Handle<pat::TriggerObjectStandAloneCollection> & triggerObjectsH, const edm::TriggerNames & trignames);
  void fillTriggerL1(const edm::Handle<l1t::EGammaBxCollection> & H_L1EG,  const edm::Handle<l1t::TauBxCollection>  & H_L1Tau,
		     const edm::Handle<l1t::JetBxCollection>    & H_L1Jet, const edm::Handle<l1t::MuonBxCollection> & H_L1Mu,
		     const edm::Handle<l1t::EtSumBxCollection>  & H_L1Sums);
  void fillAlgosL1(const edm::Event& iEvent, const edm::EventSetup & eventSetup, const edm::Handle<GlobalAlgBlkBxCollection> & H_L1Algos);


  bool readDMFromGenParticle;

  // Gen Particles and MC info
  const bool isMC;
  const bool isTriggerTree;
  const bool addTriggerObjects;
  const bool uselheweights;
  const edm::InputTag lheEventTag;
  const edm::InputTag lheRunTag;
  const bool isSignalSample;
  const bool addGenParticles;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
  edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
  edm::EDGetTokenT<LHEEventProduct>                  lheInfoToken;
  edm::EDGetTokenT<LHERunInfoProduct>                lheRunInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
  double xsec;

  // InputTags for triggers and met filters
  const edm::InputTag triggerResultsTag;
  const edm::InputTag filterResultsTag;
  const edm::InputTag prescalesTag;
  const edm::InputTag badChargedCandidateTag;
  const edm::InputTag badPFMuonTag;
  const edm::InputTag triggerObjectsTag; 
  const edm::InputTag IT_L1_EG;  
  const edm::InputTag IT_L1_Jet; 
  const edm::InputTag IT_L1_Mu;  
  const edm::InputTag IT_L1_Sums;
  const edm::InputTag IT_L1_Algos; 

  // trgger and filter tokens
  edm::EDGetTokenT<edm::TriggerResults>                    triggerResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescalesToken;
  edm::EDGetTokenT<edm::TriggerResults>                    filterResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken; 
  edm::EDGetTokenT<l1t::EGammaBxCollection>  T_L1EG;
  edm::EDGetTokenT<l1t::JetBxCollection>     T_L1Jet;
  edm::EDGetTokenT<l1t::MuonBxCollection>    T_L1Mu;
  edm::EDGetTokenT<l1t::EtSumBxCollection>   T_L1Sums;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> T_L1Algos;
  edm::EDGetTokenT<bool> badChargedCandidateToken;
  edm::EDGetTokenT<bool> badPFMuonToken;
  
  // Vertex
  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  // rho
  const edm::EDGetTokenT<double>  rhoToken;
  
  // muons
  const edm::InputTag muonsTag;
  const edm::InputTag tightmuonsTag;
  const edm::InputTag highptmuonsTag;

  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> tightmuonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> highptmuonsToken;

  // electrons
  const edm::InputTag  electronsTag;
  const edm::InputTag  looseelectronsTag;
  const edm::InputTag  tightelectronsTag;
  const edm::InputTag  heepelectronsTag;

  edm::EDGetTokenT<pat::ElectronRefVector>  electronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  looseelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  tightelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  heepelectronsToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    electronLooseIdToken;  

  // Photons
  const edm::InputTag  photonsTag;
  const edm::InputTag  mediumphotonsTag;
  const edm::InputTag  tightphotonsTag;
  const edm::InputTag  photonLooseIdTag;
  const edm::InputTag  photonHighPtIdTag;
  const bool           addPhotonPurity;

  edm::EDGetTokenT<pat::PhotonRefVector>    photonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    mediumphotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    tightphotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    photonsPurityToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    tightphotonsPurityToken;  
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonLooseIdToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonHighPtIdToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   photonsieieToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   photonPHisoToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   photonCHisoToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   photonNHisoToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   rndgammaiso04Token;
  edm::EDGetTokenT<edm::ValueMap<double> >   rndgammaiso08Token;
  edm::EDGetTokenT<edm::ValueMap<double> >   gammaisoToken;
  edm::EDGetTokenT<edm::ValueMap<double> >   rndchhadiso04Token;
  edm::EDGetTokenT<edm::ValueMap<double> >   rndchhadiso08Token;

  // Taus
  const edm::InputTag tausNewTag;
  const edm::InputTag tausOldTag;
  const edm::InputTag tausNewRawTag;
  const edm::InputTag tausOldRawTag;

  edm::EDGetTokenT<pat::TauRefVector>  tausNewToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausOldToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausNewRawToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausOldRawToken;

  //Jets AK4
  const edm::InputTag jetsTag;
  const bool addPuppiJets;
  const edm::InputTag puppijetsTag;

  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  puppijetsToken;

  // MET
  const edm::InputTag t1metTag;
  const edm::InputTag t1mumetTag;
  const edm::InputTag t1elmetTag;
  const edm::InputTag t1phmetTag;
  const edm::InputTag t1taumetTag;

  edm::EDGetTokenT<edm::View<pat::MET> >  t1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elemetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1taumetToken;

  // MET breakdown
  const bool addMETBreakDown;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetHadronHFToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetEgammaHFToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetChargedHadronToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetNeutralHadronToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetElectronsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetPhotonsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetMuonsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetUnclusteredToken;

  // Puppi MET
  const bool addPuppiMET;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1elemetToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1phmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1taumetToken;

  // MET systematics
  const bool addMETSystematics;

  // MVA met
  const bool addMVAMet;
  edm::EDGetTokenT<edm::View<reco::MET> > mvaMETToken;

  // inner bools
  const bool applyHLTFilter;
  const bool setHLTFilterFlag;
  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;   
  const double dRCleaningAK4;
  const double dRCleaningAK8;
  const std::string jetidwp;
  const std::string pileupjetidwp;
  const bool   applypileupjetid;
  const double btaggingCSVWP;
  const double btaggingMVAWP;
  const double minJetPtCountAK4;
  const double minJetPtBveto;
  const double minJetPtAK4Store;
 
  // Jet AK8
  const bool addSubstructureCHS;
  const bool addSubstructurePuppi;  
  edm::EDGetTokenT<std::vector<pat::Jet> > boostedJetsToken;
  edm::EDGetTokenT<std::vector<pat::Jet> > boostedPuppiJetsToken;
  TString boostedJetsCHSLabel;
  TString boostedJetsPuppiLabel;
  
  // Btag SF 
  bool addBTagScaleFactor;
  edm::FileInPath bTagScaleFactorFileCSV;
  edm::FileInPath bTagScaleFactorFileMVA;
  edm::FileInPath bTagScaleFactorFileSubCSV;

  BTagCalibration calibCSV;
  BTagCalibration calibMVA;
  BTagCalibration calibSubCSV;

  std::vector<BTagCalibrationReader> bMediumCSV;
  std::vector<BTagCalibrationReader> bMediumMVA;
  std::vector<BTagCalibrationReader> bMediumSubCSV;

  // additional info
  const bool addPhotonIDVariables;
  edm::EDGetTokenT<std::vector<pat::Photon> > photonIDCollectionToken;
  const bool addElectronIDVariables;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronIDCollectionToken;

  // inner vectors
  std::vector<std::string>   triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  std::vector<std::string>   filterPathsVector;
  std::map<std::string, int> filterPathsMap;

  std::unique_ptr<HLTPrescaleProvider> hltPrescaleProvider_;

  // tree
  TTree* tree;

  // pileup info
  int32_t puobs,putrue; 
  int32_t wzid,l1id,l2id;
  int32_t wzid_h,q1id,q2id;
  int32_t top_1,top_2;  
  int32_t mu1pid,mu2pid,mu1id,mu2id,mu1idm,mu2idm,mu1idt,mu2idt;
  int32_t el1pid,el2pid,el1id,el1idl,el2id,el2idl;
  int32_t tau1pid,tau2pid;
  int32_t phidl,phidm,phidt,phidh,parid,ancid; 

  // event info
  uint32_t event, run, lumi;  
  uint32_t nvtx;
  uint32_t nmuons,ntightmuons,nhighptmuons;
  uint32_t nelectrons,nlooseelectrons,ntightelectrons,nheepelectrons;
  uint32_t ntaus,ntausraw,ntausold,ntausrawold,nphotons;
  uint32_t njets,nbjets,nbjetslowpt,nbjetsMVA,nbjetsMVAlowpt;  
  uint32_t npuppijets,npuppibjets,npuppibjetsMVA,npuppibjetslowpt,npuppibjetsMVAlowpt;
  uint32_t njetsinc,npuppijetsinc;

  // trigger and met filters flags 
  uint8_t hltmet90,hltmet100,hltmet110,hltmet120;
  uint8_t hltmetwithmu90,hltmetwithmu100,hltmetwithmu110,hltmetwithmu120,hltmetwithmu170,hltmetwithmu300;
  uint8_t hltjetmet;
  uint8_t hltphoton165,hltphoton175,hltphoton120,hltphoton90,hltphoton120vbf;
  uint8_t hltdoublemu,hltsinglemu,hltdoubleel,hltsingleel,hltsingleel27,hltelnoiso;
  uint8_t hltPFHT400, hltPFHT475, hltPFHT600, hltPFHT650, hltPFHT800,hltPFHT900;
  uint8_t hltEcalHT800;
  uint8_t hltphoton90PFHT;

  //pre-scales
  double pswgt_ph120,pswgt_ph90;
  double pswgt_ht400,pswgt_ht475,pswgt_ht600,pswgt_ht650,pswgt_ht800,pswgt_ht900;
  
  //met filters
  uint8_t flagcsctight,flaghbhenoise,flaghbheiso,flageebadsc,flagecaltp,flaggoodvertices, flagglobaltighthalo, flagbadchpf, flagbadpfmu;
  
  // muon, ele, dilepton info
  double mu1pt,mu1eta,mu1phi,mu1pfpt,mu1pfeta,mu1pfphi,mu1iso,mu2pt,mu2eta,mu2phi,mu2pfpt,mu2pfeta,mu2pfphi,mu2iso;
  double el1pt,el1eta,el1phi,ele1e,el2pt,ele2e,el2eta,el2phi,phpt,pheta,phphi,phe;
  double tau1pt,tau1eta,tau1phi,tau1m,tau1iso,tau2pt,tau2eta,tau2phi,tau2m,tau2iso;
  double zmass,zpt,zeta,zphi,wmt,zeemass,zeept,zeeeta,zeephi,wemt,zttmass,zttpt,ztteta,zttphi,wtmt; 
  double emumass,emupt,emueta,emuphi,taumumass,taumupt,taumueta,taumuphi,tauemass,tauept,taueeta,tauephi;

  // photon purity studies
  uint32_t nphotonsPurity;
  double    phPHiso,  phCHiso, phNHiso, phPuritypt, phPurityeta, phPurityphi;
  double    phPurityPHiso,phPurityRND04PHiso,phPurityRND08PHiso,phPurityCHiso,phPurityRND04CHiso,phPurityRND08CHiso,phPurityNHiso;
  double    phPuritysieie, phPurityhoe, phPurityElectronVeto, phPurityEA,phPurityEAEGamma;
 
  // PF MET info (typeI and Raw)
  double rho;
  double t1pfmet,t1pfmetphi,t1mumet,t1mumetphi,t1elmet,t1elmetphi,t1phmet,t1phmetphi,t1taumet,t1taumetphi;
  double pfmet,pfmetphi,mumet,mumetphi,elmet,elmetphi,phmet,phmetphi,taumet,taumetphi;
  double calomet, calometphi;

  // MET break down
  double pfmethadronHF,pfmethadronHFphi,pfmetegammaHF,pfmetegammaHFphi,pfmetchargedhadron,pfmetchargedhadronphi;
  double pfmetneutralhadron,pfmetneutralhadronphi,pfmetelectrons,pfmetelectronsphi,pfmetmuons,pfmetmuonsphi,pfmetphotons,pfmetphotonsphi,pfmetunclustered,pfmetunclusteredphi;

  // Puppi MET info (typeI and Raw)
  double puppipfmet,puppipfmetphi,puppimumet,puppimumetphi,puppielmet,puppielmetphi,puppiphmet,puppiphmetphi,puppitaumet,puppitaumetphi;
  double puppit1pfmet,puppit1pfmetphi,puppit1mumet,puppit1mumetphi,puppit1elmet,puppit1elmetphi,puppit1phmet,puppit1phmetphi,puppit1taumet,puppit1taumetphi;

  // mva met
  double mvamet,mvametphi;

  // gen met
  double genmet,genmetphi;

  // met systematics
  double t1pfmetMuEnUp,t1pfmetMuEnDown,t1pfmetElEnUp,t1pfmetElEnDown,t1pfmetPhoEnUp,t1pfmetPhoEnDown,t1pfmetTauEnUp,t1pfmetTauEnDown;
  double t1pfmetJetEnUp,t1pfmetJetEnDown,t1pfmetJetResUp,t1pfmetJetResDown,t1pfmetUncEnUp,t1pfmetUncEnDown,t1pfmetJetSmear,t1pfmetXY;
  double t1pfmetMuEnUpPhi,t1pfmetMuEnDownPhi,t1pfmetElEnUpPhi,t1pfmetElEnDownPhi,t1pfmetPhoEnUpPhi,t1pfmetPhoEnDownPhi,t1pfmetTauEnUpPhi,t1pfmetTauEnDownPhi;
  double t1pfmetJetEnUpPhi,t1pfmetJetEnDownPhi,t1pfmetJetResUpPhi,t1pfmetJetResDownPhi,t1pfmetUncEnUpPhi,t1pfmetUncEnDownPhi,t1pfmetJetSmearPhi,t1pfmetXYPhi;

  // met systematics puppi
  double puppit1pfmetMuEnUp,puppit1pfmetMuEnDown,puppit1pfmetElEnUp,puppit1pfmetElEnDown,puppit1pfmetPhoEnUp,puppit1pfmetPhoEnDown,puppit1pfmetTauEnUp;
  double puppit1pfmetTauEnDown,puppit1pfmetJetEnUp,puppit1pfmetJetEnDown,puppit1pfmetJetResUp,puppit1pfmetJetResDown,puppit1pfmetUncEnUp,puppit1pfmetUncEnDown;
  double puppit1pfmetMuEnUpPhi,puppit1pfmetMuEnDownPhi,puppit1pfmetElEnUpPhi,puppit1pfmetElEnDownPhi;
  double puppit1pfmetPhoEnUpPhi,puppit1pfmetPhoEnDownPhi,puppit1pfmetTauEnUpPhi,puppit1pfmetTauEnDownPhi;
  double puppit1pfmetJetEnUpPhi,puppit1pfmetJetEnDownPhi,puppit1pfmetJetResUpPhi,puppit1pfmetJetResDownPhi;
  double puppit1pfmetUncEnUpPhi,puppit1pfmetUncEnDownPhi,puppit1pfmetJetSmearPhi,puppit1pfmetXYPhi;

  // AK4CHS combine jet
  std::vector<double> combinejetpt,combinejeteta,combinejetphi,combinejetm,combinejetbtag,combinejetbtagMVA;
  std::vector<double> combinejetCHfrac,combinejetNHfrac,combinejetEMfrac,combinejetCEMfrac,combinejetmetdphi;
  std::vector<double> combinejetHFlav,combinejetPFlav,combinejetQGL,combinejetPUID, combinejetPassPUID; 
  std::vector<double> combinejetGenpt,combinejetGeneta,combinejetGenphi,combinejetGenm;
  std::vector<double> combinejetBtagSF,combinejetBtagSFUp,combinejetBtagSFDown;
  std::vector<double> combinejetBtagMVASF,combinejetBtagMVASFUp,combinejetBtagMVASFDown;

  double incjetmetdphimin,incjetmumetdphimin,incjetelmetdphimin,incjetphmetdphimin,jetjetdphi,ht,htinc,ht30;
  double incjetmetdphimin4,incjetmumetdphimin4,incjetelmetdphimin4,incjetphmetdphimin4; 
  double alljetmetdphimin,alljetmetdphimin4,alljetmumetdphimin,alljetmumetdphimin4,alljetelmetdphimin,alljetelmetdphimin4,alljetphmetdphimin,alljetphmetdphimin4;

  // AK4Puppi combine jet
  std::vector<double> combinePuppijetpt,combinePuppijeteta,combinePuppijetphi,combinePuppijetm,combinePuppijetbtag,combinePuppijetbtagMVA;
  std::vector<double> combinePuppijetCHfrac,combinePuppijetNHfrac,combinePuppijetEMfrac,combinePuppijetCEMfrac,combinePuppijetmetdphi;
  std::vector<double> combinePuppijetHFlav,combinePuppijetPFlav,combinePuppijetQGL;
  std::vector<double> combinePuppijetGenpt,combinePuppijetGeneta,combinePuppijetGenphi,combinePuppijetGenm;
  std::vector<double> combinePuppijetBtagSF,combinePuppijetBtagSFUp,combinePuppijetBtagSFDown;
  std::vector<double> combinePuppijetBtagMVASF,combinePuppijetBtagMVASFUp,combinePuppijetBtagMVASFDown;

  double Puppijetmetdphimin,incPuppijetmetdphimin,Puppijetmumetdphimin,incPuppijetmumetdphimin,Puppijetelmetdphimin;
  double incPuppijetelmetdphimin,Puppijetphmetdphimin,incPuppijetphmetdphimin,PuppijetPuppijetdphi;
  double Puppijetmetdphimin4,incPuppijetmetdphimin4,Puppijetmumetdphimin4,incPuppijetmumetdphimin4;
  double Puppijetelmetdphimin4,incPuppijetelmetdphimin4,Puppijetphmetdphimin4,incPuppijetphmetdphimin4,Puppiht; 
  
  // AK8 CHS jets
  std::vector<double> boostedJetpt,boostedJeteta,boostedJetphi,boostedJetm;
  std::vector<double> boostedJetGenpt,boostedJetGenm,boostedJetGeneta,boostedJetGenphi;
  std::vector<double> boostedJettau1,boostedJettau2,boostedJettau3,boostedJettau4;
  std::vector<double> boostedJetGentau1,boostedJetGentau2,boostedJetGentau3,boostedJetGentau4;
  std::vector<double> boostedJetHFlav,boostedJetPFlav,boostedJetQGL,boostedJetBtag,boostedJetDoubleBtag;

  std::vector<double> prunedJetpt,prunedJetm,prunedJetphi,prunedJeteta;
  std::vector<double> prunedJetm_v2, prunedJetpt_v2, prunedJetphi_v2, prunedJeteta_v2;
  std::vector<double> prunedJetGenpt,prunedJetGenm,prunedJetGeneta,prunedJetGenphi;
  std::vector<double> prunedJetptraw,prunedJetmraw;
  std::vector<double> prunedJetHFlav,prunedJetPFlav,prunedJetQGL,prunedJetBtag,prunedJetDoubleBtag;

  std::vector<double> softDropJetpt,softDropJetm,softDropJeteta,softDropJetphi;
  std::vector<double> softDropJetm_v2, softDropJetpt_v2, softDropJetphi_v2, softDropJeteta_v2;
  std::vector<double> softDropJetGenpt,softDropJetGenm,softDropJetGeneta,softDropJetGenphi;
  std::vector<double> softDropJetHFlav,softDropJetPFlav,softDropJetQGL,softDropJetBtag,softDropJetDoubleBtag;
  std::vector<double> softDropJetptraw,softDropJetmraw;

  std::vector<double> prunedSubJetpt_1,prunedSubJetm_1,prunedSubJetphi_1,prunedSubJeteta_1;
  std::vector<double> prunedSubJetHFlav_1,prunedSubJetQGL_1,prunedSubJetBtag_1;
  std::vector<double> prunedSubJetGenpt_1,prunedSubJetGenm_1,prunedSubJetGeneta_1,prunedSubJetGenphi_1,prunedSubJetPFlav_1;
  std::vector<double> prunedSubJetptraw_1,prunedSubJetmraw_1;
  std::vector<double> prunedSubJetBtagSF_1,prunedSubJetBtagSFUp_1,prunedSubJetBtagSFDown_1;

  std::vector<double> prunedSubJetpt_2,prunedSubJetm_2,prunedSubJetphi_2,prunedSubJeteta_2,prunedSubJetHFlav_2,prunedSubJetQGL_2,prunedSubJetBtag_2;
  std::vector<double> prunedSubJetGenpt_2,prunedSubJetGenm_2,prunedSubJetGeneta_2,prunedSubJetGenphi_2,prunedSubJetPFlav_2;
  std::vector<double> prunedSubJetptraw_2,prunedSubJetmraw_2;
  std::vector<double> prunedSubJetBtagSF_2,prunedSubJetBtagSFUp_2,prunedSubJetBtagSFDown_2;

  std::vector<double> softDropSubJetpt_1,softDropSubJetm_1,softDropSubJetphi_1,softDropSubJeteta_1;
  std::vector<double> softDropSubJetHFlav_1,softDropSubJetQGL_1,softDropSubJetBtag_1;
  std::vector<double> softDropSubJetGenpt_1,softDropSubJetGenm_1,softDropSubJetGeneta_1,softDropSubJetGenphi_1,softDropSubJetPFlav_1;
  std::vector<double> softDropSubJetptraw_1,softDropSubJetmraw_1;
  std::vector<double> softDropSubJetBtagSF_1,softDropSubJetBtagSFUp_1,softDropSubJetBtagSFDown_1;

  std::vector<double> softDropSubJetpt_2,softDropSubJetm_2,softDropSubJetphi_2,softDropSubJeteta_2,softDropSubJetHFlav_2; 
  std::vector<double> softDropSubJetQGL_2,softDropSubJetBtag_2;
  std::vector<double> softDropSubJetGenpt_2,softDropSubJetGenm_2,softDropSubJetGeneta_2,softDropSubJetGenphi_2,softDropSubJetPFlav_2;
  std::vector<double> softDropSubJetptraw_2,softDropSubJetmraw_2;
  std::vector<double> softDropSubJetBtagSF_2,softDropSubJetBtagSFUp_2,softDropSubJetBtagSFDown_2;

  // AK8 Puppi jets
  std::vector<double> boostedPuppiJetpt,boostedPuppiJeteta,boostedPuppiJetphi,boostedPuppiJetm;
  std::vector<double> boostedPuppiJetGenpt,boostedPuppiJetGenm,boostedPuppiJetGeneta,boostedPuppiJetGenphi;
  std::vector<double> boostedPuppiJettau1,boostedPuppiJettau2,boostedPuppiJettau3,boostedPuppiJettau4;
  std::vector<double> boostedPuppiJetGentau1,boostedPuppiJetGentau2,boostedPuppiJetGentau3,boostedPuppiJetGentau4;
  std::vector<double> boostedPuppiJetHFlav,boostedPuppiJetPFlav,boostedPuppiJetQGL,boostedPuppiJetBtag,boostedPuppiJetDoubleBtag;

  std::vector<double> prunedPuppiJetpt,prunedPuppiJetm,prunedPuppiJetphi,prunedPuppiJeteta;
  std::vector<double> prunedPuppiJetm_v2, prunedPuppiJetpt_v2, prunedPuppiJetphi_v2, prunedPuppiJeteta_v2;
  std::vector<double> prunedPuppiJetGenpt,prunedPuppiJetGenm,prunedPuppiJetGeneta,prunedPuppiJetGenphi;
  std::vector<double> prunedPuppiJetptraw,prunedPuppiJetmraw;
  std::vector<double> prunedPuppiJetHFlav,prunedPuppiJetPFlav,prunedPuppiJetQGL,prunedPuppiJetBtag,prunedPuppiJetDoubleBtag;

  std::vector<double> softDropPuppiJetpt,softDropPuppiJetm,softDropPuppiJeteta,softDropPuppiJetphi;
  std::vector<double> softDropPuppiJetm_v2, softDropPuppiJetpt_v2, softDropPuppiJeteta_v2, softDropPuppiJetphi_v2;
  std::vector<double> softDropPuppiJetGenpt,softDropPuppiJetGenm,softDropPuppiJetGeneta,softDropPuppiJetGenphi;
  std::vector<double> softDropPuppiJetHFlav,softDropPuppiJetPFlav,softDropPuppiJetQGL,softDropPuppiJetBtag,softDropPuppiJetDoubleBtag;
  std::vector<double> softDropPuppiJetptraw,softDropPuppiJetmraw;

  std::vector<double> prunedPuppiSubJetpt_1,prunedPuppiSubJetm_1,prunedPuppiSubJetphi_1,prunedPuppiSubJeteta_1,prunedPuppiSubJetHFlav_1;
  std::vector<double> prunedPuppiSubJetQGL_1,prunedPuppiSubJetBtag_1;
  std::vector<double> prunedPuppiSubJetGenpt_1,prunedPuppiSubJetGenm_1,prunedPuppiSubJetGeneta_1,prunedPuppiSubJetGenphi_1,prunedPuppiSubJetPFlav_1;
  std::vector<double> prunedPuppiSubJetptraw_1,prunedPuppiSubJetmraw_1;
  std::vector<double> prunedPuppiSubJetBtagSF_1,prunedPuppiSubJetBtagSFUp_1,prunedPuppiSubJetBtagSFDown_1;

  std::vector<double> prunedPuppiSubJetpt_2,prunedPuppiSubJetm_2,prunedPuppiSubJetphi_2,prunedPuppiSubJeteta_2,prunedPuppiSubJetHFlav_2;
  std::vector<double> prunedPuppiSubJetQGL_2,prunedPuppiSubJetBtag_2;
  std::vector<double> prunedPuppiSubJetGenpt_2,prunedPuppiSubJetGenm_2,prunedPuppiSubJetGeneta_2,prunedPuppiSubJetGenphi_2,prunedPuppiSubJetPFlav_2;
  std::vector<double> prunedPuppiSubJetptraw_2,prunedPuppiSubJetmraw_2;
  std::vector<double> prunedPuppiSubJetBtagSF_2,prunedPuppiSubJetBtagSFUp_2,prunedPuppiSubJetBtagSFDown_2;

  std::vector<double> softDropPuppiSubJetpt_1,softDropPuppiSubJetm_1,softDropPuppiSubJetphi_1,softDropPuppiSubJeteta_1;
  std::vector<double> softDropPuppiSubJetHFlav_1,softDropPuppiSubJetQGL_1,softDropPuppiSubJetBtag_1;
  std::vector<double> softDropPuppiSubJetGenpt_1,softDropPuppiSubJetGenm_1,softDropPuppiSubJetGeneta_1,softDropPuppiSubJetGenphi_1,softDropPuppiSubJetPFlav_1;
  std::vector<double> softDropPuppiSubJetptraw_1,softDropPuppiSubJetmraw_1;
  std::vector<double> softDropPuppiSubJetBtagSF_1,softDropPuppiSubJetBtagSFUp_1,softDropPuppiSubJetBtagSFDown_1;

  std::vector<double> softDropPuppiSubJetpt_2,softDropPuppiSubJetm_2,softDropPuppiSubJetphi_2,softDropPuppiSubJeteta_2,softDropPuppiSubJetHFlav_2; 
  std::vector<double> softDropPuppiSubJetQGL_2,softDropPuppiSubJetBtag_2;
  std::vector<double> softDropPuppiSubJetGenpt_2,softDropPuppiSubJetGenm_2,softDropPuppiSubJetGeneta_2,softDropPuppiSubJetGenphi_2,softDropPuppiSubJetPFlav_2;
  std::vector<double> softDropPuppiSubJetptraw_2,softDropPuppiSubJetmraw_2;
  std::vector<double> softDropPuppiSubJetBtagSF_2,softDropPuppiSubJetBtagSFUp_2,softDropPuppiSubJetBtagSFDown_2;

  // special branches for photon efficiency variables dump
  std::vector<double> photonPt, photonEta, photonPhi, photonE, photonSCEta, photonSCPhi, photonSCEnergy, photonSCRawEnergy;
  std::vector<double> photonHOverE, photonSigmaIetaIeta, photonChargedIso,photonNeutralIso,photonEMIso, photonElectronVeto;
  std::vector<double> electronPt, electronEta, electronPhi, electronE, electronSCEta, electronSCPhi, electronSCEnergy;
  std::vector<double> electronSCRawEnergy, electronHOverE, electronSigmaIetaIeta, electronChargedIso, electronNeutralIso, electronEMIso, electronGsfPt;
  std::vector<double> electronEOP, electronDphi, electronDeta, electronMissHit, electronConversion, electronDxy, electronDz;

  // gen info leptoni W/Z boson (1 per event)
  double wzmass,wzmt,wzpt,wzeta,wzphi,wzmothid,l1pt,l1eta,l1phi,l2pt,l2eta,l2phi;
  // photon info
  double parpt,pareta,parphi,ancpt,anceta,ancphi,ancmass;
  int32_t ismatch, isdirect;
  // hadronic V and related quarks (1 per event)
  double wzmass_h,wzmt_h,wzpt_h,wzeta_h,wzphi_h,q1pt,q1eta,q1phi,q2pt,q2eta,q2phi;
  // one top
  double topmass,toppt,topeta,topphi;
  // second top
  double atopmass,atoppt,atopeta,atopphi;
  // DM mediator and DM particles
  double dmmass,dmpt,dmeta,dmphi,dmX1pt,dmX1eta,dmX1phi,dmX1mass,dmX2pt,dmX2eta,dmX2phi,dmX2mass;
  int   dmid,dmX1id,dmX2id;
  // for fastSIM
  double samplemedM,sampledmM;
  // weights
  double wgt,kfact,puwgt;

  // Trigger objects //ND
  uint32_t                   trig_obj_n; 
  std::vector<double>         trig_obj_pt, trig_obj_eta, trig_obj_phi; 
  std::vector< std::string > trig_obj_col; 

  int trig_L1A_check;
  int trig_L1A_n;
  std::vector< std::string > trig_L1A_list;

  std::vector<double> trig_L1EG_pt  , trig_L1EG_eta  , trig_L1EG_phi  ; 
  std::vector<double> trig_L1Jet_pt , trig_L1Jet_eta , trig_L1Jet_phi ; 
  std::vector<double> trig_L1Mu_pt  , trig_L1Mu_eta  , trig_L1Mu_phi  ; 
  double              trig_L1ETM_pt , trig_L1ETM_phi , trig_L1HTM_pt  , trig_L1HTM_phi;
  double              trig_L1ETT_pt , trig_L1ETT_phi , trig_L1HTT_pt  , trig_L1HTT_phi; 
  
 
  // sorting objects
  template<typename T> 
  class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }

  };

  PatPtSorter<pat::JetRef>      jetSorter;
  PatPtSorter<pat::MuonRef>     muonSorter;
  PatPtSorter<pat::ElectronRef> electronSorter;
  PatPtSorter<pat::PhotonRef>   photonSorter;
  PatPtSorter<pat::TauRef>      tauSorter;
  
};


MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
  ///////////// GEN INFO
  // isMC or Data --> default Data
  isMC          (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  isTriggerTree (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  addTriggerObjects (iConfig.existsAs<bool>("addTriggerObjects") ? iConfig.getParameter<bool>("addTriggerObjects") : false),
  // use lhe weights or not
  uselheweights (iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
  lheEventTag   (iConfig.getParameter<edm::InputTag>("lheinfo")),
  lheRunTag     (iConfig.getParameter<edm::InputTag>("lheRuninfo")),
  // is signal sample or not
  isSignalSample   (iConfig.existsAs<bool>("isSignalSample") ? iConfig.getParameter<bool>("isSignalSample") : false),
  addGenParticles  (iConfig.existsAs<bool>("addGenParticles") ? iConfig.getParameter<bool>("addGenParticles") : false),
  // xsec
  xsec  (iConfig.existsAs<double>("xsec") ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.),
  ///////////// TRIGGER and filter info INFO
  triggerResultsTag  (iConfig.getParameter<edm::InputTag>("triggerResults")),
  filterResultsTag   (iConfig.getParameter<edm::InputTag>("filterResults")),
  prescalesTag       (iConfig.getParameter<edm::InputTag>("prescales")),
  badChargedCandidateTag (iConfig.getParameter<edm::InputTag>("badChargedCandidate")),
  badPFMuonTag           (iConfig.getParameter<edm::InputTag>("badPFMuon")),
  triggerObjectsTag  (iConfig.existsAs<edm::InputTag>("triggerObjects") ? iConfig.getParameter<edm::InputTag>("triggerObjects") : edm::InputTag("")),
  IT_L1_EG           (iConfig.existsAs<edm::InputTag>("triggerL1EG")   ?  iConfig.getParameter<edm::InputTag>("triggerL1EG")   : edm::InputTag("")), 
  IT_L1_Jet          (iConfig.existsAs<edm::InputTag>("triggerL1Jet")  ?  iConfig.getParameter<edm::InputTag>("triggerL1Jet")  : edm::InputTag("")), 
  IT_L1_Mu           (iConfig.existsAs<edm::InputTag>("triggerL1Mu")   ?  iConfig.getParameter<edm::InputTag>("triggerL1Mu")   : edm::InputTag("")),
  IT_L1_Sums         (iConfig.existsAs<edm::InputTag>("triggerL1Sums") ?  iConfig.getParameter<edm::InputTag>("triggerL1Sums") : edm::InputTag("")), 
  IT_L1_Algos        (iConfig.existsAs<edm::InputTag>("triggerL1algos") ? iConfig.getParameter<edm::InputTag>("triggerL1algos") : edm::InputTag("gtStage2Digis")),
  // vertexes
  verticesTag (iConfig.getParameter<edm::InputTag>("vertices")),
  // rho
  rhoToken    (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  //muons
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  tightmuonsTag  (iConfig.getParameter<edm::InputTag>("tightmuons")),
  highptmuonsTag (iConfig.getParameter<edm::InputTag>("highptmuons")),
  // electrons
  electronsTag      (iConfig.getParameter<edm::InputTag>("electrons")),
  looseelectronsTag (iConfig.getParameter<edm::InputTag>("looseelectrons")),
  tightelectronsTag (iConfig.getParameter<edm::InputTag>("tightelectrons")),
  heepelectronsTag  (iConfig.getParameter<edm::InputTag>("heepelectrons")),
  // photons
  photonsTag        (iConfig.getParameter<edm::InputTag>("photons")),
  mediumphotonsTag  (iConfig.getParameter<edm::InputTag>("mediumphotons")),
  tightphotonsTag   (iConfig.getParameter<edm::InputTag>("tightphotons")),
  photonLooseIdTag  (iConfig.getParameter<edm::InputTag>("photonLooseId")),
  photonHighPtIdTag (iConfig.getParameter<edm::InputTag>("photonHighPtId")),
  // photon purity
  addPhotonPurity(iConfig.existsAs<bool>("addPhotonPurity") ? iConfig.getParameter<bool>("addPhotonPurity") : false),
  // taus
  tausNewTag(iConfig.getParameter<edm::InputTag>("taus")),
  tausOldTag(iConfig.getParameter<edm::InputTag>("tausOld")),
  tausNewRawTag(iConfig.getParameter<edm::InputTag>("tausRaw")),
  tausOldRawTag(iConfig.getParameter<edm::InputTag>("tausOldRaw")),
  // jets AK4
  jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
  addPuppiJets(iConfig.existsAs<bool>("addPuppiJets") ? iConfig.getParameter<bool>("addPuppiJets") : false),
  // met
  t1metTag(iConfig.getParameter<edm::InputTag>("t1met")),
  t1mumetTag(iConfig.getParameter<edm::InputTag>("t1mumet")),
  t1elmetTag(iConfig.getParameter<edm::InputTag>("t1elmet")),
  t1phmetTag(iConfig.getParameter<edm::InputTag>("t1phmet")),
  t1taumetTag(iConfig.getParameter<edm::InputTag>("t1taumet")),
  // add met brekdown
  addMETBreakDown(iConfig.existsAs<bool>("addMETBreakDown") ? iConfig.getParameter<bool>("addMETBreakDown") : false),
  // puppi met
  addPuppiMET(iConfig.existsAs<bool>("addPuppiMET") ? iConfig.getParameter<bool>("addPuppiMET") : false),
  // MET Systematics
  addMETSystematics(iConfig.existsAs<bool>("addMETSystematics") ? iConfig.getParameter<bool>("addMETSystematics") : false),
  // MVA met
  addMVAMet(iConfig.existsAs<bool>("addMVAMet") ? iConfig.getParameter<bool>("addMVAMet") : false),
  //filter events based on trigger
  applyHLTFilter(iConfig.existsAs<bool>("applyHLTFilter") ? iConfig.getParameter<bool>("applyHLTFilter") : false),
  setHLTFilterFlag(iConfig.existsAs<bool>("setHLTFilterFlag") ? iConfig.getParameter<bool>("setHLTFilterFlag") : false),
  // booleans
  cleanMuonJet(iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : true),
  cleanElectronJet(iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : true),
  cleanPhotonJet(iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : true),
  dRCleaningAK4(iConfig.existsAs<double>("dRCleaningAK4") ? iConfig.getParameter<double>("dRCleaningAK4") : 0.4),
  dRCleaningAK8(iConfig.existsAs<double>("dRCleaningAK8") ? iConfig.getParameter<double>("dRCleaningAK8") : 0.8),
  jetidwp(iConfig.existsAs<std::string>("jetidwp") ? iConfig.getParameter<std::string>("jetidwp") : "loose"),
  pileupjetidwp(iConfig.existsAs<std::string>("pileupjetidwp") ? iConfig.getParameter<std::string>("pileupjetidwp") : "medium"),
  applypileupjetid(iConfig.existsAs<bool>("applypileupjetid") ? iConfig.getParameter<bool>("applypileupjetid") : false),
  btaggingCSVWP(iConfig.getParameter<double>("btaggingCSVWP")),
  btaggingMVAWP(iConfig.getParameter<double>("btaggingMVAWP")),
  minJetPtCountAK4(iConfig.existsAs<double>("minJetPtCountAK4") ? iConfig.getParameter<double>("minJetPtCountAK4") : 30),
  minJetPtBveto(iConfig.existsAs<double>("minJetPtBveto") ? iConfig.getParameter<double>("minJetPtBveto") : 15),
  minJetPtAK4Store(iConfig.existsAs<double>("minJetPtAK4Store") ? iConfig.getParameter<double>("minJetPtAK4Store") : 20),
  // substructure
  addSubstructureCHS(iConfig.existsAs<bool>("addSubstructureCHS") ? iConfig.getParameter<bool>("addSubstructureCHS") : false),
  addSubstructurePuppi(iConfig.existsAs<bool>("addSubstructurePuppi") ? iConfig.getParameter<bool>("addSubstructurePuppi") : false),
  // btagging
  addBTagScaleFactor(iConfig.existsAs<bool>("addBTagScaleFactor") ? iConfig.getParameter<bool>("addBTagScaleFactor") : false),
  addPhotonIDVariables(iConfig.existsAs<bool>("addPhotonIDVariables") ? iConfig.getParameter<bool>("addPhotonIDVariables") : false),
  addElectronIDVariables(iConfig.existsAs<bool>("addElectronIDVariables") ? iConfig.getParameter<bool>("addElectronIDVariables") : false){
  
  usesResource();
  usesResource("TFileService");

  // trigger tokens
  triggerResultsToken   = consumes<edm::TriggerResults> (triggerResultsTag);
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(prescalesTag);
  filterResultsToken    = consumes<edm::TriggerResults> (filterResultsTag);
  if(isTriggerTree and addTriggerObjects){
      triggerObjectsToken   = consumes<pat::TriggerObjectStandAloneCollection>(triggerObjectsTag); 
      T_L1EG                = consumes<l1t::EGammaBxCollection>(IT_L1_EG  ); //ND
      T_L1Jet               = consumes<l1t::JetBxCollection   >(IT_L1_Jet ); //ND
      T_L1Mu                = consumes<l1t::MuonBxCollection  >(IT_L1_Mu  ); //ND
      T_L1Sums              = consumes<l1t::EtSumBxCollection >(IT_L1_Sums); //ND
      T_L1Algos             = consumes<GlobalAlgBlkBxCollection>(IT_L1_Algos); //ND
  }
  badChargedCandidateToken = consumes<bool>(badChargedCandidateTag);
  badPFMuonToken           = consumes<bool>(badPFMuonTag);

  // vertex
  verticesToken  = consumes<std::vector<reco::Vertex> > (verticesTag);

  //muons
  muonsToken       = consumes<pat::MuonRefVector> (muonsTag);
  tightmuonsToken  = consumes<pat::MuonRefVector> (tightmuonsTag);
  highptmuonsToken = consumes<pat::MuonRefVector> (highptmuonsTag);

  // electrons
  electronsToken       = consumes<pat::ElectronRefVector> (electronsTag);
  looseelectronsToken  = consumes<pat::ElectronRefVector> (looseelectronsTag);
  tightelectronsToken  = consumes<pat::ElectronRefVector> (tightelectronsTag);
  heepelectronsToken   = consumes<pat::ElectronRefVector> (heepelectronsTag);
  // photons
  photonsToken        = consumes<pat::PhotonRefVector> (photonsTag);
  mediumphotonsToken  = consumes<pat::PhotonRefVector> (mediumphotonsTag);
  tightphotonsToken   = consumes<pat::PhotonRefVector> (tightphotonsTag);
  photonLooseIdToken  = consumes<edm::ValueMap<bool> > (photonLooseIdTag);
  photonHighPtIdToken = consumes<edm::ValueMap<bool> > (photonHighPtIdTag);

  if(addPhotonPurity){
    photonsPurityToken      = consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photonsPurity"));
    tightphotonsPurityToken = consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("tightphotonsPurity")); 
    photonsieieToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonsieie")); 
    photonPHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonPHiso")); 
    photonCHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonCHiso")); 
    photonNHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonNHiso")); 
    rndgammaiso04Token      = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("rndgammaiso04")); 
    rndgammaiso08Token      = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("rndgammaiso08")); 
    gammaisoToken           = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("gammaiso"));
    rndchhadiso04Token      = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("rndchhadiso04"));  
    rndchhadiso08Token      = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("rndchhadiso08"));  
  }

  // taus
  tausNewToken = consumes<pat::TauRefVector> (tausNewTag);
  tausOldToken = consumes<pat::TauRefVector> (tausOldTag);
  tausNewRawToken = consumes<pat::TauRefVector> (tausNewRawTag);
  tausOldRawToken = consumes<pat::TauRefVector> (tausOldRawTag);

  // jets AK4
  jetsToken = consumes<std::vector<pat::Jet> > (jetsTag);
  t1metToken    = consumes<edm::View<pat::MET> > (t1metTag);
  t1mumetToken  = consumes<edm::View<pat::MET> > (t1mumetTag);
  t1elemetToken = consumes<edm::View<pat::MET> > (t1elmetTag);
  t1phmetToken  = consumes<edm::View<pat::MET> > (t1phmetTag);
  t1taumetToken  = consumes<edm::View<pat::MET> > (t1taumetTag);

  if(addMETBreakDown){
    pfMetHadronHFToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetHadronHF"));
    pfMetEgammaHFToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetEgammaHF"));
    pfMetChargedHadronToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetChargedHadron"));
    pfMetNeutralHadronToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetNeutralHadron"));
    pfMetElectronsToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetElectrons"));
    pfMetPhotonsToken   = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetPhotons"));
    pfMetMuonsToken     = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetMuons"));
    pfMetUnclusteredToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetUnclustered"));
  }
 
  // only for simulated samples
  if( isMC ){
    pileupInfoToken = consumes<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileup"));
    genevtInfoToken = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genevt"));
    lheInfoToken    = consumes<LHEEventProduct> (lheEventTag);
    lheRunInfoToken = consumes<LHERunInfoProduct,edm::InRun> (lheRunTag);
    gensToken       = consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("gens"));   
  }
  
  // consumes puppi jets
  if(addPuppiJets)
    puppijetsToken = consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("puppijets"));

  // consumes puppi met
  if(addPuppiMET){
    puppit1metToken    = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1met"));
    puppit1mumetToken  = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1mumet"));
    puppit1elemetToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1elmet"));
    puppit1phmetToken  = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1phmet"));		      
    puppit1taumetToken  = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1taumet"));		      
  }
  
  // consumes MVA met
  if(addMVAMet)
    mvaMETToken = consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("mvaMET"));

  // consumes boosted CHS jets
  if(addSubstructureCHS){
    boostedJetsToken =  consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("boostedJetsCHS"));
    boostedJetsCHSLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsCHS").label().c_str());
    boostedJetsCHSLabel.ReplaceAll("packed","");
    boostedJetsCHSLabel.ReplaceAll("PatJets","");
    boostedJetsCHSLabel.ReplaceAll("patJets","");
  }

  // consumes puppi boosted jets
  if(addSubstructurePuppi){
    boostedPuppiJetsToken =  consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("boostedJetsPuppi"));
    boostedJetsPuppiLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsPuppi").label().c_str());
    boostedJetsPuppiLabel.ReplaceAll("packed","");
    boostedJetsPuppiLabel.ReplaceAll("PatJets","");
    boostedJetsPuppiLabel.ReplaceAll("patJets","");
  }

  if(addBTagScaleFactor and isMC){
  
    bTagScaleFactorFileCSV = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileCSV");
    if ( bTagScaleFactorFileCSV.location()!=edm::FileInPath::Local)
      throw cms::Exception("MonoJetTreeMaker") << " Failed to find File = " << bTagScaleFactorFileCSV << " !!\n";

    calibCSV = BTagCalibration("CSVv2",bTagScaleFactorFileCSV.fullPath());
    bMediumCSV.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_B,"comb");
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_C,"comb");
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_UDSG,"incl");

    bTagScaleFactorFileMVA = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileMVA");
    if ( bTagScaleFactorFileMVA.location()!=edm::FileInPath::Local)
      throw cms::Exception("MonoJetTreeMaker") << " Failed to find File = " << bTagScaleFactorFileMVA << " !!\n";

    calibMVA = BTagCalibration("CMVAv2",bTagScaleFactorFileMVA.fullPath());
    bMediumMVA.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_B,"ttbar");
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_C,"ttbar");
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_UDSG,"incl");

    if(addSubstructureCHS or addSubstructurePuppi){

      bTagScaleFactorFileSubCSV = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileSubCSV");
      if ( bTagScaleFactorFileSubCSV.location()!=edm::FileInPath::Local)
	throw cms::Exception("MonoJetTreeMaker") << " Failed to find File = " << bTagScaleFactorFileSubCSV << " !!\n";
      
      calibSubCSV = BTagCalibration("CSVv2",bTagScaleFactorFileSubCSV.fullPath());
      bMediumSubCSV.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
      bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_B,"lt");
      bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_C,"lt");
      bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_UDSG,"incl");      
    }    
  }

  if(addPhotonIDVariables){
    photonIDCollectionToken =  consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonIDCollection"));
    if(not addPhotonPurity){
      photonsieieToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonsieie"));
      photonPHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonPHiso"));
      photonCHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonCHiso"));
      photonNHisoToken        = consumes<edm::ValueMap<double> > (iConfig.getParameter<edm::InputTag>("photonNHiso"));
    }
  }
  if(addElectronIDVariables){
    electronIDCollectionToken =  consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronIDCollection"));    
  }

  readDMFromGenParticle = false;

  hltPrescaleProvider_.reset(new HLTPrescaleProvider(iConfig, consumesCollector(), *this)); //ND
}


MonoJetTreeMaker::~MonoJetTreeMaker() {}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace reco;
    using namespace std;
    using namespace pat;

    // get the rho value
    Handle<double> rhoH;
    iEvent.getByToken(rhoToken, rhoH);
    rho = *rhoH;

    // Get handles to all the requisite collections
    // TRIGGER and FILTERS
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByToken(triggerResultsToken, triggerResultsH);
    Handle<pat::PackedTriggerPrescales> triggerPrescalesH;
    iEvent.getByToken(triggerPrescalesToken, triggerPrescalesH);
    Handle<TriggerResults> filterResultsH;
    iEvent.getByToken(filterResultsToken, filterResultsH);

    Handle<pat::TriggerObjectStandAloneCollection>   triggerObjectsH; 
    edm::Handle<l1t::EGammaBxCollection> H_L1EG;
    edm::Handle<l1t::TauBxCollection>    H_L1Tau;
    edm::Handle<l1t::JetBxCollection>    H_L1Jet;
    edm::Handle<l1t::MuonBxCollection>   H_L1Mu;
    edm::Handle<l1t::EtSumBxCollection>  H_L1Sums;
    edm::Handle<GlobalAlgBlkBxCollection> H_L1Algos;

    if(isTriggerTree and addTriggerObjects){
      iEvent.getByToken(triggerObjectsToken, triggerObjectsH);
      iEvent.getByToken(T_L1EG  , H_L1EG);
      iEvent.getByToken(T_L1Jet , H_L1Jet);
      iEvent.getByToken(T_L1Mu  , H_L1Mu);
      iEvent.getByToken(T_L1Sums, H_L1Sums);
      iEvent.getByToken(T_L1Algos, H_L1Algos);
    }

    edm::Handle<bool> filterBadChCandH;
    iEvent.getByToken(badChargedCandidateToken,filterBadChCandH);      
    edm::Handle<bool> filterBadPFMuonH;
    iEvent.getByToken(badPFMuonToken,filterBadPFMuonH);
   

    // GEN INFO    
    Handle<vector<PileupSummaryInfo> > pileupInfoH;
    Handle<GenEventInfoProduct>        genevtInfoH;
    Handle<LHEEventProduct>            lheInfoH;
    Handle<View<GenParticle> >         gensH;

    if(isMC){
      iEvent.getByToken(pileupInfoToken, pileupInfoH);
      if (uselheweights){
	iEvent.getByToken(genevtInfoToken, genevtInfoH);
	iEvent.getByToken(lheInfoToken, lheInfoH);
      }
      if (addGenParticles or isSignalSample)
	iEvent.getByToken(gensToken, gensH);
    }

    // VERTEX
    Handle<vector<Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);

    // MUONS
    Handle<pat::MuonRefVector> muonsH;
    iEvent.getByToken(muonsToken, muonsH);
    pat::MuonRefVector muons = *muonsH;

    Handle<pat::MuonRefVector> tightmuonsH;
    iEvent.getByToken(tightmuonsToken, tightmuonsH);
    pat::MuonRefVector tightmuons = *tightmuonsH;

    Handle<pat::MuonRefVector> highptmuonsH;
    iEvent.getByToken(highptmuonsToken, highptmuonsH);
    pat::MuonRefVector highptmuons = *highptmuonsH;

    // ELECTRONS
    Handle<pat::ElectronRefVector> electronsH;
    iEvent.getByToken(electronsToken, electronsH);
    pat::ElectronRefVector electrons = *electronsH;

    Handle<pat::ElectronRefVector> looseelectronsH;
    iEvent.getByToken(looseelectronsToken, looseelectronsH);
    pat::ElectronRefVector looseelectrons = *looseelectronsH;

    Handle<pat::ElectronRefVector> tightelectronsH;
    iEvent.getByToken(tightelectronsToken, tightelectronsH);
    pat::ElectronRefVector tightelectrons = *tightelectronsH;

    Handle<pat::ElectronRefVector> heepelectronsH;
    iEvent.getByToken(heepelectronsToken, heepelectronsH);
    pat::ElectronRefVector heepelectrons = *heepelectronsH;

    // PHOTONS
    Handle<pat::PhotonRefVector> photonsH;
    iEvent.getByToken(photonsToken, photonsH);
    pat::PhotonRefVector photons = *photonsH;

    Handle<pat::PhotonRefVector> mediumPhotonsH;
    iEvent.getByToken(mediumphotonsToken, mediumPhotonsH);
    pat::PhotonRefVector mediumphotons = *mediumPhotonsH;

    Handle<pat::PhotonRefVector> tightPhotonsH;
    iEvent.getByToken(tightphotonsToken, tightPhotonsH);
    pat::PhotonRefVector tightphotons = *tightPhotonsH;

    Handle<ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdToken, photonLooseIdH);
    Handle<ValueMap<bool> > photonHighPtIdH;
    iEvent.getByToken(photonHighPtIdToken, photonHighPtIdH);

    Handle<pat::PhotonRefVector> photonsPurityH;
    pat::PhotonRefVector photonsPurity;
    Handle<pat::PhotonRefVector> tightphotonsPurityH;
    pat::PhotonRefVector tightphotonsPurity;
    Handle<edm::ValueMap<double> > photonsieieH;
    Handle<edm::ValueMap<double> > photonPHisoH;
    Handle<edm::ValueMap<double> > photonCHisoH;
    Handle<edm::ValueMap<double> > photonNHisoH;
    Handle<edm::ValueMap<double> > rndgammaiso04H;
    Handle<edm::ValueMap<double> > rndgammaiso08H;
    Handle<edm::ValueMap<double> > gammaisoH;
    Handle<edm::ValueMap<double> > rndchhadiso04H;
    Handle<edm::ValueMap<double> > rndchhadiso08H;

    if(addPhotonPurity){
      iEvent.getByToken(photonsPurityToken, photonsPurityH);    
      photonsPurity = *photonsPurityH;
      iEvent.getByToken(tightphotonsPurityToken, tightphotonsPurityH);
      tightphotonsPurity = *tightphotonsPurityH;
      iEvent.getByToken(photonsieieToken, photonsieieH);
      iEvent.getByToken(photonPHisoToken, photonPHisoH);      
      iEvent.getByToken(photonCHisoToken, photonCHisoH);
      iEvent.getByToken(photonNHisoToken, photonNHisoH);
      iEvent.getByToken(rndgammaiso04Token, rndgammaiso04H);
      iEvent.getByToken(rndgammaiso08Token, rndgammaiso08H);      
      iEvent.getByToken(gammaisoToken, gammaisoH);
      iEvent.getByToken(rndchhadiso04Token, rndchhadiso04H);
      iEvent.getByToken(rndchhadiso08Token, rndchhadiso08H);
    }


    // TAUS
    Handle<pat::TauRefVector > tausNewH;
    iEvent.getByToken(tausNewToken, tausNewH);
    pat::TauRefVector tausNew = *tausNewH;

    Handle<pat::TauRefVector > tausOldH;
    iEvent.getByToken(tausOldToken, tausOldH);
    pat::TauRefVector tausOld = *tausOldH;

    Handle<pat::TauRefVector > tausNewRawH;
    iEvent.getByToken(tausNewRawToken, tausNewRawH);
    pat::TauRefVector tausNewRaw = *tausNewRawH;

    Handle<pat::TauRefVector > tausOldRawH;
    iEvent.getByToken(tausOldRawToken, tausOldRawH);
    pat::TauRefVector tausOldRaw = *tausOldRawH;

    // AK4 Jets
    Handle<vector<pat::Jet> > jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<vector<pat::Jet> > jetsPuppiH;
    if(addPuppiJets)
      iEvent.getByToken(puppijetsToken, jetsPuppiH);

    // MET
    Handle<View<pat::MET> > t1metH;
    iEvent.getByToken(t1metToken, t1metH);
    Handle<View<pat::MET> > puppit1metH;
    if(addPuppiMET)
      iEvent.getByToken(puppit1metToken, puppit1metH);

    Handle<View<pat::MET> > t1mumetH;
    iEvent.getByToken(t1mumetToken, t1mumetH);
    Handle<View<pat::MET> > puppit1mumetH;
    if(addPuppiMET)
      iEvent.getByToken(puppit1mumetToken, puppit1mumetH);

    Handle<View<pat::MET> > t1elmetH;
    iEvent.getByToken(t1elemetToken, t1elmetH);
    Handle<View<pat::MET> > puppit1elmetH;
    if(addPuppiMET)
      iEvent.getByToken(puppit1elemetToken, puppit1elmetH);

    Handle<View<pat::MET> > t1phmetH;
    iEvent.getByToken(t1phmetToken, t1phmetH);
    Handle<View<pat::MET> > puppit1phmetH;
    if(addPuppiMET)
      iEvent.getByToken(puppit1phmetToken, puppit1phmetH);

    Handle<View<pat::MET> > t1taumetH;
    iEvent.getByToken(t1taumetToken, t1taumetH);
    Handle<View<pat::MET> > puppit1taumetH;
    if(addPuppiMET)
      iEvent.getByToken(puppit1taumetToken, puppit1taumetH);

    // MET breakdown 
    Handle<View<pat::MET> > pfMetHadronHFH;
    Handle<View<pat::MET> > pfMetEgammaHFH;
    Handle<View<pat::MET> > pfMetChargedHadronH;
    Handle<View<pat::MET> > pfMetNeutralHadronH;
    Handle<View<pat::MET> > pfMetElectronsH;
    Handle<View<pat::MET> > pfMetPhotonsH;
    Handle<View<pat::MET> > pfMetMuonsH;
    Handle<View<pat::MET> > pfMetUnclusteredH;

    if(addMETBreakDown){
      iEvent.getByToken(pfMetHadronHFToken,pfMetHadronHFH);
      iEvent.getByToken(pfMetEgammaHFToken,pfMetEgammaHFH);
      iEvent.getByToken(pfMetChargedHadronToken,pfMetChargedHadronH);
      iEvent.getByToken(pfMetNeutralHadronToken,pfMetNeutralHadronH);
      iEvent.getByToken(pfMetElectronsToken,pfMetElectronsH);
      iEvent.getByToken(pfMetPhotonsToken,pfMetPhotonsH);
      iEvent.getByToken(pfMetMuonsToken,pfMetMuonsH);
      iEvent.getByToken(pfMetUnclusteredToken,pfMetUnclusteredH);
    }

    /// MVA MET
    Handle<View<reco::MET> > mvaMetH;
    if(addMVAMet)
      iEvent.getByToken(mvaMETToken,mvaMetH);

    // boosted jets
    Handle<vector<pat::Jet> > boostedJetsH;
    if(addSubstructureCHS)
      iEvent.getByToken(boostedJetsToken,boostedJetsH);
    
    Handle<vector<pat::Jet> > boostedPuppiJetsH;
    if(addSubstructurePuppi)
      iEvent.getByToken(boostedPuppiJetsToken,boostedPuppiJetsH);

    // Photon and electron ID
    Handle<vector<pat::Photon> > photonIDH;
    if(addPhotonIDVariables){
      iEvent.getByToken(photonIDCollectionToken,photonIDH);
      if(not addPhotonPurity){
	iEvent.getByToken(photonsieieToken, photonsieieH);
	iEvent.getByToken(photonPHisoToken, photonPHisoH);      
	iEvent.getByToken(photonCHisoToken, photonCHisoH);
	iEvent.getByToken(photonNHisoToken, photonNHisoH);
      }
    }
    Handle<vector<pat::Electron> > electronIDH;
    if(addElectronIDVariables)
      iEvent.getByToken(electronIDCollectionToken,electronIDH);
    

    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();
    
    // Trigger info
    hltmet90        = 0; hltmet100       = 0; hltmet110       = 0; hltmet120       = 0;
    hltmetwithmu90  = 0; hltmetwithmu100 = 0; hltmetwithmu110 = 0; hltmetwithmu120 = 0;
    hltmetwithmu170 = 0; hltmetwithmu300 = 0;
    hltjetmet       = 0; 
    hltphoton90     = 0; hltphoton120    = 0; hltphoton120vbf = 0;
    hltphoton165    = 0; hltphoton175    = 0;
    hltdoublemu     = 0;
    hltsinglemu     = 0;
    hltdoubleel     = 0;
    hltsingleel     = 0;
    hltsingleel27   = 0;
    hltelnoiso      = 0;
    hltPFHT400      = 0; hltPFHT475 = 0; hltPFHT600 = 0; hltPFHT650 = 0; hltPFHT800 = 0; hltPFHT900 = 0; 
    hltEcalHT800    = 0; 
    hltphoton90PFHT = 0;

    if(triggerResultsH.isValid() and setHLTFilterFlag == false){
      for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;	
	if(triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]]))

        if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger

        if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet100       = 1; // MET trigger
        if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet110       = 1; // MET trigger

        if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger

        if (i == 8  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu90   = 1; // MET trigger
        if (i == 9  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu100  = 1; // MET trigger
        if (i == 10 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu110  = 1; // MET trigger
        if (i == 11  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu120 = 1; // MET trigger

        if (i == 12  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 13  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 14 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170  = 1; // MET trigger
        if (i == 15 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170  = 1; // MET trigger

        if (i == 16 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 17 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 18 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger

        if (i == 19 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
        if (i == 20 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
        if (i == 21 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger

        if (i == 22 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
        if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger

        if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
        if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
        if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger

        if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165    = 1; // Photon trigger
        if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175    = 1; // Photon trigger
        if (i == 29 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton120    = 1; // Photon trigger
        if (i == 30 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90     = 1; // Photon trigger
        if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton120vbf = 1; // Photon trigger

        if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 32 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 33 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 34 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger

        if (i == 35 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 36 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 37 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 38 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 39 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 40 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger

        if (i == 41 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 42 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger

        if (i == 43 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 44 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 45 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 45 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel27   = 1; // Single electron trigger
        if (i == 46 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 47 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger

        if (i == 48 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger
        if (i == 49 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger

        if (i == 50 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT400      = 1; // jet ht
        if (i == 51 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT475      = 1; // jet ht
        if (i == 52 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT600      = 1; // jet ht
        if (i == 53 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT650      = 1; // jet ht
        if (i == 54 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT800      = 1; // jet ht
        if (i == 55 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT900      = 1; // jet ht
	if (i == 56 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltEcalHT800    = 1;
	if (i == 57 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90PFHT = 1;
	if (i == 58 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90PFHT = 1;

      }
    }
    else if(setHLTFilterFlag == true){
      hltmet90        = 1; hltmet100       = 1; hltmet110       = 1; hltmet120       = 1;
      hltmetwithmu90  = 1; hltmetwithmu100 = 1; hltmetwithmu110 = 1; hltmetwithmu120 = 1;
      hltmetwithmu170 = 1; hltmetwithmu300 = 1;
      hltjetmet       = 1; 
      hltphoton90     = 1; hltphoton120    = 1; hltphoton120vbf = 1;
      hltphoton165    = 1; hltphoton175    = 1;
      hltdoublemu     = 1;
      hltsinglemu     = 1;
      hltdoubleel     = 1;
      hltsingleel     = 1;
      hltsingleel27   = 1;
      hltelnoiso      = 1;
      hltPFHT400   = 1; hltPFHT475 = 1; hltPFHT600 = 1; hltPFHT650 = 1; hltPFHT800 = 1; hltPFHT900 = 1;
      hltEcalHT800 = 1; 
      hltphoton90PFHT = 1;
    }

    bool triggered = false;
    if (hltmet90        == 1) triggered = true;
    if (hltmet100       == 1) triggered = true;
    if (hltmet110       == 1) triggered = true;
    if (hltmet120       == 1) triggered = true;
    if (hltmetwithmu90  == 1) triggered = true;
    if (hltmetwithmu100 == 1) triggered = true;
    if (hltmetwithmu110 == 1) triggered = true;
    if (hltmetwithmu120 == 1) triggered = true;
    if (hltmetwithmu170 == 1) triggered = true;
    if (hltmetwithmu300 == 1) triggered = true;
    if (hltjetmet       == 1) triggered = true;
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltphoton120    == 1) triggered = true;
    if (hltphoton120vbf == 1) triggered = true;
    if (hltphoton90    == 1) triggered = true;
    if (hltdoublemu     == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltdoubleel     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (hltelnoiso      == 1) triggered = true;
    if (hltPFHT400      == 1) triggered = true;
    if (hltPFHT475      == 1) triggered = true;
    if (hltPFHT600      == 1) triggered = true;
    if (hltPFHT650      == 1) triggered = true;
    if (hltPFHT800      == 1) triggered = true;
    if (hltPFHT900      == 1) triggered = true;
    if (hltEcalHT800    == 1) triggered = true;
    if (hltphoton90PFHT == 1) triggered = true;

    if (applyHLTFilter && !triggered) return;

    pswgt_ph120 = 1.0;
    pswgt_ph90  = 1.0;
    pswgt_ht400 = 1.0; 
    pswgt_ht475 = 1.0; 
    pswgt_ht600 = 1.0; 
    pswgt_ht650 = 1.0; 
    pswgt_ht800 = 1.0; 
    pswgt_ht900 = 1.0;

    const edm::TriggerNames &trignames = iEvent.triggerNames(*triggerResultsH);
    for (size_t i = 0; i < triggerResultsH->size(); i++) {
        if (trignames.triggerName(i).find("HLT_Photon120_v") != string::npos) pswgt_ph120 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_Photon90_v") != string::npos) pswgt_ph90 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT400_v") != string::npos) pswgt_ht400 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT475_v") != string::npos) pswgt_ht475 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT600_v") != string::npos) pswgt_ht600 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT650_v") != string::npos) pswgt_ht650 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT800_v") != string::npos) pswgt_ht800 = triggerPrescalesH->getPrescaleForIndex(i);
        if (trignames.triggerName(i).find("HLT_PFHT900_v") != string::npos) pswgt_ht900 = triggerPrescalesH->getPrescaleForIndex(i);
    }
    
    // MET filter info
    flagcsctight  = 0;
    flaghbhenoise = 0;
    flaghbheiso   = 0;
    flageebadsc   = 0;
    flagecaltp    = 0;
    flaggoodvertices = 0;
    flagglobaltighthalo = 0;
    flagbadchpf = *filterBadChCandH;
    flagbadpfmu = *filterBadPFMuonH;

    if(filterResultsH.isValid()){
      // Which MET filters passed
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
      }
    }

    if(isTriggerTree and addTriggerObjects) {
      fillTriggerObjects(triggerObjectsH, trignames); //dump HL trigger objects
      fillTriggerL1(H_L1EG, H_L1Tau, H_L1Jet, H_L1Mu, H_L1Sums); //dump L1 objects
      fillAlgosL1(iEvent, iSetup, H_L1Algos); // dump L1 bits
    }

    // Pileup info -- Will need to the updated to the Run-II specifications
    if(verticesH.isValid()) nvtx = verticesH->size();
    else nvtx = 0;

    puobs  = 0;
    putrue = 0;
    puwgt  = 1.;

    // in caase the cross section is not set from outside --> fix to 1 as dummy value
    if (uselheweights && genevtInfoH.isValid())
      wgt = genevtInfoH->weight();    
    else wgt = 1.0;
    
    if (pileupInfoH.isValid()) {
      for (auto pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
	if (pileupInfo_iter->getBunchCrossing() == 0) {
	  puobs  = pileupInfo_iter->getPU_NumInteractions();
	  putrue = pileupInfo_iter->getTrueNumInteractions();
	}
      }
    }

    // MET information 
    if(t1metH.isValid()){      
      // dump gen met info
      if(t1metH->front().genMET()){
	genmet    = t1metH->front().genMET()->pt();
	genmetphi = t1metH->front().genMET()->phi();
      }
      else {
	genmet = -99.;
	genmetphi = -99.;
      }
      
      t1pfmet    = t1metH->front().corPt();
      t1pfmetphi = t1metH->front().corPhi();
      pfmet      = t1metH->front().uncorPt();
      pfmetphi   = t1metH->front().uncorPhi();
      calomet    = t1metH->front().caloMETPt();
      calometphi = t1metH->front().caloMETPhi(); 
    }
    else { 
      t1pfmet    = -99. ;
      t1pfmetphi = -99. ;
      pfmet      = -99. ;
      pfmetphi   = -99. ;
      calomet    = -99. ;
      calometphi = -99. ;
    }


    if(addMVAMet && mvaMetH.isValid()){
      mvamet    = mvaMetH->front().pt();
      mvametphi = mvaMetH->front().phi();
    }
    else{
      mvamet    = -99.;
      mvametphi = -99.;
    }

    if(t1mumetH.isValid()){
      t1mumet    = t1mumetH->front().corPt();
      t1mumetphi = t1mumetH->front().corPhi();
      mumet      = t1mumetH->front().uncorPt();
      mumetphi   = t1mumetH->front().uncorPhi();
    }
    else{
      t1mumet    = -99;
      t1mumetphi = -99;
      mumet      = -99;
      mumetphi   = -99;
    }

    if(t1elmetH.isValid()){
      t1elmet    = t1elmetH->front().corPt();
      t1elmetphi = t1elmetH->front().corPhi();
      elmet      = t1elmetH->front().uncorPt();
      elmetphi   = t1elmetH->front().uncorPhi();
    }
    else{
      t1elmet    = -99;
      t1elmetphi = -99;
      elmet      = -99;
      elmetphi   = -99;
    }

    if(t1phmetH.isValid()){
      t1phmet    = t1phmetH->front().corPt();
      t1phmetphi = t1phmetH->front().corPhi();
      phmet      = t1phmetH->front().uncorPt();
      phmetphi   = t1phmetH->front().uncorPhi();
    }
    else{
      t1phmet    = -99;
      t1phmetphi = -99;
      phmet      = -99;
      phmetphi   = -99;
    }

    if(t1taumetH.isValid()){
      t1taumet    = t1taumetH->front().corPt();
      t1taumetphi = t1taumetH->front().corPhi();
      taumet      = t1taumetH->front().uncorPt();
      taumetphi   = t1taumetH->front().uncorPhi();
    }
    else{
      t1taumet    = -99;
      t1taumetphi = -99;
      taumet      = -99;
      taumetphi   = -99;
    }

    if(addMETSystematics  and not isTriggerTree){          
      if(t1metH.isValid()){
	t1pfmetMuEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetMuEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetElEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetElEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetPhoEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetPhoEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetTauEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetTauEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetResUp   = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResUp,   pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetJetResDown = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetUncEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetUncEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetSmear   = t1metH->front().shiftedPt(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetXY         = t1metH->front().shiftedPt(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1XY);

	t1pfmetMuEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetMuEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetElEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetElEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetPhoEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetPhoEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetTauEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetTauEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetResUpPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResUp,   pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetJetResDownPhi = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetUncEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp,   pat::MET::METCorrectionLevel::Type1);
	t1pfmetUncEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);
	t1pfmetJetSmearPhi   = t1metH->front().shiftedPhi(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1Smear);
	t1pfmetXYPhi         = t1metH->front().shiftedPhi(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1XY);


      }
    }
    else{
      t1pfmetMuEnUp  = -99.; t1pfmetMuEnDown  = -99.; t1pfmetElEnUp   = -99.; t1pfmetElEnDown   = -99.;
      t1pfmetPhoEnUp = -99.; t1pfmetPhoEnDown = -99.; t1pfmetTauEnUp  = -99.; t1pfmetTauEnDown  = -99.;
      t1pfmetJetEnUp = -99.; t1pfmetJetEnDown = -99.; t1pfmetJetResUp = -99.; t1pfmetJetResDown = -99.;
      t1pfmetUncEnUp = -99.; t1pfmetUncEnDown = -99.; t1pfmetJetSmear = -99.; t1pfmetXY = -99.;     

      t1pfmetMuEnUpPhi  = -99.; t1pfmetMuEnDownPhi  = -99.; t1pfmetElEnUpPhi   = -99.; t1pfmetElEnDownPhi   = -99.;
      t1pfmetPhoEnUpPhi = -99.; t1pfmetPhoEnDownPhi = -99.; t1pfmetTauEnUpPhi  = -99.; t1pfmetTauEnDownPhi  = -99.;
      t1pfmetJetEnUpPhi = -99.; t1pfmetJetEnDownPhi = -99.; t1pfmetJetResUpPhi = -99.; t1pfmetJetResDownPhi = -99.;
      t1pfmetUncEnUpPhi = -99.; t1pfmetUncEnDownPhi = -99.; t1pfmetJetSmearPhi = -99.; t1pfmetXYPhi = -99.;     
    }
    
    // MET break down
    pfmethadronHF = -99. ; pfmethadronHFphi = -99. ; 
    pfmetegammaHF = -99. ; pfmetegammaHFphi = -99. ; 
    pfmetchargedhadron = -99. ; pfmetchargedhadronphi = -99.;
    pfmetneutralhadron = -99. ; pfmetneutralhadronphi = -99. ; 
    pfmetelectrons     = -99. ; pfmetelectronsphi   = -99. ; 
    pfmetmuons         = -99. ; pfmetmuonsphi       = -99. ; 
    pfmetphotons       = -99. ; pfmetphotonsphi     = -99. ; 
    pfmetunclustered   = -99. ; pfmetunclusteredphi = -99.;

    if(addMETBreakDown and not isTriggerTree){
      pfmethadronHF    = pfMetHadronHFH->front().pt();
      pfmethadronHFphi = pfMetHadronHFH->front().phi();
      pfmetegammaHF    = pfMetEgammaHFH->front().pt();
      pfmetegammaHFphi = pfMetEgammaHFH->front().phi();
      pfmetchargedhadron    = pfMetChargedHadronH->front().pt();
      pfmetchargedhadronphi = pfMetChargedHadronH->front().phi();
      pfmetneutralhadron    = pfMetNeutralHadronH->front().pt();
      pfmetneutralhadronphi = pfMetNeutralHadronH->front().phi();
      pfmetelectrons    = pfMetElectronsH->front().pt();
      pfmetelectronsphi = pfMetElectronsH->front().phi();
      pfmetmuons        = pfMetMuonsH->front().pt();
      pfmetmuonsphi     = pfMetMuonsH->front().phi();
      pfmetphotons      = pfMetPhotonsH->front().pt();
      pfmetphotonsphi   = pfMetPhotonsH->front().phi();
      pfmetunclustered  = pfMetUnclusteredH->front().pt();
      pfmetunclusteredphi = pfMetUnclusteredH->front().phi();	
    }
    

    // puppi met info    
    puppit1pfmet = -99.; puppit1pfmetphi = -99.;
    puppipfmet   = -99.; puppipfmetphi   = -99.;
    puppit1mumet = -99.; puppit1mumetphi = -99.;
    puppimumet   = -99.; puppimumetphi   = -99.;
    puppit1elmet = -99.; puppit1elmetphi = -99.;
    puppielmet   = -99.; puppielmetphi   = -99.;
    puppit1phmet = -99.; puppit1phmetphi = -99.;
    puppiphmet   = -99.; puppiphmetphi   = -99.;
    
    puppit1pfmetMuEnUp  = -99.; puppit1pfmetMuEnDown  = -99.; puppit1pfmetElEnUp   = -99.; puppit1pfmetElEnDown   = -99.;
    puppit1pfmetPhoEnUp = -99.; puppit1pfmetPhoEnDown = -99.; puppit1pfmetTauEnUp  = -99.; puppit1pfmetTauEnDown  = -99.;
    puppit1pfmetJetEnUp = -99.; puppit1pfmetJetEnDown = -99.; puppit1pfmetJetResUp = -99.; puppit1pfmetJetResDown = -99.;
    puppit1pfmetUncEnUp = -99.; puppit1pfmetUncEnDown = -99.;

    puppit1pfmetMuEnUpPhi  = -99.; puppit1pfmetMuEnDownPhi  = -99.; puppit1pfmetElEnUpPhi   = -99.; puppit1pfmetElEnDownPhi   = -99.;
    puppit1pfmetPhoEnUpPhi = -99.; puppit1pfmetPhoEnDownPhi = -99.; puppit1pfmetTauEnUpPhi  = -99.; puppit1pfmetTauEnDownPhi  = -99.;
    puppit1pfmetJetEnUpPhi = -99.; puppit1pfmetJetEnDownPhi = -99.; puppit1pfmetJetResUpPhi = -99.; puppit1pfmetJetResDownPhi = -99.;
    puppit1pfmetUncEnUpPhi = -99.; puppit1pfmetUncEnDownPhi = -99.;

    if(addPuppiMET and not isTriggerTree){

      if(puppit1metH.isValid()){
	puppit1pfmet        = puppit1metH->front().corPt();
	puppit1pfmetphi     = puppit1metH->front().corPhi();
	puppipfmet          = puppit1metH->front().uncorPt();
	puppipfmetphi       = puppit1metH->front().uncorPhi();
      }

      if(puppit1mumetH.isValid()){
	puppit1mumet        = puppit1mumetH->front().corPt();
	puppit1mumetphi     = puppit1mumetH->front().corPhi();
	puppimumet          = puppit1mumetH->front().uncorPt();
	puppimumetphi       = puppit1mumetH->front().uncorPhi();
      }
      
      if(puppit1elmetH.isValid()){
	puppit1elmet        = puppit1elmetH->front().corPt();
	puppit1elmetphi     = puppit1elmetH->front().corPhi();
	puppielmet          = puppit1elmetH->front().uncorPt();
	puppielmetphi       = puppit1elmetH->front().uncorPhi();
      }

      if(puppit1phmetH.isValid()){
	puppit1phmet        = puppit1phmetH->front().corPt();
	puppit1phmetphi     = puppit1phmetH->front().corPhi();
	puppiphmet          = puppit1phmetH->front().uncorPt();
	puppiphmetphi       = puppit1phmetH->front().uncorPhi();
      }

      if(puppit1taumetH.isValid()){
	puppit1taumet        = puppit1taumetH->front().corPt();
	puppit1taumetphi     = puppit1taumetH->front().corPhi();
	puppitaumet          = puppit1taumetH->front().uncorPt();
	puppitaumetphi       = puppit1taumetH->front().uncorPhi();
      }

      if(addMETSystematics and not isTriggerTree){	
	if(puppit1metH.isValid()){
	  puppit1pfmetMuEnUp     = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetMuEnDown   = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);       
	  puppit1pfmetElEnUp     = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetElEnDown   = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetPhoEnUp    = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetPhoEnDown  = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetTauEnUp    = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetTauEnDown  = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetJetEnUp    = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnDown  = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnUp    = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnDown  = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetJetResUp   = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetResUp, pat::MET::METCorrectionLevel::Type1Smear);
	  puppit1pfmetJetResDown = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);     
	  puppit1pfmetUncEnUp    = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetUncEnDown  = puppit1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);      

	  puppit1pfmetMuEnUpPhi     = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetMuEnDownPhi   = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);       
	  puppit1pfmetElEnUpPhi     = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetElEnDownPhi   = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetPhoEnUpPhi    = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetPhoEnDownPhi  = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetTauEnUpPhi    = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetTauEnDownPhi  = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetJetEnUpPhi    = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnDownPhi  = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnUpPhi    = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetJetEnDownPhi  = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);	
	  puppit1pfmetJetResUpPhi   = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResUp, pat::MET::METCorrectionLevel::Type1Smear);
	  puppit1pfmetJetResDownPhi = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);     
	  puppit1pfmetUncEnUpPhi    = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp, pat::MET::METCorrectionLevel::Type1);
	  puppit1pfmetUncEnDownPhi  = puppit1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);      
	}
      }
    }
    
    // AK4 Jet information
    vector<pat::JetRef> alljets;
    vector<pat::JetRef> incjets;
    vector<pat::JetRef> jets;
    if(jetsH.isValid()){      
      for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
	//clean from leptons
	bool skipjet = false;
	if(muonsH.isValid()){
	  for (std::size_t j = 0; j < muons.size(); j++) {
	    if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	      skipjet = true;
	  }
	}
	if(electronsH.isValid()){
	  for (std::size_t j = 0; j < electrons.size(); j++) {
	    if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	      skipjet = true;
	  }
	}
	if(photonsH.isValid()){
	  for (std::size_t j = 0; j < photons.size(); j++) {
	    if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	      skipjet = true;
	  }
	}
	
	// jet in overlap with lepton
	if (skipjet) continue;
      
	pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());	
	if(jetref.isAvailable() and jetref.isNonnull()) alljets.push_back(jetref);
	
	// apply jet id
	bool passjetid = applyJetID(*jets_iter,jetidwp);            
	if (!passjetid) 
	  continue;
	// apply pileup jet id
	bool passpuid = applyPileupJetID(*jets_iter,pileupjetidwp,false);
	if (applypileupjetid and !passpuid) continue; 
     
	if(jetref.isAvailable() and jetref.isNonnull())
	  incjets.push_back(jetref);
      }
    } 
    
    if(incjets.size() > 0) sort(incjets.begin(), incjets.end(), jetSorter);
    
    // only central jets
    for (size_t i = 0; i < incjets.size(); i++) {
      if (fabs(incjets[i]->eta()) <= 2.5) 
	jets.push_back(incjets[i]);
    }        
    
    // sort them in pt
    if(jets.size() > 0)  sort(jets.begin(), jets.end(), jetSorter);
    
    // count central jets    
    njets       = 0;
    njetsinc    = 0;
    nbjets      = 0;
    nbjetslowpt = 0;
    nbjetsMVA   = 0;
    nbjetsMVAlowpt = 0;
    
    combinejetpt        .clear(); combinejeteta       .clear(); combinejetphi       .clear(); combinejetbtag      .clear(); combinejetCHfrac    .clear();
    combinejetNHfrac    .clear(); combinejetEMfrac    .clear(); combinejetCEMfrac   .clear(); combinejetmetdphi   .clear();
    combinejetHFlav     .clear(); combinejetPFlav     .clear(); combinejetQGL       .clear(); combinejetPUID      .clear();
    combinejetGenpt     .clear(); combinejetGeneta    .clear(); combinejetGenphi    .clear(); combinejetGenm      .clear(); 
    combinejetm         .clear(); combinejetbtagMVA   .clear();
    combinejetBtagSF .clear(); combinejetBtagSFUp .clear(); combinejetBtagSFDown .clear();
    combinejetBtagMVASF .clear(); combinejetBtagMVASFUp .clear(); combinejetBtagMVASFDown .clear();
    combinejetPassPUID.clear();
    
    // all jets
    for(size_t i = 0; i < incjets.size(); i++){
      if(incjets[i]->pt() > minJetPtCountAK4) njetsinc++;
    }
    
    // only central jets
    for (size_t i = 0; i < jets.size(); i++) {
      
      if (jets[i]->pt() > minJetPtCountAK4) njets++;
      // btagging
      if (jets[i]->pt() > minJetPtCountAK4 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP) nbjets++;
      if (jets[i]->pt() > minJetPtBveto && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP) nbjetslowpt++;
      
      if (jets[i]->pt() > minJetPtCountAK4 && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP) nbjetsMVA++;
      if (jets[i]->pt() > minJetPtBveto && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP) nbjetsMVAlowpt++;
    }
    
    // fill collections
    for(size_t i = 0; i < incjets.size(); i++){
      if (incjets[i]->pt() > minJetPtAK4Store){
	
	combinejetpt.push_back(incjets[i]->pt());
	combinejeteta.push_back(incjets[i]->eta());
	combinejetphi.push_back(incjets[i]->phi());
	combinejetm.push_back(incjets[i]->mass());
	combinejetbtag.push_back(incjets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	combinejetbtagMVA.push_back(incjets[i]->bDiscriminator("pfCombinedMVAV2BJetTags"));
	combinejetCHfrac  .push_back(incjets[i]->chargedHadronEnergyFraction());
	combinejetNHfrac  .push_back(incjets[i]->neutralHadronEnergyFraction());
	combinejetEMfrac  .push_back(incjets[i]->neutralEmEnergyFraction());
	combinejetCEMfrac .push_back(incjets[i]->chargedEmEnergyFraction());
	
	if(incjets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	  combinejetQGL   .push_back(incjets[i]->userFloat("QGTagger:qgLikelihood")); 
	  // pileup jet id
	if(incjets[i]->hasUserFloat("puid:fullDiscriminant"))
	  combinejetPUID  .push_back(incjets[i]->userFloat("puid:fullDiscriminant"));	
	else
	  combinejetPUID  .push_back(incjets[i]->userFloat("pileupJetId:fullDiscriminant"));
	
	combinejetPassPUID.push_back(applyPileupJetID(*incjets[i],pileupjetidwp,false));	  
	// fill jet met dphi
	combinejetmetdphi.push_back(deltaPhi(incjets[i]->phi(), t1pfmetphi));
	
	// MC based info
	if(isMC){
	  combinejetHFlav.push_back(incjets[i]->hadronFlavour()); 
	  combinejetPFlav.push_back(incjets[i]->partonFlavour()); 
	  if(incjets[i]->genJet()){
	    combinejetGenpt.push_back(incjets[i]->genJet()->pt()); 
	    combinejetGeneta.push_back(incjets[i]->genJet()->eta()); 
	    combinejetGenphi.push_back(incjets[i]->genJet()->phi()); 
	    combinejetGenm.push_back(incjets[i]->genJet()->mass()); 
	  }
	  // b-tag SF for jets
	  if(addBTagScaleFactor){
	    calculateBtagSF(*incjets[i],"CSV",combinejetBtagSF,combinejetBtagSFUp,combinejetBtagSFDown);
	    calculateBtagSF(*incjets[i],"MVA",combinejetBtagMVASF,combinejetBtagMVASFUp,combinejetBtagMVASFDown);
	  }
	}  
      }         	
    }    
    
    // Jet met Dphi      
    jetjetdphi          = 0.0;   
    incjetmetdphimin    = 0.0; alljetmetdphimin    = 0.0;
    incjetmumetdphimin  = 0.0; alljetmumetdphimin  = 0.0;
    incjetelmetdphimin  = 0.0; alljetelmetdphimin  = 0.0;
    incjetphmetdphimin  = 0.0; alljetphmetdphimin  = 0.0;
    incjetmetdphimin4   = 0.0; alljetmetdphimin4   = 0.0;
    incjetmumetdphimin4 = 0.0; alljetmumetdphimin4 = 0.0;
    incjetelmetdphimin4 = 0.0; alljetelmetdphimin4 = 0.0;
    incjetphmetdphimin4 = 0.0; alljetphmetdphimin4 = 0.0;

    // delta phi between jets
    if (combinejetphi.size() > 1)
      jetjetdphi = deltaPhi(combinejetphi[0], combinejetphi[1]);
    
    std::vector<double> alljetmetdphiminvector;
    std::vector<double> alljetmetdphimin4vector;
    std::vector<double> alljetmumetdphiminvector;
    std::vector<double> alljetmumetdphimin4vector;
    std::vector<double> alljetelmetdphiminvector;
    std::vector<double> alljetelmetdphimin4vector;
    std::vector<double> alljetphmetdphiminvector;
    std::vector<double> alljetphmetdphimin4vector;
    
    for (size_t i = 0; i < alljets.size(); i++) {
      if (alljets[i]->pt() > minJetPtCountAK4) {	  
	double alljetphi = atan2(sin(alljets[i]->phi()), cos(alljets[i]->phi()));
	alljetmetdphiminvector  .push_back(fabs(deltaPhi(alljetphi, t1pfmetphi)));
	alljetmumetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1mumetphi)));
	alljetelmetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1elmetphi)));
	alljetphmetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1phmetphi)));
	if (i < 4) alljetmetdphimin4vector  .push_back(fabs(deltaPhi(alljetphi, t1pfmetphi)));
	if (i < 4) alljetmumetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1mumetphi)));
	if (i < 4) alljetelmetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1elmetphi)));
	if (i < 4) alljetphmetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1phmetphi)));
      }
    }
    if (alljetmetdphiminvector   .size() > 0) alljetmetdphimin    = *min_element(alljetmetdphiminvector   .begin(), alljetmetdphiminvector   .end());
    if (alljetmumetdphiminvector .size() > 0) alljetmumetdphimin  = *min_element(alljetmumetdphiminvector .begin(), alljetmumetdphiminvector .end());
    if (alljetelmetdphiminvector .size() > 0) alljetelmetdphimin  = *min_element(alljetelmetdphiminvector .begin(), alljetelmetdphiminvector .end());
    if (alljetphmetdphiminvector .size() > 0) alljetphmetdphimin  = *min_element(alljetphmetdphiminvector .begin(), alljetphmetdphiminvector .end());
    if (alljetmetdphimin4vector  .size() > 0) alljetmetdphimin4   = *min_element(alljetmetdphimin4vector  .begin(), alljetmetdphimin4vector  .end());
    if (alljetmumetdphimin4vector.size() > 0) alljetmumetdphimin4 = *min_element(alljetmumetdphimin4vector.begin(), alljetmumetdphimin4vector.end());
    if (alljetelmetdphimin4vector.size() > 0) alljetelmetdphimin4 = *min_element(alljetelmetdphimin4vector.begin(), alljetelmetdphimin4vector.end());
    if (alljetphmetdphimin4vector.size() > 0) alljetphmetdphimin4 = *min_element(alljetphmetdphimin4vector.begin(), alljetphmetdphimin4vector.end());
    
    // delta phi jet-met      
    std::vector<double> incjetmetdphiminvector;
    std::vector<double> incjetmetdphimin4vector;
    std::vector<double> incjetmumetdphiminvector;
    std::vector<double> incjetmumetdphimin4vector;
    std::vector<double> incjetelmetdphiminvector;
    std::vector<double> incjetelmetdphimin4vector;
    std::vector<double> incjetphmetdphiminvector;
    std::vector<double> incjetphmetdphimin4vector;

    for (size_t i = 0; i < incjets.size(); i++) {
      if (incjets[i]->pt() > minJetPtCountAK4) {
	double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
	incjetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1pfmetphi)));
	incjetmumetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
	incjetelmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1elmetphi)));
	incjetphmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1phmetphi)));
	if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1pfmetphi)));
	if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
	if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetphi)));
	if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetphi)));
      }
    }
    if (incjetmetdphiminvector .size() > 0) incjetmetdphimin  = *min_element(incjetmetdphiminvector .begin(), incjetmetdphiminvector .end());
    if (incjetmetdphimin4vector.size() > 0) incjetmetdphimin4 = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
    if (incjetmumetdphiminvector .size() > 0) incjetmumetdphimin  = *min_element(incjetmumetdphiminvector .begin(), incjetmumetdphiminvector .end());
    if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4 = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
    if (incjetelmetdphiminvector .size() > 0) incjetelmetdphimin  = *min_element(incjetelmetdphiminvector .begin(), incjetelmetdphiminvector .end());
    if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4 = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
    if (incjetphmetdphiminvector .size() > 0) incjetphmetdphimin  = *min_element(incjetphmetdphiminvector .begin(), incjetphmetdphiminvector .end());
    if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4 = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    
    // QCD suppression handles
    ht     = 0.;
    htinc  = 0.;
    ht30   = 0.;
    for (size_t i = 0; i < incjets.size(); i++) {
      if (incjets[i]->pt() > minJetPtCountAK4) {
	htinc += incjets[i]->pt(); 
	if (fabs(incjets[i]->eta()) < 3.0) 
	  ht30 += incjets[i]->pt();
      }
    }  
    for (size_t i = 0; i < jets.size(); i++) {
      if (jets[i]->pt() > minJetPtCountAK4)
	ht += jets[i]->pt();      
    }
    
    // PUPPI AK4 jets
    if(addPuppiJets and not isTriggerTree){
      
      vector<pat::JetRef> incPuppijets;
      vector<pat::JetRef> Puppijets;
      
      if(jetsPuppiH.isValid()){	
	for (auto jets_iter = jetsPuppiH->begin(); jets_iter != jetsPuppiH->end(); ++jets_iter) {	  
	  //clean from leptons
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
		skipjet = true;
	    }
	  }
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
		skipjet = true;
	    }
	  }
	  
	  if (skipjet) continue;
	  
	  // apply jet id
	  bool passjetid = applyJetID(*jets_iter,jetidwp);            
	  if (!passjetid) 
	    continue;
	  
	  pat::JetRef jetref(jetsPuppiH, jets_iter - jetsPuppiH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    incPuppijets.push_back(jetref);
	}
      }
      
      if(incPuppijets.size() > 0) sort(incPuppijets.begin(), incPuppijets.end(), jetSorter);
      
      // only central jets
      for (size_t i = 0; i < incPuppijets.size(); i++) {
	if (fabs(incPuppijets[i]->eta()) <= 2.5) 
	  Puppijets.push_back(incPuppijets[i]);
      }        
	
      // sort them in pt
      if(Puppijets.size() > 0)  sort(Puppijets.begin(), Puppijets.end(), jetSorter);
      
      // count central jets    
      npuppijets       = 0;
      npuppijetsinc    = 0;
      npuppibjets      = 0;
      npuppibjetslowpt = 0;
      npuppibjetsMVA   = 0;
      npuppibjetsMVAlowpt = 0;
      
      combinePuppijetpt        .clear(); combinePuppijeteta       .clear(); combinePuppijetphi       .clear(); combinePuppijetbtag      .clear(); 
      combinePuppijetCHfrac    .clear(); combinePuppijetbtagMVA   .clear();
      combinePuppijetNHfrac    .clear(); combinePuppijetEMfrac    .clear(); combinePuppijetCEMfrac   .clear(); combinePuppijetmetdphi   .clear();
      combinePuppijetHFlav     .clear(); combinePuppijetPFlav     .clear(); combinePuppijetQGL       .clear(); 
      combinePuppijetGenpt     .clear(); combinePuppijetGeneta    .clear(); combinePuppijetGenphi    .clear(); combinePuppijetGenm      .clear();
      combinePuppijetm         .clear(); 
      combinePuppijetBtagSF    .clear(); combinePuppijetBtagSFUp .clear(); combinePuppijetBtagSFDown .clear();
      combinePuppijetBtagMVASF    .clear(); combinePuppijetBtagMVASFUp .clear(); combinePuppijetBtagMVASFDown .clear();
      
      for(size_t i = 0; i < incPuppijets.size(); i++){
	if(incPuppijets[i]->pt() > minJetPtCountAK4) npuppijetsinc++;
      }
      
      for (size_t i = 0; i < Puppijets.size(); i++) {
	
	if (Puppijets[i]->pt() > minJetPtCountAK4) npuppijets++;
	if (Puppijets[i]->pt() > minJetPtCountAK4 && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP) npuppibjets++;
	if (Puppijets[i]->pt() > minJetPtBveto && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP) npuppibjetslowpt++;

	if (Puppijets[i]->pt() > minJetPtCountAK4 && Puppijets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP) npuppibjetsMVA++;
	if (Puppijets[i]->pt() > minJetPtBveto && Puppijets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP) npuppibjetsMVAlowpt++;
      }
      
      // fill collections
      for(size_t i = 0; i < incPuppijets.size(); i++){
	if (incPuppijets[i]->pt() > minJetPtCountAK4){
	  
	  combinePuppijetpt.push_back(incPuppijets[i]->pt());
	  combinePuppijeteta.push_back(incPuppijets[i]->eta());
	  combinePuppijetphi.push_back(incPuppijets[i]->phi());
	  combinePuppijetm.push_back(incPuppijets[i]->mass());
	  combinePuppijetbtag.push_back(incPuppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	  combinePuppijetbtagMVA.push_back(incPuppijets[i]->bDiscriminator("pfCombinedMVAV2BJetTags"));
	  combinePuppijetCHfrac  .push_back(incPuppijets[i]->chargedHadronEnergyFraction());
	  combinePuppijetNHfrac  .push_back(incPuppijets[i]->neutralHadronEnergyFraction());
	  combinePuppijetEMfrac  .push_back(incPuppijets[i]->neutralEmEnergyFraction());
	  combinePuppijetCEMfrac .push_back(incPuppijets[i]->chargedEmEnergyFraction());
	  combinePuppijetmetdphi.push_back(deltaPhi(incPuppijets[i]->phi(), puppit1pfmetphi));
	  
	  if(incPuppijets[i]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
	    combinePuppijetQGL   .push_back(incPuppijets[i]->userFloat("QGTaggerPuppi:qgLikelihood")); 
	  
	  // MC based info
	  if(isMC){
	    combinePuppijetHFlav.push_back(incPuppijets[i]->hadronFlavour()); 
	    combinePuppijetPFlav.push_back(incPuppijets[i]->partonFlavour()); 
	    if(incPuppijets[i]->genJet()){
	      combinePuppijetGenpt.push_back(incPuppijets[i]->genJet()->pt()); 
	      combinePuppijetGeneta.push_back(incPuppijets[i]->genJet()->eta()); 
	      combinePuppijetGenphi.push_back(incPuppijets[i]->genJet()->phi()); 
	      combinePuppijetGenm.push_back(incPuppijets[i]->genJet()->mass()); 
	    }	  
	    
	      // b-tag SF for Puppijets
	    if(addBTagScaleFactor){
	      calculateBtagSF(*incPuppijets[i],"CSV",combinePuppijetBtagSF,combinePuppijetBtagSFUp,combinePuppijetBtagSFDown);
	      calculateBtagSF(*incPuppijets[i],"MVA",combinePuppijetBtagMVASF,combinePuppijetBtagMVASFUp,combinePuppijetBtagMVASFDown);
	    }    	
	  }
	}
      }
      
      PuppijetPuppijetdphi    = 0.0;
      incPuppijetmetdphimin   = 0.0;
      incPuppijetmumetdphimin = 0.0;
      incPuppijetelmetdphimin = 0.0;
      incPuppijetphmetdphimin = 0.0;
      
      incPuppijetmetdphimin4  = 0.0;
      incPuppijetmumetdphimin4= 0.0;
      incPuppijetelmetdphimin4= 0.0;
      incPuppijetphmetdphimin4= 0.0;
      
      // delta phi between Puppijets
      if (combinePuppijetphi.size()>1) 
	PuppijetPuppijetdphi = deltaPhi(combinePuppijetphi[0], combinePuppijetphi[1]);
      
      std::vector<double> incpuppijetmetdphiminvector;
      std::vector<double> incpuppijetmetdphimin4vector;
      std::vector<double> incpuppijetmumetdphiminvector;
      std::vector<double> incpuppijetmumetdphimin4vector;
      std::vector<double> incpuppijetelmetdphiminvector;
      std::vector<double> incpuppijetelmetdphimin4vector;
      std::vector<double> incpuppijetphmetdphiminvector;
      std::vector<double> incpuppijetphmetdphimin4vector;
      
      for (size_t i = 0; i < incPuppijets.size(); i++) {
	if (incPuppijets[i]->pt() > minJetPtCountAK4) {
	  double incjetphi = atan2(sin(incPuppijets[i]->phi()), cos(incPuppijets[i]->phi()));
	  incpuppijetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1pfmetphi)));
	  incpuppijetmumetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1mumetphi)));
	  incpuppijetelmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1elmetphi)));
	  incpuppijetphmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1phmetphi)));
	  if (i < 4) incpuppijetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1pfmetphi)));
	  if (i < 4) incpuppijetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1mumetphi)));
	  if (i < 4) incpuppijetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1elmetphi)));
	  if (i < 4) incpuppijetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1phmetphi)));
	}
      }
      if (incpuppijetmetdphiminvector .size() > 0) incPuppijetmetdphimin  = *min_element(incpuppijetmetdphiminvector .begin(), incpuppijetmetdphiminvector .end());
      if (incpuppijetmetdphimin4vector.size() > 0) incPuppijetmetdphimin4 = *min_element(incpuppijetmetdphimin4vector.begin(), incpuppijetmetdphimin4vector.end());
      if (incpuppijetmumetdphiminvector .size() > 0) incPuppijetmumetdphimin  = *min_element(incpuppijetmumetdphiminvector .begin(), incpuppijetmumetdphiminvector .end());
      if (incpuppijetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4 = *min_element(incpuppijetmumetdphimin4vector.begin(), incpuppijetmumetdphimin4vector.end());
      if (incpuppijetelmetdphiminvector .size() > 0) incPuppijetelmetdphimin  = *min_element(incpuppijetelmetdphiminvector .begin(), incpuppijetelmetdphiminvector .end());
      if (incpuppijetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4 = *min_element(incpuppijetelmetdphimin4vector.begin(), incpuppijetelmetdphimin4vector.end());
      if (incpuppijetphmetdphiminvector .size() > 0) incPuppijetphmetdphimin  = *min_element(incpuppijetphmetdphiminvector .begin(), incpuppijetphmetdphiminvector .end());
      if (incpuppijetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4 = *min_element(incpuppijetphmetdphimin4vector.begin(), incpuppijetphmetdphimin4vector.end());
	
      // QCD suppression handles
      Puppiht     = 0.;
      for (size_t i = 0; i < Puppijets.size(); i++) {
	if (Puppijets[i]->pt() > minJetPtCountAK4) {
	  Puppiht += Puppijets[i]->pt();
	}
      }        
    }
    
    // Lepton part
    vector<pat::MuonRef> muonvector;
    if(muonsH.isValid() and tightmuonsH.isValid() and highptmuonsH.isValid()){
      nmuons          = muonsH->size();
      ntightmuons     = tightmuonsH->size();
      nhighptmuons    = highptmuonsH->size();
      
      for (size_t i = 0; i < muons.size(); i++) 
	muonvector.push_back(muons[i]);
    }

    vector<pat::ElectronRef> electronvector;
    if(electronsH.isValid() and looseelectronsH.isValid() and tightelectronsH.isValid() and heepelectronsH.isValid()){
      
      nelectrons      = electronsH->size();
      nlooseelectrons = looseelectronsH->size();
      ntightelectrons = tightelectronsH->size();
      nheepelectrons  = heepelectronsH->size();

      for (size_t i = 0; i < electrons.size(); i++) 
	electronvector.push_back(electrons[i]);
    }
    
    // re-apply the cleaning to be sure
    vector<pat::TauRef> tauvector;
    ntaus = 0;
    if(tausNewH.isValid()){
      for(std::size_t itau =0 ; itau < tausNew.size(); itau++){
	bool skiptau = false;
	for (std::size_t j = 0; j < muons.size(); j++) {
	  if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausNew[itau]->eta(), tausNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	for (std::size_t j = 0; j < electrons.size(); j++) {
	  if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausNew[itau]->eta(), tausNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	if(skiptau) continue;
	tauvector.push_back(tausNew[itau]);
	ntaus++;
      }
    }

 
    ntausold = 0;
    if(tausOldH.isValid()){
      for(std::size_t itau =0 ; itau < tausOld.size(); itau++){
	bool skiptau = false;
	for (std::size_t j = 0; j < muons.size(); j++) {
	  if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausOld[itau]->eta(), tausOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	for (std::size_t j = 0; j < electrons.size(); j++) {
	  if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausOld[itau]->eta(), tausOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	if(skiptau) continue;
	ntausold++;
      }
    }

    // old isolation
    ntausraw = 0;
    if(tausNewRawH.isValid()){
      for(std::size_t itau =0 ; itau < tausNewRaw.size(); itau++){
	bool skiptau = false;
	for (std::size_t j = 0; j < muons.size(); j++) {
	  if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausNewRaw[itau]->eta(), tausNewRaw[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	for (std::size_t j = 0; j < electrons.size(); j++) {
	  if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausNewRaw[itau]->eta(), tausNewRaw[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	if(skiptau) continue;
	ntausraw++;
      }
    }

    ntausrawold = 0;
    if(tausOldRawH.isValid()){
      for(std::size_t itau =0 ; itau < tausOldRaw.size(); itau++){
	bool skiptau = false;
	for (std::size_t j = 0; j < muons.size(); j++) {
	  if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausOldRaw[itau]->eta(), tausOldRaw[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	for (std::size_t j = 0; j < electrons.size(); j++) {
	  if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausOldRaw[itau]->eta(), tausOldRaw[itau]->phi()) < dRCleaningAK4) skiptau = true;
	}
	if(skiptau) continue;
	ntausrawold++;
      }
    }

    // W, Z control sample information
    zmass       = 0.0; zpt         = 0.0; zeta        = 0.0; zphi        = 0.0;
    zeemass     = 0.0; zeept       = 0.0; zeeeta      = 0.0; zeephi      = 0.0;
    zttmass     = 0.0; zttpt       = 0.0; ztteta      = 0.0; zttphi      = 0.0;
    wmt         = 0.0; wemt        = 0.0; wtmt        = 0.0; 
    emumass     = 0.0; emupt       = 0.0; emueta      = 0.0; emuphi      = 0.0;
    taumumass   = 0.0; taumupt     = 0.0; taumueta    = 0.0; taumuphi    = 0.0;
    tauemass    = 0.0; tauept      = 0.0; taueeta     = 0.0; tauephi     = 0.0;

    mu1pid      = 0;   mu1pt       = 0.0; mu1eta      = 0.0; mu1phi      = 0.0;
    mu1pfpt     = 0.0; mu1pfeta    = 0.0; mu1pfphi    = 0.0; mu1id       = 0;
    mu1idm      = 0;   mu1idt      = 0;   mu1iso      = 0.0;

    mu2pid      = 0;   mu2pt       = 0.0; mu2eta      = 0.0; mu2phi      = 0.0;
    mu2pfpt     = 0.0; mu2pfeta    = 0.0; mu2pfphi    = 0.0; mu2id       = 0;
    mu2idm      = 0;   mu2idt      = 0;   mu2iso      = 0.0;

    el1pid      = 0; el1pt       = 0.0; el1eta      = 0.0; el1phi      = 0.0; el1id       = 0;
    el2pid      = 0; el2pt       = 0.0; el2eta      = 0.0; el2phi      = 0.0; el2id       = 0;

    tau1pid     = 0;   tau1pt    = 0.0; tau1eta     = 0.0; tau1phi     = 0.0; tau1m       = 0.0; tau1iso = 0.0;
    tau2pid     = 0;   tau2pt    = 0.0; tau2eta     = 0.0; tau2phi     = 0.0; tau2m       = 0.0; tau2iso = 0.0;


    // sort electrons and muons
    sort(muonvector.begin(), muonvector.end(), muonSorter);
    sort(electronvector.begin(), electronvector.end(), electronSorter);
    sort(tauvector.begin(), tauvector.end(), tauSorter);
    
    // one or two loose muons
    if (nmuons == 1 || nmuons == 2) {
      
      pat::MuonRef muon = muonvector[0];
      mu1pid   = muon->pdgId(); 
      mu1pt    = muon->pt(); 
      mu1eta   = muon->eta(); 
      mu1phi   = muon->phi();
      mu1pfpt  = muon->pfP4().Pt();
      mu1pfeta = muon->pfP4().Eta();
      mu1pfphi = muon->pfP4().Phi();
      mu1iso   = computeMuonIso(*muon); 
      mu1idm   = (muon::isMediumMuon(*muon) ? 1 : 0);
      if (verticesH->size() > 0) 
	mu1idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);
      
      for (std::size_t i = 0; i < tightmuons.size(); i++) {
	if (muon == tightmuons[i]) 
	  mu1id = 1; // tight muon
      }
            
      // store high-pt muons that are not tight ones
      for (std::size_t i = 0; i < highptmuons.size(); i++) {
	if (muon == highptmuons[i] and mu1id != 1) 
	  mu1id = 2; // high pt muon
      }

      if (nmuons == 1) 
	wmt = sqrt(2.0 * mu1pt * t1pfmet * (1.0 - cos(deltaPhi(mu1phi, t1pfmetphi))));
    }
   
    // two loose muons
    if (nmuons == 2) {        

      pat::MuonRef muon = muonvector[1];
      mu2pid   = muon->pdgId(); 
      mu2pt    = muon->pt(); 
      mu2eta   = muon->eta(); 
      mu2phi   = muon->phi();
      mu2pfpt  = muon->pfP4().Pt();
      mu2pfeta = muon->pfP4().Eta();
      mu2pfphi = muon->pfP4().Phi();
      mu2iso   = computeMuonIso(*muon);       
      mu2idm   = (muon::isMediumMuon(*muon) ? 1 : 0);
      if (verticesH->size() > 0) 
	mu2idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);
      // check if belong to the tight / high pt collection
      for (std::size_t i = 0; i < tightmuons.size(); i++) {
	if (muon == tightmuons[i]) 
	  mu2id = 1;
      }
      
      // store high-pt muons that are not tight ones
      for (std::size_t i = 0; i < highptmuons.size(); i++) {
	if (muon == highptmuons[i] and mu2id != 1) 
	  mu2id = 2;
      }
      
      TLorentzVector mu1vec; 
      mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonvector[0]->p());
      TLorentzVector mu2vec; 
      mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muon->p());
      
      TLorentzVector zvec(mu1vec);
      zvec += mu2vec;
    
      zmass = zvec.M();
      zpt   = zvec.Pt();
      zeta  = zvec.Eta();            
      zphi  = zvec.Phi();
    }
    
    // one or two loose electrons
    if (nelectrons == 1 || nelectrons == 2) {
      pat::ElectronRef electron = electronvector[0];
      el1pid = electron->pdgId();
      el1pt  = electron->pt();
      el1eta = electron->eta();
      el1phi = electron->phi();
      for(std::size_t i = 0; i < looseelectrons.size(); i++) {
	if(electron == looseelectrons[i])
	  el1idl = 1;
      }
      
      for (std::size_t i = 0; i < tightelectrons.size(); i++) {
	if (electron == tightelectrons[i]) 
	  el1id = 1;
      }
      
      for (std::size_t i = 0; i < heepelectrons.size(); i++) {
	if (electron == heepelectrons[i] and el1id != 1) 
	  el1id = 2;
      }
              
      if (electrons.size() == 1) 
	wemt = sqrt(2.0 * el1pt * t1pfmet * (1.0 - cos(deltaPhi(el1phi, t1pfmetphi))));
    }

    // two loose electrons
    if (nelectrons == 2) {
        pat::ElectronRef electron = electronvector[1];
        el2pid = electron->pdgId();
        el2pt  = electron->pt();
        el2eta = electron->eta();
        el2phi = electron->phi();

        for (std::size_t i = 0; i < looseelectrons.size(); i++) {
            if (electron == looseelectrons[i]) el2idl = 1;
        }

        for (std::size_t i = 0; i < tightelectrons.size(); i++) {
            if (electron == tightelectrons[i]) el2id = 1;
        }

	for (std::size_t i = 0; i < heepelectrons.size(); i++) {
	  if (electron == heepelectrons[i] and el2id != 1) 
	    el2id = 2;
	}
 
        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electronvector[0]->p());
        TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electron->p());

        TLorentzVector zvec(el1vec);
        zvec += el2vec;

        zeemass = zvec.M();
        zeept   = zvec.Pt();
        zeeeta  = zvec.Eta();
        zeephi  = zvec.Phi();
    }

    ///// taus
    // one or two loose muons
    if (ntaus == 1 || ntaus == 2) {
      
      pat::TauRef tau = tauvector[0];
      tau1pid   = tau->pdgId(); 
      tau1pt    = tau->pt(); 
      tau1eta   = tau->eta(); 
      tau1phi   = tau->phi();
      tau1m     = tau->mass();
      tau1iso   = tau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");

      if (ntaus == 1) 
	wtmt = sqrt(2.0 * tau1pt * t1pfmet * (1.0 - cos(deltaPhi(tau1phi, t1pfmetphi))));
    }
   
    // two loose muons
    if (ntaus == 2) {        

      pat::TauRef tau = tauvector[1];
      tau2pid   = tau->pdgId(); 
      tau2pt    = tau->pt(); 
      tau2eta   = tau->eta(); 
      tau2phi   = tau->phi();
      tau2m     = tau->mass();
      tau2iso   = tau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT");
      
      TLorentzVector tau1vec; 
      tau1vec.SetPtEtaPhiE(tau1pt, tau1eta, tau1phi, tauvector[0]->p());
      TLorentzVector tau2vec; 
      tau2vec.SetPtEtaPhiE(tau2pt, tau2eta, tau2phi, tau->p());
      
      TLorentzVector zvec(tau1vec);
      zvec += tau2vec;
    
      zttmass = zvec.M();
      zttpt   = zvec.Pt();
      ztteta  = zvec.Eta();            
      zttphi  = zvec.Phi();
    }
    

    // one electron and one muon (ttbar DF fully leptonic)
    if (nmuons == 1 && nelectrons == 1) {
      TLorentzVector mu1vec; 
      mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonvector[0]->p());
      TLorentzVector el1vec; 
      el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electronvector[0]->p());
      
      TLorentzVector emuvec(mu1vec);
      emuvec += el1vec;
      
      emumass = emuvec.M();
      emupt   = emuvec.Pt();
      emueta  = emuvec.Eta();
      emuphi  = emuvec.Phi();
    } 

    // one electron and one tau (ttbar DF fully leptonic)
    if (ntaus == 1 && nelectrons == 1) {
      TLorentzVector tau1vec; 
      tau1vec.SetPtEtaPhiE(tau1pt, tau1eta, tau1phi, tauvector[0]->p());
      TLorentzVector el1vec; 
      el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi,  electronvector[0]->p());
      
      TLorentzVector emuvec(tau1vec);
      emuvec += el1vec;
      
      tauemass = emuvec.M();
      tauept   = emuvec.Pt();
      taueeta  = emuvec.Eta();
      tauephi  = emuvec.Phi();
    } 

    // one tau and one muon (ttbar DF fully leptonic)
    if (nmuons == 1 && ntaus == 1) {
      TLorentzVector mu1vec; 
      mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonvector[0]->p());
      TLorentzVector tau1vec; 
      tau1vec.SetPtEtaPhiE(tau1pt, tau1eta, tau1phi, tauvector[0]->p());
      
      TLorentzVector emuvec(tau1vec);
      emuvec += mu1vec;
      
      taumumass = emuvec.M();
      taumupt   = emuvec.Pt();
      taumueta  = emuvec.Eta();
      taumuphi  = emuvec.Phi();
    } 
    
    // Photon information
    phidl    = 0; phidm    = 0; phidt    = 0; phidh    = 0;
    phpt     = 0; pheta    = 0; phphi    = 0;

    int hardestPhotonIndex = -1;
    double hardestPhotonPt = 0.0;

    if(photonsH.isValid() and photonLooseIdH.isValid() and mediumPhotonsH.isValid() and tightPhotonsH.isValid() and photonHighPtIdH.isValid()){
      
      for (size_t i = 0; i < photons.size(); i++) {
        if (photons[i]->pt() > hardestPhotonPt) {
	  hardestPhotonIndex = i;
	  hardestPhotonPt = photons[i]->pt();
        }
      }

      nphotons = photons.size();
      
      if (hardestPhotonIndex >= 0) {
	phidl   = ((*photonLooseIdH )[photons[hardestPhotonIndex]] ? 1 : 0);
	phidh   = ((*photonHighPtIdH)[photons[hardestPhotonIndex]] ? 1 : 0);

	for(size_t i = 0; i < mediumphotons.size(); i++){
	  if(photons[hardestPhotonIndex] == mediumphotons[i])
	    phidm = 1;
	}
       
	for(size_t i = 0; i < tightphotons.size(); i++){
	  if(tightphotons[hardestPhotonIndex] == tightphotons[i])
	    phidt = 1;
	}

	phpt    = photons[hardestPhotonIndex]->pt();
	pheta   = photons[hardestPhotonIndex]->eta();
	phphi   = photons[hardestPhotonIndex]->phi();
      }
    }
    
    int hardestPhotonPurityIndex = -1;
    double hardestPhotonPurityPt = 0.0;

    if(addPhotonPurity and not isTriggerTree){
      phPuritypt     = 0.0;
      phPurityeta    = 0.0;
      phPurityphi    = 0.0;
      phPurityPHiso    = 0.0;
      phPurityRND04PHiso    = 0.0;
      phPurityRND08PHiso    = 0.0;
      phPurityCHiso    = 0.0;
      phPurityRND08CHiso    = 0.0;
      phPurityRND04CHiso    = 0.0;
      phNHiso          = 0.0;
      phPurityNHiso    = 0.0;
      phPuritysieie    = 0.0;
      phPurityhoe    = 0.0;
      phPurityElectronVeto = 0;
      nphotonsPurity = photonsPurityH->size();
      phPurityEA = 0.;
      phPurityEAEGamma = 0.;
      
      for (size_t i = 0; i < tightphotonsPurity.size(); i++) {
	if (tightphotonsPurity[i]->pt() > hardestPhotonPurityPt) {
	  hardestPhotonPurityIndex = i;
	  hardestPhotonPurityPt = tightphotonsPurity[i]->pt();
	}
      }
      
      if (hardestPhotonPurityIndex >= 0) {
	phPuritypt    = tightphotonsPurity[hardestPhotonPurityIndex]->pt();
	phPurityeta   = tightphotonsPurity[hardestPhotonPurityIndex]->eta();
	phPurityphi   = tightphotonsPurity[hardestPhotonPurityIndex]->phi();

	phPHiso       = tightphotonsPurity[hardestPhotonPurityIndex]->photonIso()-rho*getGammaEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta());
	phCHiso       = tightphotonsPurity[hardestPhotonPurityIndex]->chargedHadronIso()-rho*getChargedHadronEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta());
	phNHiso       = tightphotonsPurity[hardestPhotonPurityIndex]->neutralHadronIso()-rho*getNeutralHadronEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta());

	phPurityPHiso      = (*photonPHisoH)[tightphotonsPurity[hardestPhotonPurityIndex]];
	phPurityCHiso      = (*photonCHisoH)[tightphotonsPurity[hardestPhotonPurityIndex]]-rho*getChargedHadronEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta()); 
	phPurityNHiso     = (*photonNHisoH)[tightphotonsPurity[hardestPhotonPurityIndex]]-rho*getNeutralHadronEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta()); 	

	phPurityRND04CHiso = (*rndchhadiso04H)[tightphotonsPurity[hardestPhotonPurityIndex]];
	phPurityRND04PHiso = (*rndgammaiso04H)[tightphotonsPurity[hardestPhotonPurityIndex]];
	phPurityRND08PHiso = (*rndgammaiso08H)[tightphotonsPurity[hardestPhotonPurityIndex]];
	phPurityRND08CHiso = (*rndchhadiso08H)[tightphotonsPurity[hardestPhotonPurityIndex]];

	phPuritysieie     = (*photonsieieH)[tightphotonsPurity[hardestPhotonPurityIndex]];
	phPurityElectronVeto   = tightphotonsPurity[hardestPhotonPurityIndex]->passElectronVeto();
	phPurityhoe       = tightphotonsPurity[hardestPhotonPurityIndex]->hadTowOverEm();
	phPurityEAEGamma  = getGammaEAForPhotonIso(tightphotonsPurity[hardestPhotonPurityIndex]->eta());
	if(abs(tightphotonsPurity[hardestPhotonPurityIndex]->eta()) < 0.9) phPurityEA = 0.1271;
	if(abs(tightphotonsPurity[hardestPhotonPurityIndex]->eta()) > 0.9 && abs(tightphotonsPurity[hardestPhotonPurityIndex]->eta()) < 1.4442 ) phPurityEA = 0.1101;
      }
    }

    // Substructure CHS
    if(addSubstructureCHS and not isTriggerTree){      
      //sort collection to make sure it is ordered
      vector<pat::JetRef> jetsBoosted;
      if(boostedJetsH.isValid()){
	for (auto jets_iter = boostedJetsH->begin(); jets_iter != boostedJetsH->end(); ++jets_iter) {
	  //clean from leptons                                                                                                                                                  
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }
	  }
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }
	  }
	  
	  if (skipjet) continue;
	  
	  // apply jet id                                                                                                                                                 
	  bool passjetid = applyJetID(*jets_iter,jetidwp);
	  if (!passjetid)
	    continue;
	  
	  pat::JetRef jetref(boostedJetsH, jets_iter - boostedJetsH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    jetsBoosted.push_back(jetref);
	}
	
	// sort in pt
	if(jetsBoosted.size() > 0)
	  sort(jetsBoosted.begin(), jetsBoosted.end(), jetSorter);


	boostedJetpt    .clear(); boostedJeteta  .clear();  
	boostedJetphi   .clear(); boostedJetm    .clear();
	boostedJetGenpt .clear(); boostedJetGenm .clear();
	boostedJetGeneta.clear(); boostedJetGenphi .clear();

	boostedJettau1  .clear(); boostedJettau2 .clear();
	boostedJettau3  .clear(); boostedJettau4 .clear();
	boostedJetGentau1.clear(); boostedJetGentau2 .clear();
	boostedJetGentau3.clear(); boostedJetGentau4 .clear();
	boostedJetHFlav.clear(); boostedJetPFlav.clear(); boostedJetQGL.clear(); 
	boostedJetBtag .clear(), boostedJetDoubleBtag .clear();
	
	prunedJetpt     .clear(); prunedJetm    .clear(); prunedJetGenpt.clear(); prunedJetGenm .clear();
	prunedJeteta    .clear(); prunedJetphi  .clear(); prunedJetGeneta.clear(); prunedJetGenphi.clear();
	prunedJetm_v2 .clear();   prunedJetpt_v2 .clear(); prunedJeteta_v2 .clear(); prunedJetphi_v2 .clear(); 
	prunedJetHFlav  .clear(); prunedJetPFlav.clear(); prunedJetQGL  .clear(); prunedJetBtag .clear(); prunedJetDoubleBtag.clear();
	prunedJetptraw  .clear(); prunedJetmraw .clear();

	softDropJetpt   .clear(); softDropJetm.clear();    softDropJetGenpt.clear();  softDropJetGenm.clear(); 	
	softDropJetphi  .clear(); softDropJeteta.clear();  softDropJetGenphi.clear(); softDropJetGeneta.clear();
	softDropJetm_v2.clear(); softDropJetpt_v2.clear(); softDropJetphi_v2 .clear(); softDropJeteta_v2. clear();
	softDropJetHFlav.clear(); softDropJetPFlav.clear();softDropJetQGL.clear();    softDropJetBtag.clear(); softDropJetDoubleBtag.clear();
	softDropJetptraw.clear(); softDropJetmraw.clear();
	
	prunedSubJetpt_1 .clear(); prunedSubJetm_1  .clear(); prunedSubJetphi_1 .clear(); prunedSubJeteta_1 .clear();
	prunedSubJetHFlav_1 .clear(); prunedSubJetQGL_1 .clear(); prunedSubJetBtag_1 .clear();
	prunedSubJetBtagSF_1.clear(); prunedSubJetBtagSFUp_1.clear(); prunedSubJetBtagSFDown_1.clear();
	prunedSubJetGenpt_1 .clear(); prunedSubJetGenm_1 .clear(); prunedSubJetPFlav_1 .clear();
	prunedSubJetGenphi_1 .clear(); prunedSubJetGeneta_1 .clear(); 
	prunedSubJetptraw_1 .clear(); prunedSubJetmraw_1  .clear(); 
	
	prunedSubJetpt_2 .clear(); prunedSubJetm_2  .clear(); prunedSubJetphi_2 .clear(); prunedSubJeteta_2 .clear();
	prunedSubJetHFlav_2 .clear(); prunedSubJetQGL_2 .clear(); prunedSubJetBtag_2 .clear();  prunedSubJetPFlav_2 .clear();
	prunedSubJetBtagSF_2.clear(); prunedSubJetBtagSFUp_2.clear(); prunedSubJetBtagSFDown_2.clear();
	prunedSubJetGenpt_2 .clear(); prunedSubJetGenm_2 .clear();
	prunedSubJetGeneta_2 .clear(); prunedSubJetGenphi_2 .clear();
	prunedSubJetptraw_2 .clear(); prunedSubJetmraw_2  .clear(); 

	softDropSubJetpt_1 .clear(); softDropSubJetm_1  .clear(); softDropSubJetphi_1 .clear(); softDropSubJeteta_1 .clear();
	softDropSubJetHFlav_1 .clear(); softDropSubJetQGL_1 .clear(); softDropSubJetBtag_1 .clear(); softDropSubJetPFlav_1 .clear();
	softDropSubJetBtagSF_1.clear(); softDropSubJetBtagSFUp_1.clear(); softDropSubJetBtagSFDown_1.clear();
	softDropSubJetGenpt_1 .clear(); softDropSubJetGenm_1 .clear();
	softDropSubJetGeneta_1 .clear(); softDropSubJetGenphi_1 .clear();	
	softDropSubJetptraw_1 .clear(); softDropSubJetmraw_1  .clear(); 

	softDropSubJetpt_2 .clear(); softDropSubJetm_2  .clear(); softDropSubJetphi_2 .clear(); softDropSubJeteta_2 .clear();
	softDropSubJetHFlav_2 .clear(); softDropSubJetQGL_2 .clear(); softDropSubJetBtag_2 .clear(); softDropSubJetPFlav_2 .clear();
	softDropSubJetBtagSF_2.clear(); softDropSubJetBtagSFUp_2.clear(); softDropSubJetBtagSFDown_2.clear();
	softDropSubJetGenpt_2 .clear(); softDropSubJetGenm_2 .clear();
	softDropSubJetGeneta_2 .clear(); softDropSubJetGenphi_2 .clear();
	softDropSubJetptraw_2 .clear(); softDropSubJetmraw_2  .clear(); 
	
	for(size_t i = 0; i < jetsBoosted.size(); i++){
	
	  boostedJetpt  .push_back( jetsBoosted[i]->pt());
	  boostedJeteta .push_back( jetsBoosted[i]->eta());
	  boostedJetphi .push_back( jetsBoosted[i]->phi());
	  boostedJetm   .push_back( jetsBoosted[i]->mass());	
	  boostedJetBtag .push_back( jetsBoosted[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	  boostedJetDoubleBtag .push_back( jetsBoosted[i]->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	  // N-jettiness
	  if(jetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau1"))
	    boostedJettau1 .push_back( jetsBoosted[i]->userFloat("Njettiness"+boostedJetsCHSLabel+":tau1"));
	  
	  if(jetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau2"))
	    boostedJettau2 .push_back( jetsBoosted[i]->userFloat("Njettiness"+boostedJetsCHSLabel+":tau2"));
	  
	  if(jetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau3"))
	    boostedJettau3 .push_back( jetsBoosted[i]->userFloat("Njettiness"+boostedJetsCHSLabel+":tau3"));	

	  if(jetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau4"))
	    boostedJettau4 .push_back( jetsBoosted[i]->userFloat("Njettiness"+boostedJetsCHSLabel+":tau4"));

	  // Gen n-subjettiness
	  if(isMC){
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau1"))
	      boostedJetGentau1 .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau1"));

	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau2"))
	      boostedJetGentau2 .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau2"));

	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau3"))
	      boostedJetGentau3 .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau3"));

	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau4"))
	      boostedJetGentau4 .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenNjettinessMatched:tau4"));
	  }

	  /*
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf1"))
	    boostedJetecf1 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf1"));
	  
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf2"))
	    boostedJetecf2 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf2"));
	  
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf3"))
	    boostedJetecf3 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf3"));
	  */

	  // gen jets
	  if(isMC){
	    if(jetsBoosted[i]->genJet()){ // gen AK8 jet
	      boostedJetGenpt .push_back( jetsBoosted[i]->genJet()->pt());
	      boostedJetGenm  .push_back( jetsBoosted[i]->genJet()->mass());
	      boostedJetGeneta  .push_back( jetsBoosted[i]->genJet()->eta());
	      boostedJetGenphi  .push_back( jetsBoosted[i]->genJet()->phi());	      
	    }
	    boostedJetHFlav .push_back( jetsBoosted[i]->hadronFlavour());
	    boostedJetPFlav .push_back( jetsBoosted[i]->partonFlavour());	  
	  }
	  
	  // QGL 
	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"QGL:qgLikelihood"))
	    boostedJetQGL .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"QGL:qgLikelihood"));

	 	  
	  // pruned matched jet

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:mass"))
	    prunedJetm .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:mass"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:eta"))
	    prunedJeteta .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:eta"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:phi"))
	    prunedJetphi .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:phi"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pt"))
	    prunedJetpt .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:pt"));
	
	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedQGLMatched:qgLikelihood"))
	    prunedJetQGL .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedQGLMatched:qgLikelihood"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    prunedJetBtag .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    prunedJetDoubleBtag .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"))
	    prunedJetptraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass")){
	    prunedJetmraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass"));
	    // apply correction by hand from uncorrected variables
	    if(jetsBoosted[i]->availableJECSets().size() > 1 and 
	       jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:raweta")	and
	       jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawphi") and
	       jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt")){
	      
	      TLorentzVector correctedP4;
	      correctedP4.SetPtEtaPhiM(jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:raweta"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawphi"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass")
				       );
	      correctedP4 *= 1./jetsBoosted[i]->jecFactor("Uncorrected","none",jetsBoosted[i]->availableJECSets().at(1)); // apply AK8 corrections
	      prunedJetm_v2 .push_back(correctedP4.M());
	      prunedJetpt_v2 .push_back(correctedP4.Pt());
	      prunedJeteta_v2 .push_back(correctedP4.Eta());
	      prunedJetphi_v2 .push_back(correctedP4.Phi());
	    }
	    else{
	      prunedJetm_v2 .push_back(0.);
	      prunedJetpt_v2 .push_back(0.);
	      prunedJeteta_v2 .push_back(0.);
	      prunedJetphi_v2 .push_back(0.);
	    }
	  }

	  if(isMC){
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:hadronFlavour"))
	      prunedJetHFlav  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:hadronFlavour"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:partonFlavour"))
	      prunedJetPFlav .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:partonFlavour"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genMass"))
	      prunedJetGenm  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:genMass"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genPt"))
	      prunedJetGenpt  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:genPt"));	
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genEta"))
	      prunedJetGeneta  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:genEta"));	
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genPhi"))
	      prunedJetGenphi  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:genPhi"));	
	  }	

	  // soft drop matched jet
	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:mass"))
	    softDropJetm .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:mass"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pt"))
	    softDropJetpt .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pt"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:phi"))
	    softDropJetphi .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:phi"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:eta"))
	    softDropJeteta .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:eta"));
	
	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropQGLMatched:qgLikelihood"))
	    softDropJetQGL .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropQGLMatched:qgLikelihood"));
	  
	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    softDropJetBtag .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    softDropJetDoubleBtag .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw"))
	    softDropJetptraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass")){
	    softDropJetmraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass"));

            if(jetsBoosted[i]->availableJECSets().size()>1 and
               jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:raweta") and
	       jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:rawphi") and 
	       jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:rawpt")){

              TLorentzVector correctedP4;
	      correctedP4.SetPtEtaPhiM(jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawpt"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:raweta"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawphi"),
				       jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass")
                                       );
	      correctedP4 *= 1./jetsBoosted[i]->jecFactor("Uncorrected","none",jetsBoosted[i]->availableJECSets().at(1));
              softDropJetm_v2 .push_back(correctedP4.M());
              softDropJetpt_v2 .push_back(correctedP4.Pt());
              softDropJeteta_v2 .push_back(correctedP4.Eta());
              softDropJetphi_v2 .push_back(correctedP4.Phi());
            }
            else{
              softDropJetm_v2 .push_back(0.);
              softDropJetpt_v2 .push_back(0.);
              softDropJeteta_v2 .push_back(0.);
              softDropJetphi_v2 .push_back(0.);
            }
	  }

	  if(isMC){
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:hadronFlavour"))
	      softDropJetHFlav  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:hadronFlavour"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:partonFlavour"))
	      softDropJetPFlav .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:partonFlavour"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genMass"))
	      softDropJetGenm  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genMass"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genPt"))
	      softDropJetGenpt  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genPt"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genEta"))
	      softDropJetGeneta  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genEta"));
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genPhi"))
	      softDropJetGenphi  .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genPhi"));
	  }	
	  
	  // sub-jets pruned 
	  if(jetsBoosted[i]->hasSubjets("Pruned")){

	    pat::JetPtrCollection subjets = jetsBoosted[i]->subjets("Pruned");
	    if(subjets.size() > 0 ){
	      prunedSubJetpt_1  .push_back( subjets[0]->pt()); 
	      prunedSubJetm_1   .push_back( subjets[0]->mass()); 
	      prunedSubJetphi_1 .push_back( subjets[0]->phi()); 
	      prunedSubJeteta_1 .push_back( subjets[0]->eta());
	      prunedSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	      
	      prunedSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	      prunedSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	      
	      if(subjets[0]->hasUserFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"))
		prunedSubJetQGL_1 .push_back( subjets[0]->userFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"));
	      	      
	      if(isMC){
		prunedSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
		prunedSubJetPFlav_1 .push_back( subjets[0]->partonFlavour()); 
		if(subjets[0]->genJet()){
		  prunedSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		  prunedSubJetGenm_1  .push_back( subjets[0]->genJet()->mass()); 
		  prunedSubJetGeneta_1  .push_back( subjets[0]->genJet()->eta()); 
		  prunedSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi()); 
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets[0],"SubCSV",prunedSubJetBtagSF_1,prunedSubJetBtagSFUp_1,prunedSubJetBtagSFDown_1);
	      }
	    }
	  
	    if(subjets.size() > 1 ){
	      prunedSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	      prunedSubJetm_2   .push_back( subjets.at(1)->mass()); 
	      prunedSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	      prunedSubJeteta_2 .push_back( subjets.at(1)->eta());
	      prunedSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      prunedSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	      prunedSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"))
		prunedSubJetQGL_2 .push_back( subjets.at(1)->userFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"));
	      

	      if(isMC){
		prunedSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
		prunedSubJetPFlav_2 .push_back( subjets.at(1)->partonFlavour()); 
		if(subjets.at(1)->genJet()){
		  prunedSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		  prunedSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		  prunedSubJetGeneta_2 .push_back( subjets.at(1)->genJet()->eta());
		  prunedSubJetGenphi_2 .push_back( subjets.at(1)->genJet()->phi());
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(1),"SubCSV",prunedSubJetBtagSF_2,prunedSubJetBtagSFUp_2,prunedSubJetBtagSFDown_2);
	      }
	    }
	  }

	  // sub-jets soft drop 
	  if(jetsBoosted[i]->hasSubjets("SoftDrop")){
	    pat::JetPtrCollection subjets = jetsBoosted[i]->subjets("SoftDrop");
	    if(subjets.size() > 0 ){
	      softDropSubJetpt_1  .push_back( subjets[0]->pt()); 
	      softDropSubJetm_1   .push_back( subjets[0]->mass()); 
	      softDropSubJetphi_1 .push_back( subjets[0]->phi()); 
	      softDropSubJeteta_1 .push_back( subjets[0]->eta());
	      softDropSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      prunedSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	      prunedSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	    
	      if(subjets[0]->hasUserFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		softDropSubJetQGL_1 .push_back( subjets[0]->userFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"));
	    
	    
	      if(isMC){
		softDropSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
		softDropSubJetPFlav_1 .push_back( subjets[0]->partonFlavour()); 
		if(subjets[0]->genJet()){
		  softDropSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		  softDropSubJetGenm_1  .push_back( subjets[0]->genJet()->mass()); 
		  softDropSubJetGeneta_1  .push_back( subjets[0]->genJet()->eta()); 
		  softDropSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi()); 
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets[0],"SubCSV",softDropSubJetBtagSF_1,softDropSubJetBtagSFUp_1,softDropSubJetBtagSFDown_1);
	      }
	    }
	    
	    if(subjets.size() > 1 ){
	      softDropSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	      softDropSubJetm_2   .push_back( subjets.at(1)->mass()); 
	      softDropSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	      softDropSubJeteta_2 .push_back( subjets.at(1)->eta());
	      softDropSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      prunedSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	      prunedSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		softDropSubJetQGL_2 .push_back( subjets.at(1)->userFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"));
	      
	      if(isMC){
		softDropSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
		softDropSubJetPFlav_2 .push_back( subjets.at(1)->partonFlavour()); 
		if(subjets.at(1)->genJet()){
		  softDropSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		  softDropSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		  softDropSubJetGeneta_2  .push_back( subjets.at(1)->genJet()->eta());
		  softDropSubJetGenphi_2  .push_back( subjets.at(1)->genJet()->phi());
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(1),"SubCSV",softDropSubJetBtagSF_2,softDropSubJetBtagSFUp_2,softDropSubJetBtagSFDown_2);
	      }
	    }	  
	  }
	}
      }
    }

    // Substructure Puppi
    if(addSubstructurePuppi and not isTriggerTree){      
      //sort collection to make sure it is ordered
      vector<pat::JetRef> puppiJetsBoosted;
      
      if(boostedPuppiJetsH.isValid()){
	for (auto jets_iter = boostedPuppiJetsH->begin(); jets_iter != boostedPuppiJetsH->end(); ++jets_iter) {
	  
	  //clean from leptons                                                                                                                                                
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }
	  }
	  
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
		skipjet = true;
	    }	  
	  }
	  
	  if (skipjet) continue;

	  // apply jet id                                                                                                                                                    
	  bool passjetid = applyJetID(*jets_iter,jetidwp);
	  if (!passjetid)
	    continue;
	  
	  pat::JetRef jetref(boostedPuppiJetsH, jets_iter - boostedPuppiJetsH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    puppiJetsBoosted.push_back(jetref);
	}	
	// sort in pt
	if(puppiJetsBoosted.size() > 0)
	  sort(puppiJetsBoosted.begin(), puppiJetsBoosted.end(), jetSorter);
	
	boostedPuppiJetpt    .clear(); boostedPuppiJeteta  .clear();  boostedPuppiJetphi   .clear(); boostedPuppiJetm    .clear();
	boostedPuppiJetGenpt .clear(); boostedPuppiJetGenm .clear();  boostedPuppiJetGeneta .clear(); boostedPuppiJetGenphi .clear();
	boostedPuppiJettau1  .clear(); boostedPuppiJettau2 .clear();  boostedPuppiJettau3  .clear(); boostedPuppiJettau4 .clear();
	boostedPuppiJetGentau1  .clear(); boostedPuppiJetGentau2 .clear(); boostedPuppiJetGentau3  .clear(); boostedPuppiJetGentau4 .clear();
	boostedPuppiJetHFlav .clear(); boostedPuppiJetPFlav .clear(); boostedPuppiJetQGL  .clear(); 
	boostedPuppiJetBtag .clear(); boostedPuppiJetDoubleBtag .clear();
	
	prunedPuppiJetpt     .clear(); prunedPuppiJetm      .clear(); prunedPuppiJetGenpt .clear(); prunedPuppiJetGenm  .clear();
	prunedPuppiJetm_v2   .clear(); prunedPuppiJetpt_v2   .clear();prunedPuppiJeteta_v2   .clear(); prunedPuppiJetphi_v2   .clear();
	prunedPuppiJeteta     .clear(); prunedPuppiJetphi   .clear(); prunedPuppiJetGeneta .clear(); prunedPuppiJetGenphi  .clear();
	prunedPuppiJetHFlav  .clear(); prunedPuppiJetPFlav  .clear(); prunedPuppiJetQGL   .clear(); prunedPuppiJetBtag  .clear();
	prunedPuppiJetDoubleBtag  .clear(); prunedPuppiJetmraw .clear(); prunedPuppiJetptraw .clear();
	

	softDropPuppiJetpt    .clear(); softDropPuppiJetm .clear(); softDropPuppiJetGenpt .clear(); softDropPuppiJetGenm .clear(); 
	softDropPuppiJetm_v2.clear();   softDropPuppiJetpt_v2.clear(); softDropPuppiJeteta_v2.clear(); softDropPuppiJetphi_v2.clear();
	softDropPuppiJeteta    .clear(); softDropPuppiJetphi .clear(); softDropPuppiJetGeneta .clear(); softDropPuppiJetGenphi .clear(); 
	softDropPuppiJetHFlav .clear(); softDropPuppiJetPFlav .clear(); softDropPuppiJetQGL .clear(); softDropPuppiJetBtag .clear();
	softDropPuppiJetDoubleBtag .clear(); softDropPuppiJetmraw .clear(); softDropPuppiJetptraw .clear();
	
	prunedPuppiSubJetpt_1 .clear(); prunedPuppiSubJetm_1  .clear(); prunedPuppiSubJetphi_1 .clear(); prunedPuppiSubJeteta_1 .clear();
	prunedPuppiSubJetHFlav_1 .clear(); prunedPuppiSubJetPFlav_1 .clear(); prunedPuppiSubJetQGL_1 .clear(); prunedPuppiSubJetBtag_1 .clear();
	prunedPuppiSubJetGenpt_1 .clear(); prunedPuppiSubJetGenm_1 .clear(); prunedPuppiSubJetGenphi_1 .clear(); prunedPuppiSubJetGeneta_1 .clear();
	prunedPuppiSubJetptraw_1 .clear(); prunedPuppiSubJetmraw_1  .clear(); 
	prunedPuppiSubJetBtagSF_1.clear(); prunedPuppiSubJetBtagSFUp_1.clear(); prunedPuppiSubJetBtagSFDown_1.clear();
      

	prunedPuppiSubJetpt_2 .clear(); prunedPuppiSubJetm_2  .clear(); prunedPuppiSubJetphi_2 .clear(); prunedPuppiSubJeteta_2 .clear();
	prunedPuppiSubJetHFlav_2 .clear();prunedPuppiSubJetPFlav_2 .clear(); prunedPuppiSubJetQGL_2 .clear(); prunedPuppiSubJetBtag_2 .clear();
        prunedPuppiSubJetGenpt_2 .clear(); prunedPuppiSubJetGenm_2 .clear(); prunedPuppiSubJetGenphi_2 .clear(); prunedPuppiSubJetGeneta_2 .clear();
	prunedPuppiSubJetptraw_2 .clear(); prunedPuppiSubJetmraw_2  .clear(); 
	prunedPuppiSubJetBtagSF_2.clear(); prunedPuppiSubJetBtagSFUp_2.clear(); prunedPuppiSubJetBtagSFDown_2.clear();
	
	softDropPuppiSubJetpt_1 .clear(); softDropPuppiSubJetm_1  .clear(); softDropPuppiSubJetphi_1 .clear(); softDropPuppiSubJeteta_1 .clear();
	softDropPuppiSubJetHFlav_1 .clear(); softDropPuppiSubJetQGL_1 .clear(); softDropPuppiSubJetBtag_1 .clear(); softDropPuppiSubJetPFlav_1 .clear();
	softDropPuppiSubJetGenpt_1 .clear(); softDropPuppiSubJetGenm_1 .clear(); softDropPuppiSubJetGenphi_1 .clear(); softDropPuppiSubJetGeneta_1 .clear();
	softDropPuppiSubJetptraw_1 .clear(); softDropPuppiSubJetmraw_1  .clear(); 
	softDropPuppiSubJetBtagSF_1.clear(); softDropPuppiSubJetBtagSFUp_1.clear(); softDropPuppiSubJetBtagSFDown_1.clear();

	softDropPuppiSubJetpt_2 .clear(); softDropPuppiSubJetm_2  .clear(); softDropPuppiSubJetphi_2 .clear(); softDropPuppiSubJeteta_2 .clear();
	softDropPuppiSubJetHFlav_2 .clear(); softDropPuppiSubJetQGL_2 .clear(); softDropPuppiSubJetBtag_2 .clear(); softDropPuppiSubJetPFlav_2 .clear();
	softDropPuppiSubJetGenpt_2 .clear(); softDropPuppiSubJetGenm_2 .clear(); softDropPuppiSubJetGenphi_2 .clear(); softDropPuppiSubJetGeneta_2 .clear();
	softDropPuppiSubJetptraw_2 .clear(); softDropPuppiSubJetmraw_2  .clear(); 
	softDropPuppiSubJetBtagSF_2.clear(); softDropPuppiSubJetBtagSFUp_2.clear(); softDropPuppiSubJetBtagSFDown_2.clear();
	
	for(size_t i = 0; i < puppiJetsBoosted.size(); i++){
	  
	  boostedPuppiJetpt  .push_back( puppiJetsBoosted[i]->pt());
	  boostedPuppiJeteta .push_back( puppiJetsBoosted[i]->eta());
	  boostedPuppiJetphi .push_back( puppiJetsBoosted[i]->phi());
	  boostedPuppiJetm   .push_back( puppiJetsBoosted[i]->mass());
	  boostedPuppiJetBtag .push_back( puppiJetsBoosted[i]->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
	  boostedPuppiJetDoubleBtag .push_back( puppiJetsBoosted[i]->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	  
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau1"))
	    boostedPuppiJettau1 .push_back( puppiJetsBoosted[i]->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau1"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau2"))
	    boostedPuppiJettau2 .push_back( puppiJetsBoosted[i]->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau2"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau3"))
	    boostedPuppiJettau3 .push_back( puppiJetsBoosted[i]->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau3"));	
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau4"))
	    boostedPuppiJettau4 .push_back( puppiJetsBoosted[i]->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau4"));
	  /*	  	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf1")){
	    boostedPuppiJetecf1 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf1"));
	  }
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf2")){
	    boostedPuppiJetecf2 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf2"));
	  }
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf3")){
	    boostedPuppiJetecf3 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf3"));
	  }
	  */
	  if(isMC){
	    if(puppiJetsBoosted[i]->genJet()){ // gen AK8 jet
	      boostedPuppiJetGenpt .push_back( puppiJetsBoosted[i]->genJet()->pt());
	      boostedPuppiJetGenm  .push_back( puppiJetsBoosted[i]->genJet()->mass());
	      boostedPuppiJetGeneta  .push_back( puppiJetsBoosted[i]->genJet()->eta());
	      boostedPuppiJetGenphi  .push_back( puppiJetsBoosted[i]->genJet()->phi());
	    }
	    boostedPuppiJetHFlav .push_back( puppiJetsBoosted[i]->hadronFlavour());
	    boostedPuppiJetPFlav .push_back( puppiJetsBoosted[i]->partonFlavour());	  
	  }
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"QGL:qgLikelihood"))
	    boostedPuppiJetQGL .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"QGL:qgLikelihood"));
	  
	  if(isMC){
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau1"))
	      boostedPuppiJetGentau1 .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau1"));
	    
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau2"))
	      boostedPuppiJetGentau2 .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau2"));
	    
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau3"))
	      boostedPuppiJetGentau3 .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau3"));
	    
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau4"))
	      boostedPuppiJetGentau4 .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenNjettinessMatched:tau4"));
	  }
	  
	  
	  // pruned matched jet
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:mass"))
	    prunedPuppiJetm .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:mass"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pt"))
	    prunedPuppiJetpt .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pt"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:eta"))
	    prunedPuppiJeteta .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:eta"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:phi"))
	    prunedPuppiJetphi .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:phi"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedQGLMatched:qgLikelihood"))
	    prunedPuppiJetQGL .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedQGLMatched:qgLikelihood"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    prunedPuppiJetBtag .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    prunedPuppiJetDoubleBtag .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"))
	    prunedPuppiJetptraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass")){
	    prunedPuppiJetmraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass"));

            if(puppiJetsBoosted[i]->availableJECSets().size()>1 and
               puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:raweta") and
               puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawphi") and 
	       puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt")){

              TLorentzVector correctedP4;
              correctedP4.SetPtEtaPhiM(puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:raweta"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawphi"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass")
                                       );
              correctedP4 *= 1./puppiJetsBoosted[i]->jecFactor("Uncorrected","none",puppiJetsBoosted[i]->availableJECSets().at(1));
              prunedPuppiJetm_v2 .push_back(correctedP4.M());
              prunedPuppiJetpt_v2 .push_back(correctedP4.Pt());
              prunedPuppiJeteta_v2 .push_back(correctedP4.Eta());
              prunedPuppiJetphi_v2 .push_back(correctedP4.Phi());
            }
            else{
              prunedPuppiJetm_v2 .push_back(0.);
              prunedPuppiJetpt_v2 .push_back(0.);
              prunedPuppiJeteta_v2 .push_back(0.);
              prunedPuppiJetphi_v2 .push_back(0.);
            }
	  }
	  
	  if(isMC){
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:hadronFlavour"))
	      prunedPuppiJetHFlav  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:hadronFlavour"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:partonFlavour"))
	     prunedPuppiJetPFlav .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:partonFlavour"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genMass"))
	      prunedPuppiJetGenm  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genMass"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genPt"))
	      prunedPuppiJetGenpt  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genPt"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genEta"))
	      prunedPuppiJetGeneta  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genEta"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genPhi"))
	      prunedPuppiJetGenphi  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genPhi"));
	  }
      
	  
	  // soft drop matched jet
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:mass"))
	    softDropPuppiJetm .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:mass"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pt"))
	    softDropPuppiJetpt .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pt"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:eta"))
	    softDropPuppiJeteta .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:eta"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:phi"))
	    softDropPuppiJetphi .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:phi"));
	
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropQGLMatched:qgLikelihood"))
	    softDropPuppiJetQGL .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropQGLMatched:qgLikelihood"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    softDropPuppiJetBtag .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    softDropPuppiJetDoubleBtag .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"))
	    softDropPuppiJetptraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass")){
	    softDropPuppiJetmraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass"));

            if(puppiJetsBoosted[i]->availableJECSets().size()>1 and
               puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:raweta") and
               puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawphi") and
               puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt")){

              TLorentzVector correctedP4;
              correctedP4.SetPtEtaPhiM(puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:raweta"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawphi"),
                                       puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass")
                                       );
              correctedP4 *= 1./puppiJetsBoosted[i]->jecFactor("Uncorrected","none",puppiJetsBoosted[i]->availableJECSets().at(1));
              softDropPuppiJetm_v2 .push_back(correctedP4.M());
              softDropPuppiJetpt_v2 .push_back(correctedP4.Pt());
              softDropPuppiJeteta_v2 .push_back(correctedP4.Eta());
              softDropPuppiJetphi_v2 .push_back(correctedP4.Phi());
            }
            else{
              softDropPuppiJetm_v2 .push_back(0.);
              softDropPuppiJetpt_v2 .push_back(0.);
              softDropPuppiJeteta_v2 .push_back(0.);
              softDropPuppiJetphi_v2 .push_back(0.);
            }

	  }

	  if(isMC){
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:hadronFlavour"))
	      softDropPuppiJetHFlav  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:hadronFlavour"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:partonFlavour"))
	      softDropPuppiJetPFlav .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:partonFlavour"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genMass"))
	      softDropPuppiJetGenm  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genMass"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPt"))
	      softDropPuppiJetGenpt  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPt"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genEta"))
	      softDropPuppiJetGeneta  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genEta"));
	    if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPhi"))
	      softDropPuppiJetGenphi  .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPhi"));	    
	  }
	
	  // sub-jets pruned 
	  if(puppiJetsBoosted[i]->hasSubjets("Pruned")){
	    pat::JetPtrCollection subjets = puppiJetsBoosted[i]->subjets("Pruned");
	    
	    if(subjets.size() > 0){
	      prunedPuppiSubJetpt_1  .push_back( subjets[0]->pt()); 
	      prunedPuppiSubJetm_1   .push_back( subjets[0]->mass()); 
	      prunedPuppiSubJetphi_1 .push_back( subjets[0]->phi()); 
	      prunedPuppiSubJeteta_1 .push_back( subjets[0]->eta());
	      prunedPuppiSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      prunedPuppiSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	      prunedPuppiSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	      
	      if(subjets[0]->hasUserFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"))
		prunedPuppiSubJetQGL_1 .push_back( subjets[0]->userFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"));
	      
	      if(isMC){
		prunedPuppiSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
		prunedPuppiSubJetPFlav_1 .push_back( subjets[0]->partonFlavour()); 
		if(subjets[0]->genJet()){
		  prunedPuppiSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		  prunedPuppiSubJetGenm_1  .push_back( subjets[0]->genJet()->mass());
		  prunedPuppiSubJetGeneta_1 .push_back( subjets[0]->genJet()->eta());
		  prunedPuppiSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi());
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(0),"SubCSV",prunedPuppiSubJetBtagSF_1,prunedPuppiSubJetBtagSFUp_1,prunedPuppiSubJetBtagSFDown_1);
	      }
	    }
	    
	    if(subjets.size() > 1){
	      prunedPuppiSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	      prunedPuppiSubJetm_2   .push_back( subjets.at(1)->mass()); 
	      prunedPuppiSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	      prunedPuppiSubJeteta_2 .push_back( subjets.at(1)->eta());
	      prunedPuppiSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      prunedPuppiSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	      prunedPuppiSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"))
		prunedPuppiSubJetQGL_2 .push_back( subjets.at(1)->userFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"));
	      
	      if(isMC){
		prunedPuppiSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
		prunedPuppiSubJetPFlav_2 .push_back( subjets.at(1)->partonFlavour()); 
		if(subjets.at(1)->genJet()){
		  prunedPuppiSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		  prunedPuppiSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		  prunedPuppiSubJetGenphi_2  .push_back( subjets.at(1)->genJet()->phi());
		  prunedPuppiSubJetGeneta_2  .push_back( subjets.at(1)->genJet()->eta());
		}	      
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(1),"SubCSV",prunedPuppiSubJetBtagSF_2,prunedPuppiSubJetBtagSFUp_2,prunedPuppiSubJetBtagSFDown_2);
	      }
	    }
	  }
	  
	  
	  // sub-jets soft drop 
	  if(puppiJetsBoosted[i]->hasSubjets("SoftDrop")){
	    pat::JetPtrCollection subjets = puppiJetsBoosted[i]->subjets("SoftDrop");
	    if(subjets.size() > 0){
	      softDropPuppiSubJetpt_1  .push_back( subjets[0]->pt()); 
	      softDropPuppiSubJetm_1   .push_back( subjets[0]->mass()); 
	      softDropPuppiSubJetphi_1 .push_back( subjets[0]->phi()); 
	      softDropPuppiSubJeteta_1 .push_back( subjets[0]->eta());
	      softDropPuppiSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

	      softDropPuppiSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	      softDropPuppiSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	    
	      if(subjets[0]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		softDropPuppiSubJetQGL_1 .push_back( subjets[0]->userFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"));
	      
	      if(isMC){
		softDropPuppiSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
		if(subjets[0]->genJet()){
		  softDropPuppiSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		  softDropPuppiSubJetGenm_1  .push_back( subjets[0]->genJet()->mass()); 
		  softDropPuppiSubJetGeneta_1  .push_back( subjets[0]->genJet()->eta()); 
		  softDropPuppiSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi()); 
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(0),"SubCSV",softDropPuppiSubJetBtagSF_1,softDropPuppiSubJetBtagSFUp_1,softDropPuppiSubJetBtagSFDown_1);
	      }
	    }
	  
	    if(subjets.size() > 1){
	      softDropPuppiSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	      softDropPuppiSubJetm_2   .push_back( subjets.at(1)->mass()); 
	      softDropPuppiSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	      softDropPuppiSubJeteta_2 .push_back( subjets.at(1)->eta());
	      softDropPuppiSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		softDropPuppiSubJetQGL_2 .push_back( subjets.at(1)->userFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"));

	      softDropPuppiSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	      softDropPuppiSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	      
	      if(isMC){
		softDropPuppiSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
		if(subjets.at(1)->genJet()){
		  softDropPuppiSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		  softDropPuppiSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		  softDropPuppiSubJetGeneta_2  .push_back( subjets.at(1)->genJet()->eta());
		  softDropPuppiSubJetGenphi_2  .push_back( subjets.at(1)->genJet()->phi());
		}
		if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(1),"SubCSV",softDropPuppiSubJetBtagSF_2,softDropPuppiSubJetBtagSFUp_2,softDropPuppiSubJetBtagSFDown_2);
	      }
	    }
	  }
	}
      }
    }

    // phtoon and electron ID info
    if(photonIDH.isValid() and addPhotonIDVariables and not isTriggerTree){

      photonPt.clear();
      photonEta.clear();
      photonPhi.clear();
      photonE.clear(); 
      photonSCEta.clear(); 
      photonSCPhi.clear();
      photonSCEnergy.clear();
      photonSCRawEnergy.clear();
      photonHOverE.clear();
      photonSigmaIetaIeta.clear();
      photonChargedIso.clear();					      
      photonNeutralIso.clear();
      photonEMIso.clear(); 
      photonElectronVeto.clear();

      for(auto photon_iter = photonIDH->begin(); photon_iter != photonIDH->end(); ++photon_iter){
	if(photon_iter->pt() < 35 or fabs(photon_iter->superCluster()->eta()) > 2.5) continue;
	photonPt.push_back(photon_iter->pt());
	photonEta.push_back(photon_iter->eta());	
	photonPhi.push_back(photon_iter->phi());
	photonE.push_back(photon_iter->energy());
	photonElectronVeto.push_back(photon_iter->passElectronVeto());
	photonSCEta.push_back(photon_iter->superCluster()->eta());
	photonSCPhi.push_back(photon_iter->superCluster()->phi());
	photonSCEnergy.push_back(photon_iter->superCluster()->energy());
	photonSCRawEnergy.push_back(photon_iter->superCluster()->rawEnergy());
	photonHOverE.push_back(photon_iter->hadTowOverEm());

	pat::PhotonRef phoref(photonIDH, photon_iter - photonIDH->begin());
	if(phoref.isAvailable() and phoref.isNonnull()){
	  if(photonsieieH.isValid())
	    photonSigmaIetaIeta.push_back((*photonsieieH)[phoref]);
	  if(photonCHisoH.isValid()){
	    photonChargedIso.push_back((*photonCHisoH)[phoref]);
	  }
	  if(photonPHisoH.isValid()){
	    photonEMIso.push_back((*photonPHisoH)[phoref]);
	  }
	  if(photonNHisoH.isValid()){
	    photonNeutralIso.push_back((*photonNHisoH)[phoref]);
	  }
	}	
      }      
    }

    if(electronIDH.isValid() and addElectronIDVariables and not isTriggerTree){

      electronPt.clear();
      electronEta.clear();
      electronPhi.clear();
      electronE.clear();
      electronSCEta.clear();
      electronSCPhi.clear();
      electronSCEnergy.clear();
      electronSCRawEnergy.clear();
      electronHOverE.clear();
      electronSigmaIetaIeta.clear();				    
      electronChargedIso.clear(); 
      electronNeutralIso.clear(); 
      electronEMIso.clear(); 
      electronGsfPt.clear(); 
      electronEOP.clear(); 
      electronDxy.clear();
      electronDz.clear();
      electronDphi.clear(); 
      electronDeta.clear(); 
      electronMissHit.clear(); 
      electronConversion.clear();

      for(auto electron_iter = electronIDH->begin(); electron_iter != electronIDH->end(); ++electron_iter){
	if(electron_iter->pt() < 35 or fabs(electron_iter->superCluster()->eta()) > 2.5) continue;
	electronPt.push_back(electron_iter->pt());
	electronEta.push_back(electron_iter->eta());	
	electronPhi.push_back(electron_iter->phi());
	electronE.push_back(electron_iter->energy());
	electronSCEta.push_back(electron_iter->superCluster()->eta());
	electronSCPhi.push_back(electron_iter->superCluster()->phi());
	electronSCEnergy.push_back(electron_iter->superCluster()->energy());
	electronSCRawEnergy.push_back(electron_iter->superCluster()->rawEnergy());
	electronHOverE.push_back(electron_iter->hadronicOverEm());
	electronSigmaIetaIeta.push_back(electron_iter->full5x5_sigmaIetaIeta());
	electronChargedIso.push_back(electron_iter->pfIsolationVariables().sumChargedHadronPt);
	electronNeutralIso.push_back(electron_iter->pfIsolationVariables().sumNeutralHadronEt);
	electronEMIso.push_back(electron_iter->pfIsolationVariables().sumPhotonEt);
	electronGsfPt.push_back(electron_iter->gsfTrack()->pt());
	electronDxy.push_back(electron_iter->gsfTrack()->dxy(verticesH->begin()->position()));
	electronDz.push_back(electron_iter->gsfTrack()->dz(verticesH->begin()->position()));
	electronEOP.push_back(fabs(1-electron_iter->eSuperClusterOverP())*1./electron_iter->ecalEnergy()); 
	electronDphi.push_back(electron_iter->deltaPhiSuperClusterTrackAtVtx());
	electronDeta.push_back(electron_iter->deltaEtaSuperClusterTrackAtVtx());
	electronMissHit.push_back(electron_iter->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
	electronConversion.push_back(electron_iter->passConversionVeto());
      }      
    }

    // Generator-level information
    wzid          = 0; wzmass        = 0.0; wzpt          = 0.0; wzeta         = 0.0; wzphi         = 0.0;
    l1id          = 0; l1pt          = 0.0; l1eta         = 0.0; l1phi         = 0.0;
    l2id          = 0; l2pt          = 0.0; l2eta         = 0.0; l2phi         = 0.0;
    wzid_h        = 0; wzmass_h      = 0.0; wzpt_h        = 0.0; wzeta_h       = 0.0; wzphi_h       = 0.0;
    topmass       = 0; toppt         = 0.0; topeta        = 0.0; topphi        = 0.0;
    atopmass      = 0; atoppt        = 0.0; atopeta       = 0.0; atopphi       = 0.0;
    q1id          = 0; q1pt          = 0.0; q1eta         = 0.0; q1phi         = 0.0;
    q2id          = 0; q2pt          = 0.0; q2eta         = 0.0; q2phi         = 0.0;
    parid         = 0; parpt         = 0.0; pareta        = 0.0; parphi        = 0.0;
    ancid         = 0; ancpt         = 0.0; anceta        = 0.0; ancphi        = 0.0; ancmass       = 0;
    wzmothid      = 0.0;
    isdirect      = 0;
    ismatch       = 0;

    dmmass   = 0.; dmphi   = 0.; dmeta   = 0.; dmpt   = 0.; dmid   = 0;
    dmX1mass = 0.; dmX1phi = 0.; dmX1eta = 0.; dmX1pt = 0.; dmX1id = 0;
    dmX2mass = 0.; dmX2phi = 0.; dmX2eta = 0.; dmX2pt = 0.; dmX2id = 0;

    if (isSignalSample and gensH.isValid() and not isTriggerTree) {

      TLorentzVector dm1vec; 
      TLorentzVector dm2vec; 
      bool foundfirst = false;
      
      // loop on gen particles looking for the DM particles --> then to the mediator
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) { 

	bool goodParticle = false;
	if (abs(gens_iter->pdgId()) >= 1000001 and abs(gens_iter->pdgId()) <= 1000039)
	  goodParticle = true;
	else if(abs(gens_iter->pdgId()) >= 2000001 and abs(gens_iter->pdgId()) <= 2000015)
	  goodParticle = true;
	else if(abs(gens_iter->pdgId()) == 9100012)
	  goodParticle = true;
	
	if(not goodParticle)
	  continue;
	
	if(!foundfirst) { // first DM particle
	  dm1vec.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	  dmX1id = gens_iter->pdgId();
	  foundfirst = true;

	  if(readDMFromGenParticle)
	    sampledmM = gens_iter->mass();
	}
	else{
	  dm2vec.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	  dmX2id = gens_iter->pdgId();
	  break;
	}
      }
      
      dmX1pt   = dm1vec.Pt();
      dmX1eta  = dm1vec.Eta();
      dmX1phi  = dm1vec.Phi();
      dmX1mass = dm1vec.M();

      dmX2pt   = dm2vec.Pt();
      dmX2eta  = dm2vec.Eta();
      dmX2phi  = dm2vec.Phi();
      dmX2mass = dm2vec.M();
      
      TLorentzVector medvec(dm1vec);
      medvec += dm2vec;
      dmpt  = medvec.Pt();
      dmeta = medvec.Eta();
      dmphi = medvec.Phi();
      dmmass = medvec.M();


      if(foundfirst == false){ //not found the DM particles and mediator --> look for Higgs invisible
	
	for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter){
	  if(gens_iter->pdgId() != 25) continue;
	  if(gens_iter->numberOfDaughters() <= 1) continue;

	  dmpt   = gens_iter->pt();
	  dmeta  = gens_iter->eta();
	  dmphi  = gens_iter->phi();
	  dmmass = gens_iter->mass();
	  dmid   = gens_iter->pdgId();

	  dmX1pt   = gens_iter->daughter(0)->pt();
	  dmX1eta  = gens_iter->daughter(0)->eta();
	  dmX1phi  = gens_iter->daughter(0)->phi();
	  dmX1mass = gens_iter->daughter(0)->mass();
	  dmX1id   = gens_iter->daughter(0)->pdgId();
	  dmX2pt   = gens_iter->daughter(1)->pt();
	  dmX2eta  = gens_iter->daughter(1)->eta();
	  dmX2phi  = gens_iter->daughter(1)->phi();
	  dmX2mass = gens_iter->daughter(1)->mass();
	  dmX2id   = gens_iter->daughter(1)->pdgId();
	}
      }
    }

    // dump inportant gen particles
    if(addGenParticles and gensH.isValid() and not isTriggerTree){
      
      // loop on genParticles (prunedGenParticles) trying to find W/Z decying leptonically or hadronically, top and anti-top quarks
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	if ( (gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
	     gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
	     abs(gens_iter->daughter(0)->pdgId()) > 10 && 
	     abs(gens_iter->daughter(0)->pdgId()) < 17)  { // decays into leptons, neutrinos 
	  
	  wzid   = gens_iter->pdgId();
	  wzmass = gens_iter->mass();
	  wzpt   = gens_iter->pt();
	  wzeta  = gens_iter->eta();
	  wzphi  = gens_iter->phi();
	  
	  l1id   = gens_iter->daughter(0)->pdgId();
	  l1pt   = gens_iter->daughter(0)->pt();
	  l1eta  = gens_iter->daughter(0)->eta();
	  l1phi  = gens_iter->daughter(0)->phi();
	
	  l2id   = gens_iter->daughter(1)->pdgId();
	  l2pt   = gens_iter->daughter(1)->pt();
	  l2eta  = gens_iter->daughter(1)->eta();
	  l2phi  = gens_iter->daughter(1)->phi();
	  wzmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
	}
      
	else if ( (gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
		  gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
		  ( (abs(gens_iter->daughter(0)->pdgId()) > 0 && abs(gens_iter->daughter(0)->pdgId()) <= 5) or
		    (abs(gens_iter->daughter(1)->pdgId()) > 0 && abs(gens_iter->daughter(1)->pdgId()) <= 5)))  { // decays into quarks
	  
	  wzid_h   = gens_iter->pdgId();
	  wzmass_h = gens_iter->mass();
	  wzpt_h   = gens_iter->pt();
	  wzeta_h  = gens_iter->eta();
	  wzphi_h  = gens_iter->phi();
	  
	  q1id   = gens_iter->daughter(0)->pdgId();
	  q1pt   = gens_iter->daughter(0)->pt();
	  q1eta  = gens_iter->daughter(0)->eta();
	  q1phi  = gens_iter->daughter(0)->phi();
	  
	  q2id   = gens_iter->daughter(1)->pdgId();
	  q2pt   = gens_iter->daughter(1)->pt();
	  q2eta  = gens_iter->daughter(1)->eta();
	  q2phi  = gens_iter->daughter(1)->phi();
	  wzmt_h   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
	}		  
	
	else if(gens_iter->pdgId() == 6 and gens_iter->numberOfDaughters() > 1 and 
		(((abs(gens_iter->daughter(0)->pdgId()) > 0 and  abs(gens_iter->daughter(0)->pdgId()) <= 5) and abs(gens_iter->daughter(1)->pdgId()) == 24) or
		 ((abs(gens_iter->daughter(1)->pdgId()) > 0 and  abs(gens_iter->daughter(1)->pdgId()) <= 5) and abs(gens_iter->daughter(0)->pdgId()) == 24))){
	  
	  topmass = gens_iter->mass();
	  toppt   = gens_iter->pt();
	  topeta  = gens_iter->eta();
	  topphi  = gens_iter->eta();
	}
	else if(gens_iter->pdgId() == -6 and gens_iter->numberOfDaughters() > 1 and 
		(((abs(gens_iter->daughter(0)->pdgId()) > 0 and  abs(gens_iter->daughter(0)->pdgId()) <= 5) and abs(gens_iter->daughter(1)->pdgId()) == 24) or
		 ((abs(gens_iter->daughter(1)->pdgId()) > 0 and  abs(gens_iter->daughter(1)->pdgId()) <= 5) and abs(gens_iter->daughter(0)->pdgId()) == 24))){
	  
	  atopmass = gens_iter->mass();
	  atoppt   = gens_iter->pt();
	  atopeta  = gens_iter->eta();
	  atopphi  = gens_iter->eta();	  
	}	
      }
    
    
      // if a Z/W is not found look for a pair of lepton .. this way with the pdgId is not guaranteed that you catch a Z/W boson and also recover DY production
      if (wzid == 0) {
	for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	  if (gens_iter->isPromptFinalState() || gens_iter->isPromptDecayed()) {
	    if (gens_iter->pdgId() >  10 && gens_iter->pdgId() <  17) {
	      l1id   = gens_iter->pdgId();
	      l1pt   = gens_iter->pt();
	      l1eta  = gens_iter->eta();
	      l1phi  = gens_iter->phi();
	    }
	    if (gens_iter->pdgId() < -10 && gens_iter->pdgId() > -17) {
	      l2id   = gens_iter->pdgId();
	      l2pt   = gens_iter->pt();
	      l2eta  = gens_iter->eta();
	      l2phi  = gens_iter->phi();
	    }
	  }
	}
	if (l1id > 0) {
	  TLorentzVector l1vec;
	  TLorentzVector l2vec;
	  l1vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, 0.);
	  l2vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, 0.);
	  TLorentzVector wzvec(l1vec);
	  wzvec += l2vec;
	  wzmass = wzvec.M();
	  wzpt   = wzvec.Pt();
	  wzeta  = wzvec.Eta();
	  wzphi  = wzvec.Phi();
	  wzmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
	  if (l1id+l2id == 0) wzid = 23;
	  else                wzid = 24;
	}
      }
      
      // no W or Z decay leptonically
      if (wzid == 0) {

	for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) { // loop on prunedGenParticles                                                          
	  if (gens_iter->pdgId() == 22 && // photons                                                                                                                           
	      gens_iter->status() == 1 && // final state                                                                                                                       
	      gens_iter->isPromptFinalState() &&
	      gens_iter->pt() > wzpt) {

	    wzid   = gens_iter->pdgId();
	    wzpt   = gens_iter->pt();
	    wzeta  = gens_iter->eta();
	    wzphi  = gens_iter->phi();
	      
	    findFirstNonPhotonMother(&(*gens_iter), ancid, ancpt, anceta, ancphi);
	    findMother(&(*gens_iter), parid, parpt, pareta, parphi);	    	    
	  }
	  
	  if(addPhotonPurity){
	    if (gens_iter->pdgId() == 22 && // photons
		gens_iter->status() == 1 && // final state
		gens_iter->isPromptFinalState() &&
		gens_iter->pt() >= wzpt) {
	      
	      findFirstNonPhotonMother(&(*gens_iter), ancid, ancpt, anceta, ancphi);
	      findMother(&(*gens_iter), parid, parpt, pareta, parphi);
	      if( abs(ancid) <= 5 || abs(ancid)==2212){ 
		double dR = computeDR(&(*gens_iter),tightphotonsPurity[hardestPhotonPurityIndex] );
		wzid   = gens_iter->pdgId();
		wzpt   = gens_iter->pt();
		wzeta  = gens_iter->eta();
		wzphi  = gens_iter->phi();
		wzmothid = gens_iter->mother(0)->pdgId();
		if(dR<0.3 && fabs((tightphotonsPurity[hardestPhotonPurityIndex]->pt()-gens_iter->pt())/tightphotonsPurity[hardestPhotonPurityIndex]->pt()) < 0.5){
		  ismatch=1;
		  TLorentzVector genpho;
		  genpho.SetPtEtaPhiM(wzpt, wzeta,wzphi,0);
		  TLorentzVector genpart;
		  genpart.SetPtEtaPhiM(ancpt, anceta,ancphi,ancmass);
		  double dRFrag = genpho.DeltaR(genpart);
		  if(dRFrag>0.4)isdirect=1;
		}
	      }
	    }          
	  }
	}
      }
	}
    tree->Fill();    
}    

void MonoJetTreeMaker::beginJob() {

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree"       , "tree");

  // Run, Lumi, Event info
  tree->Branch("event"                , &event                , "event/I");
  tree->Branch("run"                  , &run                  , "run/I");
  tree->Branch("lumi"                 , &lumi                 , "lumi/I");
  // Event weights
  // Pileup info
  if(not isTriggerTree){
    tree->Branch("puwgt"                , &puwgt                , "puwgt/D");
    tree->Branch("puobs"                , &puobs                , "puobs/I");
    tree->Branch("xsec"                 , &xsec                 , "xsec/D");
    tree->Branch("wgt"                  , &wgt                  , "wgt/D");
  }

  tree->Branch("putrue"               , &putrue               , "putrue/I");
  tree->Branch("nvtx"                 , &nvtx                 , "nvtx/I");
  
  // Triggers
  tree->Branch("hltmet90"             , &hltmet90             , "hltmet90/B");
  tree->Branch("hltmet100"            , &hltmet100            , "hltmet100/B");
  tree->Branch("hltmet110"            , &hltmet110            , "hltmet110/B");
  tree->Branch("hltmet120"            , &hltmet120            , "hltmet120/B");
  tree->Branch("hltmetwithmu90"       , &hltmetwithmu90       , "hltmetwithmu90/B");
  tree->Branch("hltmetwithmu100"      , &hltmetwithmu100      , "hltmetwithmu100/B");
  tree->Branch("hltmetwithmu110"      , &hltmetwithmu110      , "hltmetwithmu110/B");
  tree->Branch("hltmetwithmu120"      , &hltmetwithmu120      , "hltmetwithmu120/B");
  tree->Branch("hltmetwithmu170"      , &hltmetwithmu170      , "hltmetwithmu170/B");
  tree->Branch("hltmetwithmu300"      , &hltmetwithmu300      , "hltmetwithmu300/B");
  tree->Branch("hltjetmet"            , &hltjetmet            , "hltjetmet/B");
  tree->Branch("hltphoton90"          , &hltphoton90          , "hltphoton90/B");
  tree->Branch("hltphoton120"         , &hltphoton120         , "hltphoton120/B");
  tree->Branch("hltphoton120vbf"      , &hltphoton120vbf      , "hltphoton120vbf/B");
  tree->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/B");
  tree->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/B");

  tree->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/B");
  tree->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/B");
  tree->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/B");
  tree->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/B");
  tree->Branch("hltsingleel27"        , &hltsingleel27        , "hltsingleel27/B");
  tree->Branch("hltelnoiso"           , &hltelnoiso           , "hltelnoiso/B");

  tree->Branch("hltPFHT400"           , &hltPFHT400           , "hltPFHT400/B");
  tree->Branch("hltPFHT475"           , &hltPFHT475           , "hltPFHT475/B");
  tree->Branch("hltPFHT600"           , &hltPFHT600           , "hltPFHT600/B");
  tree->Branch("hltPFHT650"           , &hltPFHT650           , "hltPFHT650/B");
  tree->Branch("hltPFHT800"           , &hltPFHT800           , "hltPFHT800/B");
  tree->Branch("hltPFHT900"           , &hltPFHT900           , "hltPFHT900/B");
  tree->Branch("hltEcalHT800"         , &hltEcalHT800         , "hltEcalHT800/B");
  tree->Branch("hltphoton90PFHT"      , &hltphoton90PFHT      , "hltphoton90PFHT/B");

  tree->Branch("pswgt_ph120"          , &pswgt_ph120          , "pswgt_ph120/D");
  tree->Branch("pswgt_ph90"           , &pswgt_ph90           , "pswgt_ph90/D");
    
  if(isTriggerTree){
    tree->Branch("pswgt_ht400"          , &pswgt_ht400          , "pswgt_ht400/D");
    tree->Branch("pswgt_ht475"          , &pswgt_ht475          , "pswgt_ht475/D");
    tree->Branch("pswgt_ht600"          , &pswgt_ht600          , "pswgt_ht600/D");
    tree->Branch("pswgt_ht650"          , &pswgt_ht650          , "pswgt_ht650/D");
    tree->Branch("pswgt_ht800"          , &pswgt_ht800          , "pswgt_ht800/D");
    tree->Branch("pswgt_ht900"          , &pswgt_ht900          , "pswgt_ht900/D");
  }

  // MET filters
  tree->Branch("flagcsctight"         , &flagcsctight         , "flagcsctight/B");
  tree->Branch("flaghbhenoise"        , &flaghbhenoise        , "flaghbhenoise/B");
  tree->Branch("flaghbheiso"          , &flaghbheiso          , "flaghbheiso/B");
  tree->Branch("flageebadsc"          , &flageebadsc          , "flageebadsc/B");
  tree->Branch("flagecaltp"           , &flagecaltp           , "flagecaltp/B");
  tree->Branch("flaggoodvertices"     , &flaggoodvertices     , "flaggoodvertices/B");
  tree->Branch("flagglobaltighthalo"  , &flagglobaltighthalo  , "flagglobaltighthalo/B");
  tree->Branch("flagbadchpf"          , &flagbadchpf          , "flagbadchpf/B");
  tree->Branch("flagbadpfmu"          , &flagbadpfmu          , "flagbadpfmu/B");

  if(isTriggerTree and addTriggerObjects){

    tree->Branch("trig_obj_n"           , &trig_obj_n           , "trig_obj_n/I"); //ND
    tree->Branch("trig_obj_pt"          , "std::vector<double>" , &trig_obj_pt);   //ND
    tree->Branch("trig_obj_eta"         , "std::vector<double>" , &trig_obj_eta);  //ND
    tree->Branch("trig_obj_phi"         , "std::vector<double>" , &trig_obj_phi);  //ND
    tree->Branch("trig_obj_col"         , "std::vector<std::string>" , &trig_obj_col); //, buffersize); //ND 

    tree->Branch("trig_L1A_check"       , &trig_L1A_check            , "trig_L1A_check/I"); //ND
    tree->Branch("trig_L1A_n"           , &trig_L1A_n                , "trig_L1A_n/I"); //ND 
    tree->Branch("trig_L1A_list"        , "std::vector<std::string>" , &trig_L1A_list); //ND
 
    tree->Branch("trig_L1EG_pt"         , "std::vector<double>" , &trig_L1EG_pt); //ND
    tree->Branch("trig_L1EG_eta"        , "std::vector<double>" , &trig_L1EG_eta); //ND
    tree->Branch("trig_L1EG_phi"        , "std::vector<double>" , &trig_L1EG_phi); //ND 
    //
    tree->Branch("trig_L1Jet_pt"        , "std::vector<double>" , &trig_L1Jet_pt); //ND
    tree->Branch("trig_L1Jet_eta"       , "std::vector<double>" , &trig_L1Jet_eta); //ND
    tree->Branch("trig_L1Jet_phi"       , "std::vector<double>" , &trig_L1Jet_phi); //ND
    //
    tree->Branch("trig_L1Mu_pt"         , "std::vector<double>" , &trig_L1Mu_pt); //ND
    tree->Branch("trig_L1Mu_eta"        , "std::vector<double>" , &trig_L1Mu_eta); //ND
    tree->Branch("trig_L1Mu_phi"        , "std::vector<double>" , &trig_L1Mu_phi); //ND
    //
    tree->Branch("trig_L1ETM_pt"        , &trig_L1ETM_pt        , "trig_L1ETM_pt/D"); //ND
    tree->Branch("trig_L1ETM_phi"       , &trig_L1ETM_phi       , "trig_L1ETM_phi/D"); //ND
    //
    tree->Branch("trig_L1ETT_pt"        , &trig_L1ETT_pt        , "trig_L1ETT_pt/D"); //ND
    tree->Branch("trig_L1ETT_phi"       , &trig_L1ETT_phi       , "trig_L1ETT_phi/D"); //ND
    //
    tree->Branch("trig_L1HTM_pt"        , &trig_L1HTM_pt        , "trig_L1HTM_pt/D"); //ND
    tree->Branch("trig_L1HTM_phi"       , &trig_L1HTM_phi       , "trig_L1HTM_phi/D"); //ND
    //
    tree->Branch("trig_L1HTT_pt"        , &trig_L1HTT_pt        , "trig_L1HTT_pt/D"); //ND
    tree->Branch("trig_L1HTT_phi"       , &trig_L1HTT_phi       , "trig_L1HTT_phi/D"); //ND


  }

  // Object counts
  tree->Branch("nmuons"               , &nmuons               , "nmuons/I");
  tree->Branch("nelectrons"           , &nelectrons           , "nelectrons/I");
  tree->Branch("nlooseelectrons"      , &nlooseelectrons      , "nlooseelectrons/I");
  tree->Branch("ntightmuons"          , &ntightmuons          , "ntightmuons/I");
  tree->Branch("nhighptmuons"         , &nhighptmuons         , "nhighptmuons/I");
  tree->Branch("ntightelectrons"      , &ntightelectrons      , "ntightelectrons/I");
  tree->Branch("nheepelectrons"       , &nheepelectrons       , "nheepelectrons/I");
  tree->Branch("ntaus"                , &ntaus                , "ntaus/I");
  tree->Branch("ntausraw"             , &ntausraw             , "ntausraw/I");
  tree->Branch("ntausold"             , &ntausold             , "ntausold/I");
  tree->Branch("ntausrawold"          , &ntausrawold          , "ntausrawold/I");
  tree->Branch("nphotons"             , &nphotons             , "nphotons/I");
  tree->Branch("njets"                , &njets                , "njets/I");
  tree->Branch("njetsinc"             , &njetsinc             , "njetsinc/I");
  tree->Branch("nbjets"               , &nbjets               , "nbjets/I");
  tree->Branch("nbjetslowpt"          , &nbjetslowpt          , "nbjetslowpt/I");
  if(not isTriggerTree){
    tree->Branch("nbjetsMVA"            , &nbjetsMVA            , "nbjetsMVA/I");
    tree->Branch("nbjetsMVAlowpt"       , &nbjetsMVAlowpt       , "nbjetsMVAlowpt/I");
  }

  if(addPuppiJets and not isTriggerTree){
    tree->Branch("npuppijets"                , &npuppijets                , "npuppijets/I");
    tree->Branch("npuppijetsinc"             , &npuppijetsinc             , "npuppijetsinc/I");
    tree->Branch("npuppibjets"               , &npuppibjets               , "npuppibjets/I");
    tree->Branch("npuppibjetslowpt"          , &npuppibjetslowpt          , "npuppibjetslowpt/I");
    tree->Branch("npuppibjetsMVA"            , &npuppibjetsMVA            , "npuppibjetsMVA/I");
    tree->Branch("npuppibjetsMVAlowpt"       , &npuppibjetsMVAlowpt       , "npuppibjetsMVAlowpt/I");
  }

  // MET
  tree->Branch("pfmet"                , &pfmet                , "pfmet/D");
  tree->Branch("pfmetphi"             , &pfmetphi             , "pfmetphi/D");
  tree->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/D");
  tree->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/D");
  tree->Branch("calomet"              , &calomet              , "calomet/D");   //ND
  tree->Branch("calometphi"           , &calometphi           , "calometphi/D");//ND
  if(not isTriggerTree){    
    tree->Branch("mumet"                , &mumet                , "mumet/D");
    tree->Branch("mumetphi"             , &mumetphi             , "mumetphi/D");
  }
  tree->Branch("t1mumet"              , &t1mumet              , "t1mumet/D");
  tree->Branch("t1mumetphi"           , &t1mumetphi           , "t1mumetphi/D");
  if(not isTriggerTree){
    tree->Branch("elmet"                , &elmet                , "elmet/D");
    tree->Branch("elmetphi"             , &elmetphi             , "elmetphi/D");
  }
  tree->Branch("t1elmet"              , &t1elmet              , "t1elmet/D");
  tree->Branch("t1elmetphi"           , &t1elmetphi           , "t1elmetphi/D");
  if(not isTriggerTree){
    tree->Branch("phmet"                , &phmet                , "phmet/D");
    tree->Branch("phmetphi"             , &phmetphi             , "phmetphi/D");
  }
  tree->Branch("t1phmet"              , &t1phmet              , "t1phmet/D");
  tree->Branch("t1phmetphi"           , &t1phmetphi           , "t1phmetphi/D");
  if(not isTriggerTree){
    tree->Branch("taumet"                , &taumet                , "taumet/D");
    tree->Branch("taumetphi"             , &taumetphi             , "taumetphi/D");
    tree->Branch("t1taumet"              , &t1taumet              , "t1taumet/D");
    tree->Branch("t1taumetphi"           , &t1taumetphi           , "t1taumetphi/D");
    tree->Branch("genmet",    &genmet,   "genmet/D");
    tree->Branch("genmetphi", &genmetphi,"genmetphi/D");
  }

  if(addMETBreakDown and not isTriggerTree){
    
    tree->Branch("pfmethadronHF",&pfmethadronHF,"pfmethadronHF/D");
    tree->Branch("pfmethadronHFphi",&pfmethadronHFphi,"pfmethadronHFphi/D");
    tree->Branch("pfmetegammaHF",&pfmetegammaHF,"pfmetegammaHF/D");
    tree->Branch("pfmetegammaHFphi",&pfmetegammaHFphi,"pfmetegammaHFphi/D");
    tree->Branch("pfmetchargedhadron",&pfmetchargedhadron,"pfmetchargedhadron/D");
    tree->Branch("pfmetchargedhadronphi",&pfmetchargedhadronphi,"pfmetchargedhadronphi/D");
    tree->Branch("pfmetneutralhadron",&pfmetneutralhadron,"pfmetneutralhadron/D");
    tree->Branch("pfmetneutralhadronphi",&pfmetneutralhadronphi,"pfmetneutralhadronphi/D");
    tree->Branch("pfmetelectrons",&pfmetelectrons,"pfmetelectrons/D");
    tree->Branch("pfmetelectronsphi",&pfmetelectronsphi,"pfmetelectronsphi/D");
    tree->Branch("pfmetmuons",&pfmetmuons,"pfmetmuons/D");
    tree->Branch("pfmetmuonsphi",&pfmetmuonsphi,"pfmetmuonsphi/D");
    tree->Branch("pfmetphotons",&pfmetphotons,"pfmetphotons/D");
    tree->Branch("pfmetphotonsphi",&pfmetphotonsphi,"pfmetphotonsphi/D");
    tree->Branch("pfmetunclustered",&pfmetunclustered,"pfmetunclustered/D");
    tree->Branch("pfmetunclusteredphi",&pfmetunclusteredphi,"pfmetunclusteredphi/D");
  }

  if(addMVAMet and not isTriggerTree){
    tree->Branch("mvamet"              , &mvamet              , "mvamet/D");
    tree->Branch("mvametphi"           , &mvametphi           , "mvametphi/D");
  }

  if(addMETSystematics and not isTriggerTree){
    tree->Branch("t1pfmetMuEnUp"       , &t1pfmetMuEnUp       , "t1pfmetMuEnUp/D");
    tree->Branch("t1pfmetMuEnDown"     , &t1pfmetMuEnDown     , "t1pfmetMuEnDown/D");
    tree->Branch("t1pfmetElEnUp"       , &t1pfmetElEnUp       , "t1pfmetElEnUp/D");
    tree->Branch("t1pfmetElEnDown"     , &t1pfmetElEnDown     , "t1pfmetElEnDown/D");
    tree->Branch("t1pfmetPhoEnUp"      , &t1pfmetPhoEnUp      , "t1pfmetPhoEnUp/D");
    tree->Branch("t1pfmetPhoEnDown"    , &t1pfmetPhoEnDown    , "t1pfmetPhoEnDown/D");
    tree->Branch("t1pfmetTauEnUp"      , &t1pfmetTauEnUp      , "t1pfmetTauEnUp/D");
    tree->Branch("t1pfmetTauEnDown"    , &t1pfmetTauEnDown    , "t1pfmetTauEnDown/D");
    tree->Branch("t1pfmetJetEnUp"      , &t1pfmetJetEnUp      , "t1pfmetJetEnUp/D");
    tree->Branch("t1pfmetJetEnDown"    , &t1pfmetJetEnDown    , "t1pfmetJetEnDown/D");
    tree->Branch("t1pfmetJetResUp"     , &t1pfmetJetResUp     , "t1pfmetJetResUp/D");
    tree->Branch("t1pfmetJetResDown"   , &t1pfmetJetResDown   , "t1pfmetJetResDown/D");
    tree->Branch("t1pfmetUncEnUp"      , &t1pfmetUncEnUp      , "t1pfmetUncEnUp/D");
    tree->Branch("t1pfmetUncEnDown"    , &t1pfmetUncEnDown    , "t1pfmetUncEnDown/D");
    tree->Branch("t1pfmetJetSmear"        , &t1pfmetJetSmear        , "t1pfmetJetSmear/D");
    tree->Branch("t1pfmetXY"           , &t1pfmetXY           , "t1pfmetXY/D");

    tree->Branch("t1pfmetMuEnUpPhi"       , &t1pfmetMuEnUpPhi       , "t1pfmetMuEnUpPhi/D");
    tree->Branch("t1pfmetMuEnDownPhi"     , &t1pfmetMuEnDownPhi     , "t1pfmetMuEnDownPhi/D");
    tree->Branch("t1pfmetElEnUpPhi"       , &t1pfmetElEnUpPhi       , "t1pfmetElEnUpPhi/D");
    tree->Branch("t1pfmetElEnDownPhi"     , &t1pfmetElEnDownPhi     , "t1pfmetElEnDownPhi/D");
    tree->Branch("t1pfmetPhoEnUpPhi"      , &t1pfmetPhoEnUpPhi      , "t1pfmetPhoEnUpPhi/D");
    tree->Branch("t1pfmetPhoEnDownPhi"    , &t1pfmetPhoEnDownPhi    , "t1pfmetPhoEnDownPhi/D");
    tree->Branch("t1pfmetTauEnUpPhi"      , &t1pfmetTauEnUpPhi      , "t1pfmetTauEnUpPhi/D");
    tree->Branch("t1pfmetTauEnDownPhi"    , &t1pfmetTauEnDownPhi    , "t1pfmetTauEnDownPhi/D");
    tree->Branch("t1pfmetJetEnUpPhi"      , &t1pfmetJetEnUpPhi      , "t1pfmetJetEnUpPhi/D");
    tree->Branch("t1pfmetJetEnDownPhi"    , &t1pfmetJetEnDownPhi    , "t1pfmetJetEnDownPhi/D");
    tree->Branch("t1pfmetJetResUpPhi"     , &t1pfmetJetResUpPhi     , "t1pfmetJetResUpPhi/D");
    tree->Branch("t1pfmetJetResDownPhi"   , &t1pfmetJetResDownPhi   , "t1pfmetJetResDownPhi/D");
    tree->Branch("t1pfmetUncEnUpPhi"      , &t1pfmetUncEnUpPhi      , "t1pfmetUncEnUpPhi/D");
    tree->Branch("t1pfmetUncEnDownPhi"    , &t1pfmetUncEnDownPhi    , "t1pfmetUncEnDownPhi/D");
    tree->Branch("t1pfmetJetSmearPhi"        , &t1pfmetJetSmearPhi        , "t1pfmetJetSmearPhi/D");
    tree->Branch("t1pfmetXYPhi"           , &t1pfmetXYPhi           , "t1pfmetXYPhi/D");

  }

  if(addPuppiMET and not isTriggerTree){
    tree->Branch("puppipfmet"                , &puppipfmet                , "puppipfmet/D");
    tree->Branch("puppipfmetphi"             , &puppipfmetphi             , "puppipfmetphi/D");
    tree->Branch("puppit1pfmet"              , &puppit1pfmet              , "puppit1pfmet/D");
    tree->Branch("puppit1pfmetphi"           , &puppit1pfmetphi           , "puppit1pfmetphi/D");
    tree->Branch("puppimumet"                , &puppimumet                , "puppimumet/D");
    tree->Branch("puppimumetphi"             , &puppimumetphi             , "puppimumetphi/D");
    tree->Branch("puppit1mumet"              , &puppit1mumet              , "puppit1mumet/D");
    tree->Branch("puppit1mumetphi"           , &puppit1mumetphi           , "puppit1mumetphi/D");
    tree->Branch("puppielmet"                , &puppielmet                , "puppielmet/D");
    tree->Branch("puppielmetphi"             , &puppielmetphi             , "elmetphi/D");
    tree->Branch("puppit1elmet"              , &puppit1elmet              , "puppit1elmet/D");
    tree->Branch("puppit1elmetphi"           , &puppit1elmetphi           , "puppit1elmetphi/D");
    tree->Branch("puppiphmet"                , &puppiphmet                , "puppiphmet/D");
    tree->Branch("puppiphmetphi"             , &puppiphmetphi             , "puppiphmetphi/D");
    tree->Branch("puppit1phmet"              , &puppit1phmet              , "puppit1phmet/D");
    tree->Branch("puppit1phmetphi"           , &puppit1phmetphi           , "puppit1phmetphi/D");

    if(addMETSystematics and not isTriggerTree){
      tree->Branch("puppit1pfmetMuEnUp"       , &puppit1pfmetMuEnUp       , "puppit1pfmetMuEnUp/D");
      tree->Branch("puppit1pfmetMuEnDown"     , &puppit1pfmetMuEnDown     , "puppit1pfmetMuEnDown/D");
      tree->Branch("puppit1pfmetElEnUp"       , &puppit1pfmetElEnUp       , "puppit1pfmetElEnUp/D");
      tree->Branch("puppit1pfmetElEnDown"     , &puppit1pfmetElEnDown     , "puppit1pfmetElEnDown/D");
      tree->Branch("puppit1pfmetPhoEnUp"      , &puppit1pfmetPhoEnUp      , "puppit1pfmetPhoEnUp/D");
      tree->Branch("puppit1pfmetPhoEnDown"    , &puppit1pfmetPhoEnDown    , "puppit1pfmetPhoEnDown/D");
      tree->Branch("puppit1pfmetTauEnUp"      , &puppit1pfmetTauEnUp      , "puppit1pfmetTauEnUp/D");
      tree->Branch("puppit1pfmetTauEnDown"    , &puppit1pfmetTauEnDown    , "puppit1pfmetTauEnDown/D");
      tree->Branch("puppit1pfmetJetEnUp"      , &puppit1pfmetJetEnUp      , "puppit1pfmetJetEnUp/D");
      tree->Branch("puppit1pfmetJetEnDown"    , &puppit1pfmetJetEnDown    , "puppit1pfmetJetEnDown/D");
      tree->Branch("puppit1pfmetJetResUp"     , &puppit1pfmetJetResUp     , "puppit1pfmetJetResUp/D");
      tree->Branch("puppit1pfmetJetResDown"   , &puppit1pfmetJetResDown   , "puppit1pfmetJetResDown/D");
      tree->Branch("puppit1pfmetUncEnUp"      , &puppit1pfmetUncEnUp      , "puppit1pfmetUncEnUp/D");
      tree->Branch("puppit1pfmetUncEnDown"    , &puppit1pfmetUncEnDown    , "puppit1pfmetUncEnDown/D");

      tree->Branch("puppit1pfmetMuEnUpPhi"       , &puppit1pfmetMuEnUpPhi       , "puppit1pfmetMuEnUpPhi/D");
      tree->Branch("puppit1pfmetMuEnDownPhi"     , &puppit1pfmetMuEnDownPhi     , "puppit1pfmetMuEnDownPhi/D");
      tree->Branch("puppit1pfmetElEnUpPhi"       , &puppit1pfmetElEnUpPhi       , "puppit1pfmetElEnUpPhi/D");
      tree->Branch("puppit1pfmetElEnDownPhi"     , &puppit1pfmetElEnDownPhi     , "puppit1pfmetElEnDownPhi/D");
      tree->Branch("puppit1pfmetPhoEnUpPhi"      , &puppit1pfmetPhoEnUpPhi      , "puppit1pfmetPhoEnUpPhi/D");
      tree->Branch("puppit1pfmetPhoEnDownPhi"    , &puppit1pfmetPhoEnDownPhi    , "puppit1pfmetPhoEnDownPhi/D");
      tree->Branch("puppit1pfmetTauEnUpPhi"      , &puppit1pfmetTauEnUpPhi      , "puppit1pfmetTauEnUpPhi/D");
      tree->Branch("puppit1pfmetTauEnDownPhi"    , &puppit1pfmetTauEnDownPhi    , "puppit1pfmetTauEnDownPhi/D");
      tree->Branch("puppit1pfmetJetEnUpPhi"      , &puppit1pfmetJetEnUpPhi      , "puppit1pfmetJetEnUpPhi/D");
      tree->Branch("puppit1pfmetJetEnDownPhi"    , &puppit1pfmetJetEnDownPhi    , "puppit1pfmetJetEnDownPhi/D");
      tree->Branch("puppit1pfmetJetResUpPhi"     , &puppit1pfmetJetResUpPhi     , "puppit1pfmetJetResUpPhi/D");
      tree->Branch("puppit1pfmetJetResDownPhi"   , &puppit1pfmetJetResDownPhi   , "puppit1pfmetJetResDownPhi/D");
      tree->Branch("puppit1pfmetUncEnUpPhi"      , &puppit1pfmetUncEnUpPhi      , "puppit1pfmetUncEnUpPhi/D");
      tree->Branch("puppit1pfmetUncEnDownPhi"    , &puppit1pfmetUncEnDownPhi    , "puppit1pfmetUncEnDownPhi/D");


    }

  }
  
  // Jet info
  tree->Branch("combinejetpt",      "std::vector<double>", &combinejetpt);
  tree->Branch("combinejeteta",     "std::vector<double>", &combinejeteta);
  tree->Branch("combinejetphi",     "std::vector<double>", &combinejetphi);
  tree->Branch("combinejetm",       "std::vector<double>", &combinejetm);
  if(not isTriggerTree){
    tree->Branch("combinejetbtag",    "std::vector<double>", &combinejetbtag);
    tree->Branch("combinejetbtagMVA", "std::vector<double>", &combinejetbtagMVA);
  }
  tree->Branch("combinejetCHfrac",  "std::vector<double>", &combinejetCHfrac);
  tree->Branch("combinejetNHfrac",  "std::vector<double>", &combinejetNHfrac);
  if(not isTriggerTree){
    tree->Branch("combinejetEMfrac",  "std::vector<double>", &combinejetEMfrac);
    tree->Branch("combinejetCEMfrac", "std::vector<double>", &combinejetCEMfrac);
    tree->Branch("combinejetmetdphi", "std::vector<double>", &combinejetmetdphi);
    tree->Branch("combinejetHFlav",   "std::vector<double>", &combinejetHFlav);
    tree->Branch("combinejetPFlav",   "std::vector<double>", &combinejetPFlav);
    tree->Branch("combinejetQGL",     "std::vector<double>", &combinejetQGL);
    tree->Branch("combinejetPUIF",    "std::vector<double>", &combinejetPUID);
    tree->Branch("combinejetPassPUIF",    "std::vector<double>", &combinejetPassPUID);
    tree->Branch("combinejetGenpt",   "std::vector<double>", &combinejetGenpt);
    tree->Branch("combinejetGeneta",  "std::vector<double>", &combinejetGeneta);
    tree->Branch("combinejetGenphi",  "std::vector<double>", &combinejetGenphi);
    tree->Branch("combinejetGenm",    "std::vector<double>", &combinejetGenm);
    tree->Branch("combinejetBtagSF",  "std::vector<double>", &combinejetBtagSF);
    tree->Branch("combinejetBtagSFUp", "std::vector<double>", &combinejetBtagSFUp);
    tree->Branch("combinejetBtagSFDown",    "std::vector<double>", &combinejetBtagSFDown);
    tree->Branch("combinejetBtagMVASF",     "std::vector<double>", &combinejetBtagMVASF);
    tree->Branch("combinejetBtagMVASFUp",   "std::vector<double>", &combinejetBtagMVASFUp);
    tree->Branch("combinejetBtagMVASFDown", "std::vector<double>", &combinejetBtagMVASFDown);    
    tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");
  

    tree->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/D");
    tree->Branch("incjetmumetdphimin"   , &incjetmumetdphimin   , "incjetmumetdphimin/D");
    tree->Branch("incjetelmetdphimin"   , &incjetelmetdphimin   , "incjetelmetdphimin/D");
    tree->Branch("incjetphmetdphimin"   , &incjetphmetdphimin   , "incjetphmetdphimin/D");
    
    tree->Branch("alljetmetdphimin"     , &alljetmetdphimin     , "alljetmetdphimin/D");
    tree->Branch("alljetmumetdphimin"   , &alljetmumetdphimin   , "alljetmumetdphimin/D");
    tree->Branch("alljetelmetdphimin"   , &alljetelmetdphimin   , "alljetelmetdphimin/D");
    tree->Branch("alljetphmetdphimin"   , &alljetphmetdphimin   , "alljetphmetdphimin/D");
  }

  tree->Branch("incjetmetdphimin4"    , &incjetmetdphimin4    , "incjetmetdphimin4/D");
  tree->Branch("incjetmumetdphimin4"  , &incjetmumetdphimin4  , "incjetmumetdphimin4/D");
  tree->Branch("incjetelmetdphimin4"  , &incjetelmetdphimin4  , "incjetelmetdphimin4/D");
  tree->Branch("incjetphmetdphimin4"  , &incjetphmetdphimin4  , "incjetphmetdphimin4/D");

  if(not isTriggerTree){
    tree->Branch("alljetmetdphimin4"    , &alljetmetdphimin4    , "alljetmetdphimin4/D");
    tree->Branch("alljetmumetdphimin4"  , &alljetmumetdphimin4  , "alljetmumetdphimin4/D");
    tree->Branch("alljetelmetdphimin4"  , &alljetelmetdphimin4  , "alljetelmetdphimin4/D");
    tree->Branch("alljetphmetdphimin4"  , &alljetphmetdphimin4  , "alljetphmetdphimin4/D");    
    tree->Branch("ht"                   , &ht                   , "ht/D");
    tree->Branch("htinc"                , &htinc                , "htinc/D");
    tree->Branch("ht30"                 , &ht30                 , "ht30/D");
  }

  if(addPuppiJets and not isTriggerTree){

    tree->Branch("combinePuppijetpt",  "std::vector<double>", &combinePuppijetpt);
    tree->Branch("combinePuppijeteta", "std::vector<double>", &combinePuppijeteta);
    tree->Branch("combinePuppijetphi", "std::vector<double>", &combinePuppijetphi);
    tree->Branch("combinePuppijetm",   "std::vector<double>", &combinePuppijetm);
    tree->Branch("combinePuppijetbtag", "std::vector<double>", &combinePuppijetbtag);
    tree->Branch("combinePuppijetbtagMVA", "std::vector<double>", &combinePuppijetbtagMVA);
    tree->Branch("combinePuppijetCHfrac", "std::vector<double>", &combinePuppijetCHfrac);
    tree->Branch("combinePuppijetNHfrac", "std::vector<double>", &combinePuppijetNHfrac);
    tree->Branch("combinePuppijetEMfrac", "std::vector<double>", &combinePuppijetEMfrac);
    tree->Branch("combinePuppijetCEMfrac", "std::vector<double>", &combinePuppijetCEMfrac);
    tree->Branch("combinePuppijetmetdphi", "std::vector<double>", &combinePuppijetmetdphi);
    tree->Branch("combinePuppijetHFlav", "std::vector<double>", &combinePuppijetHFlav);
    tree->Branch("combinePuppijetPFlav", "std::vector<double>", &combinePuppijetPFlav);
    tree->Branch("combinePuppijetQGL",   "std::vector<double>", &combinePuppijetQGL);
    tree->Branch("combinePuppijetGenpt", "std::vector<double>", &combinePuppijetGenpt);
    tree->Branch("combinePuppijetGeneta", "std::vector<double>", &combinePuppijetGeneta);
    tree->Branch("combinePuppijetGenphi", "std::vector<double>", &combinePuppijetGenphi);
    tree->Branch("combinePuppijetGenm",   "std::vector<double>", &combinePuppijetGenm);
    tree->Branch("combinePuppijetBtagSF", "std::vector<double>", &combinePuppijetBtagSF);
    tree->Branch("combinePuppijetBtagSFUp", "std::vector<double>", &combinePuppijetBtagSFUp);
    tree->Branch("combinePuppijetBtagSFDown", "std::vector<double>", &combinePuppijetBtagSFDown);
    tree->Branch("combinePuppijetBtagMVASF", "std::vector<double>", &combinePuppijetBtagMVASF);
    tree->Branch("combinePuppijetBtagMVASFUp", "std::vector<double>", &combinePuppijetBtagMVASFUp);
    tree->Branch("combinePuppijetBtagMVASFDown", "std::vector<double>", &combinePuppijetBtagMVASFDown);

    tree->Branch("PuppijetPuppijetdphi"      , &PuppijetPuppijetdphi      , "PuppijetPuppijetdphi/D");

    tree->Branch("Puppijetmetdphimin"        , &Puppijetmetdphimin        , "Puppijetmetdphimin/D");
    tree->Branch("Puppijetmumetdphimin"        , &Puppijetmumetdphimin        , "Puppijetmumetdphimin/D");
    tree->Branch("Puppijetelmetdphimin"      , &Puppijetelmetdphimin      , "Puppijetelmetdphimin/D");
    tree->Branch("Puppijetphmetdphimin"      , &Puppijetphmetdphimin      , "Puppijetphmetdphimin/D");

    tree->Branch("incPuppijetmetdphimin"     , &incPuppijetmetdphimin     , "incPuppijetmetdphimin/D");
    tree->Branch("incPuppijetmumetdphimin"     , &incPuppijetmumetdphimin     , "incPuppijetmumetdphimin/D");
    tree->Branch("incPuppijetelmetdphimin"   , &incPuppijetelmetdphimin   , "incPuppijetelmetdphimin/D");
    tree->Branch("incPuppijetphmetdphimin"   , &incPuppijetphmetdphimin   , "incPuppijetphmetdphimin/D");

    tree->Branch("Puppijetmetdphimin4"       , &Puppijetmetdphimin4       , "Puppijetmetdphimin4/D");
    tree->Branch("Puppijetmumetdphimin4"       , &Puppijetmumetdphimin4       , "Puppijetmumetdphimin4/D");
    tree->Branch("Puppijetelmetdphimin4"     , &Puppijetelmetdphimin4     , "Puppijetelmetdphimin4/D");
    tree->Branch("Puppijetphmetdphimin4"     , &Puppijetphmetdphimin4     , "Puppijetphmetdphimin4/D");

    tree->Branch("incPuppijetmetdphimin4"    , &incPuppijetmetdphimin4    , "incPuppijetmetdphimin4/D");
    tree->Branch("incPuppijetmumetdphimin4"    , &incPuppijetmumetdphimin4    , "incPuppijetmumetdphimin4/D");
    tree->Branch("incPuppijetelmetdphimin4"  , &incPuppijetelmetdphimin4  , "incPuppijetelmetdphimin4/D");
    tree->Branch("incPuppijetphmetdphimin4"  , &incPuppijetphmetdphimin4  , "incPuppijetphmetdphimin4/D");

    tree->Branch("Puppiht"                   , &Puppiht                   , "Puppiht/D");
    
  }
  
  // Lepton info
  tree->Branch("mu1pid"               , &mu1pid               , "mu1pid/I");
  tree->Branch("mu1pt"                , &mu1pt                , "mu1pt/D");
  tree->Branch("mu1eta"               , &mu1eta               , "mu1eta/D");
  tree->Branch("mu1phi"               , &mu1phi               , "mu1phi/D");
  if(not isTriggerTree){
    tree->Branch("mu1pfpt"              , &mu1pfpt              , "mu1pfpt/D");
    tree->Branch("mu1pfeta"             , &mu1pfeta             , "mu1pfeta/D");
    tree->Branch("mu1pfphi"             , &mu1pfphi             , "mu1pfphi/D");
  }
  tree->Branch("mu1id"                , &mu1id                , "mu1id/I");
  tree->Branch("mu1idm"               , &mu1idm               , "mu1idm/I");
  tree->Branch("mu1idt"               , &mu1idt               , "mu1idt/I");
  tree->Branch("mu1iso"               , &mu1iso               , "mu1iso/D");

  tree->Branch("mu2pid"               , &mu2pid               , "mu2pid/I");
  tree->Branch("mu2pt"                , &mu2pt                , "mu2pt/D");
  tree->Branch("mu2eta"               , &mu2eta               , "mu2eta/D");
  tree->Branch("mu2phi"               , &mu2phi               , "mu2phi/D");
  if(not isTriggerTree){
    tree->Branch("mu2pfpt"              , &mu2pfpt              , "mu2pfpt/D");
    tree->Branch("mu2pfeta"             , &mu2pfeta             , "mu2pfeta/D");
    tree->Branch("mu2pfphi"             , &mu2pfphi             , "mu2pfphi/D");
  }
  tree->Branch("mu2id"                , &mu2id                , "mu2id/I");
  tree->Branch("mu2idm"               , &mu2idm               , "mu2idm/I");
  tree->Branch("mu2idt"               , &mu2idt               , "mu2idt/I");
  tree->Branch("mu2iso"               , &mu2iso               , "mu2iso/D");

  tree->Branch("el1pid"               , &el1pid               , "el1pid/I");
  tree->Branch("el1pt"                , &el1pt                , "el1pt/D");
  tree->Branch("el1eta"               , &el1eta               , "el1eta/D");
  tree->Branch("el1phi"               , &el1phi               , "el1phi/D");
  tree->Branch("el1id"                , &el1id                , "el1id/I");
  tree->Branch("el1idl"               , &el1idl               , "el1idl/I");

  tree->Branch("el2pid"               , &el2pid               , "el2pid/I");
  tree->Branch("el2pt"                , &el2pt                , "el2pt/D");
  tree->Branch("el2eta"               , &el2eta               , "el2eta/D");
  tree->Branch("el2phi"               , &el2phi               , "el2phi/D");
  tree->Branch("el2id"                , &el2id                , "el2id/I");
  tree->Branch("el2idl"               , &el2idl               , "el2idl/I");

  if(not isTriggerTree){
    tree->Branch("tau1pid"               , &tau1pid               , "tau1pid/I");
    tree->Branch("tau1pt"                , &tau1pt                , "tau1pt/D");
    tree->Branch("tau1eta"               , &tau1eta               , "tau1eta/D");
    tree->Branch("tau1phi"               , &tau1phi               , "tau1phi/D");
    tree->Branch("tau1iso"               , &tau1iso               , "tau1iso/D");
    tree->Branch("tau1m"               , &tau1m               , "tau1m/D");
    
    tree->Branch("tau2pid"               , &tau2pid               , "tau2pid/I");
    tree->Branch("tau2pt"                , &tau2pt                , "tau2pt/D");
    tree->Branch("tau2eta"               , &tau2eta               , "tau2eta/D");
    tree->Branch("tau2phi"               , &tau2phi               , "tau2phi/D");
    tree->Branch("tau2iso"               , &tau2iso               , "tau2iso/D");
    tree->Branch("tau2m"               , &tau2m               , "tau2m/D");
  }

    // Dilepton info
  if(not isTriggerTree){
    tree->Branch("zmass"                , &zmass                , "zmass/D");
    tree->Branch("zpt"                  , &zpt                  , "zpt/D");
    tree->Branch("zeta"                 , &zeta                 , "zeta/D");
    tree->Branch("zphi"                 , &zphi                 , "zphi/D");
    tree->Branch("wmt"                  , &wmt                  , "wmt/D");
    tree->Branch("emumass"              , &emumass              , "emumass/D");
    tree->Branch("emupt"                , &emupt                , "emupt/D");
    tree->Branch("emueta"               , &emueta               , "emueta/D");
    tree->Branch("emuphi"               , &emuphi               , "emuphi/D");
    tree->Branch("zeemass"              , &zeemass              , "zeemass/D");
    tree->Branch("zeept"                , &zeept                , "zeept/D");
    tree->Branch("zeeeta"               , &zeeeta               , "zeeeta/D");
    tree->Branch("zeephi"               , &zeephi               , "zeephi/D");
    tree->Branch("wemt"                 , &wemt                 , "wemt/D");
    tree->Branch("zttmass"              , &zttmass              , "zttmass/D");
    tree->Branch("zttpt"                , &zttpt                , "zttept/D");
    tree->Branch("ztteta"               , &ztteta               , "ztteta/D");
    tree->Branch("zttphi"               , &zttphi               , "zttphi/D");
    tree->Branch("wtmt"                 , &wtmt                 , "wtmt/D");
    tree->Branch("taumumass"              , &taumumass              , "taumumass/D");
    tree->Branch("taumupt"                , &taumupt                , "taumupt/D");
    tree->Branch("taumueta"               , &taumueta               , "taumueta/D");
    tree->Branch("taumuphi"               , &taumuphi               , "taumuphi/D");
    tree->Branch("tauemass"              , &tauemass              , "tauemass/D");
    tree->Branch("tauept"                , &tauept                , "tauept/D");
    tree->Branch("taueeta"               , &taueeta               , "taueeta/D");
    tree->Branch("tauephi"               , &tauephi               , "tauephi/D");
  }
  // Photon info
  tree->Branch("phidl"                , &phidl                , "phidl/I");
  tree->Branch("phidm"                , &phidm                , "phidm/I");
  tree->Branch("phidt"                , &phidt                , "phidt/I");
  tree->Branch("phidh"                , &phidh                , "phidh/I");
  tree->Branch("phpt"                 , &phpt                 , "phpt/D");
  tree->Branch("pheta"                , &pheta                , "pheta/D");
  tree->Branch("phphi"                , &phphi                , "phphi/D");

  if(addPhotonPurity and not isTriggerTree){

    tree->Branch("nphotonsPurity"  , &nphotonsPurity  , "nphotonsPurity/I");
    tree->Branch("rho"             , &rho             , "rho/D");
    tree->Branch("phPHiso"         , &phPHiso         , "phPHiso/D");
    tree->Branch("phCHiso"         , &phCHiso         , "phCHiso/D");
    tree->Branch("phPuritypt"      , &phPuritypt      , "phPuritypt/D");
    tree->Branch("pPurityheta"     , &phPurityeta     , "phPurityeta/D");
    tree->Branch("phPurityphi"     , &phPurityphi     , "phPurityphi/D");
    tree->Branch("phPurityPHiso"   , &phPurityPHiso   , "phPurityPHiso/D");
    tree->Branch("phPurityRND04PHiso"   , &phPurityRND04PHiso    , "phPurityRND04PHiso/D");
    tree->Branch("phPurityRND08PHiso"   , &phPurityRND08PHiso    , "phPurityRND08PHiso/D");
    tree->Branch("phPurityCHiso"        , &phPurityCHiso         , "phPurityCHiso/D");
    tree->Branch("phPurityRND04CHiso"        , &phPurityRND04CHiso         , "phPurityRND04CHiso/D");
    tree->Branch("phPurityRND08CHiso"        , &phPurityRND08CHiso         , "phPurityRND08CHiso/D");
    tree->Branch("phNHiso"        , &phNHiso         , "phNHiso/D");
    tree->Branch("phPurityNHiso"        , &phPurityNHiso         , "phPurityNHiso/D");
    tree->Branch("phPuritysieie"        , &phPuritysieie         , "phPuritysieie/D");
    tree->Branch("phPurityhoe"          , &phPurityhoe           , "phPurityhoe/D");
    tree->Branch("phPurityEA"           , &phPurityEA            , "phPurityEA/D");
    tree->Branch("phPurityEAEGamma"     , &phPurityEAEGamma      , "phPurityEAEGamma/D");
    tree->Branch("phPurityElectronVeto" , &phPurityElectronVeto  , "phPurityElectronVeto/D");

  }
  
  // W/Z gen-level info: leptonic and hadronic
  if(not isTriggerTree){
    tree->Branch("wzid"                 , &wzid                 , "wzid/I");
    tree->Branch("wzmass"               , &wzmass               , "wzmass/D");
    tree->Branch("wzmt"                 , &wzmt                 , "wzmt/D");
    tree->Branch("wzpt"                 , &wzpt                 , "wzpt/D");
    tree->Branch("wzeta"                , &wzeta                , "wzeta/D");
    tree->Branch("wzphi"                , &wzphi                , "wzphi/D");
    tree->Branch("wzmothid"             , &wzmothid             , "wzmothid/D");
    tree->Branch("ismatch"              , &ismatch              , "ismatch/I");
    tree->Branch("isdirect"             , &isdirect             , "isdirect/I");
    
    tree->Branch("l1id"                 , &l1id                 , "l1id/I");
    tree->Branch("l1pt"                 , &l1pt                 , "l1pt/D");
    tree->Branch("l1eta"                , &l1eta                , "l1eta/D");
    tree->Branch("l1phi"                , &l1phi                , "l1phi/D");
    tree->Branch("l2id"                 , &l2id                 , "l2id/I");
    tree->Branch("l2pt"                 , &l2pt                 , "l2pt/D");
    tree->Branch("l2eta"                , &l2eta                , "l2eta/D");
    tree->Branch("l2phi"                , &l2phi                , "l2phi/D");
    
    tree->Branch("wzid_h"                 , &wzid_h                 , "wzid_h/I");
    tree->Branch("wzmass_h"               , &wzmass_h               , "wzmass_h/D");
    tree->Branch("wzmt_h"                 , &wzmt_h                 , "wzmt_h/D");
    tree->Branch("wzpt_h"                 , &wzpt_h                 , "wzpt_h/D");
    tree->Branch("wzeta_h"                , &wzeta_h                , "wzeta_h/D");
    tree->Branch("wzphi_h"                , &wzphi_h                , "wzphi_h/D");
    tree->Branch("q1id"                 , &q1id                 , "q1id/I");
    tree->Branch("q1pt"                 , &q1pt                 , "q1pt/D");
    tree->Branch("q1eta"                , &q1eta                , "q1eta/D");
    tree->Branch("q1phi"                , &q1phi                , "q1phi/D");
    tree->Branch("q2id"                 , &q2id                 , "q2id/I");
    tree->Branch("q2pt"                 , &q2pt                 , "q2pt/D");
    tree->Branch("q2eta"                , &q2eta                , "q2eta/D");
    tree->Branch("q2phi"                , &q2phi                , "q2phi/D");
    
    // Top info
    tree->Branch("topmass"               , &topmass               , "topmass/D");
    tree->Branch("toppt"                 , &toppt                 , "toppt/D");
    tree->Branch("topeta"                , &topeta                , "topeta/D");
    tree->Branch("topphi"                , &topphi                , "topphi/D");
    tree->Branch("atopmass"               , &atopmass               , "atopmass/D");
    tree->Branch("atoppt"                 , &atoppt                 , "atoppt/D");
    tree->Branch("atopeta"                , &atopeta                , "atopeta/D");
    tree->Branch("atopphi"                , &atopphi                , "atopphi/D");
    
    // photon gen info
    tree->Branch("parid"                , &parid                , "parid/I");
    tree->Branch("parpt"                , &parpt                , "parpt/D");
    tree->Branch("pareta"               , &pareta               , "pareta/D");
    tree->Branch("parphi"               , &parphi               , "parphi/D");
    tree->Branch("ancid"                , &ancid                , "ancid/I");
    tree->Branch("ancpt"                , &ancpt                , "ancpt/D");
    tree->Branch("anceta"               , &anceta               , "anceta/D");
    tree->Branch("ancphi"               , &ancphi               , "ancphi/D");
    
    // DM mediator
    tree->Branch("dmmass",&dmmass,"dmmass/D");
    tree->Branch("dmpt",&dmpt,"dmpt/D");
    tree->Branch("dmeta",&dmeta,"dmeta/D");
    tree->Branch("dmphi",&dmphi,"dmphi/D");
    tree->Branch("dmid",&dmid,"dmid/I");
    
    // DM particles
    tree->Branch("dmX1id",&dmX1id,"dmX1id/I");
    tree->Branch("dmX1pt",&dmX1pt,"dmX1pt/D");
    tree->Branch("dmX1eta",&dmX1eta,"dmX1eta/D");
    tree->Branch("dmX1phi",&dmX1phi,"dmX1phi/D");
    tree->Branch("dmX1mass",&dmX1mass,"dmX1mass/D");
    
    tree->Branch("dmX2id",&dmX2id,"dmX2id/I");
    tree->Branch("dmX2pt",&dmX2pt,"dmX2pt/D");
    tree->Branch("dmX2eta",&dmX2eta,"dmX2eta/D");
    tree->Branch("dmX2phi",&dmX2phi,"dmX2phi/D");
    tree->Branch("dmX2mass",&dmX2mass,"dmX2mass/D");
    
    // sample info: mediator and DM mass, useful for fast sim                                                                                                                     
    tree->Branch("samplemedM",   &samplemedM, "samplemedM/D");
    tree->Branch("sampledmM",    &sampledmM, "sampledmM/D");
  }
  // AK8 Puppi jets                                                                                                                                                             
  if(addSubstructureCHS and not isTriggerTree){

    tree->Branch("boostedJetpt",  "std::vector<double>", &boostedJetpt);
    tree->Branch("boostedJeteta", "std::vector<double>", &boostedJeteta);
    tree->Branch("boostedJetphi", "std::vector<double>", &boostedJetphi);
    tree->Branch("boostedJetm",   "std::vector<double>", &boostedJetm);

    tree->Branch("boostedJetGenpt",  "std::vector<double>", &boostedJetGenpt);
    tree->Branch("boostedJetGeneta", "std::vector<double>", &boostedJetGeneta);
    tree->Branch("boostedJetGenphi", "std::vector<double>", &boostedJetGenphi);
    tree->Branch("boostedJetGenm",   "std::vector<double>", &boostedJetGenm);

    tree->Branch("boostedJetHFlav", "std::vector<double>", &boostedJetHFlav);
    tree->Branch("boostedJetPFlav", "std::vector<double>", &boostedJetPFlav);
    tree->Branch("boostedJetQGL",   "std::vector<double>", &boostedJetQGL);
    tree->Branch("boostedJetBtag",  "std::vector<double>", &boostedJetBtag);
    tree->Branch("boostedJetDoubleBtag", "std::vector<double>", &boostedJetDoubleBtag);

    tree->Branch("boostedJettau1", "std::vector<double>", &boostedJettau1);
    tree->Branch("boostedJettau2", "std::vector<double>", &boostedJettau2);
    tree->Branch("boostedJettau3", "std::vector<double>", &boostedJettau3);
    tree->Branch("boostedJettau4", "std::vector<double>", &boostedJettau4);

    tree->Branch("boostedJetGentau1", "std::vector<double>", &boostedJetGentau1);
    tree->Branch("boostedJetGentau2", "std::vector<double>", &boostedJetGentau2);
    tree->Branch("boostedJetGentau3", "std::vector<double>", &boostedJetGentau3);
    tree->Branch("boostedJetGentau4", "std::vector<double>", &boostedJetGentau4);

    tree->Branch("prunedJetpt",    "std::vector<double>", &prunedJetpt);
    tree->Branch("prunedJeteta",   "std::vector<double>", &prunedJeteta);
    tree->Branch("prunedJetphi",   "std::vector<double>", &prunedJetphi);
    tree->Branch("prunedJetm",     "std::vector<double>", &prunedJetm);
    tree->Branch("prunedJetm_v2",   "std::vector<double>", &prunedJetm_v2);
    tree->Branch("prunedJetpt_v2",  "std::vector<double>", &prunedJetpt_v2);
    tree->Branch("prunedJeteta_v2", "std::vector<double>", &prunedJeteta_v2);
    tree->Branch("prunedJetphi_v2", "std::vector<double>", &prunedJetphi_v2);
    tree->Branch("prunedJetptraw", "std::vector<double>", &prunedJetptraw);
    tree->Branch("prunedJetmraw",  "std::vector<double>", &prunedJetmraw);

    tree->Branch("prunedJetGenpt",  "std::vector<double>", &prunedJetGenpt);
    tree->Branch("prunedJetGeneta", "std::vector<double>", &prunedJetGeneta);
    tree->Branch("prunedJetGenphi", "std::vector<double>", &prunedJetGenphi);
    tree->Branch("prunedJetGenm",   "std::vector<double>", &prunedJetGenm);

    tree->Branch("prunedJetHFlav", "std::vector<double>", &prunedJetHFlav);
    tree->Branch("prunedJetPFlav", "std::vector<double>", &prunedJetPFlav);
    tree->Branch("prunedJetQGL",   "std::vector<double>", &prunedJetQGL);
    tree->Branch("prunedJetBtag",  "std::vector<double>", &prunedJetBtag);
    tree->Branch("prunedJetDoubleBtag", "std::vector<double>", &prunedJetDoubleBtag);


    tree->Branch("softDropJetpt",    "std::vector<double>", &softDropJetpt);
    tree->Branch("softDropJeteta",   "std::vector<double>", &softDropJeteta);
    tree->Branch("softDropJetphi",   "std::vector<double>", &softDropJetphi);
    tree->Branch("softDropJetm",     "std::vector<double>", &softDropJetm);
    tree->Branch("softDropJetm_v2",  "std::vector<double>", &softDropJetm_v2);
    tree->Branch("softDropJetpt_v2",  "std::vector<double>", &softDropJetpt_v2);
    tree->Branch("softDropJeteta_v2",  "std::vector<double>", &softDropJeteta_v2);
    tree->Branch("softDropJetphi_v2",  "std::vector<double>", &softDropJetphi_v2);
    tree->Branch("softDropJetptraw", "std::vector<double>", &softDropJetptraw);
    tree->Branch("softDropJetmraw",  "std::vector<double>", &softDropJetmraw);

    tree->Branch("softDropJetGenpt",  "std::vector<double>", &softDropJetGenpt);
    tree->Branch("softDropJetGeneta", "std::vector<double>", &softDropJetGeneta);
    tree->Branch("softDropJetGenphi", "std::vector<double>", &softDropJetGenphi);
    tree->Branch("softDropJetGenm",   "std::vector<double>", &softDropJetGenm);

    tree->Branch("softDropJetHFlav", "std::vector<double>", &softDropJetHFlav);
    tree->Branch("softDropJetPFlav", "std::vector<double>", &softDropJetPFlav);
    tree->Branch("softDropJetQGL",   "std::vector<double>", &softDropJetQGL);
    tree->Branch("softDropJetBtag",  "std::vector<double>", &softDropJetBtag);
    tree->Branch("softDropJetDoubleBtag", "std::vector<double>", &softDropJetDoubleBtag);

    tree->Branch("prunedSubJetpt_1",  "std::vector<double>",  &prunedSubJetpt_1);
    tree->Branch("prunedSubJeteta_1", "std::vector<double>",  &prunedSubJeteta_1);
    tree->Branch("prunedSubJetphi_1", "std::vector<double>",  &prunedSubJetphi_1);
    tree->Branch("prunedSubJetm_1",   "std::vector<double>", &prunedSubJetm_1);
    tree->Branch("prunedSubJetGenpt_1","std::vector<double>",  &prunedSubJetGenpt_1);
    tree->Branch("prunedSubJetGenm_1", "std::vector<double>", &prunedSubJetGenm_1);
    tree->Branch("prunedSubJetGeneta_1", "std::vector<double>", &prunedSubJetGeneta_1);
    tree->Branch("prunedSubJetGenphi_1", "std::vector<double>", &prunedSubJetGenphi_1);
    tree->Branch("prunedSubJetHFlav_1",  "std::vector<double>", &prunedSubJetHFlav_1);
    tree->Branch("prunedSubJetPFlav_1",  "std::vector<double>", &prunedSubJetPFlav_1);
    tree->Branch("prunedSubJetQGL_1",    "std::vector<double>", &prunedSubJetQGL_1);
    tree->Branch("prunedSubJetBtag_1",   "std::vector<double>", &prunedSubJetBtag_1);
    tree->Branch("prunedSubJetptraw_1",  "std::vector<double>", &prunedSubJetptraw_1);
    tree->Branch("prunedSubJetmraw_1",   "std::vector<double>", &prunedSubJetmraw_1);
    tree->Branch("prunedSubJetBtagSF_1",   "std::vector<double>", &prunedSubJetBtagSF_1);
    tree->Branch("prunedSubJetBtagSFUp_1",   "std::vector<double>", &prunedSubJetBtagSFUp_1);
    tree->Branch("prunedSubJetBtagSFDown_1",   "std::vector<double>", &prunedSubJetBtagSFDown_1);

    tree->Branch("prunedSubJetpt_2",  "std::vector<double>",  &prunedSubJetpt_2);
    tree->Branch("prunedSubJeteta_2", "std::vector<double>",  &prunedSubJeteta_2);
    tree->Branch("prunedSubJetphi_2", "std::vector<double>",  &prunedSubJetphi_2);
    tree->Branch("prunedSubJetm_2",   "std::vector<double>", &prunedSubJetm_2);
    tree->Branch("prunedSubJetGenpt_2","std::vector<double>",  &prunedSubJetGenpt_2);
    tree->Branch("prunedSubJetGenm_2", "std::vector<double>", &prunedSubJetGenm_2);
    tree->Branch("prunedSubJetGeneta_2", "std::vector<double>", &prunedSubJetGeneta_2);
    tree->Branch("prunedSubJetGenphi_2", "std::vector<double>", &prunedSubJetGenphi_2);
    tree->Branch("prunedSubJetHFlav_2",  "std::vector<double>", &prunedSubJetHFlav_2);
    tree->Branch("prunedSubJetPFlav_2", "std::vector<double>", &prunedSubJetPFlav_2);
    tree->Branch("prunedSubJetQGL_2",   "std::vector<double>", &prunedSubJetQGL_2);
    tree->Branch("prunedSubJetBtag_2",  "std::vector<double>", &prunedSubJetBtag_2);
    tree->Branch("prunedSubJetptraw_2", "std::vector<double>", &prunedSubJetptraw_2);
    tree->Branch("prunedSubJetmraw_2",  "std::vector<double>", &prunedSubJetmraw_2);
    tree->Branch("prunedSubJetBtagSF_2",   "std::vector<double>", &prunedSubJetBtagSF_2);
    tree->Branch("prunedSubJetBtagSFUp_2",   "std::vector<double>", &prunedSubJetBtagSFUp_2);
    tree->Branch("prunedSubJetBtagSFDown_2",   "std::vector<double>", &prunedSubJetBtagSFDown_2);


    tree->Branch("softDropSubJetpt_1","std::vector<double>",  &softDropSubJetpt_1);
    tree->Branch("softDropSubJeteta_1","std::vector<double>",  &softDropSubJeteta_1);
    tree->Branch("softDropSubJetphi_1","std::vector<double>",  &softDropSubJetphi_1);
    tree->Branch("softDropSubJetm_1", "std::vector<double>", &softDropSubJetm_1);
    tree->Branch("softDropSubJetGenpt_1","std::vector<double>",  &softDropSubJetGenpt_1);
    tree->Branch("softDropSubJetGenm_1", "std::vector<double>", &softDropSubJetGenm_1);
    tree->Branch("softDropSubJetGeneta_1", "std::vector<double>", &softDropSubJetGeneta_1);
    tree->Branch("softDropSubJetGenphi_1", "std::vector<double>", &softDropSubJetGenphi_1);
    tree->Branch("softDropSubJetHFlav_1", "std::vector<double>", &softDropSubJetHFlav_1);
    tree->Branch("softDropSubJetPFlav_1", "std::vector<double>", &softDropSubJetPFlav_1);
    tree->Branch("softDropSubJetQGL_1", "std::vector<double>", &softDropSubJetQGL_1);
    tree->Branch("softDropSubJetBtag_1", "std::vector<double>", &softDropSubJetBtag_1);
    tree->Branch("softDropSubJetptraw_1", "std::vector<double>", &softDropSubJetptraw_1);
    tree->Branch("softDropSubJetmraw_1", "std::vector<double>", &softDropSubJetmraw_1);
    tree->Branch("softDropSubJetBtagSF_1",   "std::vector<double>", &softDropSubJetBtagSF_1);
    tree->Branch("softDropSubJetBtagSFUp_1",   "std::vector<double>", &softDropSubJetBtagSFUp_1);
    tree->Branch("softDropSubJetBtagSFDown_1",   "std::vector<double>", &softDropSubJetBtagSFDown_1);

    tree->Branch("softDropSubJetpt_2","std::vector<double>",  &softDropSubJetpt_2);
    tree->Branch("softDropSubJeteta_2","std::vector<double>",  &softDropSubJeteta_2);
    tree->Branch("softDropSubJetphi_2","std::vector<double>",  &softDropSubJetphi_2);
    tree->Branch("softDropSubJetm_2", "std::vector<double>", &softDropSubJetm_2);
    tree->Branch("softDropSubJetGenpt_2","std::vector<double>",  &softDropSubJetGenpt_2);
    tree->Branch("softDropSubJetGenm_2", "std::vector<double>", &softDropSubJetGenm_2);
    tree->Branch("softDropSubJetGeneta_2", "std::vector<double>", &softDropSubJetGeneta_2);
    tree->Branch("softDropSubJetGenphi_2", "std::vector<double>", &softDropSubJetGenphi_2);
    tree->Branch("softDropSubJetHFlav_2", "std::vector<double>", &softDropSubJetHFlav_2);
    tree->Branch("softDropSubJetPFlav_2", "std::vector<double>", &softDropSubJetPFlav_2);
    tree->Branch("softDropSubJetQGL_2", "std::vector<double>", &softDropSubJetQGL_2);
    tree->Branch("softDropSubJetBtag_2", "std::vector<double>", &softDropSubJetBtag_2);
    tree->Branch("softDropSubJetptraw_2", "std::vector<double>", &softDropSubJetptraw_2);
    tree->Branch("softDropSubJetmraw_2", "std::vector<double>", &softDropSubJetmraw_2);
    tree->Branch("softDropSubJetBtagSF_2",   "std::vector<double>", &softDropSubJetBtagSF_2);
    tree->Branch("softDropSubJetBtagSFUp_2",   "std::vector<double>", &softDropSubJetBtagSFUp_2);
    tree->Branch("softDropSubJetBtagSFDown_2",   "std::vector<double>", &softDropSubJetBtagSFDown_2);
  }

  if(addSubstructurePuppi and not isTriggerTree){

    tree->Branch("boostedPuppiJetpt", "std::vector<double>", &boostedPuppiJetpt);
    tree->Branch("boostedPuppiJeteta", "std::vector<double>", &boostedPuppiJeteta);
    tree->Branch("boostedPuppiJetphi", "std::vector<double>", &boostedPuppiJetphi);
    tree->Branch("boostedPuppiJetm", "std::vector<double>", &boostedPuppiJetm);

    tree->Branch("boostedPuppiJetGenpt", "std::vector<double>", &boostedPuppiJetGenpt);
    tree->Branch("boostedPuppiJetGeneta", "std::vector<double>", &boostedPuppiJetGeneta);
    tree->Branch("boostedPuppiJetGenphi", "std::vector<double>", &boostedPuppiJetGenphi);
    tree->Branch("boostedPuppiJetGenm", "std::vector<double>", &boostedPuppiJetGenm);

    tree->Branch("boostedPuppiJetHFlav", "std::vector<double>", &boostedPuppiJetHFlav);
    tree->Branch("boostedPuppiJetPFlav", "std::vector<double>", &boostedPuppiJetPFlav);
    tree->Branch("boostedPuppiJetQGL", "std::vector<double>", &boostedPuppiJetQGL);
    tree->Branch("boostedPuppiJetBtag", "std::vector<double>", &boostedPuppiJetBtag);
    tree->Branch("boostedPuppiJetDoubleBtag", "std::vector<double>", &boostedPuppiJetDoubleBtag);

    tree->Branch("boostedPuppiJettau1", "std::vector<double>", &boostedPuppiJettau1);
    tree->Branch("boostedPuppiJettau2", "std::vector<double>", &boostedPuppiJettau2);
    tree->Branch("boostedPuppiJettau3", "std::vector<double>", &boostedPuppiJettau3);
    tree->Branch("boostedPuppiJettau4", "std::vector<double>", &boostedPuppiJettau4);

    tree->Branch("boostedPuppiJetGentau1", "std::vector<double>", &boostedPuppiJetGentau1);
    tree->Branch("boostedPuppiJetGentau2", "std::vector<double>", &boostedPuppiJetGentau2);
    tree->Branch("boostedPuppiJetGentau3", "std::vector<double>", &boostedPuppiJetGentau3);
    tree->Branch("boostedPuppiJetGentau4", "std::vector<double>", &boostedPuppiJetGentau4);

    tree->Branch("prunedPuppiJetpt", "std::vector<double>", &prunedPuppiJetpt);
    tree->Branch("prunedPuppiJeteta", "std::vector<double>", &prunedPuppiJeteta);
    tree->Branch("prunedPuppiJetphi", "std::vector<double>", &prunedPuppiJetphi);
    tree->Branch("prunedPuppiJetm", "std::vector<double>", &prunedPuppiJetm);
    tree->Branch("prunedPuppiJetm_v2", "std::vector<double>", &prunedPuppiJetm_v2);
    tree->Branch("prunedPuppiJetpt_v2", "std::vector<double>", &prunedPuppiJetpt_v2);
    tree->Branch("prunedPuppiJeteta_v2", "std::vector<double>", &prunedPuppiJeteta_v2);
    tree->Branch("prunedPuppiJetphi_v2", "std::vector<double>", &prunedPuppiJetphi_v2);
    tree->Branch("prunedPuppiJetptraw", "std::vector<double>", &prunedPuppiJetptraw);
    tree->Branch("prunedPuppiJetmraw", "std::vector<double>", &prunedPuppiJetmraw);

    tree->Branch("prunedPuppiJetGenpt", "std::vector<double>", &prunedPuppiJetGenpt);
    tree->Branch("prunedPuppiJetGeneta", "std::vector<double>", &prunedPuppiJetGeneta);
    tree->Branch("prunedPuppiJetGenphi", "std::vector<double>", &prunedPuppiJetGenphi);
    tree->Branch("prunedPuppiJetGenm", "std::vector<double>", &prunedPuppiJetGenm);

    tree->Branch("prunedPuppiJetHFlav", "std::vector<double>", &prunedPuppiJetHFlav);
    tree->Branch("prunedPuppiJetPFlav", "std::vector<double>", &prunedPuppiJetPFlav);
    tree->Branch("prunedPuppiJetQGL", "std::vector<double>", &prunedPuppiJetQGL);
    tree->Branch("prunedPuppiJetBtag", "std::vector<double>", &prunedPuppiJetBtag);
    tree->Branch("prunedPuppiJetDoubleBtag", "std::vector<double>", &prunedPuppiJetDoubleBtag);


    tree->Branch("softDropPuppiJetpt", "std::vector<double>", &softDropPuppiJetpt);
    tree->Branch("softDropPuppiJeteta", "std::vector<double>", &softDropPuppiJeteta);
    tree->Branch("softDropPuppiJetphi", "std::vector<double>", &softDropPuppiJetphi);
    tree->Branch("softDropPuppiJetm", "std::vector<double>", &softDropPuppiJetm);
    tree->Branch("softDropPuppiJetm_v2", "std::vector<double>", &softDropPuppiJetm_v2);
    tree->Branch("softDropPuppiJetpt_v2", "std::vector<double>", &softDropPuppiJetpt_v2);
    tree->Branch("softDropPuppiJeteta_v2", "std::vector<double>", &softDropPuppiJeteta_v2);
    tree->Branch("softDropPuppiJetphi_v2", "std::vector<double>", &softDropPuppiJetphi_v2);
    tree->Branch("softDropPuppiJetptraw", "std::vector<double>", &softDropPuppiJetptraw);
    tree->Branch("softDropPuppiJetmraw", "std::vector<double>", &softDropPuppiJetmraw);

    tree->Branch("softDropPuppiJetGenpt", "std::vector<double>", &softDropPuppiJetGenpt);
    tree->Branch("softDropPuppiJetGeneta", "std::vector<double>", &softDropPuppiJetGeneta);
    tree->Branch("softDropPuppiJetGenphi", "std::vector<double>", &softDropPuppiJetGenphi);
    tree->Branch("softDropPuppiJetGenm", "std::vector<double>", &softDropPuppiJetGenm);

    tree->Branch("softDropPuppiJetHFlav", "std::vector<double>", &softDropPuppiJetHFlav);
    tree->Branch("softDropPuppiJetPFlav", "std::vector<double>", &softDropPuppiJetPFlav);
    tree->Branch("softDropPuppiJetQGL", "std::vector<double>", &softDropPuppiJetQGL);
    tree->Branch("softDropPuppiJetBtag", "std::vector<double>", &softDropPuppiJetBtag);
    tree->Branch("softDropPuppiJetDoubleBtag", "std::vector<double>", &softDropPuppiJetDoubleBtag);

    tree->Branch("prunedPuppiSubJetpt_1","std::vector<double>",  &prunedPuppiSubJetpt_1);
    tree->Branch("prunedPuppiSubJeteta_1","std::vector<double>",  &prunedPuppiSubJeteta_1);
    tree->Branch("prunedPuppiSubJetphi_1","std::vector<double>",  &prunedPuppiSubJetphi_1);
    tree->Branch("prunedPuppiSubJetm_1", "std::vector<double>", &prunedPuppiSubJetm_1);
    tree->Branch("prunedPuppiSubJetGenpt_1","std::vector<double>",  &prunedPuppiSubJetGenpt_1);
    tree->Branch("prunedPuppiSubJetGenm_1", "std::vector<double>", &prunedPuppiSubJetGenm_1);
    tree->Branch("prunedPuppiSubJetGeneta_1", "std::vector<double>", &prunedPuppiSubJetGeneta_1);
    tree->Branch("prunedPuppiSubJetGenphi_1", "std::vector<double>", &prunedPuppiSubJetGenphi_1);
    tree->Branch("prunedPuppiSubJetHFlav_1", "std::vector<double>", &prunedPuppiSubJetHFlav_1);
    tree->Branch("prunedPuppiSubJetPFlav_1", "std::vector<double>", &prunedPuppiSubJetPFlav_1);
    tree->Branch("prunedPuppiSubJetQGL_1", "std::vector<double>", &prunedPuppiSubJetQGL_1);
    tree->Branch("prunedPuppiSubJetBtag_1", "std::vector<double>", &prunedPuppiSubJetBtag_1);
    tree->Branch("prunedPuppiSubJetptraw_1", "std::vector<double>", &prunedPuppiSubJetptraw_1);
    tree->Branch("prunedPuppiSubJetmraw_1", "std::vector<double>", &prunedPuppiSubJetmraw_1);
    tree->Branch("prunedPuppiSubJetBtagSF_1", "std::vector<double>", &prunedPuppiSubJetBtagSF_1);
    tree->Branch("prunedPuppiSubJetBtagSFUp_1", "std::vector<double>", &prunedPuppiSubJetBtagSFUp_1);
    tree->Branch("prunedPuppiSubJetBtagSFDown_1", "std::vector<double>", &prunedPuppiSubJetBtagSFDown_1);

    tree->Branch("prunedPuppiSubJetpt_2","std::vector<double>",  &prunedPuppiSubJetpt_2);
    tree->Branch("prunedPuppiSubJeteta_2","std::vector<double>",  &prunedPuppiSubJeteta_2);
    tree->Branch("prunedPuppiSubJetphi_2","std::vector<double>",  &prunedPuppiSubJetphi_2);
    tree->Branch("prunedPuppiSubJetm_2", "std::vector<double>", &prunedPuppiSubJetm_2);
    tree->Branch("prunedPuppiSubJetGenpt_2","std::vector<double>",  &prunedPuppiSubJetGenpt_2);
    tree->Branch("prunedPuppiSubJetGenm_2", "std::vector<double>", &prunedPuppiSubJetGenm_2);
    tree->Branch("prunedPuppiSubJetGeneta_2", "std::vector<double>", &prunedPuppiSubJetGeneta_2);
    tree->Branch("prunedPuppiSubJetGenphi_2", "std::vector<double>", &prunedPuppiSubJetGenphi_2);
    tree->Branch("prunedPuppiSubJetHFlav_2", "std::vector<double>", &prunedPuppiSubJetHFlav_2);
    tree->Branch("prunedPuppiSubJetPFlav_2", "std::vector<double>", &prunedPuppiSubJetPFlav_2);
    tree->Branch("prunedPuppiSubJetQGL_2", "std::vector<double>", &prunedPuppiSubJetQGL_2);
    tree->Branch("prunedPuppiSubJetBtag_2", "std::vector<double>", &prunedPuppiSubJetBtag_2);
    tree->Branch("prunedPuppiSubJetptraw_2", "std::vector<double>", &prunedPuppiSubJetptraw_2);
    tree->Branch("prunedPuppiSubJetmraw_2", "std::vector<double>", &prunedPuppiSubJetmraw_2);
    tree->Branch("prunedPuppiSubJetBtagSF_2", "std::vector<double>", &prunedPuppiSubJetBtagSF_2);
    tree->Branch("prunedPuppiSubJetBtagSFUp_2", "std::vector<double>", &prunedPuppiSubJetBtagSFUp_2);
    tree->Branch("prunedPuppiSubJetBtagSFDown_2", "std::vector<double>", &prunedPuppiSubJetBtagSFDown_2);


    tree->Branch("softDropPuppiSubJetpt_1","std::vector<double>",  &softDropPuppiSubJetpt_1);
    tree->Branch("softDropPuppiSubJeteta_1","std::vector<double>",  &softDropPuppiSubJeteta_1);
    tree->Branch("softDropPuppiSubJetphi_1","std::vector<double>",  &softDropPuppiSubJetphi_1);
    tree->Branch("softDropPuppiSubJetm_1", "std::vector<double>", &softDropPuppiSubJetm_1);
    tree->Branch("softDropPuppiSubJetGenpt_1","std::vector<double>",  &softDropPuppiSubJetGenpt_1);
    tree->Branch("softDropPuppiSubJetGenm_1", "std::vector<double>", &softDropPuppiSubJetGenm_1);
    tree->Branch("softDropPuppiSubJetGeneta_1", "std::vector<double>", &softDropPuppiSubJetGeneta_1);
    tree->Branch("softDropPuppiSubJetGenphi_1", "std::vector<double>", &softDropPuppiSubJetGenphi_1);
    tree->Branch("softDropPuppiSubJetHFlav_1", "std::vector<double>", &softDropPuppiSubJetHFlav_1);
    tree->Branch("softDropPuppiSubJetPFlav_1", "std::vector<double>", &softDropPuppiSubJetPFlav_1);
    tree->Branch("softDropPuppiSubJetQGL_1", "std::vector<double>", &softDropPuppiSubJetQGL_1);
    tree->Branch("softDropPuppiSubJetBtag_1", "std::vector<double>", &softDropPuppiSubJetBtag_1);
    tree->Branch("softDropPuppiSubJetptraw_1", "std::vector<double>", &softDropPuppiSubJetptraw_1);
    tree->Branch("softDropPuppiSubJetmraw_1", "std::vector<double>", &softDropPuppiSubJetmraw_1);
    tree->Branch("softDropPuppiSubJetBtagSF_1", "std::vector<double>", &softDropPuppiSubJetBtagSF_1);
    tree->Branch("softDropPuppiSubJetBtagSFUp_1", "std::vector<double>", &softDropPuppiSubJetBtagSFUp_1);
    tree->Branch("softDropPuppiSubJetBtagSFDown_1", "std::vector<double>", &softDropPuppiSubJetBtagSFDown_1);

    tree->Branch("softDropPuppiSubJetpt_2","std::vector<double>",  &softDropPuppiSubJetpt_2);
    tree->Branch("softDropPuppiSubJeteta_2","std::vector<double>",  &softDropPuppiSubJeteta_2);
    tree->Branch("softDropPuppiSubJetphi_2","std::vector<double>",  &softDropPuppiSubJetphi_2);
    tree->Branch("softDropPuppiSubJetm_2", "std::vector<double>", &softDropPuppiSubJetm_2);
    tree->Branch("softDropPuppiSubJetGenpt_2","std::vector<double>",  &softDropPuppiSubJetGenpt_2);
    tree->Branch("softDropPuppiSubJetGenm_2", "std::vector<double>", &softDropPuppiSubJetGenm_2);
    tree->Branch("softDropPuppiSubJetGeneta_2", "std::vector<double>", &softDropPuppiSubJetGeneta_2);
    tree->Branch("softDropPuppiSubJetGenphi_2", "std::vector<double>", &softDropPuppiSubJetGenphi_2);
    tree->Branch("softDropPuppiSubJetHFlav_2", "std::vector<double>", &softDropPuppiSubJetHFlav_2);
    tree->Branch("softDropPuppiSubJetPFlav_2", "std::vector<double>", &softDropPuppiSubJetPFlav_2);
    tree->Branch("softDropPuppiSubJetQGL_2", "std::vector<double>", &softDropPuppiSubJetQGL_2);
    tree->Branch("softDropPuppiSubJetBtag_2", "std::vector<double>", &softDropPuppiSubJetBtag_2);
    tree->Branch("softDropPuppiSubJetptraw_2", "std::vector<double>", &softDropPuppiSubJetptraw_2);
    tree->Branch("softDropPuppiSubJetmraw_2", "std::vector<double>", &softDropPuppiSubJetmraw_2);
    tree->Branch("softDropPuppiSubJetBtagSF_2", "std::vector<double>", &softDropPuppiSubJetBtagSF_2);
    tree->Branch("softDropPuppiSubJetBtagSFUp_2", "std::vector<double>", &softDropPuppiSubJetBtagSFUp_2);
    tree->Branch("softDropPuppiSubJetBtagSFDown_2", "std::vector<double>", &softDropPuppiSubJetBtagSFDown_2);

  }

  if(addPhotonIDVariables and not isTriggerTree){
    if(not addPhotonPurity)
      tree->Branch("rho"             , &rho             , "rho/D");
    tree->Branch("photonPt", "std::vector<double>", &photonPt);
    tree->Branch("photonEta", "std::vector<double>", &photonEta);
    tree->Branch("photonPhi", "std::vector<double>", &photonPhi);
    tree->Branch("photonE", "std::vector<double>", &photonE);
    tree->Branch("photonSCEta", "std::vector<double>", &photonSCEta);
    tree->Branch("photonSCPhi", "std::vector<double>", &photonSCPhi);
    tree->Branch("photonSCEnergy", "std::vector<double>", &photonSCEnergy);
    tree->Branch("photonSCRawEnergy", "std::vector<double>", &photonSCRawEnergy);
    tree->Branch("photonHOverE", "std::vector<double>", &photonHOverE);
    tree->Branch("photonSigmaIetaIeta", "std::vector<double>", &photonSigmaIetaIeta);
    tree->Branch("photonChargedIso", "std::vector<double>", &photonChargedIso);
    tree->Branch("photonNeutralIso", "std::vector<double>", &photonNeutralIso);
    tree->Branch("photonEMIso", "std::vector<double>", &photonEMIso);
    tree->Branch("photonElectronVeto", "std::vector<double>", &photonElectronVeto);
  }
  if(addElectronIDVariables and not isTriggerTree){
    if(not addPhotonPurity)
      tree->Branch("rho"             , &rho             , "rho/D");
    tree->Branch("electronPt", "std::vector<double>", &electronPt);
    tree->Branch("electronEta", "std::vector<double>", &electronEta);
    tree->Branch("electronPhi", "std::vector<double>", &electronPhi);
    tree->Branch("electronE", "std::vector<double>", &electronE);
    tree->Branch("electronSCEta", "std::vector<double>", &electronSCEta);
    tree->Branch("electronSCPhi", "std::vector<double>", &electronSCPhi);
    tree->Branch("electronSCEnergy", "std::vector<double>", &electronSCEnergy);
    tree->Branch("electronSCRawEnergy", "std::vector<double>", &electronSCRawEnergy);
    tree->Branch("electronHOverE", "std::vector<double>", &electronHOverE);
    tree->Branch("electronSigmaIetaIeta", "std::vector<double>", &electronSigmaIetaIeta);
    tree->Branch("electronChargedIso", "std::vector<double>", &electronChargedIso);
    tree->Branch("electronNeutralIso", "std::vector<double>", &electronNeutralIso);
    tree->Branch("electronEMIso", "std::vector<double>", &electronEMIso);
    tree->Branch("electronGsfPt", "std::vector<double>", &electronGsfPt);
    tree->Branch("electronDphi", "std::vector<double>", &electronDphi);
    tree->Branch("electronDeta", "std::vector<double>", &electronDeta);
    tree->Branch("electronEOP", "std::vector<double>", &electronEOP);
    tree->Branch("electronMissHit", "std::vector<double>", &electronMissHit);
    tree->Branch("electronConversion", "std::vector<double>", &electronConversion);
    tree->Branch("electronDz", "std::vector<double>", &electronDz);
    tree->Branch("electronDxy", "std::vector<double>", &electronDxy);
  }
}

void MonoJetTreeMaker::endJob() {}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  // triggers for the Analysis
  triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight"); //0
  triggerPathsVector.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight"); //1
  triggerPathsVector.push_back("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight"); //2
  triggerPathsVector.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");//3
  triggerPathsVector.push_back("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");//4
  triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight"); //5
  triggerPathsVector.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight"); //6
  triggerPathsVector.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"); //7
  triggerPathsVector.push_back("HLT_PFMET90_PFMHT90_IDTight"); //8
  triggerPathsVector.push_back("HLT_PFMET100_PFMHT100_IDTight_v"); //9 
  triggerPathsVector.push_back("HLT_PFMET110_PFMHT110_IDTight_v");//10 
  triggerPathsVector.push_back("HLT_PFMET120_PFMHT120_IDTight"); //11
  triggerPathsVector.push_back("HLT_PFMET170_NoiseCleaned"); //12
  triggerPathsVector.push_back("HLT_PFMET170_JetIdCleaned"); //13
  triggerPathsVector.push_back("HLT_PFMET170_HBHECleaned"); //14
  triggerPathsVector.push_back("HLT_PFMET170_v"); //15
  triggerPathsVector.push_back("HLT_PFMET300_NoiseCleaned"); //16
  triggerPathsVector.push_back("HLT_PFMET300_JetIdCleaned"); //17
  triggerPathsVector.push_back("HLT_PFMET300_v"); //18
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight"); //19
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight"); //20
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight");    //21
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_v");//22
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v");//23
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight"); //24
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight"); //25
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"); //26
  triggerPathsVector.push_back("HLT_Photon165_HE10"); //27
  triggerPathsVector.push_back("HLT_Photon175");      //28
  triggerPathsVector.push_back("HLT_Photon120_v");    //29
  triggerPathsVector.push_back("HLT_Photon90_v");     //30
  triggerPathsVector.push_back("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v");     //31
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); //31
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); //32
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); //33
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); //34
  triggerPathsVector.push_back("HLT_IsoMu20_v"); //35
  triggerPathsVector.push_back("HLT_IsoMu22_v"); //36
  triggerPathsVector.push_back("HLT_IsoMu24_v"); //37
  triggerPathsVector.push_back("HLT_IsoTkMu20"); //38
  triggerPathsVector.push_back("HLT_IsoTkMu22"); //39
  triggerPathsVector.push_back("HLT_IsoTkMu24"); //40
  triggerPathsVector.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"); //41
  triggerPathsVector.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"); //42
  triggerPathsVector.push_back("HLT_Ele24_eta2p1_WPLoose_Gsf_v"); //43
  triggerPathsVector.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v"); //44
  triggerPathsVector.push_back("HLT_Ele27_WPTight_Gsf_v"); //45
  triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v"); //46 
  triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPTight_Gsf_v"); //47
  triggerPathsVector.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v"); //48
  triggerPathsVector.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v"); //49
  triggerPathsVector.push_back("HLT_PFHT400_v");//50
  triggerPathsVector.push_back("HLT_PFHT475_v");//51
  triggerPathsVector.push_back("HLT_PFHT600_v");//52
  triggerPathsVector.push_back("HLT_PFHT650_v");//53
  triggerPathsVector.push_back("HLT_PFHT800_v");//54
  triggerPathsVector.push_back("HLT_PFHT900_v");//55
  triggerPathsVector.push_back("HLT_ECALHT800_v");//56
  triggerPathsVector.push_back("HLT_Photon90_CaloIdL_PFHT500_v");//57
  triggerPathsVector.push_back("HLT_Photon90_CaloIdL_PFHT600_v");//58

  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);
  
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }

  bool changedHLTPSP = false;
  hltPrescaleProvider_->init(iRun, iSetup, triggerResultsTag.process(), changedHLTPSP); //ND 

  // MET filter Paths
  filterPathsVector.push_back("Flag_CSCTightHalo2015Filter");
  filterPathsVector.push_back("Flag_HBHENoiseFilter");
  filterPathsVector.push_back("Flag_eeBadScFilter");
  filterPathsVector.push_back("Flag_HBHENoiseIsoFilter");
  filterPathsVector.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filterPathsVector.push_back("Flag_goodVertices");
  filterPathsVector.push_back("Flag_globalTightHalo2016Filter");
  
  HLTConfigProvider fltrConfig;
  fltrConfig.init(iRun, iSetup, filterResultsTag.process(), changedConfig);
  
  for (size_t i = 0; i < filterPathsVector.size(); i++) {
    filterPathsMap[filterPathsVector[i]] = -1;
  }

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
  if(isMC and uselheweights){    
    edm::Handle<LHERunInfoProduct> run;
    iRun.getByLabel(lheRunTag,run);
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    if(xsec < 0)
      xsec = myLHERunInfoProduct.heprup().XSECUP.at(0);
    

    using namespace boost::algorithm;

    if(isSignalSample){

      for (auto iter = myLHERunInfoProduct.headers_begin(); iter != myLHERunInfoProduct.headers_end(); iter++){
	std::vector<std::string> lines = iter->lines();
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) { 
	  std::vector<std::string> tokens;
	  if(lines.at(iLine).find("DMmass") !=std::string::npos){// powheg mono-j
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    sampledmM = std::stod(tokens.at(1));
	  }
	  else if(lines.at(iLine).find("DMVmass") !=std::string::npos){// powheg mono-j
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    samplemedM = std::stod(tokens.at(1));
	  }
	  else if(lines.at(iLine).find("import model") !=std::string::npos){ // madgraph mono-V                                                                               
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    std::vector<std::string> subtokens;
	    split(subtokens,tokens.at(2),is_any_of("_"));
	    samplemedM = std::stod(subtokens.at(3));
	    sampledmM = std::stod(subtokens.at(4));
	  }
	  else if(lines.at(iLine).find("Resonance:") != std::string::npos){ // JHUGen --> only resonance mass (mediator) .. dM fixed in the event loop                     
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    samplemedM = std::stod(tokens.at(3));
	    sampledmM  = -1.;
	    readDMFromGenParticle = true;
	  }
	}
      }   
    }       
  }
}

void MonoJetTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {}

// to fill b-tag SF
void MonoJetTreeMaker::calculateBtagSF(const pat::Jet & jet, const std::string & algorithm, 
				       std::vector<double> & scalefactor, std::vector<double> & scalefactorUp, std::vector<double> & scalefactorDown){

  if(algorithm != "CSV" and algorithm != "MVA" and algorithm != "SubCSV") return;
  // bounds for CSVv2 and MVAv2
  double jetPt = jet.pt();
  double jetEta = jet.eta();

  if(algorithm == "CSV"){

    if(jet.hadronFlavour() == 5){
      scalefactor.push_back(bMediumCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));	      
    }    
    else if(jet.hadronFlavour() == 4){
      scalefactor.push_back(bMediumCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_C,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_C,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_C,jetEta,jetPt));	      
    }
    else{
      scalefactor.push_back(bMediumCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
    }	      
  }
  else if(algorithm == "MVA"){
    if(jet.hadronFlavour() == 5){            
      scalefactor.push_back(bMediumMVA.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumMVA.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumMVA.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));	      
    }
    else if(jet.hadronFlavour() == 4){
      scalefactor.push_back(bMediumMVA.back().eval_auto_bounds("central",BTagEntry::FLAV_C,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumMVA.back().eval_auto_bounds("up",BTagEntry::FLAV_C,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumMVA.back().eval_auto_bounds("down",BTagEntry::FLAV_C,jetEta,jetPt));	      
      
    }
    else{      
      scalefactor.push_back(bMediumMVA.back().eval_auto_bounds("central",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumMVA.back().eval_auto_bounds("up",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumMVA.back().eval_auto_bounds("down",BTagEntry::FLAV_UDSG,jetEta,jetPt));	      
    }
  }
  else if(algorithm == "SubCSV"){
    if(jet.hadronFlavour() == 5){            
      scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));	      
    }
    else if(jet.hadronFlavour() == 4){
      scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));	      
    }
    else{      
      scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));	      
      scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));	      
    }
  }
}


//This code is ripped off from https://github.com/krav/ElectronWork/blob/master/ElectronNtupler/plugins/PhotonNtuplerMiniAOD.cc
void MonoJetTreeMaker::findFirstNonPhotonMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {

  if (particle == 0)
    return;
  
  if (abs(particle->pdgId()) == 22) 
    findFirstNonPhotonMother(particle->mother(0), ancestorid, ancestorpt, ancestoreta, ancestorphi);
  else {
    ancestorid  = particle->pdgId();
    ancestorpt  = particle->pt();
    ancestoreta = particle->eta();
    ancestorphi = particle->phi();
  }
  return;
}

// for photons
void MonoJetTreeMaker::findMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {
  
  if (particle == 0) 
    return;

  if (abs(particle->pdgId()) == 22) {
    ancestorid  = particle->pdgId();
    ancestorpt  = particle->pt();
    ancestoreta = particle->eta();
    ancestorphi = particle->phi();
  }
  return;
}

// compute muon isolation
double MonoJetTreeMaker::computeMuonIso(const reco::Muon& mu) {

    double isoval = mu.pfIsolationR04().sumNeutralHadronEt;
    isoval += mu.pfIsolationR04().sumPhotonEt;
    isoval -= 0.5*mu.pfIsolationR04().sumPUPt;
    if (isoval < 0.) isoval = 0.;
    isoval += mu.pfIsolationR04().sumChargedHadronPt;
    isoval /= mu.pt();            

    return isoval;
}

// apply standard CMS jet ID
bool MonoJetTreeMaker::applyJetID(const pat::Jet & jet, const std::string & level){

  if(level != "loose" and level != "tight" and level != "tightLepVeto")
    return true;
   
  bool passjetid = false;

  //apply a loose jet id https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
  if(level == "loose"){ 
    if (fabs(jet.eta()) <= 2.7 and
	jet.neutralHadronEnergyFraction() < 0.99 and
	jet.neutralEmEnergyFraction()     < 0.99 and
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
     if (fabs(jet.eta()) > 2.4)
       passjetid = true;
     else if (fabs(jet.eta()) <= 2.4 and
	      jet.chargedHadronEnergyFraction() > 0. and
	      jet.chargedEmEnergyFraction()     < 0.99 and 
	      jet.chargedMultiplicity()         > 0) 
       passjetid = true;
   }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) <= 3.0 and
	     jet.neutralEmEnergyFraction() < 0.9 and
	     jet.neutralMultiplicity()     > 2) 
      passjetid = true;  
    else if(fabs(jet.eta()) > 3.0 and
	    jet.neutralEmEnergyFraction() < 0.9 and
	    jet.neutralMultiplicity()     > 10)
      passjetid = true; 
  } 
  else if(level == "tight"){
    
    if (fabs(jet.eta()) <= 2.7 && 
	jet.neutralHadronEnergyFraction() < 0.90 && 
	jet.neutralEmEnergyFraction()     < 0.90 && 
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4) 
	passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 && 
	       jet.chargedHadronEnergyFraction() > 0. && 
	       jet.chargedEmEnergyFraction() < 0.99 && 
	       jet.chargedMultiplicity() > 0) 
	passjetid = true;
    }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) < 3.0 
	     && jet.neutralEmEnergyFraction() < 0.9 
	     && jet.neutralMultiplicity() > 2) 
      passjetid = true;
    else if(fabs(jet.eta()) > 3.0 and
	    jet.neutralEmEnergyFraction() < 0.9 and
	    jet.neutralMultiplicity() > 10)
      passjetid = true;    
  }
  
  else if(level == "tightLepVeto"){
    if (fabs(jet.eta()) <= 2.7 &&
        jet.neutralHadronEnergyFraction() < 0.90 &&
        jet.neutralEmEnergyFraction() < 0.90 &&
	jet.muonEnergyFraction() < 0.80 && 
        (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4)
        passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 &&
               jet.chargedHadronEnergyFraction() > 0. &&
               jet.chargedEmEnergyFraction() < 0.90 &&
               jet.chargedMultiplicity() > 0)
	passjetid = true;
    }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) < 3.0
	     && jet.neutralEmEnergyFraction() < 0.9
	     && jet.neutralMultiplicity() > 2)
      passjetid = true;
    else if (fabs(jet.eta()) > 3.0
	     && jet.neutralEmEnergyFraction() < 0.9
	     && jet.neutralMultiplicity() > 10)
      passjetid = true;    
    
  }  
  return passjetid;  
}

bool MonoJetTreeMaker::applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi){

  bool passpuid    = false;
  double puidval   = 0;
  double jetabseta = fabs(jet.eta());
  double jetpt     = jet.pt();

  if(jet.hasUserFloat("puid:fullDiscriminant"))
    puidval = jet.userFloat("puid:fullDiscriminant");
  else if(jet.hasUserFloat("puidPuppi:fullDiscriminant"))
    puidval = jet.userFloat("puidPuppi:fullDiscriminant");
  else if(jet.hasUserFloat("pileupJetId:fullDiscriminant"))
    puidval = jet.userFloat("pileupJetId:fullDiscriminant");
  else 
    return true;

  // from twiki: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID --> to be loaded in GT soon
  if(level == "loose"){

    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > -0.89) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.52) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.38) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.30) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;

  }
  else if(level == "medium"){ 

    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > 0.61) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.23) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.17) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;

  }
  else if(level == "tight"){
    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > 0.86) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.10) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.01) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.03) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;
  }

  return passpuid;
}

double MonoJetTreeMaker::computeDR(const reco::Candidate *genPart,pat::PhotonRef phot){
  double dR = 999.;
  TLorentzVector phop4;
  phop4.SetPtEtaPhiM(phot->pt(),phot->eta(), phot->phi(),0);
  TLorentzVector p4;
  p4.SetPtEtaPhiM(genPart->pt(),genPart->eta(),genPart->phi(),0);
  dR = phop4.DeltaR(p4);
  return dR;
}

void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

double MonoJetTreeMaker::getChargedHadronEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.000000000001;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return  0.000000000001;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return  0.000000000001;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return  0.000000000001;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return  0.000000000001;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return  0.000000000001;
  else if (fabs(eta) >= 2.4) return 0.000000000001 ;
  else return 0.;
}

double MonoJetTreeMaker::getNeutralHadronEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.0599;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0819;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0696;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.036;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.036;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.0462;
  else if (fabs(eta) >= 2.4) return 0.0656;
  else return 0.;
}

double MonoJetTreeMaker::getGammaEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.1271;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.1101;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0756;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.1175;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.1498;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.1857;
  else if (fabs(eta) >= 2.4) return 0.2183;
  else return 0.;
}
double MonoJetTreeMaker::getGammaNewEAForPhotonIso(double eta) {
  if (fabs(eta) < 0.9) return 0.17;
  else if (fabs(eta) >= 0.9   && fabs(eta) < 1.4442) return 0.14;
  else if (fabs(eta) >= 1.4442 && fabs(eta) < 2.0  ) return 0.0320;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0512;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0766;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.0949;
  else if (fabs(eta) >= 2.4) return 0.116;
  else return 0.;
}

////////
void MonoJetTreeMaker::fillTriggerObjects(const edm::Handle<pat::TriggerObjectStandAloneCollection> & triggerObjectsH, 
					  const edm::TriggerNames & trignames) {

  // Initialize 
  trig_obj_n = 0;
  trig_obj_pt.clear();
  trig_obj_eta.clear();
  trig_obj_phi.clear();
  trig_obj_col.clear();
  std::string trgColl = "";
  int iObj = 0;
  
  if(triggerObjectsH.isValid()) {
    
    // Loop over trigger objects --> dump all the trigger objects
    for (pat::TriggerObjectStandAlone obj : *triggerObjectsH) {       
      iObj++ ;
      obj.unpackPathNames(trignames);
      
      // trigger objec
      trig_obj_pt.push_back( obj.pt());
      trig_obj_eta.push_back(obj.eta());
      trig_obj_phi.push_back(obj.phi());

      // collection name
      trgColl = obj.collection();
      trig_obj_col.push_back(trgColl);

      // Increment trigger object index
      trig_obj_n++ ;
    }
  }  
}


/////////
void MonoJetTreeMaker::fillAlgosL1(const edm::Event & iEvent,
				   const edm::EventSetup & eventSetup,
				   const edm::Handle<GlobalAlgBlkBxCollection> & H_L1Algos){
  
  // Initialize
  trig_L1A_check = 0;
  trig_L1A_n     = 0;
  trig_L1A_list.clear();

  // Check Handle validity
  if(!H_L1Algos.isValid()) {
    trig_L1A_check = -1;
    return;
  }

  // Get L1 Menu and utils
  edm::ESHandle<L1TUtmTriggerMenu> menu;
  eventSetup.get<L1TUtmTriggerMenuRcd>().get(menu);
  int iErrorCode = -1;
  L1GtUtils::TriggerCategory trigCategory = L1GtUtils::AlgorithmTrigger;
  // Get L1 utils for prescales
  L1GtUtils const & l1GtUtils = hltPrescaleProvider_->l1GtUtils();
  l1GtUtils.prescaleFactorSetIndex(iEvent, trigCategory, iErrorCode);
  // Get the bit/name association //
  const UInt_t nBits = 512;
  std::string algoBitToName[nBits];
  //
  for (auto const & keyval: menu->getAlgorithmMap()) { 
    std::string const & trigName  = keyval.second.getName(); 
    unsigned int index = keyval.second.getIndex(); 
    if(index<nBits) algoBitToName[index] = trigName;
  } // end algo Map

  // Get the L1 decision per algo //
  GlobalAlgBlk const &result = H_L1Algos->at(0,0);
  //
  for (unsigned int itrig = 0; itrig < result.maxPhysicsTriggers; ++itrig) {
    // Check decision for this bit
    bool myflag = result.getAlgoDecisionFinal(itrig) ; 
    if(myflag ) { 
      trig_L1A_list.push_back(algoBitToName[itrig]);
      trig_L1A_n++ ;
    }
  } // end loop: L1 trigger results  
}


void MonoJetTreeMaker::fillTriggerL1(const edm::Handle<l1t::EGammaBxCollection> & H_L1EG,  const edm::Handle<l1t::TauBxCollection>  & H_L1Tau,
				     const edm::Handle<l1t::JetBxCollection>    & H_L1Jet, const edm::Handle<l1t::MuonBxCollection> & H_L1Mu,
				     const edm::Handle<l1t::EtSumBxCollection>  & H_L1Sums) {

  // Initialize
  trig_L1EG_pt  .clear(); trig_L1EG_eta  .clear(); trig_L1EG_phi  .clear(); 
  trig_L1Jet_pt .clear(); trig_L1Jet_eta .clear(); trig_L1Jet_phi .clear(); 
  trig_L1Mu_pt  .clear(); trig_L1Mu_eta  .clear(); trig_L1Mu_phi  .clear(); 
  trig_L1ETM_pt = trig_L1ETM_phi = trig_L1HTM_pt  = trig_L1HTM_phi = 0; 
  trig_L1ETT_pt = trig_L1ETT_phi = trig_L1HTT_pt  = trig_L1HTT_phi = 0; 

  int sumType   = -1;
  double minL1EG  = 20;
  double minL1Jet = 50;
  double minL1Mu  = 15;

  // L1 EG    
  if(H_L1EG.isValid()) {
    for (l1t::EGammaBxCollection::const_iterator it=H_L1EG->begin(); it!=H_L1EG->end(); it++){
      if(it->pt() < minL1EG) continue;
      trig_L1EG_pt .push_back( it->pt()  );
      trig_L1EG_eta.push_back( it->eta() );
      trig_L1EG_phi.push_back( it->phi() );
    }
  }
  
  // L1 Jet
  if(H_L1Jet.isValid()) {
    for (l1t::JetBxCollection::const_iterator it=H_L1Jet->begin(); it!=H_L1Jet->end(); it++){
      if(it->pt() < minL1Jet) continue;
      trig_L1Jet_pt .push_back( it->pt()  );
      trig_L1Jet_eta.push_back( it->eta() );
      trig_L1Jet_phi.push_back( it->phi() );
    }
  }

  // L1 Mu
  if(H_L1Mu.isValid()) {
    for (l1t::MuonBxCollection::const_iterator it=H_L1Mu->begin(); it!=H_L1Mu->end(); it++){
      if(it->pt() < minL1Mu) continue;
      trig_L1Mu_pt .push_back( it->pt()  );
      trig_L1Mu_eta.push_back( it->eta() );
      trig_L1Mu_phi.push_back( it->phi() );
    }
  }
  
  // L1 Sums
  if(H_L1Sums.isValid()) {
    for (l1t::EtSumBxCollection::const_iterator it=H_L1Sums->begin(); it!=H_L1Sums->end(); it++){
      
      sumType = static_cast<int>( it->getType() );
      if(sumType == l1t::EtSum::kTotalEt){
	  trig_L1ETT_pt  = it->et();
	  trig_L1ETT_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kTotalHt){
	  trig_L1HTT_pt  = it->et();
	  trig_L1HTT_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kMissingEt){
	  trig_L1ETM_pt  = it->et();
	  trig_L1ETM_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kMissingHt){
	  trig_L1HTM_pt  = it->et();
	  trig_L1HTM_phi = it->phi();
      }      
    }    
  }
}
  

DEFINE_FWK_MODULE(MonoJetTreeMaker);
