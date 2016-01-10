#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <TH1F.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

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
  
  void findMother(const reco::Candidate*, int &, double &, double &, double &);
  void findFirstNonPhotonMother(const reco::Candidate*, int &, double &, double &, double &);
  double computeMuonIso(const reco::Muon&);

  bool applyJetID(const pat::Jet & jet, const std::string & level);
  bool applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi);

  // Gen Particles
  const bool isMC;
  const bool uselheweights;
  const bool isWorZMCSample;
  const bool isSignalSample;   
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
  edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
  edm::EDGetTokenT<LHEEventProduct>                  lheInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
  double xsec;
  
  // InputTags
  const edm::InputTag triggerResultsTag;
  const edm::InputTag filterResultsTag;
  
  // trgger and filter tokens
  const edm::EDGetTokenT<edm::TriggerResults>              triggerResultsToken;
  const edm::EDGetTokenT<edm::TriggerResults>              filterResultsToken;
  const edm::EDGetTokenT<bool>                             hbhelooseToken;
  const edm::EDGetTokenT<bool>                             hbhetightToken;
  const edm::EDGetTokenT<bool>                             hbheisoToken;
  
  // Vertex
  const edm::EDGetTokenT<std::vector<reco::Vertex> >       verticesToken;

  // muons
  const edm::EDGetTokenT<pat::MuonRefVector>               muonsToken;
  const edm::EDGetTokenT<pat::MuonRefVector>               tightmuonsToken;
  const edm::EDGetTokenT<pat::MuonRefVector>               highptmuonsToken;

  // electrons
  const edm::EDGetTokenT<pat::ElectronRefVector>           electronsToken;
  const edm::EDGetTokenT<pat::ElectronRefVector>           tightelectronsToken;
  const edm::EDGetTokenT<pat::ElectronRefVector>           heepelectronsToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> >             electronLooseIdToken;  

  // Photons
  const edm::EDGetTokenT<pat::PhotonRefVector>             photonsToken;
  const edm::EDGetTokenT<pat::PhotonRefVector>             tightphotonsToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> >             photonLooseIdToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> >             photonMediumIdToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> >             photonTightIdToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> >             photonHighPtIdToken;

  // Taus
  const edm::EDGetTokenT<std::vector<pat::Tau> >           tausToken;

  //Jets AK4
  const edm::EDGetTokenT<std::vector<pat::Jet> >           jetsToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >           puppijetsToken;
  const bool addPuppiJets;

  // MET
  const edm::EDGetTokenT<edm::View<pat::MET> >             t1metToken;
  const edm::EDGetTokenT<edm::View<pat::MET> >             t1mumetToken;
  const edm::EDGetTokenT<edm::View<pat::MET> >             t1elemetToken;
  const edm::EDGetTokenT<edm::View<pat::MET> >             t1phmetToken;

  // Puppi MET
  const bool addPuppiMET;
  edm::EDGetTokenT<edm::View<pat::MET> >             puppit1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> >             puppit1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >             puppit1elemetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >             puppit1phmetToken;

  // MET systematics
  const bool addMETSystematics;

  // MVA met
  const bool addMVAMet;
  edm::EDGetTokenT<edm::View<reco::MET> >           mvaMETToken;

  // inner bools
  const bool applyHLTFilter;
  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;   

  // Jet AK8
  const bool addSubstructureCHS;
  const bool addSubstructurePuppi;
  edm::EDGetTokenT<std::vector<pat::Jet> >           boostedJetsToken;
  TString boostedJetsCHSLabel;
  edm::EDGetTokenT<std::vector<pat::Jet> >           boostedPuppiJetsToken;
  TString boostedJetsPuppiLabel;
  

  // inner vectors
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  std::vector<std::string> filterPathsVector;
  std::map<std::string, int> filterPathsMap;

  // tree
  std::auto_ptr<TTree> tree;

  // pileup info
  int32_t puobs, putrue; 

  // lepton info
  int32_t wzid, l1id, l2id;
  int32_t mu1pid, mu2pid, mu1id, mu2id, mu1idm, mu2idm, mu1idt, mu2idt, mu1hptid, mu2hptid;
  int32_t el1pid, el2pid, el1id, el1idl, el2id, el2idl;
  int32_t phidl, phidm, phidt, phidh, parid, ancid; 

  // event info
  uint32_t event, run, lumi;  
  uint32_t nvtx;
  uint32_t nmuons, ntightmuons, nhighptmuons;
  uint32_t nelectrons, ntightelectrons, nheepelectrons;
  uint32_t ntaus, nphotons;
  uint32_t njets, nbjets, nbjetslowpt;
  uint32_t npuppijets, npuppibjets, npuppibjetslowpt;

  // trigger and met filters flags 
  uint8_t  hltmet90, hltmet120, hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300;
  uint8_t  hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel;
  uint8_t  flagcsctight, flaghbhenoise, flaghbheloose, flaghbhetight, flaghbheiso, flageebadsc;

  // PF MET info (typeI and Raw)
  double t1pfmet, t1pfmetphi, t1mumet, t1mumetphi, t1elmet, t1elmetphi, t1phmet, t1phmetphi;
  double pfmet, pfmetphi, mumet, mumetphi, elmet, elmetphi, phmet, phmetphi;

  // Puppi MET info (typeI and Raw)
  double  puppipfmet, puppipfmetphi, puppimumet, puppimumetphi, puppielmet, puppielmetphi, puppiphmet, puppiphmetphi;
  double  puppit1pfmet, puppit1pfmetphi, puppit1mumet, puppit1mumetphi, puppit1elmet, puppit1elmetphi, puppit1phmet, puppit1phmetphi;

  // mva met
  double mvamet, mvametphi;

  // gen met
  double genmet, genmetphi;

  // met systematics
  double t1pfmetMuEnUp, t1pfmetMuEnDown, t1pfmetElEnUp, t1pfmetElEnDown, t1pfmetPhoEnUp, t1pfmetPhoEnDown, t1pfmetTauEnUp, t1pfmetTauEnDown;
  double t1pfmetJetEnUp, t1pfmetJetEnDown, t1pfmetJetResUp, t1pfmetJetResDown, t1pfmetUncEnUp, t1pfmetUncEnDown;

  // met systematics puppi
  double puppit1pfmetMuEnUp, puppit1pfmetMuEnDown, puppit1pfmetElEnUp, puppit1pfmetElEnDown, puppit1pfmetPhoEnUp, puppit1pfmetPhoEnDown, puppit1pfmetTauEnUp;
  double puppit1pfmetTauEnDown, puppit1pfmetJetEnUp, puppit1pfmetJetEnDown, puppit1pfmetJetResUp, puppit1pfmetJetResDown, puppit1pfmetUncEnUp, puppit1pfmetUncEnDown;
  
  // AK4 CHS jets
  double   leadingjetpt, leadingjeteta, leadingjetphi, leadingjetm; 
  double   signaljetpt , signaljeteta , signaljetphi , signaljetm, signaljetbtag, signaljetRawpt;
  double   signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi;
  double   signaljetHFlav, signaljetPFlav, signaljetQGL, signaljetPUID;
  double   signaljetGenpt, signaljetGeneta, signaljetGenphi, signaljetGenm;
  //
  double   secondjetpt , secondjeteta , secondjetphi , secondjetm, secondjetbtag, secondjetRawpt;
  double   secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetCEMfrac, secondjetmetdphi;
  double   secondjetHFlav, secondjetPFlav, secondjetQGL, secondjetPUID;
  double   secondjetGenpt, secondjetGeneta, secondjetGenphi, secondjetGenm;
  //
  double   thirdjetpt , thirdjeteta , thirdjetphi , thirdjetm, thirdjetbtag, thirdjetRawpt;
  double   thirdjetCHfrac, thirdjetNHfrac, thirdjetEMfrac, thirdjetCEMfrac, thirdjetmetdphi;
  double   thirdjetHFlav, thirdjetPFlav, thirdjetQGL, thirdjetPUID, thirdjetGenpt, thirdjetGeneta, thirdjetGenphi, thirdjetGenm;
  //
  double   fourthjetpt , fourthjeteta , fourthjetphi , fourthjetm, fourthjetbtag;
  double   fourthjetCHfrac, fourthjetNHfrac, fourthjetEMfrac, fourthjetCEMfrac, fourthjetmetdphi;
  double   fourthjetHFlav, fourthjetPFlav, fourthjetQGL, fourthjetPUID, fourthjetGenpt, fourthjetGeneta, fourthjetGenphi, fourthjetGenm, fourthjetRawpt;
  //
  double   jetmetdphimin , incjetmetdphimin , jetmumetdphimin , incjetmumetdphimin, jetelmetdphimin , incjetelmetdphimin , jetphmetdphimin , incjetphmetdphimin , jetjetdphi;
  double   jetmetdphimin4, incjetmetdphimin4, jetmumetdphimin4, incjetmumetdphimin4 , jetelmetdphimin4, incjetelmetdphimin4, jetphmetdphimin4, incjetphmetdphimin4, ht; 

  // Puppijet ak4 puppi
  double   leadingPuppijetpt, leadingPuppijeteta, leadingPuppijetphi, leadingPuppijetm; 
  // leading
  double   signalPuppijetpt , signalPuppijeteta , signalPuppijetphi , signalPuppijetm, signalPuppijetbtag;
  double   signalPuppijetCHfrac, signalPuppijetNHfrac, signalPuppijetEMfrac, signalPuppijetCEMfrac, signalPuppijetmetdphi;
  double   signalPuppijetHFlav, signalPuppijetPFlav, signalPuppijetQGL, signalPuppijetPUID;
  double   signalPuppijetGenpt, signalPuppijetGeneta, signalPuppijetGenphi, signalPuppijetGenm, signalPuppijetRawpt;
  //
  double   secondPuppijetpt , secondPuppijeteta , secondPuppijetphi , secondPuppijetm, secondPuppijetbtag;
  double   secondPuppijetCHfrac, secondPuppijetNHfrac, secondPuppijetEMfrac, secondPuppijetCEMfrac, secondPuppijetmetdphi;
  double   secondPuppijetHFlav, secondPuppijetPFlav, secondPuppijetQGL, secondPuppijetPUID;
  double   secondPuppijetGenpt, secondPuppijetGeneta, secondPuppijetGenphi, secondPuppijetGenm, secondPuppijetRawpt;
  //
  double   thirdPuppijetpt , thirdPuppijeteta , thirdPuppijetphi , thirdPuppijetm, thirdPuppijetbtag;
  double   thirdPuppijetCHfrac, thirdPuppijetNHfrac, thirdPuppijetEMfrac, thirdPuppijetCEMfrac, thirdPuppijetmetdphi;
  double   thirdPuppijetHFlav, thirdPuppijetPFlav, thirdPuppijetQGL, thirdPuppijetPUID;
  double   thirdPuppijetGenpt, thirdPuppijetGeneta, thirdPuppijetGenphi, thirdPuppijetGenm, thirdPuppijetRawpt;
  //
  double   fourthPuppijetpt , fourthPuppijeteta , fourthPuppijetphi , fourthPuppijetm, fourthPuppijetbtag;
  double   fourthPuppijetCHfrac, fourthPuppijetNHfrac, fourthPuppijetEMfrac, fourthPuppijetCEMfrac, fourthPuppijetmetdphi;
  double   fourthPuppijetHFlav, fourthPuppijetPFlav, fourthPuppijetQGL, fourthPuppijetPUID;
  double   fourthPuppijetGenpt, fourthPuppijetGeneta, fourthPuppijetGenphi, fourthPuppijetGenm, fourthPuppijetRawpt;
  //
  double   Puppijetmetdphimin , incPuppijetmetdphimin , Puppijetmumetdphimin , incPuppijetmumetdphimin , Puppijetelmetdphimin , incPuppijetelmetdphimin , Puppijetphmetdphimin , incPuppijetphmetdphimin , PuppijetPuppijetdphi;
  double   Puppijetmetdphimin4, incPuppijetmetdphimin4, Puppijetmumetdphimin4, incPuppijetmumetdphimin4, Puppijetelmetdphimin4, incPuppijetelmetdphimin4, Puppijetphmetdphimin4, incPuppijetphmetdphimin4, Puppiht; 
  
  // AK8 CHS jets
  double   leadBoostedJetpt, leadBoostedJeteta, leadBoostedJetphi, leadBoostedJetm;
  double   leadBoostedJetGenpt, leadBoostedJetGenm, leadBoostedJetGeneta, leadBoostedJetGenphi;
  double   leadBoostedJettau1, leadBoostedJettau2, leadBoostedJettau3, leadBoostedJettau4, leadBoostedJetecf1, leadBoostedJetecf2, leadBoostedJetecf3;
  double   leadBoostedJetHFlav, leadBoostedJetPFlav, leadBoostedJetQGL, leadBoostedJetBtag, leadBoostedJetDoubleBtag;
  double   leadBoostedJetBosonpt, leadBoostedJetBosoneta, leadBoostedJetBosonphi, leadBoostedJetBosonm;
  double   leadPrunedJetpt, leadPrunedJetm, leadPrunedJetphi,leadPrunedJeteta,  leadPrunedJetGenpt, leadPrunedJetGenm, leadPrunedJetGeneta,leadPrunedJetGenphi;
  double   leadPrunedJetptraw, leadPrunedJetmraw;
  double   leadPrunedJetHFlav, leadPrunedJetPFlav, leadPrunedJetQGL, leadPrunedJetBtag, leadPrunedJetDoubleBtag;
  double   leadSoftDropJetpt, leadSoftDropJetm, leadSoftDropJeteta, leadSoftDropJetphi, leadSoftDropJetGenpt, leadSoftDropJetGenm,leadSoftDropJetGeneta, leadSoftDropJetGenphi;
  double   leadSoftDropJetHFlav, leadSoftDropJetPFlav, leadSoftDropJetQGL, leadSoftDropJetBtag, leadSoftDropJetDoubleBtag;
  double   leadSoftDropJetptraw, leadSoftDropJetmraw;
  double   leadPrunedSubJetpt_1, leadPrunedSubJetm_1,leadPrunedSubJetphi_1, leadPrunedSubJeteta_1, leadPrunedSubJetHFlav_1, leadPrunedSubJetQGL_1, leadPrunedSubJetBtag_1;
  double   leadPrunedSubJetGenpt_1, leadPrunedSubJetGenm_1,leadPrunedSubJetGeneta_1, leadPrunedSubJetGenphi_1, leadPrunedSubJetPFlav_1;
  double   leadPrunedSubJetptraw_1, leadPrunedSubJetmraw_1;
  double   leadPrunedSubJetpt_2, leadPrunedSubJetm_2,leadPrunedSubJetphi_2, leadPrunedSubJeteta_2, leadPrunedSubJetHFlav_2, leadPrunedSubJetQGL_2, leadPrunedSubJetBtag_2;
  double   leadPrunedSubJetGenpt_2, leadPrunedSubJetGenm_2, leadPrunedSubJetGeneta_2, leadPrunedSubJetGenphi_2, leadPrunedSubJetPFlav_2;
  double   leadPrunedSubJetptraw_2, leadPrunedSubJetmraw_2;
  double   leadSoftDropSubJetpt_1, leadSoftDropSubJetm_1,leadSoftDropSubJetphi_1, leadSoftDropSubJeteta_1;
  double   leadSoftDropSubJetHFlav_1, leadSoftDropSubJetQGL_1, leadSoftDropSubJetBtag_1;
  double   leadSoftDropSubJetGenpt_1, leadSoftDropSubJetGenm_1, leadSoftDropSubJetGeneta_1, leadSoftDropSubJetGenphi_1, leadSoftDropSubJetPFlav_1;
  double   leadSoftDropSubJetptraw_1, leadSoftDropSubJetmraw_1;
  double   leadSoftDropSubJetpt_2, leadSoftDropSubJetm_2,leadSoftDropSubJetphi_2, leadSoftDropSubJeteta_2, leadSoftDropSubJetHFlav_2; 
  double   leadSoftDropSubJetQGL_2, leadSoftDropSubJetBtag_2;
  double   leadSoftDropSubJetGenpt_2, leadSoftDropSubJetGenm_2, leadSoftDropSubJetGeneta_2, leadSoftDropSubJetGenphi_2, leadSoftDropSubJetPFlav_2;
  double   leadSoftDropSubJetptraw_2, leadSoftDropSubJetmraw_2;

  // AK8 Puppi jets
  double   leadPuppiBoostedJetpt, leadPuppiBoostedJeteta, leadPuppiBoostedJetphi, leadPuppiBoostedJetm;
  double   leadPuppiBoostedJetGenpt, leadPuppiBoostedJetGenm, leadPuppiBoostedJetGeneta, leadPuppiBoostedJetGenphi;
  double   leadPuppiBoostedJettau1, leadPuppiBoostedJettau2, leadPuppiBoostedJettau3, leadPuppiBoostedJettau4, leadPuppiBoostedJetecf1, leadPuppiBoostedJetecf2, leadPuppiBoostedJetecf3;
  double   leadPuppiBoostedJetHFlav, leadPuppiBoostedJetPFlav, leadPuppiBoostedJetQGL, leadPuppiBoostedJetBtag, leadPuppiBoostedJetDoubleBtag;
  double   leadPuppiBoostedJetBosonpt, leadPuppiBoostedJetBosoneta, leadPuppiBoostedJetBosonphi, leadPuppiBoostedJetBosonm;
  double   leadPuppiPrunedJetpt, leadPuppiPrunedJetm, leadPuppiPrunedJeteta, leadPuppiPrunedJetphi;
  double   leadPuppiPrunedJetptraw, leadPuppiPrunedJetmraw;
  double   leadPuppiPrunedJetGenpt, leadPuppiPrunedJetGenm, leadPuppiPrunedJetGenphi, leadPuppiPrunedJetGeneta;
  double   leadPuppiPrunedJetHFlav, leadPuppiPrunedJetPFlav, leadPuppiPrunedJetQGL, leadPuppiPrunedJetBtag, leadPuppiPrunedJetDoubleBtag;
  double   leadPuppiSoftDropJetpt, leadPuppiSoftDropJetm, leadPuppiSoftDropJeteta, leadPuppiSoftDropJetphi;
  double   leadPuppiSoftDropJetptraw, leadPuppiSoftDropJetmraw;
  double   leadPuppiSoftDropJetGenpt, leadPuppiSoftDropJetGenm, leadPuppiSoftDropJetGeneta, leadPuppiSoftDropJetGenphi;
  double   leadPuppiSoftDropJetHFlav, leadPuppiSoftDropJetPFlav, leadPuppiSoftDropJetQGL, leadPuppiSoftDropJetBtag, leadPuppiSoftDropJetDoubleBtag;

  double   leadPuppiPrunedSubJetpt_1, leadPuppiPrunedSubJetm_1,leadPuppiPrunedSubJetphi_1, leadPuppiPrunedSubJeteta_1, leadPuppiPrunedSubJetHFlav_1;
  double   leadPuppiPrunedSubJetQGL_1, leadPuppiPrunedSubJetBtag_1;
  double   leadPuppiPrunedSubJetptraw_1, leadPuppiPrunedSubJetmraw_1;
  double   leadPuppiPrunedSubJetGenpt_1, leadPuppiPrunedSubJetGenm_1,leadPuppiPrunedSubJetGeneta_1, leadPuppiPrunedSubJetGenphi_1, leadPuppiPrunedSubJetPFlav_1;
  double   leadPuppiPrunedSubJetpt_2, leadPuppiPrunedSubJetm_2,leadPuppiPrunedSubJetphi_2, leadPuppiPrunedSubJeteta_2, leadPuppiPrunedSubJetHFlav_2; 
  double   leadPuppiPrunedSubJetQGL_2, leadPuppiPrunedSubJetBtag_2;
  double   leadPuppiPrunedSubJetGenpt_2, leadPuppiPrunedSubJetGenm_2, leadPuppiPrunedSubJetGeneta_2, leadPuppiPrunedSubJetGenphi_2, leadPuppiPrunedSubJetPFlav_2;
  double   leadPuppiPrunedSubJetptraw_2, leadPuppiPrunedSubJetmraw_2;

  double   leadPuppiSoftDropSubJetpt_1, leadPuppiSoftDropSubJetm_1,leadPuppiSoftDropSubJetphi_1, leadPuppiSoftDropSubJeteta_1, leadPuppiSoftDropSubJetHFlav_1;
  double   leadPuppiSoftDropSubJetQGL_1, leadPuppiSoftDropSubJetBtag_1;
  double   leadPuppiSoftDropSubJetGenpt_1, leadPuppiSoftDropSubJetGenm_1, leadPuppiSoftDropSubJetGeneta_1, leadPuppiSoftDropSubJetGenphi_1, leadPuppiSoftDropSubJetPFlav_1;
  double   leadPuppiSoftDropSubJetptraw_1, leadPuppiSoftDropSubJetmraw_1;
  double   leadPuppiSoftDropSubJetpt_2, leadPuppiSoftDropSubJetm_2,leadPuppiSoftDropSubJetphi_2, leadPuppiSoftDropSubJeteta_2;
  double   leadPuppiSoftDropSubJetHFlav_2, leadPuppiSoftDropSubJetQGL_2, leadPuppiSoftDropSubJetBtag_2;
  double   leadPuppiSoftDropSubJetGenpt_2, leadPuppiSoftDropSubJetGenm_2, leadPuppiSoftDropSubJetGeneta_2, leadPuppiSoftDropSubJetGenphi_2, leadPuppiSoftDropSubJetPFlav_2;
  double   leadPuppiSoftDropSubJetptraw_2, leadPuppiSoftDropSubJetmraw_2;


  // muon info
  double   mu1pt, mu1eta, mu1phi, mu1pterr, mu1pfpt, mu1pfeta, mu1pfphi, mu1iso, mu2pt, mu2eta, mu2phi, mu2pfpt, mu2pfeta, mu2pfphi, mu2iso, mu2pterr;
  // ele info
  double   el1pt, el1eta, el1phi, ele1e, ele1eerr, el2pt, ele2e, ele2eerr, el2eta, el2phi, phpt, pheta, phphi, phe, pheerr;
  // dilepton info
  double   zmass, zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt;
  // gen info
  double   wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi, parpt, pareta, parphi, ancpt, anceta, ancphi;
  // weights
  double   wgt, kfact, puwgt;

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
  
};


MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
  ///////////// GEN INFO
  // isMC or Data
  isMC(iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  // use lhe weights or not
  uselheweights(iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
  // is signal sample or not
  isWorZMCSample(iConfig.existsAs<bool>("isWorZMCSample") ? iConfig.getParameter<bool>("isWorZMCSample") : false),
  isSignalSample(iConfig.existsAs<bool>("isSignalSample") ? iConfig.getParameter<bool>("isSignalSample") : false),
  // qcd and pdf weights
  xsec(iConfig.getParameter<double>("xsec") * 1000.0),
  ///////////// TRIGGER and filter info INFO
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
  filterResultsTag(iConfig.getParameter<edm::InputTag>("filterResults")),
  triggerResultsToken(consumes<edm::TriggerResults> (triggerResultsTag)),
  filterResultsToken(consumes<edm::TriggerResults> (filterResultsTag)),
  hbhelooseToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbheloose"))),
  hbhetightToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbhetight"))),
  hbheisoToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbheiso"))),
  // vertexes
  verticesToken(consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
  //muons
  muonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("muons"))),
  tightmuonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("tightmuons"))),
  highptmuonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("highptmuons"))),
  // electrons
  electronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("electrons"))),
  tightelectronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("tightelectrons"))),
  heepelectronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("heepelectrons"))),
  electronLooseIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("electronLooseId"))),
  // photons
  photonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photons"))),
  tightphotonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("tightphotons"))),
  photonLooseIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonLooseId"))),
  photonMediumIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonMediumId"))),
  photonTightIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonTightId"))),
  photonHighPtIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonHighPtId"))),
  // taus
  tausToken(consumes<std::vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("taus"))),
  // jets AK4
  jetsToken(consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
  addPuppiJets(iConfig.existsAs<bool>("addPuppiJets") ? iConfig.getParameter<bool>("addPuppiJets") : false),
  // met
  t1metToken(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("t1met"))),
  t1mumetToken(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("t1mumet"))),
  t1elemetToken(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("t1elmet"))),
  t1phmetToken(consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("t1phmet"))),
  // puppi met
  addPuppiMET(iConfig.existsAs<bool>("addPuppiMET") ? iConfig.getParameter<bool>("addPuppiMET") : false),
  // MET Systematics
  addMETSystematics(iConfig.existsAs<bool>("addMETSystematics") ? iConfig.getParameter<bool>("addMETSystematics") : false),
  // MVA met
  addMVAMet(iConfig.existsAs<bool>("addMVAMet") ? iConfig.getParameter<bool>("addMVAMet") : false),
  //filter events based on trigger
  applyHLTFilter(iConfig.existsAs<bool>("applyHLTFilter") ? iConfig.getParameter<bool>("applyHLTFilter") : false),
  // booleans
  cleanMuonJet(iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
  cleanElectronJet(iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
  cleanPhotonJet(iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
  // substructure
  addSubstructureCHS(iConfig.existsAs<bool>("addSubstructureCHS") ? iConfig.getParameter<bool>("addSubstructureCHS") : false),
  addSubstructurePuppi(iConfig.existsAs<bool>("addSubstructurePuppi") ? iConfig.getParameter<bool>("addSubstructurePuppi") : false){

  usesResource();
  
  // only for simulated samples

  if( isMC ){
    if(iConfig.existsAs<edm::InputTag>("pileup"))
      pileupInfoToken = consumes<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileup"));
    else
      pileupInfoToken = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));

    if(iConfig.existsAs<edm::InputTag>("genevt"))
      genevtInfoToken = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genevt"));
    else
      genevtInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    
    if(iConfig.existsAs<edm::InputTag>("lheinfo"))
      lheInfoToken = consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("lheinfo"));
    else
      lheInfoToken = consumes<LHEEventProduct> (edm::InputTag("externalLHEProducer"));

    if(iConfig.existsAs<edm::InputTag>("gens"))      
      gensToken = consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("gens"));
    else
      gensToken = consumes<edm::View<reco::GenParticle> > ( edm::InputTag("prunedGenParticles"));

  }
  
  if(addPuppiJets)
    puppijetsToken = consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("puppijets"));

  if(addPuppiMET){
    puppit1metToken    = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1met"));
    puppit1mumetToken  = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1mumet"));
    puppit1elemetToken = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1elmet"));
    puppit1phmetToken  = consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("puppit1phmet"));		      
  }
  
  if(addMVAMet)
    mvaMETToken = consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("mvaMET"));

  if(addSubstructureCHS){
    boostedJetsToken =  consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("boostedJetsCHS"));
    boostedJetsCHSLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsCHS").label().c_str());
    boostedJetsCHSLabel.ReplaceAll("packed","");
    boostedJetsCHSLabel.ReplaceAll("PatJets","");
    boostedJetsCHSLabel.ReplaceAll("patJets","");
  }

  if(addSubstructurePuppi){
    boostedPuppiJetsToken =  consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("boostedJetsPuppi"));
    boostedJetsPuppiLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsPuppi").label().c_str());
    boostedJetsPuppiLabel.ReplaceAll("packed","");
    boostedJetsPuppiLabel.ReplaceAll("PatJets","");
    boostedJetsPuppiLabel.ReplaceAll("patJets","");
  }

}


MonoJetTreeMaker::~MonoJetTreeMaker() {}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    using namespace edm;
    using namespace reco;
    using namespace std;
    using namespace pat;

    // Get handles to all the requisite collections
    // TRIGGER and FILTERS
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByToken(triggerResultsToken, triggerResultsH);
    Handle<TriggerResults> filterResultsH;
    iEvent.getByToken(filterResultsToken, filterResultsH);
    Handle<bool> hbhelooseH;
    iEvent.getByToken(hbhelooseToken, hbhelooseH);
    Handle<bool> hbhetightH;
    iEvent.getByToken(hbhetightToken, hbhetightH);
    Handle<bool> hbheisoH;
    iEvent.getByToken(hbheisoToken, hbheisoH);

    // GEN INFO    
    Handle<vector<PileupSummaryInfo> > pileupInfoH;
    Handle<GenEventInfoProduct> genevtInfoH;
    Handle<LHEEventProduct> lheInfoH;
    Handle<View<GenParticle> > gensH;

    if(isMC){
      iEvent.getByToken(pileupInfoToken, pileupInfoH);
      if (uselheweights){
	iEvent.getByToken(genevtInfoToken, genevtInfoH);
	iEvent.getByToken(lheInfoToken, lheInfoH);
      }
      if (isWorZMCSample || isSignalSample)
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

    Handle<pat::ElectronRefVector> tightelectronsH;
    iEvent.getByToken(tightelectronsToken, tightelectronsH);
    pat::ElectronRefVector tightelectrons = *tightelectronsH;

    Handle<pat::ElectronRefVector> heepelectronsH;
    iEvent.getByToken(heepelectronsToken, heepelectronsH);
    pat::ElectronRefVector heepelectrons = *heepelectronsH;

    Handle<edm::ValueMap<bool> > electronLooseIdH;
    iEvent.getByToken(electronLooseIdToken, electronLooseIdH);

    // PHOTONS
    Handle<pat::PhotonRefVector> photonsH;
    iEvent.getByToken(photonsToken, photonsH);
    pat::PhotonRefVector photons = *photonsH;

    Handle<pat::PhotonRefVector> tightphotonsH;
    iEvent.getByToken(tightphotonsToken, tightphotonsH);
    pat::PhotonRefVector tightphotons = *tightphotonsH;

    Handle<ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdToken, photonLooseIdH);
    Handle<ValueMap<bool> > photonMediumIdH;
    iEvent.getByToken(photonMediumIdToken, photonMediumIdH);
    Handle<ValueMap<bool> > photonTightIdH;
    iEvent.getByToken(photonTightIdToken, photonTightIdH);
    Handle<ValueMap<bool> > photonHighPtIdH;
    iEvent.getByToken(photonHighPtIdToken, photonHighPtIdH);

    // TAUS
    Handle<vector<pat::Tau> > tausH;
    iEvent.getByToken(tausToken, tausH);

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


    Handle<View<reco::MET> > mvaMetH;
    if(addMVAMet)
      iEvent.getByToken(mvaMETToken,mvaMetH);

    // boosted jets
    Handle<vector<pat::Jet> > boostedJetsH;
    if(addSubstructureCHS){
      iEvent.getByToken(boostedJetsToken,boostedJetsH);
    }

    Handle<vector<pat::Jet> > boostedPuppiJetsH;
    if(addSubstructurePuppi)
      iEvent.getByToken(boostedPuppiJetsToken,boostedPuppiJetsH);
    
    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();
    
    // Trigger info
    hltmet90        = 0;
    hltmet120       = 0;
    hltmetwithmu90  = 0;
    hltmetwithmu120 = 0;
    hltmetwithmu170 = 0;
    hltmetwithmu300 = 0;
    hltjetmet90     = 0;
    hltjetmet120    = 0;
    hltphoton165    = 0;
    hltphoton175    = 0;
    hltdoublemu     = 0;
    hltsinglemu     = 0;
    hltdoubleel     = 0;
    hltsingleel     = 0;

    // Which triggers fired
    if(triggerResultsH.isValid()){
      for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;	
	if(triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]]))
        if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
        if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
        if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu90  = 1; // MET trigger
        if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu120 = 1; // MET trigger
        if (i == 8  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 9  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 10 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 11 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
        if (i == 12 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 13 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 14 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
        if (i == 15 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 16 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 17 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet90     = 1; // Jet-MET trigger
        if (i == 18 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 19 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 20 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet120    = 1; // Jet-MET trigger
        if (i == 21 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165    = 1; // Photon trigger
        if (i == 22 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175    = 1; // Photon trigger
        if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 29 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 30 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 32 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      }
    }

    bool triggered = false;
    if (hltmet90        == 1) triggered = true;
    if (hltmet120       == 1) triggered = true;
    if (hltmetwithmu90  == 1) triggered = true;
    if (hltmetwithmu120 == 1) triggered = true;
    if (hltmetwithmu170 == 1) triggered = true;
    if (hltmetwithmu300 == 1) triggered = true;
    if (hltjetmet90     == 1) triggered = true;
    if (hltjetmet120    == 1) triggered = true;
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltdoublemu     == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltdoubleel     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (applyHLTFilter && !triggered) return;

    // MET filter info
    flagcsctight  = 0;
    flaghbhenoise = 0;
    flageebadsc   = 0;

    // HBHE Noise 
    flaghbheloose = (*hbhelooseH ? 1 : 0);
    flaghbhetight = (*hbhetightH ? 1 : 0);
    flaghbheiso   = (*hbheisoH   ? 1 : 0);
    if(filterResultsH.isValid()){
      // Which MET filters passed
      for (size_t i = 0; i < filterPathsVector.size(); i++) {
        if (filterPathsMap[filterPathsVector[i]] == -1) continue;
	if(filterResultsH->accept(filterPathsMap[filterPathsVector[i]]))
        if (i == 0  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagcsctight  = 1; // CSCTightHaloFilter
        if (i == 1  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbhenoise = 1; // HBHENoiseFilter
        if (i == 2  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flageebadsc   = 1; // eeBadScFilter
      }
    }

    // Pileup info -- Will need to the updated to the Run-II specifications
    if(verticesH.isValid())
      nvtx   = verticesH->size();
    else
      nvtx = 0;

    puobs  = 0;
    putrue = 0;
    puwgt  = 1.;
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
      if(t1metH->front().genMET()){
	genmet    = t1metH->front().genMET()->pt();
	genmetphi = t1metH->front().genMET()->phi();
      }
      else {
	genmet = -99.;
	genmetphi = -99.;
      }

      t1pfmet        = t1metH->front().corPt();
      t1pfmetphi     = t1metH->front().corPhi();
      pfmet          = t1metH->front().uncorPt();
      pfmetphi       = t1metH->front().uncorPhi();
    }

    if(addMVAMet && mvaMetH.isValid()){
      mvamet    = mvaMetH->front().pt();
      mvametphi = mvaMetH->front().phi();
    }
    else{
      mvamet = -99.;
      mvametphi = -99.;
    }

    if(t1mumetH.isValid()){
      t1mumet        = t1mumetH->front().corPt();
      t1mumetphi     = t1mumetH->front().corPhi();
      mumet          = t1mumetH->front().uncorPt();
      mumetphi       = t1mumetH->front().uncorPhi();
    }

    if(t1elmetH.isValid()){
      t1elmet        = t1elmetH->front().corPt();
      t1elmetphi     = t1elmetH->front().corPhi();
      elmet          = t1elmetH->front().uncorPt();
      elmetphi       = t1elmetH->front().uncorPhi();
    }

    if(t1phmetH.isValid()){
      t1phmet        = t1phmetH->front().corPt();
      t1phmetphi     = t1phmetH->front().corPhi();
      phmet          = t1phmetH->front().uncorPt();
      phmetphi       = t1phmetH->front().uncorPhi();
    }

    if(addMETSystematics){          
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
      }
    }
    else{
      t1pfmetMuEnUp = -99.; t1pfmetMuEnDown= -99.; t1pfmetElEnUp = -99.; t1pfmetElEnDown = -99.;
      t1pfmetPhoEnUp = -99.; t1pfmetPhoEnDown = -99.; t1pfmetTauEnUp = -99.; t1pfmetTauEnDown = -99.;
      t1pfmetJetEnUp = -99.; t1pfmetJetEnDown = -99.; t1pfmetJetResUp = -99.; t1pfmetJetResDown = -99.;
      t1pfmetUncEnUp = -99.; t1pfmetUncEnDown = -99.;      
    }

    // puppi met info
    
    puppit1pfmet = -99.; puppit1pfmetphi = -99.;
    puppipfmet = -99.; puppipfmetphi = -99.;
    puppit1mumet = -99.; puppit1mumetphi = -99.;
    puppimumet = -99.; puppimumetphi = -99.;
    puppit1elmet = -99.; puppit1elmetphi = -99.;
    puppielmet = -99.; puppielmetphi = -99.;
    puppit1phmet = -99.; puppit1phmetphi = -99.;
    puppiphmet = -99.; puppiphmetphi = -99.;
    
    puppit1pfmetMuEnUp = -99.; puppit1pfmetMuEnDown= -99.; puppit1pfmetElEnUp = -99.; puppit1pfmetElEnDown = -99.;
    puppit1pfmetPhoEnUp = -99.; puppit1pfmetPhoEnDown = -99.; puppit1pfmetTauEnUp = -99.; puppit1pfmetTauEnDown = -99.;
    puppit1pfmetJetEnUp = -99.; puppit1pfmetJetEnDown = -99.; puppit1pfmetJetResUp = -99.; puppit1pfmetJetResDown = -99.;
    puppit1pfmetUncEnUp = -99.; puppit1pfmetUncEnDown = -99.;


    if(addPuppiMET){

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

      if(addMETSystematics){	
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
	}
      }
    }
    
    // AK4 Jet information
    // find leading jets not overlapping with leptons
    leadingjetpt  = 0.0;
    leadingjeteta = 0.0;
    leadingjetphi = 0.0;
    leadingjetm   = 0.0;

    vector<pat::JetRef> incjets;
    vector<pat::JetRef> jets;
    if(jetsH.isValid()){
      
      for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
	
	//clean from leptons
	bool skipjet = false;
	if(muonsH.isValid()){
	  for (std::size_t j = 0; j < muons.size(); j++) {
	    if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
	      skipjet = true;
	  }
	}
	if(electronsH.isValid()){
	  for (std::size_t j = 0; j < electrons.size(); j++) {
	    if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
	      skipjet = true;
	  }
	}
	if(photonsH.isValid()){
	  for (std::size_t j = 0; j < photons.size(); j++) {
	    if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
	      skipjet = true;
	  }
	}
      
	if (skipjet) continue;
      
	if (jets_iter->pt() > leadingjetpt) {
	  leadingjetpt  = jets_iter->pt() ;
	  leadingjeteta = jets_iter->eta();
	  leadingjetphi = jets_iter->phi();
	  leadingjetm   = jets_iter->mass();
	}
	            
	// apply jet id
	bool passjetid = applyJetID(*jets_iter,"loose");            
	if (!passjetid) 
	  continue;
	
	// apply pileup jet id
	bool passpuid = applyPileupJetID(*jets_iter,"medium",false);
	if (!passpuid) 
	  continue;
	
	pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());	
	if(jetref.isAvailable() and jetref.isNonnull())
	  incjets.push_back(jetref);
      }
    
      if(incjets.size() > 0)
	sort(incjets.begin(), incjets.end(), jetSorter);
      
      // only central jets
      for (size_t i = 0; i < incjets.size(); i++) {
	if (fabs(incjets[i]->eta()) <= 2.5) 
	  jets.push_back(incjets[i]);
      }        
      
      // sort them in pt
      if(jets.size() > 0)
	sort(jets.begin(), jets.end(), jetSorter);
      
      // count central jets    
      njets       = 0;
      nbjets      = 0;
      nbjetslowpt = 0;

      for (size_t i = 0; i < jets.size(); i++) {
	if (jets[i]->pt() > 30) njets++;
	if (jets[i]->pt() > 30 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjets++;
	if (jets[i]->pt() > 15 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjetslowpt++;
      }

      signaljetpt        = 0.0; signaljeteta       = 0.0; signaljetphi       = 0.0; signaljetbtag      = 0.0; signaljetCHfrac    = 0.0;
      signaljetNHfrac    = 0.0; signaljetEMfrac    = 0.0; signaljetCEMfrac   = 0.0; signaljetmetdphi   = 0.0;
      signaljetHFlav     = 0.0; signaljetPFlav     = 0.0; signaljetQGL       = 0.0; signaljetPUID      = 0.0;
      signaljetGenpt     = 0.0; signaljetGeneta    = 0.0; signaljetGenphi    = 0.0; signaljetGenm      = 0.0; signaljetRawpt = 0.0;
      signaljetm = 0.0;

      secondjetpt        = 0.0; secondjeteta       = 0.0; secondjetphi       = 0.0; secondjetbtag      = 0.0; secondjetCHfrac    = 0.0;
      secondjetNHfrac    = 0.0; secondjetEMfrac    = 0.0; secondjetCEMfrac   = 0.0; secondjetmetdphi   = 0.0;
      secondjetHFlav     = 0.0; secondjetPFlav     = 0.0; secondjetQGL       = 0.0; secondjetPUID      = 0.0;
      secondjetGenpt     = 0.0; secondjetGeneta    = 0.0; secondjetGenphi    = 0.0; secondjetGenm      = 0.0; secondjetRawpt = 0.0;
      secondjetm = 0.0;
      
      thirdjetpt        = 0.0; thirdjeteta       = 0.0; thirdjetphi       = 0.0; thirdjetbtag      = 0.0; thirdjetCHfrac    = 0.0;
      thirdjetNHfrac    = 0.0; thirdjetEMfrac    = 0.0; thirdjetCEMfrac   = 0.0; thirdjetmetdphi   = 0.0;
      thirdjetHFlav     = 0.0; thirdjetPFlav     = 0.0; thirdjetQGL       = 0.0; thirdjetPUID      = 0.0;
      thirdjetGenpt     = 0.0; thirdjetGeneta    = 0.0; thirdjetGenphi    = 0.0; thirdjetGenm      = 0.0; thirdjetRawpt = 0.0;
      thirdjetm = 0.0;
      
      fourthjetpt        = 0.0; fourthjeteta       = 0.0; fourthjetphi       = 0.0; fourthjetbtag      = 0.0; fourthjetCHfrac    = 0.0;
      fourthjetNHfrac    = 0.0; fourthjetEMfrac    = 0.0; fourthjetCEMfrac   = 0.0; fourthjetmetdphi   = 0.0;
      fourthjetHFlav     = 0.0; fourthjetPFlav     = 0.0; fourthjetQGL       = 0.0; fourthjetPUID      = 0.0;
      fourthjetGenpt     = 0.0; fourthjetGeneta    = 0.0; fourthjetGenphi    = 0.0; fourthjetGenm      = 0.0; fourthjetRawpt = 0.0;
      fourthjetm = 0.0;
      
      // only central jets  
      if (njets > 0) {
        signaljetpt      = jets[0]->pt();
        signaljeteta     = jets[0]->eta();
        signaljetphi     = jets[0]->phi();
        signaljetm       = jets[0]->mass();
	// b-tagging
        signaljetbtag    = jets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        signaljetCHfrac  = jets[0]->chargedHadronEnergyFraction();
        signaljetNHfrac  = jets[0]->neutralHadronEnergyFraction();
        signaljetEMfrac  = jets[0]->neutralEmEnergyFraction();
        signaljetCEMfrac = jets[0]->chargedEmEnergyFraction();
	signaljetRawpt   = jets[0]->correctedP4(0).pt();
	// Quark gluon likelihood
	if(jets[0]->hasUserFloat("QGTagger:qgLikelihood"))
	  signaljetQGL   = jets[0]->userFloat("QGTagger:qgLikelihood"); 
	// pileup jet id
	if(jets[0]->hasUserFloat("puid:fullDiscriminant"))
	  signaljetPUID  = jets[0]->userFloat("puid:fullDiscriminant");	
	else
	  signaljetPUID  = jets[0]->userFloat("pileupJetId:fullDiscriminant");
	// MC based info
	if(isMC){
	  signaljetHFlav = jets[0]->hadronFlavour(); 
	  signaljetPFlav = jets[0]->partonFlavour(); 
	  if(jets[0]->genJet()){
	    signaljetGenpt   = jets[0]->genJet()->pt(); 
	    signaljetGeneta  = jets[0]->genJet()->eta(); 
	    signaljetGenphi  = jets[0]->genJet()->phi(); 
	    signaljetGenm    = jets[0]->genJet()->mass(); 
	  }
	}
      }
      

      if (njets > 1) {
        secondjetpt      = jets[1]->pt();
        secondjeteta     = jets[1]->eta();
        secondjetphi     = jets[1]->phi();
        secondjetm       = jets[1]->mass();
	// b-tagging
        secondjetbtag    = jets[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        secondjetCHfrac  = jets[1]->chargedHadronEnergyFraction();
        secondjetNHfrac  = jets[1]->neutralHadronEnergyFraction();
        secondjetEMfrac  = jets[1]->neutralEmEnergyFraction();
        secondjetCEMfrac = jets[1]->chargedEmEnergyFraction();
	secondjetRawpt   = jets[1]->correctedP4(0).pt();
	// Quark gluon likelihood
	if(jets[1]->hasUserFloat("QGTagger:qgLikelihood"))
	  secondjetQGL   = jets[1]->userFloat("QGTagger:qgLikelihood"); 
	// pileup jet id
	if(jets[1]->hasUserFloat("puid:fullDiscriminant"))
	  secondjetPUID  = jets[1]->userFloat("puid:fullDiscriminant");
	else
	  secondjetPUID  = jets[1]->userFloat("pileupJetId:fullDiscriminant");
	// MC based info	
	if(isMC){
	  secondjetHFlav   = jets[1]->hadronFlavour(); 
	  secondjetPFlav   = jets[1]->partonFlavour(); 
	  if(jets[1]->genJet()){
	    secondjetGenpt     = jets[1]->genJet()->pt(); 
	    secondjetGeneta    = jets[1]->genJet()->eta(); 
	    secondjetGenphi    = jets[1]->genJet()->phi(); 
	    secondjetGenm      = jets[1]->genJet()->mass(); 
	  }
	}
      }
      
      if (njets > 2) {
        thirdjetpt      = jets[2]->pt();
        thirdjeteta     = jets[2]->eta();
        thirdjetphi     = jets[2]->phi();
        thirdjetm       = jets[2]->mass();
	// b-tagging
        thirdjetbtag    = jets[2]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        thirdjetCHfrac  = jets[2]->chargedHadronEnergyFraction();
        thirdjetNHfrac  = jets[2]->neutralHadronEnergyFraction();
        thirdjetEMfrac  = jets[2]->neutralEmEnergyFraction();
        thirdjetCEMfrac = jets[2]->chargedEmEnergyFraction();
	thirdjetRawpt   = jets[2]->correctedP4(0).pt();
	// QGL
	if(jets[2]->hasUserFloat("QGTagger:qgLikelihood"))
	  thirdjetQGL   = jets[2]->userFloat("QGTagger:qgLikelihood"); 
	// pileup jet id
	if(jets[2]->hasUserFloat("puid:fullDiscriminant"))
	  thirdjetPUID  = jets[2]->userFloat("puid:fullDiscriminant");
	else
	  thirdjetPUID  = jets[2]->userFloat("pileupJetId:fullDiscriminant");
	// MC based info
	if(isMC){
	  thirdjetHFlav   = jets[2]->hadronFlavour(); 
	  thirdjetPFlav   = jets[2]->partonFlavour(); 
	  if(jets[2]->genJet()){
	    thirdjetGenpt   = jets[2]->genJet()->pt(); 
	    thirdjetGeneta  = jets[2]->genJet()->eta(); 
	    thirdjetGenphi  = jets[2]->genJet()->phi(); 
	    thirdjetGenm    = jets[2]->genJet()->mass(); 
	  }
	}
      }

      if (njets > 3) {
        fourthjetpt      = jets[3]->pt();
        fourthjeteta     = jets[3]->eta();
        fourthjetphi     = jets[3]->phi();
        fourthjetm       = jets[3]->mass();
        fourthjetbtag    = jets[3]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        fourthjetCHfrac  = jets[3]->chargedHadronEnergyFraction();
        fourthjetNHfrac  = jets[3]->neutralHadronEnergyFraction();
        fourthjetEMfrac  = jets[3]->neutralEmEnergyFraction();
        fourthjetCEMfrac = jets[3]->chargedEmEnergyFraction();
	fourthjetRawpt   = jets[3]->correctedP4(0).pt();

	if(jets[3]->hasUserFloat("QGTagger:qgLikelihood"))
	  fourthjetQGL   = jets[3]->userFloat("QGTagger:qgLikelihood"); 
	if(jets[3]->hasUserFloat("puid:fullDiscriminant"))
	  fourthjetPUID  = jets[3]->userFloat("puid:fullDiscriminant");
	else
	  fourthjetPUID  = jets[3]->userFloat("pileupJetId:fullDiscriminant");

	if(isMC){
	  fourthjetHFlav   = jets[3]->hadronFlavour(); 
	  fourthjetPFlav   = jets[3]->partonFlavour(); 
	  
	  if(jets[3]->genJet()){
	    fourthjetGenpt   = jets[3]->genJet()->pt(); 
	    fourthjetGeneta  = jets[3]->genJet()->eta(); 
	    fourthjetGenphi  = jets[3]->genJet()->phi(); 
	    fourthjetGenm    = jets[3]->genJet()->mass(); 
	  }
	}
       }
      
      
      jetjetdphi         = 0.0;   
      signaljetmetdphi   = 0.0;   secondjetmetdphi = 0.0; thirdjetmetdphi = 0.0; fourthjetmetdphi = 0.0; 
      jetmetdphimin      = 0.0;   incjetmetdphimin    = 0.0;
      jetmumetdphimin    = 0.0;   incjetmumetdphimin  = 0.0; 
      jetelmetdphimin    = 0.0;   incjetelmetdphimin  = 0.0;
      jetphmetdphimin    = 0.0;   incjetphmetdphimin  = 0.0;
      jetmetdphimin4     = 0.0;   incjetmetdphimin4   = 0.0;
      jetmumetdphimin4   = 0.0;   incjetmumetdphimin4 = 0.0;
      jetelmetdphimin4   = 0.0;   incjetelmetdphimin4 = 0.0;
      jetphmetdphimin4   = 0.0;   incjetphmetdphimin4 = 0.0;
      
      // delta phi between jets
      if (signaljetpt > 0.0 && secondjetpt > 0.0) 
	jetjetdphi = deltaPhi(signaljetphi, secondjetphi);
      // delta phi jet-met
      if (signaljetpt > 0.0) 
	signaljetmetdphi = deltaPhi(signaljetphi, t1pfmetphi);
      if (secondjetpt > 0.0) 
	secondjetmetdphi = deltaPhi(secondjetphi, t1pfmetphi);
      if (thirdjetpt  > 0.0) 
	thirdjetmetdphi  = deltaPhi(thirdjetphi , t1pfmetphi);
      if (fourthjetpt > 0.0) 
	fourthjetmetdphi = deltaPhi(fourthjetphi, t1pfmetphi);
      
      std::vector<double> jetmetdphiminvector;
      std::vector<double> jetmetdphimin4vector;
      // only central jets
      for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > 30) {
	  double jetphi = jets[i]->phi();
	  jetmetdphiminvector.push_back(fabs(deltaPhi(jetphi, t1pfmetphi)));
	  if (i < 4) jetmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, t1pfmetphi)));
        }
      }
      if (jetmetdphiminvector .size() > 0) jetmetdphimin  = *min_element(jetmetdphiminvector .begin(), jetmetdphiminvector .end());
      if (jetmetdphimin4vector.size() > 0) jetmetdphimin4 = *min_element(jetmetdphimin4vector.begin(), jetmetdphimin4vector.end());
      
      std::vector<double> incjetmetdphiminvector;
      std::vector<double> incjetmetdphimin4vector;
      // central and forward jet
      for (size_t i = 0; i < incjets.size(); i++) {
	if (incjets[i]->pt() > 30) {
	  double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
	  incjetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1pfmetphi)));
	  if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1pfmetphi)));
        }
      }
      if (incjetmetdphiminvector .size() > 0) incjetmetdphimin  = *min_element(incjetmetdphiminvector .begin(), incjetmetdphiminvector .end());
      if (incjetmetdphimin4vector.size() > 0) incjetmetdphimin4 = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      
      // MU met
      std::vector<double> jetmumetdphiminvector;
      std::vector<double> jetmumetdphimin4vector;
      for (size_t i = 0; i < jets.size(); i++) {
	if (jets[i]->pt() > 30) {
	  double jetphi = jets[i]->phi();
	  jetmumetdphiminvector.push_back(fabs(deltaPhi(jetphi, t1mumetphi)));
	  if (i < 4) jetmumetdphimin4vector.push_back(fabs(deltaPhi(jetphi, t1mumetphi)));
        }
    }
      if (jetmumetdphiminvector .size() > 0) jetmumetdphimin  = *min_element(jetmumetdphiminvector .begin(), jetmumetdphiminvector .end());
      if (jetmumetdphimin4vector.size() > 0) jetmumetdphimin4 = *min_element(jetmumetdphimin4vector.begin(), jetmumetdphimin4vector.end());
      
      std::vector<double> incjetmumetdphiminvector;
      std::vector<double> incjetmumetdphimin4vector;
      
      for (size_t i = 0; i < incjets.size(); i++) {
        if (incjets[i]->pt() > 30) {
	  double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
            incjetmumetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
            if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
        }
      }
      if (incjetmumetdphiminvector .size() > 0) incjetmumetdphimin  = *min_element(incjetmumetdphiminvector .begin(), incjetmumetdphiminvector .end());
      if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4 = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      
    // ELE met
      std::vector<double> jetelmetdphiminvector;
      std::vector<double> jetelmetdphimin4vector;
      for (size_t i = 0; i < jets.size(); i++) {
	if (jets[i]->pt() > 30) {
            double jetphi = jets[i]->phi();
            jetelmetdphiminvector.push_back(fabs(deltaPhi(jetphi, t1elmetphi)));
            if (i < 4) jetelmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, t1elmetphi)));
        }
      }
      if (jetelmetdphiminvector .size() > 0) jetelmetdphimin  = *min_element(jetelmetdphiminvector .begin(), jetelmetdphiminvector .end());
      if (jetelmetdphimin4vector.size() > 0) jetelmetdphimin4 = *min_element(jetelmetdphimin4vector.begin(), jetelmetdphimin4vector.end());
      
      std::vector<double> incjetelmetdphiminvector;
      std::vector<double> incjetelmetdphimin4vector;
      for (size_t i = 0; i < incjets.size(); i++) {
        if (incjets[i]->pt() > 30) {
            double incjetphi = incjets[i]->phi();
            incjetelmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1elmetphi)));
            if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetphi)));
        }
      }
      if (incjetelmetdphiminvector .size() > 0) incjetelmetdphimin  = *min_element(incjetelmetdphiminvector .begin(), incjetelmetdphiminvector .end());
      if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4 = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      
      // PH met
      std::vector<double> jetphmetdphiminvector;
      std::vector<double> jetphmetdphimin4vector;
      for (size_t i = 0; i < jets.size(); i++) {
	if (jets[i]->pt() > 30) {
	  double jetphi = jets[i]->phi();
	  jetphmetdphiminvector.push_back(fabs(deltaPhi(jetphi, t1phmetphi)));
	  if (i < 4) jetphmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, t1phmetphi)));
        }
      }
      if (jetphmetdphiminvector.size()  > 0) jetphmetdphimin  = *min_element(jetphmetdphiminvector .begin(), jetphmetdphiminvector .end());
      if (jetphmetdphimin4vector.size() > 0) jetphmetdphimin4 = *min_element(jetphmetdphimin4vector.begin(), jetphmetdphimin4vector.end());
      
      std::vector<double> incjetphmetdphiminvector;
      std::vector<double> incjetphmetdphimin4vector;
      for (size_t i = 0; i < incjets.size(); i++) {
        if (incjets[i]->pt() > 30) {
	  double incjetphi = incjets[i]->phi();
	  incjetphmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1phmetphi)));
	  if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetphi)));
        }
      }
      if (incjetphmetdphiminvector .size() > 0) incjetphmetdphimin  = *min_element(incjetphmetdphiminvector .begin(), incjetphmetdphiminvector .end());
      if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4 = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());

      // QCD suppression handles
      ht     = 0.;
      std::vector<double> jetEts;
      for (size_t i = 0; i < jets.size(); i++) {
	if (jets[i]->pt() > 30) {
	  ht += jets[i]->pt();
	  jetEts.push_back(jets[i]->pt());
	}
      }
    }
    
    
    // PUPPI AK4 jets
    if(addPuppiJets){

      leadingPuppijetpt  = 0.0;
      leadingPuppijeteta = 0.0;
      leadingPuppijetphi = 0.0;
      leadingPuppijetm   = 0.0;

      vector<pat::JetRef> incPuppijets;
      vector<pat::JetRef> Puppijets;

      if(jetsPuppiH.isValid()){

	for (auto jets_iter = jetsPuppiH->begin(); jets_iter != jetsPuppiH->end(); ++jets_iter) {
	  
	  //clean from leptons
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
		skipjet = true;
	    }
	  }
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) 
		skipjet = true;
	    }
	  }
	  
	  if (skipjet) continue;
	  
	  if (jets_iter->pt() > leadingPuppijetpt) {
	    leadingPuppijetpt  = jets_iter->pt() ;
	    leadingPuppijeteta = jets_iter->eta();
	    leadingPuppijetphi = jets_iter->phi();
	    leadingPuppijetm   = jets_iter->mass();
	  }
	  
	  // apply jet id
	  bool passjetid = applyJetID(*jets_iter,"loose");            
	  if (!passjetid) 
	    continue;
	  
	  //apply pileup jet id
	  bool passpuid = applyPileupJetID(*jets_iter,"medium",true);
	  if (!passpuid) continue;
	  
	  pat::JetRef jetref(jetsPuppiH, jets_iter - jetsPuppiH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    incPuppijets.push_back(jetref);
	}
	
	// only central jets
	if(incPuppijets.size() > 0)
	  sort(incPuppijets.begin(), incPuppijets.end(), jetSorter);

	for (size_t i = 0; i < incPuppijets.size(); i++) {
	  if (fabs(incPuppijets[i]->eta()) <= 2.5) 
	    Puppijets.push_back(incPuppijets[i]);
	}        
	
	// sort them in pt
	if(Puppijets.size() > 0)
	  sort(Puppijets.begin(), Puppijets.end(), jetSorter);
	
	// count central jets    
	npuppijets       = 0;
	npuppibjets      = 0;
	npuppibjetslowpt = 0;

	for (size_t i = 0; i < Puppijets.size(); i++) {
	  if (Puppijets[i]->pt() > 30) npuppijets++;
	  if (Puppijets[i]->pt() > 30 && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) npuppibjets++;
	  if (Puppijets[i]->pt() > 15 && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) npuppibjetslowpt++;
	}
	
	signalPuppijetpt        = 0.0; signalPuppijeteta       = 0.0; signalPuppijetphi       = 0.0; signalPuppijetbtag      = 0.0; signalPuppijetCHfrac    = 0.0;
	signalPuppijetNHfrac    = 0.0; signalPuppijetEMfrac    = 0.0; signalPuppijetCEMfrac   = 0.0; signalPuppijetmetdphi   = 0.0;
	signalPuppijetHFlav     = 0.0; signalPuppijetPFlav     = 0.0; signalPuppijetQGL       = 0.0; signalPuppijetPUID      = 0.0;
	signalPuppijetGenpt     = 0.0; signalPuppijetGeneta    = 0.0; signalPuppijetGenphi    = 0.0; signalPuppijetGenm      = 0.0; signalPuppijetRawpt = 0.0;

	secondPuppijetpt        = 0.0; secondPuppijeteta       = 0.0; secondPuppijetphi       = 0.0; secondPuppijetbtag      = 0.0; secondPuppijetCHfrac    = 0.0;
	secondPuppijetNHfrac    = 0.0; secondPuppijetEMfrac    = 0.0; secondPuppijetCEMfrac   = 0.0; secondPuppijetmetdphi   = 0.0;
	secondPuppijetHFlav     = 0.0; secondPuppijetPFlav     = 0.0; secondPuppijetQGL       = 0.0; secondPuppijetPUID      = 0.0;
	secondPuppijetGenpt     = 0.0; secondPuppijetGeneta    = 0.0; secondPuppijetGenphi    = 0.0; secondPuppijetGenm      = 0.0; secondPuppijetRawpt = 0.0;

	thirdPuppijetpt        = 0.0; thirdPuppijeteta       = 0.0; thirdPuppijetphi       = 0.0; thirdPuppijetbtag      = 0.0; thirdPuppijetCHfrac    = 0.0;
	thirdPuppijetNHfrac    = 0.0; thirdPuppijetEMfrac    = 0.0; thirdPuppijetCEMfrac   = 0.0; thirdPuppijetmetdphi   = 0.0;
	thirdPuppijetHFlav     = 0.0; thirdPuppijetPFlav     = 0.0; thirdPuppijetQGL       = 0.0; thirdPuppijetPUID      = 0.0;
	thirdPuppijetGenpt     = 0.0; thirdPuppijetGeneta    = 0.0; thirdPuppijetGenphi    = 0.0; thirdPuppijetGenm      = 0.0; thirdPuppijetRawpt = 0.0;
	
	fourthPuppijetpt        = 0.0; fourthPuppijeteta       = 0.0; fourthPuppijetphi       = 0.0; fourthPuppijetbtag      = 0.0; fourthPuppijetCHfrac    = 0.0;
	fourthPuppijetNHfrac    = 0.0; fourthPuppijetEMfrac    = 0.0; fourthPuppijetCEMfrac   = 0.0; fourthPuppijetmetdphi   = 0.0;
	fourthPuppijetHFlav     = 0.0; fourthPuppijetPFlav     = 0.0; fourthPuppijetQGL       = 0.0; fourthPuppijetPUID      = 0.0;
	fourthPuppijetGenpt     = 0.0; fourthPuppijetGeneta    = 0.0; fourthPuppijetGenphi    = 0.0; fourthPuppijetGenm      = 0.0; fourthPuppijetRawpt = 0.0;
	
	// only central jets  
	if (npuppijets > 0) {
	  signalPuppijetpt      = Puppijets[0]->pt();
	  signalPuppijeteta     = Puppijets[0]->eta();
	  signalPuppijetphi     = Puppijets[0]->phi();
	  signalPuppijetm       = Puppijets[0]->mass();
	  signalPuppijetbtag    = Puppijets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	  signalPuppijetCHfrac  = Puppijets[0]->chargedHadronEnergyFraction();
	  signalPuppijetNHfrac  = Puppijets[0]->neutralHadronEnergyFraction();
	  signalPuppijetEMfrac  = Puppijets[0]->neutralEmEnergyFraction();
	  signalPuppijetCEMfrac = Puppijets[0]->chargedEmEnergyFraction();
	  signalPuppijetRawpt   = Puppijets[0]->correctedP4(0).pt();
	  
	  if(Puppijets[0]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
	    signalPuppijetQGL     = Puppijets[0]->userFloat("QGTaggerPuppi:qgLikelihood"); 
	  
	  if(Puppijets[0]->hasUserFloat("puidPuppi:fullDiscriminant"))
	    signalPuppijetPUID    = Puppijets[0]->userFloat("puidPuppi:fullDiscriminant");
	  else if (Puppijets[0]->hasUserFloat("pileupJetIdPuppi:fullDiscriminant"))
	    signalPuppijetPUID    = Puppijets[0]->userFloat("pileupJetId:fullDiscriminant");
	  
	  if(isMC){
	    signalPuppijetHFlav   = Puppijets[0]->hadronFlavour(); 
	    signalPuppijetPFlav   = Puppijets[0]->partonFlavour(); 
	    if(Puppijets[0]->genJet()){
	      signalPuppijetGenpt     = Puppijets[0]->genJet()->pt(); 
	      signalPuppijetGeneta    = Puppijets[0]->genJet()->eta(); 
	      signalPuppijetGenphi    = Puppijets[0]->genJet()->phi(); 
	      signalPuppijetGenm      = Puppijets[0]->genJet()->mass();
	    }
	  } 
	}
      	
	if (npuppijets > 1) {
	  secondPuppijetpt      = Puppijets[1]->pt();
	  secondPuppijeteta     = Puppijets[1]->eta();
	  secondPuppijetphi     = Puppijets[1]->phi();
	  secondPuppijetm       = Puppijets[1]->mass();
	  secondPuppijetbtag    = Puppijets[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	  secondPuppijetCHfrac  = Puppijets[1]->chargedHadronEnergyFraction();
	  secondPuppijetNHfrac  = Puppijets[1]->neutralHadronEnergyFraction();
	  secondPuppijetEMfrac  = Puppijets[1]->neutralEmEnergyFraction();
	  secondPuppijetCEMfrac = Puppijets[1]->chargedEmEnergyFraction();
	  secondPuppijetRawpt   = Puppijets[1]->correctedP4(0).pt();
	  
	  if(Puppijets[1]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
	    secondPuppijetQGL     = Puppijets[1]->userFloat("QGTaggerPuppi:qgLikelihood"); 
	  
	  if(Puppijets[1]->hasUserFloat("puidPuppi:fullDiscriminant"))
	    secondPuppijetPUID    = Puppijets[1]->userFloat("puidPuppi:fullDiscriminant");
	  else if(Puppijets[1]->hasUserFloat("pileupJetIdPuppi:fullDiscriminant"))
	    secondPuppijetPUID    = Puppijets[1]->userFloat("pileupJetIdPuppi:fullDiscriminant");
	  
	  if(isMC){
	    secondPuppijetHFlav   = Puppijets[1]->hadronFlavour(); 
	    secondPuppijetPFlav   = Puppijets[1]->partonFlavour(); 
	    if(Puppijets[1]->genJet()){
	      secondPuppijetGenpt     = Puppijets[1]->genJet()->pt(); 
	      secondPuppijetGeneta    = Puppijets[1]->genJet()->eta(); 
	      secondPuppijetGenphi    = Puppijets[1]->genJet()->phi(); 
	      secondPuppijetGenm      = Puppijets[1]->genJet()->mass(); 
	    }
	  }
	}
      
	if (npuppijets > 2) {
	  thirdPuppijetpt      = Puppijets[2]->pt();
	  thirdPuppijeteta     = Puppijets[2]->eta();
	  thirdPuppijetphi     = Puppijets[2]->phi();
	  thirdPuppijetm       = Puppijets[2]->mass();
	  thirdPuppijetbtag    = Puppijets[2]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	  thirdPuppijetCHfrac  = Puppijets[2]->chargedHadronEnergyFraction();
	  thirdPuppijetNHfrac  = Puppijets[2]->neutralHadronEnergyFraction();
	  thirdPuppijetEMfrac  = Puppijets[2]->neutralEmEnergyFraction();
	  thirdPuppijetCEMfrac = Puppijets[2]->chargedEmEnergyFraction();
	  thirdPuppijetRawpt   = Puppijets[2]->correctedP4(0).pt();

	  if(Puppijets[2]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
	    thirdPuppijetQGL     = Puppijets[2]->userFloat("QGTaggerPuppi:qgLikelihood"); 
	  
	  if(Puppijets[2]->hasUserFloat("puidPuppi:fullDiscriminant"))
	    thirdPuppijetPUID    = Puppijets[2]->userFloat("puidPuppi:fullDiscriminant");
	  else if(Puppijets[2]->hasUserFloat("pileupJetIdPuppi:fullDiscriminant"))
	    thirdPuppijetPUID    = Puppijets[2]->userFloat("pileupJetIdPuppi:fullDiscriminant");

	  if(isMC){
	    thirdPuppijetHFlav   = Puppijets[2]->hadronFlavour(); 
	    thirdPuppijetPFlav   = Puppijets[2]->partonFlavour(); 
	    if(Puppijets[2]->genJet()){
	      thirdPuppijetGenpt     = Puppijets[2]->genJet()->pt(); 
	      thirdPuppijetGeneta    = Puppijets[2]->genJet()->eta(); 
	      thirdPuppijetGenphi    = Puppijets[2]->genJet()->phi(); 
	      thirdPuppijetGenm      = Puppijets[2]->genJet()->mass();
	    }
	  } 
	}
	
	if (npuppijets > 3) {
	  fourthPuppijetpt      = Puppijets[3]->pt();
	  fourthPuppijeteta     = Puppijets[3]->eta();
	  fourthPuppijetphi     = Puppijets[3]->phi();
	  fourthPuppijetm       = Puppijets[3]->mass();
	  fourthPuppijetbtag    = Puppijets[3]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	  fourthPuppijetCHfrac  = Puppijets[3]->chargedHadronEnergyFraction();
	  fourthPuppijetNHfrac  = Puppijets[3]->neutralHadronEnergyFraction();
	  fourthPuppijetEMfrac  = Puppijets[3]->neutralEmEnergyFraction();
	  fourthPuppijetCEMfrac = Puppijets[3]->chargedEmEnergyFraction();
	  fourthPuppijetRawpt   = Puppijets[3]->correctedP4(0).pt();
	  
	  if(Puppijets[3]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
	    fourthPuppijetQGL     = Puppijets[3]->userFloat("QGTaggerPuppi:qgLikelihood"); 
	  
	  if(Puppijets[3]->hasUserFloat("puidPuppi:fullDiscriminant"))
	    fourthPuppijetPUID    = Puppijets[3]->userFloat("puidPuppi:fullDiscriminant");
	  else if(Puppijets[3]->hasUserFloat("pileupJetIdPuppi:fullDiscriminant"))
	    fourthPuppijetPUID    = Puppijets[3]->userFloat("pileupJetIdPuppi:fullDiscriminant");
	  
	  if(isMC){
	    if(Puppijets[3]->genJet()){
	      fourthPuppijetGenpt     = Puppijets[3]->genJet()->pt(); 
	      fourthPuppijetGeneta    = Puppijets[3]->genJet()->eta(); 
	      fourthPuppijetGenphi    = Puppijets[3]->genJet()->phi(); 
	      fourthPuppijetGenm      = Puppijets[3]->genJet()->mass(); 
	    }
	    fourthPuppijetHFlav   = Puppijets[3]->hadronFlavour(); 
	    fourthPuppijetPFlav   = Puppijets[3]->partonFlavour(); 
	  }
	}
	
	PuppijetPuppijetdphi    = 0.0;
	signalPuppijetmetdphi   = 0.0; secondPuppijetmetdphi = 0.0; thirdPuppijetmetdphi = 0.0; fourthPuppijetmetdphi = 0.0; 
	Puppijetmetdphimin      = 0.0; incPuppijetmetdphimin   = 0.0;
	Puppijetmumetdphimin    = 0.0; incPuppijetmumetdphimin = 0.0;
	Puppijetelmetdphimin    = 0.0; incPuppijetelmetdphimin = 0.0;
	Puppijetphmetdphimin    = 0.0; incPuppijetphmetdphimin = 0.0;
	
	Puppijetmetdphimin4     = 0.0; incPuppijetmetdphimin4  = 0.0;
	Puppijetmumetdphimin4   = 0.0; incPuppijetmumetdphimin4= 0.0;
	Puppijetelmetdphimin4   = 0.0; incPuppijetelmetdphimin4= 0.0;
	Puppijetphmetdphimin4   = 0.0; incPuppijetphmetdphimin4= 0.0;
	
	// delta phi between Puppijets
	if (signalPuppijetpt > 0.0 && secondPuppijetpt > 0.0) 
	  PuppijetPuppijetdphi = deltaPhi(signalPuppijetphi, secondPuppijetphi);
	// delta phi jet-met
	if (signalPuppijetpt > 0.0) 
	  signalPuppijetmetdphi = deltaPhi(signalPuppijetphi, puppit1pfmetphi);
	if (secondPuppijetpt > 0.0) 
	  secondPuppijetmetdphi = deltaPhi(secondPuppijetphi, puppit1pfmetphi);
	if (thirdPuppijetpt  > 0.0) 
	  thirdPuppijetmetdphi  = deltaPhi(thirdPuppijetphi , puppit1pfmetphi);
	if (fourthPuppijetpt > 0.0) 
	  fourthPuppijetmetdphi = deltaPhi(fourthPuppijetphi, puppit1pfmetphi);
	
	std::vector<double> puppijetmetdphiminvector;
	std::vector<double> puppijetmetdphimin4vector;
	// only central Puppijets
	for (size_t i = 0; i < Puppijets.size(); i++) {
	  if (Puppijets[i]->pt() > 30) {
	    double jetphi = Puppijets[i]->phi();
	    puppijetmetdphiminvector.push_back(fabs(deltaPhi(jetphi, puppit1pfmetphi)));
	    if (i < 4) puppijetmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, puppit1pfmetphi)));
	  }
	}
	if (puppijetmetdphiminvector .size() > 0) Puppijetmetdphimin  = *min_element(puppijetmetdphiminvector .begin(), puppijetmetdphiminvector .end());
	if (puppijetmetdphimin4vector.size() > 0) Puppijetmetdphimin4 = *min_element(puppijetmetdphimin4vector.begin(), puppijetmetdphimin4vector.end());
	
	std::vector<double> incpuppijetmetdphiminvector;
	std::vector<double> incpuppijetmetdphimin4vector;
	// central and forward jet
	for (size_t i = 0; i < incPuppijets.size(); i++) {
	  if (incPuppijets[i]->pt() > 30) {
	    double incjetphi = atan2(sin(incPuppijets[i]->phi()), cos(incPuppijets[i]->phi()));
	    incpuppijetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1pfmetphi)));
	    if (i < 4) incpuppijetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1pfmetphi)));
	  }
	}
	if (incpuppijetmetdphiminvector .size() > 0) incPuppijetmetdphimin  = *min_element(incpuppijetmetdphiminvector .begin(), incpuppijetmetdphiminvector .end());
	if (incpuppijetmetdphimin4vector.size() > 0) incPuppijetmetdphimin4 = *min_element(incpuppijetmetdphimin4vector.begin(), incpuppijetmetdphimin4vector.end());
	
	// MU met
	std::vector<double> puppijetmumetdphiminvector;
	std::vector<double> puppijetmumetdphimin4vector;
	for (size_t i = 0; i < jets.size(); i++) {
	  if (jets[i]->pt() > 30) {
	    double jetphi = jets[i]->phi();
	    puppijetmumetdphiminvector.push_back(fabs(deltaPhi(jetphi, puppit1mumetphi)));
	    if (i < 4) puppijetmumetdphimin4vector.push_back(fabs(deltaPhi(jetphi, puppit1mumetphi)));
	  }
	}
	if (puppijetmumetdphiminvector .size() > 0) Puppijetmumetdphimin  = *min_element(puppijetmumetdphiminvector .begin(), puppijetmumetdphiminvector .end());
	if (puppijetmumetdphimin4vector.size() > 0) Puppijetmumetdphimin4 = *min_element(puppijetmumetdphimin4vector.begin(), puppijetmumetdphimin4vector.end());
	
	std::vector<double> incpuppijetmumetdphiminvector;
	std::vector<double> incpuppijetmumetdphimin4vector;
	
	for (size_t i = 0; i < incjets.size(); i++) {
	  if (incjets[i]->pt() > 30) {
	    double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
	    incpuppijetmumetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1mumetphi)));
	    if (i < 4) incpuppijetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1mumetphi)));
        }
	}
	if (incpuppijetmumetdphiminvector .size() > 0) incPuppijetmumetdphimin  = *min_element(incpuppijetmumetdphiminvector .begin(), incpuppijetmumetdphiminvector .end());
	if (incpuppijetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4 = *min_element(incpuppijetmumetdphimin4vector.begin(), incpuppijetmumetdphimin4vector.end());
	
	// ELE met
	std::vector<double> puppijetelmetdphiminvector;
	std::vector<double> puppijetelmetdphimin4vector;
	for (size_t i = 0; i < jets.size(); i++) {
	  if (jets[i]->pt() > 30) {
	    double jetphi = jets[i]->phi();
	    puppijetelmetdphiminvector.push_back(fabs(deltaPhi(jetphi, puppit1elmetphi)));
	    if (i < 4) puppijetelmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, puppit1elmetphi)));
	  }
	}
	if (puppijetelmetdphiminvector .size() > 0) Puppijetelmetdphimin  = *min_element(puppijetelmetdphiminvector .begin(), puppijetelmetdphiminvector .end());
	if (puppijetelmetdphimin4vector.size() > 0) Puppijetelmetdphimin4 = *min_element(puppijetelmetdphimin4vector.begin(), puppijetelmetdphimin4vector.end());
	
	std::vector<double> incpuppijetelmetdphiminvector;
	std::vector<double> incpuppijetelmetdphimin4vector;
	for (size_t i = 0; i < incjets.size(); i++) {
	  if (incjets[i]->pt() > 30) {
	    double incjetphi = incjets[i]->phi();
	    incpuppijetelmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, puppit1elmetphi)));
	    if (i < 4) incpuppijetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, puppit1elmetphi)));
	  }
	}
	if (incpuppijetelmetdphiminvector .size() > 0) incPuppijetelmetdphimin  = *min_element(incpuppijetelmetdphiminvector .begin(), incpuppijetelmetdphiminvector .end());
	if (incpuppijetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4 = *min_element(incpuppijetelmetdphimin4vector.begin(), incpuppijetelmetdphimin4vector.end());
	
	// PH met
	std::vector<double> puppijetphmetdphiminvector;
	std::vector<double> puppijetphmetdphimin4vector;
	for (size_t i = 0; i < jets.size(); i++) {
	  if (jets[i]->pt() > 30) {
	    double puppijetphi = jets[i]->phi();
	    puppijetphmetdphiminvector.push_back(fabs(deltaPhi(puppijetphi, puppit1phmetphi)));
	    if (i < 4) puppijetphmetdphimin4vector.push_back(fabs(deltaPhi(puppijetphi, puppit1phmetphi)));
        }
	}
	if (puppijetphmetdphiminvector.size()  > 0) Puppijetphmetdphimin  = *min_element(puppijetphmetdphiminvector .begin(), puppijetphmetdphiminvector .end());
	if (puppijetphmetdphimin4vector.size() > 0) Puppijetphmetdphimin4 = *min_element(puppijetphmetdphimin4vector.begin(), puppijetphmetdphimin4vector.end());
	
	std::vector<double> incpuppijetphmetdphiminvector;
	std::vector<double> incpuppijetphmetdphimin4vector;
	for (size_t i = 0; i < incjets.size(); i++) {
	  if (incjets[i]->pt() > 30) {
	    double incpuppijetphi = incjets[i]->phi();
	    incpuppijetphmetdphiminvector.push_back(fabs(deltaPhi(incpuppijetphi, puppit1phmetphi)));
	    if (i < 4) incpuppijetphmetdphimin4vector.push_back(fabs(deltaPhi(incpuppijetphi, puppit1phmetphi)));
	  }
	}
	if (incpuppijetphmetdphiminvector .size() > 0) incPuppijetphmetdphimin  = *min_element(incpuppijetphmetdphiminvector .begin(), incpuppijetphmetdphiminvector .end());
	if (incpuppijetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4 = *min_element(incpuppijetphmetdphimin4vector.begin(), incpuppijetphmetdphimin4vector.end());
	
	// QCD suppression handles
	Puppiht     = 0.;
	std::vector<double> PuppijetEts;
	for (size_t i = 0; i < Puppijets.size(); i++) {
	  if (Puppijets[i]->pt() > 30) {
	    Puppiht += Puppijets[i]->pt();
	    PuppijetEts.push_back(Puppijets[i]->pt());
	  }
	}        
      }
    }
    
    // Lepton counts
    ntaus  = 0;
    if(tausH.isValid()){
      for (auto taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
	bool skiptau = false;
	for (std::size_t j = 0; j < muons.size(); j++) {
	  if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
	}
	for (std::size_t j = 0; j < electrons.size(); j++) {
	  if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
	}
      if (taus_iter->pt() > 18 && 
	  fabs(taus_iter->eta()) < 2.3 && 
	  taus_iter->tauID("decayModeFinding") > 0.5 && 
	  taus_iter->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 5 && !skiptau) ntaus++;
      }
    }
    
    vector<pat::MuonRef> muonvector;
    if(muonsH.isValid() and tightmuonsH.isValid() and highptmuonsH.isValid()){
      nmuons          = muonsH->size();
      ntightmuons     = tightmuonsH->size();
      nhighptmuons    = highptmuonsH->size();
      
      for (size_t i = 0; i < muons.size(); i++) 
	muonvector.push_back(muons[i]);
    }

    vector<pat::ElectronRef> electronvector;
    if(electronsH.isValid() and tightelectronsH.isValid() and heepelectronsH.isValid()){
      
      nelectrons      = electronsH->size();
      ntightelectrons = tightelectronsH->size();
      nheepelectrons  = heepelectronsH->size();

      for (size_t i = 0; i < electrons.size(); i++) 
	electronvector.push_back(electrons[i]);
    }

    // W, Z control sample information
    zmass       = 0.0; zpt         = 0.0; zeta        = 0.0; zphi        = 0.0;
    zeemass     = 0.0; zeept       = 0.0; zeeeta      = 0.0; zeephi      = 0.0;
    wmt         = 0.0; wemt        = 0.0; emumass     = 0.0; emupt       = 0.0;
    emueta      = 0.0; emuphi      = 0.0;

    mu1pid      = 0;   mu1pt       = 0.0; mu1eta      = 0.0; mu1phi      = 0.0;
    mu1pfpt     = 0.0; mu1pfeta    = 0.0; mu1pfphi    = 0.0; mu1id       = 0;
    mu1idm      = 0;   mu1idt      = 0;   mu1iso      = 0.0;

    mu2pid      = 0;   mu2pt       = 0.0; mu2eta      = 0.0; mu2phi      = 0.0;
    mu2pfpt     = 0.0; mu2pfeta    = 0.0; mu2pfphi    = 0.0; mu2id       = 0;
    mu2idm      = 0;   mu2idt      = 0;   mu2iso      = 0.0;

    el1pid      = 0; el1pt       = 0.0; el1eta      = 0.0; el1phi      = 0.0; el1id       = 0;
    el2pid      = 0; el2pt       = 0.0; el2eta      = 0.0; el2phi      = 0.0; el2id       = 0;

    // sort electrons and muons
    sort(muonvector.begin(), muonvector.end(), muonSorter);
    sort(electronvector.begin(), electronvector.end(), electronSorter);
    
    // one or two loose muons
    if (nmuons == 1 || nmuons == 2) {
      
      pat::MuonRef muon = muons[0];
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
	  mu1id = 1;
      }
      
      // store high-pt muons that are not tight ones
      for (std::size_t i = 0; i < highptmuons.size(); i++) {
	if (muon == highptmuons[i] and mu1id != 1) 
	  mu1id = 2;
      }

      if (nmuons == 1) 
	wmt = sqrt(2.0 * mu1pt * t1pfmet * (1.0 - cos(deltaPhi(mu1phi, t1pfmetphi))));
    }
   
    // two loose muons
    if (nmuons == 2) {        

      pat::MuonRef muon = muons[1];
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
      mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
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
      pat::ElectronRef electron = electrons[0];
      el1pid = electron->pdgId();
      el1pt  = electron->pt();
      el1eta = electron->eta();
      el1phi = electron->phi();
      el1idl = ((*electronLooseIdH )[electron] ? 1 : 0);
      
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
        pat::ElectronRef electron = electrons[1];
        el2pid = electron->pdgId();
        el2pt  = electron->pt();
        el2eta = electron->eta();
        el2phi = electron->phi();
        el2idl = ((*electronLooseIdH )[electron] ? 1 : 0);

        for (std::size_t i = 0; i < tightelectrons.size(); i++) {
            if (electron == tightelectrons[i]) el2id = 1;
        }

	for (std::size_t i = 0; i < heepelectrons.size(); i++) {
	  if (electron == heepelectrons[i] and el2id != 1) 
	    el2id = 2;
	}
 
        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
        TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electron->p());

        TLorentzVector zvec(el1vec);
        zvec += el2vec;

        zeemass = zvec.M();
        zeept   = zvec.Pt();
        zeeeta  = zvec.Eta();
        zeephi  = zvec.Phi();
    }

    // one electron and one muon (ttbar DF fully leptonic)
    if (nmuons == 1 && nelectrons == 1) {
      TLorentzVector mu1vec; 
      mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
      TLorentzVector el1vec; 
      el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
      
      TLorentzVector emuvec(mu1vec);
      emuvec += el1vec;
      
      emumass = emuvec.M();
      emupt   = emuvec.Pt();
      emueta  = emuvec.Eta();
      emuphi  = emuvec.Phi();
    } 
    
    // Photon information
    phidl    = 0; phidm    = 0; phidt    = 0; phidh    = 0;
    phpt     = 0.0; pheta    = 0.0; phphi    = 0.0;

    int hardestPhotonIndex = -1;
    double hardestPhotonPt = 0.0;
    if(photonsH.isValid() and photonLooseIdH.isValid() and photonMediumIdH.isValid() and photonTightIdH.isValid() and photonHighPtIdH.isValid()){
      
      for (size_t i = 0; i < tightphotons.size(); i++) {
        if (tightphotons[i]->pt() > hardestPhotonPt) {
	  hardestPhotonIndex = i;
	  hardestPhotonPt = tightphotons[i]->pt();
        }
      }


      nphotons = photonsH->size();
      
      if (hardestPhotonIndex >= 0) {
	phidl   = ((*photonLooseIdH )[tightphotons[hardestPhotonIndex]] ? 1 : 0);
	phidm   = ((*photonMediumIdH)[tightphotons[hardestPhotonIndex]] ? 1 : 0);
	phidt   = ((*photonTightIdH )[tightphotons[hardestPhotonIndex]] ? 1 : 0);
	phidh   = ((*photonHighPtIdH)[tightphotons[hardestPhotonIndex]] ? 1 : 0);
	phpt    = tightphotons[hardestPhotonIndex]->pt();
	pheta   = tightphotons[hardestPhotonIndex]->eta();
	phphi   = tightphotons[hardestPhotonIndex]->phi();
	phpt    = tightphotons[hardestPhotonIndex]->pt();
	pheta   = tightphotons[hardestPhotonIndex]->eta();
	phphi   = tightphotons[hardestPhotonIndex]->phi();
      }
    }
    
    // Substructure CHS
    if(addSubstructureCHS){      
      //sort collection to make sure it is ordered
      vector<pat::JetRef> jetsBoosted;
      if(boostedJetsH.isValid()){
	for (auto jets_iter = boostedJetsH->begin(); jets_iter != boostedJetsH->end(); ++jets_iter) {
	  //clean from leptons                                                                                                                                                  
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }
	  }
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }
	  }
	  
	  if (skipjet) continue;
	  
	  // apply jet id                                                                                                                                                 
	  bool passjetid = applyJetID(*jets_iter,"loose");
	  if (!passjetid)
	    continue;
	  
	  pat::JetRef jetref(boostedJetsH, jets_iter - boostedJetsH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    jetsBoosted.push_back(jetref);
	}
      
	// sort in pt
	if(jetsBoosted.size() > 0)
	  sort(jetsBoosted.begin(), jetsBoosted.end(), jetSorter);


	leadBoostedJetpt    = 0.0; leadBoostedJeteta  = 0.0;  
	leadBoostedJetphi   = 0.0; leadBoostedJetm    = 0.0;
	leadBoostedJetGenpt = 0.0; leadBoostedJetGenm = 0.0;
	leadBoostedJetGeneta = 0.0; leadBoostedJetGenphi = 0.0;

	leadBoostedJettau1  = 0.0; leadBoostedJettau2 = 0.0;
	leadBoostedJettau3  = 0.0; leadBoostedJettau4 = 0.0;
	leadBoostedJetecf1  = 0.0; leadBoostedJetecf2 = 0.0; leadBoostedJetecf3  = 0.0;
	leadBoostedJetHFlav = 0.0; leadBoostedJetPFlav = 0.0; leadBoostedJetQGL  = 0.0; leadBoostedJetBtag = 0.0, leadBoostedJetDoubleBtag = 0.0;
	
	leadBoostedJetBosonpt = 0.0; leadBoostedJetBosoneta = 0.0; leadBoostedJetBosonphi = 0.0; leadBoostedJetBosonm = 0.0;
	
	leadPrunedJetpt     = 0.0; leadPrunedJetm      = 0.0; leadPrunedJetGenpt = 0.0; leadPrunedJetGenm  = 0.0;
	leadPrunedJeteta    = 0.0; leadPrunedJetphi     = 0.0; leadPrunedJetGeneta = 0.0; leadPrunedJetGenphi  = 0.0;
	leadPrunedJetHFlav  = 0.0; leadPrunedJetPFlav  = 0.0; leadPrunedJetQGL   = 0.0; leadPrunedJetBtag  = 0.0; leadPrunedJetDoubleBtag  = 0.0;
	leadPrunedJetptraw  = 0.0; leadPrunedJetmraw = 0.0;

	leadSoftDropJetpt    = 0.0; leadSoftDropJetm = 0.0; leadSoftDropJetGenpt = 0.0; leadSoftDropJetGenm = 0.0; 
	leadSoftDropJetphi   = 0.0; leadSoftDropJeteta = 0.0; leadSoftDropJetGenphi = 0.0; leadSoftDropJetGeneta = 0.0; 
	leadSoftDropJetHFlav = 0.0; leadSoftDropJetPFlav = 0.0; leadSoftDropJetQGL = 0.0; leadSoftDropJetBtag = 0.0; leadSoftDropJetDoubleBtag = 0.0;
	leadSoftDropJetptraw  = 0.0; leadSoftDropJetmraw = 0.0;
	
	leadPrunedSubJetpt_1 = 0.0; leadPrunedSubJetm_1  = 0.0; leadPrunedSubJetphi_1 = 0.0; leadPrunedSubJeteta_1 = 0.0;
	leadPrunedSubJetHFlav_1 = 0.0; leadPrunedSubJetQGL_1 = 0.0; leadPrunedSubJetBtag_1 = 0.0;
	leadPrunedSubJetGenpt_1 = 0.0; leadPrunedSubJetGenm_1 = 0.0; leadPrunedSubJetPFlav_1 = 0.0;
	leadPrunedSubJetGenphi_1 = 0.0; leadPrunedSubJetGeneta_1 = 0.0; 
	leadPrunedSubJetptraw_1 = 0.0; leadPrunedSubJetmraw_1  = 0.0; 
	
	leadPrunedSubJetpt_2 = 0.0; leadPrunedSubJetm_2  = 0.0; leadPrunedSubJetphi_2 = 0.0; leadPrunedSubJeteta_2 = 0.0;
	leadPrunedSubJetHFlav_2 = 0.0; leadPrunedSubJetQGL_2 = 0.0; leadPrunedSubJetBtag_2 = 0.0;  leadPrunedSubJetPFlav_2 = 0.0;
	leadPrunedSubJetGenpt_2 = 0.0; leadPrunedSubJetGenm_2 = 0.0;
	leadPrunedSubJetGeneta_2 = 0.0; leadPrunedSubJetGenphi_2 = 0.0;
	leadPrunedSubJetptraw_2 = 0.0; leadPrunedSubJetmraw_2  = 0.0; 

	leadSoftDropSubJetpt_1 = 0.0; leadSoftDropSubJetm_1  = 0.0; leadSoftDropSubJetphi_1 = 0.0; leadSoftDropSubJeteta_1 = 0.0;
	leadSoftDropSubJetHFlav_1 = 0.0; leadSoftDropSubJetQGL_1 = 0.0; leadSoftDropSubJetBtag_1 = 0.0; leadSoftDropSubJetPFlav_1 = 0.0;
	leadSoftDropSubJetGenpt_1 = 0.0; leadSoftDropSubJetGenm_1 = 0.0;
	leadSoftDropSubJetGeneta_1 = 0.0; leadSoftDropSubJetGenphi_1 = 0.0;	
	leadSoftDropSubJetptraw_1 = 0.0; leadSoftDropSubJetmraw_1  = 0.0; 

	leadSoftDropSubJetpt_2 = 0.0; leadSoftDropSubJetm_2  = 0.0; leadSoftDropSubJetphi_2 = 0.0; leadSoftDropSubJeteta_2 = 0.0;
	leadSoftDropSubJetHFlav_2 = 0.0; leadSoftDropSubJetQGL_2 = 0.0; leadSoftDropSubJetBtag_2 = 0.0; leadSoftDropSubJetPFlav_2 = 0.0;
	leadSoftDropSubJetGenpt_2 = 0.0; leadSoftDropSubJetGenm_2 = 0.0;
	leadSoftDropSubJetGeneta_2 = 0.0; leadSoftDropSubJetGenphi_2 = 0.0;
	leadSoftDropSubJetptraw_2 = 0.0; leadSoftDropSubJetmraw_2  = 0.0; 
	
	if(jetsBoosted.size() > 0){ // basic AK8 info
	
	  leadBoostedJetpt  = jetsBoosted.front()->pt();
	  leadBoostedJeteta = jetsBoosted.front()->eta();
	  leadBoostedJetphi = jetsBoosted.front()->phi();
	  leadBoostedJetm   = jetsBoosted.front()->mass();	
	  leadBoostedJetBtag = jetsBoosted.front()->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	  leadBoostedJetDoubleBtag = jetsBoosted.front()->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");

	  // N-jettiness
	  if(jetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau1"))
	    leadBoostedJettau1 = jetsBoosted.front()->userFloat("NjettinessAK8PFJetsCHS:tau1");
	  
	  if(jetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau2"))
	    leadBoostedJettau2 = jetsBoosted.front()->userFloat("Njettiness"+boostedJetsCHSLabel+":tau2");
	  
	  if(jetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau3"))
	    leadBoostedJettau3 = jetsBoosted.front()->userFloat("Njettiness"+boostedJetsCHSLabel+":tau3");	

	  if(jetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsCHSLabel+":tau4"))
	    leadBoostedJettau4 = jetsBoosted.front()->userFloat("Njettiness"+boostedJetsCHSLabel+":tau4");

	  // ecf
	  if(jetsBoosted.front()->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf1"))
	    leadBoostedJetecf1 = jetsBoosted.front()->userFloat("ecf"+boostedJetsCHSLabel+":ecf1");
	  
	  if(jetsBoosted.front()->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf2"))
	    leadBoostedJetecf2 = jetsBoosted.front()->userFloat("ecf"+boostedJetsCHSLabel+":ecf2");
	  
	  if(jetsBoosted.front()->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf3"))
	    leadBoostedJetecf3 = jetsBoosted.front()->userFloat("ecf"+boostedJetsCHSLabel+":ecf3");

	  // gen jets
	  if(isMC){
	    if(jetsBoosted.front()->genJet()){ // gen AK8 jet
	      leadBoostedJetGenpt = jetsBoosted.front()->genJet()->pt();
	      leadBoostedJetGenm  = jetsBoosted.front()->genJet()->mass();
	      leadBoostedJetGeneta  = jetsBoosted.front()->genJet()->eta();
	      leadBoostedJetGenphi  = jetsBoosted.front()->genJet()->phi();
	    }
	    leadBoostedJetHFlav = jetsBoosted.front()->hadronFlavour();
	    leadBoostedJetPFlav = jetsBoosted.front()->partonFlavour();	  

	  }
	  
	  // QGL 
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"QGL:qgLikelihood"))
	    leadBoostedJetQGL = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"QGL:qgLikelihood");
	  
	  
	  // matched gen boson
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:pt"))
	    leadBoostedJetBosonpt = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"GenBosonMatched:pt");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:eta"))
	    leadBoostedJetBosoneta = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"GenBosonMatched:eta");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:phi"))
	    leadBoostedJetBosonphi = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"GenBosonMatched:phi");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:mass"))
	    leadBoostedJetBosonm = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"GenBosonMatched:mass");


	  // pruned matched jet
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:mass"))
	    leadPrunedJetm = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:mass");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:eta"))
	    leadPrunedJeteta = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:eta");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:phi"))
	    leadPrunedJetphi = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:phi");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pt"))
	    leadPrunedJetpt = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:pt");
	
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedQGLMatched:qgLikelihood"))
	    leadPrunedJetQGL = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedQGLMatched:qgLikelihood");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    leadPrunedJetBtag = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    leadPrunedJetDoubleBtag = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags");


	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass"))
	    leadPrunedJetmraw = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"))
	    leadPrunedJetptraw = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt");

	  if(isMC){
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:hadronFlavour"))
	      leadPrunedJetHFlav  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:hadronFlavour");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:partonFlavour"))
	      leadPrunedJetPFlav = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:partonFlavour");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genMass"))
	      leadPrunedJetGenm  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:genMass");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genPt"))
	      leadPrunedJetGenpt  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:genPt");	
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genEta"))
	      leadPrunedJetGeneta  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:genEta");	
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:genPhi"))
	      leadPrunedJetGenphi  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"PrunedMatched:genPhi");	
	  }

	  // soft drop matched jet
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:mass"))
	    leadSoftDropJetm = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:mass");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pt"))
	    leadSoftDropJetpt = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pt");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:phi"))
	    leadSoftDropJetphi = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:phi");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:eta"))
	    leadSoftDropJeteta = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:eta");
	
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropQGLMatched:qgLikelihood"))
	    leadSoftDropJetQGL = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropQGLMatched:qgLikelihood");
	  
	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    leadSoftDropJetBtag = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    leadSoftDropJetDoubleBtag = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass"))
	    leadSoftDropJetmraw = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass");

	  if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw"))
	    leadSoftDropJetptraw = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw");

	  if(isMC){
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:hadronFlavour"))
	      leadSoftDropJetHFlav  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:hadronFlavour");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:partonFlavour"))
	      leadSoftDropJetPFlav = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:partonFlavour");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genMass"))
	      leadSoftDropJetGenm  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genMass");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genPt"))
	      leadSoftDropJetGenpt  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genPt");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genEta"))
	      leadSoftDropJetGeneta  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genEta");
	    if(jetsBoosted.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:genPhi"))
	      leadSoftDropJetGenphi  = jetsBoosted.front()->userFloat(boostedJetsCHSLabel+"SoftDropMatched:genPhi");
	  }	
	  
	  // sub-jets pruned 
	  if(jetsBoosted.front()->hasSubjets("Pruned")){

	    pat::JetPtrCollection subjets = jetsBoosted.front()->subjets("Pruned");
	    if(subjets.size() > 0 ){
	      leadPrunedSubJetpt_1  = subjets.front()->pt(); 
	      leadPrunedSubJetm_1   = subjets.front()->mass(); 
	      leadPrunedSubJetphi_1 = subjets.front()->phi(); 
	      leadPrunedSubJeteta_1 = subjets.front()->eta();
	      leadPrunedSubJetBtag_1 = subjets.front()->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPrunedSubJetptraw_1  = subjets.front()->correctedP4(0).pt(); 
	      leadPrunedSubJetmraw_1   = subjets.front()->correctedP4(0).mass(); 
	      
	      if(subjets.front()->hasUserFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"))
		leadPrunedSubJetQGL_1 = subjets.front()->userFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood");
	      	      
	      if(isMC){
		leadPrunedSubJetHFlav_1 = subjets.front()->hadronFlavour(); 
		leadPrunedSubJetPFlav_1 = subjets.front()->partonFlavour(); 
		if(subjets.front()->genJet()){
		  leadPrunedSubJetGenpt_1 = subjets.front()->genJet()->pt();
		  leadPrunedSubJetGenm_1  = subjets.front()->genJet()->mass(); 
		  leadPrunedSubJetGeneta_1  = subjets.front()->genJet()->eta(); 
		  leadPrunedSubJetGenphi_1  = subjets.front()->genJet()->phi(); 
		}
	      }
	    }
	  
	    if(subjets.size() > 1 ){
	      leadPrunedSubJetpt_2  = subjets.at(1)->pt(); 
	      leadPrunedSubJetm_2   = subjets.at(1)->mass(); 
	      leadPrunedSubJetphi_2 = subjets.at(1)->phi(); 
	      leadPrunedSubJeteta_2 = subjets.at(1)->eta();
	      leadPrunedSubJetBtag_2 = subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPrunedSubJetptraw_2  = subjets.at(1)->correctedP4(0).pt(); 
	      leadPrunedSubJetmraw_2   = subjets.at(1)->correctedP4(0).mass(); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood"))
		leadPrunedSubJetQGL_2 = subjets.at(1)->userFloat(boostedJetsCHSLabel+"PrunedSubJetsQGL:qgLikelihood");
	      

	      if(isMC){
		leadPrunedSubJetHFlav_2 = subjets.at(1)->hadronFlavour(); 
		leadPrunedSubJetPFlav_2 = subjets.at(1)->partonFlavour(); 
		if(subjets.at(1)->genJet()){
		  leadPrunedSubJetGenpt_2 = subjets.at(1)->genJet()->pt();
		  leadPrunedSubJetGenm_2  = subjets.at(1)->genJet()->mass();
		  leadPrunedSubJetGeneta_2 = subjets.at(1)->genJet()->eta();
		  leadPrunedSubJetGenphi_2 = subjets.at(1)->genJet()->phi();
		}
	      }
	    }
	  }

	  // sub-jets soft drop 
	  if(jetsBoosted.front()->hasSubjets("SoftDrop")){
	    pat::JetPtrCollection subjets = jetsBoosted.front()->subjets("SoftDrop");
	    if(subjets.size() > 0 ){
	      leadSoftDropSubJetpt_1  = subjets.front()->pt(); 
	      leadSoftDropSubJetm_1   = subjets.front()->mass(); 
	      leadSoftDropSubJetphi_1 = subjets.front()->phi(); 
	      leadSoftDropSubJeteta_1 = subjets.front()->eta();
	      leadSoftDropSubJetBtag_1 = subjets.front()->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPrunedSubJetptraw_1  = subjets.front()->correctedP4(0).pt(); 
	      leadPrunedSubJetmraw_1   = subjets.front()->correctedP4(0).mass(); 
	    
	      if(subjets.front()->hasUserFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		leadSoftDropSubJetQGL_1 = subjets.front()->userFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood");
	    
	    
	      if(isMC){
		leadSoftDropSubJetHFlav_1 = subjets.front()->hadronFlavour(); 
		leadSoftDropSubJetPFlav_1 = subjets.front()->partonFlavour(); 
		if(subjets.front()->genJet()){
		  leadSoftDropSubJetGenpt_1 = subjets.front()->genJet()->pt();
		  leadSoftDropSubJetGenm_1  = subjets.front()->genJet()->mass(); 
		  leadSoftDropSubJetGeneta_1  = subjets.front()->genJet()->eta(); 
		  leadSoftDropSubJetGenphi_1  = subjets.front()->genJet()->phi(); 
		}
	      }
	    }
	    
	    if(subjets.size() > 1 ){
	      leadSoftDropSubJetpt_2  = subjets.at(1)->pt(); 
	      leadSoftDropSubJetm_2   = subjets.at(1)->mass(); 
	      leadSoftDropSubJetphi_2 = subjets.at(1)->phi(); 
	      leadSoftDropSubJeteta_2 = subjets.at(1)->eta();
	      leadSoftDropSubJetBtag_2 = subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPrunedSubJetptraw_2  = subjets.at(1)->correctedP4(0).pt(); 
	      leadPrunedSubJetmraw_2   = subjets.at(1)->correctedP4(0).mass(); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		leadSoftDropSubJetQGL_2 = subjets.at(1)->userFloat(boostedJetsCHSLabel+"SoftDropSubJetsQGL:qgLikelihood");
	      
	      if(isMC){
		leadSoftDropSubJetHFlav_2 = subjets.at(1)->hadronFlavour(); 
		leadSoftDropSubJetPFlav_2 = subjets.at(1)->partonFlavour(); 
		if(subjets.at(1)->genJet()){
		  leadSoftDropSubJetGenpt_2 = subjets.at(1)->genJet()->pt();
		  leadSoftDropSubJetGenm_2  = subjets.at(1)->genJet()->mass();
		  leadSoftDropSubJetGeneta_2  = subjets.at(1)->genJet()->eta();
		  leadSoftDropSubJetGenphi_2  = subjets.at(1)->genJet()->phi();
		}
	      }
	    }	  
	  }
	}
      }
    }
    
    // Substructure Puppi
    if(addSubstructurePuppi){      
      //sort collection to make sure it is ordered
      vector<pat::JetRef> puppiJetsBoosted;
      
      if(boostedPuppiJetsH.isValid()){
	for (auto jets_iter = boostedPuppiJetsH->begin(); jets_iter != boostedPuppiJetsH->end(); ++jets_iter) {
	  
	  //clean from leptons                                                                                                                                                
	  bool skipjet = false;
	  if(muonsH.isValid()){
	    for (std::size_t j = 0; j < muons.size(); j++) {
	      if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }
	  }
	  
	  if(electronsH.isValid()){
	    for (std::size_t j = 0; j < electrons.size(); j++) {
	      if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }
	  }
	  if(photonsH.isValid()){
	    for (std::size_t j = 0; j < photons.size(); j++) {
	      if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.8) 
		skipjet = true;
	    }	  
	  }
	  
	  if (skipjet) continue;

	  // apply jet id                                                                                                                                                    
	  bool passjetid = applyJetID(*jets_iter,"loose");
	  if (!passjetid)
	    continue;
	  
	  pat::JetRef jetref(boostedPuppiJetsH, jets_iter - boostedPuppiJetsH->begin());
	  if(jetref.isAvailable() and jetref.isNonnull())
	    puppiJetsBoosted.push_back(jetref);
	}	
	// sort in pt
	if(puppiJetsBoosted.size() > 0)
	  sort(puppiJetsBoosted.begin(), puppiJetsBoosted.end(), jetSorter);
	
	leadPuppiBoostedJetpt    = 0.0; leadPuppiBoostedJeteta  = 0.0;  
	leadPuppiBoostedJetphi   = 0.0; leadPuppiBoostedJetm    = 0.0;
	leadPuppiBoostedJetGenpt = 0.0; leadPuppiBoostedJetGenm = 0.0;
	leadPuppiBoostedJetGeneta = 0.0; leadPuppiBoostedJetGenphi = 0.0;
	leadPuppiBoostedJettau1  = 0.0; leadPuppiBoostedJettau2 = 0.0;
	leadPuppiBoostedJettau3  = 0.0; leadPuppiBoostedJettau4 = 0.0;
	leadPuppiBoostedJetecf1  = 0.0; leadPuppiBoostedJetecf2 = 0.0; leadPuppiBoostedJetecf3  = 0.0;
	leadPuppiBoostedJetHFlav = 0.0; leadPuppiBoostedJetPFlav = 0.0; leadPuppiBoostedJetQGL  = 0.0; 
	leadPuppiBoostedJetBtag = 0.0; leadPuppiBoostedJetDoubleBtag = 0.0;
	
	leadPuppiBoostedJetBosonpt = 0.0; leadPuppiBoostedJetBosoneta = 0.0; leadPuppiBoostedJetBosonphi = 0.0; leadPuppiBoostedJetBosonm = 0.0;
	
	leadPuppiPrunedJetpt     = 0.0; leadPuppiPrunedJetm      = 0.0; leadPuppiPrunedJetGenpt = 0.0; leadPuppiPrunedJetGenm  = 0.0;
	leadPuppiPrunedJeteta     = 0.0; leadPuppiPrunedJetphi   = 0.0; leadPuppiPrunedJetGeneta = 0.0; leadPuppiPrunedJetGenphi  = 0.0;
	leadPuppiPrunedJetHFlav  = 0.0; leadPuppiPrunedJetPFlav  = 0.0; leadPuppiPrunedJetQGL   = 0.0; leadPuppiPrunedJetBtag  = 0.0;
	leadPuppiPrunedJetDoubleBtag  = 0.0;
	
	leadPuppiSoftDropJetpt    = 0.0; leadPuppiSoftDropJetm = 0.0; leadPuppiSoftDropJetGenpt = 0.0; leadPuppiSoftDropJetGenm = 0.0; 
	leadPuppiSoftDropJeteta    = 0.0; leadPuppiSoftDropJetphi = 0.0; leadPuppiSoftDropJetGeneta = 0.0; leadPuppiSoftDropJetGenphi = 0.0; 
	leadPuppiSoftDropJetHFlav = 0.0; leadPuppiSoftDropJetPFlav = 0.0; leadPuppiSoftDropJetQGL = 0.0; leadPuppiSoftDropJetBtag = 0.0;
	leadPuppiSoftDropJetDoubleBtag = 0.0;
	
	leadPuppiPrunedSubJetpt_1 = 0.0; leadPuppiPrunedSubJetm_1  = 0.0; leadPuppiPrunedSubJetphi_1 = 0.0; leadPuppiPrunedSubJeteta_1 = 0.0;
	leadPuppiPrunedSubJetHFlav_1 = 0.0; leadPuppiPrunedSubJetPFlav_1 = 0.0; leadPuppiPrunedSubJetQGL_1 = 0.0; leadPuppiPrunedSubJetBtag_1 = 0.0;
	leadPuppiPrunedSubJetGenpt_1 = 0.0; leadPuppiPrunedSubJetGenm_1 = 0.0;
	leadPuppiPrunedSubJetGenphi_1 = 0.0; leadPuppiPrunedSubJetGeneta_1 = 0.0;
	leadPuppiPrunedSubJetptraw_1 = 0.0; leadPuppiPrunedSubJetmraw_1  = 0.0; 

	leadPuppiPrunedSubJetpt_2 = 0.0; leadPuppiPrunedSubJetm_2  = 0.0; leadPuppiPrunedSubJetphi_2 = 0.0; leadPuppiPrunedSubJeteta_2 = 0.0;
	leadPuppiPrunedSubJetHFlav_2 = 0.0;leadPuppiPrunedSubJetPFlav_2 = 0.0; leadPuppiPrunedSubJetQGL_2 = 0.0; leadPuppiPrunedSubJetBtag_2 = 0.0;
	leadPuppiPrunedSubJetGenpt_2 = 0.0; leadPuppiPrunedSubJetGenm_2 = 0.0;
	leadPuppiPrunedSubJetGenphi_2 = 0.0; leadPuppiPrunedSubJetGeneta_2 = 0.0;
	leadPuppiPrunedSubJetptraw_2 = 0.0; leadPuppiPrunedSubJetmraw_2  = 0.0; 
	
	leadPuppiSoftDropSubJetpt_1 = 0.0; leadPuppiSoftDropSubJetm_1  = 0.0; leadPuppiSoftDropSubJetphi_1 = 0.0; leadPuppiSoftDropSubJeteta_1 = 0.0;
	leadPuppiSoftDropSubJetHFlav_1 = 0.0; leadPuppiSoftDropSubJetQGL_1 = 0.0; leadPuppiSoftDropSubJetBtag_1 = 0.0; leadPuppiSoftDropSubJetPFlav_1 = 0.0;
	leadPuppiSoftDropSubJetGenpt_1 = 0.0; leadPuppiSoftDropSubJetGenm_1 = 0.0;
	leadPuppiSoftDropSubJetGenphi_1 = 0.0; leadPuppiSoftDropSubJetGeneta_1 = 0.0;
	leadPuppiSoftDropSubJetptraw_1 = 0.0; leadPuppiSoftDropSubJetmraw_1  = 0.0; 

	leadPuppiSoftDropSubJetpt_2 = 0.0; leadPuppiSoftDropSubJetm_2  = 0.0; leadPuppiSoftDropSubJetphi_2 = 0.0; leadPuppiSoftDropSubJeteta_2 = 0.0;
	leadPuppiSoftDropSubJetHFlav_2 = 0.0; leadPuppiSoftDropSubJetQGL_2 = 0.0; leadPuppiSoftDropSubJetBtag_2 = 0.0; leadPuppiSoftDropSubJetPFlav_2 = 0.0;
	leadPuppiSoftDropSubJetGenpt_2 = 0.0; leadPuppiSoftDropSubJetGenm_2 = 0.0;
	leadPuppiSoftDropSubJetGenphi_2 = 0.0; leadPuppiSoftDropSubJetGeneta_2 = 0.0;
	leadPuppiSoftDropSubJetptraw_2 = 0.0; leadPuppiSoftDropSubJetmraw_2  = 0.0; 
	

	if(puppiJetsBoosted.size() > 0){ // basic AK8 info
	  
	  leadPuppiBoostedJetpt  = puppiJetsBoosted.front()->pt();
	  leadPuppiBoostedJeteta = puppiJetsBoosted.front()->eta();
	  leadPuppiBoostedJetphi = puppiJetsBoosted.front()->phi();
	  leadPuppiBoostedJetm   = puppiJetsBoosted.front()->mass();
	  leadPuppiBoostedJetBtag = puppiJetsBoosted.front()->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
	  leadPuppiBoostedJetDoubleBtag = puppiJetsBoosted.front()->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
	  
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau1"))
	    leadPuppiBoostedJettau1 = puppiJetsBoosted.front()->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau1");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau2"))
	    leadPuppiBoostedJettau2 = puppiJetsBoosted.front()->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau2");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau3"))
	    leadPuppiBoostedJettau3 = puppiJetsBoosted.front()->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau3");	
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("Njettiness"+boostedJetsPuppiLabel+":tau4"))
	    leadPuppiBoostedJettau4 = puppiJetsBoosted.front()->userFloat("Njettiness"+boostedJetsPuppiLabel+":tau4");
	  	  
	  if(puppiJetsBoosted.front()->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf1")){
	    leadPuppiBoostedJetecf1 = puppiJetsBoosted.front()->userFloat("ecf"+boostedJetsPuppiLabel+":ecf1");
	  }
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf2")){
	    leadPuppiBoostedJetecf2 = puppiJetsBoosted.front()->userFloat("ecf"+boostedJetsPuppiLabel+":ecf2");
	  }
	  
	  if(puppiJetsBoosted.front()->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf3")){
	    leadPuppiBoostedJetecf3 = puppiJetsBoosted.front()->userFloat("ecf"+boostedJetsPuppiLabel+":ecf3");
	  }
	  
	  if(isMC){
	    if(puppiJetsBoosted.front()->genJet()){ // gen AK8 jet
	      leadPuppiBoostedJetGenpt = puppiJetsBoosted.front()->genJet()->pt();
	      leadPuppiBoostedJetGenm  = puppiJetsBoosted.front()->genJet()->mass();
	      leadPuppiBoostedJetGeneta  = puppiJetsBoosted.front()->genJet()->eta();
	      leadPuppiBoostedJetGenphi  = puppiJetsBoosted.front()->genJet()->phi();
	    }
	    leadPuppiBoostedJetHFlav = puppiJetsBoosted.front()->hadronFlavour();
	    leadPuppiBoostedJetPFlav = puppiJetsBoosted.front()->partonFlavour();	  
	  }
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"QGL:qgLikelihood"))
	    leadPuppiBoostedJetQGL = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"QGL:qgLikelihood");
	  
	  // matched gen boson	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:pt"))
	    leadPuppiBoostedJetBosonpt = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:pt");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:eta"))
	    leadPuppiBoostedJetBosoneta = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:eta");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:phi"))
	    leadPuppiBoostedJetBosonphi = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:phi");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:mass"))
	    leadPuppiBoostedJetBosonm = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:mass");
	  
	  // pruned matched jet
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:mass"))
	    leadPuppiPrunedJetm = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:mass");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pt"))
	    leadPuppiPrunedJetpt = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pt");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:eta"))
	    leadPuppiPrunedJeteta = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:eta");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:phi"))
	    leadPuppiPrunedJetphi = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:phi");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedQGLMatched:qgLikelihood"))
	    leadPuppiPrunedJetQGL = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedQGLMatched:qgLikelihood");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    leadPuppiPrunedJetBtag = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    leadPuppiPrunedJetDoubleBtag = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass"))
	    leadPuppiPrunedJetmraw = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"))
	    leadPuppiPrunedJetptraw = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt");
	  
	  if(isMC){
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:hadronFlavour"))
	      leadPuppiPrunedJetHFlav  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:hadronFlavour");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:partonFlavour"))
	     leadPuppiPrunedJetPFlav = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:partonFlavour");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genMass"))
	      leadPuppiPrunedJetGenm  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genMass");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genPt"))
	      leadPuppiPrunedJetGenpt  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genPt");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genEta"))
	      leadPuppiPrunedJetGeneta  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genEta");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:genPhi"))
	      leadPuppiPrunedJetGenphi  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"PrunedMatched:genPhi");
	  }
      
	  
	  // soft drop matched jet
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:mass"))
	    leadPuppiSoftDropJetm = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:mass");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pt"))
	    leadPuppiSoftDropJetpt = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pt");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:eta"))
	    leadPuppiSoftDropJeteta = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:eta");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:phi"))
	    leadPuppiSoftDropJetphi = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:phi");
	
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropQGLMatched:qgLikelihood"))
	    leadPuppiSoftDropJetQGL = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropQGLMatched:qgLikelihood");
	  
	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	    leadPuppiSoftDropJetBtag = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	    leadPuppiSoftDropJetDoubleBtag = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass"))
	    leadPuppiSoftDropJetmraw = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass");

	  if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"))
	    leadPuppiSoftDropJetptraw = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt");
	  
	  if(isMC){
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:hadronFlavour"))
	      leadPuppiSoftDropJetHFlav  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:hadronFlavour");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:partonFlavour"))
	      leadPuppiSoftDropJetPFlav = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:partonFlavour");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genMass"))
	      leadPuppiSoftDropJetGenm  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genMass");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPt"))
	      leadPuppiSoftDropJetGenpt  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPt");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genEta"))
	      leadPuppiSoftDropJetGeneta  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genEta");
	    if(puppiJetsBoosted.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPhi"))
	      leadPuppiSoftDropJetGenphi  = puppiJetsBoosted.front()->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:genPhi");	    
	  }
	
	  // sub-jets pruned 
	  if(puppiJetsBoosted.front()->hasSubjets("Pruned")){
	    pat::JetPtrCollection subjets = puppiJetsBoosted.front()->subjets("Pruned");
	    
	    if(subjets.size() > 0){
	      leadPuppiPrunedSubJetpt_1  = subjets.front()->pt(); 
	      leadPuppiPrunedSubJetm_1   = subjets.front()->mass(); 
	      leadPuppiPrunedSubJetphi_1 = subjets.front()->phi(); 
	      leadPuppiPrunedSubJeteta_1 = subjets.front()->eta();
	      leadPuppiPrunedSubJetBtag_1 = subjets.front()->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPuppiPrunedSubJetptraw_1  = subjets.front()->correctedP4(0).pt(); 
	      leadPuppiPrunedSubJetmraw_1   = subjets.front()->correctedP4(0).mass(); 
	      
	      if(subjets.front()->hasUserFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"))
		leadPuppiPrunedSubJetQGL_1 = subjets.front()->userFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood");
	      
	      if(isMC){
		leadPuppiPrunedSubJetHFlav_1 = subjets.front()->hadronFlavour(); 
		leadPuppiPrunedSubJetPFlav_1 = subjets.front()->partonFlavour(); 
		if(subjets.front()->genJet()){
		  leadPuppiPrunedSubJetGenpt_1 = subjets.front()->genJet()->pt();
		  leadPuppiPrunedSubJetGenm_1  = subjets.front()->genJet()->mass();
		  leadPuppiPrunedSubJetGeneta_1 = subjets.front()->genJet()->eta();
		  leadPuppiPrunedSubJetGenphi_1  = subjets.front()->genJet()->phi();
		}
	      }
	    }
	    
	    if(subjets.size() > 1){
	      leadPuppiPrunedSubJetpt_2  = subjets.at(1)->pt(); 
	      leadPuppiPrunedSubJetm_2   = subjets.at(1)->mass(); 
	      leadPuppiPrunedSubJetphi_2 = subjets.at(1)->phi(); 
	      leadPuppiPrunedSubJeteta_2 = subjets.at(1)->eta();
	      leadPuppiPrunedSubJetBtag_2 = subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPuppiPrunedSubJetptraw_2  = subjets.at(1)->correctedP4(0).pt(); 
	      leadPuppiPrunedSubJetmraw_2   = subjets.at(1)->correctedP4(0).mass(); 
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood"))
		leadPuppiPrunedSubJetQGL_2 = subjets.at(1)->userFloat(boostedJetsPuppiLabel+"PrunedSubJetsQGL:qgLikelihood");
	      
	      if(isMC){
		leadPuppiPrunedSubJetHFlav_2 = subjets.at(1)->hadronFlavour(); 
		leadPuppiPrunedSubJetPFlav_2 = subjets.at(1)->partonFlavour(); 
		if(subjets.at(1)->genJet()){
		  leadPuppiPrunedSubJetGenpt_2 = subjets.at(1)->genJet()->pt();
		  leadPuppiPrunedSubJetGenm_2  = subjets.at(1)->genJet()->mass();
		  leadPuppiPrunedSubJetGenphi_2  = subjets.at(1)->genJet()->phi();
		  leadPuppiPrunedSubJetGeneta_2  = subjets.at(1)->genJet()->eta();
		}	      
	      }
	    }
	  }
	  
	  
	  // sub-jets soft drop 
	  if(puppiJetsBoosted.front()->hasSubjets("SoftDrop")){
	    pat::JetPtrCollection subjets = puppiJetsBoosted.front()->subjets("SoftDrop");
	    if(subjets.size() > 0){
	      leadPuppiSoftDropSubJetpt_1  = subjets.front()->pt(); 
	      leadPuppiSoftDropSubJetm_1   = subjets.front()->mass(); 
	      leadPuppiSoftDropSubJetphi_1 = subjets.front()->phi(); 
	      leadPuppiSoftDropSubJeteta_1 = subjets.front()->eta();
	      leadPuppiSoftDropSubJetBtag_1 = subjets.front()->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	      leadPuppiSoftDropSubJetptraw_1  = subjets.front()->correctedP4(0).pt(); 
	      leadPuppiSoftDropSubJetmraw_1   = subjets.front()->correctedP4(0).mass(); 
	    
	      if(subjets.front()->hasUserFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		leadPuppiSoftDropSubJetQGL_1 = subjets.front()->userFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood");
	      
	      if(isMC){
		leadPuppiSoftDropSubJetHFlav_1 = subjets.front()->hadronFlavour(); 
		if(subjets.front()->genJet()){
		  leadPuppiSoftDropSubJetGenpt_1 = subjets.front()->genJet()->pt();
		  leadPuppiSoftDropSubJetGenm_1  = subjets.front()->genJet()->mass(); 
		  leadPuppiSoftDropSubJetGeneta_1  = subjets.front()->genJet()->eta(); 
		  leadPuppiSoftDropSubJetGenphi_1  = subjets.front()->genJet()->phi(); 
		}
	      }
	    }
	  
	    if(subjets.size() > 1){
	      leadPuppiSoftDropSubJetpt_2  = subjets.at(1)->pt(); 
	      leadPuppiSoftDropSubJetm_2   = subjets.at(1)->mass(); 
	      leadPuppiSoftDropSubJetphi_2 = subjets.at(1)->phi(); 
	      leadPuppiSoftDropSubJeteta_2 = subjets.at(1)->eta();
	      leadPuppiSoftDropSubJetBtag_2 = subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	      
	      if(subjets.at(1)->hasUserFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood"))
		leadPuppiSoftDropSubJetQGL_2 = subjets.at(1)->userFloat(boostedJetsPuppiLabel+"SoftDropSubJetsQGL:qgLikelihood");

	      leadPuppiSoftDropSubJetptraw_2  = subjets.at(1)->correctedP4(0).pt(); 
	      leadPuppiSoftDropSubJetmraw_2   = subjets.at(1)->correctedP4(0).mass(); 
	      
	      if(isMC){
		leadPuppiSoftDropSubJetHFlav_2 = subjets.at(1)->hadronFlavour(); 
		if(subjets.at(1)->genJet()){
		  leadPuppiSoftDropSubJetGenpt_2 = subjets.at(1)->genJet()->pt();
		  leadPuppiSoftDropSubJetGenm_2  = subjets.at(1)->genJet()->mass();
		  leadPuppiSoftDropSubJetGeneta_2  = subjets.at(1)->genJet()->eta();
		  leadPuppiSoftDropSubJetGenphi_2  = subjets.at(1)->genJet()->phi();
		}
	      }
	    }
	  }
	}
      }
    }

    // Generator-level information
    wzid          = 0; wzmass        = 0.0; wzpt          = 0.0; wzeta         = 0.0; wzphi         = 0.0;
    l1id          = 0; l1pt          = 0.0; l1eta         = 0.0; l1phi         = 0.0;
    l2id          = 0; l2pt          = 0.0; l2eta         = 0.0; l2phi         = 0.0;
    parid         = 0; parpt         = 0.0; pareta        = 0.0; parphi        = 0.0;
    ancid         = 0; ancpt         = 0.0; anceta        = 0.0; ancphi        = 0.0;

    if (isWorZMCSample && gensH.isValid()) {
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) { // loop on genParticles (prunedGenParticles)
	if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
	    gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
	    abs(gens_iter->daughter(0)->pdgId()) > 10 && 
	    abs(gens_iter->daughter(0)->pdgId()) < 17)  { // leptons 

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
	}
      }
    }

    // only for signal samples
    if (isSignalSample && gensH.isValid()) {
        TLorentzVector dm1vec; 
        TLorentzVector dm2vec; 
        bool foundfirst = false;
        for(auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	  if (gens_iter->pdgId() == 1000022 && !foundfirst) { // first DM particle
	    dm1vec.SetPtEtaPhiE(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->p());
	    foundfirst = true;
	  }
	  if (gens_iter->pdgId() == 1000022 &&  foundfirst) { // second DM particle
	    dm2vec.SetPtEtaPhiE(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->p());
	    break;
	  }
        }
        TLorentzVector medvec(dm1vec);
        medvec += dm2vec;
        wzpt  = medvec.Pt();
        wzeta = medvec.Eta();
        wzphi = medvec.Phi();
    }
    tree->Fill();    
  
}


void MonoJetTreeMaker::beginJob() {

  edm::Service<TFileService> fs;
  tree = std::auto_ptr<TTree>(fs->make<TTree>("tree"       , "tree"));

  // Run, Lumi, Event info
  tree->Branch("event"                , &event                , "event/i");
  tree->Branch("run"                  , &run                  , "run/i");
  tree->Branch("lumi"                 , &lumi                 , "lumi/i");
  // Event weights
  tree->Branch("xsec"                 , &xsec                 , "xsec/D");
  tree->Branch("wgt"                  , &wgt                  , "wgt/D");
  // Pileup info
  tree->Branch("puwgt"                , &puwgt                , "puwgt/D");
  tree->Branch("puobs"                , &puobs                , "puobs/I");
  tree->Branch("putrue"               , &putrue               , "putrue/I");
  tree->Branch("nvtx"                 , &nvtx                 , "nvtx/i");
  
  // Triggers
  tree->Branch("hltmet90"             , &hltmet90             , "hltmet90/b");
  tree->Branch("hltmet120"            , &hltmet120            , "hltmet120/b");
  tree->Branch("hltmetwithmu90"       , &hltmetwithmu90       , "hltmetwithmu90/b");
  tree->Branch("hltmetwithmu120"      , &hltmetwithmu120      , "hltmetwithmu120/b");
  tree->Branch("hltmetwithmu170"      , &hltmetwithmu170      , "hltmetwithmu170/b");
  tree->Branch("hltmetwithmu300"      , &hltmetwithmu300      , "hltmetwithmu300/b");
  tree->Branch("hltjetmet90"          , &hltjetmet90          , "hltjetmet90/b");
  tree->Branch("hltjetmet120"         , &hltjetmet120         , "hltjetmet120/b");
  tree->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/b");
  tree->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/b");
  tree->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/b");
  tree->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/b");
  tree->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/b");
  tree->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/b");

  // MET filters
  tree->Branch("flagcsctight"         , &flagcsctight         , "flagcsctight/b");
  tree->Branch("flaghbhenoise"        , &flaghbhenoise        , "flaghbhenoise/b");
  tree->Branch("flaghbheloose"        , &flaghbheloose        , "flaghbheloose/b");
  tree->Branch("flaghbhetight"        , &flaghbhetight        , "flaghbhetight/b");
  tree->Branch("flaghbheiso"          , &flaghbheiso          , "flaghbheiso/b");
  tree->Branch("flageebadsc"          , &flageebadsc          , "flageebadsc/b");

  // Object counts
  tree->Branch("nmuons"               , &nmuons               , "nmuons/i");
  tree->Branch("nelectrons"           , &nelectrons           , "nelectrons/i");
  tree->Branch("ntightmuons"          , &ntightmuons          , "ntightmuons/i");
  tree->Branch("nhighptmuons"         , &nhighptmuons         , "nhighptmuons/i");
  tree->Branch("ntightelectrons"      , &ntightelectrons      , "ntightelectrons/i");
  tree->Branch("nheepelectrons"       , &nheepelectrons       , "nheepelectrons/i");
  tree->Branch("ntaus"                , &ntaus                , "ntaus/i");
  tree->Branch("nphotons"             , &nphotons             , "nphotons/i");
  tree->Branch("njets"                , &njets                , "njets/i");
  tree->Branch("nbjets"               , &nbjets               , "nbjets/i");
  tree->Branch("nbjetslowpt"          , &nbjetslowpt          , "nbjetslowpt/i");

  if(addPuppiJets){
    tree->Branch("npuppijets"                , &npuppijets                , "npuppijets/i");
    tree->Branch("npuppibjets"               , &npuppibjets               , "npuppibjets/i");
    tree->Branch("npuppibjetslowpt"          , &npuppibjetslowpt          , "npuppibjetslowpt/i");
  }

  // MET
  tree->Branch("pfmet"                , &pfmet                , "pfmet/D");
  tree->Branch("pfmetphi"             , &pfmetphi             , "pfmetphi/D");
  tree->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/D");
  tree->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/D");
  tree->Branch("mumet"                , &mumet                , "mumet/D");
  tree->Branch("mumetphi"             , &mumetphi             , "mumetphi/D");
  tree->Branch("t1mumet"              , &t1mumet              , "t1mumet/D");
  tree->Branch("t1mumetphi"           , &t1mumetphi           , "t1mumetphi/D");
  tree->Branch("elmet"                , &elmet                , "elmet/D");
  tree->Branch("elmetphi"             , &elmetphi             , "elmetphi/D");
  tree->Branch("t1elmet"              , &t1elmet              , "t1elmet/D");
  tree->Branch("t1elmetphi"           , &t1elmetphi           , "t1elmetphi/D");
  tree->Branch("phmet"                , &phmet                , "phmet/D");
  tree->Branch("phmetphi"             , &phmetphi             , "phmetphi/D");
  tree->Branch("t1phmet"              , &t1phmet              , "t1phmet/D");
  tree->Branch("t1phmetphi"           , &t1phmetphi           , "t1phmetphi/D");

  tree->Branch("genmet",&genmet,"genmet/D");

  if(addMVAMet){
    tree->Branch("mvamet"              , &mvamet              , "mvamet/D");
    tree->Branch("mvametphi"           , &mvametphi           , "mvametphi/D");
  }

  if(addMETSystematics){
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
  }

  if(addPuppiMET){

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

    if(addMETSystematics){

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
    }

  }
  
  // Jet info
  tree->Branch("leadingjetpt"         , &leadingjetpt         , "leadingjetpt/D");
  tree->Branch("leadingjeteta"        , &leadingjeteta        , "leadingjeteta/D");
  tree->Branch("leadingjetphi"        , &leadingjetphi        , "leadingjetphi/D");
  tree->Branch("leadingjetm"          , &leadingjetm          , "leadingjetm/D");

  tree->Branch("signaljetpt"          , &signaljetpt          , "signaljetpt/D");
  tree->Branch("signaljeteta"         , &signaljeteta         , "signaljeteta/D");
  tree->Branch("signaljetphi"         , &signaljetphi         , "signaljetphi/D");
  tree->Branch("signaljetbtag"        , &signaljetbtag        , "signaljetbtag/D");
  tree->Branch("signaljetCHfrac"      , &signaljetCHfrac      , "signaljetCHfrac/D");
  tree->Branch("signaljetNHfrac"      , &signaljetNHfrac      , "signaljetNHfrac/D");
  tree->Branch("signaljetEMfrac"      , &signaljetEMfrac      , "signaljetEMfrac/D");
  tree->Branch("signaljetCEMfrac"     , &signaljetCEMfrac     , "signaljetCEMfrac/D");
  tree->Branch("signaljetmetdphi"     , &signaljetmetdphi     , "signaljetmetdphi/D");
  tree->Branch("signaljetHFlav"       , &signaljetHFlav       , "signaljetHFlav/D");
  tree->Branch("signaljetPFlav"       , &signaljetPFlav       , "signaljetPFlav/D");
  tree->Branch("signaljetQGL"         , &signaljetQGL         , "signaljetQGL/D");
  tree->Branch("signaljetPUID"        , &signaljetPUID        , "signaljetPUID/D");
  tree->Branch("signaljetGenpt"       , &signaljetGenpt       , "signaljetGenpt/D");
  tree->Branch("signaljetGeneta"      , &signaljetGeneta      , "signaljetGeneta/D");
  tree->Branch("signaljetGenphi"      , &signaljetGenphi      , "signaljetGenphi/D");
  tree->Branch("signaljetGenm"        , &signaljetGenm        , "signaljetGenm/D");
  tree->Branch("signaljetRawpt"       , &signaljetRawpt       , "signaljetRawpt/D");


  tree->Branch("secondjetpt"          , &secondjetpt          , "secondjetpt/D");
  tree->Branch("secondjeteta"         , &secondjeteta         , "secondjeteta/D");
  tree->Branch("secondjetphi"         , &secondjetphi         , "secondjetphi/D");
  tree->Branch("secondjetbtag"        , &secondjetbtag        , "secondjetbtag/D");
  tree->Branch("secondjetCHfrac"      , &secondjetCHfrac      , "secondjetCHfrac/D");
  tree->Branch("secondjetNHfrac"      , &secondjetNHfrac      , "secondjetNHfrac/D");
  tree->Branch("secondjetEMfrac"      , &secondjetEMfrac      , "secondjetEMfrac/D");
  tree->Branch("secondjetCEMfrac"     , &secondjetCEMfrac     , "secondjetCEMfrac/D");
  tree->Branch("secondjetmetdphi"     , &secondjetmetdphi     , "secondjetmetdphi/D");
  tree->Branch("secondjetHFlav"       , &secondjetHFlav       , "secondjetHFlav/D");
  tree->Branch("secondjetPFlav"       , &secondjetPFlav       , "secondjetPFlav/D");
  tree->Branch("secondjetQGL"         , &secondjetQGL         , "secondjetQGL/D");
  tree->Branch("secondjetPUID"        , &secondjetPUID        , "secondjetPUID/D");
  tree->Branch("secondjetGenpt"       , &secondjetGenpt       , "secondjetGenpt/D");
  tree->Branch("secondjetGeneta"      , &secondjetGeneta      , "secondjetGeneta/D");
  tree->Branch("secondjetGenphi"      , &secondjetGenphi      , "secondjetGenphi/D");
  tree->Branch("secondjetGenm"        , &secondjetGenm        , "secondjetGenm/D");
  tree->Branch("secondjetRawpt"       , &secondjetRawpt       , "secondjetRawpt/D");

  tree->Branch("thirdjetpt"           , &thirdjetpt           , "thirdjetpt/D");
  tree->Branch("thirdjeteta"          , &thirdjeteta          , "thirdjeteta/D");
  tree->Branch("thirdjetphi"          , &thirdjetphi          , "thirdjetphi/D");
  tree->Branch("thirdjetbtag"         , &thirdjetbtag         , "thirdjetbtag/D");
  tree->Branch("thirdjetCHfrac"       , &thirdjetCHfrac       , "thirdjetCHfrac/D");
  tree->Branch("thirdjetNHfrac"       , &thirdjetNHfrac       , "thirdjetNHfrac/D");
  tree->Branch("thirdjetEMfrac"       , &thirdjetEMfrac       , "thirdjetEMfrac/D");
  tree->Branch("thirdjetCEMfrac"      , &thirdjetCEMfrac      , "thirdjetCEMfrac/D");
  tree->Branch("thirdjetmetdphi"      , &thirdjetmetdphi      , "thirdjetmetdphi/D");
  tree->Branch("thirdjetHFlav"        , &thirdjetHFlav        , "thirdjetHFlav/D");
  tree->Branch("thirdjetPFlav"        , &thirdjetPFlav        , "thirdjetPFlav/D");
  tree->Branch("thirdjetQGL"          , &thirdjetQGL          , "thirdjetQGL/D");
  tree->Branch("thirdjetPUID"         , &thirdjetPUID         , "thirdjetPUID/D");
  tree->Branch("thirdjetGenpt"        , &thirdjetGenpt        , "thirdjetGenpt/D");
  tree->Branch("thirdjetGeneta"       , &thirdjetGeneta       , "thirdjetGeneta/D");
  tree->Branch("thirdjetGenphi"       , &thirdjetGenphi       , "thirdjetGenphi/D");
  tree->Branch("thirdjetGenm"         , &thirdjetGenm         , "thirdjetGenm/D");
  tree->Branch("thirdjetRawpt"        , &thirdjetRawpt        , "thirdjetRawpt/D");

  tree->Branch("fourthjetpt"           , &fourthjetpt           , "fourthjetpt/D");
  tree->Branch("fourthjeteta"          , &fourthjeteta          , "fourthjeteta/D");
  tree->Branch("fourthjetphi"          , &fourthjetphi          , "fourthjetphi/D");
  tree->Branch("fourthjetbtag"         , &fourthjetbtag         , "fourthjetbtag/D");
  tree->Branch("fourthjetCHfrac"       , &fourthjetCHfrac       , "fourthjetCHfrac/D");
  tree->Branch("fourthjetNHfrac"       , &fourthjetNHfrac       , "fourthjetNHfrac/D");
  tree->Branch("fourthjetEMfrac"       , &fourthjetEMfrac       , "fourthjetEMfrac/D");
  tree->Branch("fourthjetCEMfrac"      , &fourthjetCEMfrac      , "fourthjetCEMfrac/D");
  tree->Branch("fourthjetmetdphi"      , &fourthjetmetdphi      , "fourthjetmetdphi/D");
  tree->Branch("fourthjetHFlav"        , &fourthjetHFlav        , "fourthjetHFlav/D");
  tree->Branch("fourthjetPFlav"        , &fourthjetPFlav        , "fourthjetPFlav/D");
  tree->Branch("fourthjetQGL"          , &fourthjetQGL          , "fourthjetQGL/D");
  tree->Branch("fourthjetPUID"         , &fourthjetPUID         , "fourthjetPUID/D");
  tree->Branch("fourthjetGenpt"        , &fourthjetGenpt        , "fourthjetGenpt/D");
  tree->Branch("fourthjetGeneta"       , &fourthjetGeneta       , "fourthjetGeneta/D");
  tree->Branch("fourthjetGenphi"       , &fourthjetGenphi       , "fourthjetGenphi/D");
  tree->Branch("fourthjetGenm"         , &fourthjetGenm         , "fourthjetGenm/D");
  tree->Branch("fourthjetRawpt"        , &fourthjetRawpt        , "fourthjetRawpt/D");


  tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");

  tree->Branch("jetmetdphimin"        , &jetmetdphimin        , "jetmetdphimin/D");
  tree->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/D");
  tree->Branch("jetmumetdphimin"      , &jetmumetdphimin      , "jetmumetdphimin/D");
  tree->Branch("incjetmumetdphimin"   , &incjetmumetdphimin   , "incjetmumetdphimin/D");
  tree->Branch("jetelmetdphimin"      , &jetelmetdphimin      , "jetelmetdphimin/D");
  tree->Branch("incjetelmetdphimin"   , &incjetelmetdphimin   , "incjetelmetdphimin/D");
  tree->Branch("jetphmetdphimin"      , &jetphmetdphimin      , "jetphmetdphimin/D");
  tree->Branch("incjetphmetdphimin"   , &incjetphmetdphimin   , "incjetphmetdphimin/D");

  tree->Branch("jetmetdphimin4"       , &jetmetdphimin4       , "jetmetdphimin4/D");
  tree->Branch("incjetmetdphimin4"    , &incjetmetdphimin4    , "incjetmetdphimin4/D");
  tree->Branch("jetmumetdphimin4"     , &jetmumetdphimin4     , "jetmumetdphimin4/D");
  tree->Branch("incjetmumetdphimin4"  , &incjetmumetdphimin4  , "incjetmumetdphimin4/D");
  tree->Branch("jetelmetdphimin4"     , &jetelmetdphimin4     , "jetelmetdphimin4/D");
  tree->Branch("incjetelmetdphimin4"  , &incjetelmetdphimin4  , "incjetelmetdphimin4/D");
  tree->Branch("jetphmetdphimin4"     , &jetphmetdphimin4     , "jetphmetdphimin4/D");
  tree->Branch("incjetphmetdphimin4"  , &incjetphmetdphimin4  , "incjetphmetdphimin4/D");

  tree->Branch("ht"                   , &ht                   , "ht/D");

  if(addPuppiJets){

    tree->Branch("leadingPuppijetpt"         , &leadingPuppijetpt         , "leadingPuppijetpt/D");
    tree->Branch("leadingPuppijeteta"        , &leadingPuppijeteta        , "leadingPuppijeteta/D");
    tree->Branch("leadingPuppijetphi"        , &leadingPuppijetphi        , "leadingPuppijetphi/D");
    tree->Branch("leadingPuppijetm"          , &leadingPuppijetm          , "leadingPuppijetm/D");
    
    tree->Branch("signalPuppijetpt"          , &signalPuppijetpt          , "signalPuppijetpt/D");
    tree->Branch("signalPuppijeteta"         , &signalPuppijeteta         , "signalPuppijeteta/D");
    tree->Branch("signalPuppijetphi"         , &signalPuppijetphi         , "signalPuppijetphi/D");
    tree->Branch("signalPuppijetbtag"        , &signalPuppijetbtag        , "signalPuppijetbtag/D");
    tree->Branch("signalPuppijetCHfrac"      , &signalPuppijetCHfrac      , "signalPuppijetCHfrac/D");
    tree->Branch("signalPuppijetNHfrac"      , &signalPuppijetNHfrac      , "signalPuppijetNHfrac/D");
    tree->Branch("signalPuppijetEMfrac"      , &signalPuppijetEMfrac      , "signalPuppijetEMfrac/D");
    tree->Branch("signalPuppijetCEMfrac"     , &signalPuppijetCEMfrac     , "signalPuppijetCEMfrac/D");
    tree->Branch("signalPuppijetmetdphi"     , &signalPuppijetmetdphi     , "signalPuppijetmetdphi/D");
    tree->Branch("signalPuppijetHFlav"       , &signalPuppijetHFlav       , "signalPuppijetHFlav/D");
    tree->Branch("signalPuppijetPFlav"       , &signalPuppijetPFlav       , "signalPuppijetPFlav/D");
    tree->Branch("signalPuppijetQGL"         , &signalPuppijetQGL         , "signalPuppijetQGL/D");
    tree->Branch("signalPuppijetPUID"        , &signalPuppijetPUID        , "signalPuppijetPUID/D");
    tree->Branch("signalPuppijetGenpt"       , &signalPuppijetGenpt       , "signalPuppijetGenpt/D");
    tree->Branch("signalPuppijetGeneta"      , &signalPuppijetGeneta      , "signalPuppijetGeneta/D");
    tree->Branch("signalPuppijetGenphi"      , &signalPuppijetGenphi      , "signalPuppijetGenphi/D");
    tree->Branch("signalPuppijetGenm"        , &signalPuppijetGenm        , "signalPuppijetGenm/D");
    tree->Branch("signalPuppijetRawpt"       , &signalPuppijetRawpt       , "signalPuppijetRawpt/D");
    
    tree->Branch("secondPuppijetpt"          , &secondPuppijetpt          , "secondPuppijetpt/D");
    tree->Branch("secondPuppijeteta"         , &secondPuppijeteta         , "secondPuppijeteta/D");
    tree->Branch("secondPuppijetphi"         , &secondPuppijetphi         , "secondPuppijetphi/D");
    tree->Branch("secondPuppijetbtag"        , &secondPuppijetbtag        , "secondPuppijetbtag/D");
    tree->Branch("secondPuppijetCHfrac"      , &secondPuppijetCHfrac      , "secondPuppijetCHfrac/D");
    tree->Branch("secondPuppijetNHfrac"      , &secondPuppijetNHfrac      , "secondPuppijetNHfrac/D");
    tree->Branch("secondPuppijetEMfrac"      , &secondPuppijetEMfrac      , "secondPuppijetEMfrac/D");
    tree->Branch("secondPuppijetCEMfrac"     , &secondPuppijetCEMfrac     , "secondPuppijetCEMfrac/D");
    tree->Branch("secondPuppijetmetdphi"     , &secondPuppijetmetdphi     , "secondPuppijetmetdphi/D");
    tree->Branch("secondPuppijetHFlav"       , &secondPuppijetHFlav       , "secondPuppijetHFlav/D");
    tree->Branch("secondPuppijetPFlav"       , &secondPuppijetPFlav       , "secondPuppijetPFlav/D");
    tree->Branch("secondPuppijetQGL"         , &secondPuppijetQGL         , "secondPuppijetQGL/D");
    tree->Branch("secondPuppijetPUID"        , &secondPuppijetPUID        , "secondPuppijetPUID/D");
    tree->Branch("secondPuppijetGenpt"       , &secondPuppijetGenpt       , "secondPuppijetGenpt/D");
    tree->Branch("secondPuppijetGeneta"      , &secondPuppijetGeneta      , "secondPuppijetGeneta/D");
    tree->Branch("secondPuppijetGenphi"      , &secondPuppijetGenphi      , "secondPuppijetGenphi/D");
    tree->Branch("secondPuppijetGenm"        , &secondPuppijetGenm        , "secondPuppijetGenm/D");
    tree->Branch("secondPuppijetRawpt"       , &secondPuppijetRawpt       , "secondPuppijetRawpt/D");

    tree->Branch("thirdPuppijetpt"           , &thirdPuppijetpt           , "thirdPuppijetpt/D");
    tree->Branch("thirdPuppijeteta"          , &thirdPuppijeteta          , "thirdPuppijeteta/D");
    tree->Branch("thirdPuppijetphi"          , &thirdPuppijetphi          , "thirdPuppijetphi/D");
    tree->Branch("thirdPuppijetbtag"         , &thirdPuppijetbtag         , "thirdPuppijetbtag/D");
    tree->Branch("thirdPuppijetCHfrac"       , &thirdPuppijetCHfrac       , "thirdPuppijetCHfrac/D");
    tree->Branch("thirdPuppijetNHfrac"       , &thirdPuppijetNHfrac       , "thirdPuppijetNHfrac/D");
    tree->Branch("thirdPuppijetEMfrac"       , &thirdPuppijetEMfrac       , "thirdPuppijetEMfrac/D");
    tree->Branch("thirdPuppijetCEMfrac"      , &thirdPuppijetCEMfrac      , "thirdPuppijetCEMfrac/D");
    tree->Branch("thirdPuppijetmetdphi"      , &thirdPuppijetmetdphi      , "thirdPuppijetmetdphi/D");
    tree->Branch("thirdPuppijetHFlav"        , &thirdPuppijetHFlav        , "thirdPuppijetHFlav/D");
    tree->Branch("thirdPuppijetPFlav"        , &thirdPuppijetPFlav        , "thirdPuppijetPFlav/D");
    tree->Branch("thirdPuppijetQGL"          , &thirdPuppijetQGL          , "thirdPuppijetQGL/D");
    tree->Branch("thirdPuppijetPUID"         , &thirdPuppijetPUID         , "thirdPuppijetPUID/D");
    tree->Branch("thirdPuppijetGenpt"        , &thirdPuppijetGenpt        , "thirdPuppijetGenpt/D");
    tree->Branch("thirdPuppijetGeneta"       , &thirdPuppijetGeneta       , "thirdPuppijetGeneta/D");
    tree->Branch("thirdPuppijetGenphi"       , &thirdPuppijetGenphi       , "thirdPuppijetGenphi/D");
    tree->Branch("thirdPuppijetGenm"         , &thirdPuppijetGenm         , "thirdPuppijetGenm/D");
    tree->Branch("thirdPuppijetRawpt"        , &thirdPuppijetRawpt        , "thirdPuppijetRawpt/D");

    tree->Branch("fourthPuppijetpt"           , &fourthPuppijetpt        , "fourthPuppijetpt/D");
    tree->Branch("fourthPuppijeteta"          , &fourthPuppijeteta       , "fourthPuppijeteta/D");
    tree->Branch("fourthPuppijetphi"          , &fourthPuppijetphi       , "fourthPuppijetphi/D");
    tree->Branch("fourthPuppijetbtag"         , &fourthPuppijetbtag      , "fourthPuppijetbtag/D");
    tree->Branch("fourthPuppijetCHfrac"       , &fourthPuppijetCHfrac    , "fourthPuppijetCHfrac/D");
    tree->Branch("fourthPuppijetNHfrac"       , &fourthPuppijetNHfrac    , "fourthPuppijetNHfrac/D");
    tree->Branch("fourthPuppijetEMfrac"       , &fourthPuppijetEMfrac    , "fourthPuppijetEMfrac/D");
    tree->Branch("fourthPuppijetCEMfrac"      , &fourthPuppijetCEMfrac   , "fourthPuppijetCEMfrac/D");
    tree->Branch("fourthPuppijetmetdphi"      , &fourthPuppijetmetdphi   , "fourthPuppijetmetdphi/D");
    tree->Branch("fourthPuppijetHFlav"        , &fourthPuppijetHFlav     , "fourthPuppijetHFlav/D");
    tree->Branch("fourthPuppijetPFlav"        , &fourthPuppijetPFlav     , "fourthPuppijetPFlav/D");
    tree->Branch("fourthPuppijetQGL"          , &fourthPuppijetQGL       , "fourthPuppijetQGL/D");
    tree->Branch("fourthPuppijetPUID"         , &fourthPuppijetPUID      , "fourthPuppijetPUID/D");
    tree->Branch("fourthPuppijetGenpt"        , &fourthPuppijetGenpt     , "fourthPuppijetGenpt/D");
    tree->Branch("fourthPuppijetGeneta"       , &fourthPuppijetGeneta    , "fourthPuppijetGeneta/D");
    tree->Branch("fourthPuppijetGenphi"       , &fourthPuppijetGenphi    , "fourthPuppijetGenphi/D");
    tree->Branch("fourthPuppijetGenm"         , &fourthPuppijetGenm      , "fourthPuppijetGenm/D");
    tree->Branch("fourthPuppijetRawpt"        , &fourthPuppijetRawpt     , "fourthPuppijetRawpt/D");
    
    tree->Branch("PuppijetPuppijetdphi"      , &PuppijetPuppijetdphi      , "PuppijetPuppijetdphi/D");

    tree->Branch("Puppijetmetdphimin"        , &Puppijetmetdphimin        , "Puppijetmetdphimin/D");
    tree->Branch("incPuppijetmetdphimin"     , &incPuppijetmetdphimin     , "incPuppijetmetdphimin/D");

    tree->Branch("Puppijetmumetdphimin"        , &Puppijetmumetdphimin        , "Puppijetmumetdphimin/D");
    tree->Branch("incPuppijetmumetdphimin"     , &incPuppijetmumetdphimin     , "incPuppijetmumetdphimin/D");
    tree->Branch("Puppijetelmetdphimin"      , &Puppijetelmetdphimin      , "Puppijetelmetdphimin/D");
    tree->Branch("incPuppijetelmetdphimin"   , &incPuppijetelmetdphimin   , "incPuppijetelmetdphimin/D");
    tree->Branch("Puppijetphmetdphimin"      , &Puppijetphmetdphimin      , "Puppijetphmetdphimin/D");
    tree->Branch("incPuppijetphmetdphimin"   , &incPuppijetphmetdphimin   , "incPuppijetphmetdphimin/D");

    tree->Branch("Puppijetmetdphimin4"       , &Puppijetmetdphimin4       , "Puppijetmetdphimin4/D");
    tree->Branch("incPuppijetmetdphimin4"    , &incPuppijetmetdphimin4    , "incPuppijetmetdphimin4/D");
    tree->Branch("Puppijetmumetdphimin4"       , &Puppijetmumetdphimin4       , "Puppijetmumetdphimin4/D");
    tree->Branch("incPuppijetmumetdphimin4"    , &incPuppijetmumetdphimin4    , "incPuppijetmumetdphimin4/D");
    tree->Branch("Puppijetelmetdphimin4"     , &Puppijetelmetdphimin4     , "Puppijetelmetdphimin4/D");
    tree->Branch("incPuppijetelmetdphimin4"  , &incPuppijetelmetdphimin4  , "incPuppijetelmetdphimin4/D");
    tree->Branch("Puppijetphmetdphimin4"     , &Puppijetphmetdphimin4     , "Puppijetphmetdphimin4/D");
    tree->Branch("incPuppijetphmetdphimin4"  , &incPuppijetphmetdphimin4  , "incPuppijetphmetdphimin4/D");
    tree->Branch("ht"                   , &ht                   , "ht/D");
    
  }
  
  // Lepton info
  tree->Branch("mu1pid"               , &mu1pid               , "mu1pid/I");
  tree->Branch("mu1pt"                , &mu1pt                , "mu1pt/D");
  tree->Branch("mu1eta"               , &mu1eta               , "mu1eta/D");
  tree->Branch("mu1phi"               , &mu1phi               , "mu1phi/D");
  tree->Branch("mu1pfpt"              , &mu1pfpt              , "mu1pfpt/D");
  tree->Branch("mu1pfeta"             , &mu1pfeta             , "mu1pfeta/D");
  tree->Branch("mu1pfphi"             , &mu1pfphi             , "mu1pfphi/D");
  tree->Branch("mu1id"                , &mu1id                , "mu1id/I");
  tree->Branch("mu1idm"               , &mu1idm               , "mu1idm/I");
  tree->Branch("mu1idt"               , &mu1idt               , "mu1idt/I");
  tree->Branch("mu1iso"               , &mu1iso               , "mu1iso/D");
  tree->Branch("mu2pid"               , &mu2pid               , "mu2pid/I");
  tree->Branch("mu2pt"                , &mu2pt                , "mu2pt/D");
  tree->Branch("mu2eta"               , &mu2eta               , "mu2eta/D");
  tree->Branch("mu2phi"               , &mu2phi               , "mu2phi/D");
  tree->Branch("mu2pfpt"              , &mu2pfpt              , "mu2pfpt/D");
  tree->Branch("mu2pfeta"             , &mu2pfeta             , "mu2pfeta/D");
  tree->Branch("mu2pfphi"             , &mu2pfphi             , "mu2pfphi/D");
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

    // Dilepton info
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
  tree->Branch("zeept"                , &zeept                , "zeeept/D");
  tree->Branch("zeeeta"               , &zeeeta               , "zeeeta/D");
  tree->Branch("zeephi"               , &zeephi               , "zeephi/D");
  tree->Branch("wemt"                 , &wemt                 , "wemt/D");

  // Photon info
  tree->Branch("phidl"                , &phidl                , "phidl/I");
  tree->Branch("phidm"                , &phidm                , "phidm/I");
  tree->Branch("phidt"                , &phidt                , "phidt/I");
  tree->Branch("phidh"                , &phidh                , "phidh/I");
  tree->Branch("phpt"                 , &phpt                 , "phpt/D");
  tree->Branch("pheta"                , &pheta                , "pheta/D");
  tree->Branch("phphi"                , &phphi                , "phphi/D");
  
  // W/Z gen-level info
  tree->Branch("wzid"                 , &wzid                 , "wzid/I");
  tree->Branch("wzmass"               , &wzmass               , "wzmass/D");
  tree->Branch("wzmt"                 , &wzmt                 , "wzmt/D");
  tree->Branch("wzpt"                 , &wzpt                 , "wzpt/D");
  tree->Branch("wzeta"                , &wzeta                , "wzeta/D");
  tree->Branch("wzphi"                , &wzphi                , "wzphi/D");
  tree->Branch("l1id"                 , &l1id                 , "l1id/I");
  tree->Branch("l1pt"                 , &l1pt                 , "l1pt/D");
  tree->Branch("l1eta"                , &l1eta                , "l1eta/D");
  tree->Branch("l1phi"                , &l1phi                , "l1phi/D");
  tree->Branch("l2id"                 , &l2id                 , "l2id/I");
  tree->Branch("l2pt"                 , &l2pt                 , "l2pt/D");
  tree->Branch("l2eta"                , &l2eta                , "l2eta/D");
  tree->Branch("l2phi"                , &l2phi                , "l2phi/D");
  tree->Branch("parid"                , &parid                , "parid/I");
  tree->Branch("parpt"                , &parpt                , "parpt/D");
  tree->Branch("pareta"               , &pareta               , "pareta/D");
  tree->Branch("parphi"               , &parphi               , "parphi/D");
  tree->Branch("ancid"                , &ancid                , "ancid/I");
  tree->Branch("ancpt"                , &ancpt                , "ancpt/D");
  tree->Branch("anceta"               , &anceta               , "anceta/D");
  tree->Branch("ancphi"               , &ancphi               , "ancphi/D");

  // AK8 Puppi jets                                                                                                                                                             
  if(addSubstructureCHS){

    tree->Branch("leadBoostedJetpt", &leadBoostedJetpt, "leadBoostedJetpt/D");
    tree->Branch("leadBoostedJeteta", &leadBoostedJeteta, "leadBoostedJeteta/D");
    tree->Branch("leadBoostedJetphi", &leadBoostedJetphi, "leadBoostedJetphi/D");
    tree->Branch("leadBoostedJetm", &leadBoostedJetm, "leadBoostedJetm/D");
    tree->Branch("leadBoostedJetGenpt", &leadBoostedJetGenpt, "leadBoostedJetGenpt/D");
    tree->Branch("leadBoostedJetGenm", &leadBoostedJetGenm, "leadBoostedJetGenm/D");
    tree->Branch("leadBoostedJetGeneta", &leadBoostedJetGeneta, "leadBoostedJetGeneta/D");
    tree->Branch("leadBoostedJetGenphi", &leadBoostedJetGenphi, "leadBoostedJetGenphi/D");
    tree->Branch("leadBoostedJetHFlav", &leadBoostedJetHFlav, "leadBoostedJetHFlav/D");
    tree->Branch("leadBoostedJetPFlav", &leadBoostedJetPFlav, "leadBoostedJetPFlav/D");
    tree->Branch("leadBoostedJetQGL", &leadBoostedJetQGL, "leadBoostedJetQGL/D");
    tree->Branch("leadBoostedJetBtag", &leadBoostedJetBtag, "leadBoostedJetBtag/D");
    tree->Branch("leadBoostedJetDoubleBtag", &leadBoostedJetDoubleBtag, "leadBoostedJetDoubleBtag/D");

    tree->Branch("leadBoostedJettau1", &leadBoostedJettau1, "leadBoostedJettau1/D");
    tree->Branch("leadBoostedJettau2", &leadBoostedJettau2, "leadBoostedJettau2/D");
    tree->Branch("leadBoostedJettau3", &leadBoostedJettau3, "leadBoostedJettau3/D");
    tree->Branch("leadBoostedJettau4", &leadBoostedJettau4, "leadBoostedJettau4/D");

    tree->Branch("leadBoostedJetecf1", &leadBoostedJetecf1, "leadBoostedJetecf1/D");
    tree->Branch("leadBoostedJetecf2", &leadBoostedJetecf2, "leadBoostedJetecf2/D");
    tree->Branch("leadBoostedJetecf3", &leadBoostedJetecf3, "leadBoostedJetecf3/D");

    tree->Branch("leadBoostedJetBosonpt", &leadBoostedJetBosonpt, "leadBoostedJetBosonpt/D");
    tree->Branch("leadBoostedJetBosoneta", &leadBoostedJetBosoneta, "leadBoostedJetBosoneta/D");
    tree->Branch("leadBoostedJetBosonphi", &leadBoostedJetBosonphi, "leadBoostedJetBosonphi/D");
    tree->Branch("leadBoostedJetBosonm", &leadBoostedJetBosonm, "leadBoostedJetBosonm/D");

    tree->Branch("leadPrunedJetpt", &leadPrunedJetpt, "leadPrunedJetpt/D");
    tree->Branch("leadPrunedJetm", &leadPrunedJetm, "leadPrunedJetm/D");
    tree->Branch("leadPrunedJeteta", &leadPrunedJeteta, "leadPrunedJeteta/D");
    tree->Branch("leadPrunedJetphi", &leadPrunedJetphi, "leadPrunedJetphi/D");
    tree->Branch("leadPrunedJetptraw", &leadPrunedJetptraw, "leadPrunedJetptraw/D");
    tree->Branch("leadPrunedJetmraw", &leadPrunedJetmraw, "leadPrunedJetmraw/D");

    tree->Branch("leadPrunedJetGenpt", &leadPrunedJetGenpt, "leadPrunedJetGenpt/D");
    tree->Branch("leadPrunedJetGenm", &leadPrunedJetGenm, "leadPrunedJetGenm/D");
    tree->Branch("leadPrunedJetGeneta", &leadPrunedJetGeneta, "leadPrunedJetGeneta/D");
    tree->Branch("leadPrunedJetGenphi", &leadPrunedJetGenphi, "leadPrunedJetGenphi/D");
    tree->Branch("leadPrunedJetHFlav", &leadPrunedJetHFlav, "leadPrunedJetHFlav/D");
    tree->Branch("leadPrunedJetPFlav", &leadPrunedJetPFlav, "leadPrunedJetPFlav/D");
    tree->Branch("leadPrunedJetQGL", &leadPrunedJetQGL, "leadPrunedJetQGL/D");
    tree->Branch("leadPrunedJetBtag", &leadPrunedJetBtag, "leadPrunedJetBtag/D");
    tree->Branch("leadPrunedJetDoubleBtag", &leadPrunedJetDoubleBtag, "leadPrunedJetDoubleBtag/D");

    tree->Branch("leadSoftDropJetpt", &leadSoftDropJetpt, "leadSoftDropJetpt/D");
    tree->Branch("leadSoftDropJetm", &leadSoftDropJetm, "leadSoftDropJetm/D");
    tree->Branch("leadSoftDropJeteta", &leadSoftDropJeteta, "leadSoftDropJeteta/D");
    tree->Branch("leadSoftDropJetphi", &leadSoftDropJetphi, "leadSoftDropJetphi/D");
    tree->Branch("leadSoftDropJetptraw", &leadSoftDropJetptraw, "leadSoftDropJetptraw/D");
    tree->Branch("leadSoftDropJetmraw", &leadSoftDropJetmraw, "leadSoftDropJetmraw/D");

    tree->Branch("leadSoftDropJetGenpt", &leadSoftDropJetGenpt, "leadSoftDropJetGenpt/D");
    tree->Branch("leadSoftDropJetGenm", &leadSoftDropJetGenm, "leadSoftDropJetGenm/D");
    tree->Branch("leadSoftDropJetGenphi", &leadSoftDropJetGenphi, "leadSoftDropJetGenphi/D");
    tree->Branch("leadSoftDropJetGeneta", &leadSoftDropJetGeneta, "leadSoftDropJetGeneta/D");
    tree->Branch("leadSoftDropJetHFlav", &leadSoftDropJetHFlav, "leadSoftDropJetHFlav/D");
    tree->Branch("leadSoftDropJetPFlav", &leadSoftDropJetPFlav, "leadSoftDropJetPFlav/D");
    tree->Branch("leadSoftDropJetQGL", &leadSoftDropJetQGL, "leadSoftDropJetQGL/D");
    tree->Branch("leadSoftDropJetBtag", &leadSoftDropJetBtag, "leadSoftDropJetBtag/D");
    tree->Branch("leadSoftDropJetDoubleBtag", &leadSoftDropJetDoubleBtag, "leadSoftDropJetDoubleBtag/D");

    tree->Branch("leadPrunedSubJetpt_1", &leadPrunedSubJetpt_1, "leadPrunedSubJetpt_1/D");
    tree->Branch("leadPrunedSubJeteta_1", &leadPrunedSubJeteta_1, "leadPrunedSubJeteta_1/D");
    tree->Branch("leadPrunedSubJetphi_1", &leadPrunedSubJetphi_1, "leadPrunedSubJetphi_1/D");
    tree->Branch("leadPrunedSubJetm_1", &leadPrunedSubJetm_1, "leadPrunedSubJetm_1/D");
    tree->Branch("leadPrunedSubJetGenpt_1", &leadPrunedSubJetGenpt_1, "leadPrunedSubJetGenpt_1/D");
    tree->Branch("leadPrunedSubJetGenm_1", &leadPrunedSubJetGenm_1, "leadPrunedSubJetGenm_1/D");
    tree->Branch("leadPrunedSubJetGeneta_1", &leadPrunedSubJetGeneta_1, "leadPrunedSubJetGeneta_1/D");
    tree->Branch("leadPrunedSubJetGenphi_1", &leadPrunedSubJetGenphi_1, "leadPrunedSubJetGenphi_1/D");
    tree->Branch("leadPrunedSubJetHFlav_1", &leadPrunedSubJetHFlav_1, "leadPrunedSubJetHFlav_1/D");
    tree->Branch("leadPrunedSubJetPFlav_1", &leadPrunedSubJetPFlav_1, "leadPrunedSubJetPFlav_1/D");
    tree->Branch("leadPrunedSubJetQGL_1", &leadPrunedSubJetQGL_1, "leadPrunedSubJetQGL_1/D");
    tree->Branch("leadPrunedSubJetBtag_1", &leadPrunedSubJetBtag_1, "leadPrunedSubJetBtag_1/D");
    tree->Branch("leadPrunedSubJetptraw_1", &leadPrunedSubJetptraw_1, "leadPrunedSubJetptraw_1/D");
    tree->Branch("leadPrunedSubJetmraw_1", &leadPrunedSubJetmraw_1, "leadPrunedSubJetmraw_1/D");

    tree->Branch("leadPrunedSubJetpt_2", &leadPrunedSubJetpt_2, "leadPrunedSubJetpt_2/D");
    tree->Branch("leadPrunedSubJeteta_2", &leadPrunedSubJeteta_2, "leadPrunedSubJeteta_2/D");
    tree->Branch("leadPrunedSubJetphi_2", &leadPrunedSubJetphi_2, "leadPrunedSubJetphi_2/D");
    tree->Branch("leadPrunedSubJetm_2", &leadPrunedSubJetm_2, "leadPrunedSubJetm_2/D");
    tree->Branch("leadPrunedSubJetGenpt_2", &leadPrunedSubJetGenpt_2, "leadPrunedSubJetGenpt_2/D");
    tree->Branch("leadPrunedSubJetGenm_2", &leadPrunedSubJetGenm_2, "leadPrunedSubJetGenm_2/D");
    tree->Branch("leadPrunedSubJetGenphi_2", &leadPrunedSubJetGenphi_2, "leadPrunedSubJetGenphi_2/D");
    tree->Branch("leadPrunedSubJetGeneta_2", &leadPrunedSubJetGeneta_2, "leadPrunedSubJetGeneta_2/D");
    tree->Branch("leadPrunedSubJetHFlav_2", &leadPrunedSubJetHFlav_2, "leadPrunedSubJetHFlav_2/D");
    tree->Branch("leadPrunedSubJetPFlav_2", &leadPrunedSubJetPFlav_2, "leadPrunedSubJetPFlav_2/D");
    tree->Branch("leadPrunedSubJetQGL_2", &leadPrunedSubJetQGL_2, "leadPrunedSubJetQGL_2/D");
    tree->Branch("leadPrunedSubJetBtag_2", &leadPrunedSubJetBtag_2, "leadPrunedSubJetBtag_2/D");
    tree->Branch("leadPrunedSubJetptraw_2", &leadPrunedSubJetptraw_2, "leadPrunedSubJetptraw_2/D");
    tree->Branch("leadPrunedSubJetmraw_2", &leadPrunedSubJetmraw_2, "leadPrunedSubJetmraw_2/D");
 
    tree->Branch("leadSoftDropSubJetpt_1", &leadSoftDropSubJetpt_1, "leadSoftDropSubJetpt_1/D");
    tree->Branch("leadSoftDropSubJeteta_1", &leadSoftDropSubJeteta_1, "leadSoftDropSubJeteta_1/D");
    tree->Branch("leadSoftDropSubJetphi_1", &leadSoftDropSubJetphi_1, "leadSoftDropSubJetphi_1/D");
    tree->Branch("leadSoftDropSubJetm_1", &leadSoftDropSubJetm_1, "leadSoftDropSubJetm_1/D");
    tree->Branch("leadSoftDropSubJetGenpt_1", &leadSoftDropSubJetGenpt_1, "leadSoftDropSubJetGenpt_1/D");
    tree->Branch("leadSoftDropSubJetGenm_1", &leadSoftDropSubJetGenm_1, "leadSoftDropSubJetGenm_1/D");
    tree->Branch("leadSoftDropSubJetGeneta_1", &leadSoftDropSubJetGeneta_1, "leadSoftDropSubJetGeneta_1/D");
    tree->Branch("leadSoftDropSubJetGenphi_1", &leadSoftDropSubJetGenphi_1, "leadSoftDropSubJetGenphi_1/D");
    tree->Branch("leadSoftDropSubJetHFlav_1", &leadSoftDropSubJetHFlav_1, "leadSoftDropSubJetHFlav_1/D");
    tree->Branch("leadSoftDropSubJetPFlav_1", &leadSoftDropSubJetPFlav_1, "leadSoftDropSubJetPFlav_1/D");
    tree->Branch("leadSoftDropSubJetQGL_1", &leadSoftDropSubJetQGL_1, "leadSoftDropSubJetQGL_1/D");
    tree->Branch("leadSoftDropSubJetBtag_1", &leadSoftDropSubJetBtag_1, "leadSoftDropSubJetBtag_1/D");
    tree->Branch("leadSoftDropSubJetptraw_1", &leadSoftDropSubJetptraw_1, "leadSoftDropSubJetptraw_1/D");
    tree->Branch("leadSoftDropSubJetmraw_1", &leadSoftDropSubJetmraw_1, "leadSoftDropSubJetmraw_1/D");

    tree->Branch("leadSoftDropSubJetpt_2", &leadSoftDropSubJetpt_2, "leadSoftDropSubJetpt_2/D");
    tree->Branch("leadSoftDropSubJeteta_2", &leadSoftDropSubJeteta_2, "leadSoftDropSubJeteta_2/D");
    tree->Branch("leadSoftDropSubJetphi_2", &leadSoftDropSubJetphi_2, "leadSoftDropSubJetphi_2/D");
    tree->Branch("leadSoftDropSubJetm_2", &leadSoftDropSubJetm_2, "leadSoftDropSubJetm_2/D");
    tree->Branch("leadSoftDropSubJetGenpt_2", &leadSoftDropSubJetGenpt_2, "leadSoftDropSubJetGenpt_2/D");
    tree->Branch("leadSoftDropSubJetGenm_2", &leadSoftDropSubJetGenm_2, "leadSoftDropSubJetGenm_2/D");
    tree->Branch("leadSoftDropSubJetGenphi_2", &leadSoftDropSubJetGenphi_2, "leadSoftDropSubJetGenphi_2/D");
    tree->Branch("leadSoftDropSubJetGeneta_2", &leadSoftDropSubJetGeneta_2, "leadSoftDropSubJetGeneta_2/D");
    tree->Branch("leadSoftDropSubJetHFlav_2", &leadSoftDropSubJetHFlav_2, "leadSoftDropSubJetHFlav_2/D");
    tree->Branch("leadSoftDropSubJetPFlav_2", &leadSoftDropSubJetPFlav_2, "leadSoftDropSubJetPFlav_2/D");
    tree->Branch("leadSoftDropSubJetQGL_2", &leadSoftDropSubJetQGL_2, "leadSoftDropSubJetQGL_2/D");
    tree->Branch("leadSoftDropSubJetBtag_2", &leadSoftDropSubJetBtag_2, "leadSoftDropSubJetBtag_2/D");
    tree->Branch("leadSoftDropSubJetptraw_2", &leadSoftDropSubJetptraw_2, "leadSoftDropSubJetptraw_2/D");
    tree->Branch("leadSoftDropSubJetmraw_2", &leadSoftDropSubJetmraw_2, "leadSoftDropSubJetmraw_2/D");

  }

  if(addSubstructurePuppi){

    tree->Branch("leadPuppiBoostedJetpt", &leadPuppiBoostedJetpt, "leadPuppiBoostedJetpt/D");
    tree->Branch("leadPuppiBoostedJeteta", &leadPuppiBoostedJeteta, "leadPuppiBoostedJeteta/D");
    tree->Branch("leadPuppiBoostedJetphi", &leadPuppiBoostedJetphi, "leadPuppiBoostedJetphi/D");
    tree->Branch("leadPuppiBoostedJetm", &leadPuppiBoostedJetm, "leadPuppiBoostedJetm/D");
    tree->Branch("leadPuppiBoostedJetGenpt", &leadPuppiBoostedJetGenpt, "leadPuppiBoostedJetGenpt/D");
    tree->Branch("leadPuppiBoostedJetGenm", &leadPuppiBoostedJetGenm, "leadPuppiBoostedJetGenm/D");
    tree->Branch("leadPuppiBoostedJetGeneta", &leadPuppiBoostedJetGeneta, "leadPuppiBoostedJetGeneta/D");
    tree->Branch("leadPuppiBoostedJetGenphi", &leadPuppiBoostedJetGenphi, "leadPuppiBoostedJetGenphi/D");
    tree->Branch("leadPuppiBoostedJetHFlav", &leadPuppiBoostedJetHFlav, "leadPuppiBoostedJetHFlav/D");
    tree->Branch("leadPuppiBoostedJetPFlav", &leadPuppiBoostedJetPFlav, "leadPuppiBoostedJetPFlav/D");
    tree->Branch("leadPuppiBoostedJetQGL", &leadPuppiBoostedJetQGL, "leadPuppiBoostedJetQGL/D");
    tree->Branch("leadPuppiBoostedJetBtag", &leadPuppiBoostedJetBtag, "leadPuppiBoostedJetBtag/D");
    tree->Branch("leadPuppiBoostedJetBtag", &leadPuppiBoostedJetDoubleBtag, "leadPuppiBoostedJetDoubleBtag/D");

    tree->Branch("leadPuppiBoostedJettau1", &leadPuppiBoostedJettau1, "leadPuppiBoostedJettau1/D");
    tree->Branch("leadPuppiBoostedJettau2", &leadPuppiBoostedJettau2, "leadPuppiBoostedJettau2/D");
    tree->Branch("leadPuppiBoostedJettau3", &leadPuppiBoostedJettau3, "leadPuppiBoostedJettau3/D");
    tree->Branch("leadPuppiBoostedJettau4", &leadPuppiBoostedJettau4, "leadPuppiBoostedJettau4/D");

    tree->Branch("leadPuppiBoostedJetecf1", &leadPuppiBoostedJetecf1, "leadPuppiBoostedJetecf1/D");
    tree->Branch("leadPuppiBoostedJetecf2", &leadPuppiBoostedJetecf2, "leadPuppiBoostedJetecf2/D");
    tree->Branch("leadPuppiBoostedJetecf3", &leadPuppiBoostedJetecf3, "leadPuppiBoostedJetecf3/D");

    tree->Branch("leadPuppiBoostedJetBosonpt", &leadPuppiBoostedJetBosonpt, "leadPuppiBoostedJetBosonpt/D");
    tree->Branch("leadPuppiBoostedJetBosoneta", &leadPuppiBoostedJetBosoneta, "leadPuppiBoostedJetBosoneta/D");
    tree->Branch("leadPuppiBoostedJetBosonphi", &leadPuppiBoostedJetBosonphi, "leadPuppiBoostedJetBosonphi/D");
    tree->Branch("leadPuppiBoostedJetBosonm", &leadPuppiBoostedJetBosonm, "leadPuppiBoostedJetBosonm/D");

    tree->Branch("leadPuppiPrunedJetpt", &leadPuppiPrunedJetpt, "leadPuppiPrunedJetpt/D");
    tree->Branch("leadPuppiPrunedJetm", &leadPuppiPrunedJetm, "leadPuppiPrunedJetm/D");
    tree->Branch("leadPuppiPrunedJeteta", &leadPuppiPrunedJeteta, "leadPuppiPrunedJeteta/D");
    tree->Branch("leadPuppiPrunedJetphi", &leadPuppiPrunedJetphi, "leadPuppiPrunedJetphi/D");
    tree->Branch("leadPuppiPrunedJetGenpt", &leadPuppiPrunedJetGenpt, "leadPuppiPrunedJetGenpt/D");
    tree->Branch("leadPuppiPrunedJetGenm", &leadPuppiPrunedJetGenm, "leadPuppiPrunedJetGenm/D");
    tree->Branch("leadPuppiPrunedJetGenphi", &leadPuppiPrunedJetGenphi, "leadPuppiPrunedJetGenphi/D");
    tree->Branch("leadPuppiPrunedJetGeneta", &leadPuppiPrunedJetGeneta, "leadPuppiPrunedJetGeneta/D");
    tree->Branch("leadPuppiPrunedJetHFlav", &leadPuppiPrunedJetHFlav, "leadPuppiPrunedJetHFlav/D");
    tree->Branch("leadPuppiPrunedJetPFlav", &leadPuppiPrunedJetPFlav, "leadPuppiPrunedJetPFlav/D");
    tree->Branch("leadPuppiPrunedJetQGL", &leadPuppiPrunedJetQGL, "leadPuppiPrunedJetQGL/D");
    tree->Branch("leadPuppiPrunedJetBtag", &leadPuppiPrunedJetBtag, "leadPuppiPrunedJetBtag/D");
    tree->Branch("leadPuppiPrunedJetDoubleBtag", &leadPuppiPrunedJetDoubleBtag, "leadPuppiPrunedJetDoubleBtag/D");

    tree->Branch("leadPuppiSoftDropJetpt", &leadPuppiSoftDropJetpt, "leadPuppiSoftDropJetpt/D");
    tree->Branch("leadPuppiSoftDropJetm", &leadPuppiSoftDropJetm, "leadPuppiSoftDropJetm/D");
    tree->Branch("leadPuppiSoftDropJeteta", &leadPuppiSoftDropJeteta, "leadPuppiSoftDropJeteta/D");
    tree->Branch("leadPuppiSoftDropJetphi", &leadPuppiSoftDropJetphi, "leadPuppiSoftDropJetphi/D");
    tree->Branch("leadPuppiSoftDropJetGenpt", &leadPuppiSoftDropJetGenpt, "leadPuppiSoftDropJetGenpt/D");
    tree->Branch("leadPuppiSoftDropJetGenm", &leadPuppiSoftDropJetGenm, "leadPuppiSoftDropJetGenm/D");
    tree->Branch("leadPuppiSoftDropJetGeneta", &leadPuppiSoftDropJetGeneta, "leadPuppiSoftDropJetGeneta/D");
    tree->Branch("leadPuppiSoftDropJetGenphi", &leadPuppiSoftDropJetGenphi, "leadPuppiSoftDropJetGenphi/D");
    tree->Branch("leadPuppiSoftDropJetHFlav", &leadPuppiSoftDropJetHFlav, "leadPuppiSoftDropJetHFlav/D");
    tree->Branch("leadPuppiSoftDropJetPFlav", &leadPuppiSoftDropJetPFlav, "leadPuppiSoftDropJetPFlav/D");
    tree->Branch("leadPuppiSoftDropJetQGL", &leadPuppiSoftDropJetQGL, "leadPuppiSoftDropJetQGL/D");
    tree->Branch("leadPuppiSoftDropJetBtag", &leadPuppiSoftDropJetBtag, "leadPuppiSoftDropJetBtag/D");
    tree->Branch("leadPuppiSoftDropJetDoubleBtag", &leadPuppiSoftDropJetDoubleBtag, "leadPuppiSoftDropJetDoubleBtag/D");

    tree->Branch("leadPuppiPrunedSubJetpt_1", &leadPuppiPrunedSubJetpt_1, "leadPuppiPrunedSubJetpt_1/D");
    tree->Branch("leadPuppiPrunedSubJeteta_1", &leadPuppiPrunedSubJeteta_1, "leadPuppiPrunedSubJeteta_1/D");
    tree->Branch("leadPuppiPrunedSubJetphi_1", &leadPuppiPrunedSubJetphi_1, "leadPuppiPrunedSubJetphi_1/D");
    tree->Branch("leadPuppiPrunedSubJetm_1", &leadPuppiPrunedSubJetm_1, "leadPuppiPrunedSubJetm_1/D");
    tree->Branch("leadPuppiPrunedSubJetGenpt_1", &leadPuppiPrunedSubJetGenpt_1, "leadPuppiPrunedSubJetGenpt_1/D");
    tree->Branch("leadPuppiPrunedSubJetGenm_1", &leadPuppiPrunedSubJetGenm_1, "leadPuppiPrunedSubJetGenm_1/D");
    tree->Branch("leadPuppiPrunedSubJetGenphi_1", &leadPuppiPrunedSubJetGenphi_1, "leadPuppiPrunedSubJetGenphi_1/D");
    tree->Branch("leadPuppiPrunedSubJetGeneta_1", &leadPuppiPrunedSubJetGeneta_1, "leadPuppiPrunedSubJetGeneta_1/D");
    tree->Branch("leadPuppiPrunedSubJetHFlav_1", &leadPuppiPrunedSubJetHFlav_1, "leadPuppiPrunedSubJetHFlav_1/D");
    tree->Branch("leadPuppiPrunedSubJetPFlav_1", &leadPuppiPrunedSubJetPFlav_1, "leadPuppiPrunedSubJetPFlav_1/D");
    tree->Branch("leadPuppiPrunedSubJetQGL_1", &leadPuppiPrunedSubJetQGL_1, "leadPuppiPrunedSubJetQGL_1/D");
    tree->Branch("leadPuppiPrunedSubJetBtag_1", &leadPuppiPrunedSubJetBtag_1, "leadPuppiPrunedSubJetBtag_1/D");
    tree->Branch("leadPuppiPrunedSubJetptraw_1", &leadPuppiPrunedSubJetptraw_1, "leadPuppiPrunedSubJetptraw_1/D");
    tree->Branch("leadPuppiPrunedSubJetmraw_1", &leadPuppiPrunedSubJetmraw_1, "leadPuppiPrunedSubJetmraw_1/D");

    tree->Branch("leadPuppiPrunedSubJetpt_2", &leadPuppiPrunedSubJetpt_2, "leadPuppiPrunedSubJetpt_2/D");
    tree->Branch("leadPuppiPrunedSubJeteta_2", &leadPuppiPrunedSubJeteta_2, "leadPuppiPrunedSubJeteta_2/D");
    tree->Branch("leadPuppiPrunedSubJetphi_2", &leadPuppiPrunedSubJetphi_2, "leadPuppiPrunedSubJetphi_2/D");
    tree->Branch("leadPuppiPrunedSubJetm_2", &leadPuppiPrunedSubJetm_2, "leadPuppiPrunedSubJetm_2/D");
    tree->Branch("leadPuppiPrunedSubJetGenpt_2", &leadPuppiPrunedSubJetGenpt_2, "leadPuppiPrunedSubJetGenpt_2/D");
    tree->Branch("leadPuppiPrunedSubJetGenm_2", &leadPuppiPrunedSubJetGenm_2, "leadPuppiPrunedSubJetGenm_2/D");
    tree->Branch("leadPuppiPrunedSubJetGeneta_2", &leadPuppiPrunedSubJetGeneta_2, "leadPuppiPrunedSubJetGeneta_2/D");
    tree->Branch("leadPuppiPrunedSubJetGenphi_2", &leadPuppiPrunedSubJetGenphi_2, "leadPuppiPrunedSubJetGenphi_2/D");
    tree->Branch("leadPuppiPrunedSubJetHFlav_2", &leadPuppiPrunedSubJetHFlav_2, "leadPuppiPrunedSubJetHFlav_2/D");
    tree->Branch("leadPuppiPrunedSubJetPFlav_2", &leadPuppiPrunedSubJetPFlav_2, "leadPuppiPrunedSubJetPFlav_2/D");
    tree->Branch("leadPuppiPrunedSubJetQGL_2", &leadPuppiPrunedSubJetQGL_2, "leadPuppiPrunedSubJetQGL_2/D");
    tree->Branch("leadPuppiPrunedSubJetBtag_2", &leadPuppiPrunedSubJetBtag_2, "leadPuppiPrunedSubJetBtag_2/D");
    tree->Branch("leadPuppiPrunedSubJetptraw_2", &leadPuppiPrunedSubJetptraw_2, "leadPuppiPrunedSubJetptraw_2/D");
    tree->Branch("leadPuppiPrunedSubJetmraw_2", &leadPuppiPrunedSubJetmraw_2, "leadPuppiPrunedSubJetmraw_2/D");
 
    tree->Branch("leadPuppiSoftDropSubJetpt_1", &leadPuppiSoftDropSubJetpt_1, "leadPuppiSoftDropSubJetpt_1/D");
    tree->Branch("leadPuppiSoftDropSubJeteta_1", &leadPuppiSoftDropSubJeteta_1, "leadPuppiSoftDropSubJeteta_1/D");
    tree->Branch("leadPuppiSoftDropSubJetphi_1", &leadPuppiSoftDropSubJetphi_1, "leadPuppiSoftDropSubJetphi_1/D");
    tree->Branch("leadPuppiSoftDropSubJetm_1", &leadPuppiSoftDropSubJetm_1, "leadPuppiSoftDropSubJetm_1/D");
    tree->Branch("leadPuppiSoftDropSubJetGenpt_1", &leadPuppiSoftDropSubJetGenpt_1, "leadPuppiSoftDropSubJetGenpt_1/D");
    tree->Branch("leadPuppiSoftDropSubJetGenm_1", &leadPuppiSoftDropSubJetGenm_1, "leadPuppiSoftDropSubJetGenm_1/D");
    tree->Branch("leadPuppiSoftDropSubJetGeneta_1", &leadPuppiSoftDropSubJetGeneta_1, "leadPuppiSoftDropSubJetGeneta_1/D");
    tree->Branch("leadPuppiSoftDropSubJetGenphi_1", &leadPuppiSoftDropSubJetGenphi_1, "leadPuppiSoftDropSubJetGenphi_1/D");
    tree->Branch("leadPuppiSoftDropSubJetHFlav_1", &leadPuppiSoftDropSubJetHFlav_1, "leadPuppiSoftDropSubJetHFlav_1/D");
    tree->Branch("leadPuppiSoftDropSubJetPFlav_1", &leadPuppiSoftDropSubJetPFlav_1, "leadPuppiSoftDropSubJetPFlav_1/D");
    tree->Branch("leadPuppiSoftDropSubJetQGL_1", &leadPuppiSoftDropSubJetQGL_1, "leadPuppiSoftDropSubJetQGL_1/D");
    tree->Branch("leadPuppiSoftDropSubJetBtag_1", &leadPuppiSoftDropSubJetBtag_1, "leadPuppiSoftDropSubJetBtag_1/D");
    tree->Branch("leadPuppiSoftDropSubJetptraw_1", &leadPuppiSoftDropSubJetptraw_1, "leadPuppiSoftDropSubJetptraw_1/D");
    tree->Branch("leadPuppiSoftDropSubJetmraw_1", &leadPuppiSoftDropSubJetmraw_1, "leadPuppiSoftDropSubJetmraw_1/D");

    tree->Branch("leadPuppiSoftDropSubJetpt_2", &leadPuppiSoftDropSubJetpt_2, "leadPuppiSoftDropSubJetpt_2/D");
    tree->Branch("leadPuppiSoftDropSubJeteta_2", &leadPuppiSoftDropSubJeteta_2, "leadPuppiSoftDropSubJeteta_2/D");
    tree->Branch("leadPuppiSoftDropSubJetphi_2", &leadPuppiSoftDropSubJetphi_2, "leadPuppiSoftDropSubJetphi_2/D");
    tree->Branch("leadPuppiSoftDropSubJetm_2", &leadPuppiSoftDropSubJetm_2, "leadPuppiSoftDropSubJetm_2/D");
    tree->Branch("leadPuppiSoftDropSubJetGenpt_2", &leadPuppiSoftDropSubJetGenpt_2, "leadPuppiSoftDropSubJetGenpt_2/D");
    tree->Branch("leadPuppiSoftDropSubJetGenm_2", &leadPuppiSoftDropSubJetGenm_2, "leadPuppiSoftDropSubJetGenm_2/D");
    tree->Branch("leadPuppiSoftDropSubJetGenphi_2", &leadPuppiSoftDropSubJetGenphi_2, "leadPuppiSoftDropSubJetGenphi_2/D");
    tree->Branch("leadPuppiSoftDropSubJetGeneta_2", &leadPuppiSoftDropSubJetGeneta_2, "leadPuppiSoftDropSubJetGeneta_2/D");
    tree->Branch("leadPuppiSoftDropSubJetHFlav_2", &leadPuppiSoftDropSubJetHFlav_2, "leadPuppiSoftDropSubJetHFlav_2/D");
    tree->Branch("leadPuppiSoftDropSubJetPFlav_2", &leadPuppiSoftDropSubJetPFlav_2, "leadPuppiSoftDropSubJetPFlav_2/D");
    tree->Branch("leadPuppiSoftDropSubJetQGL_2", &leadPuppiSoftDropSubJetQGL_2, "leadPuppiSoftDropSubJetQGL_2/D");
    tree->Branch("leadPuppiSoftDropSubJetBtag_2", &leadPuppiSoftDropSubJetBtag_2, "leadPuppiSoftDropSubJetBtag_2/D");
    tree->Branch("leadPuppiSoftDropSubJetptraw_2", &leadPuppiSoftDropSubJetptraw_2, "leadPuppiSoftDropSubJetptraw_2/D");
    tree->Branch("leadPuppiSoftDropSubJetmraw_2", &leadPuppiSoftDropSubJetmraw_2, "leadPuppiSoftDropSubJetmraw_2/D");
  }
}

void MonoJetTreeMaker::endJob() {}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_PFMET90_PFMHT90_IDTight");
  triggerPathsVector.push_back("HLT_PFMET120_PFMHT120_IDTight");
  triggerPathsVector.push_back("HLT_PFMET170_NoiseCleaned");
  triggerPathsVector.push_back("HLT_PFMET170_JetIdCleaned");
  triggerPathsVector.push_back("HLT_PFMET170_HBHECleaned");
  triggerPathsVector.push_back("HLT_PFMET170_v");
  triggerPathsVector.push_back("HLT_PFMET300_NoiseCleaned");
  triggerPathsVector.push_back("HLT_PFMET300_JetIdCleaned");
  triggerPathsVector.push_back("HLT_PFMET300_v");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight");
  triggerPathsVector.push_back("HLT_Photon165_HE10");
  triggerPathsVector.push_back("HLT_Photon175");
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  triggerPathsVector.push_back("HLT_IsoMu20");
  triggerPathsVector.push_back("HLT_IsoTkMu20");
  triggerPathsVector.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggerPathsVector.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  triggerPathsVector.push_back("HLT_Ele23_WPLoose_Gsf_v");
  triggerPathsVector.push_back("HLT_Ele27_WPLoose_Gsf_v");
  triggerPathsVector.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerPathsVector.push_back("HLT_Ele27_WP85_Gsf_v");
  
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
  
  // MET filter Paths
  filterPathsVector.push_back("Flag_CSCTightHaloFilter");
  filterPathsVector.push_back("Flag_HBHENoiseFilter");
  filterPathsVector.push_back("Flag_eeBadScFilter");
  
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
}

void MonoJetTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}


//This code is ripped off from https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/PhotonNtuplerMiniAOD.cc
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

double MonoJetTreeMaker::computeMuonIso(const reco::Muon& mu) {

    double isoval = mu.pfIsolationR04().sumNeutralHadronEt;
    isoval += mu.pfIsolationR04().sumPhotonEt;
    isoval -= 0.5*mu.pfIsolationR04().sumPUPt;
    if (isoval < 0.) isoval = 0.;
    isoval += mu.pfIsolationR04().sumChargedHadronPt;
    isoval /= mu.pt();            

    return isoval;
}

bool MonoJetTreeMaker::applyJetID(const pat::Jet & jet, const std::string & level){

  if(level != "loose" and level != "tight" and level != "tightLepVeto")
    return true;
   
  bool passjetid = false;

  //apply a loose jet id https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
  if(level == "loose"){
 
   if (fabs(jet.eta()) <= 3.0 && 
	jet.neutralHadronEnergyFraction() < 0.99 && 
	jet.neutralEmEnergyFraction() < 0.99 && 
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4) 
	passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 && 
	       jet.chargedHadronEnergyFraction() > 0. && 
	       jet.chargedEmEnergyFraction() < 0.99 && 
	       jet.chargedMultiplicity() > 0) 
      passjetid = true;
    }
    if (fabs(jet.eta()) > 3.0 
	&& jet.neutralEmEnergyFraction() < 0.9 
      && jet.neutralMultiplicity() > 10) 
      passjetid = true;
  }
  else if(level == "tight"){

   if (fabs(jet.eta()) <= 3.0 && 
	jet.neutralHadronEnergyFraction() < 0.90 && 
	jet.neutralEmEnergyFraction() < 0.90 && 
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4) 
	passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 && 
	       jet.chargedHadronEnergyFraction() > 0. && 
	       jet.chargedEmEnergyFraction() < 0.99 && 
	       jet.chargedMultiplicity() > 0) 
      passjetid = true;
    }
    if (fabs(jet.eta()) > 3.0 
	&& jet.neutralEmEnergyFraction() < 0.9 
      && jet.neutralMultiplicity() > 10) 
      passjetid = true;

  }

  else if(level == "tightLepVeto"){
    if (fabs(jet.eta()) <= 3.0 &&
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
    if (fabs(jet.eta()) > 3.0
        && jet.neutralEmEnergyFraction() < 0.9
	&& jet.neutralMultiplicity() > 10)
      passjetid = true;


  }

  return passjetid;

}

bool MonoJetTreeMaker::applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi){

  bool passpuid = false;
  double puidval = 0;
  double jetabseta = fabs(jet.eta());
  if(jet.hasUserFloat("puid:fullDiscriminant"))
    puidval = jet.userFloat("puid:fullDiscriminant");
  else if(jet.hasUserFloat("puidPuppi:fullDiscriminant"))
    puidval = jet.userFloat("puidPuppi:fullDiscriminant");
  else if(jet.hasUserFloat("pileupJetId:fullDiscriminant"))
    puidval = jet.userFloat("pileupJetId:fullDiscriminant");
  else 
    return true;

  // https://indico.cern.ch/event/450785/contribution/2/attachments/1167545/1683858/151008_JMAR_pileupJetIDtraining.pdf
  if(level == "loose"){
    if (jetabseta >= 0.00 && jetabseta < 2.00 && puidval > -0.82) passpuid = true;
    if (jetabseta >= 2.00 && jetabseta < 2.50 && puidval > -0.81) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 3.00 && puidval > -0.57) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.36) passpuid = true;
  }
  else if(level == "medium"){
    if (jetabseta >= 0.00 && jetabseta < 2.00 && puidval > -0.48) passpuid = true;
    if (jetabseta >= 2.00 && jetabseta < 2.50 && puidval > -0.66) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 3.00 && puidval > -0.44) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.29) passpuid = true;    
  }
  else if(level == "tight"){
    if (jetabseta >= 0.00 && jetabseta < 2.00 && puidval >  0.29) passpuid = true;
    if (jetabseta >= 2.00 && jetabseta < 2.50 && puidval > -0.30) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 3.00 && puidval > -0.37) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.25) passpuid = true;    

  }
  return passpuid;
}



void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MonoJetTreeMaker);

//  LocalWords:  TypeI
