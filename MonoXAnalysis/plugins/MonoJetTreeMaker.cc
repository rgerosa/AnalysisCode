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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// HLT info
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

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
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

// Jet corrections
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// b-tagging SF
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


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
  
  void findMother(const reco::Candidate*, int &, double &, double &, double &);
  void findFirstNonPhotonMother(const reco::Candidate*, int &, double &, double &, double &);
  double computeMuonIso(const reco::Muon&);

  bool applyJetID(const pat::Jet & jet, const std::string & level);
  bool applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi);

  bool readDMFromGenParticle;

  // Gen Particles
  const bool isMC;
  const bool uselheweights;
  const edm::InputTag lheEventTag;
  const edm::InputTag lheRunTag;
  const bool isWorZorSignalMCSample;
  
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
  edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
  edm::EDGetTokenT<LHEEventProduct>                  lheInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
  double xsec;

  // InputTags
  const edm::InputTag triggerResultsTag;
  const edm::InputTag filterResultsTag;
  const edm::InputTag prescalesTag;

  const edm::InputTag hbhelooseTag;
  const edm::InputTag hbhetightTag;
  const edm::InputTag hbheisoTag;
  
  // trgger and filter tokens
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
  edm::EDGetTokenT<edm::TriggerResults> filterResultsToken;
  edm::EDGetTokenT<bool> hbhelooseToken;
  edm::EDGetTokenT<bool> hbhetightToken;
  edm::EDGetTokenT<bool> hbheisoToken;
  
  // Vertex
  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  // muons
  const edm::InputTag muonsTag;
  const edm::InputTag tightmuonsTag;
  const edm::InputTag highptmuonsTag;

  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> tightmuonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> highptmuonsToken;

  // electrons
  const edm::InputTag  electronsTag;
  const edm::InputTag  tightelectronsTag;
  const edm::InputTag  heepelectronsTag;
  const edm::InputTag  electronLooseIdTag;  

  edm::EDGetTokenT<pat::ElectronRefVector>  electronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  tightelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  heepelectronsToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    electronLooseIdToken;  

  // Photons
  const edm::InputTag  photonsTag;
  const edm::InputTag  tightphotonsTag;
  const edm::InputTag  photonLooseIdTag;
  const edm::InputTag  photonMediumIdTag;
  const edm::InputTag  photonTightIdTag;
  const edm::InputTag  photonHighPtIdTag;

  edm::EDGetTokenT<pat::PhotonRefVector>    photonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    tightphotonsToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonLooseIdToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonMediumIdToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonTightIdToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonHighPtIdToken;

  // Taus
  const edm::InputTag tausTag;
  edm::EDGetTokenT<std::vector<pat::Tau> >  tausToken;

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

  edm::EDGetTokenT<edm::View<pat::MET> >  t1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elemetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetToken;

  // Puppi MET
  const bool addPuppiMET;

  const edm::InputTag puppit1metTag;
  const edm::InputTag puppit1mumetTag;
  const edm::InputTag puppit1elemetTag;
  const edm::InputTag puppit1phmetTag;

  edm::EDGetTokenT<edm::View<pat::MET> > puppit1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1elemetToken;
  edm::EDGetTokenT<edm::View<pat::MET> > puppit1phmetToken;

  // MET systematics
  const bool addMETSystematics;

  // MVA met
  const bool addMVAMet;

  const edm::InputTag mvaMETTag;
  edm::EDGetTokenT<edm::View<reco::MET> > mvaMETToken;

  // inner bools
  const bool applyHLTFilter;
  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;   

  // Jet AK8
  const bool addSubstructureCHS;
  const bool addSubstructurePuppi;  
  edm::EDGetTokenT<std::vector<pat::Jet> > boostedJetsToken;
  TString boostedJetsCHSLabel;
  edm::EDGetTokenT<std::vector<pat::Jet> > boostedPuppiJetsToken;
  TString boostedJetsPuppiLabel;
  
  // inner vectors
  std::vector<std::string>   triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  std::vector<std::string>   filterPathsVector;
  std::map<std::string, int> filterPathsMap;

  // Btag SF 
  bool addBTagScaleFactor;
  edm::FileInPath bTagScaleFactorFile;
  std::auto_ptr<BTagCalibration> calib;
  std::vector<BTagCalibrationReader> bMedium;
  std::vector<BTagCalibrationReader> bMediumUp;
  std::vector<BTagCalibrationReader> bMediumDown;

  // tree
  TTree* tree;

  // pileup info
  int32_t puobs, putrue; 
  int32_t wzid, l1id, l2id;
  int32_t wzid_h, q1id, q2id;
  int32_t top_1, top_2;  
  int32_t mu1pid, mu2pid, mu1id, mu2id, mu1idm, mu2idm, mu1idt, mu2idt;
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
  uint32_t njetsinc, npuppijetsinc;
  // trigger and met filters flags 
  uint8_t hltmet90,    hltmet120,    hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300;
  uint8_t hltjetmet90, hltjetmet120, hltphoton165,   hltphoton175,    hltphoton120,    hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel, hltelnoiso;
  uint8_t flagcsctight, flaghbhenoise, flaghbheloose, flaghbhetight, flaghbheiso, flageebadsc;
  // muon, ele, dilepton info
  double mu1pt,mu1eta,mu1phi,mu1pfpt,mu1pfeta, mu1pfphi, mu1iso, mu2pt, mu2eta, mu2phi, mu2pfpt, mu2pfeta, mu2pfphi, mu2iso;
  double el1pt,el1eta,el1phi,ele1e,el2pt, ele2e, el2eta, el2phi, phpt, pheta, phphi, phe;
  double zmass,zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt; 
  // PF MET info (typeI and Raw)
  double t1pfmet, t1pfmetphi, t1mumet, t1mumetphi, t1elmet, t1elmetphi, t1phmet, t1phmetphi;
  double pfmet, pfmetphi, mumet, mumetphi, elmet, elmetphi, phmet, phmetphi;
  // Puppi MET info (typeI and Raw)
  double puppipfmet, puppipfmetphi, puppimumet, puppimumetphi, puppielmet, puppielmetphi, puppiphmet, puppiphmetphi;
  double puppit1pfmet, puppit1pfmetphi, puppit1mumet, puppit1mumetphi, puppit1elmet, puppit1elmetphi, puppit1phmet, puppit1phmetphi;
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
  double leadingjetpt, leadingjeteta, leadingjetphi, leadingjetm; 

  // AK4CHS central jet
  std::vector<double> centraljetpt, centraljeteta, centraljetphi, centraljetm, centraljetbtag;
  std::vector<double> centraljetCHfrac, centraljetNHfrac, centraljetEMfrac, centraljetCEMfrac, centraljetmetdphi;
  std::vector<double> centraljetHFlav, centraljetPFlav, centraljetQGL, centraljetPUID;
  std::vector<double> centraljetGenpt, centraljetGeneta, centraljetGenphi, centraljetGenm;
  std::vector<double> centraljetBtagSF, centraljetBtagSFUp, centraljetBtagSFDown;
  // AK4CHS forward jet
  std::vector<double> forwardjetpt, forwardjeteta, forwardjetphi, forwardjetm, forwardjetbtag;
  std::vector<double> forwardjetCHfrac, forwardjetNHfrac, forwardjetEMfrac, forwardjetCEMfrac, forwardjetmetdphi;
  std::vector<double> forwardjetHFlav, forwardjetPFlav, forwardjetQGL, forwardjetPUID;
  std::vector<double> forwardjetGenpt, forwardjetGeneta, forwardjetGenphi, forwardjetGenm;

  double jetmetdphimin , incjetmetdphimin , jetmumetdphimin , incjetmumetdphimin, jetelmetdphimin , incjetelmetdphimin , jetphmetdphimin , incjetphmetdphimin , jetjetdphi;
  double jetmetdphimin4, incjetmetdphimin4, jetmumetdphimin4, incjetmumetdphimin4 , jetelmetdphimin4, incjetelmetdphimin4, jetphmetdphimin4, incjetphmetdphimin4, ht; 
  double alljetmetdphimin, alljetmetdphimin4, alljetmumetdphimin, alljetmumetdphimin4, alljetelmetdphimin, alljetelmetdphimin4, alljetphmetdphimin, alljetphmetdphimin4;

  // Puppijet ak4 puppi
  double leadingPuppijetpt, leadingPuppijeteta, leadingPuppijetphi, leadingPuppijetm; 

  // AK4Puppi central jet
  std::vector<double> centralPuppijetpt,     centralPuppijeteta,      centralPuppijetphi,    centralPuppijetm,      centralPuppijetbtag;
  std::vector<double> centralPuppijetCHfrac, centralPuppijetNHfrac,   centralPuppijetEMfrac, centralPuppijetCEMfrac,centralPuppijetmetdphi;
  std::vector<double> centralPuppijetHFlav,  centralPuppijetPFlav,    centralPuppijetQGL,    centralPuppijetPUID;
  std::vector<double> centralPuppijetGenpt,  centralPuppijetGeneta,   centralPuppijetGenphi, centralPuppijetGenm;
  std::vector<double> centralPuppijetBtagSF, centralPuppijetBtagSFUp, centralPuppijetBtagSFDown;
  // AK4Puppi forward Puppijet
  std::vector<double> forwardPuppijetpt,     forwardPuppijeteta,    forwardPuppijetphi,    forwardPuppijetm,       forwardPuppijetbtag;
  std::vector<double> forwardPuppijetCHfrac, forwardPuppijetNHfrac, forwardPuppijetEMfrac, forwardPuppijetCEMfrac, forwardPuppijetmetdphi;
  std::vector<double> forwardPuppijetHFlav,  forwardPuppijetPFlav,  forwardPuppijetQGL,    forwardPuppijetPUID;
  std::vector<double> forwardPuppijetGenpt,  forwardPuppijetGeneta, forwardPuppijetGenphi, forwardPuppijetGenm;
  //
  double Puppijetmetdphimin , incPuppijetmetdphimin , Puppijetmumetdphimin , incPuppijetmumetdphimin , Puppijetelmetdphimin , incPuppijetelmetdphimin , Puppijetphmetdphimin , incPuppijetphmetdphimin , PuppijetPuppijetdphi;
  double Puppijetmetdphimin4, incPuppijetmetdphimin4, Puppijetmumetdphimin4, incPuppijetmumetdphimin4, Puppijetelmetdphimin4, incPuppijetelmetdphimin4, Puppijetphmetdphimin4, incPuppijetphmetdphimin4, Puppiht; 
  
  // AK8 CHS jets
  std::vector<double> boostedJetpt, boostedJeteta, boostedJetphi, boostedJetm;
  std::vector<double> boostedJetGenpt, boostedJetGenm, boostedJetGeneta, boostedJetGenphi;
  std::vector<double> boostedJettau1, boostedJettau2, boostedJettau3, boostedJettau4, boostedJetecf1, boostedJetecf2, boostedJetecf3;  
  std::vector<double> boostedJetGentau1,boostedJetGentau2,boostedJetGentau3, boostedJetGentau4;
  std::vector<double> boostedJetHFlav, boostedJetPFlav, boostedJetQGL, boostedJetBtag, boostedJetDoubleBtag;
  std::vector<double> boostedJetBosonpt, boostedJetBosoneta, boostedJetBosonphi, boostedJetBosonm;
  std::vector<double> prunedJetpt, prunedJetm, prunedJetphi,prunedJeteta,  prunedJetGenpt, prunedJetGenm, prunedJetGeneta,prunedJetGenphi;
  std::vector<double> prunedJetptraw, prunedJetmraw;
  std::vector<double> prunedJetHFlav, prunedJetPFlav, prunedJetQGL, prunedJetBtag, prunedJetDoubleBtag;
  std::vector<double> softDropJetpt, softDropJetm, softDropJeteta, softDropJetphi, softDropJetGenpt, softDropJetGenm,softDropJetGeneta, softDropJetGenphi;
  std::vector<double> softDropJetHFlav, softDropJetPFlav, softDropJetQGL, softDropJetBtag, softDropJetDoubleBtag;
  std::vector<double> softDropJetptraw, softDropJetmraw;
  std::vector<double> prunedSubJetpt_1, prunedSubJetm_1,prunedSubJetphi_1, prunedSubJeteta_1, prunedSubJetHFlav_1, prunedSubJetQGL_1, prunedSubJetBtag_1;
  std::vector<double> prunedSubJetGenpt_1, prunedSubJetGenm_1,prunedSubJetGeneta_1, prunedSubJetGenphi_1, prunedSubJetPFlav_1;
  std::vector<double> prunedSubJetptraw_1, prunedSubJetmraw_1;
  std::vector<double> prunedSubJetpt_2, prunedSubJetm_2,prunedSubJetphi_2, prunedSubJeteta_2, prunedSubJetHFlav_2, prunedSubJetQGL_2, prunedSubJetBtag_2;
  std::vector<double> prunedSubJetGenpt_2, prunedSubJetGenm_2, prunedSubJetGeneta_2, prunedSubJetGenphi_2, prunedSubJetPFlav_2;
  std::vector<double> prunedSubJetptraw_2, prunedSubJetmraw_2;
  std::vector<double> softDropSubJetpt_1, softDropSubJetm_1,softDropSubJetphi_1, softDropSubJeteta_1;
  std::vector<double> softDropSubJetHFlav_1, softDropSubJetQGL_1, softDropSubJetBtag_1;
  std::vector<double> softDropSubJetGenpt_1, softDropSubJetGenm_1, softDropSubJetGeneta_1, softDropSubJetGenphi_1, softDropSubJetPFlav_1;
  std::vector<double> softDropSubJetptraw_1, softDropSubJetmraw_1;
  std::vector<double> softDropSubJetpt_2, softDropSubJetm_2,softDropSubJetphi_2, softDropSubJeteta_2, softDropSubJetHFlav_2; 
  std::vector<double> softDropSubJetQGL_2, softDropSubJetBtag_2;
  std::vector<double> softDropSubJetGenpt_2, softDropSubJetGenm_2, softDropSubJetGeneta_2, softDropSubJetGenphi_2, softDropSubJetPFlav_2;
  std::vector<double> softDropSubJetptraw_2, softDropSubJetmraw_2;

  // AK8 Puppi jets
  std::vector<double> boostedPuppiJetpt, boostedPuppiJeteta, boostedPuppiJetphi, boostedPuppiJetm;
  std::vector<double> boostedPuppiJetGenpt, boostedPuppiJetGenm, boostedPuppiJetGeneta, boostedPuppiJetGenphi;
  std::vector<double> boostedPuppiJettau1, boostedPuppiJettau2, boostedPuppiJettau3, boostedPuppiJettau4, boostedPuppiJetecf1, boostedPuppiJetecf2, boostedPuppiJetecf3;
  std::vector<double> boostedPuppiJetGentau1, boostedPuppiJetGentau2, boostedPuppiJetGentau3, boostedPuppiJetGentau4;
  std::vector<double> boostedPuppiJetHFlav, boostedPuppiJetPFlav, boostedPuppiJetQGL, boostedPuppiJetBtag, boostedPuppiJetDoubleBtag;
  std::vector<double> boostedPuppiJetBosonpt, boostedPuppiJetBosoneta, boostedPuppiJetBosonphi, boostedPuppiJetBosonm;
  std::vector<double> prunedPuppiJetpt, prunedPuppiJetm, prunedPuppiJetphi,prunedPuppiJeteta,  prunedPuppiJetGenpt, prunedPuppiJetGenm, prunedPuppiJetGeneta,prunedPuppiJetGenphi;
  std::vector<double> prunedPuppiJetptraw, prunedPuppiJetmraw;
  std::vector<double> prunedPuppiJetHFlav, prunedPuppiJetPFlav, prunedPuppiJetQGL, prunedPuppiJetBtag, prunedPuppiJetDoubleBtag;
  std::vector<double> softDropPuppiJetpt, softDropPuppiJetm, softDropPuppiJeteta, softDropPuppiJetphi, softDropPuppiJetGenpt, softDropPuppiJetGenm,softDropPuppiJetGeneta, softDropPuppiJetGenphi;
  std::vector<double> softDropPuppiJetHFlav, softDropPuppiJetPFlav, softDropPuppiJetQGL, softDropPuppiJetBtag, softDropPuppiJetDoubleBtag;
  std::vector<double> softDropPuppiJetptraw, softDropPuppiJetmraw;
  std::vector<double> prunedPuppiSubJetpt_1, prunedPuppiSubJetm_1,prunedPuppiSubJetphi_1, prunedPuppiSubJeteta_1, prunedPuppiSubJetHFlav_1, prunedPuppiSubJetQGL_1, prunedPuppiSubJetBtag_1;
  std::vector<double> prunedPuppiSubJetGenpt_1, prunedPuppiSubJetGenm_1,prunedPuppiSubJetGeneta_1, prunedPuppiSubJetGenphi_1, prunedPuppiSubJetPFlav_1;
  std::vector<double> prunedPuppiSubJetptraw_1, prunedPuppiSubJetmraw_1;
  std::vector<double> prunedPuppiSubJetpt_2, prunedPuppiSubJetm_2,prunedPuppiSubJetphi_2, prunedPuppiSubJeteta_2, prunedPuppiSubJetHFlav_2, prunedPuppiSubJetQGL_2, prunedPuppiSubJetBtag_2;
  std::vector<double> prunedPuppiSubJetGenpt_2, prunedPuppiSubJetGenm_2, prunedPuppiSubJetGeneta_2, prunedPuppiSubJetGenphi_2, prunedPuppiSubJetPFlav_2;
  std::vector<double> prunedPuppiSubJetptraw_2, prunedPuppiSubJetmraw_2;
  std::vector<double> softDropPuppiSubJetpt_1, softDropPuppiSubJetm_1,softDropPuppiSubJetphi_1, softDropPuppiSubJeteta_1;
  std::vector<double> softDropPuppiSubJetHFlav_1, softDropPuppiSubJetQGL_1, softDropPuppiSubJetBtag_1;
  std::vector<double> softDropPuppiSubJetGenpt_1, softDropPuppiSubJetGenm_1, softDropPuppiSubJetGeneta_1, softDropPuppiSubJetGenphi_1, softDropPuppiSubJetPFlav_1;
  std::vector<double> softDropPuppiSubJetptraw_1, softDropPuppiSubJetmraw_1;
  std::vector<double> softDropPuppiSubJetpt_2, softDropPuppiSubJetm_2,softDropPuppiSubJetphi_2, softDropPuppiSubJeteta_2, softDropPuppiSubJetHFlav_2; 
  std::vector<double> softDropPuppiSubJetQGL_2, softDropPuppiSubJetBtag_2;
  std::vector<double> softDropPuppiSubJetGenpt_2, softDropPuppiSubJetGenm_2, softDropPuppiSubJetGeneta_2, softDropPuppiSubJetGenphi_2, softDropPuppiSubJetPFlav_2;
  std::vector<double> softDropPuppiSubJetptraw_2, softDropPuppiSubJetmraw_2;

  // gen info leptoni W/Z boson (1 per event)
  double wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
  // photon info
  double parpt, pareta, parphi, ancpt, anceta, ancphi;
  // hadronic V and related quarks (1 per event)
  double wzmass_h, wzmt_h, wzpt_h, wzeta_h, wzphi_h, q1pt, q1eta, q1phi, q2pt, q2eta, q2phi;
  // one top
  double topmass, toppt, topeta, topphi;
  // second top
  double atopmass, atoppt, atopeta, atopphi;
  // DM mediator and DM particles
  double dmmass, dmpt, dmeta, dmphi, dmX1pt, dmX1eta, dmX1phi, dmX1mass, dmX2pt, dmX2eta, dmX2phi, dmX2mass;
  int    dmid, dmX1id, dmX2id;
  // for fastSIM
  double samplemedM, sampledmM;

  // weights
  double wgt, kfact, puwgt, pswgt;

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
  // isMC or Data --> default Data
  isMC(iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  // use lhe weights or not
  uselheweights(iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
  lheEventTag(iConfig.getParameter<edm::InputTag>("lheinfo")),
  lheRunTag(iConfig.getParameter<edm::InputTag>("lheRuninfo")),
  // is signal sample or not
  isWorZorSignalMCSample(iConfig.existsAs<bool>("isWorZorSignalMCSample") ? iConfig.getParameter<bool>("isWorZorSignalMCSample") : false),
  // xsec
  xsec(iConfig.existsAs<double>("xsec") ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.),
  ///////////// TRIGGER and filter info INFO
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
  filterResultsTag(iConfig.getParameter<edm::InputTag>("filterResults")),
  prescalesTag(iConfig.getParameter<edm::InputTag>("prescales")),
  hbhelooseTag(iConfig.getParameter<edm::InputTag>("hbheloose")),
  hbhetightTag(iConfig.getParameter<edm::InputTag>("hbhetight")),
  hbheisoTag(iConfig.getParameter<edm::InputTag>("hbheiso")),
  // vertexes
  verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
  //muons
  muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
  tightmuonsTag(iConfig.getParameter<edm::InputTag>("tightmuons")),
  highptmuonsTag(iConfig.getParameter<edm::InputTag>("highptmuons")),
  // electrons
  electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
  tightelectronsTag(iConfig.getParameter<edm::InputTag>("tightelectrons")),
  heepelectronsTag(iConfig.getParameter<edm::InputTag>("heepelectrons")),
  electronLooseIdTag(iConfig.getParameter<edm::InputTag>("electronLooseId")),
  // photons
  photonsTag(iConfig.getParameter<edm::InputTag>("photons")),
  tightphotonsTag(iConfig.getParameter<edm::InputTag>("tightphotons")),
  photonLooseIdTag(iConfig.getParameter<edm::InputTag>("photonLooseId")),
  photonMediumIdTag(iConfig.getParameter<edm::InputTag>("photonMediumId")),
  photonTightIdTag(iConfig.getParameter<edm::InputTag>("photonTightId")),
  photonHighPtIdTag(iConfig.getParameter<edm::InputTag>("photonHighPtId")),
  // taus
  tausTag(iConfig.getParameter<edm::InputTag>("taus")),
  // jets AK4
  jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
  addPuppiJets(iConfig.existsAs<bool>("addPuppiJets") ? iConfig.getParameter<bool>("addPuppiJets") : false),
  // met
  t1metTag(iConfig.getParameter<edm::InputTag>("t1met")),
  t1mumetTag(iConfig.getParameter<edm::InputTag>("t1mumet")),
  t1elmetTag(iConfig.getParameter<edm::InputTag>("t1elmet")),
  t1phmetTag(iConfig.getParameter<edm::InputTag>("t1phmet")),
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
  addSubstructurePuppi(iConfig.existsAs<bool>("addSubstructurePuppi") ? iConfig.getParameter<bool>("addSubstructurePuppi") : false),
  addBTagScaleFactor(iConfig.existsAs<bool>("addBTagScaleFactor") ? iConfig.getParameter<bool>("addBTagScaleFactor") : false){

  usesResource();
  usesResource("TFileService");

  // trigger tokens
  triggerResultsToken   = consumes<edm::TriggerResults> (triggerResultsTag);
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(prescalesTag);
  filterResultsToken    = consumes<edm::TriggerResults> (filterResultsTag);

  //HBHE Noise
  hbhelooseToken = consumes<bool> (hbhelooseTag);
  hbhetightToken = consumes<bool> (hbhetightTag);
  hbheisoToken   = consumes<bool> (hbheisoTag);
  verticesToken  = consumes<std::vector<reco::Vertex> > (verticesTag);

  //muons
  muonsToken       = consumes<pat::MuonRefVector> (muonsTag);
  tightmuonsToken  = consumes<pat::MuonRefVector> (tightmuonsTag);
  highptmuonsToken = consumes<pat::MuonRefVector> (highptmuonsTag);

  // electrons
  electronsToken       = consumes<pat::ElectronRefVector> (electronsTag);
  tightelectronsToken  = consumes<pat::ElectronRefVector>(tightelectronsTag);
  heepelectronsToken   = consumes<pat::ElectronRefVector> (heepelectronsTag);
  electronLooseIdToken = consumes<edm::ValueMap<bool> > (electronLooseIdTag);
  // photons
  photonsToken        = consumes<pat::PhotonRefVector> (photonsTag);
  tightphotonsToken   = consumes<pat::PhotonRefVector> (tightphotonsTag);
  photonLooseIdToken  = consumes<edm::ValueMap<bool> > (photonLooseIdTag);
  photonMediumIdToken = consumes<edm::ValueMap<bool> > (photonMediumIdTag);
  photonTightIdToken  = consumes<edm::ValueMap<bool> > (photonTightIdTag);
  photonHighPtIdToken = consumes<edm::ValueMap<bool> > (photonHighPtIdTag);
  // taus
  tausToken = consumes<std::vector<pat::Tau> > (tausTag);
  // jets AK4
  jetsToken = consumes<std::vector<pat::Jet> > (jetsTag);

  t1metToken    = consumes<edm::View<pat::MET> > (t1metTag);
  t1mumetToken  = consumes<edm::View<pat::MET> > (t1mumetTag);
  t1elemetToken = consumes<edm::View<pat::MET> > (t1elmetTag);
  t1phmetToken  = consumes<edm::View<pat::MET> > (t1phmetTag);
   
  // only for simulated samples
  if( isMC ){
    pileupInfoToken = consumes<std::vector<PileupSummaryInfo> > (iConfig.getParameter<edm::InputTag>("pileup"));
    genevtInfoToken = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genevt"));
    lheInfoToken    = consumes<LHEEventProduct> (lheEventTag);
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

  if(addBTagScaleFactor){
    bTagScaleFactorFile = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFile");
    if ( bTagScaleFactorFile.location()!=edm::FileInPath::Local)
      throw cms::Exception("MonoJetTreeMaker") << " Failed to find File = " << bTagScaleFactorFile << " !!\n";

    calib = std::auto_ptr<BTagCalibration>(new BTagCalibration("CSVv2",bTagScaleFactorFile.fullPath()));

    bMedium.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"comb","central")); // for light flavor
    bMedium.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"mujets","central")); // for b and c-jets
    bMediumUp.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"comb","up"));
    bMediumUp.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"mujets","up"));
    bMediumDown.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"comb","down"));
    bMediumDown.push_back(BTagCalibrationReader(calib.get(),BTagEntry::OP_MEDIUM,"mujets","down"));
  }

  readDMFromGenParticle = false;

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
    Handle<pat::PackedTriggerPrescales> triggerPrescalesH;
    iEvent.getByToken(triggerPrescalesToken, triggerPrescalesH);
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
      if (isWorZorSignalMCSample)
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
    hltphoton120    = 0;
    hltdoublemu     = 0;
    hltsinglemu     = 0;
    hltdoubleel     = 0;
    hltsingleel     = 0;
    hltelnoiso      = 0;

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
        if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton120    = 1; // Photon trigger
        if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
        if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
        if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 29 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
        if (i == 30 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 32 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 33 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
        if (i == 34 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger
        if (i == 35 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger
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
    if (hltphoton120    == 1) triggered = true;
    if (hltdoublemu     == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltdoubleel     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (hltelnoiso      == 1) triggered = true;
    if (applyHLTFilter && !triggered) return;

    pswgt = 1.0;
    const edm::TriggerNames &trignames = iEvent.triggerNames(*triggerResultsH);
    for (size_t i = 0; i < triggerResultsH->size(); i++) {
        if (trignames.triggerName(i).find("HLT_Photon120_v") != string::npos) pswgt = triggerPrescalesH->getPrescaleForIndex(i);
    }

    // MET filter info
    flagcsctight  = 0;
    flaghbhenoise = 0;
    flageebadsc   = 0;

    // HBHE Noise 
    if(hbhelooseH.isValid())
      flaghbheloose = (*hbhelooseH ? 1 : 0);
    if(hbhetightH.isValid())
      flaghbhetight = (*hbhetightH ? 1 : 0);
    if(hbheisoH.isValid())
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

    // in caase the cross section is not set from outside --> fix to 1 as dummy value
    if (uselheweights && genevtInfoH.isValid()){
      wgt = genevtInfoH->weight();
    }
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

    vector<pat::JetRef> alljets;
    vector<pat::JetRef> incjets;
    vector<pat::JetRef> forwardjets;
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
      
	pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());	
	if(jetref.isAvailable() and jetref.isNonnull()) alljets.push_back(jetref);

	// apply jet id
	bool passjetid = applyJetID(*jets_iter,"loose");            
	if (!passjetid) 
	  continue;
	
	// apply pileup jet id
	bool passpuid = applyPileupJetID(*jets_iter,"medium",false);
	if (!passpuid) 
	  continue;

	if (jets_iter->pt() > leadingjetpt) {
	  leadingjetpt  = jets_iter->pt() ;
	  leadingjeteta = jets_iter->eta();
	  leadingjetphi = jets_iter->phi();
	  leadingjetm   = jets_iter->mass();
	}
	            	
	if(jetref.isAvailable() and jetref.isNonnull())
	  incjets.push_back(jetref);
      }
    
      if(incjets.size() > 0)
	sort(incjets.begin(), incjets.end(), jetSorter);
      
      // only central jets
      for (size_t i = 0; i < incjets.size(); i++) {
	if (fabs(incjets[i]->eta()) <= 2.5) 
	  jets.push_back(incjets[i]);
	else
	  forwardjets.push_back(incjets[i]);
      }        
      
      // sort them in pt
      if(jets.size() > 0)
	sort(jets.begin(), jets.end(), jetSorter);
      if(forwardjets.size() > 0)
	sort(forwardjets.begin(), forwardjets.end(), jetSorter);

      // count central jets    
      njets       = 0;
      njetsinc    = 0;
      nbjets      = 0;
      nbjetslowpt = 0;
      
      float minJetPtCentral  = 20.;
      float minJetPtForward  = 30.;
      float minJetPtCount = 30.;

      centraljetpt        .clear(); centraljeteta       .clear(); centraljetphi       .clear(); centraljetbtag      .clear(); centraljetCHfrac    .clear();
      centraljetNHfrac    .clear(); centraljetEMfrac    .clear(); centraljetCEMfrac   .clear(); centraljetmetdphi   .clear();
      centraljetNHfrac    .clear(); centraljetEMfrac    .clear(); centraljetCEMfrac   .clear(); centraljetmetdphi   .clear();
      centraljetHFlav     .clear(); centraljetPFlav     .clear(); centraljetQGL       .clear(); centraljetPUID      .clear();
      centraljetGenpt     .clear(); centraljetGeneta    .clear(); centraljetGenphi    .clear(); centraljetGenm      .clear(); 
      centraljetm         .clear(); 
      centraljetBtagSF .clear(); centraljetBtagSFUp .clear(); centraljetBtagSFDown .clear();

      for(size_t i = 0; i < incjets.size(); i++){
	if(incjets[i]->pt() > 30) njetsinc++;
      }
  
      for (size_t i = 0; i < jets.size(); i++) {

	if (jets[i]->pt() > minJetPtCount) njets++;
	if (jets[i]->pt() > minJetPtCount && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjets++;
	if (jets[i]->pt() > 15 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjetslowpt++;

	// fill collections
	if (jets[i]->pt() > minJetPtCentral){

	  centraljetpt.push_back(jets[i]->pt());
	  centraljeteta.push_back(jets[i]->eta());
	  centraljetphi.push_back(jets[i]->phi());
	  centraljetm.push_back(jets[i]->mass());
	  centraljetbtag.push_back(jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	  centraljetCHfrac  .push_back(jets[i]->chargedHadronEnergyFraction());
	  centraljetNHfrac  .push_back(jets[i]->neutralHadronEnergyFraction());
	  centraljetEMfrac  .push_back(jets[i]->neutralEmEnergyFraction());
	  centraljetCEMfrac .push_back(jets[i]->chargedEmEnergyFraction());
	  
	  if(jets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	    centraljetQGL   .push_back(jets[i]->userFloat("QGTagger:qgLikelihood")); 
	  // pileup jet id
	  if(jets[i]->hasUserFloat("puid:fullDiscriminant"))
	    centraljetPUID  .push_back(jets[i]->userFloat("puid:fullDiscriminant"));	
	  else
	    centraljetPUID  .push_back(jets[i]->userFloat("pileupJetId:fullDiscriminant"));
	  // MC based info
	  if(isMC){
	    centraljetHFlav.push_back(jets[i]->hadronFlavour()); 
	    centraljetPFlav.push_back(jets[i]->partonFlavour()); 
	    if(jets[i]->genJet()){
	      centraljetGenpt.push_back(jets[i]->genJet()->pt()); 
	      centraljetGeneta.push_back(jets[i]->genJet()->eta()); 
	      centraljetGenphi.push_back(jets[i]->genJet()->phi()); 
	      centraljetGenm.push_back(jets[i]->genJet()->mass()); 
	    }	  
	    
	    // b-tag SF for jets
	    if(addBTagScaleFactor){
	      
	      float MaxBJetPt = 670.;
	      float minBJetPt = 30.;
	      float maxEta    = 2.4;
	      	      
	      if(jets[i]->hadronFlavour() == 0){
		MaxBJetPt = 1000.;
		minBJetPt = 20.;
	      }
	      
	      float jetPt = jets[i]->pt();
	      bool  doubleUncertainty = false;
	      if(jetPt > MaxBJetPt){
		jetPt = MaxBJetPt;
		doubleUncertainty = true;
	      }
	      if(jetPt < minBJetPt){
		jetPt = minBJetPt;
		doubleUncertainty = true;
	      }
	      
	      float jetEta = jets[i]->eta();
	      if(fabs(jetEta) > maxEta){
		jetEta = maxEta;
		doubleUncertainty = true;
	      }
	      
	      if(jets[i]->hadronFlavour() == 5){
		centraljetBtagSF.push_back(bMedium[1].eval(BTagEntry::FLAV_B, jetEta, jetPt));	      
		centraljetBtagSFUp.push_back(bMediumUp[1].eval(BTagEntry::FLAV_B, jetEta, jetPt));
		centraljetBtagSFDown.push_back(bMediumDown[1].eval(BTagEntry::FLAV_B, jetEta, jetPt));		
	      }
	      if(jets[i]->hadronFlavour() == 4){
		centraljetBtagSF.push_back(bMedium[1].eval(BTagEntry::FLAV_C, jetEta, jetPt));
		centraljetBtagSFUp.push_back(bMediumUp[1].eval(BTagEntry::FLAV_C, jetEta, jetPt));
		centraljetBtagSFDown.push_back(bMediumDown[1].eval(BTagEntry::FLAV_C, jetEta, jetPt));
	      }
	      else{
		centraljetBtagSF.push_back(bMedium[0].eval(BTagEntry::FLAV_UDSG, jetEta, jetPt));
		centraljetBtagSFUp.push_back(bMediumUp[0].eval(BTagEntry::FLAV_UDSG, jetEta, jetPt));
		centraljetBtagSFDown.push_back(bMediumDown[0].eval(BTagEntry::FLAV_UDSG, jetEta, jetPt));	      
	      }	      
	      if(doubleUncertainty){
		centraljetBtagSFUp.back() = 2*(centraljetBtagSFUp.back()-centraljetBtagSF.back())+centraljetBtagSF.back();
		centraljetBtagSFDown.back() = 2*(centraljetBtagSFDown.back()-centraljetBtagSF.back())+centraljetBtagSF.back();
	      }
	    }
	  }
	  // fill jet met dphi
	  centraljetmetdphi.push_back(deltaPhi(jets[i]->phi(), t1pfmetphi));
	}
      }
          
      forwardjetpt        .clear(); forwardjeteta       .clear(); forwardjetphi       .clear(); forwardjetbtag      .clear(); forwardjetCHfrac    .clear();
      forwardjetNHfrac    .clear(); forwardjetEMfrac    .clear(); forwardjetCEMfrac   .clear(); forwardjetmetdphi   .clear();
      forwardjetHFlav     .clear(); forwardjetPFlav     .clear(); forwardjetQGL       .clear(); forwardjetPUID      .clear();
      forwardjetGenpt     .clear(); forwardjetGeneta    .clear(); forwardjetGenphi    .clear(); forwardjetGenm      .clear(); 
      forwardjetm         .clear(); 

      // forward jets
      for (size_t i = 0; i < forwardjets.size(); i++) {

	// fill collections
	if (forwardjets[i]->pt() > minJetPtForward){

	  forwardjetpt.push_back(forwardjets[i]->pt());
	  forwardjeteta.push_back(forwardjets[i]->eta());
	  forwardjetphi.push_back(forwardjets[i]->phi());
	  forwardjetm.push_back(forwardjets[i]->mass());
	  forwardjetbtag.push_back(forwardjets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	  forwardjetCHfrac  .push_back(forwardjets[i]->chargedHadronEnergyFraction());
	  forwardjetNHfrac  .push_back(forwardjets[i]->neutralHadronEnergyFraction());
	  forwardjetEMfrac  .push_back(forwardjets[i]->neutralEmEnergyFraction());
	  forwardjetCEMfrac .push_back(forwardjets[i]->chargedEmEnergyFraction());
	
	  if(forwardjets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	    forwardjetQGL   .push_back(forwardjets[i]->userFloat("QGTagger:qgLikelihood")); 
	  // pileup jet id
	  if(forwardjets[i]->hasUserFloat("puid:fullDiscriminant"))
	    forwardjetPUID  .push_back(forwardjets[i]->userFloat("puid:fullDiscriminant"));	
	  else
	  forwardjetPUID  .push_back(forwardjets[i]->userFloat("pileupJetId:fullDiscriminant"));
	  // MC based info
	  if(isMC){
	    forwardjetHFlav.push_back(forwardjets[i]->hadronFlavour()); 
	    forwardjetPFlav.push_back(forwardjets[i]->partonFlavour()); 
	    if(forwardjets[i]->genJet()){
	      forwardjetGenpt.push_back(forwardjets[i]->genJet()->pt()); 
	      forwardjetGeneta.push_back(forwardjets[i]->genJet()->eta()); 
	      forwardjetGenphi.push_back(forwardjets[i]->genJet()->phi()); 
	      forwardjetGenm.push_back(forwardjets[i]->genJet()->mass()); 
	    }	  	    
	  }
	}
      }
      
      
      jetjetdphi         = 0.0;   
      jetmetdphimin      = 0.0;   incjetmetdphimin    = 0.0; alljetmetdphimin    = 0.0;
      jetmumetdphimin    = 0.0;   incjetmumetdphimin  = 0.0; alljetmumetdphimin  = 0.0;
      jetelmetdphimin    = 0.0;   incjetelmetdphimin  = 0.0; alljetelmetdphimin  = 0.0;
      jetphmetdphimin    = 0.0;   incjetphmetdphimin  = 0.0; alljetphmetdphimin  = 0.0;
      jetmetdphimin4     = 0.0;   incjetmetdphimin4   = 0.0; alljetmetdphimin4   = 0.0;
      jetmumetdphimin4   = 0.0;   incjetmumetdphimin4 = 0.0; alljetmumetdphimin4 = 0.0;
      jetelmetdphimin4   = 0.0;   incjetelmetdphimin4 = 0.0; alljetelmetdphimin4 = 0.0;
      jetphmetdphimin4   = 0.0;   incjetphmetdphimin4 = 0.0; alljetphmetdphimin4 = 0.0;

      // delta phi between jets
      if (centraljetphi.size() > 1)
	jetjetdphi = deltaPhi(centraljetphi[0], centraljetphi[1]);

      std::vector<double> alljetmetdphiminvector;
      std::vector<double> alljetmetdphimin4vector;
      std::vector<double> alljetmumetdphiminvector;
      std::vector<double> alljetmumetdphimin4vector;
      std::vector<double> alljetelmetdphiminvector;
      std::vector<double> alljetelmetdphimin4vector;
      std::vector<double> alljetphmetdphiminvector;
      std::vector<double> alljetphmetdphimin4vector;
      for (size_t i = 0; i < alljets.size(); i++) {
    if (alljets[i]->pt() > minJetPtCount) {
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
      std::vector<double> jetmetdphiminvector;
      std::vector<double> jetmetdphimin4vector;
      // only central jets
      for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > minJetPtCount) {
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
	if (incjets[i]->pt() > minJetPtCount) {
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
	if (jets[i]->pt() >minJetPtCount) {
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
        if (incjets[i]->pt() > minJetPtCount) {
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
	if (jets[i]->pt() > minJetPtCount) {
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
        if (incjets[i]->pt() > minJetPtCount) {
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
	if (jets[i]->pt() > minJetPtCount) {
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
        if (incjets[i]->pt() > minJetPtCount) {
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
	if (jets[i]->pt() > minJetPtCount) {
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
      vector<pat::JetRef> forwardPuppijets;

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
	  
	  // apply jet id
	  bool passjetid = applyJetID(*jets_iter,"loose");            
	  if (!passjetid) 
	    continue;
	  
	  //apply pileup jet id
	  bool passpuid = applyPileupJetID(*jets_iter,"medium",true);
	  if (!passpuid) continue;

	  if (jets_iter->pt() > leadingPuppijetpt) {
	    leadingPuppijetpt  = jets_iter->pt() ;
	    leadingPuppijeteta = jets_iter->eta();
	    leadingPuppijetphi = jets_iter->phi();
	    leadingPuppijetm   = jets_iter->mass();
	  }
	  	  
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
	  else 
	    forwardPuppijets.push_back(incPuppijets[i]);
	}        
	
	// sort them in pt
	if(Puppijets.size() > 0)
	  sort(Puppijets.begin(), Puppijets.end(), jetSorter);
	if(forwardPuppijets.size() > 0)
	  sort(forwardPuppijets.begin(), forwardPuppijets.end(), jetSorter);
	
	// count central jets    
	npuppijets       = 0;
	npuppijetsinc    = 0;
	npuppibjets      = 0;
	npuppibjetslowpt = 0;
	
	float minPuppiJetPtCentral  = 20.;
	float minPuppiJetPtForward  = 30.;
	float minPuppiJetPtCount = 30.;
	
	centralPuppijetpt        .clear(); centralPuppijeteta       .clear(); centralPuppijetphi       .clear(); centralPuppijetbtag      .clear(); 
	centralPuppijetCHfrac    .clear();
	centralPuppijetNHfrac    .clear(); centralPuppijetEMfrac    .clear(); centralPuppijetCEMfrac   .clear(); centralPuppijetmetdphi   .clear();
	centralPuppijetHFlav     .clear(); centralPuppijetPFlav     .clear(); centralPuppijetQGL       .clear(); centralPuppijetPUID      .clear();
	centralPuppijetGenpt     .clear(); centralPuppijetGeneta    .clear(); centralPuppijetGenphi    .clear(); centralPuppijetGenm      .clear();
	centralPuppijetm         .clear(); 
	centralPuppijetBtagSF    .clear(); centralPuppijetBtagSFUp .clear(); centralPuppijetBtagSFDown .clear();
	
	for(size_t i = 0; i < incPuppijets.size(); i++){
	  if(incPuppijets[i]->pt() > 30) npuppijetsinc++;
	}
	
	for (size_t i = 0; i < Puppijets.size(); i++) {
	  
	  if (Puppijets[i]->pt() > minPuppiJetPtCount) npuppijets++;
	  if (Puppijets[i]->pt() > minPuppiJetPtCount && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BPuppijetTags") > 0.89) npuppibjets++;
	  if (Puppijets[i]->pt() > 15 && Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BPuppijetTags") > 0.89) npuppibjetslowpt++;
	  
	  // fill collections
	  if (Puppijets[i]->pt() > minPuppiJetPtCentral){
	    
	    centralPuppijetpt.push_back(Puppijets[i]->pt());
	    centralPuppijeteta.push_back(Puppijets[i]->eta());
	    centralPuppijetphi.push_back(Puppijets[i]->phi());
	    centralPuppijetm.push_back(Puppijets[i]->mass());
	    centralPuppijetbtag.push_back(Puppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BPuppijetTags"));
	    centralPuppijetCHfrac  .push_back(Puppijets[i]->chargedHadronEnergyFraction());
	    centralPuppijetNHfrac  .push_back(Puppijets[i]->neutralHadronEnergyFraction());
	    centralPuppijetEMfrac  .push_back(Puppijets[i]->neutralEmEnergyFraction());
	    centralPuppijetCEMfrac .push_back(Puppijets[i]->chargedEmEnergyFraction());
	    
	    if(Puppijets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	      centralPuppijetQGL   .push_back(Puppijets[i]->userFloat("QGTagger:qgLikelihood")); 
	    // pileup Puppijet id
	    if(Puppijets[i]->hasUserFloat("puid:fullDiscriminant"))
	      centralPuppijetPUID  .push_back(Puppijets[i]->userFloat("puid:fullDiscriminant"));	
	    else
	      centralPuppijetPUID  .push_back(Puppijets[i]->userFloat("pileupPuppijetId:fullDiscriminant"));
	    // MC based info
	    if(isMC){
	      centralPuppijetHFlav.push_back(Puppijets[i]->hadronFlavour()); 
	      centralPuppijetPFlav.push_back(Puppijets[i]->partonFlavour()); 
	      if(Puppijets[i]->genJet()){
		centralPuppijetGenpt.push_back(Puppijets[i]->genJet()->pt()); 
		centralPuppijetGeneta.push_back(Puppijets[i]->genJet()->eta()); 
		centralPuppijetGenphi.push_back(Puppijets[i]->genJet()->phi()); 
		centralPuppijetGenm.push_back(Puppijets[i]->genJet()->mass()); 
	      }	  
	    
	      // b-tag SF for Puppijets
	      if(addBTagScaleFactor){
		
		float MaxBPuppijetPt = 670.;
		float minBPuppijetPt = 30.;
		float maxEta    = 2.4;
	      	      
		if(Puppijets[i]->hadronFlavour() == 0){
		  MaxBPuppijetPt = 1000.;
		  minBPuppijetPt = 20.;
		}
		
		float PuppijetPt = Puppijets[i]->pt();
		bool  doubleUncertainty = false;
		if(PuppijetPt > MaxBPuppijetPt){
		  PuppijetPt = MaxBPuppijetPt;
		  doubleUncertainty = true;
		}
		if(PuppijetPt < minBPuppijetPt){
		  PuppijetPt = minBPuppijetPt;
		  doubleUncertainty = true;
		}
		
		float PuppijetEta = Puppijets[0]->eta();
		if(fabs(PuppijetEta) > maxEta){
		  PuppijetEta = maxEta;
		  doubleUncertainty = true;
		}
		
		if(Puppijets[i]->hadronFlavour() == 5){
		  centralPuppijetBtagSF.push_back(bMedium[1].eval(BTagEntry::FLAV_B, PuppijetEta, PuppijetPt));	      
		  centralPuppijetBtagSFUp.push_back(bMediumUp[1].eval(BTagEntry::FLAV_B, PuppijetEta, PuppijetPt));
		  centralPuppijetBtagSFDown.push_back(bMediumDown[1].eval(BTagEntry::FLAV_B, PuppijetEta, PuppijetPt));		
		}
		if(Puppijets[i]->hadronFlavour() == 4){
		  centralPuppijetBtagSF.push_back(bMedium[1].eval(BTagEntry::FLAV_C, PuppijetEta, PuppijetPt));
		  centralPuppijetBtagSFUp.push_back(bMediumUp[1].eval(BTagEntry::FLAV_C, PuppijetEta, PuppijetPt));
		  centralPuppijetBtagSFDown.push_back(bMediumDown[1].eval(BTagEntry::FLAV_C, PuppijetEta, PuppijetPt));
		}
		else{
		  centralPuppijetBtagSF.push_back(bMedium[0].eval(BTagEntry::FLAV_UDSG, PuppijetEta, PuppijetPt));
		  centralPuppijetBtagSFUp.push_back(bMediumUp[0].eval(BTagEntry::FLAV_UDSG, PuppijetEta, PuppijetPt));
		  centralPuppijetBtagSFDown.push_back(bMediumDown[0].eval(BTagEntry::FLAV_UDSG, PuppijetEta, PuppijetPt));	      
		}	      
		if(doubleUncertainty){
		  centralPuppijetBtagSFUp.back() = 2*(centralPuppijetBtagSFUp.back()-centralPuppijetBtagSF.back())+centralPuppijetBtagSF.back();
		  centralPuppijetBtagSFDown.back() = 2*(centralPuppijetBtagSFDown.back()-centralPuppijetBtagSF.back())+centralPuppijetBtagSF.back();
		}
	      }
	    }
	    centralPuppijetmetdphi.push_back(deltaPhi(Puppijets[i]->phi(),t1pfmetphi));
	  }	  
	}
	
	forwardPuppijetpt        .clear(); forwardPuppijeteta       .clear(); forwardPuppijetphi       .clear(); forwardPuppijetbtag      .clear(); forwardPuppijetCHfrac    .clear();
	forwardPuppijetNHfrac    .clear(); forwardPuppijetEMfrac    .clear(); forwardPuppijetCEMfrac   .clear(); forwardPuppijetmetdphi   .clear();
	forwardPuppijetHFlav     .clear(); forwardPuppijetPFlav     .clear(); forwardPuppijetQGL       .clear(); forwardPuppijetPUID      .clear();
	forwardPuppijetGenpt     .clear(); forwardPuppijetGeneta    .clear(); forwardPuppijetGenphi    .clear(); forwardPuppijetGenm      .clear(); 
	forwardPuppijetm         .clear(); 

	// forward Puppijets
	for (size_t i = 0; i < forwardPuppijets.size(); i++) {
	  
	  // fill collections
	  if (forwardPuppijets[i]->pt() > minPuppiJetPtForward){
	    
	    forwardPuppijetpt.push_back(forwardPuppijets[i]->pt());
	    forwardPuppijeteta.push_back(forwardPuppijets[i]->eta());
	    forwardPuppijetphi.push_back(forwardPuppijets[i]->phi());
	    forwardPuppijetm.push_back(forwardPuppijets[i]->mass());
	    forwardPuppijetbtag.push_back(forwardPuppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BPuppijetTags"));
	    forwardPuppijetCHfrac  .push_back(forwardPuppijets[i]->chargedHadronEnergyFraction());
	    forwardPuppijetNHfrac  .push_back(forwardPuppijets[i]->neutralHadronEnergyFraction());
	    forwardPuppijetEMfrac  .push_back(forwardPuppijets[i]->neutralEmEnergyFraction());
	    forwardPuppijetCEMfrac .push_back(forwardPuppijets[i]->chargedEmEnergyFraction());
	    
	    if(forwardPuppijets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	      forwardPuppijetQGL   .push_back(forwardPuppijets[i]->userFloat("QGTagger:qgLikelihood")); 
	    // pileup Puppijet id
	    if(forwardPuppijets[i]->hasUserFloat("puid:fullDiscriminant"))
	      forwardPuppijetPUID  .push_back(forwardPuppijets[i]->userFloat("puid:fullDiscriminant"));	
	    else
	      forwardPuppijetPUID  .push_back(forwardPuppijets[i]->userFloat("pileupPuppijetId:fullDiscriminant"));
	    // MC based info
	    if(isMC){
	      forwardPuppijetHFlav.push_back(forwardPuppijets[i]->hadronFlavour()); 
	      forwardPuppijetPFlav.push_back(forwardPuppijets[i]->partonFlavour()); 
	      if(forwardPuppijets[i]->genJet()){
		forwardPuppijetGenpt.push_back(forwardPuppijets[i]->genJet()->pt()); 
		forwardPuppijetGeneta.push_back(forwardPuppijets[i]->genJet()->eta()); 
		forwardPuppijetGenphi.push_back(forwardPuppijets[i]->genJet()->phi()); 
		forwardPuppijetGenm.push_back(forwardPuppijets[i]->genJet()->mass()); 
	      }	  	      
	    }
	  }	  
	}      
	
	PuppijetPuppijetdphi    = 0.0;
	Puppijetmetdphimin      = 0.0; incPuppijetmetdphimin   = 0.0;
	Puppijetmumetdphimin    = 0.0; incPuppijetmumetdphimin = 0.0;
	Puppijetelmetdphimin    = 0.0; incPuppijetelmetdphimin = 0.0;
	Puppijetphmetdphimin    = 0.0; incPuppijetphmetdphimin = 0.0;
	
	Puppijetmetdphimin4     = 0.0; incPuppijetmetdphimin4  = 0.0;
	Puppijetmumetdphimin4   = 0.0; incPuppijetmumetdphimin4= 0.0;
	Puppijetelmetdphimin4   = 0.0; incPuppijetelmetdphimin4= 0.0;
	Puppijetphmetdphimin4   = 0.0; incPuppijetphmetdphimin4= 0.0;
	
	// delta phi between Puppijets
	if (centralPuppijetphi.size()>1) 
	  PuppijetPuppijetdphi = deltaPhi(centralPuppijetphi[0], centralPuppijetphi[1]);
	
	std::vector<double> puppijetmetdphiminvector;
	std::vector<double> puppijetmetdphimin4vector;
	// only central Puppijets
	for (size_t i = 0; i < Puppijets.size(); i++) {
	  if (Puppijets[i]->pt() > minPuppiJetPtCount) {
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
	  if (incPuppijets[i]->pt() > minPuppiJetPtCount) {
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
	  if (jets[i]->pt() > minPuppiJetPtCount) {
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
	  if (incjets[i]->pt() > minPuppiJetPtCount) {
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
	  if (jets[i]->pt() > minPuppiJetPtCount) {
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
	  if (incjets[i]->pt() > minPuppiJetPtCount) {
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
	  if (jets[i]->pt() > minPuppiJetPtCount) {
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
	  if (incjets[i]->pt() > minPuppiJetPtCount) {
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
	  if (Puppijets[i]->pt() > minPuppiJetPtCount) {
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


	boostedJetpt    .clear(); boostedJeteta  .clear();  
	boostedJetphi   .clear(); boostedJetm    .clear();
	boostedJetGenpt .clear(); boostedJetGenm .clear();
	boostedJetGeneta .clear(); boostedJetGenphi .clear();

	boostedJettau1  .clear(); boostedJettau2 .clear();
	boostedJettau3  .clear(); boostedJettau4 .clear();
	boostedJetGentau1  .clear(); boostedJetGentau2 .clear();
	boostedJetGentau3  .clear(); boostedJetGentau4 .clear();
	boostedJetecf1  .clear(); boostedJetecf2 .clear(); boostedJetecf3  .clear();
	boostedJetHFlav .clear(); boostedJetPFlav .clear(); boostedJetQGL  .clear(); boostedJetBtag .clear(), boostedJetDoubleBtag .clear();
	
	boostedJetBosonpt .clear(); boostedJetBosoneta .clear(); boostedJetBosonphi .clear(); boostedJetBosonm .clear();
	
	prunedJetpt     .clear(); prunedJetm      .clear(); prunedJetGenpt .clear(); prunedJetGenm  .clear();
	prunedJeteta    .clear(); prunedJetphi     .clear(); prunedJetGeneta .clear(); prunedJetGenphi  .clear();
	prunedJetHFlav  .clear(); prunedJetPFlav  .clear(); prunedJetQGL   .clear(); prunedJetBtag  .clear(); prunedJetDoubleBtag  .clear();
	prunedJetptraw  .clear(); prunedJetmraw .clear();

	softDropJetpt    .clear(); softDropJetm .clear(); softDropJetGenpt .clear(); softDropJetGenm .clear(); 
	softDropJetphi   .clear(); softDropJeteta .clear(); softDropJetGenphi .clear(); softDropJetGeneta .clear(); 
	softDropJetHFlav .clear(); softDropJetPFlav .clear(); softDropJetQGL .clear(); softDropJetBtag .clear(); softDropJetDoubleBtag .clear();
	softDropJetptraw  .clear(); softDropJetmraw .clear();
	
	prunedSubJetpt_1 .clear(); prunedSubJetm_1  .clear(); prunedSubJetphi_1 .clear(); prunedSubJeteta_1 .clear();
	prunedSubJetHFlav_1 .clear(); prunedSubJetQGL_1 .clear(); prunedSubJetBtag_1 .clear();
	prunedSubJetGenpt_1 .clear(); prunedSubJetGenm_1 .clear(); prunedSubJetPFlav_1 .clear();
	prunedSubJetGenphi_1 .clear(); prunedSubJetGeneta_1 .clear(); 
	prunedSubJetptraw_1 .clear(); prunedSubJetmraw_1  .clear(); 
	
	prunedSubJetpt_2 .clear(); prunedSubJetm_2  .clear(); prunedSubJetphi_2 .clear(); prunedSubJeteta_2 .clear();
	prunedSubJetHFlav_2 .clear(); prunedSubJetQGL_2 .clear(); prunedSubJetBtag_2 .clear();  prunedSubJetPFlav_2 .clear();
	prunedSubJetGenpt_2 .clear(); prunedSubJetGenm_2 .clear();
	prunedSubJetGeneta_2 .clear(); prunedSubJetGenphi_2 .clear();
	prunedSubJetptraw_2 .clear(); prunedSubJetmraw_2  .clear(); 

	softDropSubJetpt_1 .clear(); softDropSubJetm_1  .clear(); softDropSubJetphi_1 .clear(); softDropSubJeteta_1 .clear();
	softDropSubJetHFlav_1 .clear(); softDropSubJetQGL_1 .clear(); softDropSubJetBtag_1 .clear(); softDropSubJetPFlav_1 .clear();
	softDropSubJetGenpt_1 .clear(); softDropSubJetGenm_1 .clear();
	softDropSubJetGeneta_1 .clear(); softDropSubJetGenphi_1 .clear();	
	softDropSubJetptraw_1 .clear(); softDropSubJetmraw_1  .clear(); 

	softDropSubJetpt_2 .clear(); softDropSubJetm_2  .clear(); softDropSubJetphi_2 .clear(); softDropSubJeteta_2 .clear();
	softDropSubJetHFlav_2 .clear(); softDropSubJetQGL_2 .clear(); softDropSubJetBtag_2 .clear(); softDropSubJetPFlav_2 .clear();
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

	  // ecf
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf1"))
	    boostedJetecf1 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf1"));
	  
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf2"))
	    boostedJetecf2 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf2"));
	  
	  if(jetsBoosted[i]->hasUserFloat("ecf"+boostedJetsCHSLabel+":ecf3"))
	    boostedJetecf3 .push_back( jetsBoosted[i]->userFloat("ecf"+boostedJetsCHSLabel+":ecf3"));

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

	 	  
	  // matched gen boson
	  if(isMC){
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:pt"))
	      boostedJetBosonpt .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenBosonMatched:pt"));
	    
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:eta"))
	      boostedJetBosoneta .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenBosonMatched:eta"));
	    
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:phi"))
	      boostedJetBosonphi .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenBosonMatched:phi"));
	    
	    if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"GenBosonMatched:mass"))
	      boostedJetBosonm .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"GenBosonMatched:mass"));
	  }

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

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass"))
	    prunedJetmraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawmass"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"))
	    prunedJetptraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"PrunedMatched:rawpt"));

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

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass"))
	    softDropJetmraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:rawmass"));

	  if(jetsBoosted[i]->hasUserFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw"))
	    softDropJetptraw .push_back( jetsBoosted[i]->userFloat(boostedJetsCHSLabel+"SoftDropMatched:ptraw"));

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
	
	boostedPuppiJetpt    .clear(); boostedPuppiJeteta  .clear();  
	boostedPuppiJetphi   .clear(); boostedPuppiJetm    .clear();
	boostedPuppiJetGenpt .clear(); boostedPuppiJetGenm .clear();
	boostedPuppiJetGeneta .clear(); boostedPuppiJetGenphi .clear();
	boostedPuppiJettau1  .clear(); boostedPuppiJettau2 .clear();
	boostedPuppiJettau3  .clear(); boostedPuppiJettau4 .clear();
	boostedPuppiJetGentau1  .clear(); boostedPuppiJetGentau2 .clear();
	boostedPuppiJetGentau3  .clear(); boostedPuppiJetGentau4 .clear();
	boostedPuppiJetecf1  .clear(); boostedPuppiJetecf2 .clear(); boostedPuppiJetecf3  .clear();
	boostedPuppiJetHFlav .clear(); boostedPuppiJetPFlav .clear(); boostedPuppiJetQGL  .clear(); 
	boostedPuppiJetBtag .clear(); boostedPuppiJetDoubleBtag .clear();
	
	boostedPuppiJetBosonpt .clear(); boostedPuppiJetBosoneta .clear(); boostedPuppiJetBosonphi .clear(); boostedPuppiJetBosonm .clear();
	
	prunedPuppiJetpt     .clear(); prunedPuppiJetm      .clear(); prunedPuppiJetGenpt .clear(); prunedPuppiJetGenm  .clear();
	prunedPuppiJeteta     .clear(); prunedPuppiJetphi   .clear(); prunedPuppiJetGeneta .clear(); prunedPuppiJetGenphi  .clear();
	prunedPuppiJetHFlav  .clear(); prunedPuppiJetPFlav  .clear(); prunedPuppiJetQGL   .clear(); prunedPuppiJetBtag  .clear();
	prunedPuppiJetDoubleBtag  .clear(); prunedPuppiJetmraw .clear(); prunedPuppiJetptraw .clear();
	

	softDropPuppiJetpt    .clear(); softDropPuppiJetm .clear(); softDropPuppiJetGenpt .clear(); softDropPuppiJetGenm .clear(); 
	softDropPuppiJeteta    .clear(); softDropPuppiJetphi .clear(); softDropPuppiJetGeneta .clear(); softDropPuppiJetGenphi .clear(); 
	softDropPuppiJetHFlav .clear(); softDropPuppiJetPFlav .clear(); softDropPuppiJetQGL .clear(); softDropPuppiJetBtag .clear();
	softDropPuppiJetDoubleBtag .clear(); softDropPuppiJetmraw .clear(); softDropPuppiJetptraw .clear();
	
	prunedPuppiSubJetpt_1 .clear(); prunedPuppiSubJetm_1  .clear(); prunedPuppiSubJetphi_1 .clear(); prunedPuppiSubJeteta_1 .clear();
	prunedPuppiSubJetHFlav_1 .clear(); prunedPuppiSubJetPFlav_1 .clear(); prunedPuppiSubJetQGL_1 .clear(); prunedPuppiSubJetBtag_1 .clear();
	prunedPuppiSubJetGenpt_1 .clear(); prunedPuppiSubJetGenm_1 .clear();
	prunedPuppiSubJetGenphi_1 .clear(); prunedPuppiSubJetGeneta_1 .clear();
	prunedPuppiSubJetptraw_1 .clear(); prunedPuppiSubJetmraw_1  .clear(); 

	prunedPuppiSubJetpt_2 .clear(); prunedPuppiSubJetm_2  .clear(); prunedPuppiSubJetphi_2 .clear(); prunedPuppiSubJeteta_2 .clear();
	prunedPuppiSubJetHFlav_2 .clear();prunedPuppiSubJetPFlav_2 .clear(); prunedPuppiSubJetQGL_2 .clear(); prunedPuppiSubJetBtag_2 .clear();
	prunedPuppiSubJetGenpt_2 .clear(); prunedPuppiSubJetGenm_2 .clear();
	prunedPuppiSubJetGenphi_2 .clear(); prunedPuppiSubJetGeneta_2 .clear();
	prunedPuppiSubJetptraw_2 .clear(); prunedPuppiSubJetmraw_2  .clear(); 
	
	softDropPuppiSubJetpt_1 .clear(); softDropPuppiSubJetm_1  .clear(); softDropPuppiSubJetphi_1 .clear(); softDropPuppiSubJeteta_1 .clear();
	softDropPuppiSubJetHFlav_1 .clear(); softDropPuppiSubJetQGL_1 .clear(); softDropPuppiSubJetBtag_1 .clear(); softDropPuppiSubJetPFlav_1 .clear();
	softDropPuppiSubJetGenpt_1 .clear(); softDropPuppiSubJetGenm_1 .clear();
	softDropPuppiSubJetGenphi_1 .clear(); softDropPuppiSubJetGeneta_1 .clear();
	softDropPuppiSubJetptraw_1 .clear(); softDropPuppiSubJetmraw_1  .clear(); 

	softDropPuppiSubJetpt_2 .clear(); softDropPuppiSubJetm_2  .clear(); softDropPuppiSubJetphi_2 .clear(); softDropPuppiSubJeteta_2 .clear();
	softDropPuppiSubJetHFlav_2 .clear(); softDropPuppiSubJetQGL_2 .clear(); softDropPuppiSubJetBtag_2 .clear(); softDropPuppiSubJetPFlav_2 .clear();
	softDropPuppiSubJetGenpt_2 .clear(); softDropPuppiSubJetGenm_2 .clear();
	softDropPuppiSubJetGenphi_2 .clear(); softDropPuppiSubJetGeneta_2 .clear();
	softDropPuppiSubJetptraw_2 .clear(); softDropPuppiSubJetmraw_2  .clear(); 
	
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
	  	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf1")){
	    boostedPuppiJetecf1 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf1"));
	  }
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf2")){
	    boostedPuppiJetecf2 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf2"));
	  }
	  
	  if(puppiJetsBoosted[i]->hasUserFloat("ecf"+boostedJetsPuppiLabel+":ecf3")){
	    boostedPuppiJetecf3 .push_back( puppiJetsBoosted[i]->userFloat("ecf"+boostedJetsPuppiLabel+":ecf3"));
	  }
	  
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
	  
	  // matched gen boson	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:pt"))
	    boostedPuppiJetBosonpt .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:pt"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:eta"))
	    boostedPuppiJetBosoneta .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:eta"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:phi"))
	    boostedPuppiJetBosonphi .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:phi"));
	  
	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"GenBosonMatched:mass"))
	    boostedPuppiJetBosonm .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"GenBosonMatched:mass"));
	  
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

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass"))
	    prunedPuppiJetmraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawmass"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"))
	    prunedPuppiJetptraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"PrunedMatched:rawpt"));
	  
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

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass"))
	    softDropPuppiJetmraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawmass"));

	  if(puppiJetsBoosted[i]->hasUserFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"))
	    softDropPuppiJetptraw .push_back( puppiJetsBoosted[i]->userFloat(boostedJetsPuppiLabel+"SoftDropMatched:rawpt"));
	  
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
    wzid_h        = 0; wzmass_h      = 0.0; wzpt_h        = 0.0; wzeta_h       = 0.0; wzphi_h       = 0.0;
    topmass       = 0; toppt         = 0.0; topeta        = 0.0; topphi        = 0.0;
    atopmass      = 0; atoppt        = 0.0; atopeta       = 0.0; atopphi       = 0.0;
    q1id          = 0; q1pt          = 0.0; q1eta         = 0.0; q1phi         = 0.0;
    q2id          = 0; q2pt          = 0.0; q2eta         = 0.0; q2phi         = 0.0;
    parid         = 0; parpt         = 0.0; pareta        = 0.0; parphi        = 0.0;
    ancid         = 0; ancpt         = 0.0; anceta        = 0.0; ancphi        = 0.0;

    dmmass   = 0.; dmphi   = 0.; dmeta   = 0.; dmpt   = 0.; dmid   = 0;
    dmX1mass = 0.; dmX1phi = 0.; dmX1eta = 0.; dmX1pt = 0.; dmX1id = 0;
    dmX2mass = 0.; dmX2phi = 0.; dmX2eta = 0.; dmX2pt = 0.; dmX2id = 0;

    if (isWorZorSignalMCSample && gensH.isValid()) {

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

      // loop on genParticles (prunedGenParticles) trying to find W/Z decying leptonically or hadronically, top and anti-top quarks
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	if ( (gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
	     gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
	     abs(gens_iter->daughter(0)->pdgId()) > 10 && 
	     abs(gens_iter->daughter(0)->pdgId()) < 17)  { // decays into leptons, neutrinos and quarks
	  
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
		    (abs(gens_iter->daughter(1)->pdgId()) > 0 && abs(gens_iter->daughter(1)->pdgId()) <= 5)))  { // decays into leptons, neutrinos and quarks
	  
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
	}
      }
    }

    tree->Fill();    
  
} 


void MonoJetTreeMaker::beginJob() {

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree"       , "tree");

  // Run, Lumi, Event info
  tree->Branch("event"                , &event                , "event/i");
  tree->Branch("run"                  , &run                  , "run/i");
  tree->Branch("lumi"                 , &lumi                 , "lumi/i");
  // Event weights
  tree->Branch("xsec"                 , &xsec                 , "xsec/D");
  tree->Branch("wgt"                  , &wgt                  , "wgt/D");
  tree->Branch("pswgt"                , &pswgt                , "pswgt/D");
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
  tree->Branch("hltphoton120"         , &hltphoton120         , "hltphoton120/b");
  tree->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/b");
  tree->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/b");
  tree->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/b");
  tree->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/b");
  tree->Branch("hltelnoiso"           , &hltelnoiso           , "hltelnoiso/b");

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
  tree->Branch("njetsinc"             , &njetsinc             , "njetsinc/i");
  tree->Branch("nbjets"               , &nbjets               , "nbjets/i");
  tree->Branch("nbjetslowpt"          , &nbjetslowpt          , "nbjetslowpt/i");

  if(addPuppiJets){
    tree->Branch("npuppijets"                , &npuppijets                , "npuppijets/i");
    tree->Branch("npuppijetsinc"             , &npuppijetsinc             , "npuppijetsinc/i");
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
  tree->Branch("genmetphi",&genmetphi,"genmetphi/D");

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

  tree->Branch("centraljetpt", "std::vector<double>", &centraljetpt);
  tree->Branch("centraljeteta", "std::vector<double>", &centraljeteta);
  tree->Branch("centraljetphi", "std::vector<double>", &centraljetphi);
  tree->Branch("centraljetm", "std::vector<double>", &centraljetm);
  tree->Branch("centraljetbtag", "std::vector<double>", &centraljetbtag);
  tree->Branch("centraljetCHfrac", "std::vector<double>", &centraljetCHfrac);
  tree->Branch("centraljetNHfrac", "std::vector<double>", &centraljetNHfrac);
  tree->Branch("centraljetEMfrac", "std::vector<double>", &centraljetEMfrac);
  tree->Branch("centraljetCEMfrac", "std::vector<double>", &centraljetCEMfrac);
  tree->Branch("centraljetmetdphi", "std::vector<double>", &centraljetmetdphi);
  tree->Branch("centraljetHFlav", "std::vector<double>", &centraljetHFlav);
  tree->Branch("centraljetPFlav", "std::vector<double>", &centraljetPFlav);
  tree->Branch("centraljetQGL", "std::vector<double>", &centraljetQGL);
  tree->Branch("centraljetPUID", "std::vector<double>", &centraljetPUID);
  tree->Branch("centraljetGenpt", "std::vector<double>", &centraljetGenpt);
  tree->Branch("centraljetGeneta", "std::vector<double>", &centraljetGeneta);
  tree->Branch("centraljetGenphi", "std::vector<double>", &centraljetGenphi);
  tree->Branch("centraljetGenm", "std::vector<double>", &centraljetGenm);
  tree->Branch("centraljetBtagSF", "std::vector<double>", &centraljetBtagSF);
  tree->Branch("centraljetBtagSFUp", "std::vector<double>", &centraljetBtagSFUp);
  tree->Branch("centraljetBtagSFDown", "std::vector<double>", &centraljetBtagSFDown);

  tree->Branch("forwardjetpt", "std::vector<double>", &forwardjetpt);
  tree->Branch("forwardjeteta", "std::vector<double>", &forwardjeteta);
  tree->Branch("forwardjetphi", "std::vector<double>", &forwardjetphi);
  tree->Branch("forwardjetm", "std::vector<double>", &forwardjetm);
  tree->Branch("forwardjetbtag", "std::vector<double>", &forwardjetbtag);
  tree->Branch("forwardjetCHfrac", "std::vector<double>", &forwardjetCHfrac);
  tree->Branch("forwardjetNHfrac", "std::vector<double>", &forwardjetNHfrac);
  tree->Branch("forwardjetEMfrac", "std::vector<double>", &forwardjetEMfrac);
  tree->Branch("forwardjetCEMfrac", "std::vector<double>", &forwardjetCEMfrac);
  tree->Branch("forwardjetmetdphi", "std::vector<double>", &forwardjetmetdphi);
  tree->Branch("forwardjetHFlav", "std::vector<double>", &forwardjetHFlav);
  tree->Branch("forwardjetPFlav", "std::vector<double>", &forwardjetPFlav);
  tree->Branch("forwardjetQGL", "std::vector<double>", &forwardjetQGL);
  tree->Branch("forwardjetPUID", "std::vector<double>", &forwardjetPUID);
  tree->Branch("forwardjetGenpt", "std::vector<double>", &forwardjetGenpt);
  tree->Branch("forwardjetGeneta", "std::vector<double>", &forwardjetGeneta);
  tree->Branch("forwardjetGenphi", "std::vector<double>", &forwardjetGenphi);
  tree->Branch("forwardjetGenm", "std::vector<double>", &forwardjetGenm);

  tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");

  tree->Branch("jetmetdphimin"        , &jetmetdphimin        , "jetmetdphimin/D");
  tree->Branch("jetmumetdphimin"      , &jetmumetdphimin      , "jetmumetdphimin/D");
  tree->Branch("jetelmetdphimin"      , &jetelmetdphimin      , "jetelmetdphimin/D");
  tree->Branch("jetphmetdphimin"      , &jetphmetdphimin      , "jetphmetdphimin/D");

  tree->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/D");
  tree->Branch("incjetmumetdphimin"   , &incjetmumetdphimin   , "incjetmumetdphimin/D");
  tree->Branch("incjetelmetdphimin"   , &incjetelmetdphimin   , "incjetelmetdphimin/D");
  tree->Branch("incjetphmetdphimin"   , &incjetphmetdphimin   , "incjetphmetdphimin/D");

  tree->Branch("alljetmetdphimin"     , &alljetmetdphimin     , "alljetmetdphimin/D");
  tree->Branch("alljetmumetdphimin"   , &alljetmumetdphimin   , "alljetmumetdphimin/D");
  tree->Branch("alljetelmetdphimin"   , &alljetelmetdphimin   , "alljetelmetdphimin/D");
  tree->Branch("alljetphmetdphimin"   , &alljetphmetdphimin   , "alljetphmetdphimin/D");

  tree->Branch("jetmetdphimin4"       , &jetmetdphimin4       , "jetmetdphimin4/D");
  tree->Branch("jetmumetdphimin4"     , &jetmumetdphimin4     , "jetmumetdphimin4/D");
  tree->Branch("jetelmetdphimin4"     , &jetelmetdphimin4     , "jetelmetdphimin4/D");
  tree->Branch("jetphmetdphimin4"     , &jetphmetdphimin4     , "jetphmetdphimin4/D");

  tree->Branch("incjetmetdphimin4"    , &incjetmetdphimin4    , "incjetmetdphimin4/D");
  tree->Branch("incjetmumetdphimin4"  , &incjetmumetdphimin4  , "incjetmumetdphimin4/D");
  tree->Branch("incjetelmetdphimin4"  , &incjetelmetdphimin4  , "incjetelmetdphimin4/D");
  tree->Branch("incjetphmetdphimin4"  , &incjetphmetdphimin4  , "incjetphmetdphimin4/D");

  tree->Branch("alljetmetdphimin4"    , &alljetmetdphimin4    , "alljetmetdphimin4/D");
  tree->Branch("alljetmumetdphimin4"  , &alljetmumetdphimin4  , "alljetmumetdphimin4/D");
  tree->Branch("alljetelmetdphimin4"  , &alljetelmetdphimin4  , "alljetelmetdphimin4/D");
  tree->Branch("alljetphmetdphimin4"  , &alljetphmetdphimin4  , "alljetphmetdphimin4/D");

  tree->Branch("ht"                   , &ht                   , "ht/D");

  if(addPuppiJets){

    tree->Branch("centralPuppijetpt", "std::vector<double>", &centralPuppijetpt);
    tree->Branch("centralPuppijeteta", "std::vector<double>", &centralPuppijeteta);
    tree->Branch("centralPuppijetphi", "std::vector<double>", &centralPuppijetphi);
    tree->Branch("centralPuppijetm", "std::vector<double>", &centralPuppijetm);
    tree->Branch("centralPuppijetbtag", "std::vector<double>", &centralPuppijetbtag);
    tree->Branch("centralPuppijetCHfrac", "std::vector<double>", &centralPuppijetCHfrac);
    tree->Branch("centralPuppijetNHfrac", "std::vector<double>", &centralPuppijetNHfrac);
    tree->Branch("centralPuppijetEMfrac", "std::vector<double>", &centralPuppijetEMfrac);
    tree->Branch("centralPuppijetCEMfrac", "std::vector<double>", &centralPuppijetCEMfrac);
    tree->Branch("centralPuppijetmetdphi", "std::vector<double>", &centralPuppijetmetdphi);
    tree->Branch("centralPuppijetHFlav", "std::vector<double>", &centralPuppijetHFlav);
    tree->Branch("centralPuppijetPFlav", "std::vector<double>", &centralPuppijetPFlav);
    tree->Branch("centralPuppijetQGL", "std::vector<double>", &centralPuppijetQGL);
    tree->Branch("centralPuppijetPUID", "std::vector<double>", &centralPuppijetPUID);
    tree->Branch("centralPuppijetGenpt", "std::vector<double>", &centralPuppijetGenpt);
    tree->Branch("centralPuppijetGeneta", "std::vector<double>", &centralPuppijetGeneta);
    tree->Branch("centralPuppijetGenphi", "std::vector<double>", &centralPuppijetGenphi);
    tree->Branch("centralPuppijetGenm", "std::vector<double>", &centralPuppijetGenm);
    tree->Branch("centralPuppijetBtagSF", "std::vector<double>", &centralPuppijetBtagSF);
    tree->Branch("centralPuppijetBtagSFUp", "std::vector<double>", &centralPuppijetBtagSFUp);
    tree->Branch("centralPuppijetBtagSFDown", "std::vector<double>", &centralPuppijetBtagSFDown);

    tree->Branch("forwardPuppijetpt", "std::vector<double>", &forwardPuppijetpt);
    tree->Branch("forwardPuppijeteta", "std::vector<double>", &forwardPuppijeteta);
    tree->Branch("forwardPuppijetphi", "std::vector<double>", &forwardPuppijetphi);
    tree->Branch("forwardPuppijetm", "std::vector<double>", &forwardPuppijetm);
    tree->Branch("forwardPuppijetbtag", "std::vector<double>", &forwardPuppijetbtag);
    tree->Branch("forwardPuppijetCHfrac", "std::vector<double>", &forwardPuppijetCHfrac);
    tree->Branch("forwardPuppijetNHfrac", "std::vector<double>", &forwardPuppijetNHfrac);
    tree->Branch("forwardPuppijetEMfrac", "std::vector<double>", &forwardPuppijetEMfrac);
    tree->Branch("forwardPuppijetCEMfrac", "std::vector<double>", &forwardPuppijetCEMfrac);
    tree->Branch("forwardPuppijetmetdphi", "std::vector<double>", &forwardPuppijetmetdphi);
    tree->Branch("forwardPuppijetHFlav", "std::vector<double>", &forwardPuppijetHFlav);
    tree->Branch("forwardPuppijetPFlav", "std::vector<double>", &forwardPuppijetPFlav);
    tree->Branch("forwardPuppijetQGL", "std::vector<double>", &forwardPuppijetQGL);
    tree->Branch("forwardPuppijetPUID", "std::vector<double>", &forwardPuppijetPUID);
    tree->Branch("forwardPuppijetGenpt", "std::vector<double>", &forwardPuppijetGenpt);
    tree->Branch("forwardPuppijetGeneta", "std::vector<double>", &forwardPuppijetGeneta);
    tree->Branch("forwardPuppijetGenphi", "std::vector<double>", &forwardPuppijetGenphi);
    tree->Branch("forwardPuppijetGenm", "std::vector<double>", &forwardPuppijetGenm);
    
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

  // Top infor
  tree->Branch("topmass"               , &topmass               , "topmass/D");
  tree->Branch("toppt"                 , &toppt                 , "toppt/D");
  tree->Branch("topeta"                , &topeta                , "topeta/D");
  tree->Branch("topphi"                , &topphi                , "topphi/D");

  tree->Branch("atopmass"               , &atopmass               , "atopmass/D");
  tree->Branch("atoppt"                 , &atoppt                 , "atoppt/D");
  tree->Branch("atopeta"                , &atopeta                , "atopeta/D");
  tree->Branch("atopphi"                , &atopphi                , "atopphi/D");



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

  // AK8 Puppi jets                                                                                                                                                             
  if(addSubstructureCHS){

    tree->Branch("boostedJetpt", "std::vector<double>", &boostedJetpt);
    tree->Branch("boostedJeteta", "std::vector<double>", &boostedJeteta);
    tree->Branch("boostedJetphi", "std::vector<double>", &boostedJetphi);
    tree->Branch("boostedJetm", "std::vector<double>", &boostedJetm);

    tree->Branch("boostedJetGenpt", "std::vector<double>", &boostedJetGenpt);
    tree->Branch("boostedJetGeneta", "std::vector<double>", &boostedJetGeneta);
    tree->Branch("boostedJetGenphi", "std::vector<double>", &boostedJetGenphi);
    tree->Branch("boostedJetGenm", "std::vector<double>", &boostedJetGenm);

    tree->Branch("boostedJetHFlav", "std::vector<double>", &boostedJetHFlav);
    tree->Branch("boostedJetPFlav", "std::vector<double>", &boostedJetPFlav);
    tree->Branch("boostedJetQGL", "std::vector<double>", &boostedJetQGL);
    tree->Branch("boostedJetBtag", "std::vector<double>", &boostedJetBtag);
    tree->Branch("boostedJetDoubleBtag", "std::vector<double>", &boostedJetDoubleBtag);

    tree->Branch("boostedJettau1", "std::vector<double>", &boostedJettau1);
    tree->Branch("boostedJettau2", "std::vector<double>", &boostedJettau2);
    tree->Branch("boostedJettau3", "std::vector<double>", &boostedJettau3);
    tree->Branch("boostedJettau4", "std::vector<double>", &boostedJettau4);

    tree->Branch("boostedJetGentau1", "std::vector<double>", &boostedJetGentau1);
    tree->Branch("boostedJetGentau2", "std::vector<double>", &boostedJetGentau2);
    tree->Branch("boostedJetGentau3", "std::vector<double>", &boostedJetGentau3);
    tree->Branch("boostedJetGentau4", "std::vector<double>", &boostedJetGentau4);

    tree->Branch("boostedJetecf1", "std::vector<double>", &boostedJetecf1);
    tree->Branch("boostedJetecf2", "std::vector<double>", &boostedJetecf2);
    tree->Branch("boostedJetecf3", "std::vector<double>", &boostedJetecf3);

    tree->Branch("boostedJetBosonpt", "std::vector<double>", &boostedJetBosonpt);
    tree->Branch("boostedJetBosoneta", "std::vector<double>", &boostedJetBosoneta);
    tree->Branch("boostedJetBosonphi", "std::vector<double>", &boostedJetBosonphi);
    tree->Branch("boostedJetBosonm", "std::vector<double>", &boostedJetBosonm);


    tree->Branch("prunedJetpt", "std::vector<double>", &prunedJetpt);
    tree->Branch("prunedJeteta", "std::vector<double>", &prunedJeteta);
    tree->Branch("prunedJetphi", "std::vector<double>", &prunedJetphi);
    tree->Branch("prunedJetm", "std::vector<double>", &prunedJetm);
    tree->Branch("prunedJetptraw", "std::vector<double>", &prunedJetptraw);
    tree->Branch("prunedJetmraw", "std::vector<double>", &prunedJetmraw);

    tree->Branch("prunedJetGenpt", "std::vector<double>", &prunedJetGenpt);
    tree->Branch("prunedJetGeneta", "std::vector<double>", &prunedJetGeneta);
    tree->Branch("prunedJetGenphi", "std::vector<double>", &prunedJetGenphi);
    tree->Branch("prunedJetGenm", "std::vector<double>", &prunedJetGenm);

    tree->Branch("prunedJetHFlav", "std::vector<double>", &prunedJetHFlav);
    tree->Branch("prunedJetPFlav", "std::vector<double>", &prunedJetPFlav);
    tree->Branch("prunedJetQGL", "std::vector<double>", &prunedJetQGL);
    tree->Branch("prunedJetBtag", "std::vector<double>", &prunedJetBtag);
    tree->Branch("prunedJetDoubleBtag", "std::vector<double>", &prunedJetDoubleBtag);


    tree->Branch("softDropJetpt", "std::vector<double>", &softDropJetpt);
    tree->Branch("softDropJeteta", "std::vector<double>", &softDropJeteta);
    tree->Branch("softDropJetphi", "std::vector<double>", &softDropJetphi);
    tree->Branch("softDropJetm", "std::vector<double>", &softDropJetm);
    tree->Branch("softDropJetptraw", "std::vector<double>", &softDropJetptraw);
    tree->Branch("softDropJetmraw", "std::vector<double>", &softDropJetmraw);

    tree->Branch("softDropJetGenpt", "std::vector<double>", &softDropJetGenpt);
    tree->Branch("softDropJetGeneta", "std::vector<double>", &softDropJetGeneta);
    tree->Branch("softDropJetGenphi", "std::vector<double>", &softDropJetGenphi);
    tree->Branch("softDropJetGenm", "std::vector<double>", &softDropJetGenm);

    tree->Branch("softDropJetHFlav", "std::vector<double>", &softDropJetHFlav);
    tree->Branch("softDropJetPFlav", "std::vector<double>", &softDropJetPFlav);
    tree->Branch("softDropJetQGL", "std::vector<double>", &softDropJetQGL);
    tree->Branch("softDropJetBtag", "std::vector<double>", &softDropJetBtag);
    tree->Branch("softDropJetDoubleBtag", "std::vector<double>", &softDropJetDoubleBtag);

    tree->Branch("prunedSubJetpt_1","std::vector<double>",  &prunedSubJetpt_1);
    tree->Branch("prunedSubJeteta_1","std::vector<double>",  &prunedSubJeteta_1);
    tree->Branch("prunedSubJetphi_1","std::vector<double>",  &prunedSubJetphi_1);
    tree->Branch("prunedSubJetm_1", "std::vector<double>", &prunedSubJetm_1);
    tree->Branch("prunedSubJetGenpt_1","std::vector<double>",  &prunedSubJetGenpt_1);
    tree->Branch("prunedSubJetGenm_1", "std::vector<double>", &prunedSubJetGenm_1);
    tree->Branch("prunedSubJetGeneta_1", "std::vector<double>", &prunedSubJetGeneta_1);
    tree->Branch("prunedSubJetGenphi_1", "std::vector<double>", &prunedSubJetGenphi_1);
    tree->Branch("prunedSubJetHFlav_1", "std::vector<double>", &prunedSubJetHFlav_1);
    tree->Branch("prunedSubJetPFlav_1", "std::vector<double>", &prunedSubJetPFlav_1);
    tree->Branch("prunedSubJetQGL_1", "std::vector<double>", &prunedSubJetQGL_1);
    tree->Branch("prunedSubJetBtag_1", "std::vector<double>", &prunedSubJetBtag_1);
    tree->Branch("prunedSubJetptraw_1", "std::vector<double>", &prunedSubJetptraw_1);
    tree->Branch("prunedSubJetmraw_1", "std::vector<double>", &prunedSubJetmraw_1);

    tree->Branch("prunedSubJetpt_2","std::vector<double>",  &prunedSubJetpt_2);
    tree->Branch("prunedSubJeteta_2","std::vector<double>",  &prunedSubJeteta_2);
    tree->Branch("prunedSubJetphi_2","std::vector<double>",  &prunedSubJetphi_2);
    tree->Branch("prunedSubJetm_2", "std::vector<double>", &prunedSubJetm_2);
    tree->Branch("prunedSubJetGenpt_2","std::vector<double>",  &prunedSubJetGenpt_2);
    tree->Branch("prunedSubJetGenm_2", "std::vector<double>", &prunedSubJetGenm_2);
    tree->Branch("prunedSubJetGeneta_2", "std::vector<double>", &prunedSubJetGeneta_2);
    tree->Branch("prunedSubJetGenphi_2", "std::vector<double>", &prunedSubJetGenphi_2);
    tree->Branch("prunedSubJetHFlav_2", "std::vector<double>", &prunedSubJetHFlav_2);
    tree->Branch("prunedSubJetPFlav_2", "std::vector<double>", &prunedSubJetPFlav_2);
    tree->Branch("prunedSubJetQGL_2", "std::vector<double>", &prunedSubJetQGL_2);
    tree->Branch("prunedSubJetBtag_2", "std::vector<double>", &prunedSubJetBtag_2);
    tree->Branch("prunedSubJetptraw_2", "std::vector<double>", &prunedSubJetptraw_2);
    tree->Branch("prunedSubJetmraw_2", "std::vector<double>", &prunedSubJetmraw_2);


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
  }

  if(addSubstructurePuppi){

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

    tree->Branch("boostedPuppiJetecf1", "std::vector<double>", &boostedPuppiJetecf1);
    tree->Branch("boostedPuppiJetecf2", "std::vector<double>", &boostedPuppiJetecf2);
    tree->Branch("boostedPuppiJetecf3", "std::vector<double>", &boostedPuppiJetecf3);

    tree->Branch("boostedPuppiJetBosonpt", "std::vector<double>", &boostedPuppiJetBosonpt);
    tree->Branch("boostedPuppiJetBosoneta", "std::vector<double>", &boostedPuppiJetBosoneta);
    tree->Branch("boostedPuppiJetBosonphi", "std::vector<double>", &boostedPuppiJetBosonphi);
    tree->Branch("boostedPuppiJetBosonm", "std::vector<double>", &boostedPuppiJetBosonm);


    tree->Branch("prunedPuppiJetpt", "std::vector<double>", &prunedPuppiJetpt);
    tree->Branch("prunedPuppiJeteta", "std::vector<double>", &prunedPuppiJeteta);
    tree->Branch("prunedPuppiJetphi", "std::vector<double>", &prunedPuppiJetphi);
    tree->Branch("prunedPuppiJetm", "std::vector<double>", &prunedPuppiJetm);
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
  triggerPathsVector.push_back("HLT_Photon120_v");
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
  triggerPathsVector.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  triggerPathsVector.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  
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


  // info about MC cross section in case the xsec parsed has a dummy value
  if(isMC){    
    edm::Handle<LHERunInfoProduct> run;
    iRun.getByLabel(lheRunTag,run);
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    if(xsec < 0)
      xsec = myLHERunInfoProduct.heprup().XSECUP.at(0);
    

    using namespace boost::algorithm;

    if(isWorZorSignalMCSample){

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
    if (jetabseta >= 0.00 && jetabseta < 2.00 && puidval > -0.826) passpuid = true;
    if (jetabseta >= 2.00 && jetabseta < 2.50 && puidval > -0.813) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 3.00 && puidval > -0.574) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.363) passpuid = true;
  }
  else if(level == "medium"){ // left as done in EXO-15-003 for the time being
    if (jetabseta >= 0.00 && jetabseta < 2.50 && puidval > -0.63) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 2.75 && puidval > -0.60) passpuid = true;
    if (jetabseta >= 2.75 && jetabseta < 3.00 && puidval > -0.55) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.45) passpuid = true;
  }
  else if(level == "tight"){
    if (jetabseta >= 0.00 && jetabseta < 2.00 && puidval >  0.291) passpuid = true;
    if (jetabseta >= 2.00 && jetabseta < 2.50 && puidval > -0.306) passpuid = true;
    if (jetabseta >= 2.50 && jetabseta < 3.00 && puidval > -0.369) passpuid = true;
    if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.247) passpuid = true;    

  }
  return passpuid;
}



void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MonoJetTreeMaker);

//  LocalWords:  TypeI jetSorter
