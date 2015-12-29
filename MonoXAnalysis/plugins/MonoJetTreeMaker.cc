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

class MonoJetTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit MonoJetTreeMaker(const edm::ParameterSet&);
        ~MonoJetTreeMaker();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

        void findMother(const reco::Candidate*, int &, double &, double &, double &);
        void findFirstNonPhotonMother(const reco::Candidate*, int &, double &, double &, double &);
        double computeMuonIso(const reco::Muon&);

        // InputTags
        edm::InputTag triggerResultsTag;
        edm::InputTag filterResultsTag;

        // Tokens
        edm::EDGetTokenT<edm::TriggerResults>              triggerResultsToken;
        edm::EDGetTokenT<edm::TriggerResults>              filterResultsToken;
        edm::EDGetTokenT<bool>                             hbhelooseToken;
        edm::EDGetTokenT<bool>                             hbhetightToken;
        edm::EDGetTokenT<bool>                             hbheisoToken;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo> >  pileupInfoToken;
        edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
        edm::EDGetTokenT<LHEEventProduct>                  lheInfoToken;
        edm::EDGetTokenT<std::vector<reco::Vertex> >       verticesToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
        edm::EDGetTokenT<pat::MuonRefVector>               muonsToken;
        edm::EDGetTokenT<pat::ElectronRefVector>           electronsToken;
        edm::EDGetTokenT<pat::PhotonRefVector>             photonsToken;
        edm::EDGetTokenT<pat::MuonRefVector>               tightmuonsToken;
        edm::EDGetTokenT<pat::ElectronRefVector>           tightelectronsToken;
        edm::EDGetTokenT<pat::PhotonRefVector>             tightphotonsToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             electronLooseIdToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonLooseIdToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonMediumIdToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonTightIdToken;
        edm::EDGetTokenT<edm::ValueMap<bool> >             photonHighPtIdToken;
        edm::EDGetTokenT<std::vector<pat::Tau> >           tausToken;
        edm::EDGetTokenT<std::vector<pat::Jet> >           jetsToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            pfmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            t1pfmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            mumetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            t1mumetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            elmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            t1elmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            phmetToken;
        edm::EDGetTokenT<edm::View<reco::MET> >            t1phmetToken;

        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        std::vector<std::string> filterPathsVector;
        std::map<std::string, int> filterPathsMap;

        bool applyHLTFilter;
        bool isWorZMCSample, isSignalSample;   
        bool cleanMuonJet, cleanElectronJet, cleanPhotonJet;   
        bool uselheweights;   
        bool addqcdpdfweights;   
        TTree* tree;

        int32_t  puobs, putrue; 
        int32_t  wzid, l1id, l2id, mu1pid, mu2pid, mu1id, mu2id, mu1idm, mu2idm, mu1idt, mu2idt, el1pid, el2pid, el1id, el1idl, el2id, el2idl, phidl, phidm, phidt, phidh, parid, ancid; 
        uint32_t event, run, lumi;
        uint32_t nvtx, nmuons, nelectrons, ntaus, ntightmuons, ntightelectrons, nphotons, njets, nbjets, nbjetslowpt;
        uint8_t  hltmet90, hltmet120, hltmetwithmu90, hltmetwithmu120, hltmetwithmu170, hltmetwithmu300, hltjetmet90, hltjetmet120, hltphoton165, hltphoton175, hltdoublemu, hltsinglemu, hltdoubleel, hltsingleel;
        uint8_t  flagcsctight, flaghbhenoise, flaghbheloose, flaghbhetight, flaghbheiso, flageebadsc;
        double   pfmet, pfmetphi, t1pfmet, t1pfmetphi, mumet, mumetphi, t1mumet, t1mumetphi, elmet, elmetphi, t1elmet, t1elmetphi, phmet, phmetphi, t1phmet, t1phmetphi;
        double   hmet, hmetphi, amet, ametphi, bmet, bmetphi, cmet, cmetphi, emet, emetphi, mmet, mmetphi, pmet, pmetphi, omet, ometphi;
        double   leadingjetpt, leadingjeteta, leadingjetphi;
        double   signaljetpt , signaljeteta , signaljetphi , signaljetbtag, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetCEMfrac, signaljetmetdphi;
        double   secondjetpt , secondjeteta , secondjetphi , secondjetbtag, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetCEMfrac, secondjetmetdphi;
        double   thirdjetpt  , thirdjeteta  , thirdjetphi  , thirdjetbtag , thirdjetCHfrac , thirdjetNHfrac , thirdjetEMfrac , thirdjetCEMfrac , thirdjetmetdphi ;
        double   fourthjetpt , fourthjeteta , fourthjetphi , fourthjetbtag, fourthjetCHfrac, fourthjetNHfrac, fourthjetEMfrac, fourthjetCEMfrac, fourthjetmetdphi;
        double   jetmetdphimin , incjetmetdphimin , jetelmetdphimin , incjetelmetdphimin , jetphmetdphimin , incjetphmetdphimin , jetjetdphi;
        double   jetmetdphimin4, incjetmetdphimin4, jetelmetdphimin4, incjetelmetdphimin4, jetphmetdphimin4, incjetphmetdphimin4, ht; 
        double   wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi, parpt, pareta, parphi, ancpt, anceta, ancphi;
        double   mu1pt, mu1eta, mu1phi, mu1pfpt, mu1pfeta, mu1pfphi, mu1iso, mu2pt, mu2eta, mu2phi, mu2pfpt, mu2pfeta, mu2pfphi, mu2iso;
        double   el1pt, el1eta, el1phi, el2pt, el2eta, el2phi, phpt, pheta, phphi;
        double   zmass, zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt;
        double   xsec, wgt, kfact, puwgt;
        double*  wgtpdf;
        double*  wgtqcd;

        struct PatJetPtSorter {
            bool operator() (pat::JetRef i, pat::JetRef j) {
                return (i->pt() > j->pt());
            }
        } jetsorter;
        
        struct PatMuonPtSorter {
            bool operator() (pat::MuonRef i, pat::MuonRef j) {
                return (i->pt() > j->pt());
            }
        } muonsorter;
        
        struct PatElectronPtSorter {
            bool operator() (pat::ElectronRef i, pat::ElectronRef j) {
                return (i->pt() > j->pt());
            }
        } electronsorter;

        struct PatPhotonPtSorter {
            bool operator() (pat::PhotonRef i, pat::PhotonRef j) {
                return (i->pt() > j->pt());
            }
        } photonsorter;

};

MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    filterResultsTag(iConfig.getParameter<edm::InputTag>("filterResults")),
    triggerResultsToken(consumes<edm::TriggerResults> (triggerResultsTag)),
    filterResultsToken(consumes<edm::TriggerResults> (filterResultsTag)),
    hbhelooseToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbheloose"))),
    hbhetightToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbhetight"))),
    hbheisoToken(consumes<bool> (iConfig.getParameter<edm::InputTag>("hbheiso"))),
    pileupInfoToken(consumes<std::vector<PileupSummaryInfo> > ((iConfig.existsAs<edm::InputTag>("pileup") ? iConfig.getParameter<edm::InputTag>("pileup") : edm::InputTag("addPileupInfo")))),
    genevtInfoToken(consumes<GenEventInfoProduct> ((iConfig.existsAs<edm::InputTag>("genevt") ? iConfig.getParameter<edm::InputTag>("genevt") : edm::InputTag("generator")))),
    lheInfoToken(consumes<LHEEventProduct> ((iConfig.existsAs<edm::InputTag>("lheinfo") ? iConfig.getParameter<edm::InputTag>("lheinfo") : edm::InputTag("externalLHEProducer")))),
    verticesToken(consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
    gensToken(consumes<edm::View<reco::GenParticle> > ((iConfig.existsAs<edm::InputTag>("gens") ? iConfig.getParameter<edm::InputTag>("gens") : edm::InputTag("prunedGenParticles")))),
    muonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("electrons"))),
    photonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photons"))),
    tightmuonsToken(consumes<pat::MuonRefVector> (iConfig.getParameter<edm::InputTag>("tightmuons"))),
    tightelectronsToken(consumes<pat::ElectronRefVector> (iConfig.getParameter<edm::InputTag>("tightelectrons"))),
    tightphotonsToken(consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("tightphotons"))),
    electronLooseIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("electronLooseId"))),
    photonLooseIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonLooseId"))),
    photonMediumIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonMediumId"))),
    photonTightIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonTightId"))),
    photonHighPtIdToken(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("photonHighPtId"))),
    tausToken(consumes<std::vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("taus"))),
    jetsToken(consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
    pfmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("pfmet"))),
    t1pfmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("t1pfmet"))),
    mumetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("mumet"))),
    t1mumetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("t1mumet"))),
    elmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("elmet"))),
    t1elmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("t1elmet"))),
    phmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("phmet"))),
    t1phmetToken(consumes<edm::View<reco::MET> > (iConfig.getParameter<edm::InputTag>("t1phmet"))),
    applyHLTFilter(iConfig.existsAs<bool>("applyHLTFilter") ? iConfig.getParameter<bool>("applyHLTFilter") : false),
    isWorZMCSample(iConfig.existsAs<bool>("isWorZMCSample") ? iConfig.getParameter<bool>("isWorZMCSample") : false),
    isSignalSample(iConfig.existsAs<bool>("isSignalSample") ? iConfig.getParameter<bool>("isSignalSample") : false),
    cleanMuonJet(iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
    cleanElectronJet(iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
    cleanPhotonJet(iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
    uselheweights(iConfig.existsAs<bool>("uselheweights") ? iConfig.getParameter<bool>("uselheweights") : false),
    addqcdpdfweights(iConfig.existsAs<bool>("addqcdpdfweights") ? iConfig.getParameter<bool>("addqcdpdfweights") : false),
    xsec(iConfig.getParameter<double>("xsec") * 1000.0)
{

  usesResource("TFileService");
  
  wgtqcd = new double[8];
  wgtpdf = new double[100];
}


MonoJetTreeMaker::~MonoJetTreeMaker() {
    delete wgtqcd;
    delete wgtpdf;
}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get handles to all the requisite collections
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

    Handle<vector<PileupSummaryInfo> > pileupInfoH;
    iEvent.getByToken(pileupInfoToken, pileupInfoH);

    Handle<GenEventInfoProduct> genevtInfoH;
    if (uselheweights) iEvent.getByToken(genevtInfoToken, genevtInfoH);

    Handle<LHEEventProduct> lheInfoH;
    if (uselheweights) iEvent.getByToken(lheInfoToken, lheInfoH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByToken(verticesToken, verticesH);

    Handle<View<GenParticle> > gensH;
    if (isWorZMCSample || isSignalSample) iEvent.getByToken(gensToken, gensH);

    Handle<pat::MuonRefVector> muonsH;
    iEvent.getByToken(muonsToken, muonsH);
    pat::MuonRefVector muons = *muonsH;

    Handle<pat::ElectronRefVector> electronsH;
    iEvent.getByToken(electronsToken, electronsH);
    pat::ElectronRefVector electrons = *electronsH;

    Handle<pat::PhotonRefVector> photonsH;
    iEvent.getByToken(photonsToken, photonsH);
    pat::PhotonRefVector photons = *photonsH;

    Handle<pat::MuonRefVector> tightmuonsH;
    iEvent.getByToken(tightmuonsToken, tightmuonsH);
    pat::MuonRefVector tightmuons = *tightmuonsH;

    Handle<pat::ElectronRefVector> tightelectronsH;
    iEvent.getByToken(tightelectronsToken, tightelectronsH);
    pat::ElectronRefVector tightelectrons = *tightelectronsH;

    Handle<pat::PhotonRefVector> tightphotonsH;
    iEvent.getByToken(tightphotonsToken, tightphotonsH);
    pat::PhotonRefVector tightphotons = *tightphotonsH;

    Handle<edm::ValueMap<bool> > electronLooseIdH;
    iEvent.getByToken(electronLooseIdToken, electronLooseIdH);

    Handle<edm::ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdToken, photonLooseIdH);

    Handle<edm::ValueMap<bool> > photonMediumIdH;
    iEvent.getByToken(photonMediumIdToken, photonMediumIdH);

    Handle<edm::ValueMap<bool> > photonTightIdH;
    iEvent.getByToken(photonTightIdToken, photonTightIdH);

    Handle<edm::ValueMap<bool> > photonHighPtIdH;
    iEvent.getByToken(photonHighPtIdToken, photonHighPtIdH);

    Handle<std::vector<pat::Tau> > tausH;
    iEvent.getByToken(tausToken, tausH);

    Handle<std::vector<pat::Jet> > jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<View<reco::MET> > pfmetH;
    iEvent.getByToken(pfmetToken, pfmetH);

    Handle<View<reco::MET> > t1pfmetH;
    iEvent.getByToken(t1pfmetToken, t1pfmetH);

    Handle<View<reco::MET> > mumetH;
    iEvent.getByToken(mumetToken, mumetH);

    Handle<View<reco::MET> > t1mumetH;
    iEvent.getByToken(t1mumetToken, t1mumetH);

    Handle<View<reco::MET> > elmetH;
    iEvent.getByToken(elmetToken, elmetH);

    Handle<View<reco::MET> > t1elmetH;
    iEvent.getByToken(t1elmetToken, t1elmetH);

    Handle<View<reco::MET> > phmetH;
    iEvent.getByToken(phmetToken, phmetH);

    Handle<View<reco::MET> > t1phmetH;
    iEvent.getByToken(t1phmetToken, t1phmetH);

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
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
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

    // Which MET filters passed
    for (size_t i = 0; i < filterPathsVector.size(); i++) {
        if (filterPathsMap[filterPathsVector[i]] == -1) continue;
        if (i == 0  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagcsctight  = 1; // CSCTightHaloFilter
        if (i == 1  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbhenoise = 1; // HBHENoiseFilter
        if (i == 2  && filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flageebadsc   = 1; // eeBadScFilter
    }

    // Pileup info -- Will need to the updated to the Run-II specifications
    nvtx   = verticesH->size();
    puobs  = 0;
    putrue = 0;
    puwgt  = 1.;
    if (uselheweights && genevtInfoH.isValid()) wgt = genevtInfoH->weight();
    else wgt = 1.0;

    for (size_t i = 0; i < 8  ; i++) wgtqcd[i] = 0.;
    for (size_t i = 0; i < 100; i++) wgtpdf[i] = 0.;

    if (addqcdpdfweights && lheInfoH.isValid()) {
        vector<gen::WeightsInfo> weights = lheInfoH->weights();
        for (size_t i = 0; i < weights.size(); i++) {
            for (size_t j = 2; j <= 9; j++) {
                stringstream ss;
                ss << j;
                if (weights[i].id == ss.str()) wgtqcd[j-2]  = weights[i].wgt;
            }
            for (size_t j = 11; j <= 110; j++) {
                stringstream ss;
                ss << j;
                if (weights[i].id == ss.str()) wgtpdf[j-11] = weights[i].wgt;
            }
        }
    }


    if (pileupInfoH.isValid()) {
        for (vector<PileupSummaryInfo>::const_iterator pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
            if (pileupInfo_iter->getBunchCrossing() == 0) {
                puobs  = pileupInfo_iter->getPU_NumInteractions();
                putrue = pileupInfo_iter->getTrueNumInteractions();
            }
        }
    }

    // MET information 
    pfmet          = pfmetH->front().et();
    pfmetphi       = pfmetH->front().phi();

    t1pfmet        = t1pfmetH->front().et();
    t1pfmetphi     = t1pfmetH->front().phi();

    mumet          = mumetH->front().et();
    mumetphi       = mumetH->front().phi();
 
    t1mumet        = t1mumetH->front().et();
    t1mumetphi     = t1mumetH->front().phi();
 
    elmet          = elmetH->front().et();
    elmetphi       = elmetH->front().phi();
 
    t1elmet        = t1elmetH->front().et();
    t1elmetphi     = t1elmetH->front().phi();
 
    phmet          = phmetH->front().et();
    phmetphi       = phmetH->front().phi();
 
    t1phmet        = t1phmetH->front().et();
    t1phmetphi     = t1phmetH->front().phi();
 
    // Jet information
    int hardestPhotonIndex = -1;
    double hardestPhotonPt = 0.0;
    for (size_t i = 0; i < tightphotons.size(); i++) {
        if (tightphotons[i]->pt() > hardestPhotonPt) {
            hardestPhotonIndex = i;
            hardestPhotonPt = tightphotons[i]->pt();
        }
    }

    leadingjetpt  = 0.0;
    leadingjeteta = 0.0;
    leadingjetphi = 0.0;

    vector<pat::JetRef> incjets;
    vector<pat::JetRef> jets;
    for (vector<pat::Jet>::const_iterator jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
        bool skipjet = false;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        for (std::size_t j = 0; j < photons.size(); j++) {
            if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < 0.4) skipjet = true;
        }
        if (skipjet) continue;
        if (jets_iter->pt() > leadingjetpt) {
            leadingjetpt  = jets_iter->pt() ;
            leadingjeteta = jets_iter->eta();
            leadingjetphi = jets_iter->phi();
        }
        bool passjetid = false;
        if (fabs(jets_iter->eta()) <= 3.0 && jets_iter->neutralHadronEnergyFraction() < 0.99 && jets_iter->neutralEmEnergyFraction() < 0.99 && (jets_iter->chargedMultiplicity() + jets_iter->neutralMultiplicity()) > 1) {
            if (fabs(jets_iter->eta()) > 2.4) passjetid = true;
            else if (fabs(jets_iter->eta()) <= 2.4 && jets_iter->chargedHadronEnergyFraction() > 0. && jets_iter->chargedEmEnergyFraction() < 0.99 && jets_iter->chargedMultiplicity() > 0) passjetid = true;
        }
        if (fabs(jets_iter->eta()) > 3.0 && jets_iter->neutralEmEnergyFraction() < 0.9 && jets_iter->neutralMultiplicity() > 10) passjetid = true;
            
        if (!passjetid) continue;
        bool passpuid = false;
        double puidval = jets_iter->userFloat("pileupJetId:fullDiscriminant");
        double jetabseta = fabs(jets_iter->eta());
        if (jetabseta >= 0.00 && jetabseta < 2.50 && puidval > -0.63) passpuid = true;
        if (jetabseta >= 2.50 && jetabseta < 2.75 && puidval > -0.60) passpuid = true;
        if (jetabseta >= 2.75 && jetabseta < 3.00 && puidval > -0.55) passpuid = true;
        if (jetabseta >= 3.00 && jetabseta < 5.00 && puidval > -0.45) passpuid = true;
        if (!passpuid) continue;
        pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());
        incjets.push_back(jetref);
    }

    for (size_t i = 0; i < incjets.size(); i++) {
        if (fabs(incjets[i]->eta()) <= 2.5) jets.push_back(incjets[i]);
    }        

    sort(jets.begin(), jets.end(), jetsorter);
    sort(incjets.begin(), incjets.end(), jetsorter);

    // AK4 Jets
    njets       = 0;
    nbjets      = 0;
    nbjetslowpt = 0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > 30) njets++;
        if (jets[i]->pt() > 30 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjets++;
        if (jets[i]->pt() > 15 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) nbjetslowpt++;
    }

    signaljetpt        = 0.0;
    signaljeteta       = 0.0;
    signaljetphi       = 0.0;
    signaljetbtag      = 0.0;
    signaljetCHfrac    = 0.0;
    signaljetNHfrac    = 0.0;
    signaljetEMfrac    = 0.0;
    signaljetCEMfrac   = 0.0;
    signaljetmetdphi   = 0.0;
    secondjetpt        = 0.0;
    secondjeteta       = 0.0;
    secondjetphi       = 0.0;
    secondjetbtag      = 0.0;
    secondjetCHfrac    = 0.0;
    secondjetNHfrac    = 0.0;
    secondjetEMfrac    = 0.0;
    secondjetCEMfrac   = 0.0;
    secondjetmetdphi   = 0.0;
    thirdjetpt         = 0.0;
    thirdjeteta        = 0.0;
    thirdjetphi        = 0.0;
    thirdjetbtag       = 0.0;
    thirdjetCHfrac     = 0.0;
    thirdjetNHfrac     = 0.0;
    thirdjetEMfrac     = 0.0;
    thirdjetCEMfrac    = 0.0;
    thirdjetmetdphi    = 0.0;
    fourthjetpt        = 0.0;
    fourthjeteta       = 0.0;
    fourthjetphi       = 0.0;
    fourthjetbtag      = 0.0;
    fourthjetCHfrac    = 0.0;
    fourthjetNHfrac    = 0.0;
    fourthjetEMfrac    = 0.0;
    fourthjetCEMfrac   = 0.0;
    fourthjetmetdphi   = 0.0;
    jetjetdphi         = 0.0;
    jetmetdphimin      = 0.0;
    incjetmetdphimin   = 0.0;
    jetelmetdphimin    = 0.0;
    incjetelmetdphimin = 0.0;
    jetphmetdphimin    = 0.0;
    incjetphmetdphimin = 0.0;
    jetmetdphimin4     = 0.0;
    incjetmetdphimin4  = 0.0;
    jetelmetdphimin4   = 0.0;
    incjetelmetdphimin4= 0.0;
    jetphmetdphimin4   = 0.0;
    incjetphmetdphimin4= 0.0;


    if (njets > 0) {
        signaljetpt      = jets[0]->pt();
        signaljeteta     = jets[0]->eta();
        signaljetphi     = jets[0]->phi();
        signaljetbtag    = jets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        signaljetCHfrac  = jets[0]->chargedHadronEnergyFraction();
        signaljetNHfrac  = jets[0]->neutralHadronEnergyFraction();
        signaljetEMfrac  = jets[0]->neutralEmEnergyFraction();
        signaljetCEMfrac = jets[0]->chargedEmEnergyFraction();
    }

    if (njets > 1) {
        secondjetpt      = jets[1]->pt();
        secondjeteta     = jets[1]->eta();
        secondjetphi     = jets[1]->phi();
        secondjetbtag    = jets[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        secondjetCHfrac  = jets[1]->chargedHadronEnergyFraction();
        secondjetNHfrac  = jets[1]->neutralHadronEnergyFraction();
        secondjetEMfrac  = jets[1]->neutralEmEnergyFraction();
        secondjetCEMfrac = jets[1]->chargedEmEnergyFraction();
    }

    if (njets > 2) {
        thirdjetpt       = jets[2]->pt();
        thirdjeteta      = jets[2]->eta();
        thirdjetphi      = jets[2]->phi();
        thirdjetbtag     = jets[2]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        thirdjetCHfrac   = jets[2]->chargedHadronEnergyFraction();
        thirdjetNHfrac   = jets[2]->neutralHadronEnergyFraction();
        thirdjetEMfrac   = jets[2]->neutralEmEnergyFraction();
        thirdjetCEMfrac  = jets[2]->chargedEmEnergyFraction();
    }

    if (njets > 3) {
        fourthjetpt      = jets[0]->pt();
        fourthjeteta     = jets[0]->eta();
        fourthjetphi     = jets[0]->phi();
        fourthjetbtag    = jets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        fourthjetCHfrac  = jets[0]->chargedHadronEnergyFraction();
        fourthjetNHfrac  = jets[0]->neutralHadronEnergyFraction();
        fourthjetEMfrac  = jets[0]->neutralEmEnergyFraction();
        fourthjetCEMfrac = jets[0]->chargedEmEnergyFraction();
    }

    if (signaljetpt > 0.0 && secondjetpt > 0.0) jetjetdphi = deltaPhi(signaljetphi, secondjetphi);
    if (signaljetpt > 0.0) signaljetmetdphi = deltaPhi(signaljetphi, t1pfmetphi);
    if (secondjetpt > 0.0) secondjetmetdphi = deltaPhi(secondjetphi, t1pfmetphi);
    if (thirdjetpt  > 0.0) thirdjetmetdphi  = deltaPhi(thirdjetphi , t1pfmetphi);
    if (fourthjetpt > 0.0) fourthjetmetdphi = deltaPhi(fourthjetphi, t1pfmetphi);

    std::vector<double> jetmetdphiminvector;
    std::vector<double> jetmetdphimin4vector;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i]->pt() > 30) {
            double jetphi = jets[i]->phi();
            jetmetdphiminvector.push_back(fabs(deltaPhi(jetphi, t1mumetphi)));
            if (i < 4) jetmetdphimin4vector.push_back(fabs(deltaPhi(jetphi, t1mumetphi)));
        }
    }
    if (jetmetdphiminvector .size() > 0) jetmetdphimin  = *min_element(jetmetdphiminvector .begin(), jetmetdphiminvector .end());
    if (jetmetdphimin4vector.size() > 0) jetmetdphimin4 = *min_element(jetmetdphimin4vector.begin(), jetmetdphimin4vector.end());

    std::vector<double> incjetmetdphiminvector;
    std::vector<double> incjetmetdphimin4vector;
    for (size_t i = 0; i < incjets.size(); i++) {
        if (incjets[i]->pt() > 30) {
            //double incjetphi = incjets[i]->phi();
            double incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
            incjetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
            if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetphi)));
        }
    }
    if (incjetmetdphiminvector .size() > 0) incjetmetdphimin  = *min_element(incjetmetdphiminvector .begin(), incjetmetdphiminvector .end());
    if (incjetmetdphimin4vector.size() > 0) incjetmetdphimin4 = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());

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

    // Lepton counts
    nmuons          = muonsH->size();
    nelectrons      = electronsH->size();
    ntightmuons     = tightmuonsH->size();
    ntightelectrons = tightelectronsH->size();
    ntaus           = 0;

    for (vector<pat::Tau>::const_iterator taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
        bool skiptau = false;
        for (std::size_t j = 0; j < muons.size(); j++) {
            if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
        }
        for (std::size_t j = 0; j < electrons.size(); j++) {
            if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < 0.4) skiptau = true;
        }
        if (taus_iter->pt() > 18 && fabs(taus_iter->eta()) < 2.3 && taus_iter->tauID("decayModeFinding") > 0.5 && taus_iter->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 5 && !skiptau) ntaus++;
    }

    vector<pat::MuonRef> muonvector;
    for (size_t i = 0; i < muons.size(); i++) muonvector.push_back(muons[i]);

    vector<pat::ElectronRef> electronvector;
    for (size_t i = 0; i < electrons.size(); i++) electronvector.push_back(electrons[i]);

    // W, Z control sample information
    zmass       = 0.0; 
    zpt         = 0.0;
    zeta        = 0.0; 
    zphi        = 0.0;
    zeemass     = 0.0; 
    zeept       = 0.0;
    zeeeta      = 0.0; 
    zeephi      = 0.0;
    wmt         = 0.0;
    wemt        = 0.0;
    emumass     = 0.0;
    emupt       = 0.0;
    emueta      = 0.0;
    emuphi      = 0.0;
    mu1pid      = 0;
    mu1pt       = 0.0;
    mu1eta      = 0.0;
    mu1phi      = 0.0;
    mu1pfpt     = 0.0;
    mu1pfeta    = 0.0;
    mu1pfphi    = 0.0;
    mu1id       = 0;
    mu1idm      = 0;
    mu1idt      = 0;
    mu1iso      = 0.0;
    mu2pid      = 0;
    mu2pt       = 0.0;
    mu2eta      = 0.0; 
    mu2phi      = 0.0;
    mu2pfpt     = 0.0;
    mu2pfeta    = 0.0;
    mu2pfphi    = 0.0;
    mu2id       = 0;
    mu2idm      = 0;
    mu2idt      = 0;
    mu2iso      = 0.0;
    el1pid      = 0;
    el1pt       = 0.0;
    el1eta      = 0.0;
    el1phi      = 0.0;
    el1id       = 0;
    el2pid      = 0;
    el2pt       = 0.0;
    el2eta      = 0.0;
    el2phi      = 0.0;
    el2id       = 0;

    sort(muonvector.begin(), muonvector.end(), muonsorter);
    sort(electronvector.begin(), electronvector.end(), electronsorter);

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
        if (verticesH->size() > 0) mu1idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);

        for (std::size_t i = 0; i < tightmuons.size(); i++) {
            if (muon == tightmuons[i]) mu1id = 1;
        }

        if (nmuons == 1) wmt = sqrt(2.0 * mu1pt * t1pfmet * (1.0 - cos(deltaPhi(mu1phi, t1pfmetphi))));
    }
   

 
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
        if (verticesH->size() > 0) mu2idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);
    
        for (std::size_t i = 0; i < tightmuons.size(); i++) {
            if (muon == tightmuons[i]) mu2id = 1;
        }
    
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
        TLorentzVector mu2vec; mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muon->p());
    
        TLorentzVector zvec(mu1vec);
        zvec += mu2vec;
    
        zmass = zvec.M();
        zpt   = zvec.Pt();
        zeta  = zvec.Eta();            
        zphi  = zvec.Phi();
    }

    if (nelectrons == 1 || nelectrons == 2) {
        pat::ElectronRef electron = electrons[0];
        el1pid = electron->pdgId();
        el1pt  = electron->pt();
        el1eta = electron->eta();
        el1phi = electron->phi();
        el1idl = ((*electronLooseIdH )[electron] ? 1 : 0);
        
        for (std::size_t i = 0; i < tightelectrons.size(); i++) {
            if (electron == tightelectrons[i]) el1id = 1;
        }
        
        if (electrons.size() == 1) wemt = sqrt(2.0 * el1pt * t1pfmet * (1.0 - cos(deltaPhi(el1phi, t1pfmetphi))));
    }

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

        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
        TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electron->p());

        TLorentzVector zvec(el1vec);
        zvec += el2vec;

        zeemass = zvec.M();
        zeept   = zvec.Pt();
        zeeeta  = zvec.Eta();
        zeephi  = zvec.Phi();
    }

    if (nmuons == 1 && nelectrons == 1) {
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[0]->p());
        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electrons[0]->p());
        
        TLorentzVector emuvec(mu1vec);
        emuvec += el1vec;
        
        emumass = emuvec.M();
        emupt   = emuvec.Pt();
        emueta  = emuvec.Eta();
        emuphi  = emuvec.Phi();
    } 

    // Photon information
    phidl    = 0;
    phidm    = 0;
    phidt    = 0;
    phidh    = 0;
    phpt     = 0.0;
    pheta    = 0.0;
    phphi    = 0.0;
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

    // Generator-level information
    wzid          = 0;
    wzmass        = 0.0;
    wzpt          = 0.0;
    wzeta         = 0.0;
    wzphi         = 0.0;
    l1id          = 0;
    l1pt          = 0.0;
    l1eta         = 0.0;
    l1phi         = 0.0;
    l2id          = 0;
    l2pt          = 0.0;
    l2eta         = 0.0;
    l2phi         = 0.0;
    parid         = 0;
    parpt         = 0.0;
    pareta        = 0.0;
    parphi        = 0.0;
    ancid         = 0;
    ancpt         = 0.0;
    anceta        = 0.0;
    ancphi        = 0.0;

    if (isWorZMCSample && gensH.isValid()) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() > 1 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) {
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

        if (wzid == 0) {
            for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
                if (gens_iter->pdgId() == 22 && gens_iter->status() == 1 && gens_iter->isPromptFinalState() && gens_iter->pt() > wzpt) {
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

    if (isSignalSample && gensH.isValid()) {
        TLorentzVector dm1vec; 
        TLorentzVector dm2vec; 
        bool foundfirst = false;
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if (gens_iter->pdgId() == 1000022 && !foundfirst) {
                dm1vec.SetPtEtaPhiE(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->p());
                foundfirst = true;
            }
            if (gens_iter->pdgId() == 1000022 &&  foundfirst) {
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
    tree = fs->make<TTree>("tree"       , "tree");
    // Run, Lumi, Event info
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
    // Event weights
    tree->Branch("xsec"                 , &xsec                 , "xsec/D");
    tree->Branch("wgt"                  , &wgt                  , "wgt/D");
    if (addqcdpdfweights) {
    tree->Branch("wgtpdf"               ,  wgtpdf               , "wgtpdf[100]/D");
    tree->Branch("wgtqcd"               ,  wgtqcd               , "wgtqcd[8]/D");
    }
    tree->Branch("puwgt"                , &puwgt                , "puwgt/D");
    // Pileup info
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
    tree->Branch("ntightelectrons"      , &ntightelectrons      , "ntightelectrons/i");
    tree->Branch("ntaus"                , &ntaus                , "ntaus/i");
    tree->Branch("njets"                , &njets                , "njets/i");
    tree->Branch("nbjets"               , &nbjets               , "nbjets/i");
    tree->Branch("nbjetslowpt"          , &nbjetslowpt          , "nbjetslowpt/i");
    tree->Branch("nphotons"             , &nphotons             , "nphotons/i");
    // MET info
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
    tree->Branch("t1phmet"              , &t1phmet              , "t1phmet/D");
    tree->Branch("t1phmetphi"           , &t1phmetphi           , "t1phmetphi/D");
    // Jet info
    tree->Branch("leadingjetpt"         , &leadingjetpt         , "leadingjetpt/D");
    tree->Branch("leadingjeteta"        , &leadingjeteta        , "leadingjeteta/D");
    tree->Branch("leadingjetphi"        , &leadingjetphi        , "leadingjetphi/D");
    tree->Branch("signaljetpt"          , &signaljetpt          , "signaljetpt/D");
    tree->Branch("signaljeteta"         , &signaljeteta         , "signaljeteta/D");
    tree->Branch("signaljetphi"         , &signaljetphi         , "signaljetphi/D");
    tree->Branch("signaljetbtag"        , &signaljetbtag        , "signaljetbtag/D");
    tree->Branch("signaljetCHfrac"      , &signaljetCHfrac      , "signaljetCHfrac/D");
    tree->Branch("signaljetNHfrac"      , &signaljetNHfrac      , "signaljetNHfrac/D");
    tree->Branch("signaljetEMfrac"      , &signaljetEMfrac      , "signaljetEMfrac/D");
    tree->Branch("signaljetCEMfrac"     , &signaljetCEMfrac     , "signaljetCEMfrac/D");
    tree->Branch("signaljetmetdphi"     , &signaljetmetdphi     , "signaljetmetdphi/D");
    tree->Branch("secondjetpt"          , &secondjetpt          , "secondjetpt/D");
    tree->Branch("secondjeteta"         , &secondjeteta         , "secondjeteta/D");
    tree->Branch("secondjetphi"         , &secondjetphi         , "secondjetphi/D");
    tree->Branch("secondjetbtag"        , &secondjetbtag        , "secondjetbtag/D");
    tree->Branch("secondjetCHfrac"      , &secondjetCHfrac      , "secondjetCHfrac/D");
    tree->Branch("secondjetNHfrac"      , &secondjetNHfrac      , "secondjetNHfrac/D");
    tree->Branch("secondjetEMfrac"      , &secondjetEMfrac      , "secondjetEMfrac/D");
    tree->Branch("secondjetCEMfrac"     , &secondjetCEMfrac     , "secondjetCEMfrac/D");
    tree->Branch("secondjetmetdphi"     , &secondjetmetdphi     , "secondjetmetdphi/D");
    tree->Branch("thirdjetpt"           , &thirdjetpt           , "thirdjetpt/D");
    tree->Branch("thirdjeteta"          , &thirdjeteta          , "thirdjeteta/D");
    tree->Branch("thirdjetphi"          , &thirdjetphi          , "thirdjetphi/D");
    tree->Branch("thirdjetbtag"         , &thirdjetbtag         , "thirdjetbtag/D");
    tree->Branch("thirdjetCHfrac"       , &thirdjetCHfrac       , "thirdjetCHfrac/D");
    tree->Branch("thirdjetNHfrac"       , &thirdjetNHfrac       , "thirdjetNHfrac/D");
    tree->Branch("thirdjetEMfrac"       , &thirdjetEMfrac       , "thirdjetEMfrac/D");
    tree->Branch("thirdjetCEMfrac"      , &thirdjetCEMfrac      , "thirdjetCEMfrac/D");
    tree->Branch("thirdjetmetdphi"      , &thirdjetmetdphi      , "thirdjetmetdphi/D");
    tree->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/D");
    tree->Branch("jetmetdphimin"        , &jetmetdphimin        , "jetmetdphimin/D");
    tree->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/D");
    tree->Branch("jetelmetdphimin"      , &jetelmetdphimin      , "jetelmetdphimin/D");
    tree->Branch("incjetelmetdphimin"   , &incjetelmetdphimin   , "incjetelmetdphimin/D");
    tree->Branch("jetphmetdphimin"      , &jetphmetdphimin      , "jetphmetdphimin/D");
    tree->Branch("incjetphmetdphimin"   , &incjetphmetdphimin   , "incjetphmetdphimin/D");
    tree->Branch("jetmetdphimin4"       , &jetmetdphimin4       , "jetmetdphimin4/D");
    tree->Branch("incjetmetdphimin4"    , &incjetmetdphimin4    , "incjetmetdphimin4/D");
    tree->Branch("jetelmetdphimin4"     , &jetelmetdphimin4     , "jetelmetdphimin4/D");
    tree->Branch("incjetelmetdphimin4"  , &incjetelmetdphimin4  , "incjetelmetdphimin4/D");
    tree->Branch("jetphmetdphimin4"     , &jetphmetdphimin4     , "jetphmetdphimin4/D");
    tree->Branch("incjetphmetdphimin4"  , &incjetphmetdphimin4  , "incjetphmetdphimin4/D");
    tree->Branch("ht"                   , &ht                   , "ht/D");
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
}

void MonoJetTreeMaker::endJob() {
}

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

void MonoJetTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

/*
This code is ripped off from https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/PhotonNtuplerMiniAOD.cc
*/
void MonoJetTreeMaker::findFirstNonPhotonMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {
    if (particle == 0) {
        return;
    }
    if (abs(particle->pdgId()) == 22) {
        findFirstNonPhotonMother(particle->mother(0), ancestorid, ancestorpt, ancestoreta, ancestorphi);
    }
    else {
        ancestorid  = particle->pdgId();
        ancestorpt  = particle->pt();
        ancestoreta = particle->eta();
        ancestorphi = particle->phi();
    }
    return;
}

void MonoJetTreeMaker::findMother(const reco::Candidate *particle, int& ancestorid, double& ancestorpt, double& ancestoreta, double& ancestorphi) {
    if (particle == 0) {
        return;
    }
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

void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MonoJetTreeMaker);

