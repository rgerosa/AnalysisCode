#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include <TH1F.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

class MonoJetTreeMaker : public edm::EDAnalyzer {
    public:
        explicit MonoJetTreeMaker(const edm::ParameterSet&);
        ~MonoJetTreeMaker();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

        void initPileupWeights();            

        edm::InputTag pileupInfoTag;
        edm::InputTag verticesTag;
        edm::InputTag gensTag;
        edm::InputTag muonsTag;
        edm::InputTag electronsTag;
        edm::InputTag tausTag;
        edm::InputTag jetsTag;
        edm::InputTag c0metTag;
        edm::InputTag triggerResultsTag;
        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        std::string jetCorrectionService;
        bool isControlSample;            
        TTree* tree;
        TH1F* puhist;

        int32_t puobs, putrue; 
        int32_t wzid, l1id, l2id, mu1id, mu2id; 
        uint32_t event, run, lumi;
        uint32_t nvtx, ngoodmuons, nmuons, nelectrons, ntaus, njets, nsignaljets;
        uint32_t hltmet120, hltmet95jet80, hltmet105jet80, hltdoublemu, hltsinglemu;
        double pfmet, pfmetphi, c0met, c0metphi, c1met, c1metphi, c2met, c2metphi;
        double signaljetpt, signaljeteta, signaljetphi, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac;
        double secondjetpt, secondjeteta, secondjetphi, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, jetjetdphi;
        double wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
        double zmass, zpt, zeta, zphi, wmt, mu1pt, mu1eta, mu1phi, mu2pt, mu2eta, mu2phi;
        double wgt, puwgt;
};

MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
    pileupInfoTag(iConfig.getParameter<edm::InputTag>("pileup")),
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
    gensTag((iConfig.existsAs<edm::InputTag>("gens") ? iConfig.getParameter<edm::InputTag>("gens") : edm::InputTag("genParticles"))),
    muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
    electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
    tausTag(iConfig.getParameter<edm::InputTag>("taus")),
    jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
    c0metTag(iConfig.getParameter<edm::InputTag>("met")),
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    jetCorrectionService(iConfig.getParameter<std::string>("jec")),
    isControlSample(iConfig.existsAs<bool>("isControlSample") ? iConfig.getParameter<bool>("isControlSample") : false),
    wgt(iConfig.existsAs<double>("weight") ? iConfig.getParameter<double>("weight") : 1.0)
{
    initPileupWeights();
}


MonoJetTreeMaker::~MonoJetTreeMaker() {
    if (puhist) delete puhist;
}

void MonoJetTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;


    // Get handles to all the requisite collections
    Handle<TriggerResults> triggerResultsH;
    iEvent.getByLabel(triggerResultsTag, triggerResultsH);

    Handle<vector<PileupSummaryInfo> > pileupInfoH;
    iEvent.getByLabel(pileupInfoTag, pileupInfoH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByLabel(verticesTag, verticesH);

    Handle<View<GenParticle> > gensH;
    if (isControlSample) iEvent.getByLabel(gensTag, gensH);

    Handle<pat::MuonCollection> muonsH;
    iEvent.getByLabel(muonsTag, muonsH);
    pat::MuonCollection muons = *muonsH;

    Handle<pat::ElectronCollection> electronsH;
    iEvent.getByLabel(electronsTag, electronsH);
    pat::ElectronCollection electrons = *electronsH;

    Handle<pat::TauCollection> tausH;
    iEvent.getByLabel(tausTag, tausH);

    Handle<pat::JetCollection> jetsH;
    iEvent.getByLabel(jetsTag, jetsH);
    pat::JetCollection inputjets = *jetsH;

    Handle<View<MET> > c0metH;
    iEvent.getByLabel(c0metTag, c0metH);

    Handle<View<MET> > pfMetH;
    iEvent.getByLabel("pfMet", pfMetH);

    // Event, lumi, run info
    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();

    // Trigger info
    hltdoublemu    = 0;
    hltsinglemu    = 0;
    hltmet120      = 0;
    hltmet95jet80  = 0;
    hltmet105jet80 = 0;
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
        if (i == 0 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu    = 1; // Double muon triggers
        if (i == 1 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu    = 1; // Double muon triggers
        if (i == 2 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu    = 1; // Single muon triggers
        if (i == 3 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120      = 1; // MET trigger trigger
        if (i == 4 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet95jet80  = 1; // Jet-MET trigger trigger
        if (i == 5 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet105jet80 = 1; // Jet-MET trigger trigger
    }

    // Pileup info
    nvtx = verticesH->size();
    if (pileupInfoH.isValid()) {
        for (vector<PileupSummaryInfo>::const_iterator pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
            if (pileupInfo_iter->getBunchCrossing() == 0) {
                puobs  = pileupInfo_iter->getPU_NumInteractions();
                putrue = pileupInfo_iter->getTrueNumInteractions();
            }
        }
    }

    // MET information 
    pfmet    = pfMetH->front().et();
    pfmetphi = pfMetH->front().phi();
 
    c0met    = c0metH->front().et();
    c0metphi = c0metH->front().phi();
  
    double c1metx = pfmet * cos(pfmetphi);
    double c1mety = pfmet * sin(pfmetphi);
 
    double c2metx = pfmet * cos(pfmetphi);
    double c2mety = pfmet * sin(pfmetphi);
 
    for (size_t i = 0; i < muons.size(); i++) {
        if (muons[i].isPFMuon()) {
            c1metx += muons[i].pfP4().Pt() * cos(muons[i].pfP4().Phi());
            c1mety += muons[i].pfP4().Pt() * sin(muons[i].pfP4().Phi());
        }
        if (muons[i].pt() > 10 && muons[i].isPFMuon() && (muons[i].isGlobalMuon() || muons[i].isTrackerMuon())) {
            c2metx += muons[i].pfP4().Pt() * cos(muons[i].pfP4().Phi());
            c2mety += muons[i].pfP4().Pt() * sin(muons[i].pfP4().Phi());
        }
    }
   
    c1met    = sqrt(c1metx*c1metx + c1mety*c1mety);  
    c1metphi = atan2(c1mety, c1metx);
  
    c2met    = sqrt(c2metx*c2metx + c2mety*c2mety);  
    c2metphi = atan2(c2mety, c2metx);
 
    // Jet selection
    const JetCorrector* corrector = JetCorrector::getJetCorrector(jetCorrectionService, iSetup); 
    pat::JetCollection jets;
    for (size_t i = 0; i < inputjets.size(); i++) {
        if (inputjets[i].electronEnergyFraction() > 0.5 || inputjets[i].muonEnergyFraction() > 0.5) continue;
        pat::Jet jet = inputjets[i].correctedJet("Uncorrected");
        double jec = corrector->correction(jet, iEvent, iSetup);
        jet.scaleEnergy(jec);
        if (jet.pt() < 30 || fabs(jet.eta()) > 4.5) continue;
        jets.push_back(jet);
    }

    // Jet information
    njets = jets.size();

    int hardestJetIndex = -1;
    float hardestJetPt = 0.0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > hardestJetPt) {
            hardestJetIndex = i;
            hardestJetPt = jets[i].pt();
        }                        
    }

    int secondJetIndex = -1;
    float secondJetPt = 0.0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (hardestJetIndex >= 0 && i != unsigned(hardestJetIndex) && jets[i].pt() > secondJetPt) {
            secondJetIndex = i;
            secondJetPt = jets[i].pt();
        }                        
    }

    nsignaljets = 0;
    if (hardestJetIndex >= 0) {
        if (jets[hardestJetIndex].chargedHadronEnergyFraction() > 0.2 && jets[hardestJetIndex].neutralHadronEnergyFraction() < 0.7 && jets[hardestJetIndex].neutralEmEnergyFraction() < 0.7) nsignaljets = 1;
        signaljetpt     = jets[hardestJetIndex].pt();
        signaljeteta    = jets[hardestJetIndex].eta();
        signaljetphi    = jets[hardestJetIndex].phi();
        signaljetCHfrac = jets[hardestJetIndex].chargedHadronEnergyFraction();
        signaljetNHfrac = jets[hardestJetIndex].neutralHadronEnergyFraction();
        signaljetEMfrac = jets[hardestJetIndex].neutralEmEnergyFraction();
    }
    else {
        signaljetpt     = 0.0;
        signaljeteta    = 0.0;
        signaljetphi    = 0.0;
        signaljetCHfrac = 0.0;
        signaljetNHfrac = 0.0;
        signaljetEMfrac = 0.0;
    }

    if (secondJetIndex >= 0) {
        secondjetpt     = jets[secondJetIndex].pt();
        secondjeteta    = jets[secondJetIndex].eta();
        secondjetphi    = jets[secondJetIndex].phi();
        secondjetCHfrac = jets[secondJetIndex].chargedHadronEnergyFraction();
        secondjetNHfrac = jets[secondJetIndex].neutralHadronEnergyFraction();
        secondjetEMfrac = jets[secondJetIndex].neutralEmEnergyFraction();
    }
    else {
        secondjetpt     = 0.0;
        secondjeteta    = 0.0;
        secondjetphi    = 0.0;
        secondjetCHfrac = 0.0;
        secondjetNHfrac = 0.0;
        secondjetEMfrac = 0.0;
    }

    jetjetdphi = 0.0;
    if (signaljetpt > 0.0 && secondjetpt > 0.0) jetjetdphi = deltaPhi(signaljetphi, secondjetphi);

    // Lepton counts
    nmuons      = 0;
    nelectrons  = 0;
    ntaus       = tausH->size();

    for (size_t i = 0; i < muons.size(); i++) {
        if (muons[i].pt() > 10 && fabs(muons[i].eta()) < 2.4 && muons[i].isPFMuon() && (muons[i].isGlobalMuon() || muons[i].isTrackerMuon())) nmuons++;
    }

    for (size_t i = 0; i < electrons.size(); i++) {
        if (electrons[i].pt() > 10 && 
            (fabs(electrons[i].superCluster()->eta()) < 1.4442 || fabs(electrons[i].superCluster()->eta()) > 1.566) && fabs(electrons[i].superCluster()->eta()) < 2.5 && 
            (electrons[i].electronID("eidVBTFCom95") == 5 || electrons[i].electronID("eidVBTFCom95") == 7) &&
            (electrons[i].pfIsolationVariables().chargedHadronIso + electrons[i].pfIsolationVariables().neutralHadronIso + electrons[i].pfIsolationVariables().photonIso)/electrons[i].pt() < 0.2) nelectrons++;
    }

    // Generator-level information
    wzid   = 0;
    wzmass = 0.0;
    wzpt   = 0.0;
    wzeta  = 0.0;
    wzphi  = 0.0;
    l1id   = 0;
    l1pt   = 0.0;
    l1eta  = 0.0;
    l1phi  = 0.0;
    l2id   = 0;
    l2pt   = 0.0;
    l2eta  = 0.0;
    l2phi  = 0.0;

    if (isControlSample && gensH.isValid()) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() == 3) {
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
                break;
            }
        }            
    }

    // W, Z control sample information
    ngoodmuons  = 0;
    zmass       = 0.0; 
    zpt         = 0.0;
    zeta        = 0.0; 
    zphi        = 0.0;
    wmt         = 0.0;
    mu1id       = 0;
    mu1pt       = 0.0;
    mu1eta      = 0.0;
    mu1phi      = 0.0;
    mu2id       = 0;
    mu2pt       = 0.0;
    mu2eta      = 0.0; 
    mu2phi      = 0.0;
    if (isControlSample) {
        std::vector<unsigned> goodMuonIndices;
        for (size_t i = 0; i < muons.size(); i++) {
            if (muons[i].pt() > 20 && muons[i].eta() < 2.4 && muons[i].isGlobalMuon() && muons[i].isTrackerMuon() && muons[i].isPFMuon() && 
                muons[i].globalTrack()->normalizedChi2() < 10 && 
                muons[i].globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
                muons[i].innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && 
                muons[i].innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && 
                muons[i].numberOfMatchedStations() > 1 && 
                muons[i].dB(pat::Muon::PV2D) < 0.02 && 
                (muons[i].pfIsolationR04().sumChargedHadronPt + muons[i].pfIsolationR04().sumNeutralHadronEt + muons[i].pfIsolationR04().sumPhotonEt)/muons[i].pt() < 0.2) goodMuonIndices.push_back(i);
        } 
        ngoodmuons = goodMuonIndices.size();

        if (ngoodmuons == 2 && muons[goodMuonIndices[0]].pdgId() == -muons[goodMuonIndices[1]].pdgId()) {
            mu1id  = muons[goodMuonIndices[0]].pdgId(); 
            mu1pt  = muons[goodMuonIndices[0]].pt(); 
            mu1eta = muons[goodMuonIndices[0]].eta(); 
            mu1phi = muons[goodMuonIndices[0]].phi();
 
            mu2id  = muons[goodMuonIndices[1]].pdgId(); 
            mu2pt  = muons[goodMuonIndices[1]].pt(); 
            mu2eta = muons[goodMuonIndices[1]].eta(); 
            mu2phi = muons[goodMuonIndices[1]].phi();

            TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muons[goodMuonIndices[0]].p());
            TLorentzVector mu2vec; mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muons[goodMuonIndices[1]].p());

            TLorentzVector zvec(mu1vec);
            zvec += mu2vec;
        
            zmass = zvec.M();
            zpt   = zvec.Pt();
            zeta  = zvec.Eta();            
            zphi  = zvec.Phi();
        }
        if (ngoodmuons == 1) {
            mu1id  = muons[goodMuonIndices[0]].pdgId(); 
            mu1pt  = muons[goodMuonIndices[0]].pt(); 
            mu1eta = muons[goodMuonIndices[0]].eta(); 
            mu1phi = muons[goodMuonIndices[0]].phi();
            wmt    = sqrt(2.0 * mu1pt * pfmet * (1.0 - cos(deltaPhi(mu1phi, pfmetphi))));
        } 
    }

    // Adding pileup reweighting information
    puwgt = puhist->GetBinContent(puhist->FindBin(putrue));

    tree->Fill();
}


void MonoJetTreeMaker::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"  , "tree");
    tree->Branch("puobs"           , &puobs           , "puobs/I");
    tree->Branch("putrue"          , &putrue          , "putrue/I");
    tree->Branch("event"           , &event           , "event/i");
    tree->Branch("run"             , &run             , "run/i");
    tree->Branch("lumi"            , &lumi            , "lumi/i");
    tree->Branch("nvtx"            , &nvtx            , "nvtx/i");
    tree->Branch("nmuons"          , &nmuons          , "nmuons/i");
    tree->Branch("nelectrons"      , &nelectrons      , "nelectrons/i");
    tree->Branch("ntaus"           , &ntaus           , "ntaus/i");
    tree->Branch("njets"           , &njets           , "njets/i");
    tree->Branch("nsignaljets"     , &nsignaljets     , "nsignaljets/i");
    tree->Branch("pfmet"           , &pfmet           , "pfmet/D");
    tree->Branch("pfmetphi"        , &pfmetphi        , "pfmetphi/D");
    tree->Branch("c0met"           , &c0met           , "c0met/D");
    tree->Branch("c0metphi"        , &c0metphi        , "c0metphi/D");
    tree->Branch("c1met"           , &c1met           , "c1met/D");
    tree->Branch("c1metphi"        , &c1metphi        , "c1metphi/D");
    tree->Branch("c2met"           , &c2met           , "c2met/D");
    tree->Branch("c2metphi"        , &c2metphi        , "c2metphi/D");
    tree->Branch("signaljetpt"     , &signaljetpt     , "signaljetpt/D");
    tree->Branch("signaljeteta"    , &signaljeteta    , "signaljeteta/D");
    tree->Branch("signaljetphi"    , &signaljetphi    , "signaljetphi/D");
    tree->Branch("signaljetCHfrac" , &signaljetCHfrac , "signaljetCHfrac/D");
    tree->Branch("signaljetNHfrac" , &signaljetNHfrac , "signaljetNHfrac/D");
    tree->Branch("signaljetEMfrac" , &signaljetEMfrac , "signaljetEMfrac/D");
    tree->Branch("secondjetpt"     , &secondjetpt     , "secondjetpt/D");
    tree->Branch("secondjeteta"    , &secondjeteta    , "secondjeteta/D");
    tree->Branch("secondjetphi"    , &secondjetphi    , "secondjetphi/D");
    tree->Branch("secondjetCHfrac" , &secondjetCHfrac , "secondjetCHfrac/D");
    tree->Branch("secondjetNHfrac" , &secondjetNHfrac , "secondjetNHfrac/D");
    tree->Branch("secondjetEMfrac" , &secondjetEMfrac , "secondjetEMfrac/D");
    tree->Branch("jetjetdphi"      , &jetjetdphi      , "jetjetdphi/D");
    tree->Branch("wgt"             , &wgt             , "wgt/D");
    tree->Branch("puwgt"           , &puwgt           , "puwgt/D");
    tree->Branch("hltmet120"       , &hltmet120       , "hltmet120/i");
    tree->Branch("hltmet95jet80"   , &hltmet95jet80   , "hltmet95jet80/i");
    tree->Branch("hltmet105jet80"  , &hltmet105jet80  , "hltmet105jet80/i");
    
    if (isControlSample) {
    tree->Branch("wzid"            , &wzid            , "wzid/I");
    tree->Branch("wzmass"          , &wzmass          , "wzmass/D");
    tree->Branch("wzmt"            , &wzmt            , "wzmt/D");
    tree->Branch("wzpt"            , &wzpt            , "wzpt/D");
    tree->Branch("wzeta"           , &wzeta           , "wzeta/D");
    tree->Branch("wzphi"           , &wzphi           , "wzphi/D");
    tree->Branch("l1id"            , &l1id            , "l1id/I");
    tree->Branch("l1pt"            , &l1pt            , "l1pt/D");
    tree->Branch("l1eta"           , &l1eta           , "l1eta/D");
    tree->Branch("l1phi"           , &l1phi           , "l1phi/D");
    tree->Branch("l2id"            , &l2id            , "l2id/I");
    tree->Branch("l2pt"            , &l2pt            , "l2pt/D");
    tree->Branch("l2eta"           , &l2eta           , "l2eta/D");
    tree->Branch("l2phi"           , &l2phi           , "l2phi/D");
    tree->Branch("ngoodmuons"      , &ngoodmuons      , "ngoodmuons/i");
    tree->Branch("zmass"           , &zmass           , "zmass/D");
    tree->Branch("zpt"             , &zpt             , "zpt/D");
    tree->Branch("zeta"            , &zeta            , "zeta/D");
    tree->Branch("zphi"            , &zphi            , "zphi/D");
    tree->Branch("wmt"             , &wmt             , "wmt/D");
    tree->Branch("mu1id"           , &mu1id           , "mu1id/I");
    tree->Branch("mu1pt"           , &mu1pt           , "mu1pt/D");
    tree->Branch("mu1eta"          , &mu1eta          , "mu1eta/D");
    tree->Branch("mu1phi"          , &mu1phi          , "mu1phi/D");
    tree->Branch("mu2id"           , &mu2id           , "mu2id/I");
    tree->Branch("mu2pt"           , &mu2pt           , "mu2pt/D");
    tree->Branch("mu2eta"          , &mu2eta          , "mu2eta/D");
    tree->Branch("mu2phi"          , &mu2phi          , "mu2phi/D");
    tree->Branch("hltdoublemu"     , &hltdoublemu     , "hltdoublemu/i");
    tree->Branch("hltsinglemu"     , &hltsinglemu     , "hltsinglemu/i");
    }
}

void MonoJetTreeMaker::endJob() {
}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    triggerPathsVector.push_back("HLT_Mu17_Mu8");
    triggerPathsVector.push_back("HLT_Mu17_TkMu8");
    triggerPathsVector.push_back("HLT_IsoMu24_eta2p1");
    triggerPathsVector.push_back("HLT_MET120_HBHENoiseCleaned");
    triggerPathsVector.push_back("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95");
    triggerPathsVector.push_back("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95");

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
}

void MonoJetTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void MonoJetTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MonoJetTreeMaker);

void MonoJetTreeMaker::initPileupWeights() {
    float puweights[60];
    puweights[0]  = 0.242421;
    puweights[1]  = 0.314567;
    puweights[2]  = 0.328551;
    puweights[3]  = 0.341325;
    puweights[4]  = 0.311977;
    puweights[5]  = 0.560934;
    puweights[6]  = 0.443155;
    puweights[7]  = 0.444942;
    puweights[8]  = 0.625725;
    puweights[9]  = 0.940482;
    puweights[10] = 1.33114;
    puweights[11] = 1.67986;
    puweights[12] = 1.7303;
    puweights[13] = 1.54809;
    puweights[14] = 1.32353;
    puweights[15] = 1.15668;
    puweights[16] = 1.07105;
    puweights[17] = 1.04548;
    puweights[18] = 1.06319;
    puweights[19] = 1.10733;
    puweights[20] = 1.14954;
    puweights[21] = 1.17563;
    puweights[22] = 1.188;
    puweights[23] = 1.18718;
    puweights[24] = 1.16671;
    puweights[25] = 1.121;
    puweights[26] = 1.04977;
    puweights[27] = 0.956974;
    puweights[28] = 0.84729;
    puweights[29] = 0.727003;
    puweights[30] = 0.603974;
    puweights[31] = 0.485796;
    puweights[32] = 0.377733;
    puweights[33] = 0.283343;
    puweights[34] = 0.204364;
    puweights[35] = 0.14118;
    puweights[36] = 0.0934506;
    puweights[37] = 0.059445;
    puweights[38] = 0.0365081;
    puweights[39] = 0.0218306;
    puweights[40] = 0.012844;
    puweights[41] = 0.00753269;
    puweights[42] = 0.00447223;
    puweights[43] = 0.00273386;
    puweights[44] = 0.00175157;
    puweights[45] = 0.00118879;
    puweights[46] = 0.000857334;
    puweights[47] = 0.000653996;
    puweights[48] = 0.000522478;
    puweights[49] = 0.000432433;
    puweights[50] = 0.000367567;
    puweights[51] = 0.000318451;
    puweights[52] = 0.000279865;
    puweights[53] = 0.000248423;
    puweights[54] = 0.000221711;
    puweights[55] = 0.000198398;
    puweights[56] = 0.000177509;
    puweights[57] = 0.000158456;
    puweights[58] = 0.000140801;
    puweights[59] = 0.000261544;

    puhist = new TH1F("puhist", "", 60, 0., 60.);

    for(int k = 0; k < 60; k++) puhist->SetBinContent(k+1, puweights[k]);
}
