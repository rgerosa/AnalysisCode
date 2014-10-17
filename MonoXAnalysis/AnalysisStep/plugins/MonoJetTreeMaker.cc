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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

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
        edm::InputTag pfmetTag;
        edm::InputTag t1pfmetTag;
        edm::InputTag calometTag;
        edm::InputTag pfpileupTag;
        edm::InputTag photonsTag;
        edm::InputTag triggerResultsTag;
        std::vector<std::string> triggerPathsVector;
        std::map<std::string, int> triggerPathsMap;
        bool isWorZMCSample;            
        bool isPhotonSample;            
        TTree* tree;
        TH1F* puhist;

        int32_t  puobs, putrue; 
        int32_t  wzid, l1id, l2id, mu1pid, mu2pid, mu1id, mu2id, el1pid, el2pid, el1id, el2id; 
        uint32_t event, run, lumi;
        uint32_t nvtx, nmuons, nelectrons, ntaus, njets;
        uint32_t hltmet120, hltmet95jet80, hltmet105jet80, hltdoublemu, hltsinglemu, hltdoubleel;
        double   pfmet, pfmetphi, t1pfmet, t1pfmetphi, calomet, calometphi, mumet, mumetphi, elmet, elmetphi;
        double   signaljetpt, signaljeteta, signaljetphi, signaljetCHfrac, signaljetNHfrac, signaljetEMfrac, signaljetmetdphi;
        double   secondjetpt, secondjeteta, secondjetphi, secondjetCHfrac, secondjetNHfrac, secondjetEMfrac, secondjetmetdphi;
        double   jetjetdphi;
        double   thirdjetpt, thirdjeteta, thirdjetphi;
        double   wzmass, wzmt, wzpt, wzeta, wzphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
        double   mu1pt, mu1eta, mu1phi, mu1iso, mu2pt, mu2eta, mu2phi, mu2iso, el1pt, el1eta, el1phi, el1iso, el2pt, el2eta, el2phi, el2iso, phpt, pheta, phphi, phmet, phmetphi;
        double   zmass, zpt, zeta, zphi, wmt, emumass, emupt, emueta, emuphi, zeemass, zeept, zeeeta, zeephi, wemt;
        double   wgt, puwgt, weight;
};

MonoJetTreeMaker::MonoJetTreeMaker(const edm::ParameterSet& iConfig): 
    pileupInfoTag(iConfig.getParameter<edm::InputTag>("pileup")),
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
    gensTag((iConfig.existsAs<edm::InputTag>("gens") ? iConfig.getParameter<edm::InputTag>("gens") : edm::InputTag("genParticles"))),
    muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
    electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),
    tausTag(iConfig.getParameter<edm::InputTag>("taus")),
    jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
    pfmetTag(iConfig.getParameter<edm::InputTag>("pfmet")),
    t1pfmetTag(iConfig.getParameter<edm::InputTag>("t1pfmet")),
    calometTag(iConfig.getParameter<edm::InputTag>("calomet")),
    pfpileupTag(iConfig.getParameter<edm::InputTag>("pfpileup")),
    photonsTag((iConfig.existsAs<edm::InputTag>("photons") ? iConfig.getParameter<edm::InputTag>("photons") : edm::InputTag("GAMMA"))),
    triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
    isWorZMCSample((iConfig.existsAs<bool>("isWorZMCSample") ? iConfig.getParameter<bool>("isWorZMCSample") : false)),
    isPhotonSample((iConfig.existsAs<bool>("isPhotonSample") ? iConfig.getParameter<bool>("isPhotonSample") : false)),
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
    if (isWorZMCSample) iEvent.getByLabel(gensTag, gensH);

    Handle<View<Muon> > muonsH;
    iEvent.getByLabel(muonsTag, muonsH);

    Handle<View<GsfElectron> > electronsH;
    iEvent.getByLabel(electronsTag, electronsH);

    Handle<View<pat::Tau> > tausH;
    iEvent.getByLabel(tausTag, tausH);

    Handle<View<pat::Jet> > jetsH;
    iEvent.getByLabel(jetsTag, jetsH);

    Handle<View<MET> > pfmetH;
    iEvent.getByLabel(pfmetTag, pfmetH);

    Handle<View<MET> > t1pfmetH;
    iEvent.getByLabel(t1pfmetTag, t1pfmetH);

    Handle<View<MET> > calometH;
    iEvent.getByLabel(calometTag, calometH);

    Handle<View<PFCandidate> > pfpileupH;
    iEvent.getByLabel(pfpileupTag, pfpileupH);

    Handle<View<Photon> > photonsH;
    if (isPhotonSample) iEvent.getByLabel(photonsTag, photonsH);

    ESHandle<TransientTrackBuilder> trackBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

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
        if (i == 3 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel    = 1; // Double electron triggers
        if (i == 4 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel    = 1; // Double electron triggers
        if (i == 5 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120      = 1; // MET trigger trigger
        if (i == 6 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet95jet80  = 1; // Jet-MET trigger trigger
        if (i == 7 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet105jet80 = 1; // Jet-MET trigger trigger
    }

    // Pileup info
    nvtx   = verticesH->size();
    puobs  = 0;
    putrue = 0;
    puwgt  = 1.;
    if (pileupInfoH.isValid()) {
        for (vector<PileupSummaryInfo>::const_iterator pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
            if (pileupInfo_iter->getBunchCrossing() == 0) {
                puobs  = pileupInfo_iter->getPU_NumInteractions();
                putrue = pileupInfo_iter->getTrueNumInteractions();
                puwgt = puhist->GetBinContent(puhist->FindBin(putrue));
            }
        }
    }
    weight = wgt * puwgt;

    // MET information 
    calomet    = calometH->front().et();
    calometphi = calometH->front().phi();
 
    pfmet      = pfmetH->front().et();
    pfmetphi   = pfmetH->front().phi();
 
    t1pfmet    = t1pfmetH->front().et();
    t1pfmetphi = t1pfmetH->front().phi();
  
    double mumetx = t1pfmet * cos(t1pfmetphi);
    double mumety = t1pfmet * sin(t1pfmetphi);

    for (View<Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        mumetx += muons_iter->pt() * cos(muons_iter->phi());
        mumety += muons_iter->pt() * sin(muons_iter->phi());
    }

    mumet    = sqrt(mumetx*mumetx + mumety*mumety);
    mumetphi = atan2(mumety, mumetx);
  
    double elmetx = t1pfmet * cos(t1pfmetphi);
    double elmety = t1pfmet * sin(t1pfmetphi);

    for (View<GsfElectron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        elmetx += electrons_iter->pt() * cos(electrons_iter->phi());
        elmety += electrons_iter->pt() * sin(electrons_iter->phi());
    }

    elmet    = sqrt(elmetx*elmetx + elmety*elmety);
    elmetphi = atan2(elmety, elmetx);
  
    // Jet information
    pat::JetCollection jets;
    for (View<pat::Jet>::const_iterator jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
        if (jets_iter->electronEnergyFraction() > 0.5 || jets_iter->muonEnergyFraction() > 0.5) continue;
        if (isPhotonSample && jets_iter->photonEnergyFraction() > 0.5) continue;
        if (fabs(jets_iter->eta()) > 4.5) continue;
        bool passjetid = false;
        if (jets_iter->neutralHadronEnergyFraction() < 0.99 && jets_iter->neutralEmEnergyFraction() < 0.99 && jets_iter->getPFConstituents().size() > 1) {
            if (fabs(jets_iter->eta()) > 2.4) passjetid = true;
            else if (fabs(jets_iter->eta()) <= 2.4 && jets_iter->chargedHadronEnergyFraction() > 0. && jets_iter->chargedEmEnergyFraction() < 0.99 && jets_iter->chargedMultiplicity() > 0) passjetid = true;
        }
        if (!passjetid) continue;
        jets.push_back(*jets_iter);
    }

    // Jet information
    njets = 0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > 30) njets++;
    }

    int hardestJetIndex = -1;
    double hardestJetPt = 0.0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > hardestJetPt) {
            hardestJetIndex = i;
            hardestJetPt = jets[i].pt();
        }                        
    }

    int secondJetIndex = -1;
    double secondJetPt = 0.0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (hardestJetIndex >= 0 && i != unsigned(hardestJetIndex) && jets[i].pt() > secondJetPt) {
            secondJetIndex = i;
            secondJetPt = jets[i].pt();
        }                        
    }

    int thirdJetIndex = -1;
    double thirdJetPt = 0.0;
    for (size_t i = 0; i < jets.size(); i++) {
        if (hardestJetIndex >= 0 && secondJetIndex >= 0 && i != unsigned(hardestJetIndex) && i != unsigned(secondJetIndex) && jets[i].pt() > thirdJetPt) {
            thirdJetIndex = i;
            thirdJetPt = jets[i].pt();
        }                        
    }

    signaljetpt      = 0.0;
    signaljeteta     = 0.0;
    signaljetphi     = 0.0;
    signaljetCHfrac  = 0.0;
    signaljetNHfrac  = 0.0;
    signaljetEMfrac  = 0.0;
    signaljetmetdphi = 0.0;
    secondjetpt      = 0.0;
    secondjeteta     = 0.0;
    secondjetphi     = 0.0;
    secondjetCHfrac  = 0.0;
    secondjetNHfrac  = 0.0;
    secondjetEMfrac  = 0.0;
    secondjetmetdphi = 0.0;
    jetjetdphi       = 0.0;
    thirdjetpt       = 0.0;
    thirdjeteta      = 0.0;
    thirdjetphi      = 0.0;

    if (hardestJetIndex >= 0) {
        signaljetpt     = jets[hardestJetIndex].pt();
        signaljeteta    = jets[hardestJetIndex].eta();
        signaljetphi    = jets[hardestJetIndex].phi();
        signaljetCHfrac = jets[hardestJetIndex].chargedHadronEnergyFraction();
        signaljetNHfrac = jets[hardestJetIndex].neutralHadronEnergyFraction();
        signaljetEMfrac = jets[hardestJetIndex].neutralEmEnergyFraction();
    }

    if (secondJetIndex >= 0) {
        secondjetpt     = jets[secondJetIndex].pt();
        secondjeteta    = jets[secondJetIndex].eta();
        secondjetphi    = jets[secondJetIndex].phi();
        secondjetCHfrac = jets[secondJetIndex].chargedHadronEnergyFraction();
        secondjetNHfrac = jets[secondJetIndex].neutralHadronEnergyFraction();
        secondjetEMfrac = jets[secondJetIndex].neutralEmEnergyFraction();
    }

    if (thirdJetIndex >= 0) {
        thirdjetpt     = jets[thirdJetIndex].pt();
        thirdjeteta    = jets[thirdJetIndex].eta();
        thirdjetphi    = jets[thirdJetIndex].phi();
    }

    if (signaljetpt > 0.0 && secondjetpt > 0.0) jetjetdphi       = deltaPhi(signaljetphi, secondjetphi);
    if (signaljetpt > 0.0 && secondjetpt > 0.0) signaljetmetdphi = deltaPhi(signaljetphi, mumetphi);
    if (signaljetpt > 0.0 && secondjetpt > 0.0) secondjetmetdphi = deltaPhi(secondjetphi, mumetphi);

    // Lepton counts
    nmuons      = muonsH->size();
    nelectrons  = electronsH->size();
    ntaus       = tausH->size();

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

    if (isWorZMCSample && gensH.isValid()) {
        for (View<GenParticle>::const_iterator gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
            if ((gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && gens_iter->numberOfDaughters() == 3 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) {
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
    mu1id       = 0;
    mu1iso      = 0.0;
    mu2pid      = 0;
    mu2pt       = 0.0;
    mu2eta      = 0.0; 
    mu2phi      = 0.0;
    mu2id       = 0;
    mu2iso      = 0.0;
    el1pid      = 0;
    el1pt       = 0.0;
    el1eta      = 0.0;
    el1phi      = 0.0;
    el1id       = 0;
    el1iso      = 0.0;
    el2pid      = 0;
    el2pt       = 0.0;
    el2eta      = 0.0;
    el2phi      = 0.0;
    el2id       = 0;
    el2iso      = 0.0;

    if (muonsH->size() == 1 || muonsH->size() == 2) {
        View<Muon>::const_iterator muons_iter = muonsH->begin();
        mu1pid = muons_iter->pdgId(); 
        mu1pt  = muons_iter->pt(); 
        mu1eta = muons_iter->eta(); 
        mu1phi = muons_iter->phi();

        double mu1tip = 1000.0;
        if (muons_iter->innerTrack().isAvailable() && verticesH->size() > 0) {
            TransientTrack ttrack = (*trackBuilder).build(*(muons_iter->innerTrack()));
            mu1tip = IPTools::absoluteTransverseImpactParameter(ttrack, verticesH->front()).second.value();
        }
    
        if (muons_iter->pt() > 20 && fabs(muons_iter->eta()) < 2.4 && muons_iter->isGlobalMuon() && muons_iter->isTrackerMuon() &&
            muons_iter->globalTrack()->normalizedChi2() < 10 &&
            muons_iter->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
            muons_iter->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
            muons_iter->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
            muons_iter->numberOfMatchedStations() > 1 &&
            fabs(mu1tip) < 0.02) mu1id = 1;    
    
        mu1iso = (muons_iter->pfIsolationR04().sumChargedHadronPt + muons_iter->pfIsolationR04().sumNeutralHadronEt + muons_iter->pfIsolationR04().sumPhotonEt - 0.5*muons_iter->pfIsolationR04().sumPUPt)/muons_iter->pt();
    
        if (muonsH->size() == 1) wmt = sqrt(2.0 * mu1pt * t1pfmet * (1.0 - cos(deltaPhi(mu1phi, t1pfmetphi))));
    }
   

 
    if (muonsH->size() == 2) {        
        View<Muon>::const_iterator muons_iter = muonsH->begin();
        ++muons_iter;
        mu2pid = muons_iter->pdgId(); 
        mu2pt  = muons_iter->pt(); 
        mu2eta = muons_iter->eta(); 
        mu2phi = muons_iter->phi();
    
        double mu2tip = 1000.0;
        if (muons_iter->innerTrack().isAvailable() && verticesH->size() > 0) {
            TransientTrack ttrack = (*trackBuilder).build(*(muons_iter->innerTrack()));
            mu2tip = IPTools::absoluteTransverseImpactParameter(ttrack, verticesH->front()).second.value();
        }

        if (muons_iter->pt() > 20 && fabs(muons_iter->eta()) < 2.4 && muons_iter->isGlobalMuon() && muons_iter->isTrackerMuon() &&
            muons_iter->globalTrack()->normalizedChi2() < 10 &&
            muons_iter->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
            muons_iter->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
            muons_iter->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
            muons_iter->numberOfMatchedStations() > 1 &&
            fabs(mu2tip) < 0.02) mu2id = 1;
    
        mu2iso = (muons_iter->pfIsolationR04().sumChargedHadronPt + muons_iter->pfIsolationR04().sumNeutralHadronEt + muons_iter->pfIsolationR04().sumPhotonEt - 0.5*muons_iter->pfIsolationR04().sumPUPt)/muons_iter->pt();
    
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonsH->begin()->p());
        TLorentzVector mu2vec; mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muons_iter->p());
    
        TLorentzVector zvec(mu1vec);
        zvec += mu2vec;
    
        zmass = zvec.M();
        zpt   = zvec.Pt();
        zeta  = zvec.Eta();            
        zphi  = zvec.Phi();
    }

    if (nelectrons == 1 || nelectrons == 2) {
        View<GsfElectron>::const_iterator electrons_iter = electronsH->begin();
        el1pid = electrons_iter->pdgId();
        el1pt  = electrons_iter->pt();
        el1eta = electrons_iter->eta();
        el1phi = electrons_iter->phi();
        
        double el1tip = 1000.0;
        if (electrons_iter->gsfTrack().isAvailable() && verticesH->size() > 0) {
            TransientTrack ttrack = (*trackBuilder).build(electrons_iter->gsfTrack());
            el1tip = IPTools::absoluteTransverseImpactParameter(ttrack, verticesH->front()).second.value();
        }

        if (electrons_iter->pt() > 20 && fabs(electrons_iter->superCluster()->eta()) < 2.5 && electrons_iter->gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 && 
            (fabs(electrons_iter->convDist()) > 0.02 || fabs(electrons_iter->convDcot()) > 0.02) &&
            ((electrons_iter->isEB() && electrons_iter->sigmaIetaIeta() < 0.01 && electrons_iter->hadronicOverEm() < 0.040 &&
              fabs(electrons_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.06 && fabs(electrons_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.004) || 
             (electrons_iter->isEE() && electrons_iter->sigmaIetaIeta() < 0.03 && electrons_iter->hadronicOverEm() < 0.025 &&
              fabs(electrons_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.03 && fabs(electrons_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.007)) && 
            fabs(el1tip) < 0.02) el1id = 1;

        el1iso = electrons_iter->pfIsolationVariables().neutralHadronIso;
        el1iso += electrons_iter->pfIsolationVariables().photonIso;
        for (View<PFCandidate>::const_iterator pfpileup_iter = pfpileupH->begin(); pfpileup_iter != pfpileupH->end(); ++pfpileup_iter) {
            if (pfpileup_iter->charge() == 0) continue;
            if (deltaR(pfpileup_iter->eta(), pfpileup_iter->phi(), electrons_iter->eta(), electrons_iter->phi()) > 0.4) continue;
            el1iso -= 0.5*(pfpileup_iter->pt());
        }
        if (el1iso < 0.) el1iso = 0.;
        el1iso += electrons_iter->pfIsolationVariables().chargedHadronIso;
        el1iso /= electrons_iter->pt();   
 
        if (electronsH->size() == 1) wemt = sqrt(2.0 * el1pt * t1pfmet * (1.0 - cos(deltaPhi(el1phi, t1pfmetphi))));
    }

    if (nelectrons == 2) {
        View<GsfElectron>::const_iterator electrons_iter = electronsH->begin();
        ++electrons_iter;
        el2pid = electrons_iter->pdgId();
        el2pt  = electrons_iter->pt();
        el2eta = electrons_iter->eta();
        el2phi = electrons_iter->phi();

        double el2tip = 1000.0;
        if (electrons_iter->gsfTrack().isAvailable() && verticesH->size() > 0) {
            TransientTrack ttrack = (*trackBuilder).build(*(electrons_iter->gsfTrack()));
            el2tip = IPTools::absoluteTransverseImpactParameter(ttrack, verticesH->front()).second.value();
        }

        if (electrons_iter->pt() > 20 && fabs(electrons_iter->superCluster()->eta()) < 2.5 && electrons_iter->gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 && 
            (fabs(electrons_iter->convDist()) > 0.02 || fabs(electrons_iter->convDcot()) > 0.02) &&
            ((electrons_iter->isEB() && electrons_iter->sigmaIetaIeta() < 0.01 && electrons_iter->hadronicOverEm() < 0.040 &&
              fabs(electrons_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.06 && fabs(electrons_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.004) || 
             (electrons_iter->isEE() && electrons_iter->sigmaIetaIeta() < 0.03 && electrons_iter->hadronicOverEm() < 0.025 &&
              fabs(electrons_iter->deltaPhiSuperClusterTrackAtVtx()) < 0.03 && fabs(electrons_iter->deltaEtaSuperClusterTrackAtVtx()) < 0.007)) && 
            fabs(el2tip) < 0.02) el2id = 1;

        el2iso = electrons_iter->pfIsolationVariables().neutralHadronIso;
        el2iso += electrons_iter->pfIsolationVariables().photonIso;
        for (View<PFCandidate>::const_iterator pfpileup_iter = pfpileupH->begin(); pfpileup_iter != pfpileupH->end(); ++pfpileup_iter) {
            if (pfpileup_iter->charge() == 0) continue;
            if (deltaR(pfpileup_iter->eta(), pfpileup_iter->phi(), electrons_iter->eta(), electrons_iter->phi()) > 0.4) continue;
            el2iso -= 0.5*(pfpileup_iter->pt());
        }
        if (el2iso < 0.) el2iso = 0.;
        el2iso += electrons_iter->pfIsolationVariables().chargedHadronIso;
        el2iso /= electrons_iter->pt();   

        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electronsH->begin()->p());
        TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electrons_iter->p());

        TLorentzVector zvec(el1vec);
        zvec += el2vec;

        zeemass = zvec.M();
        zeept   = zvec.Pt();
        zeeeta  = zvec.Eta();
        zeephi  = zvec.Phi();
    }

    if (muonsH->size() == 1 && nelectrons == 1) {
        TLorentzVector mu1vec; mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonsH->begin()->p());
        TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electronsH->begin()->p());
        
        TLorentzVector emuvec(mu1vec);
        emuvec += el1vec;
        
        emumass = emuvec.M();
        emupt   = emuvec.Pt();
        emueta  = emuvec.Eta();
        emuphi  = emuvec.Phi();
    } 

    // Photon information
    phpt     = 0.0;
    pheta    = 0.0;
    phphi    = 0.0;
    phmet    = 0.0;
    phmetphi = 0.0;

    if (isPhotonSample && photonsH.isValid() && photonsH->size() == 1) {
        phpt    = photonsH->begin()->pt();
        pheta   = photonsH->begin()->eta();
        phphi   = photonsH->begin()->phi();

        double phmetx = t1pfmet * cos(t1pfmetphi);
        double phmety = t1pfmet * sin(t1pfmetphi);
        phmetx += photonsH->begin()->et() * cos(photonsH->begin()->phi());
        phmety += photonsH->begin()->et() * sin(photonsH->begin()->phi());

        phmet    = sqrt(phmetx*phmetx + phmety*phmety);
        phmetphi = atan2(phmety, phmetx);
    }

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
    tree->Branch("pfmet"           , &pfmet           , "pfmet/D");
    tree->Branch("pfmetphi"        , &pfmetphi        , "pfmetphi/D");
    tree->Branch("t1pfmet"         , &t1pfmet         , "t1pfmet/D");
    tree->Branch("t1pfmetphi"      , &t1pfmetphi      , "t1pfmetphi/D");
    tree->Branch("calomet"         , &calomet         , "calomet/D");
    tree->Branch("calometphi"      , &calometphi      , "calometphi/D");
    tree->Branch("mumet"           , &mumet           , "mumet/D");
    tree->Branch("mumetphi"        , &mumetphi        , "mumetphi/D");
    tree->Branch("elmet"           , &elmet           , "elmet/D");
    tree->Branch("elmetphi"        , &elmetphi        , "elmetphi/D");
    tree->Branch("signaljetpt"     , &signaljetpt     , "signaljetpt/D");
    tree->Branch("signaljeteta"    , &signaljeteta    , "signaljeteta/D");
    tree->Branch("signaljetphi"    , &signaljetphi    , "signaljetphi/D");
    tree->Branch("signaljetCHfrac" , &signaljetCHfrac , "signaljetCHfrac/D");
    tree->Branch("signaljetNHfrac" , &signaljetNHfrac , "signaljetNHfrac/D");
    tree->Branch("signaljetEMfrac" , &signaljetEMfrac , "signaljetEMfrac/D");
    tree->Branch("signaljetmetdphi", &signaljetmetdphi, "signaljetmetdphi/D");
    tree->Branch("secondjetpt"     , &secondjetpt     , "secondjetpt/D");
    tree->Branch("secondjeteta"    , &secondjeteta    , "secondjeteta/D");
    tree->Branch("secondjetphi"    , &secondjetphi    , "secondjetphi/D");
    tree->Branch("secondjetCHfrac" , &secondjetCHfrac , "secondjetCHfrac/D");
    tree->Branch("secondjetNHfrac" , &secondjetNHfrac , "secondjetNHfrac/D");
    tree->Branch("secondjetEMfrac" , &secondjetEMfrac , "secondjetEMfrac/D");
    tree->Branch("secondjetmetdphi", &secondjetmetdphi, "secondjetmetdphi/D");
    tree->Branch("jetjetdphi"      , &jetjetdphi      , "jetjetdphi/D");
    tree->Branch("thirdjetpt"      , &thirdjetpt      , "thirdjetpt/D");
    tree->Branch("thirdjeteta"     , &thirdjeteta     , "thirdjeteta/D");
    tree->Branch("thirdjetphi"     , &thirdjetphi     , "thirdjetphi/D");
    tree->Branch("wgt"             , &wgt             , "wgt/D");
    tree->Branch("puwgt"           , &puwgt           , "puwgt/D");
    tree->Branch("weight"          , &weight          , "weight/D");
    tree->Branch("hltmet120"       , &hltmet120       , "hltmet120/i");
    tree->Branch("hltmet95jet80"   , &hltmet95jet80   , "hltmet95jet80/i");
    tree->Branch("hltmet105jet80"  , &hltmet105jet80  , "hltmet105jet80/i");
    tree->Branch("hltdoublemu"     , &hltdoublemu     , "hltdoublemu/i");
    tree->Branch("hltsinglemu"     , &hltsinglemu     , "hltsinglemu/i");
    tree->Branch("hltdoubleel"     , &hltdoubleel     , "hltdoubleel/i");
    tree->Branch("zmass"           , &zmass           , "zmass/D");
    tree->Branch("zpt"             , &zpt             , "zpt/D");
    tree->Branch("zeta"            , &zeta            , "zeta/D");
    tree->Branch("zphi"            , &zphi            , "zphi/D");
    tree->Branch("wmt"             , &wmt             , "wmt/D");
    tree->Branch("emumass"         , &emumass         , "emumass/D");
    tree->Branch("emupt"           , &emupt           , "emupt/D");
    tree->Branch("emueta"          , &emueta          , "emueta/D");
    tree->Branch("emuphi"          , &emuphi          , "emuphi/D");
    tree->Branch("zeemass"         , &zeemass         , "zeemass/D");
    tree->Branch("zeept"           , &zeept           , "zeeept/D");
    tree->Branch("zeeeta"          , &zeeeta          , "zeeeta/D");
    tree->Branch("zeephi"          , &zeephi          , "zeephi/D");
    tree->Branch("wemt"            , &wemt            , "wemt/D");
    tree->Branch("mu1pid"          , &mu1pid          , "mu1pid/I");
    tree->Branch("mu1pt"           , &mu1pt           , "mu1pt/D");
    tree->Branch("mu1eta"          , &mu1eta          , "mu1eta/D");
    tree->Branch("mu1phi"          , &mu1phi          , "mu1phi/D");
    tree->Branch("mu1id"           , &mu1id           , "mu1id/I");
    tree->Branch("mu1iso"          , &mu1iso          , "mu1iso/D");
    tree->Branch("mu2pid"          , &mu2pid          , "mu2pid/I");
    tree->Branch("mu2pt"           , &mu2pt           , "mu2pt/D");
    tree->Branch("mu2eta"          , &mu2eta          , "mu2eta/D");
    tree->Branch("mu2phi"          , &mu2phi          , "mu2phi/D");
    tree->Branch("mu2id"           , &mu2id           , "mu2id/I");
    tree->Branch("mu2iso"          , &mu2iso          , "mu2iso/D");
    tree->Branch("el1pid"          , &el1pid          , "el1pid/I");
    tree->Branch("el1pt"           , &el1pt           , "el1pt/D");
    tree->Branch("el1eta"          , &el1eta          , "el1eta/D");
    tree->Branch("el1phi"          , &el1phi          , "el1phi/D");
    tree->Branch("el1id"           , &el1id           , "el1id/I");
    tree->Branch("el1iso"          , &el1iso          , "el1iso/D");
    tree->Branch("el2pid"          , &el2pid          , "el2pid/I");
    tree->Branch("el2pt"           , &el2pt           , "el2pt/D");
    tree->Branch("el2eta"          , &el2eta          , "el2eta/D");
    tree->Branch("el2phi"          , &el2phi          , "el2phi/D");
    tree->Branch("el2id"           , &el2id           , "el2id/I");
    tree->Branch("el2iso"          , &el2iso          , "el2iso/D");
    
    if (isWorZMCSample) {
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
    }

    if (isPhotonSample) {
    tree->Branch("phpt"            , &phpt            , "phpt/D");
    tree->Branch("pheta"           , &pheta           , "pheta/D");
    tree->Branch("phphi"           , &phphi           , "phphi/D");
    tree->Branch("phmet"           , &phmet           , "phmet/D");
    tree->Branch("phmetphi"        , &phmetphi        , "phmetphi/D");
    }
}

void MonoJetTreeMaker::endJob() {
}

void MonoJetTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    triggerPathsVector.push_back("HLT_Mu17_Mu8");
    triggerPathsVector.push_back("HLT_Mu17_TkMu8");
    triggerPathsVector.push_back("HLT_IsoMu24_eta2p1");
    triggerPathsVector.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL"),
    triggerPathsVector.push_back("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL"),
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
