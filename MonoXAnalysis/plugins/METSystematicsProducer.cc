#ifndef METSystematicsProducer_h
#define METSystematicsProducer_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/CorrMETData.h"
#include "DataFormats/METReco/interface/MET.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "PhysicsTools/PatUtils/interface/PATJetCorrExtractor.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TFormula.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"


class METSystematicsProducer : public edm::stream::EDProducer<> {
public:
  explicit METSystematicsProducer(const edm::ParameterSet&);
  virtual ~METSystematicsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&) override;

  template<typename T>
  reco::Candidate::LorentzVector findParticle(const T & particle, const edm::View<reco::Candidate> & pfCandCollection);

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // input collections
  const edm::EDGetTokenT<pat::METCollection>    metToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> >      jetToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> >     muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
  edm::EDGetTokenT<std::vector<pat::Photon> >   photonToken_;
  edm::EDGetTokenT<std::vector<pat::Tau> >      tauToken_;
  const edm::EDGetTokenT<edm::View<reco::Candidate> > pfCandToken_;
  const edm::EDGetTokenT<double>  rhoToken_;
  edm::EDGetTokenT<reco::JetCorrector> jetCorrToken_; 
  edm::EDGetTokenT<reco::JetCorrector> jetCorrTokenRes_; 

  // PSet for the object uncertainties
  const edm::ParameterSet electronPSet_;
  const edm::ParameterSet muonPSet_;
  const edm::ParameterSet photonPSet_;
  const edm::ParameterSet tauPSet_;
  const edm::ParameterSet jetPSet_;
  const edm::ParameterSet unclusteredPSet_;
  const edm::InputTag     inputMET_;

  // options
  bool storeSmearedShiftedCollections_;
  bool skipMuon_;
  bool skipElectron_;
  bool skipPhoton_;
  bool skipTau_;
  bool skipJet_;

  // in case of binning and object selections --> from PSET
  std::vector<StringCutObjectSelector<pat::Muon> > muonSelection_;
  std::vector<float> muonScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Electron> > electronSelection_;
  std::vector<float> electronScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Photon> > photonSelection_;
  std::vector<float> photonScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Tau> > tauSelection_;
  std::vector<float> tauScaleUnc_;

  // JES external file
  edm::FileInPath jetJECUncFile_;
  StringCutObjectSelector<pat::Jet> jetSelection_;

  // JER and JERSF external files
  edm::FileInPath jetJERFile_;
  edm::FileInPath jetJERSFFile_;
  // Random numbers
  std::auto_ptr<TRandom3> rand_;

  // unclustered
  std::vector<StringCutObjectSelector<reco::Candidate> > unclusteredSelection_;
  std::vector<std::string> unclusteredUnc_;

  // sorting
  template<typename T>
  class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i.pt() > j.pt());
    }

  };

  PatPtSorter<pat::Jet>      jetSorter;
  PatPtSorter<pat::Muon>     muonSorter;
  PatPtSorter<pat::Electron> electronSorter;
  PatPtSorter<pat::Photon>   photonSorter;
  PatPtSorter<pat::Tau>      tauSorter;

};

#endif

METSystematicsProducer::~METSystematicsProducer(){}

METSystematicsProducer::METSystematicsProducer(const edm::ParameterSet& iConfig):
  metToken_                 (consumes<pat::METCollection> (iConfig.getParameter<edm::InputTag>("inputMET"))),
  pfCandToken_              (consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("pfCandidate"))),
  rhoToken_                 (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  electronPSet_             (iConfig.getParameter<edm::ParameterSet>("electron")),
  muonPSet_                 (iConfig.getParameter<edm::ParameterSet>("muon")),
  photonPSet_               (iConfig.getParameter<edm::ParameterSet>("photon")),
  tauPSet_                  (iConfig.getParameter<edm::ParameterSet>("tau")),
  jetPSet_                  (iConfig.getParameter<edm::ParameterSet>("jet")),
  unclusteredPSet_          (iConfig.getParameter<edm::ParameterSet>("unclustered")),
  inputMET_                 (iConfig.getParameter<edm::InputTag>("inputMET")),
  storeSmearedShiftedCollections_(iConfig.existsAs<bool>("storeSmearedShiftedCollections") ? iConfig.getParameter<bool>("storeSmearedShiftedCollections") : false),
  skipMuon_(iConfig.existsAs<bool>("Muon") ? iConfig.getParameter<bool>("skipMuon") : false),
  skipElectron_(iConfig.existsAs<bool>("Electron") ? iConfig.getParameter<bool>("skiElectron") : false),
  skipPhoton_(iConfig.existsAs<bool>("skipPhoton") ? iConfig.getParameter<bool>("skipPhoton") : false),
  skipTau_(iConfig.existsAs<bool>("skipTau") ? iConfig.getParameter<bool>("skipTau") : false),
  skipJet_(iConfig.existsAs<bool>("skipJet") ? iConfig.getParameter<bool>("skipJet") : false),
  jetSelection_             (StringCutObjectSelector<pat::Jet>(jetPSet_.getParameter<std::string>("selection")))
 {

   jetToken_      = consumes<std::vector<pat::Jet> >  (jetPSet_.getParameter<edm::InputTag>("src"));
   tauToken_      = consumes<std::vector<pat::Tau> >  (tauPSet_.getParameter<edm::InputTag>("src"));
   photonToken_   = consumes<std::vector<pat::Photon> >  (photonPSet_.getParameter<edm::InputTag>("src"));
   muonToken_     = consumes<std::vector<pat::Muon> > (muonPSet_.getParameter<edm::InputTag>("src"));
   electronToken_ = consumes<std::vector<pat::Electron> >  (electronPSet_.getParameter<edm::InputTag>("src"));
   // take correction if exist
   if(jetPSet_.existsAs<edm::InputTag>("jetCorrLabel"))
     jetCorrToken_ = mayConsume<reco::JetCorrector>(jetPSet_.getParameter<edm::InputTag>("jetCorrLabel"));
   if(jetPSet_.existsAs<edm::InputTag>("jetCorrLabelRes"))
     jetCorrTokenRes_ = mayConsume<reco::JetCorrector>(jetPSet_.getParameter<edm::InputTag>("jetCorrLabelRes")); 

  // selector definitions
  for(auto bin : muonPSet_.getParameter<std::vector<edm::ParameterSet> >("binning")){
    muonSelection_.push_back(StringCutObjectSelector<pat::Muon>(bin.getParameter<std::string>("binSelection")));
    muonScaleUnc_.push_back(bin.getParameter<double>("uncertainty"));
  }

  for(auto bin : electronPSet_.getParameter<std::vector<edm::ParameterSet> >("binning")){
    electronSelection_.push_back(StringCutObjectSelector<pat::Electron>(bin.getParameter<std::string>("binSelection")));
    electronScaleUnc_.push_back(bin.getParameter<double>("uncertainty"));
  }

  for(auto bin : photonPSet_.getParameter<std::vector<edm::ParameterSet> >("binning")){
    photonSelection_.push_back(StringCutObjectSelector<pat::Photon>(bin.getParameter<std::string>("binSelection")));
    photonScaleUnc_.push_back(bin.getParameter<double>("uncertainty"));
  }

  for(auto bin : tauPSet_.getParameter<std::vector<edm::ParameterSet> >("binning")){
    tauSelection_.push_back(StringCutObjectSelector<pat::Tau>(bin.getParameter<std::string>("binSelection")));
    tauScaleUnc_.push_back(bin.getParameter<double>("uncertainty"));
  }

  for(auto bin : unclusteredPSet_.getParameter<std::vector<edm::ParameterSet> >("binning")){
    unclusteredSelection_.push_back(StringCutObjectSelector<reco::Candidate>(bin.getParameter<std::string>("binSelection")));
    unclusteredUnc_.push_back(bin.getParameter<std::string>("binUncertainty"));
  }

  // file for JEC unc
  jetJECUncFile_ = jetPSet_.getParameter<edm::FileInPath>("JECUncFile");
  if(not jetJECUncFile_.location() and jetPSet_.getParameter<bool> ("useExternalJECUncertainty") and not skipJet_)
    throw cms::Exception("METSystematicsProducer") << " Failed to find File = " << jetJERFile_ << " !!\n";
  
  // file for JER truth
  jetJERFile_    = jetPSet_.getParameter<edm::FileInPath>("JERFile");
  jetJERSFFile_  = jetPSet_.getParameter<edm::FileInPath>("JERSFFile");
  
  if(not jetJERFile_.location() and jetPSet_.getParameter<bool> ("useExternalJER") and not skipJet_)
    throw cms::Exception("METSystematicsProducer") << " Failed to find File = " << jetJERFile_ << " !!\n";

  if(not jetJERSFFile_.location() and jetPSet_.getParameter<bool> ("useExternalJERSF") and not skipJet_)
    throw cms::Exception("METSystematicsProducer") << " Failed to find File = " << jetJERSFFile_ << " !!\n";

  

  // produce shifed and smearead particles
  produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label()));

  if (skipMuon_ == false){
    if(storeSmearedShiftedCollections_){
      produces<pat::MuonCollection>(std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::MuonCollection>(std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"MuonEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"MuonEnDown");
  }

  if(not skipElectron_){
    if(storeSmearedShiftedCollections_){
      produces<pat::ElectronCollection>(std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::ElectronCollection>(std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"ElectronEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"ElectronEnDown");
  }
  
  if(not skipPhoton_){
    if(storeSmearedShiftedCollections_){
      produces<pat::PhotonCollection>(std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::PhotonCollection>(std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"PhotonEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"PhotonEnDown");
  }

  if(not skipTau_){
    if(storeSmearedShiftedCollections_){
      produces<pat::TauCollection>(std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::TauCollection>(std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"TauEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"TauEnDown");
  }

  if(not skipJet_){
    if(storeSmearedShiftedCollections_){
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");    
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"Smear");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResUp");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetEnDown");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetResUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetResDown");    
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"Smear");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"SmearJetResUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"SmearJetResDown");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"UnclusteredEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"UnclusteredEnDown");
  }

  rand_ = std::auto_ptr<TRandom3>(new TRandom3());
  rand_->SetSeed(0);
}


void METSystematicsProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup){

  // get rho
  edm::Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);
  double rho = *rhoH;

  //get original MET collection
  edm::Handle<pat::METCollection> metCollection;
  iEvent.getByToken(metToken_, metCollection);  

  //get pfCands
  edm::Handle<edm::View<reco::Candidate> > pfCandCollection;
  iEvent.getByToken(pfCandToken_, pfCandCollection);
  
  //get objects
  edm::Handle<std::vector<pat::Jet> > jetColl;
  iEvent.getByToken(jetToken_, jetColl);
  edm::Handle<reco::JetCorrector> jetCorr;
  iEvent.getByToken(jetCorrToken_, jetCorr);
  edm::Handle<reco::JetCorrector> jetCorrRes;
  if(iEvent.isRealData())
    iEvent.getByToken(jetCorrTokenRes_, jetCorrRes);

  edm::Handle<std::vector<pat::Electron> > eleColl;
  iEvent.getByToken(electronToken_, eleColl);

  edm::Handle<std::vector<pat::Muon> > muColl;
  iEvent.getByToken(muonToken_, muColl);

  edm::Handle<std::vector<pat::Photon> > phColl;
  iEvent.getByToken(photonToken_, phColl);

  edm::Handle<std::vector<pat::Tau> > tauColl;
  iEvent.getByToken(tauToken_, tauColl);

  // create output autoPtr
  std::auto_ptr<pat::METCollection> MET (new pat::METCollection);
  MET->push_back((*metCollection)[0]);
  iEvent.put(MET,inputMET_.label());

  std::auto_ptr<pat::MuonCollection> muonEnUp(new pat::MuonCollection);
  std::auto_ptr<pat::MuonCollection> muonEnDown(new pat::MuonCollection);
  CorrMETData corrMuonEnUp;
  CorrMETData corrMuonEnDown;
  std::auto_ptr<pat::METCollection> metMuonEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metMuonEnDown(new pat::METCollection);

  std::auto_ptr<pat::ElectronCollection> electronEnUp(new pat::ElectronCollection);
  std::auto_ptr<pat::ElectronCollection> electronEnDown(new pat::ElectronCollection);
  CorrMETData corrElectronEnUp;
  CorrMETData corrElectronEnDown;
  std::auto_ptr<pat::METCollection> metElectronEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metElectronEnDown(new pat::METCollection);

  std::auto_ptr<pat::PhotonCollection> photonEnUp(new pat::PhotonCollection);
  std::auto_ptr<pat::PhotonCollection> photonEnDown(new pat::PhotonCollection);
  CorrMETData corrPhotonEnUp;
  CorrMETData corrPhotonEnDown;
  std::auto_ptr<pat::METCollection> metPhotonEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metPhotonEnDown(new pat::METCollection);

  std::auto_ptr<pat::TauCollection> tauEnUp(new pat::TauCollection);
  std::auto_ptr<pat::TauCollection> tauEnDown(new pat::TauCollection);
  CorrMETData corrTauEnUp;
  CorrMETData corrTauEnDown;
  std::auto_ptr<pat::METCollection> metTauEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metTauEnDown(new pat::METCollection);

  std::auto_ptr<pat::JetCollection> jetEnUp(new pat::JetCollection);
  std::auto_ptr<pat::JetCollection> jetEnDown(new pat::JetCollection);
  std::auto_ptr<pat::JetCollection> jetSmear(new pat::JetCollection);
  std::auto_ptr<pat::JetCollection> jetSmearResUp(new pat::JetCollection);
  std::auto_ptr<pat::JetCollection> jetSmearResDown(new pat::JetCollection);
  CorrMETData corrJetEnUp;
  CorrMETData corrJetEnDown;
  CorrMETData corrJetResUp;
  CorrMETData corrJetResDown;
  CorrMETData corrJetSmear;
  CorrMETData corrJetSmearResUp;
  CorrMETData corrJetSmearResDown;
  std::auto_ptr<pat::METCollection> metJetEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetEnDown(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetResUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetResDown(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmearResUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmearResDown(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmear(new pat::METCollection);

  CorrMETData corrUnclusteredEnUp;
  CorrMETData corrUnclusteredEnDown;
  std::auto_ptr<pat::METCollection> metUnclusteredEnUp (new pat::METCollection);
  std::auto_ptr<pat::METCollection> metUnclusteredEnDown (new pat::METCollection);

  // boolean for puppi MET
  bool isPuppiMET = false;
  if(TString(inputMET_.label()).Contains("Puppi") or TString(inputMET_.label()).Contains("puppi") or TString(inputMET_.label()).Contains("PUPPI")) 
    isPuppiMET = true;

  // link from high level objects to pfCandidates
  std::vector<reco::CandidatePtr> particlesInMuon;
  std::vector<reco::CandidatePtr> particlesInElectron;
  std::vector<reco::CandidatePtr> particlesInPhoton;
  std::vector<reco::CandidatePtr> particlesInTau;
  std::vector<reco::CandidatePtr> particlesInJet;

  // loop on muons
  if(not skipMuon_){
    for(auto muon : *muColl){    
      if(isPuppiMET){
	// find the particle in the puppi pfCandidate rescaling correctly the momentum
	reco::Candidate::LorentzVector p4 = findParticle(muon,*pfCandCollection);
	if(p4.pt() == 0)
	  continue;
	muon.setP4(p4);
      }

      // store PF candidates connected
      for(size_t ipart = 0 ; ipart < muon.numberOfSourceCandidatePtrs(); ipart++){
	if(muon.sourceCandidatePtr(ipart).isNonnull() and muon.sourceCandidatePtr(ipart).isAvailable()){
	  particlesInMuon.push_back(muon.sourceCandidatePtr(ipart));
	}
      }

      float ptError = 0.;
      // used numeric uncertainties provided at python level and cloned from: http://cmslxr.fnal.gov/source/PhysicsTools/PatUtils/python/tools/runMETCorrsAndUncertainties.py
      if(muonPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( muonScaleUnc_.size() == muonSelection_.size()){
	  for(auto selection : muonSelection_){	    
	    if(selection(muon)){
	      ptError = muonScaleUnc_.at(iBin);	    	      
	      break;
	    }
	    else iBin++;
	  }	
	}
      }
      else{	  
	//use stadard ptError or the tunePMuonBestTrack
	reco::Track muonTrack;
	if(muon.tunePMuonBestTrack().isNonnull())
	   muonTrack = *(muon.tunePMuonBestTrack().get());      
	else if(muon.muonBestTrack().isNonnull())
	  muonTrack = *(muon.muonBestTrack().get());     
	else
	  muonTrack = *(muon.globalTrack().get());	     
	
	ptError = muonTrack.ptError()/(muon.pt());	  
      }

      // build vector scaled on the transverse plane
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(muon.pt()*(1+ptError),muon.phi());
      p2Down.SetMagPhi(muon.pt()*(1-ptError),muon.phi());
      corrMuonEnUp.mex += (-p2Up.Px()+muon.px());
      corrMuonEnUp.mey += (-p2Up.Py()+muon.py());
      corrMuonEnUp.sumet += fabs(p2Up.Mod()-muon.pt());
      corrMuonEnDown.mex += (-p2Down.Px()+muon.px());
      corrMuonEnDown.mey += (-p2Down.Py()+muon.py());
      corrMuonEnDown.sumet += fabs(p2Down.Mod()-muon.pt());

      // appoximation scaling the whole p4
      pat::Muon newMuonEnUp(muon);
      newMuonEnUp.setP4(muon.p4()*(1+ptError));
      pat::Muon newMuonEnDown(muon);
      newMuonEnDown.setP4(muon.p4()*(1-ptError));

      muonEnUp->push_back(newMuonEnUp);
      muonEnDown->push_back(newMuonEnDown);      
    }

    // correct met --> loop on met collection of input
    for(auto met : *metCollection){
      //build corrected MET
      reco::MET corrMETMuonEnUp(met.sumEt()+corrMuonEnUp.sumet,
		       reco::Candidate::LorentzVector(met.px()+corrMuonEnUp.mex, met.py()+corrMuonEnUp.mey, 0., 
		       sqrt((met.px()+corrMuonEnUp.mex)*(met.px()+corrMuonEnUp.mex)+(met.py()+corrMuonEnUp.mey)*(met.py()+corrMuonEnUp.mey))),met.vertex());

      reco::MET corrMETMuonEnDown(met.sumEt()-corrMuonEnDown.sumet,
	       reco::Candidate::LorentzVector(met.px()+corrMuonEnDown.mex, met.py()+corrMuonEnDown.mey, 0., 
	       sqrt((met.px()+corrMuonEnDown.mex)*(met.px()+corrMuonEnDown.mex)+(met.py()+corrMuonEnDown.mey)*(met.py()+corrMuonEnDown.mey))),met.vertex());

      pat::MET metMuonUp(corrMETMuonEnUp,met);
      pat::MET metMuonDown(corrMETMuonEnDown,met);
      metMuonEnUp->push_back(metMuonUp);
      metMuonEnDown->push_back(metMuonDown);
    }

    if(storeSmearedShiftedCollections_){
      std::sort(muonEnUp->begin(),muonEnUp->end(),muonSorter);
      iEvent.put(muonEnUp,std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      std::sort(muonEnDown->begin(),muonEnDown->end(),muonSorter);
      iEvent.put(muonEnDown,std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metMuonEnUp,std::string(inputMET_.label())+"MuonEnUp");
    iEvent.put(metMuonEnDown,std::string(inputMET_.label())+"MuonEnDown");    
  }


  // loop on electrons
  if(not skipElectron_){
    for(auto ele : *eleColl){    
      if(isPuppiMET){
	reco::Candidate::LorentzVector p4 = findParticle(ele,*pfCandCollection);
	if(p4.pt() == 0)
	  continue;
	ele.setP4(p4);
      }

      for(size_t ipart = 0 ; ipart < ele.numberOfSourceCandidatePtrs(); ipart++){
	if(ele.sourceCandidatePtr(ipart).isNonnull() and ele.sourceCandidatePtr(ipart).isAvailable()){
          particlesInElectron.push_back(ele.sourceCandidatePtr(ipart));
	}
      }

      float p4Error = 0.;
      if(electronPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( electronScaleUnc_.size() == electronSelection_.size()){
	  for(auto selection : electronSelection_){
	    if(selection(ele)){
	      p4Error = electronScaleUnc_.at(iBin);	    
	      break;
	    }
	    else iBin++;
	  }	
	}
      }
      else{
	//uncertainty on the 4 Vector
	if(ele.correctedEcalEnergyError() == 999. and ele.p4Error(reco::GsfElectron::P4_COMBINATION) == 999.){
	  int iBin = 0;
	  if( electronScaleUnc_.size() == electronSelection_.size()){
	    for(auto selection : electronSelection_){
	      if(selection(ele)){
		p4Error = electronScaleUnc_.at(iBin);
		break;
	      }
	      else iBin++;
	    }
	  }
	}
	else{
	  if(ele.correctedEcalEnergyError() != 999.) 
	    p4Error = ele.correctedEcalEnergyError()/ele.energy();	
	  else 
	    p4Error = ele.p4Error(reco::GsfElectron::P4_COMBINATION)/ele.energy();
	}
      }

      // build vector scaled on the transverse plane
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(ele.pt()*(1+p4Error),ele.phi());
      p2Down.SetMagPhi(ele.pt()*(1-p4Error),ele.phi());
      corrElectronEnUp.mex += (-p2Up.Px()+ele.px());
      corrElectronEnUp.mey += (-p2Up.Py()+ele.py());
      corrElectronEnUp.sumet += fabs(p2Up.Mod()-ele.pt());
      corrElectronEnDown.mex += (-p2Down.Px()+ele.px());
      corrElectronEnDown.mey += (-p2Down.Py()+ele.py());
      corrElectronEnDown.sumet += fabs(p2Down.Mod()-ele.pt());

      pat::Electron newElectronEnUp(ele);
      newElectronEnUp.setP4(ele.p4()*(1+p4Error));
      pat::Electron newElectronEnDown(ele);
      newElectronEnDown.setP4(ele.p4()*(1-p4Error));

      electronEnUp->push_back(newElectronEnUp);
      electronEnDown->push_back(newElectronEnDown);
      
    }

    // correct met --> loop on met collection of input
    for(auto met : *metCollection){

      //build corrected MET
      reco::MET corrMETElectronEnUp(met.sumEt()+corrElectronEnUp.sumet,
		       reco::Candidate::LorentzVector(met.px()+corrElectronEnUp.mex, met.py()+corrElectronEnUp.mey, 0., 
		       sqrt((met.px()+corrElectronEnUp.mex)*(met.px()+corrElectronEnUp.mex)+(met.py()+corrElectronEnUp.mey)*(met.py()+corrElectronEnUp.mey))),met.vertex());
      reco::MET corrMETElectronEnDown(met.sumEt()-corrElectronEnDown.sumet,
		reco::Candidate::LorentzVector(met.px()+corrElectronEnDown.mex, met.py()+corrElectronEnDown.mey, 0., 
		sqrt((met.px()+corrElectronEnDown.mex)*(met.px()+corrElectronEnDown.mex)+(met.py()+corrElectronEnDown.mey)*(met.py()+corrElectronEnDown.mey))),met.vertex());

      pat::MET metElectronUp(corrMETElectronEnUp,met);
      pat::MET metElectronDown(corrMETElectronEnDown,met);
      metElectronEnUp->push_back(metElectronUp);
      metElectronEnDown->push_back(metElectronDown);
    }

    if(storeSmearedShiftedCollections_){
      std::sort(electronEnUp->begin(),electronEnUp->end(),electronSorter);
      iEvent.put(electronEnUp,std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      std::sort(electronEnDown->begin(),electronEnDown->end(),electronSorter);
      iEvent.put(electronEnDown,std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metElectronEnUp,std::string(inputMET_.label())+"ElectronEnUp");
    iEvent.put(metElectronEnDown,std::string(inputMET_.label())+"ElectronEnDown");
  }
  
  // loop on photons
  if(not skipPhoton_){
    for(auto pho : *phColl){    
      if(isPuppiMET){
	reco::Candidate::LorentzVector p4 = findParticle(pho,*pfCandCollection);
	if(p4.pt() == 0)
	  continue;
	pho.setP4(p4);
      }
      
      for(size_t ipart = 0 ; ipart < pho.numberOfSourceCandidatePtrs(); ipart++){
	if(pho.sourceCandidatePtr(ipart).isNonnull() and pho.sourceCandidatePtr(ipart).isAvailable()){
          particlesInPhoton.push_back(pho.sourceCandidatePtr(ipart));
	}
      }

      float p4Error = 0.;
      if(photonPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( photonScaleUnc_.size() == photonSelection_.size()){
	  for(auto selection : photonSelection_){
	    if(selection(pho)){
	      p4Error = photonScaleUnc_.at(iBin);	    
	      break;
	    }	
	    else iBin++;
	  }
	}
      }
      else{
	if(pho.getCorrectedEnergyError(reco::Photon::regression2) != 999.)
	  p4Error = pho.getCorrectedEnergyError(reco::Photon::regression2)/pho.energy();
	else{
	  int iBin = 0;
	  if( photonScaleUnc_.size() == photonSelection_.size()){
	    for(auto selection : photonSelection_){
	      if(selection(pho)){
		p4Error = photonScaleUnc_.at(iBin);
		break;
	      }
	      else iBin++;
	    }
	  }	  
	}
      }


      // build vector scaled on the transverse plane
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(pho.pt()*(1+p4Error),pho.phi());
      p2Down.SetMagPhi(pho.pt()*(1-p4Error),pho.phi());
      corrPhotonEnUp.mex += (-p2Up.Px()+pho.px());
      corrPhotonEnUp.mey += (-p2Up.Py()+pho.py());
      corrPhotonEnUp.sumet += fabs(p2Up.Mod()-pho.pt());
      corrPhotonEnDown.mex += (-p2Down.Px()+pho.px());
      corrPhotonEnDown.mey += (-p2Down.Py()+pho.py());
      corrPhotonEnDown.sumet += fabs(p2Down.Mod()-pho.pt());

      pat::Photon newPhotonEnUp(pho);
      newPhotonEnUp.setP4(pho.p4()*(1+p4Error));
      pat::Photon newPhotonEnDown(pho);
      newPhotonEnDown.setP4(pho.p4()*(1-p4Error));

      photonEnUp->push_back(newPhotonEnUp);
      photonEnDown->push_back(newPhotonEnDown);

    }

    // correct met --> loop on met collection of input
    for(auto met : *metCollection){
      //build corrected MET
      reco::MET corrMETPhotonEnUp(met.sumEt()+corrPhotonEnUp.sumet,
		       reco::Candidate::LorentzVector(met.px()+corrPhotonEnUp.mex, met.py()+corrPhotonEnUp.mey, 0., 
		       sqrt((met.px()+corrPhotonEnUp.mex)*(met.px()+corrPhotonEnUp.mex)+(met.py()+corrPhotonEnUp.mey)*(met.py()+corrPhotonEnUp.mey))),met.vertex());
      reco::MET corrMETPhotonEnDown(met.sumEt()-corrPhotonEnDown.sumet,
	       reco::Candidate::LorentzVector(met.px()+corrPhotonEnDown.mex, met.py()+corrPhotonEnDown.mey, 0., 
	       sqrt((met.px()+corrPhotonEnDown.mex)*(met.px()+corrPhotonEnDown.mex)+(met.py()+corrPhotonEnDown.mey)*(met.py()+corrPhotonEnDown.mey))),met.vertex());

      pat::MET metPhotonUp(corrMETPhotonEnUp,met);
      pat::MET metPhotonDown(corrMETPhotonEnDown,met);

      metPhotonEnUp->push_back(metPhotonUp);
      metPhotonEnDown->push_back(metPhotonDown);
    }

    if(storeSmearedShiftedCollections_){
      std::sort(photonEnUp->begin(),photonEnUp->end(),photonSorter);
      iEvent.put(photonEnUp,std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      std::sort(photonEnDown->begin(),photonEnDown->end(),photonSorter);
      iEvent.put(photonEnDown,std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metPhotonEnUp,std::string(inputMET_.label())+"PhotonEnUp");
    iEvent.put(metPhotonEnDown,std::string(inputMET_.label())+"PhotonEnDown");
  }
  
  // taus
  if(not skipTau_){
    for(auto tau : *tauColl){    
      if(isPuppiMET){
	reco::Candidate::LorentzVector p4 = findParticle(tau,*pfCandCollection);
	if(p4.pt() == 0)
	  continue;
	tau.setP4(p4);
      }

      for(size_t ipart = 0 ; ipart < tau.numberOfSourceCandidatePtrs(); ipart++){
	if(tau.sourceCandidatePtr(ipart).isNonnull() and tau.sourceCandidatePtr(ipart).isAvailable()){
          particlesInTau.push_back(tau.sourceCandidatePtr(ipart));
	}
      }

      float p4Error = 0.;
      if(tauPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( tauScaleUnc_.size() == tauSelection_.size()){
	  for(auto selection : tauSelection_){
	    if(selection(tau)){
	      p4Error = tauScaleUnc_.at(iBin);	    
	      break;
	    }	
	    else iBin++;
	  }
	}
      }
      else{
	// do nothing, set to 0
	p4Error = 0;
      }
      
      // build vector scaled on the transverse plane
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(tau.pt()*(1+p4Error),tau.phi());
      p2Down.SetMagPhi(tau.pt()*(1-p4Error),tau.phi());
      corrTauEnUp.mex += (-p2Up.Px()+tau.px());
      corrTauEnUp.mey += (-p2Up.Py()+tau.py());
      corrTauEnUp.sumet += fabs(p2Up.Mod()-tau.pt());
      corrTauEnDown.mex += (-p2Down.Px()+tau.px());
      corrTauEnDown.mey += (-p2Down.Py()+tau.py());
      corrTauEnDown.sumet += fabs(p2Down.Mod()-tau.pt());

      pat::Tau newTauEnUp(tau);
      newTauEnUp.setP4(tau.p4()*(1+p4Error));
      pat::Tau newTauEnDown(tau);
      newTauEnDown.setP4(tau.p4()*(1-p4Error));

      tauEnUp->push_back(newTauEnUp);
      tauEnDown->push_back(newTauEnDown);
 
    }

    // correct met --> loop on met collection of input
    for(auto met : *metCollection){

      //build corrected MET
      reco::MET corrMETTauEnUp(met.sumEt()+corrTauEnUp.sumet,
		       reco::Candidate::LorentzVector(met.px()+corrTauEnUp.mex, met.py()+corrTauEnUp.mey, 0., 
		       sqrt((met.px()+corrTauEnUp.mex)*(met.px()+corrTauEnUp.mex)+(met.py()+corrTauEnUp.mey)*(met.py()+corrTauEnUp.mey))),met.vertex());
      reco::MET corrMETTauEnDown(met.sumEt()-corrTauEnDown.sumet,
		       reco::Candidate::LorentzVector(met.px()+corrTauEnDown.mex, met.py()+corrTauEnDown.mey, 0., 
		       sqrt((met.px()+corrTauEnDown.mex)*(met.px()+corrTauEnDown.mex)+(met.py()+corrTauEnDown.mey)*(met.py()+corrTauEnDown.mey))),met.vertex());
      
      pat::MET metTauUp(corrMETTauEnUp,met);
      pat::MET metTauDown(corrMETTauEnDown,met);
      metTauEnUp->push_back(metTauUp);
      metTauEnDown->push_back(metTauDown);
    }

    if(storeSmearedShiftedCollections_){
      std::sort(tauEnUp->begin(),tauEnUp->end(),tauSorter);
      iEvent.put(tauEnUp,std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      std::sort(tauEnDown->begin(),tauEnDown->end(),tauSorter);
      iEvent.put(tauEnDown,std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metTauEnUp,std::string(inputMET_.label())+"TauEnUp");
    iEvent.put(metTauEnDown,std::string(inputMET_.label())+"TauEnDown");
  }

  if(not skipJet_ ){
    
    std::auto_ptr<JetCorrectionUncertainty> jecUnc;
    JME::JetResolution Resolution;
    JME::JetResolutionScaleFactor ScalarFactor;

    // get the proper file for JEC uncertintay --> jets already calibrated outside this module
    if(jetPSet_.getParameter<bool>("useExternalJECUncertainty") == false){
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get(jetPSet_.getParameter<std::string>("payloadName"),JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      jecUnc = std::auto_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(JetCorPar));
    }
    else{
      if(not jetJECUncFile_.location())
	throw cms::Exception("METSystematicsProducer") << " Failed to find File = " << jetJECUncFile_ << " !!\n";

      jecUnc = std::auto_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jetJECUncFile_.fullPath()));
    
    }

    // take jet energy resolution
    if(jetPSet_.getParameter<bool>("useExternalJER") == false)
      Resolution = JME::JetResolution::get(iSetup,jetPSet_.getParameter<std::string>("payloadName")+"_pt");
    else
      Resolution = JME::JetResolution(jetJERFile_.fullPath());
    
    if(jetPSet_.getParameter<bool>("useExternalJERSF") == false)
      ScalarFactor = JME::JetResolutionScaleFactor::get(iSetup,jetPSet_.getParameter<std::string>("payloadName"));
    else
      ScalarFactor = JME::JetResolutionScaleFactor(jetJERSFFile_.fullPath());

    // loop on jets
    for(auto jet : *jetColl){
      // skip jets
      if(not jetSelection_(jet)) continue;

      // subtract muons that are global or standAlone
      reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4("Uncorrected");

      reco::Candidate::LorentzVector muonP4;
      const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
      for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin(); cand != cands.end(); ++cand ) {
	const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
	const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
	if ( mu != 0 && (mu->isGlobalMuon() || mu->isStandAloneMuon() )) {
	  muonP4 += (*cand)->p4();
	}
      }      
      
      rawJetP4 -= muonP4;

      // re-correct jets
      reco::Candidate::LorentzVector correctedP4;
      PATJetCorrExtractor jetCorrExtractor;
      const pat::Jet & jet_temp = jet;
      if(iEvent.isRealData())
	correctedP4 = jetCorrExtractor(jet_temp, jetPSet_.getParameter<edm::InputTag>("jetCorrLabelRes").label(),9.9, &rawJetP4);      
      else
	correctedP4 = jetCorrExtractor(jet_temp, jetPSet_.getParameter<edm::InputTag>("jetCorrLabel").label(),9.9, &rawJetP4);

      // residual just in case of data
      float residualJES = 0.;
      if(iEvent.isRealData()){
	if (rawJetP4.pt() > 1.e-1 ) {
	  reco::Candidate::LorentzVector corrJetP4upToL3 = jetCorrExtractor(jet_temp, jetPSet_.getParameter<edm::InputTag>("jetCorrLabel").label(),9.9, &rawJetP4);
	  reco::Candidate::LorentzVector corrJetP4upToL3Res = jetCorrExtractor(jet_temp, jetPSet_.getParameter<edm::InputTag>("jetCorrLabelRes").label(),9.9, &rawJetP4);
	  if ( corrJetP4upToL3.pt() > 1.e-1 && corrJetP4upToL3Res.pt() > 1.e-1 ) 
	    residualJES =  (corrJetP4upToL3Res.pt()/corrJetP4upToL3.pt()) - 1.;
	}
      }

      //set eta and pt
      jecUnc->setJetEta(correctedP4.eta());
      jecUnc->setJetPt(correctedP4.pt());
      float uncertainty = jecUnc->getUncertainty(true);
      float p4Error     = sqrt(uncertainty*uncertainty+residualJES*residualJES);
      
           
      // build vector scaled on the transverse plane                                                                                                                            
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(correctedP4.pt()*(1+p4Error),correctedP4.phi());
      p2Down.SetMagPhi(correctedP4.pt()*(1-p4Error),correctedP4.phi());
      corrJetEnUp.mex += (-p2Up.Px()+correctedP4.px());
      corrJetEnUp.mey += (-p2Up.Py()+correctedP4.py());
      corrJetEnUp.sumet += fabs(p2Up.Mod()-correctedP4.pt());
      corrJetEnDown.mex += (-p2Down.Px()+correctedP4.px());
      corrJetEnDown.mey += (-p2Down.Py()+correctedP4.py());
      corrJetEnDown.sumet += fabs(p2Down.Mod()-correctedP4.pt());
     
      pat::Jet newJetEnUp(jet);
      newJetEnUp.setP4(jet.p4()*(1+p4Error));
      pat::Jet newJetEnDown(jet);
      newJetEnDown.setP4(jet.p4()*(1-p4Error));

      jetEnUp->push_back(newJetEnUp);
      jetEnDown->push_back(newJetEnDown);

      //take scale factor value for JER
      float jerSF   = 1.;
      float jerSFUp = 1.;
      float jerSFDw = 1.;

      // for jet energy resolution
      JME::JetParameters jetParam;
      jetParam.setJetPt(correctedP4.pt());
      jetParam.setJetEta(correctedP4.eta());
      jetParam.setRho(rho);

      jerSF   = ScalarFactor.getScaleFactor(jetParam, Variation::NOMINAL);
      jerSFUp = ScalarFactor.getScaleFactor(jetParam, Variation::UP);
      jerSFDw = ScalarFactor.getScaleFactor(jetParam, Variation::DOWN);
      
      //take the truth resolution
      float resCorrection     = 0;
      float resCorrectionUp   = 0;
      float resCorrectionDown = 0;

      bool useRandomSmear = true;
      float sigmaPt = Resolution.getResolution(jetParam);

      if(jet.genJet() !=0){ // gen matched jet found
	if(reco::deltaR(correctedP4,jet.genJet()->p4()) < 0.2 and fabs((correctedP4.pt()-jet.genJet()->pt()))/correctedP4.pt() < 3*sigmaPt){ // mathing withing dR = 0.2 and deviates less than 3 sigma in E
	  resCorrection     = 1+(jerSF-1)*(correctedP4.pt()-jet.genJet()->pt())/correctedP4.pt(); // 1+correction*
	  resCorrectionUp   = 1+(jerSFUp-1)*(correctedP4.pt()-jet.genJet()->pt())/correctedP4.pt();
	  resCorrectionDown = 1+(jerSFDw-1)*(correctedP4.pt()-jet.genJet()->pt())/correctedP4.pt();
	  useRandomSmear    = false;
	}
      }
	
      if(useRandomSmear){
	    
	// take jet resolution from db/text file in case of random smearing
	// apply smearing to jet 4V
	float smearingWidth     = sqrt(jerSF*jerSF-1)*sigmaPt;
	float smearingWidthUp   = sqrt(jerSFUp*jerSFUp-1)*sigmaPt;
	float smearingWidthDown = sqrt(jerSFDw*jerSFDw-1)*sigmaPt;
	resCorrection     = 1+rand_->Gaus(0,smearingWidth)/correctedP4.pt();       
	resCorrectionUp   = 1+rand_->Gaus(0,smearingWidthUp)/correctedP4.pt(); 
        resCorrectionDown = 1+rand_->Gaus(0,smearingWidthDown)/correctedP4.pt();
      }
      
      pat::Jet newJetSmear(jet);
      pat::Jet newJetSmearResUp(jet);
      pat::Jet newJetSmearResDown(jet);
      newJetSmear.setP4(jet.p4()*resCorrection);
      newJetSmearResUp.setP4(jet.p4()*resCorrectionUp);
      newJetSmearResDown.setP4(jet.p4()*resCorrectionDown);           

      jetSmear->push_back(newJetSmear);
      jetSmearResUp->push_back(newJetSmearResUp);
      jetSmearResDown->push_back(newJetSmearResDown);
      

      // build vector scaled on the transverse plane                                                                                                                            
      TVector2 p2Smear;
      pat::Jet jetResMET   = jet;
      pat::Jet jetResUpMET = jet;
      pat::Jet jetResDwMET = jet;
      jetResMET.setP4(correctedP4*resCorrection);
      jetResUpMET.setP4(correctedP4*resCorrectionUp);
      jetResDwMET.setP4(correctedP4*resCorrectionDown);

      // re-correct jets after smearing, only the central value for smeared met                                                                                                                  
      const pat::Jet & jet_temp_2 = jetResMET;
      reco::Candidate::LorentzVector rawP4Smear = jet_temp_2.correctedP4("Uncorrected");
      reco::Candidate::LorentzVector correctedP4Smear;
      if(iEvent.isRealData())
        correctedP4Smear = jetCorrExtractor(jet_temp_2, jetPSet_.getParameter<edm::InputTag>("jetCorrLabelRes").label(),9.9, &rawP4Smear);
      else
        correctedP4Smear = jetCorrExtractor(jet_temp_2, jetPSet_.getParameter<edm::InputTag>("jetCorrLabel").label(),9.9, &rawP4Smear);

      jetResMET.setP4(correctedP4Smear);
      

      // clone correction level from imput jet      
      p2Smear.SetMagPhi(jetResMET.pt(),jetResMET.phi());
      p2Up.SetMagPhi(jetResUpMET.pt(),jetResUpMET.phi());
      p2Down.SetMagPhi(jetResDwMET.pt(),jetResDwMET.phi());

      corrJetResUp.mex += (-p2Up.Px()+correctedP4.px());
      corrJetResUp.mey += (-p2Up.Py()+correctedP4.py());
      corrJetResUp.sumet += p2Up.Mod()-correctedP4.pt();

      corrJetResDown.mex += (-p2Down.Px()+correctedP4.px());
      corrJetResDown.mey += (-p2Down.Py()+correctedP4.py());
      corrJetResDown.sumet += p2Down.Mod()-correctedP4.pt();

      if(jetSelection_(jetResMET)){
	corrJetSmear.mex += (-p2Smear.Px()+correctedP4.px());
	corrJetSmear.mey += (-p2Smear.Py()+correctedP4.py());
	corrJetSmear.sumet += p2Smear.Mod()-correctedP4.pt();
      }
      if(jetSelection_(jetResUpMET)){
	corrJetSmearResUp.mex += (-p2Up.Px()+jetResMET.px());
	corrJetSmearResUp.mey += (-p2Up.Py()+jetResMET.py());
	corrJetSmearResUp.sumet += p2Up.Mod()-jetResMET.pt();
      }
      if(jetSelection_(jetResDwMET)){
	corrJetSmearResDown.mex += (-p2Down.Px()+jetResMET.px());
	corrJetSmearResDown.mey += (-p2Down.Py()+jetResMET.py());
	corrJetSmearResDown.sumet += p2Down.Mod()-jetResMET.pt();
      }
      
      //store all the candidates inside the jet
      particlesInJet.insert(particlesInJet.end(), jet.daughterPtrVector().begin(), jet.daughterPtrVector().end());
    }

    //sort Ptr according to the object key
    std::sort(particlesInJet.begin(),particlesInJet.end());

  
    // check all the particles not belonging to a jet over threshold
    for(unsigned int icand = 0; icand < pfCandCollection->size(); icand++){

      reco::CandidatePtr ptrCand = pfCandCollection->ptrAt(icand);
      float weightp4 =  1.;

      if(isPuppiMET){
	const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&(pfCandCollection->at(icand)));
	weightp4 = lPack->puppiWeight();
	if(lPack->puppiWeight() == 0)
	  continue;
      }

      // check overlaps with other particles
      bool isfound = false;
      for(auto part : particlesInJet){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}      	
      }

      if(isfound) continue;
       
      for(auto part : particlesInTau){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}      	
      }

      if(isfound) continue;

      for(auto part : particlesInMuon){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}      	
      }

      if(isfound) continue;

      for(auto part : particlesInElectron){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}      	
      }

      if(isfound) continue;

      for(auto part : particlesInPhoton){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}      	
      }

      if(isfound) continue;


      float uncertaintyP4 = 0;
      int iBin = 0;
      for(auto selection : unclusteredSelection_){ // use the formula for single particle PF resolution	   
	if(selection(pfCandCollection->at(icand))){
	  TString formula (unclusteredUnc_.at(iBin).c_str());
	  double  varX = 0;
	  double  varY = 0;
	  if(formula.Contains("pt")){
	    formula.ReplaceAll("pt","x");
	    varX = pfCandCollection->at(icand).pt();
	  }
	  else if(formula.Contains("energy")){
	    formula.ReplaceAll("energy","x");
	    varX = pfCandCollection->at(icand).energy();
	  }
	  if(formula.Contains("eta")){
	    formula.ReplaceAll("eta","y");
	    varY = pfCandCollection->at(icand).eta();
	  }
	  TFormula form ("form",formula.Data());
	  uncertaintyP4 = form.Eval(varX,varY);
	  iBin++;
	  break;
	}
	else
	  iBin++;
      }
      
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(pfCandCollection->at(icand).pt()*weightp4*(1+uncertaintyP4),pfCandCollection->at(icand).phi());
      p2Down.SetMagPhi(pfCandCollection->at(icand).pt()*weightp4*(1-uncertaintyP4),pfCandCollection->at(icand).phi());
      corrUnclusteredEnUp.mex += (-p2Up.Px()+pfCandCollection->at(icand).px());
      corrUnclusteredEnUp.mey += (-p2Up.Py()+pfCandCollection->at(icand).py());
      corrUnclusteredEnUp.sumet += fabs(p2Up.Mod()-pfCandCollection->at(icand).pt());
      corrUnclusteredEnDown.mex += (-p2Down.Px()+pfCandCollection->at(icand).px());
      corrUnclusteredEnDown.mey += (-p2Down.Py()+pfCandCollection->at(icand).py());
      corrUnclusteredEnDown.sumet += fabs(p2Down.Mod()-pfCandCollection->at(icand).pt());	  
    }

    // correct met --> loop on met collection of input                                                                                                                          
    for(auto met : *metCollection){

      //build corrected MET                                                                                                                                                     
      reco::MET corrMETJetEnUp(met.sumEt()+corrJetEnUp.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetEnUp.mex, met.py()+corrJetEnUp.mey, 0.,
		sqrt((met.px()+corrJetEnUp.mex)*(met.px()+corrJetEnUp.mex)+(met.py()+corrJetEnUp.mey)*(met.py()+corrJetEnUp.mey))),met.vertex());

      reco::MET corrMETJetEnDown(met.sumEt()-corrJetEnDown.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetEnDown.mex, met.py()+corrJetEnDown.mey, 0.,
		sqrt((met.px()+corrJetEnDown.mex)*(met.px()+corrJetEnDown.mex)+(met.py()+corrJetEnDown.mey)*(met.py()+corrJetEnDown.mey))),met.vertex());

      pat::MET metJetUp(corrMETJetEnUp,met);
      pat::MET metJetDown(corrMETJetEnDown,met);
      metJetEnUp->push_back(metJetUp);
      metJetEnDown->push_back(metJetDown);

      reco::MET corrMETJetResUp(met.sumEt()+corrJetResUp.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetResUp.mex, met.py()+corrJetResUp.mey, 0.,
		sqrt((met.px()+corrJetResUp.mex)*(met.px()+corrJetResUp.mex)+(met.py()+corrJetResUp.mey)*(met.py()+corrJetResUp.mey))),met.vertex());

      reco::MET corrMETJetResDown(met.sumEt()-corrJetResDown.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetResDown.mex, met.py()+corrJetResDown.mey, 0.,
		sqrt((met.px()+corrJetResDown.mex)*(met.px()+corrJetResDown.mex)+(met.py()+corrJetResDown.mey)*(met.py()+corrJetResDown.mey))),met.vertex());

      pat::MET metJetRUp(corrMETJetResUp,met);
      pat::MET metJetRDown(corrMETJetResDown,met);
      metJetResUp->push_back(metJetRUp);
      metJetResDown->push_back(metJetRDown);

      
      reco::MET corrMETJetSmear(met.sumEt()+corrJetSmear.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmear.mex, met.py()+corrJetSmear.mey, 0.,
		sqrt((met.px()+corrJetSmear.mex)*(met.px()+corrJetSmear.mex)+(met.py()+corrJetSmear.mey)*(met.py()+corrJetSmear.mey))),met.vertex());

      pat::MET metJetSme(corrMETJetSmear,met);
      metJetSmear->push_back(metJetSme);

      reco::MET corrMETJetSmearResUp(met.sumEt()+corrJetSmearResUp.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmearResUp.mex, met.py()+corrJetSmearResUp.mey, 0.,
					       sqrt((met.px()+corrJetSmearResUp.mex)*(met.px()+corrJetSmearResUp.mex)+(met.py()+corrJetSmearResUp.mey)*(met.py()+corrJetSmearResUp.mey))),met.vertex());

      reco::MET corrMETJetSmearResDown(met.sumEt()+corrJetSmearResDown.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmearResDown.mex, met.py()+corrJetSmearResDown.mey, 0.,
		sqrt((met.px()+corrJetSmearResDown.mex)*(met.px()+corrJetSmearResDown.mex)+(met.py()+corrJetSmearResDown.mey)*(met.py()+corrJetSmearResDown.mey))),met.vertex());

      pat::MET metJetSmeUp(corrMETJetSmearResUp,met);
      pat::MET metJetSmeDown(corrMETJetSmearResDown,met);


      metJetSmearResUp->push_back(metJetSmeUp);
      metJetSmearResDown->push_back(metJetSmeDown);

      // unclustered
      reco::MET corrMETUnclusteredEnUp(met.sumEt()+corrUnclusteredEnUp.sumet,
	 reco::Candidate::LorentzVector(met.px()+corrUnclusteredEnUp.mex, met.py()+corrUnclusteredEnUp.mey, 0.,
	 sqrt((met.px()+corrUnclusteredEnUp.mex)*(met.px()+corrUnclusteredEnUp.mex)+(met.py()+corrUnclusteredEnUp.mey)*(met.py()+corrUnclusteredEnUp.mey))),met.vertex());

      reco::MET corrMETUnclusteredEnDown(met.sumEt()-corrUnclusteredEnDown.sumet,
       reco::Candidate::LorentzVector(met.px()+corrUnclusteredEnDown.mex, met.py()+corrUnclusteredEnDown.mey, 0.,
       sqrt((met.px()+corrUnclusteredEnDown.mex)*(met.px()+corrUnclusteredEnDown.mex)+(met.py()+corrUnclusteredEnDown.mey)*(met.py()+corrUnclusteredEnDown.mey))),met.vertex());

      pat::MET metUncEnUp(corrMETUnclusteredEnUp,met);
      pat::MET metUncEnDown(corrMETUnclusteredEnDown,met);
      metUnclusteredEnUp->push_back(metUncEnUp);
      metUnclusteredEnDown->push_back(metUncEnDown);

    }

    if(storeSmearedShiftedCollections_){
      std::sort(jetEnUp->begin(),jetEnUp->end(),jetSorter); // sort again collections
      iEvent.put(jetEnUp,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      std::sort(jetEnDown->begin(),jetEnDown->end(),jetSorter); // sort again collections
      iEvent.put(jetEnDown,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metJetEnUp,std::string(inputMET_.label())+"JetEnUp");
    iEvent.put(metJetEnDown,std::string(inputMET_.label())+"JetEnDown");    
   
    if(storeSmearedShiftedCollections_){
      std::sort(jetSmear->begin(),jetSmear->end(),jetSorter);
      iEvent.put(jetSmear,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"Smear");
      std::sort(jetSmearResUp->begin(),jetSmearResUp->end(),jetSorter);
      iEvent.put(jetSmearResUp,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResUp");
      std::sort(jetSmearResDown->begin(),jetSmearResDown->end(),jetSorter);
      iEvent.put(jetSmearResDown,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResDown");
    }
    iEvent.put(metJetResUp,std::string(inputMET_.label())+"JetResUp");
    iEvent.put(metJetResDown,std::string(inputMET_.label())+"JetResDown");    
    iEvent.put(metJetSmear,std::string(inputMET_.label())+"Smear");
    iEvent.put(metJetSmearResUp,std::string(inputMET_.label())+"SmearJetResUp");
    iEvent.put(metJetSmearResDown,std::string(inputMET_.label())+"SmearJetResDown");    

    iEvent.put(metUnclusteredEnUp,std::string(inputMET_.label())+"UnclusteredEnUp");
    iEvent.put(metUnclusteredEnDown,std::string(inputMET_.label())+"UnclusteredEnDown");    
  }
}


template<typename T>
reco::Candidate::LorentzVector METSystematicsProducer::findParticle(const T & particle, const edm::View<reco::Candidate> & pfCandCollection){


  reco::Candidate::LorentzVector total4V;
  std::vector<reco::CandidatePtr> particles;

  // take all the PF candidate linked to the particles
  for(size_t ipart = 0 ; ipart < particle.numberOfSourceCandidatePtrs(); ipart++){      
    if(particle.sourceCandidatePtr(ipart).isNonnull() and particle.sourceCandidatePtr(ipart).isAvailable()){
      particles.push_back(particle.sourceCandidatePtr(ipart));
    }
  }

  // no-pfCandidates return 0
  if(particles.size() == 0)
    return particle.p4();  

  for(unsigned int icand = 0; icand <  pfCandCollection.size(); icand++){
    
    reco::CandidatePtr ptrCand = pfCandCollection.ptrAt(icand);
    const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>(&(pfCandCollection.at(icand)));	
    // check the puppiWeight used in the met calculation
    if(lPack->puppiWeightNoLep() == 0) // in case of puppi, only the one used for MET calculations are useful
      continue;
    
    // find and return the resclae 4V
    for(auto ipart : particles){
      if(ipart == ptrCand){
	total4V += ptrCand->p4()*lPack->puppiWeightNoLep();
	break;
      }      
    }
  }
     
  return total4V;

}


void METSystematicsProducer::beginJob() {
}

void METSystematicsProducer::endJob() {
}

void METSystematicsProducer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void METSystematicsProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void METSystematicsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void METSystematicsProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void METSystematicsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(METSystematicsProducer);
