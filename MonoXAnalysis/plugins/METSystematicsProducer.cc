#ifndef METSystematicsProducer_h
#define METSystematicsProducer_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/METReco/interface/CorrMETData.h"
#include "DataFormats/METReco/interface/MET.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TFormula.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

class JetResolutionTruth {

public:
  JetResolutionTruth(const float & etaMin, const float & etaMax,
		     const float & rhoMin, const float & rhoMax,
		     const std::vector<float> & par){

    etaMin_ = etaMin;
    etaMax_ = etaMax;
    rhoMin_ = rhoMin;
    rhoMax_ = rhoMax;
    par_ = par;
  }

  bool goodJet (const pat::Jet & jet, const float & rho){
    if(jet.eta() > etaMin_ and jet.eta() < etaMax_){
	if(rho > rhoMin_ and rho < rhoMax_)
	  return true;	
    }    
    return false;
  }

  float getEtaMin(){return etaMin_;}

  float getEtaMax(){return etaMax_;}

  float getRhoMin(){return rhoMin_;}

  float getRhoMax(){return rhoMax_;}

  float getPar(const int & nPar){
    if(nPar <= int(par_.size()))
      return par_.at(nPar);
    else
      return 0.;
  }

  std::vector<float> getPar(){return par_;}

  ~JetResolutionTruth(){};

private:

  float etaMin_;
  float etaMax_;
  float rhoMin_;
  float rhoMax_;
  std::vector<float> par_;
};


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

  const edm::EDGetTokenT<pat::METCollection> metToken_;
  const edm::EDGetTokenT<edm::View<reco::Candidate> > pfCandToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> >      jetToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> >     muonToken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronToken_;
  edm::EDGetTokenT<std::vector<pat::Photon> >   photonToken_;
  edm::EDGetTokenT<std::vector<pat::Tau> >      tauToken_;
  const edm::EDGetTokenT<double>  rhoToken_;

  const edm::ParameterSet electronPSet_;
  const edm::ParameterSet muonPSet_;
  const edm::ParameterSet photonPSet_;
  const edm::ParameterSet tauPSet_;
  const edm::ParameterSet jetPSet_;
  const edm::ParameterSet unclusteredPSet_;


  const edm::InputTag inputMET_;

  bool skipJet_;
  bool skipMuon_;
  bool skipElectron_;
  bool skipTau_;
  bool skipPhoton_;
  bool storeSmearedShiftedCollections_;

  // in case of binning
  std::vector<StringCutObjectSelector<pat::Muon> > muonSelection_;
  std::vector<float> muonScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Electron> > electronSelection_;
  std::vector<float> electronScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Photon> > photonSelection_;
  std::vector<float> photonScaleUnc_;
  std::vector<StringCutObjectSelector<pat::Tau> > tauSelection_;
  std::vector<float> tauScaleUnc_;

  std::string jetJECUncFile_;
  std::string jetJERFile_;
  std::auto_ptr<TFormula>   jetJERFormula_;

  std::vector<StringCutObjectSelector<pat::Jet> > jetJERSelection_;
  std::vector<float> jetJERSF_;
  std::vector<float> jetJERSFUnc_;
  std::vector<JetResolutionTruth> jetResolution_;
  std::auto_ptr<TRandom3> rand_;

  std::vector<StringCutObjectSelector<reco::Candidate> > unclusteredSelection_;
  std::vector<float> unclusteredUnc_;

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
  storeSmearedShiftedCollections_(iConfig.existsAs<bool>("storeSmearedShiftedCollections") ? iConfig.getParameter<bool>("storeSmearedShiftedCollections") : false){

  skipJet_  = false;
  skipTau_  = false;
  skipPhoton_ = false;
  skipMuon_ = false;
  skipElectron_ = false;
   
  if(!(jetPSet_.getParameter<edm::InputTag>("src") == edm::InputTag("")))
    jetToken_     = consumes<std::vector<pat::Jet> >       (jetPSet_.getParameter<edm::InputTag>("src"));
  else{
    skipJet_ = true;
    std::cout<<"METSystematicsProducer::Jet collection will be skipped"<<std::endl;
  }

  if(!(tauPSet_.getParameter<edm::InputTag>("src") == edm::InputTag("")))
    tauToken_      = consumes<std::vector<pat::Tau> >       (tauPSet_.getParameter<edm::InputTag>("src"));
  else{
    skipTau_ = true;
    std::cout<<"METSystematicsProducer::Tau collection will be skipped"<<std::endl;
  }

  if(!(photonPSet_.getParameter<edm::InputTag>("src") == edm::InputTag("")))
    photonToken_   = consumes<std::vector<pat::Photon> >    (photonPSet_.getParameter<edm::InputTag>("src"));
  else{
    skipPhoton_ = true;
    std::cout<<"METSystematicsProducer::Photon collection will be skipped"<<std::endl;
  }
  
  if(!(muonPSet_.getParameter<edm::InputTag>("src") == edm::InputTag("")))
    muonToken_   = consumes<std::vector<pat::Muon> >      (muonPSet_.getParameter<edm::InputTag>("src"));
  else{
    skipMuon_ = true;
    std::cout<<"METSystematicsProducer::Muon collection will be skipped"<<std::endl;
  }

  if(!(electronPSet_.getParameter<edm::InputTag>("src") == edm::InputTag("")))
    electronToken_ = consumes<std::vector<pat::Electron> >  (electronPSet_.getParameter<edm::InputTag>("src"));
  else{
    skipElectron_ = true;
    std::cout<<"METSystematicsProducer::Electron collection will be skipped"<<std::endl;
  }


  // produce shifed and smearead particles
  produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label()));

  if (skipMuon_ == false){
    if(storeSmearedShiftedCollections_ == false){
      produces<pat::MuonCollection>(std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::MuonCollection>(std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"MuonEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"MuonEnDown");
  }

  if(not skipElectron_){
    if(storeSmearedShiftedCollections_ == false){
      produces<pat::ElectronCollection>(std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::ElectronCollection>(std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"ElectronEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"ElectronEnDown");
  }
  
  if(not skipPhoton_){
    if(storeSmearedShiftedCollections_ == false){
      produces<pat::PhotonCollection>(std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::PhotonCollection>(std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"PhotonEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"PhotonEnDown");
  }

  if(not skipTau_){
    if(storeSmearedShiftedCollections_ == false){
      produces<pat::TauCollection>(std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::TauCollection>(std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"TauEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"TauEnDown");
  }

  if(not skipJet_){
    if(storeSmearedShiftedCollections_ == false){
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");    
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"Smear");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResUp");
      produces<pat::JetCollection>(std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResDown");
    }
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"JetEnDown");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"Smear");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"SmearJetResUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"SmearJetResDown");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"UnclusteredEnUp");
    produces<pat::METCollection>(std::string(iConfig.getParameter<edm::InputTag>("inputMET").label())+"UnclusteredEnDown");
  }

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
    unclusteredUnc_.push_back(bin.getParameter<double>("uncertainty"));
  }

  // file for JEC unc
  jetJECUncFile_ = jetPSet_.getParameter<std::string>("JECUncFile");

  // JER SF
  for(auto bin : jetPSet_.getParameter<std::vector<edm::ParameterSet> >("binningJERSF")){
    jetJERSelection_.push_back(StringCutObjectSelector<pat::Jet>(bin.getParameter<std::string>("binSelection")));
    jetJERSF_.push_back(bin.getParameter<double>("scaleFactor"));
    jetJERSFUnc_.push_back(bin.getParameter<double>("scaleFactorUnc"));
  }
  
  // file for JER truth
  jetJERFile_    = jetPSet_.getParameter<std::string>("JERFile");
  jetJERFormula_ = std::auto_ptr<TFormula>(new TFormula("jerFormula",jetPSet_.getParameter<std::string>("JERFormula").c_str()));

  // this setup really depends on the format coded in https://twiki.cern.ch/twiki/pub/CMS/JetResolution/Summer15_25nsV6_MC_PtResolution_AK4PFchs.txt
  std::string line;
  std::ifstream file(jetJERFile_);
  if(file.is_open()){
    while(getline(file,line)){
      if(line.at(0)=='{')
	continue;
      std::vector<std::string> tokens;
      std::stringstream ss(line);
      std::string buf;
      while (ss >> buf)
        tokens.push_back(buf);
      std::vector<float> param;
      for(size_t itok = tokens.size()-4; itok < tokens.size(); itok++){
	param.push_back(std::stof(tokens[itok]));
      }
      jetResolution_.push_back(JetResolutionTruth(std::stof(tokens[0]),std::stof(tokens[1]),std::stof(tokens[2]),std::stof(tokens[3]), param));
    }
    file.close();
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
  CorrMETData corrJetSmear;
  CorrMETData corrJetSmearResUp;
  CorrMETData corrJetSmearResDown;
  std::auto_ptr<pat::METCollection> metJetEnUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetEnDown(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmearResUp(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmearResDown(new pat::METCollection);
  std::auto_ptr<pat::METCollection> metJetSmear(new pat::METCollection);

  CorrMETData corrUnclusteredEnUp;
  CorrMETData corrUnclusteredEnDown;
  std::auto_ptr<pat::METCollection> metUnclusteredEnUp (new pat::METCollection);
  std::auto_ptr<pat::METCollection> metUnclusteredEnDown (new pat::METCollection);

  bool isPuppiMET = false;
  if(TString(inputMET_.label()).Contains("Puppi") or TString(inputMET_.label()).Contains("puppi") or TString(inputMET_.label()).Contains("PUPPI")) 
    isPuppiMET = true;

  // loop on muons
  if(not skipMuon_){
    for(auto muon : *muColl){    

      if(isPuppiMET){
	reco::Candidate::LorentzVector p4 = findParticle(muon,*pfCandCollection);
	if(p4.pt() == 0)
	  continue;
	muon.setP4(p4);
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
	    else
	      iBin++;
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

    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(muonEnUp,std::string(muonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
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

      float p4Error = 0.;
      if(electronPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( electronScaleUnc_.size() == electronSelection_.size()){
	  for(auto selection : electronSelection_){
	    if(selection(ele)){
	      p4Error = electronScaleUnc_.at(iBin);	    
	      break;
	    }
	    else
	      iBin++;
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
	      else
		iBin++;
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

    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(electronEnUp,std::string(electronPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
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

      float p4Error = 0.;
      if(photonPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( photonScaleUnc_.size() == photonSelection_.size()){
	  for(auto selection : photonSelection_){
	    if(selection(pho)){
	      p4Error = photonScaleUnc_.at(iBin);	    
	      break;
	    }	
	    else
	      iBin++;
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
	      else
		iBin++;
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

    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(photonEnUp,std::string(photonPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
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

      float p4Error = 0.;
      if(tauPSet_.getParameter<bool>("useExternalUncertainty")){	
	int iBin = 0;
	if( tauScaleUnc_.size() == tauSelection_.size()){
	  for(auto selection : tauSelection_){
	    if(selection(tau)){
	      p4Error = tauScaleUnc_.at(iBin);	    
	      break;
	    }	
	    else
	      iBin++;
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

    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(tauEnUp,std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      iEvent.put(tauEnDown,std::string(tauPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metTauEnUp,std::string(inputMET_.label())+"TauEnUp");
    iEvent.put(metTauEnDown,std::string(inputMET_.label())+"TauEnDown");
  }

  if(not skipJet_ ){
    
    std::auto_ptr<JetCorrectionUncertainty> jecUnc;
    std::vector<reco::CandidatePtr> particlesInJet;

    if(jetPSet_.getParameter<bool>("useExternalJECUncertainty") == false){
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      iSetup.get<JetCorrectionsRecord>().get(jetPSet_.getParameter<std::string>("payloadName"),JetCorParColl); 
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      jecUnc = std::auto_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(JetCorPar));
    }
    else
      jecUnc = std::auto_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jetJECUncFile_));
    

    for(auto jet : *jetColl){
      //set eta and pt
      jecUnc->setJetEta(jet.eta());
      jecUnc->setJetPt(jet.pt()); 
      float p4Error = jecUnc->getUncertainty(true);

      // build vector scaled on the transverse plane                                                                                                                            
      TVector2 p2Up, p2Down;
      p2Up.SetMagPhi(jet.pt()*(1+p4Error),jet.phi());
      p2Down.SetMagPhi(jet.pt()*(1-p4Error),jet.phi());
      corrJetEnUp.mex += (-p2Up.Px()+jet.px());
      corrJetEnUp.mey += (-p2Up.Py()+jet.py());
      corrJetEnUp.sumet += fabs(p2Up.Mod()-jet.pt());
      corrJetEnDown.mex += (-p2Down.Px()+jet.px());
      corrJetEnDown.mey += (-p2Down.Py()+jet.py());
      corrJetEnDown.sumet += fabs(p2Down.Mod()-jet.pt());

      pat::Jet newJetEnUp(jet);
      newJetEnUp.setP4(jet.p4()*(1+p4Error));
      pat::Jet newJetEnDown(jet);
      newJetEnDown.setP4(jet.p4()*(1-p4Error));

      jetEnUp->push_back(newJetEnUp);
      jetEnDown->push_back(newJetEnDown);
      
      //take scale factor value for JER
      float jerSF = 1.;
      float jerSFUnc = 0.;
      int iBin = 0;
      if(jetJERSF_.size() == jetJERSelection_.size() and jetJERSFUnc_.size() == jetJERSelection_.size()){
	for(auto selection : jetJERSelection_){
	  if(selection(jet)){
	    jerSF = jetJERSF_.at(iBin);	    
	    jerSFUnc = jetJERSFUnc_.at(iBin);
	    break;
	  }
	  else
	    iBin++;	
	}
      }

      //take the truth resolution
      float resCorrection = 1 ;
      float resCorrectionUp = 1 ;
      float resCorrectionDown = 1 ;

      bool useRandomSmear = true;
      if(jet.genJet() !=0){ // gen matched jet found
	if(reco::deltaR(jet.p4(),jet.genJet()->p4()) < std::min(0.4, 0.1 + 0.3*exp(-0.05*(jet.genJet()->pt()-10.)))){
	  resCorrection = std::max(0.,(jerSF-1)*(jet.pt()-jet.genJet()->pt())/jet.pt()); // 1+correction*
	  resCorrectionUp = std::max(0.,(jerSF+jerSFUnc-1)*(jet.pt()-jet.genJet()->pt())/jet.pt());
	  resCorrectionDown = std::max(0.,(jerSF-jerSFUnc-1)*(jet.pt()-jet.genJet()->pt())/jet.pt());
	  useRandomSmear = false;
	}
      }
	
      if(useRandomSmear){
	float sigmaPt = 0;
	for(auto jetRes : jetResolution_){	
	    if(jetRes.goodJet(jet,rho)){
	      std::vector<float> par = jetRes.getPar();
	      for(size_t ipar = 0; ipar < par.size(); ipar++){
		jetJERFormula_->SetParameter(ipar,par[ipar]);
	      }
	    }
	}
	    
	sigmaPt = jetJERFormula_->Eval(jet.pt());
	
	// apply smearing to jet 4V
	float smearingWidth     = sqrt(jerSF*jerSF-1)*sigmaPt;
	float smearingWidthUp   = sqrt((jerSF+jerSFUnc)*(jerSF+jerSFUnc)-1)*sigmaPt;
	float smearingWidthDown = sqrt((jerSF-jerSFUnc)*(jerSF-jerSFUnc)-1)*sigmaPt;
	do{ resCorrection   = rand_->Gaus(0,smearingWidth); } while(resCorrection <= -1.); 
	do{ resCorrectionUp = rand_->Gaus(0,smearingWidthUp); }while(resCorrection <= -1.);
	do{ resCorrectionDown = rand_->Gaus(0,smearingWidthDown); }while(resCorrection <= -1.);
      }
      
      float smearedPt = std::max(1.e-2,jet.pt()*(1+resCorrection));
      float smearedPtUp = std::max(1.e-2,jet.pt()*(1+resCorrectionUp));
      float smearedPtDown = std::max(1.e-2,jet.pt()*(1+resCorrectionDown));
	

      // build vector scaled on the transverse plane                                                                                                                            
      TVector2 p2Smear;
      p2Smear.SetMagPhi(smearedPt,jet.phi());
      p2Up.SetMagPhi(smearedPtUp,jet.phi());
      p2Down.SetMagPhi(smearedPtDown,jet.phi());
      corrJetSmear.mex += (-p2Smear.Px()+jet.px());
      corrJetSmear.mey += (-p2Smear.Py()+jet.py());
      corrJetSmear.sumet += p2Smear.Mod()-jet.pt();
      corrJetSmearResUp.mex += (-p2Up.Px()+jet.px());
      corrJetSmearResUp.mey += (-p2Up.Py()+jet.py());
      corrJetSmearResUp.sumet += p2Up.Mod()-jet.pt();
      corrJetSmearResDown.mex += (-p2Down.Px()+jet.px());
      corrJetSmearResDown.mey += (-p2Down.Py()+jet.py());
      corrJetSmearResDown.sumet += p2Down.Mod()-jet.pt();

      pat::Jet newJetSmear(jet);
      newJetSmear.setP4(jet.p4()*(1+resCorrection));
      pat::Jet newJetSmearResUp(jet);
      newJetSmearResUp.setP4(jet.p4()*(1+resCorrectionUp));
      pat::Jet newJetSmearResDown(jet);
      newJetSmearResDown.setP4(jet.p4()*(1+resCorrectionDown));           

      if(newJetSmear.pt() > 1.e-2) 
	jetSmear->push_back(newJetSmear);
      if(newJetSmearResUp.pt() > 1.e-2)     
	jetSmearResUp->push_back(newJetSmearResUp);
      if(newJetSmearResDown.pt() > 1.e-2)
	jetSmearResDown->push_back(newJetSmearResDown);
      
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

      bool isfound = false;
      for(auto part : particlesInJet){  	
	if(part == ptrCand){
	  isfound = true;
	  break;
	}
      }

      if(isfound)
	continue;
       
      float uncertaintyP4 = 0;
      int iBin = 0;
      for(auto selection : unclusteredSelection_){	  
	if(selection(pfCandCollection->at(icand))){
	  uncertaintyP4 = unclusteredUnc_.at(iBin);
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
      
      reco::MET corrMETJetSmear(met.sumEt()+corrJetSmear.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmear.mex, met.py()+corrJetSmear.mey, 0.,
		sqrt((met.px()+corrJetSmear.mex)*(met.px()+corrJetSmear.mex)+(met.py()+corrJetSmear.mey)*(met.py()+corrJetSmear.mey))),met.vertex());

      reco::MET corrMETJetSmearResUp(met.sumEt()+corrJetSmearResUp.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmearResUp.mex, met.py()+corrJetSmearResUp.mey, 0.,
		sqrt((met.px()+corrJetSmearResUp.mex)*(met.px()+corrJetSmearResUp.mex)+(met.py()+corrJetSmearResUp.mey)*(met.py()+corrJetSmearResUp.mey))),met.vertex());

      reco::MET corrMETJetSmearResDown(met.sumEt()+corrJetSmearResDown.sumet,
		reco::Candidate::LorentzVector(met.px()+corrJetSmearResDown.mex, met.py()+corrJetSmearResDown.mey, 0.,
		sqrt((met.px()+corrJetSmearResDown.mex)*(met.px()+corrJetSmearResDown.mex)+(met.py()+corrJetSmearResDown.mey)*(met.py()+corrJetSmearResDown.mey))),met.vertex());

      pat::MET metJetSme(corrMETJetSmear,met);
      pat::MET metJetRUp(corrMETJetSmearResUp,met);
      pat::MET metJetRDown(corrMETJetSmearResDown,met);
      metJetSmear->push_back(metJetSme);
      metJetSmearResUp->push_back(metJetRUp);
      metJetSmearResDown->push_back(metJetRDown);

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

    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(jetEnUp,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnUp");
      iEvent.put(jetEnDown,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"EnDown");
    }
    iEvent.put(metJetEnUp,std::string(inputMET_.label())+"JetEnUp");
    iEvent.put(metJetEnDown,std::string(inputMET_.label())+"JetEnDown");    
   
    if(storeSmearedShiftedCollections_ == false){
      iEvent.put(jetSmear,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"Smear");
      iEvent.put(jetSmearResUp,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResUp");
      iEvent.put(jetSmearResDown,std::string(jetPSet_.getParameter<edm::InputTag>("src").label())+"SmearJetResDown");
    }
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

  for(size_t ipart = 0 ; ipart < particle.numberOfSourceCandidatePtrs(); ipart++){      
    if(particle.sourceCandidatePtr(ipart).isNonnull() and particle.sourceCandidatePtr(ipart).isAvailable()){
      particles.push_back(particle.sourceCandidatePtr(ipart));
    }
  }

  if(particles.size() == 0)
    return particle.p4();  

  for(unsigned int icand = 0; icand <  pfCandCollection.size(); icand++){
    
    reco::CandidatePtr ptrCand = pfCandCollection.ptrAt(icand);
    const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>(&(pfCandCollection.at(icand)));	
    if(lPack->puppiWeightNoLep() == 0) // in case of puppi, only the one used for MET calculations are useful
      continue;

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
