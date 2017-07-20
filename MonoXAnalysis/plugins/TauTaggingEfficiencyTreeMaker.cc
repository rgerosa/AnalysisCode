// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"



// Dataformats
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// additional
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

//ROOT
#include "TH2F.h"


class TauTaggingEfficiencyTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
public:

  explicit TauTaggingEfficiencyTreeMaker(const edm::ParameterSet&);
  ~TauTaggingEfficiencyTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;


  // ----------member data ---------------------------
  const edm::InputTag srcTaus;		       
  const edm::EDGetTokenT<std::vector<pat::Tau> > tausToken;
  const double dRClean;
  const bool cleanMuonJet;
  const edm::EDGetTokenT<pat::MuonRefVector> muonsToken;
  const bool cleanElectronJet;
  const edm::EDGetTokenT<pat::ElectronRefVector> electronsToken;
  const bool cleanPhotonJet;
  const edm::EDGetTokenT<pat::PhotonRefVector> photonsToken;

  const std::string   selection;
  const std::vector<edm::ParameterSet>  tauDiscriminatorInfo;
  const std::vector<double> ptBins;
  const std::vector<double> etaBins;
  
  std::map<std::string,TH2F*> eff_histo_Denom_tau;
  std::map<std::string,TH2F*> eff_histo_Num_tau;
  std::map<std::string,TH2F*> eff_histo_Denom_fake;
  std::map<std::string,TH2F*> eff_histo_Num_fake;
};

TauTaggingEfficiencyTreeMaker::TauTaggingEfficiencyTreeMaker(const edm::ParameterSet& iConfig) :
  srcTaus(iConfig.getParameter<edm::InputTag>("srcTaus")),
  tausToken(consumes<std::vector<pat::Tau> >(srcTaus)),	  
  dRClean(iConfig.getParameter<double>("dRClean")),
  cleanMuonJet(iConfig.getParameter<bool>("cleanMuonJet")),
  muonsToken(consumes<pat::MuonRefVector>(iConfig.getParameter<edm::InputTag>("srcMuons"))),
  cleanElectronJet(iConfig.getParameter<bool>("cleanElectronJet")),
  electronsToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("srcElectrons"))),
  cleanPhotonJet(iConfig.getParameter<bool>("cleanPhotonJet")),
  photonsToken(consumes<pat::PhotonRefVector>(iConfig.getParameter<edm::InputTag>("srcPhotons"))),
  selection(iConfig.getParameter<std::string>("selection")),
  tauDiscriminatorInfo(iConfig.getParameter<std::vector<edm::ParameterSet> >("tauDiscriminatorInfo")),
  ptBins(iConfig.getParameter<std::vector<double> >("ptBins")),
  etaBins(iConfig.getParameter<std::vector<double> >("etaBins")){
  usesResource();
  usesResource("TFileService");
}


TauTaggingEfficiencyTreeMaker::~TauTaggingEfficiencyTreeMaker(){}
 
void TauTaggingEfficiencyTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // take jets
  edm::Handle<std::vector<pat::Tau> >tausH;
  iEvent.getByToken(tausToken,tausH);

  // take muons
  edm::Handle<pat::MuonRefVector> muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  pat::MuonRefVector muons = *muonsH;

  // take electrons
  edm::Handle<pat::ElectronRefVector> electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  pat::ElectronRefVector electrons = *electronsH;

  // take photons
  edm::Handle<pat::PhotonRefVector> photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  pat::PhotonRefVector photons = *photonsH;

  // loop over jets
  StringCutObjectSelector<pat::Tau> tauSelection(selection);

  for(auto itTau = tausH->begin(); itTau != tausH->end(); ++itTau){      
    // apply selection
    if(not tauSelection(*itTau)) continue;
    
    //clean from identified and isolated leptons                                                                                                                                                   
    bool skipjet = false;
    if(muonsH.isValid()){
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), itTau->eta(), itTau->phi()) < dRClean)
	  skipjet = true;
      }
    }
    if(electronsH.isValid()){
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), itTau->eta(), itTau->phi()) < dRClean)
	  skipjet = true;
      }
    }
    if(photonsH.isValid()){
      for (std::size_t j = 0; j < photons.size(); j++) {
	if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), itTau->eta(), itTau->phi()) < dRClean)
	  skipjet = true;
      }
    }
    
    if (skipjet) continue;
    
    for(auto iPSet : tauDiscriminatorInfo){       
      std::string discriminatorName = iPSet.getParameter<std::string>("discriminatorName");
      std::string wpLabel = iPSet.getParameter<std::string>("wpLabel");
      std::string decayMode = iPSet.getParameter<std::string>("decayModeFinding");
      
      // to handling correctly the overflow
      double jetpt = itTau->pt();
      if(jetpt > ptBins.back()) jetpt = ptBins.back();
      
      if(not itTau->genJet()){	
	eff_histo_Denom_fake ["eff_"+discriminatorName+"_Denom_fake"]->Fill(jetpt, fabs(itTau->eta()));
	if(itTau->tauID(decayMode) > 0.5 and itTau->tauID(wpLabel) > 0.5)
	  eff_histo_Num_fake ["eff_"+discriminatorName+"_Num_fake"]->Fill(jetpt, fabs(itTau->eta())); // passing the ID
      }
      
      else{
	eff_histo_Denom_tau ["eff_"+discriminatorName+"_Denom_tau"]->Fill(jetpt, fabs(itTau->eta()));
	if(itTau->tauID(decayMode) > 0.5 and itTau->tauID(wpLabel) > 0.5)
	  eff_histo_Num_tau ["eff_"+discriminatorName+"_Num_tau"]->Fill(jetpt, fabs(itTau->eta())); // passing the ID   
      } 
      
    }
  }  
}



void TauTaggingEfficiencyTreeMaker::beginJob(){
  
  edm::Service<TFileService>  fs;
  
  //make histograms in the outputFile
  for(auto iPSet : tauDiscriminatorInfo){

    std::string discriminatorName = iPSet.getParameter<std::string>("discriminatorName");

    std::string name = "eff_"+discriminatorName+"_Denom_tau";
    eff_histo_Denom_tau[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);
    eff_histo_Denom_tau[name] -> Sumw2();

    name = "eff_"+discriminatorName+"_Denom_fake";
    eff_histo_Denom_fake[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);
    eff_histo_Denom_fake[name] -> Sumw2();

    name = "eff_"+discriminatorName+"_Num_tau";
    eff_histo_Num_tau[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);
    eff_histo_Num_tau[name] -> Sumw2();

    name = "eff_"+discriminatorName+"_Num_fake";
    eff_histo_Num_fake[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);
    eff_histo_Num_fake[name] -> Sumw2();

  }  
}

void TauTaggingEfficiencyTreeMaker::endJob() {}

void TauTaggingEfficiencyTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void TauTaggingEfficiencyTreeMaker::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void TauTaggingEfficiencyTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TauTaggingEfficiencyTreeMaker);
