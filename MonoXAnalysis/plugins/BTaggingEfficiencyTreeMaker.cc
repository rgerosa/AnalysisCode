// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TH2F.h"


class BTaggingEfficiencyTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
public:

  explicit BTaggingEfficiencyTreeMaker(const edm::ParameterSet&);
  ~BTaggingEfficiencyTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // to apply jet id
  bool applyJetID(const pat::Jet & jet, const std::string & level);
  // to apply pu jet id
  bool applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi);

  // ----------member data ---------------------------
  const std::string   dirName;
  const edm::InputTag srcJets;		       
  const edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
  const double dRClean;
  const bool cleanMuonJet;
  const edm::EDGetTokenT<pat::MuonRefVector> muonsToken;
  const bool cleanElectronJet;
  const edm::EDGetTokenT<pat::ElectronRefVector> electronsToken;
  const bool cleanPhotonJet;
  const edm::EDGetTokenT<pat::PhotonRefVector> photonsToken;

  const std::string   selection;
  const std::vector<edm::ParameterSet>  bDiscriminatorInfo;
  const std::vector<double> ptBins;
  const std::vector<double> etaBins;

  std::map<std::string,TH2F*> eff_histo_Denom_b;
  std::map<std::string,TH2F*> eff_histo_Denom_c;
  std::map<std::string,TH2F*> eff_histo_Denom_ucsdg;

  std::map<std::string,TH2F*> eff_histo_Num_b;
  std::map<std::string,TH2F*> eff_histo_Num_c;
  std::map<std::string,TH2F*> eff_histo_Num_ucsdg;
};

BTaggingEfficiencyTreeMaker::BTaggingEfficiencyTreeMaker(const edm::ParameterSet& iConfig) :
  dirName(iConfig.getParameter<std::string>("directoryName")),
  srcJets(iConfig.getParameter<edm::InputTag>("srcJets")),
  jetsToken(consumes<std::vector<pat::Jet> >(srcJets)),	  
  dRClean(iConfig.getParameter<double>("dRClean")),
  cleanMuonJet(iConfig.getParameter<bool>("cleanMuonJet")),
  muonsToken(consumes<pat::MuonRefVector>(iConfig.getParameter<edm::InputTag>("srcMuons"))),
  cleanElectronJet(iConfig.getParameter<bool>("cleanElectronJet")),
  electronsToken(consumes<pat::ElectronRefVector>(iConfig.getParameter<edm::InputTag>("srcElectrons"))),
  cleanPhotonJet(iConfig.getParameter<bool>("cleanPhotonJet")),
  photonsToken(consumes<pat::PhotonRefVector>(iConfig.getParameter<edm::InputTag>("srcPhotons"))),
  selection(iConfig.getParameter<std::string>("selection")),
  bDiscriminatorInfo(iConfig.getParameter<std::vector<edm::ParameterSet> >("bDiscriminatorInfo")),
  ptBins(iConfig.getParameter<std::vector<double> >("ptBins")),
  etaBins(iConfig.getParameter<std::vector<double> >("etaBins")){

  usesResource();
  usesResource("TFileService");
}


BTaggingEfficiencyTreeMaker::~BTaggingEfficiencyTreeMaker(){}
 
void BTaggingEfficiencyTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // take jets
  edm::Handle<std::vector<pat::Jet> >jetsH;
  iEvent.getByToken(jetsToken,jetsH);

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
  StringCutObjectSelector<pat::Jet> jetSelection(selection);

  for(auto itJet = jetsH->begin(); itJet != jetsH->end(); ++itJet){

    // apply selection
    if(not jetSelection(*itJet)) continue;

    // make sure is from ME
    bool jetfromME = false;      
    if(itJet->genJet()) 
      jetfromME = true;
    if(not jetfromME)
      continue;

    //clean from leptons                                                                                                                                                     
    bool skipjet = false;
    if(muonsH.isValid()){
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), itJet->eta(), itJet->phi()) < dRClean)
	  skipjet = true;
      }
    }
    if(electronsH.isValid()){
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), itJet->eta(), itJet->phi()) < dRClean)
	  skipjet = true;
      }
    }
    if(photonsH.isValid()){
      for (std::size_t j = 0; j < photons.size(); j++) {
	if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), itJet->eta(), itJet->phi()) < dRClean)
	  skipjet = true;
      }
    }

    if (skipjet) continue;

    // require jet id
    bool passjetid = applyJetID(*itJet,"loose");
    if (!passjetid)
      continue;

    // apply pileup jet id                                                                                                                                                   
    bool isPuppi = false;
    if(TString(srcJets.label()).Contains("Puppi"))
      isPuppi = true;

    bool passpuid = applyPileupJetID(*itJet,"medium",isPuppi);
    if (!passpuid)
      continue;

    int hadronFlavor = itJet->hadronFlavour();      
    
    for(auto iPSet : bDiscriminatorInfo){ 
      
      std::string discriminatorName = iPSet.getParameter<std::string>("discriminatorName");
      std::string wpLabel = iPSet.getParameter<std::string>("wpLabel");
      double wpValue = iPSet.getParameter<double>("wpValue");

      if( abs(hadronFlavor)==5 ){
	// fill denominator
	eff_histo_Denom_b["eff_"+discriminatorName+"_"+wpLabel+"_Denom_b"]->Fill(itJet->pt(), fabs(itJet->eta()));
	// selection
	if(itJet->bDiscriminator(discriminatorName) > wpValue)
	  eff_histo_Num_b["eff_"+discriminatorName+"_"+wpLabel+"_Num_b"]->Fill(itJet->pt(), fabs(itJet->eta()));
      }
      else if( abs(hadronFlavor)==4 ){
	// fill denominator
	eff_histo_Denom_c["eff_"+discriminatorName+"_"+wpLabel+"_Denom_c"]->Fill(itJet->pt(), fabs(itJet->eta()));
	// selection
	if(itJet->bDiscriminator(discriminatorName) > wpValue)
	  eff_histo_Num_c["eff_"+discriminatorName+"_"+wpLabel+"_Num_c"]->Fill(itJet->pt(), fabs(itJet->eta()));
      }
      else{
	// fill denominator
	eff_histo_Denom_ucsdg["eff_"+discriminatorName+"_"+wpLabel+"_Denom_ucsdg"]->Fill(itJet->pt(), fabs(itJet->eta()));
	// selection
	if(itJet->bDiscriminator(discriminatorName) > wpValue)
	  eff_histo_Num_ucsdg["eff_"+discriminatorName+"_"+wpLabel+"_Num_ucsdg"]->Fill(itJet->pt(), fabs(itJet->eta()));
      }
    }
  }
}



void BTaggingEfficiencyTreeMaker::beginJob(){
  
  edm::Service<TFileService>  fs;
  
  //make histograms in the outputFile
  for(auto iPSet : bDiscriminatorInfo){

    std::string discriminatorName = iPSet.getParameter<std::string>("discriminatorName");
    std::string wpLabel = iPSet.getParameter<std::string>("wpLabel");

    std::string name = "eff_"+discriminatorName+"_"+wpLabel+"_Denom_b";
    eff_histo_Denom_b[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);

    name = "eff_"+discriminatorName+"_"+wpLabel+"_Denom_c";
    eff_histo_Denom_c[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);

    name = "eff_"+discriminatorName+"_"+wpLabel+"_Denom_ucsdg";
    eff_histo_Denom_ucsdg[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);

    name = "eff_"+discriminatorName+"_"+wpLabel+"_Num_b";
    eff_histo_Num_b[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);

    name = "eff_"+discriminatorName+"_"+wpLabel+"_Num_c";
    eff_histo_Num_c[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);

    name = "eff_"+discriminatorName+"_"+wpLabel+"_Num_ucsdg";
    eff_histo_Num_ucsdg[name] = fs->make<TH2F>(name.c_str(),"", ptBins.size()-1, &ptBins[0], etaBins.size()-1,&etaBins[0]);
  }
  
}

void BTaggingEfficiencyTreeMaker::endJob() {}

void BTaggingEfficiencyTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void BTaggingEfficiencyTreeMaker::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void BTaggingEfficiencyTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool BTaggingEfficiencyTreeMaker::applyJetID(const pat::Jet & jet, const std::string & level){
  
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


bool BTaggingEfficiencyTreeMaker::applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi){
  
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

DEFINE_FWK_MODULE(BTaggingEfficiencyTreeMaker);
