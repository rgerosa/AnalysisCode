#ifndef XConeJetProducer_h
#define XConeJetProducer_h

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/contrib/Nsubjettiness.hh" 
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/XConePlugin.hh"

#include "TMath.h"

using namespace std;

class XConeJetProducer : public edm::stream::EDProducer<> {
  
public:

  XConeJetProducer(){};
  XConeJetProducer(edm::ParameterSet const & params );  
  virtual ~XConeJetProducer() {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
private:
  
  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;      

  template <typename T>
  void  writeJets(edm::Event & iEvent,
		  const edm::EventSetup & iSetup,
		  const reco::CandidateView & pfCand,
		  const reco::VertexCollection & vertex,
		  const bool & isXCone = true);

  const edm::EDGetTokenT<edm::View<reco::Candidate> >  srcToken_;
  const edm::EDGetTokenT<reco::VertexCollection>  srcVertexToken_;
  const edm::ParameterSet xConeParameters_;
  bool addArea_;
  bool addNsubjettnessJets_;
  std::string jetType_;
  edm::ParameterSet nSubjettinessParameters_;
  double GhostEtaMax_;
  std::string jetSelection_;

  // fastjet shared ptr
  std::shared_ptr<fastjet::contrib::XConePlugin > xCone_;
  std::shared_ptr<fastjet::JetDefinition > xConeJetDef_;
  std::shared_ptr<fastjet::AreaDefinition > area_def_;
  std::shared_ptr<fastjet::ClusterSequence> xConeSeq_;
  std::shared_ptr<fastjet::ClusterSequenceArea> xConeSeqArea_;

  std::shared_ptr<fastjet::contrib::NjettinessPlugin > Njettiness_;
  std::shared_ptr<fastjet::contrib::OnePass_GenET_GenKT_Axes> jetAxis_;
  std::shared_ptr<fastjet::JetDefinition> NJettinessDef_;
  std::shared_ptr<fastjet::contrib::UnnormalizedCutoffMeasure> jetMeasureUnnormalized_;
  std::shared_ptr<fastjet::contrib::ConicalMeasure> jetMeasureConical_;
  std::shared_ptr<fastjet::contrib::OriginalGeometricMeasure> jetMeasureOriginalGeometric_;
  std::shared_ptr<fastjet::contrib::ModifiedGeometricMeasure> jetMeasureModifiedGeometric_;
  std::shared_ptr<fastjet::contrib::ConicalGeometricMeasure> jetMeasureConicalGeometric_;
  std::shared_ptr<fastjet::ClusterSequence> NJettinessSeq_;
  std::shared_ptr<fastjet::ClusterSequenceArea> NJettinessAreaSeq_;

  // candidates
  std::vector<fastjet::PseudoJet> fjcandidate_;
  std::vector<fastjet::PseudoJet> xConeJets_;
  std::vector<fastjet::PseudoJet> xConeJetsConstituent_;
  std::vector<fastjet::PseudoJet> xConeJetsGhosts_;
  std::vector<fastjet::PseudoJet> xConeJetsParticles_;
  std::vector<fastjet::PseudoJet> NjettinssJets_;
  std::vector<reco::CandidatePtr> pfjcand_;
};

#endif

XConeJetProducer::XConeJetProducer(edm::ParameterSet const & params ):
  srcToken_( consumes<edm::View<reco::Candidate> >( params.getParameter<edm::InputTag>("src") ) ), // input particles to be clustered
  srcVertexToken_( consumes<reco::VertexCollection>( params.getParameter<edm::InputTag>("srcVertex") ) ), // input vertex
  xConeParameters_(params.getParameter<edm::ParameterSet>("xConeParameters")),
  addArea_(params.existsAs<bool> ("addArea") ? params.getParameter<bool>("addArea") : false),
  addNsubjettnessJets_(params.existsAs<bool> ("addNsubjettnessJets") ? params.getParameter<bool>("addNsubjettnessJets") : false),
  jetType_(params.getParameter<std::string>("jetType")),
  GhostEtaMax_(params.existsAs<double>("GhostEtaMax") ? params.getParameter<double>("GhostEtaMax") : 5.0), 
  jetSelection_(params.existsAs<std::string>("jetSelection") ? params.getParameter<std::string>("jetSelection") : ""){
		   
  // Define XCone properties
  double beta = xConeParameters_.existsAs<double>("beta") ? xConeParameters_.getParameter<double>("beta") : -99. ;
  if(beta  == -99)
    throw cms::Exception("Configuration","xConeParameters_ beta value not found");
  double Rcut = xConeParameters_.existsAs<double>("Rcut") ? xConeParameters_.getParameter<double>("Rcut") : -99. ;
  if(Rcut  == -99)
    throw cms::Exception("Configuration","xConeParameters_ Rcut value not found");
  int N = xConeParameters_.existsAs<unsigned int>("N") ? xConeParameters_.getParameter<unsigned int>("N") : -99. ;
  if(N  == -99)
    throw cms::Exception("Configuration","xConeParameters_ N value not found");

  // define properties of XCone jets
  xCone_ = std::shared_ptr<fastjet::contrib::XConePlugin > ( new fastjet::contrib::XConePlugin (N,Rcut,beta));
  xConeJetDef_ = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition (xCone_.get()));
  if(addArea_)
    area_def_  = std::shared_ptr<fastjet::AreaDefinition>(new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(GhostEtaMax_)));  
  
  // check and define Nsubjettiness properties
  if(addNsubjettnessJets_){
    nSubjettinessParameters_ = params.getParameter<edm::ParameterSet>("nSubjettinessParameters");
    double delta = nSubjettinessParameters_.existsAs<double>("delta") ? nSubjettinessParameters_.getParameter<double>("delta") : -99. ;
    if(delta  == -99)
      throw cms::Exception("Configuration","nSubjettinessParameters_ delta value not found");	  
    double p = nSubjettinessParameters_.existsAs<double>("p") ? nSubjettinessParameters_.getParameter<double>("delta") : -99. ;
    if(p  == -99)
      throw cms::Exception("Configuration","nSubjettinessParameters_ p value not found");	  
    Rcut = nSubjettinessParameters_.existsAs<double>("Rcut") ? nSubjettinessParameters_.getParameter<double>("Rcut") : -99. ;
    if(Rcut  == -99)
      throw cms::Exception("Configuration","nSubjettinessParameters_ Rcut value not found");	  
    N = nSubjettinessParameters_.existsAs<int>("N") ? nSubjettinessParameters_.getParameter<int>("N") : -99. ;
    if(N  == -99)
      throw cms::Exception("Configuration","nSubjettinessParameters_ N value not found");
    std::string measureType = nSubjettinessParameters_.existsAs<std::string>("measureType") ? nSubjettinessParameters_.getParameter<std::string>("measureType") : "" ;
    if(measureType == "")
      throw cms::Exception("Configuration","nSubjettinessParameters_ measureType string not found");

    jetAxis_ = std::shared_ptr<fastjet::contrib::OnePass_GenET_GenKT_Axes>(new fastjet::contrib::OnePass_GenET_GenKT_Axes(delta,p,Rcut));
    
    // measure definition
    if(measureType == "UnnormalizedCutoffMeasure"){
      jetMeasureUnnormalized_ = std::shared_ptr<fastjet::contrib::UnnormalizedCutoffMeasure> (new fastjet::contrib::UnnormalizedCutoffMeasure(delta,Rcut,fastjet::contrib::pt_R));
      Njettiness_    = std::shared_ptr<fastjet::contrib::NjettinessPlugin> (new fastjet::contrib::NjettinessPlugin(N,*jetAxis_,*jetMeasureUnnormalized_));
      NJettinessDef_ = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(Njettiness_.get()));
    }
    else if(measureType == "ConicalMeasure"){
      jetMeasureConical_ = std::shared_ptr<fastjet::contrib::ConicalMeasure> (new fastjet::contrib::ConicalMeasure(delta,Rcut));
      Njettiness_        = std::shared_ptr<fastjet::contrib::NjettinessPlugin> (new fastjet::contrib::NjettinessPlugin(N,*jetAxis_,*jetMeasureConical_));
      NJettinessDef_     = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(Njettiness_.get()));
    }
    else if(measureType == "OriginalGeometricMeasure"){
      jetMeasureOriginalGeometric_ = std::shared_ptr<fastjet::contrib::OriginalGeometricMeasure> (new fastjet::contrib::OriginalGeometricMeasure(Rcut));
      Njettiness_                  = std::shared_ptr<fastjet::contrib::NjettinessPlugin> (new fastjet::contrib::NjettinessPlugin(N,*jetAxis_,*jetMeasureOriginalGeometric_));
      NJettinessDef_               = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(Njettiness_.get()));
    }	  
    else if(measureType == "ModifiedGeometricMeasure"){
      jetMeasureModifiedGeometric_ = std::shared_ptr<fastjet::contrib::ModifiedGeometricMeasure>(new fastjet::contrib::ModifiedGeometricMeasure(Rcut));
      Njettiness_                  = std::shared_ptr<fastjet::contrib::NjettinessPlugin> (new fastjet::contrib::NjettinessPlugin(N,*jetAxis_,*jetMeasureModifiedGeometric_));
      NJettinessDef_ = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(Njettiness_.get()));
    }
    else{
      jetMeasureConicalGeometric_ = std::shared_ptr<fastjet::contrib::ConicalGeometricMeasure> (new fastjet::contrib::ConicalGeometricMeasure(delta,p,Rcut));
      Njettiness_                 = std::shared_ptr<fastjet::contrib::NjettinessPlugin> (new fastjet::contrib::NjettinessPlugin(N,*jetAxis_,*jetMeasureConicalGeometric_));
      NJettinessDef_ = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(Njettiness_.get()));
    }      
  


  }
  // loop over parmeter list to checkout name
  if(jetType_ == "PFJet" or jetType_ == "RecoJet"){
    produces<reco::PFJetCollection> ();
    if(addNsubjettnessJets_)
      produces<reco::PFJetCollection> ("nSub");
  }
  else if(jetType_ == "GenJet"){
    produces<reco::GenJetCollection>();
    if(addNsubjettnessJets_)
      produces<reco::GenJetCollection> ("nSub");
  }
  else
    throw cms::Exception("Configuration","requested jetType not supported : "+jetType_+" ... please check");		
  
}
  

void XConeJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::CandidateView> candidate;
  iEvent.getByToken( srcToken_, candidate);
  if(not candidate.isValid())
    throw cms::Exception("Configuration","Invalid list of input candidates --> please check ");

  edm::Handle<reco::VertexCollection> vertex;
  iEvent.getByToken( srcVertexToken_, vertex);
  if(not vertex.isValid())
    throw cms::Exception("Configuration","Invalid list of primary vertexes --> please check ");
    
  // convert input candidates into fastjet pseudoket for the clustering
  fjcandidate_.clear();
  for(auto cand_iter = candidate->begin(); cand_iter != candidate->end(); cand_iter++){
    fjcandidate_.push_back(fastjet::PseudoJet(cand_iter->px(),cand_iter->py(),cand_iter->pz(),cand_iter->energy()));
    fjcandidate_.back().set_user_index(cand_iter-candidate->begin());
  }

  //define xCone cluster sequence
  if(not addArea_){  
    xConeSeq_  = std::shared_ptr<fastjet::ClusterSequence> (new fastjet::ClusterSequence(fjcandidate_,*xConeJetDef_));  
    // cluster particles and take jets
    xConeJets_.clear();
    xConeJets_ = fastjet::sorted_by_pt(xConeSeq_->inclusive_jets());
  }
  else{
    xConeSeqArea_ = std::shared_ptr<fastjet::ClusterSequenceArea>(new fastjet::ClusterSequenceArea(fjcandidate_,*xConeJetDef_,*area_def_));
    xConeJets_.clear();
    xConeJets_ = fastjet::sorted_by_pt(xConeSeqArea_->inclusive_jets());
  }
  

  //write jets
  if(jetType_ == "RecoJet" or jetType_ == "PFJet")
    writeJets<reco::PFJet>(iEvent,iSetup,*candidate.product(),*vertex.product());
  else if(jetType_ == "GenJet")
    writeJets<reco::GenJet>(iEvent,iSetup,*candidate.product(),*vertex.product());

  if(addNsubjettnessJets_){  
    // run more generic N-subjettiness finder
    if(not addArea_){
      NJettinessSeq_ = std::shared_ptr<fastjet::ClusterSequence> (new fastjet::ClusterSequence(fjcandidate_,*NJettinessDef_));
      NjettinssJets_.clear();
      NjettinssJets_ = fastjet::sorted_by_pt(NJettinessSeq_->inclusive_jets());
    }
    else{
      NJettinessAreaSeq_ = std::shared_ptr<fastjet::ClusterSequenceArea>( new fastjet::ClusterSequenceArea(fjcandidate_,*NJettinessDef_,*area_def_));
      NjettinssJets_.clear();
      NjettinssJets_ = fastjet::sorted_by_pt(NJettinessAreaSeq_->inclusive_jets());
    }
        
    if(jetType_ == "RecoJet" or jetType_ == "PFJet")
      writeJets<reco::PFJet>(iEvent,iSetup,*candidate.product(),*vertex.product(),false);
    else if(jetType_ == "GenJet")
      writeJets<reco::GenJet>(iEvent,iSetup,*candidate.product(),*vertex.product(),false);
  }   
}


template <typename T>
void  XConeJetProducer::writeJets(edm::Event & iEvent, 
				  const edm::EventSetup & iSetup,
				  const reco::CandidateView & pfCand, 
				  const reco::VertexCollection & vertex,
				  const bool & isXCone){

  StringCutObjectSelector<T> jetSelctor (jetSelection_); // in order to apply selections, typically eta and pt
  std::auto_ptr<std::vector<T> > outputJets(new std::vector<T>() );
  outputJets->clear();

  std::vector<fastjet::PseudoJet> jets;
  if(isXCone)
    jets = xConeJets_;
  else
    jets = NjettinssJets_;

  // loop over the clustered jets
  for(auto fjet : jets){
    // get constituent from fastJet
    xConeJetsConstituent_.clear();
    xConeJetsConstituent_ = fastjet::sorted_by_pt(fjet.constituents());
    // filter ghosts if any
    xConeJetsGhosts_.clear();
    xConeJetsParticles_.clear();
    fastjet::SelectorIsPureGhost().sift(xConeJetsConstituent_, xConeJetsGhosts_, xConeJetsParticles_);    
    // loop over constituent
    pfjcand_.clear();
    for(auto fjconst : xConeJetsParticles_){
      int index = fjconst.user_index();
      pfjcand_.push_back(pfCand.ptrAt(index));
    }

    if(pfjcand_.size() != xConeJetsParticles_.size())
      throw cms::Exception("writeJets","PFCandidate and fastjet particle collection have a different size");

    // write jet info JetSpecific.h
    T jet;
    if(jetType_ != "GenJet")
      writeSpecific(jet,reco::Particle::LorentzVector(fjet.px(),fjet.py(),fjet.pz(),fjet.E()),vertex.begin()->position(),pfjcand_,iSetup);
    else
      writeSpecific(jet,reco::Particle::LorentzVector(fjet.px(),fjet.py(),fjet.pz(),fjet.E()),reco::Jet::Point(0,0,0),pfjcand_,iSetup);

    double jetArea = 0;
    if(fjet.has_area())
      jetArea = fjet.area();
    else{
      if(isXCone)
	jetArea = TMath::Pi()*xConeParameters_.getParameter<double>("Rcut");
      else
	jetArea = TMath::Pi()*nSubjettinessParameters_.getParameter<double>("Rcut");
    }
    jet.setJetArea(jetArea);

    // apply selection
    if(not jetSelctor(jet)) continue;
    outputJets->push_back(jet);    
  }

  if(isXCone)
    iEvent.put(outputJets);
  else
    iEvent.put(outputJets,"nSub");  
}
 
void XConeJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(XConeJetProducer);
