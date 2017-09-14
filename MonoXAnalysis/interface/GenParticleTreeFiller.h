#ifndef AnalysisCode_MonoXAnalysis_GenParticleTreeFiller_h
#define AnalysisCode_MonoXAnalysis_GenParticleTreeFiller_h

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
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"


#include "TTree.h"
#include "TLorentzVector.h"

class GenParticleTreeFiller {

 public:

  ~GenParticleTreeFiller(){};

  ////// -------
  GenParticleTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);
  ////// -------
  void ReadLHERunProduct(edm::Run const& iRun, float & xsec);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  void findMother(const reco::Candidate*, int &, float &, float &, float &, float &);
  void findFirstNonPhotonMother(const reco::Candidate*, int &, float &, float &, float &, float &);
  float computeDR(const reco::Candidate *genPart,pat::PhotonRef phot);

  const bool isMC;
  const bool useLHEWeights;
  const bool isSignalSample;
  const bool addGenParticles;
  const edm::InputTag lheEventTag;
  const edm::InputTag lheRunTag;
  const bool isQCDTree;
  const bool applyDiMuonFilter;
  const bool applyDiElectronFilter;
  const bool applyPhotonJetsFilter; 
  const bool isTriggerTree; 
  const bool isPhotonPurity;

  edm::EDGetTokenT<GenEventInfoProduct>              genevtInfoToken;
  edm::EDGetTokenT<LHEEventProduct>                  lheInfoToken;
  edm::EDGetTokenT<LHERunInfoProduct>                lheRunInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >    gensToken;
  edm::EDGetTokenT<pat::PhotonRefVector>             photonsPurityToken;

  TTree* tree_;
  bool readDMFromGenParticle_;

  // W/Z boson ifo
  int32_t wzid,l1id,l2id;
  int32_t wzid_h,q1id,q2id;
  int32_t top_1,top_2;  
  int32_t parid,ancid; 

  // gen info leptoni W/Z boson (1 per event)
  float wzmass,wzmt,wzpt,wzeta,wzphi,wzmothid,l1pt,l1eta,l1phi,l2pt,l2eta,l2phi;

  // photon info
  float parpt,pareta,parphi,parmass,ancpt,anceta,ancphi,ancmass;
  int32_t ismatch, isdirect;

  // hadronic V and related quarks (1 per event)
  float wzmass_h,wzmt_h,wzpt_h,wzeta_h,wzphi_h,q1pt,q1eta,q1phi,q2pt,q2eta,q2phi;

  // one top
  float topmass,toppt,topeta,topphi;
  float atopmass,atoppt,atopeta,atopphi;

  // DM mediator and DM particles
  float dmmass,dmpt,dmeta,dmphi,dmX1pt,dmX1eta,dmX1phi,dmX1mass,dmX2pt,dmX2eta,dmX2phi,dmX2mass;
  int   dmid,dmX1id,dmX2id;

  // for fastSIM
  float samplemedM,sampledmM;

  // weights
  float wgt;
  std::vector<float> qcdscalewgt;
  std::vector<int>   qcdscale;
  std::vector<float> gTheta;
  std::vector<float> gDMV;
  std::vector<float> gDMA;
  std::vector<float> gV;
  std::vector<float> gA;
  std::vector<float> couplingwgt;

};

#endif

