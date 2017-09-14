#include "AnalysisCode/MonoXAnalysis/interface/TreeFillerUtils.h"


// to apply jet ID                                                                                                                                                                                  
bool applyJetID(const pat::Jet & jet, const std::string & level){
  
  if(level != "loose" and level != "tight" and level != "tightLepVeto")
    return true;
  
  bool passjetid = false;
  
  //apply a loose jet id https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
  if(level == "loose"){ 
    if (fabs(jet.eta()) <= 2.7 and
	jet.neutralHadronEnergyFraction() < 0.99 and
	jet.neutralEmEnergyFraction()     < 0.99 and
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4)
	passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 and
	       jet.chargedHadronEnergyFraction() > 0. and
	       jet.chargedEmEnergyFraction()     < 0.99 and 
	       jet.chargedMultiplicity()         > 0) 
	passjetid = true;
    }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) <= 3.0 and
	     jet.neutralHadronEnergyFraction() < 0.98 and
             jet.neutralEmEnergyFraction() > 0.01 and
             jet.neutralMultiplicity()     > 2)
      passjetid = true;  
    else if(fabs(jet.eta()) > 3.0 and
	    jet.neutralEmEnergyFraction() < 0.9 and
	    jet.neutralMultiplicity()     > 10)
      passjetid = true; 
  } 
  else if(level == "tight"){
    
    if (fabs(jet.eta()) <= 2.7 and 
	jet.neutralHadronEnergyFraction() < 0.90 and 
	jet.neutralEmEnergyFraction()     < 0.90 and 
	(jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4) 
	passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 and 
	       jet.chargedHadronEnergyFraction() > 0. and 
	       jet.chargedEmEnergyFraction() < 0.99 and 
	       jet.chargedMultiplicity() > 0) 
	passjetid = true;
    }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) < 3.0  and
	     jet.neutralHadronEnergyFraction() < 0.98 and
             jet.neutralEmEnergyFraction() > 0.01 and
             jet.neutralMultiplicity()     > 2)
      passjetid = true;
    else if(fabs(jet.eta()) > 3.0 and
	    jet.neutralEmEnergyFraction() < 0.9 and
	    jet.neutralMultiplicity() > 10)
      passjetid = true;    
  }
  
  else if(level == "tightLepVeto"){
    if (fabs(jet.eta()) <= 2.7 and
        jet.neutralHadronEnergyFraction() < 0.90 and
        jet.neutralEmEnergyFraction() < 0.90 and
	jet.muonEnergyFraction() < 0.80 and 
        (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1) {
      
      if (fabs(jet.eta()) > 2.4)
        passjetid = true;
      else if (fabs(jet.eta()) <= 2.4 and
               jet.chargedHadronEnergyFraction() > 0. and
               jet.chargedEmEnergyFraction() < 0.90 and
               jet.chargedMultiplicity() > 0)
	passjetid = true;
    }
    else if (fabs(jet.eta()) > 2.7 and fabs(jet.eta()) < 3.0 and
	     jet.neutralHadronEnergyFraction() < 0.98 and
             jet.neutralEmEnergyFraction() > 0.01 and
             jet.neutralMultiplicity()     > 2)
      passjetid = true;
    else if (fabs(jet.eta()) > 3.0 and 
	     jet.neutralEmEnergyFraction() < 0.9 and
	     jet.neutralMultiplicity() > 10)
      passjetid = true;    
    
  }  
  return passjetid;  

}

// to apply pileup-jet id                                                                                                                                                                            
bool applyPileupJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi){

  bool passpuid    = false;
  float puidval   = 0;
  float jetabseta = fabs(jet.eta());
  float jetpt     = jet.pt();

  if(jet.hasUserFloat("puid:fullDiscriminant"))
    puidval = jet.userFloat("puid:fullDiscriminant");
  else if(jet.hasUserFloat("puidPuppi:fullDiscriminant"))
    puidval = jet.userFloat("puidPuppi:fullDiscriminant");
  else if(jet.hasUserFloat("pileupJetId:fullDiscriminant"))
    puidval = jet.userFloat("pileupJetId:fullDiscriminant");
  else 
    return true;

  // from twiki: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID --> to be loaded in GT soon
  if(level == "loose"){

    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > -0.97) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > -0.89) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.68) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.52) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.53) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.38) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.47) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.30) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;

  }
  else if(level == "medium"){ 

    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > 0.18) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > 0.61) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.55) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.42) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.23) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.36) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.17) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;

  }

  else if(level == "tight"){
    if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 10 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 20 and jetpt > 10 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 30 and jetpt > 20 and puidval > 0.69) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt < 50 and jetpt > 30 and puidval > 0.86) passpuid = true;
    else if (jetabseta >= 0.00 and jetabseta < 2.50 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 10 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 20 and jetpt > 10 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 30 and jetpt > 20 and puidval > -0.35) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt < 50 and jetpt > 30 and puidval > -0.10) passpuid = true;
    else if (jetabseta >= 2.50 and jetabseta < 2.75 and jetpt > 50) passpuid = true;

    if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 10 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 20 and jetpt > 10 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 30 and jetpt > 20 and puidval > -0.21) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt < 50 and jetpt > 30 and puidval > -0.01) passpuid = true;
    else if (jetabseta >= 2.75 and jetabseta < 3.00 and jetpt > 50) passpuid = true;

    if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 10 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 20 and jetpt > 10 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 30 and jetpt > 20 and puidval > -0.26) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt < 50 and jetpt > 30 and puidval > -0.03) passpuid = true;
    else if (jetabseta >= 3.00 and jetabseta < 5.00 and jetpt > 50) passpuid = true;
  }

  return passpuid;

}
