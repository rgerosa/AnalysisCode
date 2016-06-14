#include "plot.h"

void showPlots() {

  gROOT->SetBatch(kTRUE);
      
  //  plot(Sample::zmm, "zmass"     , ""       , "M_{#mu#mu} [GeV]"        , "Events", "outputPlots","zmm_zmass" , 2.1, 30, 60, 120, 0.05, 7e6,false);
  //  plot(Sample::zmm, "nvtx"      , ""       , "N_{PV}"                  , "Events", "outputPlots","zmm_nvtx"  , 2.1, 25, 5.5, 30.5, 0.05, 7e6,false);
  //  plot(Sample::zmm, "zpt"       , ""       , "p_{T}^{Z} [GeV]"         , "Events", "outputPlots","zmm_pt"    , 2.1, 35, 0,  250, 0.05, 7e6,false);
  /*
  plot(Sample::zmm, "zeta"       , ""       , "#eta_{Z}"         , "Events", "outputPlots","zmm_eta"    , 0.85, 35, -4,  4, 0.05, 7e6,true);
  plot(Sample::zmm, "mu1pt"     , ""       , "p_{T}^{#mu} [GeV]"       , "Events", "outputPlots","zmm_mu1pt" , 0.85, 35, 20,  200, 0.05, 7e6,true);
  plot(Sample::zmm, "mu1eta"     , ""       , "#eta_{#mu}"       , "Events", "outputPlots","zmm_mu1eta" , 0.85, 35, -2.4,  2.4, 0.05, 7e6,true);
  plot(Sample::zmm, "mu2pt"     , ""       , "p_{T}^{#mu} [GeV]"       , "Events", "outputPlots","zmm_mu2pt" , 0.85, 35, 0,  150, 0.05, 7e6,true);
  plot(Sample::zmm, "mu2eta"     , ""       , "#eta_{#mu}"       , "Events", "outputPlots","zmm_mu2eta" , 0.85, 35, -2.4,  2.4, 0.05, 7e6,true);
  plot(Sample::zmm, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","zmm_t1recoil"  , 0.85, 25, 0, 400, 0.05, 7e6,true);
  plot(Sample::zmm, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","zmm_t1met"  ,    0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","zmm_met",        0.85, 25, 0, 150, 0.05, 1e5,true);  
  /*
  plot(Sample::zmm, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","zmm_met_hadronHF"  ,   0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","zmm_met_egammaHF", 0.85, 20, 0, 50, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","zmm_met_chargeHadron",  0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","zmm_met_neutralHadron", 0.85, 25, 0, 120, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","zmm_met_electrons", 0.85, 20, 0, 100, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","zmm_met_muons",     0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zmm, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","zmm_met_photons",   0.85, 25, 0, 150, 0.05, 7e6,true);
  */  
  /*  
  plot(Sample::zee, "zeemass"     , ""       , "M_{ee} [GeV]"        , "Events", "outputPlots","zee_zmass" , 0.84, 30, 60, 120, 0.05, 7e6,true);
  plot(Sample::zee, "nvtx"        , ""       , "N_{PV}"                  , "Events", "outputPlots","zee_nvtx"  , 0.84, 25, 5.5, 30.5, 0.05, 7e6,true);
  plot(Sample::zee, "zeept"       , ""       , "p_{T}^{Z} [GeV]"         , "Events", "outputPlots","zee_pt"    , 0.84, 35, 0,  250, 0.05, 7e6,true);
  plot(Sample::zee, "zeeeta"      , ""       , "#eta_{Z}"         , "Events", "outputPlots","zee_eta"    , 0.84, 35, -4,  4, 0.05, 7e6,true);
  plot(Sample::zee, "el1pt"     , ""       , "p_{T}^{e} [GeV]"       , "Events", "outputPlots","zee_el1pt" , 0.84, 35, 40,  200, 0.05, 7e6,true);
  plot(Sample::zee, "el1eta"     , ""       , "#eta_{e}"       , "Events", "outputPlots","zee_el1eta" , 0.84, 35, -2.5,  2.5, 0.05, 7e6,true);
  plot(Sample::zee, "el2pt"     , ""       , "p_{T}^{e} [GeV]"       , "Events", "outputPlots","zee_el2pt" , 0.84, 35, 20,  200, 0.05, 7e6,true);
  plot(Sample::zee, "el2eta"     , ""       , "#eta_{e}"       , "Events", "outputPlots","zee_el2eta" , 0.84, 35, -2.5,  2.5, 0.05, 7e6,true);
  plot(Sample::zee, "t1elmet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","zee_t1recoil"  , 0.84, 25, 0, 400, 0.05, 7e6,true);
  plot(Sample::zee, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","zee_t1met"  ,    0.84, 25, 0, 150, 0.05, 7e6,true);
  /*
  plot(Sample::zee, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","zee_met",        0.84, 25, 0, 150, 0.05, 1e5,true);  
  plot(Sample::zee, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","zee_met_hadronHF"  ,   0.84, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","zee_met_egammaHF", 0.84, 20, 0, 50, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","zee_met_chargeHadron",  0.84, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","zee_met_neutralHadron", 0.84, 25, 0, 120, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","zee_met_electrons", 0.84, 20, 0, 100, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","zee_met_muons",     0.84, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::zee, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","zee_met_photons",   0.84, 25, 0, 150, 0.05, 7e6,true);

  */
  //plot(Sample::wmn, "nvtx"      , ""       , "N_{PV}"                  , "Events", "outputPlots","wmn_nvtx"  , 0.85, 25, 5.5, 30.5, 0.05, 7e6,true);
  //plot(Sample::wmn, "mu1pt"     , ""       , "p_{T}^{#mu} [GeV]"       , "Events", "outputPlots","wmn_mu1pt" , 0.85, 35, 20,  200, 0.05, 7e6,true);
  //plot(Sample::wmn, "mu1eta"     , ""       , "#eta_{#mu}"       , "Events", "outputPlots","wmn_mu1eta" , 0.85, 35, -2.4,  2.4, 0.05, 7e6,true);
  //plot(Sample::wmn, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","wmn_t1recoil"  , 0.86, 25, 200, 1000, 0.05, 7e6,false);
  //  plot(Sample::wen, "t1elmet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","wen_t1recoil"  , 0.86, 25, 200, 1000, 0.05, 7e6,false);
  plot(Sample::gam, "t1phmet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","gam_t1recoil"  , 0.816, 25, 200, 1200, 0.05, 7e6,false);
  //plot(Sample::wen, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","wmn_t1met"  ,    0.85, 25, 0, 150, 0.05, 7e6,true);
  /*  
  plot(Sample::wmn, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","wmn_met",        0.85, 25, 0, 150, 0.05, 1e5,true);  
  plot(Sample::wmn, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","wmn_met_hadronHF"  ,   0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","wmn_met_egammaHF", 0.85, 20, 0, 50, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","wmn_met_chargeHadron",  0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","wmn_met_neutralHadron", 0.85, 25, 0, 120, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","wmn_met_electrons", 0.85, 20, 0, 100, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","wmn_met_muons",     0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wmn, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","wmn_met_photons",   0.85, 25, 0, 150, 0.05, 7e6,true);
  */
  /*  
  plot(Sample::wen, "nvtx"      , ""       , "N_{PV}"                 , "Events", "outputPlots","wen_nvtx"  , 0.85, 25, 5.5, 30.5, 0.05, 7e6,true);
  plot(Sample::wen, "el1pt"     , ""       , "p_{T}^{el} [GeV]"       , "Events", "outputPlots","wen_el1pt" , 0.85, 35, 40,  200, 0.05, 7e6,true);
  plot(Sample::wen, "el1eta"     , ""       , "#eta_{el}"             , "Events", "outputPlots","wen_el1eta" , 0.85, 35, -2.5,  2.5, 0.05, 7e6,true);
  plot(Sample::wen, "t1elmet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","wen_t1recoil"  , 0.85, 25, 0, 400, 0.05, 7e6,true);
  plot(Sample::wen, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","wen_t1met"  ,    0.85, 25, 50, 200, 0.05, 7e6,true);
  /*
  plot(Sample::wen, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","wen_met",        0.85, 25, 0, 150, 0.05, 1e5,true);  
  plot(Sample::wen, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","wen_met_hadronHF"  ,   0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","wen_met_egammaHF", 0.85, 20, 0, 50, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","wen_met_chargeHadron",  0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","wen_met_neutralHadron", 0.85, 25, 0, 120, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","wen_met_electrons", 0.85, 20, 0, 100, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","wen_met_muons",     0.85, 25, 0, 150, 0.05, 7e6,true);
  plot(Sample::wen, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","wen_met_photons",   0.85, 25, 0, 150, 0.05, 7e6,true);
  */
  //  plot(Sample::gam, "phpt"     , ""       , "p_{T}^{#gamma} [GeV]"       , "Events", "outputPlots","gam_ph1pt" , 0.864, 25, 175, 700, 0.05, 7e6,false);
    //  plot(Sample::gam, "pheta"     , ""       , "#eta_{#gamma}"       , "Events", "outputPlots","gam_ph1eta" , 0.864, 20, -1.5,  1.5, 0.05, 7e6,false);
  //plot(Sample::gam, "nvtx"      , ""       , "N_{PV}"                  , "Events", "outputPlots","gam_nvtx"  , 0.864, 25, 5.5, 30.5, 0.05, 7e6,false);
  //  plot(Sample::gam, "t1phmet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","gam_t1recoil"  , 0.864, 25, 175, 700, 0.05, 7e6,false);
  //  plot(Sample::gam, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","gam_t1met"  ,    0.864, 25, 0, 150, 0.05, 7e6,false);
  //  plot(Sample::gam, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","gam_met",        0.864, 25, 0, 150, 0.05, 1e5,false);  
  /*
  plot(Sample::gam, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","gam_met_hadronHF"  ,   0.864, 25, 0, 150, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","gam_met_egammaHF", 0.864, 20, 0, 50, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","gam_met_chargeHadron",  0.864, 25, 0, 150, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","gam_met_neutralHadron", 0.864, 25, 0, 120, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","gam_met_electrons", 0.864, 20, 0, 100, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","gam_met_muons",     0.864, 25, 0, 150, 0.05, 7e6,false);
  plot(Sample::gam, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","gam_met_photons",   0.864, 25, 0, 150, 0.05, 7e6,false);
  */
}
