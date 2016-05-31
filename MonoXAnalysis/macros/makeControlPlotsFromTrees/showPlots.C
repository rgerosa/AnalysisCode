#include "plot.h"

void showPlots() {

  gROOT->SetBatch(kTRUE);

  plot(Sample::zmm, "zmass"        , ""       , "M_{#mu#mu} [GeV]"        , "Events", "outputPlots","zmm_zmass" , 0.59, 30, 60, 120, 0.05, 1e6,true);
  plot(Sample::zee, "zeemass"        , ""       , "M_{ee} [GeV]"        , "Events", "outputPlots","zee_zmass" , 0.59, 30, 60, 120, 0.05, 1e6,true);
  plot(Sample::zmm, "nvtx"         , ""       , "N_{PV}"        , "Events", "outputPlots","zmm_nvtx"  , 0.59, 25, 5.5, 30.5, 0.05, 1e6,true);
  plot(Sample::zmm, "zpt"          , ""       , "N_{PV}"        , "Events", "outputPlots","zmm_pt"    , 0.59, 35, 0,  250, 0.05, 1e6,true);
  plot(Sample::zmm, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","zmm_t1recoil"  , 0.59, 25, 0, 400, 0.05, 1e6,true);
  plot(Sample::zmm, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","zmm_t1met"  ,    0.59, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","zmm_met",        0.59, 25, 0, 150, 0.05, 1e5,true);  
  plot(Sample::zmm, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","zmm_met_hadronHF"  ,   0.59, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","zmm_met_egammaHF", 0.59, 20, 0, 50, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","zmm_met_chargeHadron",  0.59, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","zmm_met_neutralHadron", 0.59, 25, 0, 120, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","zmm_met_electrons", 0.59, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","zmm_met_muons",     0.59, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::zmm, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","zmm_met_photons",   0.59, 25, 0, 150, 0.05, 1e6,true);

  /*
    plot(Sample::zee, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","zee_t1met"  ,    0.218, 25, 0, 150, 0.05, 1e5,true);
    plot(Sample::zee, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","zee_met",        0.218, 25, 0, 150, 0.05, 1e5,true);
    plot(Sample::zee, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","zee_met_hadronHF"  ,   0.218, 25, 0, 150, 0.05, 1e5,true);
    plot(Sample::zee, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","zee_met_egammaHF", 0.218, 20, 0, 40, 0.05, 1e5,true);
  plot(Sample::zee, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","zee_met_chargeHadron",  0.218, 25, 0, 150, 0.05, 1e5,true);
  plot(Sample::zee, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","zee_met_neutralHadron", 0.218, 25, 0, 120, 0.05, 1e5,true);
  plot(Sample::zee, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","zee_met_electrons", 0.218, 25, 0, 150, 0.05, 1e5,true);
  plot(Sample::zee, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","zee_met_muons",     0.218, 20, 0, 100, 0.05, 1e5,true);
  plot(Sample::zee, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","zee_met_photons",   0.218, 25, 0, 150, 0.05, 1e5,true);

  plot(Sample::wmn, "nvtx"         , ""       , "N_{PV}"        , "Events", "outputPlots","wmn_nvtx"  , 0.218, 25, 5.5, 30.5, 0.05, 1e6,true);
  plot(Sample::wmn, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","wmn_t1recoil"  , 0.218, 25, 0, 400, 0.05, 1e6,true);
  plot(Sample::wmn, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","wmn_t1met"  ,    0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","wmn_met",        0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","wmn_met_hadronHF"  ,   0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","wmn_met_egammaHF", 0.218, 20, 0, 50, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","wmn_met_chargeHadron",  0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","wmn_met_neutralHadron", 0.218, 20, 0, 125, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","wmn_met_electrons", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","wmn_met_muons",     0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wmn, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","wmn_met_photons",   0.218, 25, 0, 150, 0.05, 1e6,true);

  plot(Sample::wen, "nvtx"         , ""       , "N_{PV}"        , "Events", "outputPlots","wen_nvtx"  , 0.218, 25, 5.5, 30.5, 0.05, 1e6,true);
  plot(Sample::wen, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","wen_t1recoil"  , 0.218, 25, 0, 400, 0.05, 1e6,true);
  plot(Sample::wen, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","wen_t1met"  ,    0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wen, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","wen_met",        0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wen, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","wen_met_hadronHF"  ,   0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","wen_met_egammaHF", 0.218, 20, 0, 50, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","wen_met_chargeHadron",  0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","wen_met_neutralHadron", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","wen_met_electrons", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","wen_met_muons",     0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::wen, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","wen_met_photons",   0.218, 25, 0, 150, 0.05, 1e6,true);


  plot(Sample::topmu, "nvtx"         , ""       , "N_{PV}"        , "Events", "outputPlots","topmu_nvtx"  , 0.218, 25, 5.5, 30.5, 0.05, 1e6,true);
  plot(Sample::topmu, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","topmu_t1recoil"  , 0.218, 25, 0, 400, 0.05, 1e6,true);
  plot(Sample::topmu, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","topmu_t1met"  ,    0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","topmu_met",        0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","topmu_met_hadronHF"  ,   0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","topmu_met_egammaHF", 0.218, 20, 0, 50, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","topmu_met_chargeHadron",  0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","topmu_met_neutralHadron", 0.218, 20, 0, 125, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","topmu_met_electrons", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","topmu_met_muons",     0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topmu, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","topmu_met_photons",   0.218, 25, 0, 150, 0.05, 1e6,true);

  plot(Sample::topel, "nvtx"         , ""       , "N_{PV}"        , "Events", "outputPlots","topel_nvtx"  , 0.218, 25, 5.5, 30.5, 0.05, 1e6,true);
  plot(Sample::topel, "t1mumet"      , ""       , "Recoil [GeV]"        , "Events", "outputPlots","topel_t1recoil"  , 0.218, 25, 0, 400, 0.05, 1e6,true);
  plot(Sample::topel, "t1pfmet"      , ""       , "Met [GeV]"           , "Events", "outputPlots","topel_t1met"  ,    0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topel, "pfmet"        , ""       , "Raw Met [GeV]"       , "Events", "outputPlots","topel_met",        0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topel, "pfmethadronHF", ""       , "Met hadron HF [GeV]" , "Events", "outputPlots","topel_met_hadronHF"  ,   0.218, 25, 0, 150, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetegammaHF", ""       , "Met Egamma HF [GeV]"     , "Events", "outputPlots","topel_met_egammaHF", 0.218, 20, 0, 50, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetchargedhadron", ""  , "Met Charge Hadron [GeV]" , "Events", "outputPlots","topel_met_chargeHadron",  0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetneutralhadron", ""  , "Met Neutral Hadron [GeV]", "Events", "outputPlots","topel_met_neutralHadron", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetelectrons", ""      , "Met Electrons [GeV]"     , "Events", "outputPlots","topel_met_electrons", 0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetmuons"    , ""      , "Met Muons [GeV]"         , "Events", "outputPlots","topel_met_muons",     0.218, 20, 0, 100, 0.05, 1e6,true);
  plot(Sample::topel, "pfmetphotons"  , ""      , "Met Photons [GeV]"       , "Events", "outputPlots","topel_met_photons",   0.218, 25, 0, 150, 0.05, 1e6,true);
  */

}
