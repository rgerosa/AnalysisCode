#include "plot.h"

void showPlots() {
    plot(Sample::gam, "phpt"           , ""         , "Photon p_{T} [GeV]"        , "Events / 30 GeV", "gam_phpt"  , 0.218, 30, 200, 1100, 0.05, 1e5);
    plot(Sample::gam, "t1phmet"        , ""         , "Recoil [GeV]"              , "Events / 30 GeV", "gam_recoil", 0.218, 30, 200, 1100, 0.05, 1e5);
    plot(Sample::gam, "njets"          , ""         , "Number of Jets"            , "Events"         , "gam_njets" , 0.218, 10,   0,   10, 0.05, 1e5);
    plot(Sample::gam, "centraljetpt[0]", ""         , "Leading Jet p_{T} [GeV]"   , "Events / 50 GeV", "gam_jet1pt", 0.218, 20,   0, 1000, 0.05, 1e5);
    plot(Sample::gam, "centraljetpt[1]", "njets > 1", "Subleading Jet p_{T} [GeV]", "Events / 30 GeV", "gam_jet2pt", 0.218, 30,   0,  900, 0.05, 1e5);
}
