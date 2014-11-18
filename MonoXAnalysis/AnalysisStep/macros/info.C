#include <iostream>
#include <sstream>
#include "yield.h"

void info() {
    bool isData      = false;
    double lumi      = isData ? 1.0 : 19.7;
    //double rescale   = 1./0.995;
    double rescale   = 1.0;
    const char* path = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/dmAVM10/reducedtree.root";

    std::stringstream str_cut;

    /*
    str_cut << (!isData ? "puwgt * wgt * (1 " : "(1 "); 
    str_cut << " && (hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
    str_cut << " && signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9))";
    str_cut << " && signaljetpt > 110 && abs(signaljeteta) < 2.4";
    str_cut << " && njets <= 2";
    str_cut << " && (njets == 1 || abs(jetjetdphi) < 2.5)";
    str_cut << " && nmuons == 0";
    str_cut << " && nelectrons == 0";
    str_cut << " && ntaus == 0";
    str_cut << " && abs(pfmet - calomet) < 2*calomet && mumet > 500";
    str_cut << ")";

    std::pair<double, double> yld  = yieldwitherror(path, str_cut.str().c_str(), lumi, rescale);
    std::cout << "Yield : " << yld.first  << " +/- " << yld.second << std::endl;
    */

    str_cut << (!isData ? "puwgt * wgt * " : "") << "((hltmet120 > 0 || hltmet95jet80 > 0 || hltmet105jet80 > 0)";
    std::string cut1 = str_cut.str() + ")";
    str_cut << " && abs(pfmet - calomet) < 2*calomet && mumet > 200";
    std::string cut2 = str_cut.str() + ")";
    str_cut << " && signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9))";
    std::string cut3 = str_cut.str() + ")";
    str_cut << " && signaljetpt > 110 && abs(signaljeteta) < 2.4";
    std::string cut4 = str_cut.str() + ")";
    str_cut << " && njets <= 2";
    std::string cut5 = str_cut.str() + ")";
    str_cut << " && (njets == 1 || abs(jetjetdphi) < 2.5)";
    std::string cut6 = str_cut.str() + ")";
    str_cut << " && nmuons == 0";
    std::string cut7 = str_cut.str() + ")";
    str_cut << " && nelectrons == 0";
    std::string cut8 = str_cut.str() + ")";
    str_cut << " && ntaus == 0";
    std::string cut9 = str_cut.str() + ")";
    str_cut << " && mumet > 250";
    std::string cut10 = str_cut.str() + ")";
    str_cut << " && mumet > 300";
    std::string cut11 = str_cut.str() + ")";
    str_cut << " && mumet > 350";
    std::string cut12 = str_cut.str() + ")";
    str_cut << " && mumet > 400";
    std::string cut13 = str_cut.str() + ")";
    str_cut << " && mumet > 450";
    std::string cut14 = str_cut.str() + ")";
    str_cut << " && mumet > 500";
    std::string cut15 = str_cut.str() + ")";
    str_cut << " && mumet > 550";
    std::string cut16 = str_cut.str() + ")";

    std::pair<double, double> yld2  = yieldwitherror(path, cut2 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut2    : " << yld2 .first  << " +/- " << yld2 .second << std::endl;
    std::pair<double, double> yld3  = yieldwitherror(path, cut3 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut3    : " << yld3 .first  << " +/- " << yld3 .second << std::endl;
    std::pair<double, double> yld4  = yieldwitherror(path, cut4 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut4    : " << yld4 .first  << " +/- " << yld4 .second << std::endl;
    std::pair<double, double> yld5  = yieldwitherror(path, cut5 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut5    : " << yld5 .first  << " +/- " << yld5 .second << std::endl;
    std::pair<double, double> yld6  = yieldwitherror(path, cut6 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut6    : " << yld6 .first  << " +/- " << yld6 .second << std::endl;
    std::pair<double, double> yld7  = yieldwitherror(path, cut7 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut7    : " << yld7 .first  << " +/- " << yld7 .second << std::endl;
    std::pair<double, double> yld8  = yieldwitherror(path, cut8 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut8    : " << yld8 .first  << " +/- " << yld8 .second << std::endl;
    std::pair<double, double> yld9  = yieldwitherror(path, cut9 .c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut9    : " << yld9 .first  << " +/- " << yld9 .second << std::endl;
    std::pair<double, double> yld10 = yieldwitherror(path, cut10.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut10   : " << yld10.first  << " +/- " << yld10.second << std::endl;
    std::pair<double, double> yld11 = yieldwitherror(path, cut11.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut11   : " << yld11.first  << " +/- " << yld11.second << std::endl;
    std::pair<double, double> yld12 = yieldwitherror(path, cut12.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut12   : " << yld12.first  << " +/- " << yld12.second << std::endl;
    std::pair<double, double> yld13 = yieldwitherror(path, cut13.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut13   : " << yld13.first  << " +/- " << yld13.second << std::endl;
    std::pair<double, double> yld14 = yieldwitherror(path, cut14.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut14   : " << yld14.first  << " +/- " << yld14.second << std::endl;
    std::pair<double, double> yld15 = yieldwitherror(path, cut15.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut15   : " << yld15.first  << " +/- " << yld15.second << std::endl;
    std::pair<double, double> yld16 = yieldwitherror(path, cut16.c_str(), lumi, rescale, "tree/tree");
    std::cout << "Cut16   : " << yld16.first  << " +/- " << yld16.second << std::endl;
}
