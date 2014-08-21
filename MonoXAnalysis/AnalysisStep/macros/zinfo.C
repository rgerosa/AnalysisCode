#include <iostream>
#include <sstream>
#include "yield.h"

void zinfo() {
    double lumi       = 19.7;
    double metcut     = 200.;
    const char* pathz = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/zll/tree.root";
    const char* path1 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/znn50to100/tree.root";
    const char* path2 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/znn100to200/tree.root";
    const char* path3 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/znn200to400/tree.root";
    const char* path4 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/znn400toinf/tree.root";
   
    const char* chr_den_cut = "puwgt * wgt * (abs(l1id) == 13 && abs(l2id) == 13 && met > "; 
    const char* chr_num_cut = "puwgt * wgt * (abs(l1id) == 13 && abs(l2id) == 13 && wzmass > 60 && wzmass < 120 && abs(l1eta) < 2.4 && abs(l2eta) < 2.4 && l1pt > 20 && l2pt > 20 && met > "; 
    const char* chr_rec_cut = "puwgt * wgt * (abs(l1id) == 13 && abs(l2id) == 13 && wzmass > 60 && wzmass < 120 && abs(l1eta) < 2.4 && abs(l2eta) < 2.4 && l1pt > 20 && l2pt > 20 && zmass > 60 && zmass < 120 && met > "; 
    const char* chr_yld_cut = "puwgt * wgt * (zmass > 60 && zmass < 120 && nsignaljets == 1 && signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && abs(jetjetdphi) < 2.5 && nelectrons == 0 && ntaus == 0 && met > "; 
    const char* chr_znn_cut = "puwgt * wgt * (nsignaljets == 1 && signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && abs(jetjetdphi) < 2.5 && nmuons == 0 && nelectrons == 0 && ntaus == 0 && met > "; 
 
    std::stringstream str_den_cut;
    std::stringstream str_num_cut;
    std::stringstream str_rec_cut;
    std::stringstream str_yld_cut;
    std::stringstream str_znn_cut;

    str_den_cut << chr_den_cut << metcut << ")";
    str_num_cut << chr_num_cut << metcut << ")";
    str_rec_cut << chr_rec_cut << metcut << ")";
    str_yld_cut << chr_yld_cut << metcut << ")";
    str_znn_cut << chr_znn_cut << metcut << ")";

    double val_den = yield(pathz, str_den_cut.str().c_str(), lumi);
    double val_num = yield(pathz, str_num_cut.str().c_str(), lumi);
    double val_rec = yield(pathz, str_rec_cut.str().c_str(), lumi);
    double val_yld = yield(pathz, str_yld_cut.str().c_str(), lumi);
    double val_znn = yield(path1, str_znn_cut.str().c_str(), lumi);
           val_znn+= yield(path2, str_znn_cut.str().c_str(), lumi);
           val_znn+= yield(path3, str_znn_cut.str().c_str(), lumi);              
           val_znn+= yield(path4, str_znn_cut.str().c_str(), lumi);              

    std::cout << "Acceptance      : " << val_num / val_den                     << std::endl;
    std::cout << "Efficiency      : " << val_rec / val_num                     << std::endl;
    std::cout << "Z->mumu Yield   : " << val_yld                               << std::endl;
    std::cout << "Predicted Yield : " << val_yld * 5.942 / (val_rec / val_den) << std::endl;
    std::cout << "Actual Yield    : " << val_znn                               << std::endl;
}
