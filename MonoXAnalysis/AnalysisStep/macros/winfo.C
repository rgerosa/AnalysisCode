#include <iostream>
#include <sstream>
#include "yield.h"

void winfo() {
    double lumi       = 19.7;
    double metcut     = 250.;
    const char* pathw = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/wjets/tree.root";
    const char* pathd = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/met/tree.root";
    const char* path1 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/wjets/tree.root";
    const char* path2 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgw/tree.root";
    const char* path3 = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgnowz/tree.root";
  
    std::string nomucut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    nomucut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nelectrons == 0 && ntaus == 0      && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string noelcut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    noelcut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nmuons == 0     && ntaus == 0      && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string notacut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    notacut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nmuons == 0     && nelectrons == 0 && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string chr_den_cut0 = " (abs(l1id) == 13 || abs(l2id) == 13) && " + nomucut; 
    std::string chr_num_cut0 = "((abs(l1id) == 13 && abs(l1eta) < 2.4 && l1pt > 20) || (abs(l2id) == 13 && abs(l2eta) < 2.4 && l2pt > 20)) && wzmt > 50 && wzmt < 100 && " + nomucut; 
    std::string chr_rec_cut0 = "(abs(l1id) == 13 || abs(l2id) == 13) && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1 && wmt > 50 && wmt < 100 && " + nomucut; 

    std::string chr_den_cut1 = " (abs(l1id) == 13 || abs(l2id) == 13) && " + nomucut; 
    std::string chr_num_cut1 = "((abs(l1id) == 13 && abs(l1eta) < 2.4 && l1pt > 10) || (abs(l2id) == 13 && abs(l2eta) < 2.4 && l2pt > 10)) && " + nomucut; 
    std::string chr_rec_cut1 = "((abs(l1id) == 13 && abs(l1eta) < 2.4 && l1pt > 10) || (abs(l2id) == 13 && abs(l2eta) < 2.4 && l2pt > 10)) && nmuons > 0 && " + nomucut; 

    std::string chr_den_cut2 = " (abs(l1id) == 11 || abs(l2id) == 11) && " + noelcut; 
    std::string chr_num_cut2 = "((abs(l1id) == 11 && abs(l1eta) < 2.5 && l1pt > 10) || (abs(l2id) == 11 && abs(l2eta) < 2.5 && l2pt > 10)) && " + noelcut; 
    std::string chr_rec_cut2 = "((abs(l1id) == 11 && abs(l1eta) < 2.5 && l1pt > 10) || (abs(l2id) == 11 && abs(l2eta) < 2.5 && l2pt > 10)) && nelectrons > 0 && " + noelcut; 

    std::string chr_den_cut3 = " (abs(l1id) == 15 || abs(l2id) == 15) && " + notacut; 
    std::string chr_num_cut3 = "((abs(l1id) == 15 && abs(l1eta) < 2.3 && l1pt > 20) || (abs(l2id) == 15 && abs(l2eta) < 2.3 && l2pt > 20)) && " + notacut; 
    std::string chr_rec_cut3 = "((abs(l1id) == 15 && abs(l1eta) < 2.3 && l1pt > 20) || (abs(l2id) == 15 && abs(l2eta) < 2.3 && l2pt > 20)) && ntaus > 0 && " + notacut; 

    std::string chr_yld_cut  = "wmt > 50 && wmt < 100 && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1 && " + nomucut; 
    std::string chr_wln_cut  = "nmuons == 0 && " + nomucut; 
 
    std::stringstream str_den_cut0;
    std::stringstream str_num_cut0;
    std::stringstream str_rec_cut0;
    std::stringstream str_den_cut1;
    std::stringstream str_num_cut1;
    std::stringstream str_rec_cut1;
    std::stringstream str_den_cut2;
    std::stringstream str_num_cut2;
    std::stringstream str_rec_cut2;
    std::stringstream str_den_cut3;
    std::stringstream str_num_cut3;
    std::stringstream str_rec_cut3;
    std::stringstream str_yld_cut;
    std::stringstream str_wln_cut;

    str_den_cut0 << "puwgt * wgt * (" << chr_den_cut0 << metcut << ")";
    str_num_cut0 << "puwgt * wgt * (" << chr_num_cut0 << metcut << ")";
    str_rec_cut0 << "puwgt * wgt * (" << chr_rec_cut0 << metcut << ")";
    str_den_cut1 << "puwgt * wgt * (" << chr_den_cut1 << metcut << ")";
    str_num_cut1 << "puwgt * wgt * (" << chr_num_cut1 << metcut << ")";
    str_rec_cut1 << "puwgt * wgt * (" << chr_rec_cut1 << metcut << ")";
    str_den_cut2 << "puwgt * wgt * (" << chr_den_cut2 << metcut << ")";
    str_num_cut2 << "puwgt * wgt * (" << chr_num_cut2 << metcut << ")";
    str_rec_cut2 << "puwgt * wgt * (" << chr_rec_cut2 << metcut << ")";
    str_den_cut3 << "puwgt * wgt * (" << chr_den_cut3 << metcut << ")";
    str_num_cut3 << "puwgt * wgt * (" << chr_num_cut3 << metcut << ")";
    str_rec_cut3 << "puwgt * wgt * (" << chr_rec_cut3 << metcut << ")";
    str_yld_cut  << "puwgt * wgt * (" << chr_yld_cut  << metcut << ")";
    str_wln_cut  << "puwgt * wgt * (" << chr_wln_cut  << metcut << ")";

    double val_den0 = yield(pathw, str_den_cut0.str().c_str(), lumi);
    double val_num0 = yield(pathw, str_num_cut0.str().c_str(), lumi);
    double val_rec0 = yield(pathw, str_rec_cut0.str().c_str(), lumi);
    double val_den1 = yield(pathw, str_den_cut1.str().c_str(), lumi);
    double val_num1 = yield(pathw, str_num_cut1.str().c_str(), lumi);
    double val_rec1 = yield(pathw, str_rec_cut1.str().c_str(), lumi);
    double val_den2 = yield(pathw, str_den_cut2.str().c_str(), lumi);
    double val_num2 = yield(pathw, str_num_cut2.str().c_str(), lumi);
    double val_rec2 = yield(pathw, str_rec_cut2.str().c_str(), lumi);
    double val_den3 = yield(pathw, str_den_cut3.str().c_str(), lumi);
    double val_num3 = yield(pathw, str_num_cut3.str().c_str(), lumi);
    double val_rec3 = yield(pathw, str_rec_cut3.str().c_str(), lumi);
    double val_wln  = yield(path1, str_wln_cut .str().c_str(), lumi);
    double val_ymc  = yield(pathw, str_yld_cut .str().c_str(), lumi);
    double val_ydt  = yield(pathd, str_yld_cut .str().c_str(), 1.0 );
    double val_ybk  = yield(path2, str_yld_cut .str().c_str(), lumi);
           val_ybk += yield(path3, str_yld_cut .str().c_str(), lumi);

    double fe = val_den2  / val_den0;
    double ft = val_den3  / val_den0;

    std::cout << "Acceptance             : " << val_num0 / val_den0                   << std::endl;
    std::cout << "Efficiency             : " << val_rec0 / val_num0                   << std::endl;
    std::cout << "Acceptance (mu)        : " << val_num1 / val_den1                   << std::endl;
    std::cout << "Efficiency (mu)        : " << val_rec1 / val_num1                   << std::endl;
    std::cout << "Acceptance (el)        : " << val_num2 / val_den2                   << std::endl;
    std::cout << "Efficiency (el)        : " << val_rec2 / val_num2                   << std::endl;
    std::cout << "Acceptance (ta)        : " << val_num3 / val_den3                   << std::endl;
    std::cout << "Efficiency (ta)        : " << val_rec3 / val_num3                   << std::endl;
    std::cout << "fel                    : " << fe                                    << std::endl;
    std::cout << "fta                    : " << ft                                    << std::endl;
    std::cout << "W->lnu Yield (MC)      : " << val_ymc                               << std::endl;
    std::cout << "W->lnu Yield (Data)    : " << val_ydt                               << std::endl;
    std::cout << "Yield (MC)             : " << val_wln                               << std::endl;
    std::cout << "Predicted Yield (MC)   : " << val_ymc           * ((1. - (val_rec1 / val_den1)) + fe * (1. - (val_rec2 / val_den2)) + ft * (1. - (val_rec3 / val_den3))) / (val_rec0 / val_den0) << std::endl;
    std::cout << "Predicted Yield (Data) : " << (val_ydt-val_ybk) * ((1. - (val_rec1 / val_den1)) + fe * (1. - (val_rec2 / val_den2)) + ft * (1. - (val_rec3 / val_den3))) / (val_rec0 / val_den0) << std::endl;
    std::cout << "Yield Bkg Control      : " << val_ybk           * ((1. - (val_rec1 / val_den1)) + fe * (1. - (val_rec2 / val_den2)) + ft * (1. - (val_rec3 / val_den3))) / (val_rec0 / val_den0) << std::endl;
}
