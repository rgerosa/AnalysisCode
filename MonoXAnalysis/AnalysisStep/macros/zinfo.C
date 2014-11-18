#include <iostream>
#include <sstream>
#include "yield.h"

void zinfo() {
    bool useNuNuAcceptance = false;
    double lumi            = 19.7;
    double metcut          = 250.;
    const char* pathz      = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/zjets/tree.root";
    const char* pathd      = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/met/tree.root";
    const char* path1      = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/znunu/tree.root";
    const char* path2      = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgz/tree.root";
    const char* path3      = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgnowz/tree.root";

    std::string nomucut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    nomucut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nelectrons == 0 && ntaus == 0 && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string nucut   = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    nucut   += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nelectrons == 0 && ntaus == 0 && nmuons == 0 && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string chr_de0_cut = "(abs(l1id) == 12 || abs(l1id) == 14 || abs(l1id) == 16) && " + nucut; 
    std::string chr_nu0_cut = "(abs(l1id) == 12 || abs(l1id) == 14 || abs(l1id) == 16) && wzmass > 60 && wzmass < 120 && abs(l1eta) < 2.4 && abs(l2eta) < 2.4 && l1pt > 20 && l2pt > 20 && " + nucut; 
    std::string chr_den_cut = " abs(l1id) == 13 && " + nomucut; 
    std::string chr_num_cut = " abs(l1id) == 13 && wzmass > 60 && wzmass < 120 &&  l1pt > 20 &&  l2pt > 20 && abs(l1eta) < 2.4 && abs(l2eta) < 2.4 && " + nomucut; 
    std::string chr_rec_cut = " abs(l1id) == 13 &&  zmass > 60 &&  zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1) && " + nomucut; 
    std::string chr_yld_cut = "                     zmass > 60 &&  zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1) && " + nomucut; 
    std::string chr_znn_cut = "nmuons == 0 && " + nomucut; 
  
    std::stringstream str_de0_cut;
    std::stringstream str_nu0_cut;
    std::stringstream str_den_cut;
    std::stringstream str_num_cut;
    std::stringstream str_rec_cut;
    std::stringstream str_yld_cut;
    std::stringstream str_znn_cut;

    str_de0_cut << "puwgt * wgt * (" << chr_de0_cut << metcut << ")";
    str_nu0_cut << "puwgt * wgt * (" << chr_nu0_cut << metcut << ")";
    str_den_cut << "puwgt * wgt * (" << chr_den_cut << metcut << " && wzpt > " << metcut - 50. << ")";
    str_num_cut << "puwgt * wgt * (" << chr_num_cut << metcut << " && wzpt > " << metcut - 50. << ")";
    str_rec_cut << "puwgt * wgt * (" << chr_rec_cut << metcut << " && wzpt > " << metcut - 50. << ")";
    str_yld_cut << "puwgt * wgt * (" << chr_yld_cut << metcut << ")";
    str_znn_cut << "puwgt * wgt * (" << chr_znn_cut << metcut << ")";

    double val_de0 = yield(path1, str_de0_cut.str().c_str(), lumi);
    double val_nu0 = yield(path1, str_nu0_cut.str().c_str(), lumi);
    double val_den = yield(pathz, str_den_cut.str().c_str(), lumi);
    double val_num = yield(pathz, str_num_cut.str().c_str(), lumi);
    double val_rec = yield(pathz, str_rec_cut.str().c_str(), lumi);
    double val_znn = yield(path1, str_znn_cut.str().c_str(), lumi);
    double val_ymc = yield(pathz, str_yld_cut.str().c_str(), lumi);
    double val_ydt = yield(pathd, str_yld_cut.str().c_str(), 1.0 );
    double val_ybk = yield(path2, str_yld_cut.str().c_str(), lumi);
           val_ybk+= yield(path3, str_yld_cut.str().c_str(), lumi);

    double Axe = (val_rec / val_den);
    if (useNuNuAcceptance) Axe = (val_rec / val_num) * (val_nu0 / val_de0);

    std::cout << "MET cut                : " << metcut                                              << std::endl;
    std::cout << "Acceptance(nunu)       : " << val_nu0 / val_de0                                   << std::endl;
    std::cout << "Acceptance             : " << val_num / val_den                                   << std::endl;
    std::cout << "Efficiency             : " << val_rec / val_num                                   << std::endl;
    std::cout << "Z->mumu Yield (MC)     : " << val_ymc                                             << std::endl;
    std::cout << "Z->mumu Yield (Data)   : " << val_ydt                                             << std::endl;
    std::cout << "Bkg Yield (MC)         : " << val_ybk                                             << std::endl;
    std::cout << "Yield (MC)             : " << val_znn                                             << std::endl;
    std::cout << "Predicted Yield (MC)   : " << val_ymc * 5.942 / Axe                               << std::endl;
    std::cout << "Predicted Yield (Data) : " << (val_ydt - val_ybk) * 5.942 * 1.023 / Axe           << std::endl;
    std::cout << std::endl                   << std::endl                                           << std::endl;
}

