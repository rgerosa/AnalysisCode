#include <iostream>
#include <sstream>
#include "yield.h"

void zparam() {
    double lumi       = 19.7;
    const char* pathz = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/zjets/tree.root";

    std::string nomucut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    nomucut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nelectrons == 0 && ntaus == 0 && abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string chr_den_cut = " abs(l1id) == 13 && " + nomucut; 
    std::string chr_rec_cut = " abs(l1id) == 13 &&  zmass > 60 &&  zmass < 120 && mu1pid == -mu2pid && mu1pt > 20 && mu2pt > 20 && (mu1id == 1 || mu2id == 1) && " + nomucut; 

    float bins[] = {200., 220., 240., 260., 300., 360., 460., 640., 1000.};
    int nbins = 8;
    TH1F* hnum = new TH1F("hnum", "", nbins, bins);
    TH1F* hden = new TH1F("hden", "", nbins, bins);

    for (int i = 0; i < nbins; i++) {
        std::stringstream str_den_cut;
        std::stringstream str_rec_cut;
        std::stringstream histnum_ss;
        std::stringstream histden_ss;
        
        str_den_cut << "puwgt * wgt * (" << chr_den_cut << bins[i] << " && mumet < " << bins[i+1] << " && wzpt > " << bins[i] - 50. << ")";
        str_rec_cut << "puwgt * wgt * (" << chr_rec_cut << bins[i] << " && mumet < " << bins[i+1] << " && wzpt > " << bins[i] - 50. << ")";
        histnum_ss  << "histnum" << i;        
        histden_ss  << "histden" << i;        

        std::pair<double, double> val_den = yieldwitherror(pathz, str_den_cut.str().c_str(), lumi, 1.0, "tree/tree", histden_ss.str().c_str());
        std::pair<double, double> val_rec = yieldwitherror(pathz, str_rec_cut.str().c_str(), lumi, 1.0, "tree/tree", histnum_ss.str().c_str());

        hnum->SetBinContent(i+1, val_rec.first);
        hden->SetBinContent(i+1, val_den.first);

        hnum->SetBinError  (i+1, val_rec.second);
        hden->SetBinError  (i+1, val_den.second);
    }

    hnum->Sumw2();
    hden->Sumw2();
    
    hnum->Scale(1, "width");
    hden->Scale(1, "width");

    TEfficiency* teff = new TEfficiency(*hnum, *hden);
    TEfficiency* deff = new TEfficiency(*hnum, *hden);

    TF1* fitfunc = new TF1("fitfunc", "0.79*(1.0-exp([0]*(x-[1])))", 200., 1000.);
    //TF1* fitfunc = new TF1("fitfunc", "[2]*(1.0-exp([0]*(x-[1])))", 200., 1000.);
    fitfunc->SetParameter(0, -0.01);
    fitfunc->SetParameter(1, 30);
    //fitfunc->SetParameter(2, 0.8);
    fitfunc->SetLineColor(kBlack);
    teff->Fit(fitfunc);
    fitfunc->SetLineColor(kBlue);
    fitfunc->SetLineStyle(2);
    fitfunc->SetLineWidth(2);

    deff->Draw();
    fitfunc->Draw("SAME");
}

