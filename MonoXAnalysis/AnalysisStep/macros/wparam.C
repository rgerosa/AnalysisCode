#include <iostream>
#include <sstream>
#include "yield.h"
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

void wparam() {
    double lumi       = 19.7;
    const char* path = "/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/wjets/tree.root";

    std::string numcut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    numcut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nmuons == 0 && nelectrons == 0 && ntaus == 0 && ";
    numcut += "abs(pfmet - calomet) < 2*calomet && mumet > ";

    std::string dencut = "signaljetNHfrac < 0.7 && signaljetEMfrac < 0.7 && signaljetCHfrac > 0.2 && (njets < 2 || (secondjetNHfrac < 0.7 && secondjetEMfrac < 0.9)) && ";
    dencut += "signaljetpt > 110 && abs(signaljeteta) < 2.4 && njets <= 2 && (njets == 1 || abs(jetjetdphi) < 2.5) && nelectrons == 0 && ntaus == 0 && ";
    dencut += "(abs(l1id) == 13 || abs(l2id) == 13) && wmt > 50 && wmt < 100 && abs(mu1eta) < 2.4 && mu1pt > 20 && mu1id == 1 && abs(pfmet - calomet) < 2*calomet && mumet > ";

    float bins[] = {200., 210., 230., 240., 260., 280., 300., 350., 400., 450., 500., 600., 800., 1000.};
    int nbins = 13;

    TH1F* hnum = new TH1F("hnum", "", nbins, bins);
    TH1F* hden = new TH1F("hden", "", nbins, bins);

    hnum->Sumw2();
    hden->Sumw2();

    for (int i = 0; i < nbins; i++) {
        std::stringstream str_den_cut;
        std::stringstream str_num_cut;
        std::stringstream histpas_ss;
        std::stringstream histall_ss;
        
        str_den_cut << "weight * (" << dencut.c_str() << bins[i] << " && mumet < " << bins[i+1] << ")";
        str_num_cut << "weight * (" << numcut.c_str() << bins[i] << " && mumet < " << bins[i+1] << ")";
        histpas_ss  << "histpas" << i;        
        histall_ss  << "histall" << i;        

        std::pair<double, double> val_den = yieldwitherror(path, str_den_cut.str().c_str(), lumi, 1.0, "tree/tree", histall_ss.str().c_str());
        std::pair<double, double> val_rec = yieldwitherror(path, str_num_cut.str().c_str(), lumi, 1.0, "tree/tree", histpas_ss.str().c_str());

        hnum->SetBinContent(i+1, val_rec.first);
        hden->SetBinContent(i+1, val_den.first);

        hnum->SetBinError  (i+1, val_rec.second);
        hden->SetBinError  (i+1, val_den.second);
    }

    hnum->Scale(1, "width");
    hden->Scale(1, "width");

    TF1* fitfunc = new TF1("fitfunc", "[0]/pow([1] + x, [2])", bins[0], bins[nbins]);
    fitfunc->SetParameter(0, 70);
    fitfunc->SetParameter(1, -60);
    fitfunc->SetParameter(2, 1);

    hnum->Divide(hden);
    hnum->Fit(fitfunc, "I");

    fitfunc->SetLineStyle(7);
    fitfunc->SetLineColor(kBlue);
    fitfunc->SetLineWidth(2);
    TH1F* h2 = new TH1F("h2", "h2", nbins, bins);
    for (int i = 0; i <= nbins; i++) {
        h2->SetBinContent(i, hnum->GetBinContent(i));
        h2->SetBinError(i, hnum->GetBinError(i));
    }
    h2->Draw();
    fitfunc->Draw("SAME");
}

