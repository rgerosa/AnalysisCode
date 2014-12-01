#include <iostream>
#include <sstream>
#include "yield.h"

void fillhist(const char* filename, const char* treename, int process, double scale, TH1* hist) {
    TFile* file = new TFile(filename);
    TTree* tree = (TTree*)file->Get(treename);

    TBranch *bmet    = tree->GetBranch("mumet");
    TBranch *bweight = tree->GetBranch("weight");
    TBranch *bnjets  = tree->GetBranch("njets");
    TBranch *bjjdphi = tree->GetBranch("jetjetdphi");

    double vmet     = 0.0;
    double vweight  = 0.0;
    unsigned vnjets = 0;
    double vjjdphi  = 0.0;

    bmet   ->SetAddress(&vmet);
    bweight->SetAddress(&vweight);
    bnjets ->SetAddress(&vnjets);
    bjjdphi->SetAddress(&vjjdphi);

    for (int i = 0; i < tree->GetEntries(); i++) {
        bmet   ->GetEvent(i);
        bweight->GetEvent(i);
        bnjets ->GetEvent(i);
        bjjdphi->GetEvent(i);

        if (process == 0) vweight *= scale;
        if (process == 1) vweight *= scale * (5.942 * 1.023) / (0.79*(1.0-exp(-0.00910276*(vmet-36.1669))));
        if (process == 2) vweight *= scale * 38.7823/pow(-85.7023 + vmet, 0.667232);
        hist->Fill(vmet, vweight);
    }
}

void shapehist() {
    bool makeplot  = true;

    double bins[] = {250., 266, 285., 305., 330., 360., 400., 450., 500., 550., 1000.};
    int nbins = 10;

    TH1F* hisdt = new TH1F("data_obs"                 , "", nbins, bins);
    TH1F* hisdm = new TH1F("DM"                       , "", nbins, bins);
    TH1F* hiszj = new TH1F("ZJets"                    , "", nbins, bins);
    TH1F* hisdi = new TH1F("Dibosons"                 , "", nbins, bins);
    TH1F* histt = new TH1F("ttbar"                    , "", nbins, bins);
    TH1F* hisst = new TH1F("singletop"                , "", nbins, bins);
    TH1F* hisqc = new TH1F("QCD"                      , "", nbins, bins);
    TH1F* hiswj = new TH1F("WJets"                    , "", nbins, bins);
    TH1F* hiszn = new TH1F("Znunu"                    , "", nbins, bins);
    TH1F* hist1 = new TH1F("hist1"                    , "", nbins, bins);
    TH1F* hist2 = new TH1F("hist2"                    , "", nbins, bins);
    TH1F* hiswu = new TH1F("WJets_Stat_WJETSUp"       , "", nbins, bins);
    TH1F* hiswd = new TH1F("WJets_Stat_WJETSDown"     , "", nbins, bins);
    TH1F* hiszu = new TH1F("Znunu_Stat_ZNUNUUp"       , "", nbins, bins);
    TH1F* hiszd = new TH1F("Znunu_Stat_ZNUNUDown"     , "", nbins, bins);

    hisdt->Sumw2();
    hisdm->Sumw2();
    hiszn->Sumw2();
    hiswj->Sumw2();
    hiszj->Sumw2();
    hisdi->Sumw2();
    histt->Sumw2();
    hisst->Sumw2();
    hisqc->Sumw2();
    hiszu->Sumw2();
    hiszd->Sumw2();
    hiswu->Sumw2();
    hiswd->Sumw2();
    hist1->Sumw2();
    hist2->Sumw2();

    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/met/ztree.root",    "tree/tree", 1, 1.0 , hiszn);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgz/ztree.root",   "tree/tree", 1, 19.7, hist1);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgnowz/ztree.root","tree/tree", 1, 19.7, hist1);

    hiszn->Add(hist1, -1.0);

    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/met/wtree.root",    "tree/tree", 2, 1.0 , hiswj);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgw/wtree.root",   "tree/tree", 2, 19.7, hist2);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/bkgnowz/wtree.root","tree/tree", 2, 19.7, hist2);

    hiswj->Add(hist2, -1.0);

    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/dmVM100/reducedtree.root"  ,  "tree/tree", 0, 19.7     , hisdm);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/met/reducedtree.root"      ,  "tree/tree", 0, 1.0      , hisdt);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/zjets/reducedtree.root"    ,  "tree/tree", 0, 19.7     , hiszj);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/dibosons/reducedtree.root" ,  "tree/tree", 0, 19.7     , hisdi);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/ttbar/reducedtree.root"    ,  "tree/tree", 0, 19.7     , histt);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/singletop/reducedtree.root",  "tree/tree", 0, 19.7     , hisst);
    fillhist("/home/avartak/CMS/MonoX/CMSSW_5_3_12/src/MonoXAnalysis/AnalysisStep/trees/qcd/reducedtree.root"      ,  "tree/tree", 0, 19.7*1.3 , hisqc);

    for (int i = 1; i <= nbins; i++) {
        hiszu->SetBinContent(i, hiszn->GetBinContent(i) + hiszn->GetBinError(i));
        hiszd->SetBinContent(i, hiszn->GetBinContent(i) - hiszn->GetBinError(i));
        hiswu->SetBinContent(i, hiswj->GetBinContent(i) + hiswj->GetBinError(i));
        hiswd->SetBinContent(i, hiswj->GetBinContent(i) - hiswj->GetBinError(i));
    }

    for (int i = 1; i <= nbins; i++) {
        if (hiszj->GetBinContent(i) == 0) hiszj->SetBinContent(i, 0.001);
        if (hisdi->GetBinContent(i) == 0) hisdi->SetBinContent(i, 0.001);
        if (histt->GetBinContent(i) == 0) histt->SetBinContent(i, 0.001);
        if (hisst->GetBinContent(i) == 0) hisst->SetBinContent(i, 0.001);
        if (hisqc->GetBinContent(i) == 0) hisqc->SetBinContent(i, 0.001);

        hiszj->SetBinError(i, 0.5*hiszj->GetBinContent(i));
        hisdi->SetBinError(i, 0.5*hisdi->GetBinContent(i));
        histt->SetBinError(i, 0.5*histt->GetBinContent(i));
        hisst->SetBinError(i, 0.5*hisst->GetBinContent(i));
        hisqc->SetBinError(i, 0.5*hisqc->GetBinContent(i));
    }


    std::cout << "data_obs  : " << hisdt->Integral() << std::endl;
    std::cout << "DM        : " << hisdm->Integral() << std::endl;
    std::cout << "Znunu     : " << hiszn->Integral() << std::endl;
    std::cout << "WJets     : " << hiswj->Integral() << std::endl;
    std::cout << "ZJets     : " << hiszj->Integral() << std::endl;
    std::cout << "Dibosons  : " << hisdi->Integral() << std::endl;
    std::cout << "ttbar     : " << histt->Integral() << std::endl;
    std::cout << "singletop : " << hisst->Integral() << std::endl;
    std::cout << "QCD       : " << hisqc->Integral() << std::endl;

    if (makeplot) {
    hiszj->Add(hisdi);
    hiszj->Add(histt);
    hiszj->Add(hisst);
    hiszj->Add(hisqc);

    hiswj->Add(hiszj);
    hiszn->Add(hiswj);
    hiszu->Add(hiswj);
    hiszd->Add(hiswj);

    hiszj->SetLineColor(kGreen+2);
    hiswj->SetLineColor(kRed);
    hiszn->SetLineColor(kBlue);
    hiszj->SetLineWidth(2);
    hiswj->SetLineWidth(2);
    hiszn->SetLineWidth(2);

    hiszu->SetMarkerSize(0);
    hiszu->SetMarkerColor(0);
    hiszu->SetLineColor(kOrange);
    hiszu->SetFillColor(kOrange);

    hiszd->SetMarkerSize(0);
    hiszd->SetMarkerColor(0);
    hiszd->SetLineColor(10);
    hiszd->SetFillColor(10);

    hiszj->SetMarkerSize(0);
    hiswj->SetMarkerSize(0);
    hiszn->SetMarkerSize(0);

    hisdm->SetMarkerSize(0);
    hisdm->SetLineStyle(7);
    hisdm->SetLineWidth(2);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
    canvas->SetRightMargin(0.075);
    canvas->SetTopMargin(0.075);
    TH1* frame = canvas->DrawFrame(250., 2.0, 1000., 100000.0, "");
    frame->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    frame->GetXaxis()->SetLabelSize(0.9*frame->GetXaxis()->GetLabelSize());
    frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
    canvas->SetLogy();
    frame->Draw();
    hiszu->Draw("HIST SAME");
    hiszd->Draw("HIST SAME");
    hiszn->Draw("HIST SAME");
    hiswj->Draw("HIST SAME");
    hiszj->Draw("HIST SAME");
    hisdm->Draw("HIST SAME");
    hisdt->Draw("PE SAME");

    canvas->RedrawAxis();
    }

    else {
    TFile* outfile = new TFile("monojet_8TeV_shape_histogram.root", "RECREATE");
    hisdt->Write();
    hisdm->Write();
    hiszn->Write();
    hiswj->Write();
    hiszj->Write();
    hisdi->Write();
    histt->Write();
    hisst->Write();
    hisqc->Write();
    hiszu->Write();
    hiszd->Write();
    hiswu->Write();
    hiswd->Write();
    outfile->Close();
    }
}

