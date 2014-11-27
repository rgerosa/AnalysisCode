#ifndef YIELD_H
#define YIELD_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <utility>
#include <string>

double yield(const char* path, const char* cut, double lumi, double rescale = 1.0, const char* treename = "tree/tree") {
    TFile* file = new TFile(path);
    TTree* tree = (TTree*)file->Get(treename);
    if (!tree) tree = (TTree*)file->Get("tree");
    tree->Draw("njets>>hist(1, 0, 1000)", cut);
    TH1* hist = (TH1F*)gDirectory->Get("hist");
    hist->Scale(lumi * rescale);
    double yld = hist->GetBinContent(1);
    file->Close();
    return yld;
}

std::pair<double, double> yieldwitherror(const char* path, const char* cut, double lumi, double rescale = 1.0, const char* treename = "tree/tree", const char* histname = "hist") {
    TFile* file = new TFile(path);
    TTree* tree = (TTree*)file->Get(treename);
    if (!tree) tree = (TTree*)file->Get("tree");
    std::string drawstr = "njets>>";
    drawstr += histname;
    TH1F* hist = new TH1F(histname, "", 1, 0, 1000);
    hist->Sumw2();
    tree->Draw(drawstr.c_str(), cut);
    hist->SetName(histname);
    hist->Scale(lumi * rescale);
    double yld = hist->GetBinContent(1);
    double ylderr = hist->GetBinError(1);
    file->Close();
    return std::pair<double, double>(yld, ylderr);
}


#endif
