#ifndef YIELD_H
#define YIELD_H

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

double yield(const char* path, const char* cut, double lumi) {
    TFile* file = new TFile(path);
    TTree* tree = (TTree*)file->Get("tree/tree");
    tree->Draw("nsignaljets>>hist(1, 0, 10)", cut);
    TH1* hist = (TH1F*)gDirectory->Get("hist");
    hist->Scale(lumi);
    double yld = hist->Integral();
    file->Close();
    return yld;
}

#endif
