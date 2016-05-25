#include <iostream>
#include <sstream>
#include <cmath>

Double_t ApproxErf(Double_t arg) {
    static const double erflim = 5.0;
    if( arg >  erflim) return 1.0;
    if( arg < -erflim) return -1.0;
    
    return TMath::Erf(arg);
}

Double_t ErfCB(double *x, double *par) { 
    double m = x[0];
    double m0 = par[0];
    double sigma = par[1];
    double alpha = par[2];
    double n = par[3];
    double norm = par[4];
    
    const double sqrtPiOver2 = 1.2533141373; // sqrt(pi/2)
    const double sqrt2 = 1.4142135624;
    
    Double_t sig = fabs((Double_t) sigma);
    Double_t t = (m - m0)/sig ;
    
    if (alpha < 0) t = -t;
    
    Double_t absAlpha = fabs(alpha / sig);
    Double_t a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t b = absAlpha - n/absAlpha;
    
    Double_t leftArea = (1 + ApproxErf( absAlpha / sqrt2 )) * sqrtPiOver2 ;
    Double_t rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
    Double_t area = leftArea + rightArea ;
    
    if ( t <= absAlpha )return norm * (1 + ApproxErf( t / sqrt2 )) * sqrtPiOver2 / area ;
    else return norm * (leftArea +  a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area ;
    
} 

TLatex *CMSPreliminary() {

    TLatex* CP = new TLatex(0.18,0.96, "#font[22]{CMS Preliminary        #sqrt{s} = 13 TeV, L = 0.21 fb^{-1}}");
    CP->SetNDC(kTRUE);
    CP->SetTextSize(0.040);

    return CP;

}

void trigeff() {

    TLatex *CP = CMSPreliminary();
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->cd();

    TF1 *fitfunc = new TF1("fitfunc", ErfCB, 100, 1000, 5);
    fitfunc->SetParameters(1., 100., 0.1, 1., 1.);

    Float_t bins[] = {100., 105., 110., 115., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 500., 600., 1000.};
    Int_t   nbins  = 19;

    TFile* file1 = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SingleMuon_noMetCut/wmnfilter/wmn_SingleMuon.root");
    //TFile* file2 = new TFile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/wmntree.root");
    TTree* tree1 = (TTree*)file1->Get("tree/tree");
    //TTree* tree2 = (TTree*)file2->Get("tree/tree");

    TH1F* hnum1 = new TH1F("hnum1", "", nbins, bins);
    TH1F* hden1 = new TH1F("hden1", "", nbins, bins);

    //TH1F* hnum2 = new TH1F("hnum2", "", nbins, bins);
    //TH1F* hden2 = new TH1F("hden2", "", nbins, bins);

    tree1->Draw("t1mumet>>hnum1","          (hltmet90==1 || hltmetwithmu170==1) && hltsinglemu==1 && mu1pt>30 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    
    tree1->Draw("t1mumet>>hden1","                                                 hltsinglemu==1 && mu1pt>30 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    

    //tree2->Draw("t1mumet>>hnum2","xsec*wgt*(hltmet90==1 && flaghbheloose==1 && flagcsctight==1 && hltsinglemu==1 && mu1pt>30 && njets>0 && signaljetCHfrac>0.1)/wgtsum");    
    //tree2->Draw("t1mumet>>hden2","xsec*wgt*(                                       flaghbheloose==1 && flagcsctight==1 && hltsinglemu==1 && mu1pt>30 && njets>0 && signaljetCHfrac>0.1)/wgtsum");    

    TEfficiency* teff1 = new TEfficiency(*hnum1, *hden1);
    //TEfficiency* teff2 = new TEfficiency(*hnum2, *hden2);

    //teff1->Fit(fitfunc);

    //teff2->SetMarkerColor(kRed);
    //teff2->SetLineColor(kRed);

    TH1* frame = canvas->DrawFrame(bins[0], 0.2, bins[nbins], 1.1, "");
    frame->GetXaxis()->SetTitle("E_{T#mu}^{miss} [GeV]    ");
    frame->GetYaxis()->SetTitle("Trigger Efficiency");
    frame->GetYaxis()->CenterTitle();
    frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
    frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());

    canvas->SetRightMargin(0.075);
    canvas->SetTopMargin(0.06);
    canvas->Draw();
    canvas->cd();
    frame->Draw();
    teff1->Draw("SAME");
    //teff2->Draw("SAME");

    teff1->SetFillColor(0);
    //teff2->SetFillColor(0);

    TLegend* leg = new TLegend(0.6, 0.55, 0.9, 0.80);
    leg->SetFillColor(0);
    leg->AddEntry(teff1, "Data");
    //leg->AddEntry(teff2, "W(l#nu) MC");
    //leg->Draw("SAME");

    canvas->RedrawAxis();

    fitfunc->Draw();

    CP->Draw("SAME");

}

