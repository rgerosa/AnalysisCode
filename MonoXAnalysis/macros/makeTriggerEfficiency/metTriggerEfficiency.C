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

void trigeff() {

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->cd();

    TF1 *fitfunc = new TF1("fitfunc", ErfCB, 100, 1000, 5);
    fitfunc->SetParameters(117., 25., 30., 4., 1.);

    Float_t bins[] = {100., 105., 110., 115., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 500., 600., 1000.};
    Int_t   nbins  = 19;

    TFile* file = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SingleMuon_noMetCut/wmnfilter/wmn_SingleMuon.root");
    TTree* tree = (TTree*)file->Get("tree/tree");

    TH1F* hnum = new TH1F("hnum", "", nbins, bins);
    TH1F* hden = new TH1F("hden", "", nbins, bins);

    tree->Draw("t1mumet>>hnum","          (hltmet90==1 || hltmetwithmu170==1) && hltsinglemu==1 && mu1pt>30 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    
    tree->Draw("t1mumet>>hden","                                                 hltsinglemu==1 && mu1pt>30 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    

    TEfficiency* teff = new TEfficiency(*hnum, *hden);
    teff->SetFillColor(0);

    TEfficiency* deff = (TEfficiency*)teff->Clone("deff1");
    deff->Fit(fitfunc);
    fitfunc->SetLineColor(kBlue);
    fitfunc->SetLineWidth(2);

    TH1* frame = canvas->DrawFrame(bins[0], 0.2, bins[nbins], 1.1, "");
    frame->GetXaxis()->SetTitle("E_{T#mu}^{miss} [GeV]    ");
    frame->GetYaxis()->SetTitle("Trigger Efficiency");
    frame->GetYaxis()->CenterTitle();
    frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
    frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
    frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
    frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
    frame->GetXaxis()->SetTitleOffset(1.0);

    canvas->SetRightMargin(0.075);
    canvas->SetTopMargin(0.06);
    canvas->Draw();
    canvas->cd();
    frame->Draw();
    teff->Draw("SAME");
    fitfunc->Draw("SAME");

    canvas->RedrawAxis();

    CMS_lumi(canvas, 4, 0);
}

