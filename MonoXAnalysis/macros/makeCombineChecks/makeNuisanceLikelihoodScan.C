#include "../CMS_lumi.h"

class likelihoodPoint {

public:
  
  likelihoodPoint(float r, float nll){
    r_ = r;
    nll_ = nll;
  }
  bool operator < (const likelihoodPoint & a) const{
    if(r_ < a.r_) return true;
    else return false;
  }
  float r_;
  float nll_;
};

void makeNuisanceLikelihoodScan(string inputFileName, string nameNuisance, string outputPlots){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputPlots).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* limit = (TTree*) inputFile->Get("limit");

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();

  TGraph* likelihoodScan = new TGraph();
  TGraph* likelihoodScanAbove = new TGraph();
  TGraph* likelihoodScanBelow = new TGraph();

  TF1* line = new TF1("line","1",-10,10);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  float r, deltaNLL;
  limit->SetBranchAddress(nameNuisance.c_str(),&r);
  limit->SetBranchAddress("deltaNLL",&deltaNLL);
  vector<likelihoodPoint> point;
  for(int ientry = 0; ientry < limit->GetEntries(); ientry++){
    limit->GetEntry(ientry);
    point.push_back(likelihoodPoint(r,2*deltaNLL));      
  }
  sort(point.begin(),point.end());

  int ipoint = 0;
  likelihoodScan->Set(0);
  for(auto p: point){
    if(p.nll_ > 50) continue;
    likelihoodScan->SetPoint(ipoint,p.r_,p.nll_);
    ipoint++;
  }


  limit->GetEntry(0);
  float muref = r;

  vector<likelihoodPoint> pointBelow;
  vector<likelihoodPoint> pointAbove;
  for(int ientry = 0; ientry < limit->GetEntries(); ientry++){
    limit->GetEntry(ientry);
    if(r < muref)
      pointBelow.push_back(likelihoodPoint(2*deltaNLL,r));
    else
      pointAbove.push_back(likelihoodPoint(2*deltaNLL,r));    
  }

  sort(pointBelow.begin(),pointBelow.end());
  sort(pointAbove.begin(),pointAbove.end());

  ipoint = 0;
  likelihoodScanAbove->Set(0);
  for(auto p: pointAbove){
    if(p.nll_ > 50) continue;
    likelihoodScanAbove->SetPoint(ipoint,p.r_,p.nll_);
    ipoint++;
  }
  ipoint = 0;
  likelihoodScanBelow->Set(0);
  for(auto p: pointBelow){
    if(p.nll_ > 50) continue;
    likelihoodScanBelow->SetPoint(ipoint,p.r_,p.nll_);
    ipoint++;
  }
  
  
  likelihoodScan->SetLineColor(kBlack);
  likelihoodScan->SetLineWidth(2);
  likelihoodScan->SetMarkerColor(kBlack);
  likelihoodScan->SetMarkerStyle(20);
  likelihoodScan->SetMarkerSize(1);
  likelihoodScan->GetXaxis()->SetTitle(nameNuisance.c_str());
  likelihoodScan->GetYaxis()->SetTitle("-2*#Delta Log(L)");
  likelihoodScan->SetMinimum(0);
  line->SetRange(point.front().r_,point.back().r_);
  
  likelihoodScan->Draw("APL");
  line->Draw("Lsame");
  CMS_lumi(canvas,"35.9");

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.4,0.8,Form("#theta = %.2f + (%.2f) - (%.2f)",muref,fabs(muref-likelihoodScanAbove->Eval(1)),fabs(muref-likelihoodScanBelow->Eval(1))));

  canvas->SaveAs((outputPlots+"/scan_"+nameNuisance+".png").c_str(),"png");
  canvas->SaveAs((outputPlots+"/scan_"+nameNuisance+".pdf").c_str(),"pdf");
}
