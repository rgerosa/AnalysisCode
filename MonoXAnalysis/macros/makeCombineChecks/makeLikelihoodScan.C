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

void makeLikelihoodScan(string inputDirectory, string nameToGrep, string outputPlots, string postfix){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputPlots).c_str());

  system(("ls "+inputDirectory+" | grep root | grep "+nameToGrep+" > list.temp ").c_str());
  vector<TFile*> inputFiles;

  ifstream infile ("list.temp");
  if(infile.is_open()){
    string line;
    while(!infile.eof()){
      getline(infile,line);
      if(TString(line).Contains(".root"))
	inputFiles.push_back(TFile::Open((inputDirectory+"/"+line).c_str(),"READ"));
    }
  }
  infile.close();
  system("rm list.temp");

  int nbadfit_bonly = 0;
  int invalid_bonly = 0;
  int nbadfit_sb = 0;
  int invalid_sb = 0;

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();

  TGraph* likelihoodScan = new TGraph();
  TGraph* likelihoodScanAbove = new TGraph();
  TGraph* likelihoodScanBelow = new TGraph();

  TF1* line = new TF1("line","1",-10,10);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  TGraph* positiveMu  = new TGraph();
  TGraph* negativeMu  = new TGraph();
  

  int ifile = 0;
  for(auto file: inputFiles){
    if(ifile >= 1 )continue;
    TTree* limit = (TTree*) file->Get("limit");
    double mh;
    float r, deltaNLL;
    limit->SetBranchAddress("mh",&mh);
    limit->SetBranchAddress("r",&r);
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
    likelihoodScan->GetXaxis()->SetTitle("signal strenght");
    likelihoodScan->GetYaxis()->SetTitle("-2*#Delta Log(L)");
    line->SetRange(point.front().r_,point.back().r_);
    
    likelihoodScan->Draw("APL");
    line->Draw("Lsame");
    CMS_lumi(canvas,"35.9");

    cout<<"best mu "<<muref<<" 1sigma-up "<<fabs(muref-likelihoodScanAbove->Eval(1))<<" 1 sigma-dw "<<fabs(muref-likelihoodScanBelow->Eval(1))<<endl;

    canvas->SaveAs((outputPlots+"/scan_"+to_string(mh)+"_"+nameToGrep+"_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputPlots+"/scan_"+to_string(mh)+"_"+nameToGrep+"_"+postfix+".pdf").c_str(),"pdf");
  }
  
}
