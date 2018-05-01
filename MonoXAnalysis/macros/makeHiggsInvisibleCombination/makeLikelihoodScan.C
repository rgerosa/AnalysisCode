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

//// code to plot likelihood scans
void makeLikelihoodScan(string outputPlots, string postfix){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputPlots).c_str());

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();

  TFile* file_scan_monoj_obs = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoJ/higgsCombine_scan_mu_obs.MultiDimFit.mH120.root","READ");
  TFile* file_scan_monov_obs = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoV/higgsCombine_scan_mu_obs.MultiDimFit.mH120.root","READ");
  TFile* file_scan_monoz_obs = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-052/higgsCombine_scan_mu_obs.MultiDimFit.mH120.root","READ");
  TFile* file_scan_vbf_obs   = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/HIG-17-023/higgsCombine_scan_mu_obs.MultiDimFit.mH125.root","READ");
  TFile* file_scan_combined_obs = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/Combination/higgsCombine_scan_mu_obs.MultiDimFit.mH125.root","READ");

  TFile* file_scan_monoj_exp = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoJ/higgsCombine_scan_mu_exp.MultiDimFit.mH120.root","READ");
  TFile* file_scan_monov_exp = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoV/higgsCombine_scan_mu_exp.MultiDimFit.mH120.root","READ");
  TFile* file_scan_monoz_exp = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-052/higgsCombine_scan_mu_exp.MultiDimFit.mH120.root","READ");
  TFile* file_scan_vbf_exp   = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/HIG-17-023/higgsCombine_scan_mu_exp.MultiDimFit.mH120.root","READ");
  TFile* file_scan_combined_exp = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/Combination/higgsCombine_scan_mu_exp.MultiDimFit.mH125.root","READ");

  //////
  TGraph* scan_monoj_obs = new TGraph();
  TGraph* scan_monoz_obs = new TGraph();
  TGraph* scan_monov_obs = new TGraph();
  TGraph* scan_vbf_obs = new TGraph();
  TGraph* scan_combined_obs = new TGraph();

  TGraph* scan_monoj_exp = new TGraph();
  TGraph* scan_monoz_exp = new TGraph();
  TGraph* scan_monov_exp = new TGraph();
  TGraph* scan_vbf_exp = new TGraph();
  TGraph* scan_combined_exp = new TGraph();

  /////
  TTree* limit_scan_monoj_obs = (TTree*) file_scan_monoj_obs->Get("limit");
  TTree* limit_scan_monoz_obs = (TTree*) file_scan_monoz_obs->Get("limit");
  TTree* limit_scan_monov_obs = (TTree*) file_scan_monov_obs->Get("limit");
  TTree* limit_scan_vbf_obs = (TTree*) file_scan_vbf_obs->Get("limit");
  TTree* limit_scan_combined_obs = (TTree*) file_scan_combined_obs->Get("limit");

  TTree* limit_scan_monoj_exp = (TTree*) file_scan_monoj_exp->Get("limit");
  TTree* limit_scan_monoz_exp = (TTree*) file_scan_monoz_exp->Get("limit");
  TTree* limit_scan_monov_exp = (TTree*) file_scan_monov_exp->Get("limit");
  TTree* limit_scan_vbf_exp = (TTree*) file_scan_vbf_exp->Get("limit");
  TTree* limit_scan_combined_exp = (TTree*) file_scan_combined_exp->Get("limit");

  TTreeReader reader(limit_scan_monoj_obs);
  TTreeReaderValue<float> r (reader,"r");
  TTreeReaderValue<float> deltaNLL (reader,"deltaNLL");

  int ipoint = 0;
  float min  = 0;
  
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  reader.SetEntry(0);
  
  while(reader.Next()){
    if(*r < 0) continue;
    scan_monoj_obs->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monoj_obs->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_monoz_obs);
  reader.SetEntry(0);

  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  reader.SetEntry(0);

  while(reader.Next()){
    if(*r < 0) continue;
    scan_monoz_obs->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monoz_obs->Sort();
  ipoint = 0;
  min    = 0;

  ////
  reader.SetTree(limit_scan_monov_obs);
  reader.SetEntry(0);

  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  reader.SetEntry(0);

  while(reader.Next()){
    if(*r < 0) continue;
    scan_monov_obs->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monov_obs->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_vbf_obs);
  reader.SetEntry(0);

  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_vbf_obs->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_vbf_obs->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_combined_obs);
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_combined_obs->SetPoint(ipoint,*r,2*(*deltaNLL));
    ipoint++;
  }
  scan_combined_obs->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_monoj_exp);

  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_monoj_exp->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monoj_exp->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_monoz_exp);

  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }

  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_monoz_exp->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monoz_exp->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_monov_exp);

  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }

  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_monov_exp->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_monov_exp->Sort();
  ipoint = 0;
  min = 0;

  ////
  reader.SetTree(limit_scan_vbf_exp);
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_vbf_exp->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_vbf_exp->Sort();
  ipoint = 0;
  min = 0;
  ////
  reader.SetTree(limit_scan_combined_exp);
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    if(2*(*deltaNLL) < min) min = 2*(*deltaNLL);
  }
  reader.SetEntry(0);
  while(reader.Next()){
    if(*r < 0) continue;
    scan_combined_exp->SetPoint(ipoint,*r,2*(*deltaNLL)+fabs(min));
    ipoint++;
  }
  scan_combined_exp->Sort();

  //// Produce the final plot
  scan_combined_exp->GetXaxis()->SetTitle("BR(H #rightarrow inv)");
  scan_combined_exp->GetXaxis()->SetTitleOffset(1.1);
  scan_combined_exp->GetYaxis()->SetTitle("-2 #Delta Log(L)");
  scan_combined_exp->GetYaxis()->SetTitleOffset(1.1);
  scan_combined_exp->GetYaxis()->SetRangeUser(0,8);
  scan_combined_exp->GetXaxis()->SetRangeUser(0,1);
  scan_combined_exp->SetLineColor(kBlack);
  scan_combined_exp->SetLineWidth(2);
  scan_combined_exp->SetLineStyle(7);
  scan_combined_obs->SetLineColor(kBlack);
  scan_combined_obs->SetLineWidth(2);
  scan_combined_exp->Draw("AL");
  scan_combined_obs->Draw("Lsame");

  scan_monoj_obs->SetLineColor(kOrange);
  scan_monoj_obs->SetLineWidth(2);
  scan_monoj_obs->Draw("Lsame");
  scan_monoj_exp->SetLineColor(kOrange);
  scan_monoj_exp->SetLineWidth(2);
  scan_monoj_exp->SetLineStyle(7);
  scan_monoj_exp->Draw("Lsame");

  scan_monov_obs->SetLineColor(kGreen+2);
  scan_monov_obs->SetLineWidth(2);
  scan_monov_obs->Draw("Lsame");
  scan_monov_exp->SetLineColor(kGreen+2);
  scan_monov_exp->SetLineWidth(2);
  scan_monov_exp->SetLineStyle(7);
  scan_monov_exp->Draw("Lsame");

  scan_monoz_obs->SetLineColor(kBlue);
  scan_monoz_obs->SetLineWidth(2);
  scan_monoz_obs->Draw("Lsame");
  scan_monoz_exp->SetLineColor(kBlue);
  scan_monoz_exp->SetLineWidth(2);
  scan_monoz_exp->SetLineStyle(7);
  scan_monoz_exp->Draw("Lsame");

  scan_vbf_obs->SetLineColor(kRed);
  scan_vbf_obs->SetLineWidth(2);
  scan_vbf_obs->Draw("Lsame");
  scan_vbf_exp->SetLineColor(kRed);
  scan_vbf_exp->SetLineWidth(2);
  scan_vbf_exp->SetLineStyle(7);
  scan_vbf_exp->Draw("Lsame");

  scan_combined_exp->Draw("Lsame");
  scan_combined_obs->Draw("Lsame");


  TLegend leg (0.67,0.17,0.93,0.40);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(scan_combined_obs,"Combined","L");
  leg.AddEntry(scan_vbf_obs,"qqH-tagged","L");
  leg.AddEntry(scan_monoz_obs,"Z(ll)H-tagged","L");
  leg.AddEntry(scan_monov_obs,"V(qq')H-tagged","L");
  leg.AddEntry(scan_monoj_obs,"ggH-tagged","L");
  leg.Draw("same");

  TLegend leg2 (0.48,0.17,0.65,0.26);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(scan_combined_obs,"Observed","L");
  leg2.AddEntry(scan_combined_exp,"Expected","L");
  leg2.Draw("same");


  CMS_lumi(canvas,"35.9");
  canvas->SaveAs((outputPlots+"/scan_profile_likelihood.png").c_str(),"png");
  canvas->SaveAs((outputPlots+"/scan_profile_likelihood.pdf").c_str(),"pdf");  
  canvas->SaveAs((outputPlots+"/scan_profile_likelihood.C").c_str(),"C");  
}
