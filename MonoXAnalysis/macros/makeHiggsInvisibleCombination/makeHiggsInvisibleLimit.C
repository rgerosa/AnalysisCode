#include "../CMS_lumi.h"

void  fillPoint(TTree* tree, TGraphAsymmErrors* limit_obs, TGraphAsymmErrors* limit_exp, TGraphAsymmErrors* limit_1s, TGraphAsymmErrors* limit_2s, const int & ipoint){

  TTreeReader reader(tree);
  TTreeReaderValue<double> limit (reader,"limit");
  TTreeReaderValue<float> quantile (reader,"quantileExpected");

  float err_dw_2s = 0;
  float err_dw_1s = 0;
  float err_up_2s = 0;
  float err_up_1s = 0;
  float val_expected = 0;

  while(reader.Next()){
    if(*quantile == -1){
      limit_obs->SetPoint(ipoint,ipoint+0.5,*limit);
      limit_obs->SetPointError(ipoint,0.5,0.5,0,0);
    }      
    else if(*quantile == 0.5){
      limit_exp->SetPoint(ipoint,ipoint+0.5,*limit);
      limit_1s->SetPoint(ipoint,ipoint+0.5,*limit);
      limit_2s->SetPoint(ipoint,ipoint+0.5,*limit);
      limit_exp->SetPointError(ipoint,0.5,0.5,0,0);
      val_expected = *limit;
    }
    else if(*quantile < 0.03) err_dw_2s = *limit;
    else if(*quantile > 0.14 and *quantile < 0.18) err_dw_1s= *limit;
    else if(*quantile > 0.96) err_up_2s = *limit;
    else if(*quantile > 0.83 and *quantile < 0.86) err_up_1s= *limit;
  }
  limit_1s->SetPointError(ipoint,0.5,0.5,fabs(err_dw_1s-val_expected),fabs(err_up_1s-val_expected));
  limit_2s->SetPointError(ipoint,0.5,0.5,fabs(err_dw_2s-val_expected),fabs(err_up_2s-val_expected));
  

}  

void makeHiggsInvisibleLimit(string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);


  TFile* file_monoj = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoJ/higgsCombine_limit_unblind.AsymptoticLimits.mH125.root","READ");
  TFile* file_monov = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-048/MonoV/higgsCombine_limit_unblind.AsymptoticLimits.mH125.root","READ");
  TFile* file_monoz = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/EXO-16-052/higgsCombine_limit_unblind.AsymptoticLimits.mH125.root","READ");
  TFile* file_vbf   = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/HIG-17-023/higgsCombine_limit_unblind.AsymptoticLimits.mH125.root","READ");
  TFile* file_combined = new TFile("../makeWorkspace/HiggsInvisible/HiggsInvisibleCombination/Combination/higgsCombine_limit_unblind.AsymptoticLimits.mH125.root","READ");

  /////                                                                                                                                                                                                
  TTree* limit_monoj = (TTree*) file_monoj->Get("limit");
  TTree* limit_monoz = (TTree*) file_monoz->Get("limit");
  TTree* limit_monov = (TTree*) file_monov->Get("limit");
  TTree* limit_vbf = (TTree*) file_vbf->Get("limit");
  TTree* limit_combined = (TTree*) file_combined->Get("limit");

  ///
  TGraphAsymmErrors* limit_obs = new TGraphAsymmErrors();
  TGraphAsymmErrors* limit_exp = new TGraphAsymmErrors();
  TGraphAsymmErrors* limit_1s = new TGraphAsymmErrors();
  TGraphAsymmErrors* limit_2s = new TGraphAsymmErrors();

  
  int ipoint = 0;
  fillPoint(limit_combined,limit_obs,limit_exp,limit_1s,limit_2s,ipoint);
  ipoint++;

  fillPoint(limit_vbf,limit_obs,limit_exp,limit_1s,limit_2s,ipoint);
  ipoint++;

  fillPoint(limit_monoz,limit_obs,limit_exp,limit_1s,limit_2s,ipoint);
  ipoint++;

  fillPoint(limit_monov,limit_obs,limit_exp,limit_1s,limit_2s,ipoint);
  ipoint++;

  fillPoint(limit_monoj,limit_obs,limit_exp,limit_1s,limit_2s,ipoint);
  ipoint++;

  TH1* frame = (TH1*) canvas->DrawFrame(0.,0.,5,1.5);
  frame->SetBins(5,0,5);
  frame->GetXaxis()->SetBinLabel(1,"Combined");
  frame->GetXaxis()->SetBinLabel(2,"VBF-tag");
  frame->GetXaxis()->SetBinLabel(3,"Z(ll)H-tag");
  frame->GetXaxis()->SetBinLabel(4,"V(qq')H-tag");
  frame->GetXaxis()->SetBinLabel(5,"ggH-tag");
  frame->GetXaxis()->LabelsOption("h");
  frame->GetYaxis()->SetRangeUser(0.,1.5);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitle("95% CL upper limit on #sigma x B(H #rightarrow inv.)/#sigma_{SM}");
  frame->GetYaxis()->SetTitleSize(0.038);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->Draw();

  limit_exp->SetLineColor(kBlack);
  limit_exp->SetLineWidth(2);
  limit_exp->SetLineStyle(2);
  limit_exp->Draw("PE0same");
  limit_1s->SetFillColor(kGreen+1);
  limit_1s->SetLineColor(kGreen+1);
  limit_2s->SetFillColor(kOrange);
  limit_2s->SetLineColor(kOrange);
  limit_2s->Draw("2same");
  limit_1s->Draw("2same");
  limit_exp->SetMarkerSize(1.0);
  limit_exp->SetMarkerStyle(24);
  limit_exp->SetMarkerColor(kBlack);
  limit_exp->Draw("PE0 same");
  limit_obs->SetMarkerSize(1.0);
  limit_obs->SetMarkerStyle(20);
  limit_obs->SetMarkerColor(kBlack);
  limit_obs->Draw("PE0 same");
  
  TLegend* leg = new TLegend(0.16,0.57,0.45,0.84);
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(limit_obs,"Observed","PL");
  leg->AddEntry(limit_exp,"Median expected","PL");
  leg->AddEntry(limit_1s,"68% expected","F");
  leg->AddEntry(limit_2s,"95% expected","F");
  leg->Draw("same");
  CMS_lumi(canvas,"35.9",false);

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/higgsinvisble_summary.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/higgsinvisble_summary.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/higgsinvisble_summary.C").c_str(),"C");


}
