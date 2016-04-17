#include <cmath>
#include "CMS_lumi.h"
#include "histoUtils.h"

using namespace std;

void rzmm(string fileName, int category, string observable, bool alongX = false) {

  TFile* file = new TFile(fileName.c_str());  
  TH1F*  hist = (TH1F*)file->FindObjectAny(("zmmcorhist_"+observable).c_str());
  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);

  TCanvas* canvas = new TCanvas("czmm", "czmm", 600, 600);
  
  for(vector<TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);
  
    TH1F* ehist = (TH1F*) (*hist)->Clone("ehist");

    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitle("R_{Z(#mu#mu)}");
    (*hist)->GetYaxis()->CenterTitle();
    (*hist)->GetYaxis()->SetTitleOffset(1.15);
    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);
    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);
    (*hist)->GetYaxis()->SetRangeUser(5.,15.);
    
    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= (*hist)->GetNbinsX(); i++) {
      double err = 0.0;
      err += (*hist)->GetBinError(i)*(*hist)->GetBinError(i);
      err += pow((*hist)->GetBinContent(i)*0.02, 2);      
      ehist->SetBinError(i, sqrt(err));
    }

    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist)->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist)->Draw("PE SAME");
    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(*hist,"R(Z(#mu#mu)) Stat. Unc.","pl");
    leg->AddEntry(ehist,"Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");
    
    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%2.f <= %s < %.2f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.2f ",text.second.c_str(),bin.at(ihisto)));
   
    canvas->SaveAs(Form("rzmm_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rzmm_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;
  }
}


void rzee(string fileName, int category, string observable, bool alongX = false) {

  TFile* file = new TFile(fileName.c_str());  
  TH1F*  hist = (TH1F*)file->FindObjectAny(("zeecorhist_"+observable).c_str());
  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);

  TCanvas* canvas = new TCanvas("czee", "czee", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetRightMargin(0.075);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);  
    
    TH1F* ehist = (TH1F*)(*hist)->Clone("ehist");
    (*hist)->GetYaxis()->SetTitle("R_{Z(ee)}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->CenterTitle();
    (*hist)->GetYaxis()->SetTitleOffset(1.15);
    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);
    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);
    (*hist)->GetYaxis()->SetRangeUser(5.,20.);
    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
      double err = 0.0;
      err +=    (*hist)->GetBinError(i)*(*hist)->GetBinError(i);
      err += pow((*hist)->GetBinContent(i)*0.04, 2);
      err += pow((*hist)->GetBinContent(i)*0.01, 2);      
      ehist->SetBinError(i, sqrt(err));
    }
    
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");

    canvas->RedrawAxis();
    
    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist), "R(Z(ee)) Stat. Unc.","pl");
    leg->AddEntry(ehist, "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rzee_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rzee_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;
  }
}


void rwmn(string fileName, int category, string observable, bool alongX = false) {

  TFile* file = new TFile(fileName.c_str());  
  TH1F*  hist = (TH1F*)file->FindObjectAny(("wmncorhist_"+observable).c_str());
  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);
  TCanvas* canvas = new TCanvas("cwmn", "cwmn", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehist = (TH1F*)(*hist)->Clone("ehist");

    (*hist)->GetYaxis()->SetTitle("R_{W(#mu#nu)}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);
    (*hist)->GetXaxis()->SetLabelSize(0.040);

    (*hist)->GetYaxis()->CenterTitle();

    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);

    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);
    (*hist)->GetYaxis()->SetRangeUser(0.,1.);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
      double err = 0.0;
      err +=    (*hist)->GetBinError(i)*(*hist)->GetBinError(i);
      err += pow((*hist)->GetBinContent(i)*0.01, 2);
      err += pow((*hist)->GetBinContent(i)*0.03, 2);
      
      ehist->SetBinError(i, sqrt(err));
    }
    
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R(W(#mu#nu)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rwmn_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rwmn_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;
  }
}

void rwen(string fileName, int category, string observable, bool alongX = false) {

  TFile* file = new TFile(fileName.c_str());  
  TH1F*  hist = (TH1F*)file->FindObjectAny(("wencorhist_"+observable).c_str());
  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);  
  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);

  TCanvas* canvas = new TCanvas("cwen", "cwen", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehist = (TH1F*)(*hist)->Clone("ehist");

    (*hist)->GetYaxis()->SetTitle("R_{W(e#nu)}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);
    (*hist)->GetYaxis()->CenterTitle();
    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);
    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);
    
    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
      double err = 0.0;
      err +=    (*hist)->GetBinError(i)*(*hist)->GetBinError(i);
      err += pow((*hist)->GetBinContent(i)*0.05, 2);
      err += pow((*hist)->GetBinContent(i)*0.01, 2);
      err += pow((*hist)->GetBinContent(i)*0.03, 2);      
      ehist->SetBinError(i, sqrt(err));
    }
    
    (*hist)->GetYaxis()->SetRangeUser(0.,2.0);
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");
    
    canvas->RedrawAxis();
    
    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R(W(e#nu)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);
 
    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rwen_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rwen_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;

  }
}

void rgam(string fileName, int category, string observable, bool alongX = false) {

  TFile* file    = new TFile(fileName.c_str());  
  TH1F*  hist    = (TH1F*)file->FindObjectAny(("gamcorewkhist_"+observable).c_str());
  TH1F*  ewkhist = (TH1F*)file->FindObjectAny(("ZG_EWK_"+observable).c_str());
  TH1F*  re1hist = (TH1F*)file->FindObjectAny(("ZG_RenScale1_"+observable).c_str());
  TH1F*  re2hist = (TH1F*)file->FindObjectAny(("ZG_RenScale2_"+observable).c_str());
  TH1F*  fa1hist = (TH1F*)file->FindObjectAny(("ZG_FactScale1_"+observable).c_str());
  TH1F*  fa2hist = (TH1F*)file->FindObjectAny(("ZG_FactScale2_"+observable).c_str());
  TH1F*  pdfhist = (TH1F*)file->FindObjectAny(("ZG_PDF_"+observable).c_str());
  TH1F*  fophist = (TH1F*)file->FindObjectAny(("ZG_Footprint_"+observable).c_str());
  
  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  vector<TH1F*> histograms_ewk = transformUnrolledHistogram(ewkhist,observable,category,alongX);
  vector<TH1F*> histograms_re1 = transformUnrolledHistogram(re1hist,observable,category,alongX);
  vector<TH1F*> histograms_re2 = transformUnrolledHistogram(re2hist,observable,category,alongX);
  vector<TH1F*> histograms_fa1 = transformUnrolledHistogram(fa1hist,observable,category,alongX);
  vector<TH1F*> histograms_fa2 = transformUnrolledHistogram(fa2hist,observable,category,alongX);
  vector<TH1F*> histograms_pdf = transformUnrolledHistogram(pdfhist,observable,category,alongX);
  vector<TH1F*> histograms_fop = transformUnrolledHistogram(fophist,observable,category,alongX);

  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);
  TCanvas* canvas = new TCanvas("cgam", "cgam", 600, 600);
    
  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehistEWK    = (TH1F*)(*hist)->Clone("ehistEWK");
    TH1F* ehist       = (TH1F*)(*hist)->Clone("ehist");

    (*hist)->GetYaxis()->SetTitle("R_{#gamma}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->CenterTitle();
    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);
    (*hist)->SetLineColor(kBlack);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);

    ehistEWK->SetFillColor(kOrange+1);
    ehistEWK->SetLineColor(kBlack);
    ehist->SetFillColor(kGreen+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
      double err = 0.0;
      err +=    (*hist)->GetBinError(i)*(*hist)->GetBinError  (i);
      err += pow(histograms_ewk.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      ehistEWK->SetBinError(i, sqrt(err));
      err += pow(histograms_re1.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_re2.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_fa1.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_fa2.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_pdf.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_fop.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      ehist->SetBinError(i, sqrt(err));
    }

    (*hist)->GetYaxis()->SetRangeUser(0.,1.);
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    ehistEWK->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");
    
    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R_{#gamma} Stat. Unc.","pl");
    leg->AddEntry(ehistEWK, "Stat. + Syst. (EWK)", "F");
    leg->AddEntry(ehist, "Stat. + Syst. (QCD+EWK+PDF)", "F");
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);
 
    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rgam_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rgam_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;

  }    
}

void rzwj(string fileName, int category, string observable, bool alongX = false) {

  TFile* file = new TFile(fileName.c_str());
  TH1F*  hist = (TH1F*)file->FindObjectAny(("zwjcorewkhist_"+observable).c_str());
  TH1F*  ewkhist = (TH1F*)file->FindObjectAny(("ZW_EWK_"+observable).c_str());
  TH1F*  re1hist = (TH1F*)file->FindObjectAny(("ZW_RenScale1_"+observable).c_str());
  TH1F*  re2hist = (TH1F*)file->FindObjectAny(("ZW_RenScale2_"+observable).c_str());
  TH1F*  fa1hist = (TH1F*)file->FindObjectAny(("ZW_FactScale1_"+observable).c_str());
  TH1F*  fa2hist = (TH1F*)file->FindObjectAny(("ZW_FactScale2_"+observable).c_str());
  TH1F*  pdfhist = (TH1F*)file->FindObjectAny(("ZW_PDF_"+observable).c_str());

  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  vector<TH1F*> histograms_ewk = transformUnrolledHistogram(ewkhist,observable,category,alongX);
  vector<TH1F*> histograms_re1 = transformUnrolledHistogram(re1hist,observable,category,alongX);
  vector<TH1F*> histograms_re2 = transformUnrolledHistogram(re2hist,observable,category,alongX);
  vector<TH1F*> histograms_fa1 = transformUnrolledHistogram(fa1hist,observable,category,alongX);
  vector<TH1F*> histograms_fa2 = transformUnrolledHistogram(fa2hist,observable,category,alongX);
  vector<TH1F*> histograms_pdf = transformUnrolledHistogram(pdfhist,observable,category,alongX);

  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);
  TCanvas* canvas = new TCanvas("czwj", "czwj", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehist    = (TH1F*)(*hist)->Clone("ehist");
    TH1F* ehistEWK = (TH1F*)(*hist)->Clone("ehistEWK");

    (*hist)->GetYaxis()->SetTitle("R_{Z/W}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);

    (*hist)->GetYaxis()->CenterTitle();

    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);

    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);

    ehist->SetFillColor(kGreen+1);
    ehist->SetLineColor(kBlack);
    ehistEWK->SetFillColor(kOrange+1);
    ehistEWK->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
      double err = 0.0;
      err +=    (*hist)->GetBinError  (i)*   (*hist)->GetBinError  (i);
      err += pow(histograms_ewk.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      ehistEWK->SetBinError(i, sqrt(err));
      err += pow(histograms_re1.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_re2.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_fa1.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_fa2.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      err += pow(histograms_pdf.at(ihisto)->GetBinContent(i)*(*hist)->GetBinContent(i), 2);
      ehist->SetBinError(i, sqrt(err));
    }

    (*hist)->GetYaxis()->SetRangeUser(0.,15.);
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    ehistEWK->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");
    
    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R(Z/W) Stat. Unc.","pl");
    leg->AddEntry(ehistEWK, "Stat. + Syst. (EWK)", "F");
    leg->AddEntry(ehist, "Stat. + Syst. (QCD+EWK+PDF)", "F");
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);
 
    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rwzj_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rwzj_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;

  }
}

void rtopmu(string fileName, int category, string observable, bool alongX) {

  TFile* file = new TFile(fileName.c_str());
  TH1F*  hist  = (TH1F*)file->FindObjectAny(("topmucorhist_"+observable).c_str());
  TH1F*  histb = (TH1F*)file->FindObjectAny(("TOP_MU_B_"+observable).c_str());

  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  vector<TH1F*> histograms_b = transformUnrolledHistogram(histb,observable,category,alongX);

  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);

  TCanvas* canvas = new TCanvas("ctopmu", "ctopmu", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehist = (TH1F*)(*hist)->Clone("ehist");

    (*hist)->GetYaxis()->SetTitle("R_{top,#mu}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);

    (*hist)->GetYaxis()->CenterTitle();

    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);

    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    (*hist)->GetBinError  (i)*   (*hist)->GetBinError  (i);
        err += pow((*hist)->GetBinContent(i)*0.02,2);
        err += pow(fabs(histograms_b.at(ihisto)->GetBinContent(i))*(*hist)->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    (*hist)->GetYaxis()->SetRangeUser(0.,0.4);
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R(top,#mu) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);
 
    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rtopmu_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rtopmu_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;
  }
}



void rtopel(string fileName, int category, string observable, bool alongX) {

  TFile* file = new TFile(fileName.c_str());
  TH1F*  hist  = (TH1F*)file->FindObjectAny(("topelcorhist_"+observable).c_str());
  TH1F*  histb = (TH1F*)file->FindObjectAny(("TOP_EL_B_"+observable).c_str());

  vector<TH1F*> histograms = transformUnrolledHistogram(hist,observable,category,alongX);
  vector<TH1F*> histograms_b = transformUnrolledHistogram(histb,observable,category,alongX);

  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != histograms.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  int ihisto = 0;
  pair<string,string> text = observableName(observable,alongX);
  TCanvas* canvas = new TCanvas("ctopel", "ctopel", 600, 600);

  for(vector<TH1F*>::iterator hist = histograms.begin(); hist != histograms.end(); hist++){

    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TH1F* ehist = (TH1F*)(*hist)->Clone("ehist");

    (*hist)->GetYaxis()->SetTitle("R_{top,el}");
    (*hist)->GetXaxis()->SetTitle(text.first.c_str());
    (*hist)->GetXaxis()->SetTitleSize(0.045);
    (*hist)->GetXaxis()->SetLabelSize(0.040);

    (*hist)->GetYaxis()->CenterTitle();

    (*hist)->GetYaxis()->SetLabelSize(0.040);
    (*hist)->GetYaxis()->SetTitleSize(0.045);
    (*hist)->GetYaxis()->SetTitleOffset(1.15);

    (*hist)->SetLineColor(kBlack);
    (*hist)->SetLineWidth(1);
    (*hist)->SetMarkerSize(1.2);
    (*hist)->SetMarkerStyle(20);
    (*hist)->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    (*hist)->GetBinError  (i)*   (*hist)->GetBinError  (i);
        err += pow((*hist)->GetBinContent(i)*0.02,2);
        err += pow(fabs(histograms_b.at(ihisto)->GetBinContent(i))*(*hist)->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    (*hist)->GetYaxis()->SetRangeUser(0.,0.4);
    (*hist)->Draw();
    CMS_lumi(canvas,"2.30",true);
    (*hist) ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    (*hist) ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.5, 0.7, 0.85, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry((*hist)  , "R(top,#mu) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    
    leg->Draw("SAME");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%1.f <= %s < %.1f ",bin.at(ihisto),text.second.c_str(),bin.at(ihisto+1)));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihisto)));

    canvas->SaveAs(Form("rtopel_%s_bin_%d.pdf",observable.c_str(),ihisto));
    canvas->SaveAs(Form("rtopel_%s_bin_%d.png",observable.c_str(),ihisto));
    ihisto++;


  }
}



void makeTransferFactorPlot(string fileName, int category, string observable, bool addtop = false, bool alongX = false) {

  gROOT->SetBatch(kTRUE);

  setTDRStyle();

  rzmm(fileName,category,observable,alongX);
  rzee(fileName,category,observable,alongX);
  rwmn(fileName,category,observable,alongX);
  rwen(fileName,category,observable,alongX);
  rgam(fileName,category,observable,alongX);
  rzwj(fileName,category,observable,alongX);
  
  if(addtop){
    rtopmu(fileName,category,observable,alongX);
    rtopel(fileName,category,observable,alongX);
  }  
}


