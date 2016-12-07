#include <cmath>
#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

using namespace std;

float musf = 0.02;
float elsf = 0.02;
float phsf = 0.02;
float mutrack = 0.01;
float eltrack = 0.01;
float mettrig = 0.01;
float eltrig = 0.02;
float phtrig = 0.02;
float lepveto = 0.03;

/////////////////////////////
void rzmm(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("czmm", "czmm", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    vector<double> bins = selectBinning(observable,category);
    TH1F*  hist = (TH1F*)file->FindObjectAny(("zmmcorhist_"+observable).c_str());
    if(isEWK)
      hist = (TH1F*)file->FindObjectAny(("zewkmmcorhist_"+observable).c_str());

    TH1F* ehist = (TH1F*)hist->Clone("ehist_");
    TH1* frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 18., "");
    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 18., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 2.0, bins.back(), 25., "");
	frame->GetXaxis()->SetTitle("Recoil [GeV]");
      }
      else if(TString(observable).Contains("mjj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 6.0, bins.back(), 18., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 30., "");

	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 6.0, bins.back(), 18., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 30., "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 6.0, bins.back(), 18., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 30., "");

	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    
    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{Z(#mu#mu)}");
    else
      frame->GetYaxis()->SetTitle("R_{Z(#mu#mu)-EWK}");

    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);


    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);
    
    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError(i)*hist->GetBinError(i);
        err += pow(hist->GetBinContent(i)*musf*2, 2);
        err += pow(hist->GetBinContent(i)*mutrack*2, 2);	
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK)
      leg->AddEntry(hist  , "R(Z(#mu#mu)) Stat. Unc.","pl");
    else
      leg->AddEntry(hist  , "R(Z(#mu#mu)-EWK) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rzmm.pdf");
      canvas->SaveAs("rzmm.png");
    }
    else{
      canvas->SaveAs("rzewkmm.pdf");
      canvas->SaveAs("rzewkmm.png");
    }

}


/////////////////////////////
void rzee(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("czee", "czee", 600, 600);
    canvas->SetRightMargin(0.075);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    vector<double> bins = selectBinning(observable,category);

    TH1F*  hist = (TH1F*)file->FindObjectAny(("zeecorhist_"+observable).c_str());
    if(isEWK)
      hist = (TH1F*)file->FindObjectAny(("zewkeecorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    TH1* frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 17., "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 22., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 35., "");

	frame->GetXaxis()->SetTitle("Recoil [GeV]");
      }
      else if(TString(observable).Contains("mjj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 22., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 35., "");
	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 4.0, bins.back(), 22., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 35., "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 6.0, bins.back(), 22., "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0.0, bins.back(), 35., "");
	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");
    
    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{Z(ee)}");
    else
      frame->GetYaxis()->SetTitle("R_{Z(ee)-EWK}");
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);


    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError(i)*hist->GetBinError(i);
        err += pow(hist->GetBinContent(i)*elsf*2, 2);
        err += pow(hist->GetBinContent(i)*eltrack*2, 2);
        err += pow(hist->GetBinContent(i)*mettrig, 2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK)
      leg->AddEntry(hist  , "R(Z(ee)) Stat. Unc.","pl");
    else
      leg->AddEntry(hist  , "R(Z(ee)-EWK) Stat. Unc.","pl");

    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rzee.pdf");
      canvas->SaveAs("rzee.png");
    }
    else{
      canvas->SaveAs("rzewkee.pdf");
      canvas->SaveAs("rzewkee.png");
    }
}

/////////////////////////////
void rwmn(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("cwmn", "cwmn", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());
    TH1F*  hist = (TH1F*)file->FindObjectAny(("wmncorhist_"+observable).c_str());
    if(isEWK)
      hist = (TH1F*)file->FindObjectAny(("wewkmncorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.5, "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.5, "");
	frame->GetXaxis()->SetTitle("Recoil [GeV]");
      }
      else if(TString(observable).Contains("mjj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 1.5, "");
	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 1.5, "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 1.5, "");
	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{W(#mu#nu)}");
    else
      frame->GetYaxis()->SetTitle("R_{W(#mu#nu)-EWK}");
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError(i)*hist->GetBinError(i);
        err += pow(hist->GetBinContent(i)*musf, 2);
        err += pow(hist->GetBinContent(i)*mutrack, 2);
        err += pow(hist->GetBinContent(i)*lepveto, 2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK)
      leg->AddEntry(hist  , "R(W(#mu#nu)) Stat. Unc.","pl");
    else
      leg->AddEntry(hist  , "R(W(#mu#nu)-EWK) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");
    
    if(not isEWK){
      canvas->SaveAs("rwmn.pdf");
      canvas->SaveAs("rwmn.png");
    }
    else{
      canvas->SaveAs("rwewkmn.pdf");
      canvas->SaveAs("rwewkmn.png");
    }
}

/////////////////////////////
void rwen(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("cwen", "cwen", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->FindObjectAny(("wencorhist_"+observable).c_str());
    if(isEWK)
      hist = (TH1F*)file->FindObjectAny(("wewkencorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 2.5, "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 2.5, "");
	frame->GetXaxis()->SetTitle("Recoil [GeV]");
      }
      else if(TString(observable).Contains("mjj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.5, "");
	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.5, "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.5, "");
	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else 
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{W(e#nu)}");
    else
      frame->GetYaxis()->SetTitle("R_{W(e#nu)}-EWK");

    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError(i)*hist->GetBinError(i);
        err += pow(hist->GetBinContent(i)*elsf, 2);
        err += pow(hist->GetBinContent(i)*eltrack, 2);
        err += pow(hist->GetBinContent(i)*eltrig, 2);
        err += pow(hist->GetBinContent(i)*mettrig, 2);
        err += pow(hist->GetBinContent(i)*lepveto, 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK)
      leg->AddEntry(hist  , "R(W(e#nu)) Stat. Unc.","pl");
    else
      leg->AddEntry(hist  , "R(W(e#nu)-EWK) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rwen.pdf");
      canvas->SaveAs("rwen.png");
    }
    else{
      canvas->SaveAs("rwewken.pdf");
      canvas->SaveAs("rwewken.png");
    }
}

/////////////////////////////
void rgam(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("cgam", "cgam", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist    = (TH1F*)file->FindObjectAny(("gamcorewkhist_"+observable).c_str());
    if(isEWK)
      hist    = (TH1F*)file->FindObjectAny(("gamewkcorewkhist_"+observable).c_str());    
    TH1F*  ewkhist = (TH1F*)file->FindObjectAny(("ZG_EWK_"+observable).c_str());
    TH1F*  re1hist = (TH1F*)file->FindObjectAny(("ZG_RenScale1_"+observable).c_str());
    TH1F*  re2hist = (TH1F*)file->FindObjectAny(("ZG_RenScale2_"+observable).c_str());
    TH1F*  fa1hist = (TH1F*)file->FindObjectAny(("ZG_FactScale1_"+observable).c_str());
    TH1F*  fa2hist = (TH1F*)file->FindObjectAny(("ZG_FactScale2_"+observable).c_str());
    TH1F*  pdfhist = (TH1F*)file->FindObjectAny(("ZG_PDF_"+observable).c_str());
    TH1F*  fophist = (TH1F*)file->FindObjectAny(("ZG_Footprint_"+observable).c_str());

    TH1F* ehistEWK    = (TH1F*)hist->Clone("ehistEWK");
    TH1F* ehist       = (TH1F*)hist->Clone("ehist");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0.2, bins.back(), 1.0, "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi"))
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
      else if(TString(observable).Contains("mjj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{#gamma}");
    else
      frame->GetYaxis()->SetTitle("R_{#gamma}-EWK");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    if(not isEWK){
      ehistEWK->SetFillColor(kOrange+1);
      ehistEWK->SetLineColor(kBlack);
      ehist->SetFillColor(kGreen+1);
      ehist->SetLineColor(kBlack);
    }
    else{
      ehist->SetFillColor(kOrange+1);
      ehist->SetLineColor(kBlack);
    }
    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
	if(not isEWK){
	  err += pow(ewkhist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  ehistEWK->SetBinError(i, sqrt(err));
	  err += pow(re1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(re2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fa1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fa2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(pdfhist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fophist->GetBinContent(i)*hist->GetBinContent(i), 2);       
	}
        err += pow(hist->GetBinContent(i)*phtrig, 2);
        err += pow(hist->GetBinContent(i)*phsf, 2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    if(not isEWK)
      ehistEWK->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK){
      leg->AddEntry(hist  , "R_{#gamma} Stat. Unc.","pl");
      leg->AddEntry(ehistEWK, "Stat. + Syst. (EWK)", "F");
      leg->AddEntry(ehist, "Stat. + Syst. (QCD+EWK+PDF)", "F");
    }
    else{
      leg->AddEntry(hist  , "R_{#gamma}-EWK Stat. Unc.","pl");
      leg->AddEntry(ehist, "Stat. + Syst.", "F");
    }

    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rgam.pdf");
      canvas->SaveAs("rgam.png");
    }
    else{
      canvas->SaveAs("rgamewk.pdf");
      canvas->SaveAs("rgamewk.png");
    }
}

/////////////////////////////
void rzwj(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("czwj", "czwj", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->FindObjectAny(("zwjcorewkhist_"+observable).c_str());
    if(isEWK)
      hist = (TH1F*)file->FindObjectAny(("zwjewkcorhist_"+observable).c_str());

    TH1F*  ewkhist = (TH1F*)file->FindObjectAny(("ZW_EWK_"+observable).c_str());
    TH1F*  re1hist = (TH1F*)file->FindObjectAny(("ZW_RenScale1_"+observable).c_str());
    TH1F*  re2hist = (TH1F*)file->FindObjectAny(("ZW_RenScale2_"+observable).c_str());
    TH1F*  fa1hist = (TH1F*)file->FindObjectAny(("ZW_FactScale1_"+observable).c_str());
    TH1F*  fa2hist = (TH1F*)file->FindObjectAny(("ZW_FactScale2_"+observable).c_str());
    TH1F*  pdfhist = (TH1F*)file->FindObjectAny(("ZW_PDF_"+observable).c_str());

    TH1F* ehist    = (TH1F*)hist->Clone("ehist");
    TH1F* ehistEWK = (TH1F*)hist->Clone("ehistEWK");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 12.0, "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 3.0, "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 4.0, "");

	frame->GetXaxis()->SetTitle("Recoil [GeV]");
      }
      else if(TString(observable).Contains("mjj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 3.0, "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 4.0, "");

	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 3.0, "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 4.0, "");

	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
	if(not isEWK)
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 3.0, "");
	else
	  frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 4.0, "");

	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{Z/W}");
    else
      frame->GetYaxis()->SetTitle("R_{Z/W}-EWK");

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);
    ehist->SetFillColor(kGreen+1);
    ehist->SetLineColor(kBlack);
    ehistEWK->SetFillColor(kOrange+1);
    ehistEWK->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
	err += pow(ewkhist->GetBinContent(i)*hist->GetBinContent(i), 2);
	ehistEWK->SetBinError(i, sqrt(err));
	err += pow(re1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	err += pow(re2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	err += pow(fa1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	err += pow(fa2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	err += pow(pdfhist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(hist->GetBinContent(i)*lepveto, 2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    ehistEWK->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK){
      leg->AddEntry(hist  , "R(Z/W) Stat. Unc.","pl");
      leg->AddEntry(ehistEWK, "Stat. + Syst. (EWK)", "F");
      leg->AddEntry(ehist, "Stat. + Syst. (QCD+EWK+PDF)", "F");
    }
    else{
      leg->AddEntry(hist  , "R(Z/W)-EWK Stat. Unc.","pl");
      leg->AddEntry(ehistEWK, "Stat. + Syst. (QCD+EWK+PDF)", "F");
      leg->AddEntry(ehist, "Stat. + Syst.", "F");
    }
    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rzwj.pdf");
      canvas->SaveAs("rzwj.png");
    }
    else{
      canvas->SaveAs("rzwjewk.pdf");
      canvas->SaveAs("rzwjewk.png");
    }
}


/////////////////////////////
void rwgam(string fileName, Category category, string observable, bool isEWK) {

    TCanvas* canvas = new TCanvas("cwgam", "cwgam", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist    = (TH1F*)file->FindObjectAny(("wgamcorewkhist_"+observable).c_str());
    if(isEWK)
      hist    = (TH1F*)file->FindObjectAny(("wgamewkcorhist_"+observable).c_str());

    TH1F*  ewkhist = (TH1F*)file->FindObjectAny(("WG_EWK_"+observable).c_str());
    TH1F*  re1hist = (TH1F*)file->FindObjectAny(("WG_RenScale1_"+observable).c_str());
    TH1F*  re2hist = (TH1F*)file->FindObjectAny(("WG_RenScale2_"+observable).c_str());
    TH1F*  fa1hist = (TH1F*)file->FindObjectAny(("WG_FactScale1_"+observable).c_str());
    TH1F*  fa2hist = (TH1F*)file->FindObjectAny(("WG_FactScale2_"+observable).c_str());
    TH1F*  pdfhist = (TH1F*)file->FindObjectAny(("WG_PDF_"+observable).c_str());
    TH1F*  fophist = (TH1F*)file->FindObjectAny(("WG_Footprint_"+observable).c_str());

    TH1F* ehistEWK    = (TH1F*)hist->Clone("ehistEWK");
    TH1F* ehist       = (TH1F*)hist->Clone("ehist");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.0, "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(category == Category::VBF){
      if(TString(observable).Contains("met") and TString(observable).Contains("jetmetdphi"))
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
      else if(TString(observable).Contains("mjj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("M_{jj} [GeV]");
      }
      else if(TString(observable).Contains("jetmetdphi")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("#Delta#phi(jet,met)");
      }
      else if(TString(observable).Contains("detajj")){
        frame = canvas->DrawFrame(bins.front(), 0, bins.back(), 2.0, "");
	frame->GetXaxis()->SetTitle("#Delta#eta_{jj}");
      }
    }
    else
      frame->GetXaxis()->SetTitle("Recoil [GeV]");

    if(not isEWK)
      frame->GetYaxis()->SetTitle("R_{W#gamma}");
    else
      frame->GetYaxis()->SetTitle("R_{W#gamma}-EWK");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);
    
    if(not isEWK){
      ehistEWK->SetFillColor(kOrange+1);
      ehistEWK->SetLineColor(kBlack);
      ehist->SetFillColor(kGreen+1);
      ehist->SetLineColor(kBlack);
    }
    else{
      ehist->SetFillColor(kOrange+1);
      ehist->SetLineColor(kBlack);
    }

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
	if(not isEWK){
	  err += pow(ewkhist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  ehistEWK->SetBinError(i, sqrt(err));
	  err += pow(re1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(re2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fa1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fa2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(pdfhist->GetBinContent(i)*hist->GetBinContent(i), 2);
	  err += pow(fophist->GetBinContent(i)*hist->GetBinContent(i), 2);
	}
	err += pow(hist->GetBinContent(i)*phsf, 2);
	err += pow(hist->GetBinContent(i)*phtrig, 2);
	err += pow(hist->GetBinContent(i)*mettrig, 2);
        err += pow(hist->GetBinContent(i)*lepveto, 2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    if(not isEWK)
      ehistEWK->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(not isEWK){
      leg->AddEntry(hist  , "R_{W#gamma} Stat. Unc.","pl");
      leg->AddEntry(ehistEWK, "Stat. + Syst. (EWK)", "F");
      leg->AddEntry(ehist, "Stat. + Syst. (QCD+EWK+PDF)", "F");
    }
    else{
      leg->AddEntry(hist  , "R_{W#gamma}-EWK Stat. Unc.","pl");
      leg->AddEntry(ehistEWK, "Stat. + Syst.", "F");
    }
    leg->Draw("SAME");

    if(not isEWK){
      canvas->SaveAs("rwgam.pdf");
      canvas->SaveAs("rwgam.png");
    }
    else{
      canvas->SaveAs("rwgamewk.pdf");
      canvas->SaveAs("rwgamewk.png");
    }
}


/////////////////////////////
void rtopmu(string fileName, Category category, string observable) {

    TCanvas* canvas = new TCanvas("ctopmu", "ctopmu", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist  = (TH1F*)file->FindObjectAny(("topmucorhist_"+observable).c_str());
    TH1F*  histb = (TH1F*)file->FindObjectAny(("TOP_MU_B_"+observable).c_str());

    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1., "");
    frame->GetYaxis()->SetTitle("R_{top,#mu}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.02,2);
        err += pow(fabs(histb->GetBinContent(i))*hist->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(top,#mu) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    ;
    leg->Draw("SAME");

    canvas->SaveAs("rtopmu.pdf");
    canvas->SaveAs("rtopmu.png");
}


/////////////////////////////
void rtopel(string fileName, Category category, string observable) {

    TCanvas* canvas = new TCanvas("ctopel", "ctopel", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist  = (TH1F*)file->FindObjectAny(("topelcorhist_"+observable).c_str());
    TH1F*  histb = (TH1F*)file->FindObjectAny(("TOP_EL_B_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<double> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1., "");
    frame->GetYaxis()->SetTitle("R_{top,el}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->CenterTitle();

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(1.2);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);

    ehist->SetFillColor(kOrange+1);
    ehist->SetLineColor(kBlack);

    for (int i = 1; i <= ehist->GetNbinsX(); i++) {
        double err = 0.0;
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.02,2);
        err += pow(histb->GetBinContent(i)*hist->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas,"35.9",true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(top,el) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rtopel.pdf");
    canvas->SaveAs("rtopel.png");
}

/////////////////////////////
void makeTransferFactorPlot(string fileName, Category category, string observable, bool isEWK = false, bool addwgamma = false, bool addtop = false) {

  gROOT->SetBatch(kTRUE);

  initializeBinning();

  rzmm(fileName,category,observable,isEWK);
  rzee(fileName,category,observable,isEWK);
  rwmn(fileName,category,observable,isEWK);
  rwen(fileName,category,observable,isEWK);
  rzwj(fileName,category,observable,isEWK);

  if(category != Category::VBF)
    rgam(fileName,category,observable,false);
  
  if(addtop){
    rtopmu(fileName,category,observable);
    rtopel(fileName,category,observable);
  }
  
  if(addwgamma and category != Category::VBF)
    rwgam(fileName,category,observable,false);
}


