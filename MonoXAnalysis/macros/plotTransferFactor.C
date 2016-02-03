#include <cmath>
#include "CMS_lumi.h"
#include "makehist.h"

using namespace std;

void rzmm(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("czmm", "czmm", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    vector<float> bins = selectBinning(observable,category);

    TH1F*  hist = (TH1F*)file->Get(("zmmcorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    TH1* frame = canvas->DrawFrame(bins.front(), 6.0, bins.back(), 12., "");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);

    frame->GetYaxis()->SetTitle("R_{Z(#mu#mu)}");
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
        err += pow(hist->GetBinContent(i)*0.02, 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(Z(#mu#mu)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rzmm.pdf");
    canvas->SaveAs("rzmm.png");

    //file->Close();
}


void rzee(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("czee", "czee", 600, 600);
    canvas->SetRightMargin(0.075);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    vector<float> bins = selectBinning(observable,category);

    TH1F*  hist = (TH1F*)file->Get(("zeecorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");
    TH1* frame = canvas->DrawFrame(bins.front(), 5.0, bins.back(), 30., "");
    frame->GetYaxis()->SetTitle("R_{Z(ee)}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err += pow(hist->GetBinContent(i)*0.04, 2);
        err += pow(hist->GetBinContent(i)*0.01, 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(Z(ee)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rzee.pdf");
    canvas->SaveAs("rzee.png");

    //file->Close();
}

void rwmn(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("cwmn", "cwmn", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());
    TH1F*  hist = (TH1F*)file->Get(("wmncorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.0, "");
    frame->GetYaxis()->SetTitle("R_{W(#mu#nu)}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err += pow(hist->GetBinContent(i)*0.01, 2);
        err += pow(hist->GetBinContent(i)*0.03, 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(W(#mu#nu)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rwmn.pdf");
    canvas->SaveAs("rwmn.png");

    //file->Close();
}

void rwen(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("cwen", "cwen", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->Get(("wencorhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist_");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 2., "");
    frame->GetYaxis()->SetTitle("R_{W(e#nu)}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err += pow(hist->GetBinContent(i)*0.05, 2);
        err += pow(hist->GetBinContent(i)*0.01, 2);
        err += pow(hist->GetBinContent(i)*0.03, 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(W(e#nu)) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rwen.pdf");
    canvas->SaveAs("rwen.png");

    //file->Close();
}

void rgam(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("cgam", "cgam", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist    = (TH1F*)file->Get(("gamcorewkhist_"+observable).c_str());
    TH1F*  ewkhist = (TH1F*)file->Get(("ZG_EWK_"+observable).c_str());
    TH1F*  re1hist = (TH1F*)file->Get(("ZG_RenScale1_"+observable).c_str());
    TH1F*  re2hist = (TH1F*)file->Get(("ZG_RenScale2_"+observable).c_str());
    TH1F*  fa1hist = (TH1F*)file->Get(("ZG_FactScale1_"+observable).c_str());
    TH1F*  fa2hist = (TH1F*)file->Get(("ZG_FactScale2_"+observable).c_str());
    TH1F*  pdfhist = (TH1F*)file->Get(("ZG_PDF_"+observable).c_str());
    TH1F*  fophist = (TH1F*)file->Get(("ZG_Footprint_"+observable).c_str());

    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0.2, bins.back(), 1.0, "");
    frame->GetYaxis()->SetTitle("R_{#gamma}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(ewkhist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(re1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(re2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(fa1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(fa2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(pdfhist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(fophist->GetBinContent(i)*hist->GetBinContent(i), 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();

    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R_{#gamma} Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rgam.pdf");
    canvas->SaveAs("rgam.png");

    //file->Close();
}

void rzwj(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("czwj", "czwj", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->Get(("zwjcorewkhist_"+observable).c_str());
    TH1F*  ewkhist = (TH1F*)file->Get(("ZW_EWK_"+observable).c_str());
    TH1F*  re1hist = (TH1F*)file->Get(("ZW_RenScale1_"+observable).c_str());
    TH1F*  re2hist = (TH1F*)file->Get(("ZW_RenScale2_"+observable).c_str());
    TH1F*  fa1hist = (TH1F*)file->Get(("ZW_FactScale1_"+observable).c_str());
    TH1F*  fa2hist = (TH1F*)file->Get(("ZW_FactScale2_"+observable).c_str());
    TH1F*  pdfhist = (TH1F*)file->Get(("ZW_PDF_"+observable).c_str());

    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 12.0, "");
    frame->GetYaxis()->SetTitle("R_{Z/W}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(ewkhist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(re1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(re2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(fa1hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(fa2hist->GetBinContent(i)*hist->GetBinContent(i), 2);
        err += pow(pdfhist->GetBinContent(i)*hist->GetBinContent(i), 2);

        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(Z/W) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rzwj.pdf");
    canvas->SaveAs("rzwj.png");

    //file->Close();
}



void rtopmu(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("ctopmu", "ctopmu", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist  = (TH1F*)file->Get(("topmucorhist_"+observable).c_str());
    TH1F*  histb = (TH1F*)file->Get(("TOP_MU_B_"+observable).c_str());

    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1., "");
    frame->GetYaxis()->SetTitle("R_{top,#mu}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.02,2);
        err += pow(fabs(histb->GetBinContent(i))*hist->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
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

    //file->Close();
}


void rtopel(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("ctopel", "ctopel", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist  = (TH1F*)file->Get(("topelcorhist_"+observable).c_str());
    TH1F*  histb = (TH1F*)file->Get(("TOP_EL_B_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1., "");
    frame->GetYaxis()->SetTitle("R_{top,el}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.02,2);
        err += pow(histb->GetBinContent(i)*hist->GetBinContent(i),2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
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

    //file->Close();
}

void rsidebandZ(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("csidebandZ", "csidebandZ", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->Get(("sidebandcorZhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 4.0, "");
    frame->GetYaxis()->SetTitle("R_{sideband,Z}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.10,2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(sideband,Z) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rsidebandZ.pdf");
    canvas->SaveAs("rsidebandZ.png");

    //file->Close();
}


void rsidebandW(string fileName, int category, string observable) {

    TCanvas* canvas = new TCanvas("csidebandW", "csidebandW", 600, 600);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.075);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.06);

    TFile* file = new TFile(fileName.c_str());

    TH1F*  hist = (TH1F*)file->Get(("sidebandcorWhist_"+observable).c_str());
    TH1F* ehist = (TH1F*)hist->Clone("ehist");

    vector<float> bins = selectBinning(observable,category);

    TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 4.0, "");
    frame->GetYaxis()->SetTitle("R_{sideband,W}");
    frame->GetXaxis()->SetTitle("Recoil [GeV]");
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
        err +=    hist->GetBinError  (i)*   hist->GetBinError  (i);
        err += pow(hist->GetBinContent(i)*0.10,2);
        ehist->SetBinError(i, sqrt(err));
    }

    frame->Draw();
    CMS_lumi(canvas, 4, 0, true);
    hist ->Draw("PE SAME");
    ehist->Draw("E2 SAME");
    hist ->Draw("PE SAME");

    canvas->RedrawAxis();
    TLegend* leg = new TLegend(0.3, 0.7, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist  , "R(sideband,W) Stat. Unc.","pl");
    leg->AddEntry(ehist , "Stat. + Syst. Uncertainty", "F");
    leg->Draw("SAME");

    canvas->SaveAs("rsidebandW.pdf");
    canvas->SaveAs("rsidebandW.png");

    //file->Close();
}




void plotTransferFactor(string fileName, int category, string observable, bool addtop = false, bool addsideband = false) {

  gROOT->SetBatch(kTRUE);

  //  rzmm(fileName,category,observable);
  rzee(fileName,category,observable);
  //  rwmn(fileName,category,observable);
  //  rwen(fileName,category,observable);
  //  rgam(fileName,category,observable);
  //  rzwj(fileName,category,observable);

  //  if(addtop){
    //    rtopmu(fileName,category,observable);
    //    rtopel(fileName,category,observable);
    //  }

  //  if(addsideband){
  //    rsidebandZ(fileName,category,observable);
  //    rsidebandW(fileName,category,observable);
  //  }

}


