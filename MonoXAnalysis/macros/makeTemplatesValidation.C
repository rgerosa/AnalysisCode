#include "CMS_lumi.h"

#include <vector>

void makeTemplatesValidation(string fileFullSIM, string fileInterpolation, string interaction, string signalType, int category, string outputDIR){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TFile* inputFileFullSIM = TFile::Open(fileFullSIM.c_str());
  TFile* inputFileInterpolation = TFile::Open(fileInterpolation.c_str());

  system(("mkdir -p "+outputDIR).c_str());

  vector<string> medMassFullSIM;
  vector<string> dmMassFullSIM;
  vector<string> medMassInterpolation;
  vector<string> dmMassInterpolation;

  if(interaction == "Vector"){
    medMassFullSIM.push_back("200"); medMassInterpolation.push_back("0200");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("150"); dmMassInterpolation.push_back("0150");
    dmMassFullSIM.push_back("50");  dmMassInterpolation.push_back("0050");
    dmMassFullSIM.push_back("150"); dmMassInterpolation.push_back("0150");
  }
  else if(interaction == "Axial"){
    medMassFullSIM.push_back("200"); medMassInterpolation.push_back("0200");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("2000"); medMassInterpolation.push_back("2000");

    dmMassFullSIM.push_back("10"); dmMassInterpolation.push_back("0010");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("150"); dmMassInterpolation.push_back("0150");
    dmMassFullSIM.push_back("10");  dmMassInterpolation.push_back("0010");
    dmMassFullSIM.push_back("100"); dmMassInterpolation.push_back("0100");
			     
  }
  else if(interaction == "Scalar"){
    medMassFullSIM.push_back("100");  medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("2000"); medMassInterpolation.push_back("2000");
    dmMassFullSIM.push_back("50"); dmMassInterpolation.push_back("0050");
    dmMassFullSIM.push_back("10"); dmMassInterpolation.push_back("0010");
    dmMassFullSIM.push_back("150"); dmMassInterpolation.push_back("0150");
    dmMassFullSIM.push_back("10");  dmMassInterpolation.push_back("0010");
  }

  TCanvas* canvas = new TCanvas("canvas","",500,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.27);
  pad2->SetTickx();
  pad2->SetTicky();

  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();

  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();

  string cat;
  if(category <= 1)
    cat = "category_monojet";
  else
    cat = "category_monov";

  string model;
  if(interaction == "Vector")
    model = "800";
  else if(interaction == "Axial")
    model = "801";
  else if(interaction == "Scalar")
    model = "805";
  else if(interaction == "Pseduoscalar")
    model = "806";

  for(size_t ipoint = 0; ipoint < medMassFullSIM.size(); ipoint++){

    canvas->cd();
    pad1->cd();

   
    TH1F* histoFullSIM = (TH1F*) inputFileFullSIM->Get((signalType+"hist_"+interaction+"_"+medMassFullSIM.at(ipoint)+"_"+dmMassFullSIM.at(ipoint)+"_met").c_str());      
    TH1F* histoInterpolation = (TH1F*) inputFileInterpolation->Get((cat+"/signal_signal_"+model+medMassInterpolation.at(ipoint)+dmMassInterpolation.at(ipoint)).c_str());
    
    histoFullSIM->SetTitle("");
    histoInterpolation->SetTitle("");
    
    histoFullSIM->SetLineColor(kBlack);
    histoFullSIM->SetLineWidth(2);

    histoInterpolation->SetLineColor(kRed);
    histoInterpolation->SetMarkerColor(kRed);
    histoInterpolation->SetMarkerSize(0.8);
    histoInterpolation->SetMarkerStyle(20);    
    histoInterpolation->Draw("PE");

    CMS_lumi(pad1,"2.30",true);

    TH1* histoFullSIM_band = (TH1*) histoFullSIM->Clone("histoFullSIM_band");
    histoFullSIM_band->SetFillColor(kBlack);
    histoFullSIM_band->SetFillStyle(3001);
    histoFullSIM_band->Draw("E2 same");
    histoFullSIM->Draw("hist same");


    TLegend* leg = new TLegend(0.35,0.6,0.90,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(histoFullSIM_band,Form("full SIM, g_{SM} = 1, m_{MED}=%s, m_{DM}=%s",medMassFullSIM.at(ipoint).c_str(),dmMassFullSIM.at(ipoint).c_str()),"FL");
    leg->AddEntry(histoInterpolation,Form("Interpolation"),"P");
    leg->Draw("same");
    pad1->RedrawAxis("sameaxis");
    pad1->SetLogy();

    canvas->cd();

    pad2->cd();

    TH1* frame2 = pad2->DrawFrame(histoFullSIM->GetXaxis()->GetXmin(),0.,histoFullSIM->GetXaxis()->GetXmax(),2.0);
    frame2->GetXaxis()->SetLabelSize(0.10);
    frame2->GetXaxis()->SetLabelOffset(0.03);
    frame2->GetXaxis()->SetTitleSize(0.13);
    frame2->GetXaxis()->SetTitleOffset(1.05);
    frame2->GetYaxis()->SetLabelSize(0.08);
    frame2->GetYaxis()->SetTitleSize(0.10);
    frame2->GetXaxis()->SetTitle("MET [GeV]");
    frame2->GetYaxis()->SetTitle("FullSIM/Interpolation");
    frame2->GetYaxis()->SetTitleOffset(0.5);
    frame2->GetYaxis()->SetNdivisions(504, false);
    frame2->Draw();
    
    TH1F* ratio = (TH1F*) histoFullSIM->Clone("ratio");
    ratio->Divide(histoInterpolation);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerSize(0.8);
    ratio->SetMarkerStyle(20);

    ratio->Draw("PEsame");

    pad2->RedrawAxis("sameaxis");
    
    canvas->SaveAs((outputDIR+"/template_"+interaction+"_"+cat+"_"+string(medMassFullSIM.at(ipoint))+"_"+string(dmMassFullSIM.at(ipoint))+".pdf").c_str(),"pdf");
    canvas->SaveAs((outputDIR+"/template_"+interaction+"_"+cat+"_"+string(medMassFullSIM.at(ipoint))+"_"+string(dmMassFullSIM.at(ipoint))+".png").c_str(),"png");
  }
  
  return;
  
}
