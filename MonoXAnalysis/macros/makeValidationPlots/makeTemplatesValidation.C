#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

#include <vector>

void makeTemplatesValidation(string fileFullSIM, // file with templates from full sim
			     string fileInterpolation, // file with templates from interpolation
			     string interaction, //Vector, Axial, Scalar or Pseudoscalar
			     string signalType, // monoJ, monoW, monoZ
			     Category category, 
			     string outputDIR, 
			     double lumiScaleFullSIM = 1,
			     bool shapeComparison = false){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  setTDRStyle();

  TFile* inputFileFullSIM = TFile::Open(fileFullSIM.c_str());
  TFile* inputFileInterpolation = TFile::Open(fileInterpolation.c_str());

  system(("mkdir -p "+outputDIR).c_str());

  vector<string> medMassFullSIM;
  vector<string> dmMassFullSIM;
  vector<string> medMassInterpolation;
  vector<string> dmMassInterpolation;

  // mkae by had comparison list
  if(interaction == "Vector"){
    medMassFullSIM.push_back("100"); medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0500");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1200");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1200");
    medMassFullSIM.push_back("1500"); medMassInterpolation.push_back("1525");
    medMassFullSIM.push_back("1500"); medMassInterpolation.push_back("1525");
    dmMassFullSIM.push_back("1");    dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("200");  dmMassInterpolation.push_back("0200");
    dmMassFullSIM.push_back("350");  dmMassInterpolation.push_back("0300");
    dmMassFullSIM.push_back("1");    dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("400");  dmMassInterpolation.push_back("0400");
    dmMassFullSIM.push_back("10");   dmMassInterpolation.push_back("0010");
    dmMassFullSIM.push_back("200");  dmMassInterpolation.push_back("0200");
  }
  else if(interaction == "Axial"){
    medMassFullSIM.push_back("100"); medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("300"); medMassInterpolation.push_back("0325");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1200");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1200");
    medMassFullSIM.push_back("1500"); medMassInterpolation.push_back("1525");
    medMassFullSIM.push_back("1500"); medMassInterpolation.push_back("1525");

    dmMassFullSIM.push_back("50"); dmMassInterpolation.push_back("0050");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("10");  dmMassInterpolation.push_back("0010");
    dmMassFullSIM.push_back("350"); dmMassInterpolation.push_back("0300");
    dmMassFullSIM.push_back("50");  dmMassInterpolation.push_back("0050");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
			     
  }
  else if(interaction == "Scalar"){
    medMassFullSIM.push_back("50");  medMassInterpolation.push_back("0050");
    medMassFullSIM.push_back("100");  medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("100");  medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("300");  medMassInterpolation.push_back("0300");
    medMassFullSIM.push_back("500");  medMassInterpolation.push_back("0500");
    medMassFullSIM.push_back("500");  medMassInterpolation.push_back("0500");
    medMassFullSIM.push_back("1000"); medMassInterpolation.push_back("1000");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1250");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("50"); dmMassInterpolation.push_back("0050");
    dmMassFullSIM.push_back("100");  dmMassInterpolation.push_back("0100");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("1");  dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("200");  dmMassInterpolation.push_back("0200");
  }
  else if(interaction == "Pseudoscalar"){
    medMassFullSIM.push_back("100"); medMassInterpolation.push_back("0100");
    medMassFullSIM.push_back("300"); medMassInterpolation.push_back("0300");
    medMassFullSIM.push_back("500"); medMassInterpolation.push_back("0525");
    medMassFullSIM.push_back("1250"); medMassInterpolation.push_back("1200");
    medMassFullSIM.push_back("1500"); medMassInterpolation.push_back("1525");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("1"); dmMassInterpolation.push_back("0001");
    dmMassFullSIM.push_back("150"); dmMassInterpolation.push_back("0150");
    dmMassFullSIM.push_back("400"); dmMassInterpolation.push_back("0400");
    dmMassFullSIM.push_back("350"); dmMassInterpolation.push_back("0300");

  }

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  TPad *pad1 = new TPad("pad1","pad1",0,0.20,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.25);
  pad2->SetTickx();
  pad2->SetTicky();
  pad1->Draw();
  canvas->cd();
  pad2->SetGridy();
  pad2->Draw();


  string cat;
  if(category == Category::monojet)
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
  else if(interaction == "Pseudoscalar")
    model = "806";

  for(size_t ipoint = 0; ipoint < medMassFullSIM.size(); ipoint++){

    canvas->cd();
    pad1->cd();

  
    TH1F* histoFullSIM = (TH1F*) inputFileFullSIM->FindObjectAny((signalType+"hist_"+interaction+"_"+medMassFullSIM.at(ipoint)+"_"+dmMassFullSIM.at(ipoint)+"_met").c_str());   
    TH1F* histoInterpolation = (TH1F*) inputFileInterpolation->Get((cat+"/signal_signal_"+model+medMassInterpolation.at(ipoint)+dmMassInterpolation.at(ipoint)).c_str());
    cout<<signalType+"hist_"+interaction+"_"+medMassFullSIM.at(ipoint)+"_"+dmMassFullSIM.at(ipoint)+"_met"<<endl;
    cout<<histoFullSIM<<" "<<histoInterpolation<<" integral "<<histoFullSIM->Integral()<<" "<<histoInterpolation->Integral()<<endl;    
    histoFullSIM->SetTitle("");
    histoFullSIM->Scale(lumiScaleFullSIM);
    histoInterpolation->GetXaxis()->SetTitle("");    
    histoFullSIM->GetYaxis()->SetTitle("");
    if(shapeComparison)
      histoInterpolation->GetYaxis()->SetTitle("a.u.");
    else
      histoInterpolation->GetYaxis()->SetTitle("Entries");
      
    histoInterpolation->SetTitle("");    
    histoFullSIM->SetLineColor(kRed);
    histoFullSIM->SetLineWidth(2);

    histoInterpolation->SetLineColor(kBlack);
    histoInterpolation->SetMarkerColor(kBlack);
    histoInterpolation->SetMarkerSize(0.8);
    histoInterpolation->SetMarkerStyle(20);    
    if(shapeComparison){
      histoInterpolation->Scale(1./histoInterpolation->Integral());
      histoFullSIM->Scale(1./histoFullSIM->Integral());
    }

    histoInterpolation->GetYaxis()->SetRangeUser(max(0.00001,min(histoInterpolation->GetMinimum(),histoFullSIM->GetMinimum())*0.5),max(histoInterpolation->GetMaximum(),histoFullSIM->GetMaximum())*10);
    histoInterpolation->Draw("PE");

    CMS_lumi(pad1,"12.9",true);

    TH1* histoFullSIM_band = (TH1*) histoFullSIM->Clone("histoFullSIM_band");
    histoFullSIM_band->SetFillColor(kGray);
    histoFullSIM_band->Draw("E2 same");
    histoFullSIM->Draw("hist same");
    histoInterpolation->Draw("PEsame");


    TLegend* leg = new TLegend(0.5,0.70,0.90,0.85);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(histoFullSIM_band, Form("full SIM, m_{MED} = %s, m_{DM} = %s",medMassFullSIM.at(ipoint).c_str(),dmMassFullSIM.at(ipoint).c_str()),"FL");
    leg->AddEntry(histoInterpolation,Form("Interpolation, m_{MED} = %s, m_{DM} = %s",medMassInterpolation.at(ipoint).c_str(),dmMassInterpolation.at(ipoint).c_str()),"PEL");
    leg->Draw("same");
    pad1->RedrawAxis("sameaxis");
    pad1->SetLogy();

    canvas->cd();
    pad2->cd();

    TH1* frame2 = pad2->DrawFrame(histoFullSIM->GetXaxis()->GetXmin(),0.5,histoFullSIM->GetXaxis()->GetXmax(),1.5);
    frame2->GetXaxis()->SetLabelSize(0.10);
    frame2->GetXaxis()->SetLabelOffset(0.03);
    frame2->GetXaxis()->SetTitleSize(0.13);
    frame2->GetXaxis()->SetTitleOffset(1.05);
    frame2->GetYaxis()->SetLabelSize(0.08);
    frame2->GetYaxis()->SetTitleSize(0.10);
    frame2->GetXaxis()->SetTitle("MET [GeV]");
    frame2->GetYaxis()->SetTitle("Ratio");
    frame2->GetYaxis()->SetTitleOffset(0.5);
    frame2->GetYaxis()->SetNdivisions(504, false);
    frame2->Draw();
    
    TH1F* ratio = (TH1F*) histoInterpolation->Clone("ratio");
    ratio->Divide(histoFullSIM);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerSize(0.8);
    ratio->SetMarkerStyle(20);
    
    ratio->Draw("PEsame");

    pad2->RedrawAxis("sameaxis");
    
    canvas->SaveAs((outputDIR+"/template_"+interaction+"_"+cat+"_"+string(medMassFullSIM.at(ipoint))+"_"+string(dmMassFullSIM.at(ipoint))+".pdf").c_str(),"pdf");
    canvas->SaveAs((outputDIR+"/template_"+interaction+"_"+cat+"_"+string(medMassFullSIM.at(ipoint))+"_"+string(dmMassFullSIM.at(ipoint))+".png").c_str(),"png");
  }
  
  return;
  
}
