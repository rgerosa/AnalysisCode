#include "../CMS_lumi.h"

void makeDMSimpDMFComparison(string file_DMF,
			     string file_DMSimp,
			     string outputDIR,
			     string postfix){

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TFile* file_dmf    = TFile::Open(file_DMF.c_str());
  TFile* file_dmsimp = TFile::Open(file_DMSimp.c_str());

  TGraph* monojet_dmf = (TGraph*) file_dmf->Get("monojet_xsec");
  TGraph* monoW_dmf = (TGraph*) file_dmf->Get("monoW_xsec");
  TGraph* monoZ_dmf = (TGraph*) file_dmf->Get("monoZ_xsec");
  
  TGraph* monojet_dmsimp = (TGraph*) file_dmsimp->Get("monojet_xsec");
  TGraph* monoW_dmsimp = (TGraph*) file_dmsimp->Get("monoW_xsec");
  TGraph* monoZ_dmsimp = (TGraph*) file_dmsimp->Get("monoZ_xsec");

  double min_monojet = min(TMath::MinElement(monojet_dmf->GetN(),monojet_dmf->GetX()),TMath::MinElement(monojet_dmsimp->GetN(),monojet_dmsimp->GetX()));
  double max_monojet = max(TMath::MaxElement(monojet_dmf->GetN(),monojet_dmf->GetX()),TMath::MaxElement(monojet_dmsimp->GetN(),monojet_dmsimp->GetX()));

  double min_monoW = 0;
  double max_monoW = 0;
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    min(TMath::MinElement(monoW_dmf->GetN(),monoW_dmf->GetX()),TMath::MinElement(monoW_dmsimp->GetN(),monoW_dmsimp->GetX()));
    max(TMath::MaxElement(monoW_dmf->GetN(),monoW_dmf->GetX()),TMath::MaxElement(monoW_dmsimp->GetN(),monoW_dmsimp->GetX()));
  }

  double min_monoZ = 0;
  double max_monoZ = 0;
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    min(TMath::MinElement(monoZ_dmf->GetN(),monoZ_dmf->GetX()),TMath::MinElement(monoZ_dmsimp->GetN(),monoZ_dmsimp->GetX()));
    max(TMath::MaxElement(monoZ_dmf->GetN(),monoZ_dmf->GetX()),TMath::MaxElement(monoZ_dmsimp->GetN(),monoZ_dmsimp->GetX()));
  }
  
  double minimum = min_monojet;
  double maximum = max_monojet;
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    minimum = min(min_monojet,min(min_monoW,min_monoW));
    maximum = max(max_monojet,max(max_monoW,max_monoW));
  }

  TGraph* ratio_monojet = new TGraph();
  int ratioPoint = 0;
  for(int iPoint = 0; iPoint < monojet_dmf->GetN(); iPoint++){
    double x1,y1;
    monojet_dmf->GetPoint(iPoint,x1,y1);
    for(int jPoint = 0; jPoint < monojet_dmsimp->GetN(); jPoint++){
      double x2,y2;
      monojet_dmsimp->GetPoint(jPoint,x2,y2);
      if(x1 == x2){
	cout<<"Monojet common point "<<x1<<" DMF "<<y1<<" DMSimp "<<y2<<" Ratio: "<<y2/y1<<endl;
	ratio_monojet->SetPoint(ratioPoint,x1,y2/y1);
	ratioPoint++;
      }
    }
  }

  TGraph* ratio_monoW = new TGraph();
  ratioPoint = 0;
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    for(int iPoint = 0; iPoint < monoW_dmf->GetN(); iPoint++){
      double x1,y1;
      monoW_dmf->GetPoint(iPoint,x1,y1);
      for(int jPoint = 0; jPoint < monoW_dmsimp->GetN(); jPoint++){
	double x2,y2;
	monoW_dmsimp->GetPoint(jPoint,x2,y2);
	if(x1 == x2){
	  cout<<"MonoW common point "<<x1<<" DMF "<<y1<<" DMSimp "<<y2<<" Ratio: "<<y2/y1<<endl;
	  ratio_monoW->SetPoint(ratioPoint,x1,y2/y1);
	  ratioPoint++;
      }
      }
    }
  }
  
  TGraph* ratio_monoZ = new TGraph();
  ratioPoint = 0;
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    for(int iPoint = 0; iPoint < monoZ_dmf->GetN(); iPoint++){
      double x1,y1;
      monoZ_dmf->GetPoint(iPoint,x1,y1);
      for(int jPoint = 0; jPoint < monoZ_dmsimp->GetN(); jPoint++){
	double x2,y2;
	monoZ_dmsimp->GetPoint(jPoint,x2,y2);
	if(x1 == x2){
	  cout<<"MonoZ common point "<<x1<<" DMF "<<y1<<" DMSimp "<<y2<<" Ratio: "<<y2/y1<<endl;
	  ratio_monoZ->SetPoint(ratioPoint,x1,y2/y1);
	  ratioPoint++;
	}
      }
    }
  }

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  TH1F* frame = new TH1F("frame","",100,minimum,maximum); 
  frame->GetXaxis()->SetTitle("Mediator mass [GeV]");
  frame->GetYaxis()->SetTitle("DMSimp/DMF");
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(2);
  frame->Draw("");

  if(postfix != "scalar" and postfix !="pseudoscalar")
    frame->GetYaxis()->SetRangeUser(min(TMath::MinElement(ratio_monojet->GetN(),ratio_monojet->GetY()),min(TMath::MinElement(ratio_monoW->GetN(),ratio_monoW->GetY()),TMath::MinElement(ratio_monoZ->GetN(),ratio_monoZ->GetY())))*0.75,max(TMath::MaxElement(ratio_monojet->GetN(),ratio_monojet->GetY()),max(TMath::MaxElement(ratio_monoW->GetN(),ratio_monoW->GetY()),TMath::MaxElement(ratio_monoZ->GetN(),ratio_monoZ->GetY())))*1.5);  
  else
    frame->GetYaxis()->SetRangeUser(TMath::MinElement(ratio_monojet->GetN(),ratio_monojet->GetY())*0.75,TMath::MaxElement(ratio_monojet->GetN(),ratio_monojet->GetY())*1.5);
  
  CMS_lumi(canvas,"");

  ratio_monojet->SetMarkerColor(kBlack);
  ratio_monojet->SetLineColor(kBlack);
  ratio_monojet->SetLineWidth(2);
  ratio_monojet->SetMarkerStyle(20);
  ratio_monojet->Draw("EPLsame");
  ratio_monoW->SetMarkerColor(kRed);
  ratio_monoW->SetLineColor(kRed);
  ratio_monoW->SetLineWidth(2);
  ratio_monoW->SetMarkerStyle(20);
  if(postfix != "scalar" and postfix !="pseudoscalar")
    ratio_monoW->Draw("EPLsame");
  ratio_monoZ->SetMarkerColor(kBlue);
  ratio_monoZ->SetMarkerStyle(20);
  ratio_monoZ->SetLineColor(kBlue);
  ratio_monoZ->SetLineWidth(2);
  if(postfix != "scalar" and postfix !="pseudoscalar")
    ratio_monoZ->Draw("EPLsame");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry(ratio_monojet,"mono-jet, m_{DM}=1 GeV","PL");
  if(postfix != "scalar" and postfix !="pseudoscalar"){
    leg.AddEntry(ratio_monoW,"mono-W, m_{DM}=1 GeV","PL");
    leg.AddEntry(ratio_monoZ,"mono-Z, m_{DM}=1 GeV","PL");
  }
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".pdf").c_str(),"pdf");
}
