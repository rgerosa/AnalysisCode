#include "CMS_lumi.h"
#include "histoUtils.h"

void drawCorrelationPlots(string inputFileName, string controlRegion, vector<string> observables2D, vector<string> observablesLatex2D_X, vector<string> observablesLatex2D_Y,
			  string mediatorMass = "1000", string DMMass = "50",bool isLog = true){


  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TFile* templateFile = TFile::Open(inputFileName.c_str());

  vector<TH2*> datahist;
  vector<TH2*> bkghist;
  vector<TH2*> monoJhist;
  vector<TH2*> monoWhist;
  vector<TH2*> monoZhist;
  
  for(auto obs : observables2D){
    if(controlRegion != "SR"){
      datahist.push_back((TH2*) templateFile->Get(("datahist"+controlRegion+"_2D_"+obs).c_str()));

      if(controlRegion == "topmu" or controlRegion == "topel"){
	bkghist.push_back((TH2*) templateFile->Get(("tbkghist"+controlRegion+"_2D_"+obs).c_str()));
	bkghist.back()->Add((TH2*) templateFile->Get(("dbkghist"+controlRegion+"_2D_"+obs).c_str()));
	bkghist.back()->Add((TH2*) templateFile->Get(("qbkghist"+controlRegion+"_2D_"+obs).c_str()));
	bkghist.back()->Add((TH2*) templateFile->Get(("vlbkghist"+controlRegion+"_2D_"+obs).c_str()));
	bkghist.back()->Add((TH2*) templateFile->Get(("vllbkghist"+controlRegion+"_2D_"+obs).c_str()));
      }
      else if(controlRegion == "wmn" or controlRegion == "wen"){
	bkghist.push_back((TH2*) templateFile->Get(("vlbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("dbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("qbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("tbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("vllbkghist"+controlRegion+"_2D_"+obs).c_str()));	
      }
      else if(controlRegion == "zmm" or controlRegion == "zee"){
	bkghist.push_back((TH2*) templateFile->Get(("vllbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("dbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("qbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("tbkghist"+controlRegion+"_2D_"+obs).c_str()));
        bkghist.back()->Add((TH2*) templateFile->Get(("vlbkghist"+controlRegion+"_2D_"+obs).c_str()));
      }
    }
    if(controlRegion == "SR"){
      monoJhist.push_back((TH2*) templateFile->Get(("monoJhist_"+mediatorMass+"_"+DMMass+"_2D_"+obs).c_str()));
      monoWhist.push_back((TH2*) templateFile->Get(("monoWhist_"+mediatorMass+"_"+DMMass+"_2D_"+obs).c_str()));
      monoZhist.push_back((TH2*) templateFile->Get(("monoZhist_"+mediatorMass+"_"+DMMass+"_2D_"+obs).c_str()));
    }
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.10);

  if(isLog)
    canvas->SetLogz();

  for(auto hist : datahist)
    hist->Scale(1.0,"width");
  for(auto hist : bkghist)
    hist->Scale(1.0,"width");
  for(auto hist : monoJhist)
    hist->Scale(1.0,"width");
  for(auto hist : monoWhist)
    hist->Scale(1.0,"width");
  for(auto hist : monoZhist)
    hist->Scale(1.0,"width");

  int iObs = 0;
  for(auto obs : observables2D){    

    bin2D bins = selectBinning2D (obs,0);
    TH2F* frame = new TH2F("","",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    frame->GetXaxis()->SetTitle(observablesLatex2D_X.at(iObs).c_str());
    frame->GetYaxis()->SetTitle(observablesLatex2D_Y.at(iObs).c_str());
    frame->GetZaxis()->SetTitle("Events/GeV");
    frame->GetZaxis()->SetTitleOffset(0.1);
    frame->GetZaxis()->CenterTitle();
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.045);

    if(controlRegion == "SR"){ // just plot correlation matrxi + profile
      frame->Draw();
      CMS_lumi(canvas,4, 0, true);
      TProfile* profileX = monoJhist.at(iObs)->ProfileX("_pfx",1,monoJhist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      TGraph2D* graph2D_1 = new TGraph2D (monoJhist.at(iObs));
      for(int iBinX = 0; iBinX < monoJhist.at(iObs)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < monoJhist.at(iObs)->GetNbinsY(); iBinY++){
	  graph2D_1->Interpolate(monoJhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1),monoJhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
	}
      }
	

      graph2D_1->GetHistogram()->Draw("colz same");
      profileX->Draw("EP same");

      TString corr = Form("#rho = %0.2f",monoJhist.at(iObs)->GetCorrelationFactor());
      TLatex* correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");

      canvas->RedrawAxis("sameaxis");
      canvas->SaveAs(("monoJ_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoJ_"+obs+".png").c_str());

      frame->Draw();
      CMS_lumi(canvas,4, 0, true);
      profileX = monoWhist.at(iObs)->ProfileX("_pfx",1,monoWhist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);

      TGraph2D* graph2D_2 = new TGraph2D (monoWhist.at(iObs));
      for(int iBinX = 0; iBinX < monoWhist.at(iObs)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < monoWhist.at(iObs)->GetNbinsY(); iBinY++){
	  graph2D_2->Interpolate(monoWhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1),monoWhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
	}
      }

      graph2D_2->GetHistogram()->Draw("colz same");
      profileX->Draw("EP same");
      corr = Form("#rho = %0.2f",monoWhist.at(iObs)->GetCorrelationFactor());
      correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs(("monoW_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoW_"+obs+".png").c_str());

      frame->Draw();
      CMS_lumi(canvas,4, 0, true);
      profileX = monoZhist.at(iObs)->ProfileX("_pfx",1,monoZhist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      TGraph2D* graph2D_3 = new TGraph2D (monoZhist.at(iObs));
      for(int iBinX = 0; iBinX < monoZhist.at(iObs)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < monoZhist.at(iObs)->GetNbinsY(); iBinY++){
	  graph2D_3->Interpolate(monoZhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1),monoZhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
	}
      }

      graph2D_3->GetHistogram()->Draw("colz same");
      profileX->Draw("EP same");
      corr = Form("#rho = %0.2f",monoZhist.at(iObs)->GetCorrelationFactor());
      correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs(("monoZ_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoW_"+obs+".png").c_str());
    }
    else{

      // plot MC
      frame->Draw();
      CMS_lumi(canvas,4, 0, true);
      TGraph2D* graph2D_1  = new TGraph2D(bkghist.at(iObs));
      for(int iBinX = 0; iBinX < bkghist.at(iObs)->GetNbinsX(); iBinX++){
        for(int iBinY = 0; iBinY < bkghist.at(iObs)->GetNbinsY(); iBinY++){
          graph2D_1->Interpolate(bkghist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1),bkghist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
        }
      }
      graph2D_1->GetHistogram()->Draw("colz same");
      TProfile* profileX = bkghist.at(iObs)->ProfileX("_pfx",1,bkghist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      profileX->Draw("EP same");
      TString corr = Form("#rho = %0.2f",bkghist.at(iObs)->GetCorrelationFactor());
      TLatex* correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs((controlRegion+"_"+obs+"_MC.pdf").c_str());
      canvas->SaveAs((controlRegion+"_"+obs+"_MC.png").c_str());
 
      // plot data
      frame->Draw();
      CMS_lumi(canvas,4, 0, true);
      TGraph2D* graph2D_2 = new TGraph2D (datahist.at(iObs));
      for(int iBinX = 0; iBinX < datahist.at(iObs)->GetNbinsX(); iBinX++){
        for(int iBinY = 0; iBinY < datahist.at(iObs)->GetNbinsY(); iBinY++){
          graph2D_2->Interpolate(datahist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1),datahist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
        }
      }
      graph2D_2->GetHistogram()->Draw("colz same");
      profileX = datahist.at(iObs)->ProfileX("_pfx",1,datahist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      profileX->Draw("EP same");
      corr = Form("#rho = %0.2f",datahist.at(iObs)->GetCorrelationFactor());
      correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs((controlRegion+"_"+obs+"_Data.pdf").c_str());
      canvas->SaveAs((controlRegion+"_"+obs+"_Data.png").c_str());

    }      
    iObs++;
  }
}
