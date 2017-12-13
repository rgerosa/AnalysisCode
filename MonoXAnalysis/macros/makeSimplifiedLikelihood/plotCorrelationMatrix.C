static bool debug = false;
static bool skipCorrelations = false;

#include "../makeTemplates/histoUtils.h"
#include "simplifiedLikelihoodUtils.h"
#include "../CMS_lumi.h"

static bool addPreliminary = false;

void plotCorrelationMatrix(string inputFile, Category category, bool isZeynep, string outputDIR, bool addLabel = false, bool suppressDiagonal = false){

  system(("mkdir -p "+outputDIR).c_str());
  initializeBinning();
  gROOT->SetBatch(kTRUE);
  setTDRStyle();  
  //gStyle->SetPalette(kBird);
  gStyle->SetPalette(kThermometer);
  gStyle->SetNumberContours(999);

  TFile *file = new TFile(inputFile.c_str(),"READ");
  string dir = "ch1_ch1";
  if(category == Category::monoV and not isZeynep)
    dir = "ch2_ch1";
  else if(category == Category::monojet and isZeynep)
    dir = "monojet_signal";
  else if(category == Category::monoV and isZeynep)
    dir = "monov_signal";
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    dir = "ch1";
    
    
  TCanvas* canvas = NULL;
  if(category == Category::monoV)
    canvas = new TCanvas("canvas","",900,900);
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    canvas = new TCanvas("canvas","",900,900);
  else if(category == Category::monojet)
    canvas = new TCanvas("canvas","",1300,900);
  else if(category == Category::total)
    canvas = new TCanvas("canvas","",1600,1000);
  
  canvas->cd();  
  canvas->SetGrid();
  canvas->SetRightMargin(0.12);
  if(addLabel){
    canvas->SetBottomMargin(0.24);
    canvas->SetLeftMargin(0.2);
  }

  TH1F* bkg   = NULL;
  TH2F* covar = NULL;
  TH2F* covar_total = NULL;

  if(category == Category::total){
    covar_total = (TH2F*) file->Get("shapes_fit_b/overall_total_covar");
  }
  else{
    bkg   = (TH1F*)file->Get(("shapes_fit_b/"+dir+"/total_background").c_str());
    covar = (TH2F*)file->Get(("shapes_fit_b/"+dir+"/total_covar").c_str());
  }

  TH1F* bkg_monojet = NULL;
  TH1F* bkg_monov   = NULL;
  if(category == Category::total and not isZeynep){
    bkg_monojet = (TH1F*) file->Get("shapes_fit_b/ch1_ch1/total_background");
    bkg_monov   = (TH1F*) file->Get("shapes_fit_b/ch1_ch2/total_background");
  }
  else if(category == Category::total and isZeynep){
    bkg_monojet = (TH1F*) file->Get("shapes_fit_b/monojet_signal/total_background");
    bkg_monov   = (TH1F*) file->Get("shapes_fit_b/monov_signal/total_background");
  }

  if(category == Category::total){
    // merge for the background total
    bkg = new TH1F("bkg_total","",bkg_monojet->GetNbinsX()+bkg_monov->GetNbinsX(),0,bkg_monojet->GetNbinsX()+bkg_monov->GetNbinsX()+1);
    bkg->Sumw2();
    
    int iBin = 0;
    for(; iBin < bkg_monojet->GetNbinsX(); iBin++){
      bkg->SetBinContent(iBin+1,bkg_monojet->GetBinContent(iBin+1));
      bkg->SetBinError(iBin+1,bkg_monojet->GetBinError(iBin+1));
    }
    for(int iBinX = 0; iBinX < bkg_monov->GetNbinsX(); iBinX++,iBin++){
      bkg->SetBinContent(iBin+1,bkg_monov->GetBinContent(iBinX+1));
      bkg->SetBinError(iBin+1,bkg_monojet->GetBinError(iBinX+1));
    }
    // merge for the total covariance
    if(not isZeynep)
      covar = importCorrelationMatrix(covar_total,"ch1_ch1","ch2_ch1");
    else
      covar = importCorrelationMatrix(covar_total,"monojet_signal","monov_signal");	
  }
  

  TH2F* corr = (TH2F*) covar->Clone("test");
  corr->Reset();
  int nbins  = bkg->GetNbinsX();
  // loop on the bin --> square matrix --> multiply by the bin width to get real yields
  for (int b=1 ; b<=nbins; b++){
    double bw = bkg->GetBinWidth(b);
    for (int j=1 ; j<=nbins; j++){
      double bj = bkg->GetBinWidth(j);
      covar->SetBinContent(b,j,covar->GetBinContent(b,j)*bw*bj);
    }
  }

  for (int b=1;b<=nbins;b++){
    for (int j=1;j<=nbins;j++){
      // calculate the standard deviation (diagonal term)
      double sigb = TMath::Sqrt(covar->GetBinContent(b,b));
      double sigj = TMath::Sqrt(covar->GetBinContent(j,j));
      corr->SetBinContent(b,j,covar->GetBinContent(b,j)/(sigb*sigj));      
      
      if(addLabel){
	if(category != Category::total and category != Category::VBF and category != Category::VBFrelaxed){
	  if(b == 1)
	    corr->GetYaxis()->SetBinLabel(j,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(j)),int(bkg->GetXaxis()->GetBinLowEdge(j+1))));
	  if(j == 1) 
	    corr->GetXaxis()->SetBinLabel(b,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(b)),int(bkg->GetXaxis()->GetBinLowEdge(b+1))));
	}
	else if(category == Category::total and category != Category::VBF and category != Category::VBFrelaxed){
	  if(b == 1 and j <= bkg_monojet->GetNbinsX())
	    corr->GetYaxis()->SetBinLabel(j,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg_monojet->GetXaxis()->GetBinLowEdge(j)),int(bkg_monojet->GetXaxis()->GetBinLowEdge(j+1))));
	  if(b == 1 and j > bkg_monojet->GetNbinsX())
	    corr->GetYaxis()->SetBinLabel(j,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg_monov->GetXaxis()->GetBinLowEdge(j-bkg_monojet->GetNbinsX())),int(bkg_monov->GetXaxis()->GetBinLowEdge(j+1-bkg_monojet->GetNbinsX()))));
	  if(j == 1 and b <= bkg_monojet->GetNbinsX())
	    corr->GetXaxis()->SetBinLabel(b,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg_monojet->GetXaxis()->GetBinLowEdge(b)),int(bkg_monojet->GetXaxis()->GetBinLowEdge(b+1))));
	  if(j == 1 and b > bkg_monojet->GetNbinsX())
	    corr->GetXaxis()->SetBinLabel(b,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg_monov->GetXaxis()->GetBinLowEdge(b-bkg_monojet->GetNbinsX())),int(bkg_monov->GetXaxis()->GetBinLowEdge(b+1-bkg_monojet->GetNbinsX()))));
	}
	else if(category != Category::total and (category == Category::VBF or category == Category::VBFrelaxed)){
	  if(b == 1)
	    corr->GetYaxis()->SetBinLabel(j,Form("%d < m_{jj} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(j)),int(bkg->GetXaxis()->GetBinLowEdge(j+1))));
	  if(j == 1) 
	    corr->GetXaxis()->SetBinLabel(b,Form("%d < m_{jj} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(b)),int(bkg->GetXaxis()->GetBinLowEdge(b+1))));
	}
      }
      else{
	if(b == 1)
	  corr->GetYaxis()->SetBinLabel(j,Form("%d",int(bkg->GetXaxis()->GetBinLowEdge(j))));
	if(j == 1)
	  corr->GetXaxis()->SetBinLabel(b,Form("%d",int(bkg->GetXaxis()->GetBinLowEdge(b))));
	corr->GetXaxis()->CenterLabels(kFALSE);
	corr->GetYaxis()->CenterLabels(kFALSE);
      }
    }
  }

  if(corr->GetMinimum() >= 0)
    corr->GetZaxis()->SetRangeUser(max(-1.,corr->GetMinimum()*0.9),min(corr->GetMaximum()*1.1,1.));
  else
    corr->GetZaxis()->SetRangeUser(max(-1.,corr->GetMinimum()*1.1),min(corr->GetMaximum()*1.1,1.));

  if(suppressDiagonal){
    for(int i = 0; i <= corr->GetNbinsX(); i++){
      for(int j = 0; j <= corr->GetNbinsY(); j++){
	if(j < i)
	  corr->SetBinContent(i,j,-100);
      }
    }
  }

  corr->GetXaxis()->SetTitle("");
  corr->GetYaxis()->SetTitle("");
  corr->GetZaxis()->SetTitle("Correlation");
  if(addLabel)
    corr->GetXaxis()->LabelsOption("v");
  corr->GetYaxis()->LabelsOption("v");
  corr->GetXaxis()->SetLabelSize(0.027);
  corr->GetYaxis()->SetLabelSize(0.027);  
  corr->GetZaxis()->SetLabelSize(0.029);
  corr->GetZaxis()->SetTitleSize(0.035);
  corr->GetZaxis()->SetTitleOffset(1.1);

  gStyle->SetPaintTextFormat("1.2f");
  corr->Draw("COLZTEXT");
  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",true,true,true,0.05,-0.06);
  else
    CMS_lumi(canvas,"35.9",true,false,true,0.05,-0.06);

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)corr->GetListOfFunctions()->FindObject("palette");
  
  // the following lines moe the paletter. Choose the values you need for the position.
  palette->SetX1NDC(0.885);
  palette->SetX2NDC(0.92);
  gPad->Modified();
  gPad->Update();

  if(category == Category::monojet){
    corr->SetMarkerColor(kBlack);
    canvas->SaveAs((outputDIR+"/correlation_monojet.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monojet.png").c_str());
    //    canvas->SetLogz();
    //    canvas->SaveAs((outputDIR+"/correlation_monojet_log.pdf").c_str());
    //    canvas->SaveAs((outputDIR+"/correlation_monojet_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    corr->SetMarkerColor(kBlack);
    CMS_lumi(canvas,"35.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_monojet_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monojet_text.png").c_str());
    TFile* outputFile = new TFile((outputDIR+"/correlation_monojet.root").c_str(),"RECREATE");
    outputFile->cd();
    corr->Write("correlation_monojet");
    outputFile->Close();
  }
  else if(category == Category::monoV){
    corr->SetMarkerColor(kBlack);
    canvas->SaveAs((outputDIR+"/correlation_monov.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monov.png").c_str());
    //    canvas->SetLogz();
    //    canvas->SaveAs((outputDIR+"/correlation_monov_log.pdf").c_str());
    //    canvas->SaveAs((outputDIR+"/correlation_monov_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    corr->SetMarkerColor(kBlack);
    CMS_lumi(canvas,"35.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_monov_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monov_text.png").c_str());
    TFile* outputFile = new TFile((outputDIR+"/correlation_monov.root").c_str(),"RECREATE");
    outputFile->cd();
    corr->Write("correlation_monov");
    outputFile->Close();
  }
  else if(category == Category::VBF or category == Category::VBFrelaxed){
    corr->SetMarkerColor(kBlack);
    canvas->SaveAs((outputDIR+"/correlation_VBF.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_VBF.png").c_str());
    //    canvas->SetLogz();
    //    canvas->SaveAs((outputDIR+"/correlation_monov_log.pdf").c_str());
    //    canvas->SaveAs((outputDIR+"/correlation_monov_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    corr->SetMarkerColor(kBlack);
    CMS_lumi(canvas,"35.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_VBF_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_VBF_text.png").c_str());
    TFile* outputFile = new TFile((outputDIR+"/correlation_monov.root").c_str(),"RECREATE");
    outputFile->cd();
    corr->Write("correlation_VBF");
    outputFile->Close();
  }
  else if(category == Category::total){
    corr->SetMarkerColor(kBlack);
    canvas->SaveAs((outputDIR+"/correlation_total.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_total.png").c_str());
    //    canvas->SetLogz();
    //    canvas->SaveAs((outputDIR+"/correlation_total_log.pdf").c_str());
    //    canvas->SaveAs((outputDIR+"/correlation_total_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    corr->SetMarkerColor(kBlack);
    CMS_lumi(canvas,"35.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_total_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_total_text.png").c_str());
    canvas->SaveAs((outputDIR+"/correlation_total_text.root").c_str());
    TFile* outputFile = new TFile((outputDIR+"/correlation_total.root").c_str(),"RECREATE");
    outputFile->cd();
    corr->Write("correlation_total");
    outputFile->Close();
  }
}
