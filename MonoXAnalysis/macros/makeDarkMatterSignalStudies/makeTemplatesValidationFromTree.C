#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static string luminosity = "35.9";

void makeTemplatesValidationFromTree(string inputTreeInterpolationName,  // tree used for the interpoaltion
				     string inputTemplateName, // template file
				     string outputDirectory, 
				     string interpolationPoint, //string with the mass point to extract the templated  
				     Category category, 
				     string medMass, string dmMass){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());
  
  TFile* inputTreeInterpFile = TFile::Open(inputTreeInterpolationName.c_str());
  TFile* inputFileTemplate   = TFile::Open(inputTemplateName.c_str());
  TH1F* histo_interpolation  = (TH1F*) inputFileTemplate->FindObjectAny(("signal_signal_"+interpolationPoint).c_str());
  
  TTree* tree = (TTree*) inputTreeInterpFile->Get("tree");
  TH1F* histo_fullSIM = (TH1F*) histo_interpolation->Clone("histo_fullSIM");
  histo_fullSIM->Reset();
  if(category == Category::monojet)
    tree->Draw("pfMetPt >> histo_fullSIM",("("+luminosity+"*weightPU*weightTurnOn*genWeight*xsec*(1./sumwgt))*(id == 1 && genMediatorMass == "+medMass+" && genX1Mass == "+dmMass+")").c_str(),"goff");
  else if(category == Category::monoV)
    tree->Draw("pfMetPt >> histo_fullSIM",("("+luminosity+"*weightPU*weightTurnOn*genWeight*xsec*(1./sumwgt))*(id == 2 && genMediatorMass == "+medMass+" && genX1Mass == "+dmMass+")").c_str(),"goff");
  
  // fix the overflow
  histo_fullSIM->SetBinContent(histo_fullSIM->GetNbinsX(),histo_fullSIM->GetBinContent(histo_fullSIM->GetNbinsX())+histo_fullSIM->GetBinContent(histo_fullSIM->GetNbinsX()+1));
  histo_fullSIM->SetBinError(histo_fullSIM->GetNbinsX(),sqrt(TMath::Power(histo_fullSIM->GetBinError(histo_fullSIM->GetNbinsX()),2)+
							     TMath::Power(histo_fullSIM->GetBinError(histo_fullSIM->GetNbinsX()+1),2)));
 
  
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.05);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  canvas->cd();
  histo_interpolation->Scale(1,"width");
  histo_fullSIM->Scale(1,"width");
  histo_interpolation->GetYaxis()->SetTitle("Events / GeV");
  histo_interpolation->GetYaxis()->SetTitleOffset(1.35);
  histo_interpolation->SetMarkerStyle(20);
  histo_interpolation->SetMarkerSize(1);
  histo_interpolation->SetMarkerColor(kBlack);
  histo_interpolation->SetLineColor(kBlack);
  histo_interpolation->GetXaxis()->SetLabelSize(0);
  histo_interpolation->GetXaxis()->SetTitleSize(0);
  histo_interpolation->Draw("EP");
 
  histo_fullSIM->SetLineColor(kRed);
  histo_fullSIM->SetLineWidth(2);

  TH1* histo_fullSIM_band = (TH1*) histo_fullSIM->Clone("histo_fullSIM_band");
  histo_fullSIM_band->SetFillColor(kGray);
  histo_fullSIM_band->Draw("E2 same");
  histo_fullSIM->Draw("hist same");
  histo_interpolation->Draw("EPsame");

  histo_interpolation->GetYaxis()->SetRangeUser(min(histo_interpolation->GetMinimum(),histo_fullSIM->GetMinimum())*0.01,
						max(histo_interpolation->GetMaximum(),histo_fullSIM->GetMaximum())*100);

  if(min(histo_interpolation->GetMinimum(),histo_fullSIM->GetMinimum()) == 0)
    histo_interpolation->GetYaxis()->SetRangeUser(0.0001,max(histo_interpolation->GetMaximum(),histo_fullSIM->GetMaximum())*100);

  TLegend leg (0.6,0.65,0.92,0.92);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.AddEntry((TObject*)(0),("m_{MED}="+medMass+" GeV "+"m_{DM}="+dmMass+" GeV").c_str(),"");
  leg.AddEntry(histo_interpolation,"Interpolated template","EP");
  leg.AddEntry(histo_fullSIM_band,"Full SIM template","FL");
  leg.Draw("same");

  CMS_lumi(canvas,luminosity);

  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  TH1* ratio = (TH1*) histo_fullSIM->Clone("ratio");
  ratio->Divide(histo_interpolation);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);
  ratio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ratio->GetYaxis()->SetTitle("FullSIM/Interp.");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.3,1.7);
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->GetYaxis()->SetTitleOffset(1.5);
  ratio->GetYaxis()->SetLabelSize(0.04);
  ratio->GetYaxis()->SetTitleSize(0.04);
  ratio->GetXaxis()->SetLabelSize(0.04);
  ratio->GetXaxis()->SetTitleSize(0.05);
  ratio->GetXaxis()->SetTitleOffset(1.1);
  ratio->Draw("EP");

  if(category == Category::monojet){
    canvas->SaveAs((outputDirectory+"/comparison_monojet_"+interpolationPoint+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/comparison_monojet_"+interpolationPoint+".pdf").c_str(),"pdf");
  }
  else if(category == Category::monoV){
    canvas->SaveAs((outputDirectory+"/comparison_monoV_"+interpolationPoint+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/comparison_monoV_"+interpolationPoint+".pdf").c_str(),"pdf");
  }

}
