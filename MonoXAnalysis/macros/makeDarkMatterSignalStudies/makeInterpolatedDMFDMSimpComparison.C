#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float luminosity = 35.9;

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) <= 0)
      histo->SetBinContent(iBin+1,0.0001);
  }
}


void makeInterpolatedDMFDMSimpComparison(string inputTemplateDMFName, string inputTemplateDMSimpName, string outputDirectory, string interpolationPoint, Category category){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());
  
  TFile* inputTemplateDMF = TFile::Open(inputTemplateDMFName.c_str());
  TFile* inputTemplateDMSimp = TFile::Open(inputTemplateDMSimpName.c_str());

  TH1F* histo_interpolation_DMF    = (TH1F*) inputTemplateDMF->FindObjectAny(("signal_signal_"+interpolationPoint).c_str());
  TH1F* histo_interpolation_DMSimp = (TH1F*) inputTemplateDMSimp->FindObjectAny(("signal_signal_"+interpolationPoint).c_str());
  
  // compare correctly the luminosity
  histo_interpolation_DMF->Scale(luminosity/12.9);

  checkNegativeBin(histo_interpolation_DMF);
  checkNegativeBin(histo_interpolation_DMSimp);
  
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetBottomMargin(0.3);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.05);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  canvas->cd();
  histo_interpolation_DMF->Scale(1,"width");
  histo_interpolation_DMSimp->Scale(1,"width");

  histo_interpolation_DMF->GetYaxis()->SetTitle("Events / GeV ");
  histo_interpolation_DMF->GetYaxis()->SetTitleOffset(1.25);
  histo_interpolation_DMF->SetMarkerStyle(20);
  histo_interpolation_DMF->SetMarkerColor(kBlack);
  histo_interpolation_DMF->SetLineColor(kBlack);
  histo_interpolation_DMF->GetXaxis()->SetLabelSize(0);
  histo_interpolation_DMF->GetXaxis()->SetTitleSize(0);
  histo_interpolation_DMF->Draw("EP");

  histo_interpolation_DMSimp->SetLineColor(kRed);
  histo_interpolation_DMSimp->SetLineWidth(2);

  TH1* histo_interpolation_DMSimp_band = (TH1*) histo_interpolation_DMSimp->Clone("histo_interpolation_DMSimp_band");
  histo_interpolation_DMSimp_band->SetFillColor(kGray);
  histo_interpolation_DMSimp_band->Draw("E2 same");
  histo_interpolation_DMSimp->Draw("hist same");
  histo_interpolation_DMF->Draw("EPsame");

  histo_interpolation_DMF->GetYaxis()->SetRangeUser(min(histo_interpolation_DMF->GetMinimum(),histo_interpolation_DMSimp->GetMinimum())*0.01,
						    max(histo_interpolation_DMF->GetMaximum(),histo_interpolation_DMSimp->GetMaximum())*100);

  if(min(histo_interpolation_DMF->GetMinimum(),histo_interpolation_DMSimp->GetMinimum()) == 0)
    histo_interpolation_DMF->GetYaxis()->SetRangeUser(0.0001,max(histo_interpolation_DMF->GetMaximum(),histo_interpolation_DMSimp->GetMaximum())*100);
  
  TLegend leg (0.65,0.65,0.92,0.92);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.AddEntry((TObject*)(0),("Sample ID "+interpolationPoint).c_str(),"");
  leg.AddEntry(histo_interpolation_DMF,"DMF Interpolated","EP");
  leg.AddEntry(histo_interpolation_DMSimp_band,"DMSimp template","FL");
  leg.Draw("same");

  CMS_lumi(canvas,Form("%.1f",luminosity));
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* ratio = (TH1*) histo_interpolation_DMSimp->Clone("ratio");
  ratio->Divide(histo_interpolation_DMF);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);

  ratio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ratio->GetYaxis()->SetTitle("DMSimp/DMF");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->GetYaxis()->SetRangeUser(0,2);
  ratio->GetYaxis()->SetTitleOffset(1.25);
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
