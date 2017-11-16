#include "../CMS_lumi.h"

void plotHistogram(TCanvas* canvas, const vector<TH1*> & histogramToPlot, const string & outputDIR, const bool & scaleToOne, const string & postfix){

  if(scaleToOne){
    for(auto hist : histogramToPlot)
      hist->Scale(1./hist->Integral());
  }
  else{
    for(auto hist : histogramToPlot)
      hist->Scale(1,"width");
  }    
    

  canvas->cd();
  int icolor = 1;
  double minimum = 99999;
  double maximum = -9999;

  TLegend leg (0.45,0.53,0.97,0.91);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  

  for(auto hist : histogramToPlot){

    hist->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    hist->GetYaxis()->SetTitle("Events/GeV");
    minimum = min(minimum,hist->GetMinimum());
    maximum = max(maximum,hist->GetMaximum());

    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    
    hist->SetLineColor(icolor);
    hist->SetLineWidth(2);
    if(icolor == 1) hist->Draw("hist");
    else hist->Draw("hist same");

    leg.AddEntry(hist,hist->GetName(),"L");

    icolor++;
  }
  
  leg.Draw("same");
  CMS_lumi(canvas,"35.9");

  if(not scaleToOne){
    if(minimum > 0)
      histogramToPlot.at(0)->GetYaxis()->SetRangeUser(minimum*0.1,maximum*100);
    else
      histogramToPlot.at(0)->GetYaxis()->SetRangeUser(0.001,maximum*100);
  }
  else{
    if(minimum > 0)
      histogramToPlot.at(0)->GetYaxis()->SetRangeUser(minimum*0.01,maximum*10);
    else
      histogramToPlot.at(0)->GetYaxis()->SetRangeUser(0.001,maximum*10);
  }

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");

}


void makeTemplatesLLPlots(string inputFileName, string outputDIR, string toGrep){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();

  vector<TH1*> histogramToPlot;
  TIter next(inputFile->GetListOfKeys());
  TKey *key;

  std::vector<std::string> Mmed;
  std::vector<std::string> Mchi2;
  std::vector<std::string> ctau;

  while ((key = (TKey*)next())) {
    // break the name
    TString histoName(key->GetName());
    if(not histoName.Contains(toGrep.c_str())) continue;
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    histogramToPlot.push_back((TH1*) key->ReadObj());
  }
  
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();  
  plotHistogram(canvas,histogramToPlot,outputDIR,false,"templates_"+toGrep);
  plotHistogram(canvas,histogramToPlot,outputDIR,true,"templates_"+toGrep+"_shape");
  
}
