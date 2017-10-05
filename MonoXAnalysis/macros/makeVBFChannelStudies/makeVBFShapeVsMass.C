#include "../CMS_lumi.h"

static vector<int> vcolor = {602,600,855,858,867,866,852};

void plotDistribution(TCanvas* canvas, vector<pair<string,TH1F*> > & histo, const string & outputDIR, const string & obs){

  canvas->cd();
  int icolor = 0;

  double min = 1;
  double max = 0;

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  double integral = 0;

  //// -----
  for(auto hist : histo){

    //// -----
    if(icolor == 0){
      hist.second->Scale(1.,"width");
      integral = hist.second->Integral();
    }
    else{
      hist.second->Scale(1.,"width");
      hist.second->Scale(integral/hist.second->Integral());
    }      

    hist.second->SetLineColor(vcolor.at(icolor));
    hist.second->SetMarkerColor(vcolor.at(icolor));
    hist.second->SetLineWidth(2);
    hist.second->SetMarkerSize(1);
    hist.second->SetMarkerStyle(20);

    //// -----
    if(obs == "mjj")
      hist.second->GetXaxis()->SetTitle("M_{jj} [GeV]");
    else if(obs == "met")
      hist.second->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    else if(obs == "dphijj")
      hist.second->GetXaxis()->SetTitle("#Delta#phi_{jj}");
    else if(obs == "detajj")
      hist.second->GetXaxis()->SetTitle("#Delta#eta_{jj}");

    hist.second->GetYaxis()->SetTitle("a.u.");
    hist.second->GetYaxis()->SetTitleOffset(1.25);
    hist.second->GetXaxis()->SetTitleOffset(1.25);

    //// -----
    icolor++;
    if(icolor == 1){
      hist.second->Draw("hist");
      hist.second->Draw("P same");
      leg.AddEntry(hist.second,Form("VBF m_{H} = %s GeV",hist.first.c_str()),"LP");
    }
    else{
      hist.second->Draw("hist same");      
      hist.second->Draw("P same");      
      leg.AddEntry(hist.second,Form("VBF m_{H} = %s GeV",hist.first.c_str()),"LP");
    }

  }

  histo.at(0).second->GetYaxis()->SetRangeUser(0,histo.at(0).second->GetMaximum()*1.5);

  CMS_lumi(canvas,"");
  leg.Draw();

  //canvas->SetLogy();
  
  canvas->SaveAs((outputDIR+"/shapeComparison_vsMass_"+obs+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/shapeComparison_vsMass_"+obs+".pdf").c_str(),"pdf");

}

void makeVBFShapeVsMass(string inputTemplateFileName, string outputDIR, vector<string> observables){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  
  TFile* inputTemplateFile = TFile::Open(inputTemplateFileName.c_str(),"READ");

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  
  for(auto obs: observables){

    vector<pair<string,TH1F*> > histoToPlot;
    histoToPlot.push_back(make_pair<string,TH1F*>("110",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_110_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("125",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_125_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("200",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_200_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("400",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_400_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("600",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_600_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("800",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_800_"+obs).c_str())));
    histoToPlot.push_back(make_pair<string,TH1F*>("1000",(TH1F*) inputTemplateFile->Get(("vbfH/vbfHhist_1000_"+obs).c_str())));   
    plotDistribution(canvas,histoToPlot,outputDIR,obs);
    
  }

}
