#include "../CMS_lumi.h"

static float medMin = 0;
static float medMax = 500;
static float yMin   = 0;
static float yMax   = 10;


void makeBrazilianPlots (string inputDIR, string category, string mediatorType, string outputPlotDIR, double DMmass = 0){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetPalette(57);

  system(("mkdir -p "+outputPlotDIR).c_str());

  TChain* limit = new TChain("limit","limit");
  limit->Add((inputDIR+"/*"+category+"*.root").c_str());

  TGraphAsymmErrors* observedLimit = new TGraphAsymmErrors();
  TGraphAsymmErrors* expectedLimit = new TGraphAsymmErrors();
  TGraphAsymmErrors* expectedLimit1s = new TGraphAsymmErrors();
  TGraphAsymmErrors* expectedLimit2s = new TGraphAsymmErrors();
  
  TTreeReader* reader = new TTreeReader(limit);

  TTreeReaderValue<float>* limitObs = new TTreeReaderValue<float>(*reader,"limitObs");
  TTreeReaderValue<float>* limitExp = new TTreeReaderValue<float>(*reader,"limitExp");
  TTreeReaderValue<float>* limitExp1sUp = new TTreeReaderValue<float>(*reader,"limitExpUp1s");
  TTreeReaderValue<float>* limitExp2sUp = new TTreeReaderValue<float>(*reader,"limitExpUp2s");
  TTreeReaderValue<float>* limitExp1sDw = new TTreeReaderValue<float>(*reader,"limitExpDw1s");
  TTreeReaderValue<float>* limitExp2sDw = new TTreeReaderValue<float>(*reader,"limitExpDw2s");
  TTreeReaderValue<double>* mh = new TTreeReaderValue<double>(*reader,"mh");

  long int nPoint = 0;

  string mass ;
  string mediatorMass;
  string dmMass;
  double MH;

  // set branches
  while(reader->Next()){
    
    MH = **mh;
    mass = to_string(MH);
    
    mediatorMass = mass[3];
    mediatorMass += mass[4];
    mediatorMass += mass[5];
    mediatorMass += mass[6];
    dmMass       = mass[7];
    dmMass      += mass[8];
    dmMass      += mass[9];
    dmMass      += mass[10];

    if(atof(dmMass.c_str()) != DMmass) continue;

    observedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),**limitObs);
    observedLimit->SetPointError(nPoint,0,0,0,0);

    expectedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()), **limitExp);
    expectedLimit->SetPointError(nPoint,0,0,0,0);

    expectedLimit1s->SetPoint(nPoint,atof(mediatorMass.c_str()),**limitExp);
    expectedLimit1s->SetPointError(nPoint,50,50,fabs(**limitExp1sDw-**limitExp),fabs(**limitExp-**limitExp1sUp));

    expectedLimit2s->SetPoint(nPoint,atof(mediatorMass.c_str()), **limitExp);
    expectedLimit2s->SetPointError(nPoint,50,50,fabs(**limitExp2sDw-**limitExp),fabs(**limitExp-**limitExp2sUp));

    nPoint++;				  
  }


  TCanvas* canvas = new TCanvas("canvas","",650,650);
  canvas->cd();
  canvas->SetTickx();
  canvas->SetTicky();
  TH1F* frame = canvas->DrawFrame(medMin,yMin,medMax,yMax,"");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetYaxis()->SetTitle("95% upper limit on #sigma/#sigma_{th}");
  frame->GetYaxis()->CenterTitle();
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetLabelSize(0.035,"Z");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);  

  expectedLimit2s->SetFillColor(kYellow);
  expectedLimit2s->SetLineColor(kBlack);
  expectedLimit2s->Draw("3same");
  expectedLimit1s->SetFillColor(kGreen);
  expectedLimit1s->SetLineColor(kBlack);
  expectedLimit1s->Draw("3same");

  expectedLimit->SetLineColor(kBlack);
  expectedLimit->SetLineWidth(2);
  expectedLimit->SetLineStyle(7);
  expectedLimit->Draw("Csame");

  observedLimit->SetLineColor(kBlack);
  observedLimit->SetLineWidth(2);
  observedLimit->Draw("Csame");
 
  TLegend leg (0.2,0.6,0.65,0.9);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  if(mediatorType == "scalar")
    leg.AddEntry((TObject*)0,"Scalar med., Dirac DM, g_{q}=1, g_{DM}=1","");
  else if(mediatorType == "pseudoscalar")
    leg.AddEntry((TObject*)0,"Pseudoscalar med., Dirac DM, g_{q}=1, g_{DM}=1","");

  TF1* f1 = new TF1("f1","1",medMin,medMax);
  f1->SetLineColor(kRed);
  f1->SetLineWidth(2);
  f1->Draw("Lsame");
  expectedLimit->Draw("Csame");
  observedLimit->Draw("Csame");

  leg.AddEntry(observedLimit,"Observed 95% CL","L"); 
  leg.AddEntry(expectedLimit,"Median Expected 95% CL","L"); 
  leg.AddEntry(expectedLimit1s,"Expected #pm 1#sigma","F"); 
  leg.AddEntry(expectedLimit2s,"Expected #pm 2#sigma","F"); 
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputPlotDIR+"/brazilian_"+category+"_"+mediatorType+".png").c_str());
  canvas->SaveAs((outputPlotDIR+"/brazilian_"+category+"_"+mediatorType+".pdf").c_str());

}  
