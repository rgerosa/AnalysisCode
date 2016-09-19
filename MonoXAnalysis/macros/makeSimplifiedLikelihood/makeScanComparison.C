#include "../CMS_lumi.h"

static float medMin = 50;
static float medMax = 2500;
static float dmMin  = 5;
static float dmMax  = 1200;
static int   nPointInterpolation = 500000;

static float bugLimitObs = 0.001;
static float bugLimitExp = 280.49469;
static float minYaxis = 0.01;
static float maxYaxis = 10;

TList* getContourList(TH2* inputHisto, const string & postfix, double* contLevel){
  
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",750,650);
  canvas->cd();
  TH2* cont = (TH2*) inputHisto->Clone(("cont_"+postfix).c_str());
  cont->SetContour(1,contLevel);
  cont->Draw("CONT,list");
  canvas->Update();  
  TList* obsContour = new TList();
  TObjArray *contList = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  for(int i = 0; i < contList->GetSize(); i++){
    TList*  contLevel = (TList*) contList->At(i);
    for(int j = 0; j < contLevel->GetSize(); j++){
      TGraph* obs = new TGraph(*((TGraph*)(contLevel->At(j))));      
      if(TString(postfix).Contains("observed")){
	obs->SetLineColor(kRed+1);
	obs->SetLineWidth(3);
      }
      else if(TString(postfix).Contains("expected")){
	obs->SetLineColor(kBlack);
	obs->SetLineWidth(3);
      }
      else{
	obs->SetLineColor(kGray+2);
	obs->SetLineWidth(2);
      }      
      obsContour->Add((TGraph*) obs->Clone());
    }
  }
  
  return obsContour;
  
}

// compare observed contours
void makeScanComparison (string category, string mediatorType, string outputPlotDIR){

  string inputDIRSL = "ResultsSimplifiedLikelihood/result_vector_g25/";
  string inputDIRSLNoCorr = "ResultsSimplifiedLikelihood/result_vector_g25_nocorr/";
  string inputDIRSLSR = "ResultsSimplifiedLikelihood/result_vector_g25_SR/";
  string inputDIRSLSRNoCorr = "ResultsSimplifiedLikelihood/result_vector_g25_SR_nocorr/";
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetPalette(57);

  system(("mkdir -p "+outputPlotDIR).c_str());

  if(mediatorType == "scalar" or mediatorType == "pseudoscalar"){
    medMin = 0;
    medMax = 500;
    dmMin  = 5;
    dmMax  = 250;
  }

  if(mediatorType == "scalar" or mediatorType == "pseudoscalar"){
    minYaxis = minYaxis*20;
    maxYaxis = maxYaxis*0.7;
  }


  TChain* limitSL = new TChain("limit","limit");
  TChain* limitSLNoCorr = new TChain("limit","limit");
  TChain* limitSLSR = new TChain("limit","limit");
  TChain* limitSLSRNoCorr = new TChain("limit","limit");

  limitSL->Add((inputDIRSL+"/*"+category+"*.root").c_str());
  limitSLNoCorr->Add((inputDIRSLNoCorr+"/*"+category+"*.root").c_str());
  limitSLSR->Add((inputDIRSLSR+"/*"+category+"*.root").c_str());
  limitSLSRNoCorr->Add((inputDIRSLSRNoCorr+"/*"+category+"*.root").c_str());

  TGraph2D* observedLimitSL = new TGraph2D();
  observedLimitSL->SetNpx(500);
  observedLimitSL->SetNpy(500);
  TGraph2D* observedLimitSLNoCorr = new TGraph2D();
  observedLimitSLNoCorr->SetNpx(500);
  observedLimitSLNoCorr->SetNpy(500);
  TGraph2D* observedLimitSLSR = new TGraph2D();
  observedLimitSLSR->SetNpx(500);
  observedLimitSLSR->SetNpy(500);
  TGraph2D* observedLimitSLSRNoCorr = new TGraph2D();
  observedLimitSLSRNoCorr->SetNpx(500);
  observedLimitSLSRNoCorr->SetNpy(500);

  TGraph2D* expectedLimitSL = new TGraph2D();
  expectedLimitSL->SetNpx(500);
  expectedLimitSL->SetNpy(500);
  TGraph2D* expectedLimitSLNoCorr = new TGraph2D();
  expectedLimitSLNoCorr->SetNpx(500);
  expectedLimitSLNoCorr->SetNpy(500);
  TGraph2D* expectedLimitSLSR = new TGraph2D();
  expectedLimitSLSR->SetNpx(500);
  expectedLimitSLSR->SetNpy(500);
  TGraph2D* expectedLimitSLSRNoCorr = new TGraph2D();
  expectedLimitSLSRNoCorr->SetNpx(500);
  expectedLimitSLSRNoCorr->SetNpy(500);

  string mass ;
  string mediatorMass;
  string dmMass;
  double MH;
  double minObs = 1e06;
  double minExp = 1e06;
  
  TTreeReader* readerSL = new TTreeReader(limitSL);
  TTreeReaderValue<float>* limitObs = new TTreeReaderValue<float>(*readerSL,"limitObs");
  TTreeReaderValue<float>* limitExp = new TTreeReaderValue<float>(*readerSL,"limitExp");
  TTreeReaderValue<double>* mh = new TTreeReaderValue<double>(*readerSL,"mh");
  long int nPoint = 0;
  // set branches
  while(readerSL->Next()){    
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
    if(**limitObs < minObs and (**limitObs) != bugLimitObs) minObs = **limitObs;    
    if(**limitExp < minExp and **limitExp   != bugLimitExp) minExp = **limitExp;    
    observedLimitSL->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
    expectedLimitSL->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
    nPoint++;				  
  }

  ////////
  minObs = 1e06;
  minExp = 1e06;
  TTreeReader* readerSLNoCorr = new TTreeReader(limitSLNoCorr);
  limitObs = new TTreeReaderValue<float>(*readerSLNoCorr,"limitObs");
  limitExp = new TTreeReaderValue<float>(*readerSLNoCorr,"limitExp");
  mh = new TTreeReaderValue<double>(*readerSLNoCorr,"mh");
  nPoint = 0;
  while(readerSLNoCorr->Next()){
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
    if(**limitObs < minObs and (**limitObs) != bugLimitObs) minObs = **limitObs;
    if(**limitExp < minExp and **limitExp   != bugLimitExp) minExp = **limitExp;
    observedLimitSLNoCorr->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
    expectedLimitSLNoCorr->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
    nPoint++;
  }

  ////////
  minObs = 1e06;
  minExp = 1e06;
  TTreeReader* readerSLSR = new TTreeReader(limitSLSR);
  limitObs = new TTreeReaderValue<float>(*readerSLSR,"limitObs");
  limitExp = new TTreeReaderValue<float>(*readerSLSR,"limitExp");
  mh = new TTreeReaderValue<double>(*readerSLSR,"mh");
  nPoint = 0;
  while(readerSLSR->Next()){
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
    if(**limitObs < minObs and (**limitObs) != bugLimitObs) minObs = **limitObs;
    if(**limitExp < minExp and **limitExp   != bugLimitExp) minExp = **limitExp;
    observedLimitSLSR->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
    expectedLimitSLSR->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
    nPoint++;
  }

  ////////
  minObs = 1e06;
  minExp = 1e06;
  TTreeReader* readerSLSRNoCorr = new TTreeReader(limitSLSRNoCorr);
  limitObs = new TTreeReaderValue<float>(*readerSLSRNoCorr,"limitObs");
  limitExp = new TTreeReaderValue<float>(*readerSLSRNoCorr,"limitExp");
  mh = new TTreeReaderValue<double>(*readerSLSRNoCorr,"mh");
  nPoint = 0;
  while(readerSLSRNoCorr->Next()){
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
    if(**limitObs < minObs and (**limitObs) != bugLimitObs) minObs = **limitObs;
    if(**limitExp < minExp and **limitExp   != bugLimitExp) minExp = **limitExp;
    observedLimitSLSRNoCorr->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
    expectedLimitSLSRNoCorr->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
    nPoint++;
  }
    
  // avoid problem with low mass points
  double* z = observedLimitSL->GetZ();
  double* x = observedLimitSL->GetX();
  double* y = observedLimitSL->GetY();

  for(int N = 0; N < observedLimitSL->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimitSL->SetPoint(N,x[N],y[N],z[N]);
  }
  
  z = expectedLimitSL->GetZ();
  x = expectedLimitSL->GetX();
  y = expectedLimitSL->GetY();

  for(int N = 0; N < expectedLimitSL->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimitSL->SetPoint(N,x[N],y[N],z[N]);
  }

  z = observedLimitSLNoCorr->GetZ();
  x = observedLimitSLNoCorr->GetX();
  y = observedLimitSLNoCorr->GetY();

  for(int N = 0; N < observedLimitSLNoCorr->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimitSLNoCorr->SetPoint(N,x[N],y[N],z[N]);
  }
  
  z = expectedLimitSLNoCorr->GetZ();
  x = expectedLimitSLNoCorr->GetX();
  y = expectedLimitSLNoCorr->GetY();

  for(int N = 0; N < expectedLimitSLNoCorr->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimitSLNoCorr->SetPoint(N,x[N],y[N],z[N]);
  }

  z = observedLimitSLSR->GetZ();
  x = observedLimitSLSR->GetX();
  y = observedLimitSLSR->GetY();

  for(int N = 0; N < observedLimitSLSR->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimitSLSR->SetPoint(N,x[N],y[N],z[N]);
  }
  
  z = expectedLimitSLSR->GetZ();
  x = expectedLimitSLSR->GetX();
  y = expectedLimitSLSR->GetY();

  for(int N = 0; N < expectedLimitSLSR->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimitSLSR->SetPoint(N,x[N],y[N],z[N]);
  }

  z = observedLimitSLSRNoCorr->GetZ();
  x = observedLimitSLSRNoCorr->GetX();
  y = observedLimitSLSRNoCorr->GetY();

  for(int N = 0; N < observedLimitSLSRNoCorr->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimitSLSRNoCorr->SetPoint(N,x[N],y[N],z[N]);
  }
  
  z = expectedLimitSLSRNoCorr->GetZ();
  x = expectedLimitSLSRNoCorr->GetX();
  y = expectedLimitSLSRNoCorr->GetY();

  for(int N = 0; N < expectedLimitSLSRNoCorr->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimitSLSRNoCorr->SetPoint(N,x[N],y[N],z[N]);
  }

  ///////////
  TCanvas* canvas = new TCanvas("canvas","",750,650);
  canvas->cd();
  canvas->SetLogz();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.10);
  TH1F* frame = canvas->DrawFrame(medMin,dmMin,medMax,dmMax,"");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetYaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetZaxis()->SetTitleOffset(0.1);
  frame->GetYaxis()->SetNdivisions(510);
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetLabelSize(0.035,"Z");

  
  //make interpolation of graphs
  TRandom3 rand;
  for (Int_t N=0; N<nPointInterpolation; N++) {
    double x,y;
    x = rand.Uniform(medMin,medMax);
    y = rand.Uniform(dmMin,dmMax);
    observedLimitSL->Interpolate(x,y);
    expectedLimitSL->Interpolate(x,y);
    observedLimitSLNoCorr->Interpolate(x,y);
    expectedLimitSLNoCorr->Interpolate(x,y);
    observedLimitSLSR->Interpolate(x,y);
    expectedLimitSLSR->Interpolate(x,y);
    observedLimitSLSRNoCorr->Interpolate(x,y);
    expectedLimitSLSRNoCorr->Interpolate(x,y);
  }

  TH2* observedSL = observedLimitSL->GetHistogram();
  TH2* expectedSL = expectedLimitSL->GetHistogram();
  TH2* observedSLNoCorr = observedLimitSLNoCorr->GetHistogram();
  TH2* expectedSLNoCorr = expectedLimitSLNoCorr->GetHistogram();
  TH2* observedSLSR = observedLimitSLSR->GetHistogram();
  TH2* expectedSLSR = expectedLimitSLSR->GetHistogram();
  TH2* observedSLSRNoCorr = observedLimitSLSRNoCorr->GetHistogram();
  TH2* expectedSLSRNoCorr = expectedLimitSLSRNoCorr->GetHistogram();

  observedSL->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expectedSL->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  observedSLNoCorr->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expectedSLNoCorr->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  observedSLSR->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expectedSLSR->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  observedSLSRNoCorr->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expectedSLSRNoCorr->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
		   
  // contour lines  
  double contLevel[1];
  contLevel[0] = 1;
  
  TH2* observedContSL = (TH2*) observedSL->Clone("observedContSL");
  TH2* expectedContSL = (TH2*) expectedSL->Clone("expectedContSL");
  TH2* observedContSLNoCorr = (TH2*) observedSLNoCorr->Clone("observedContSLNoCorr");
  TH2* expectedContSLNoCorr = (TH2*) expectedSLNoCorr->Clone("expectedContSLNoCorr");
  TH2* observedContSLSR = (TH2*) observedSLSR->Clone("observedContSLSR");
  TH2* expectedContSLSR = (TH2*) expectedSLSR->Clone("expectedContSLSR");
  TH2* observedContSLSRNoCorr = (TH2*) observedSLSR->Clone("observedContSLSRNoCorr");
  TH2* expectedContSLSRNoCorr = (TH2*) expectedSLSR->Clone("expectedContSLSRNoCorr");

  TList* obsContourSL = getContourList(observedContSL,"observedSL",contLevel);
  TList* expContourSL = getContourList(expectedContSL,"expectedSL",contLevel);
  TList* obsContourSLNoCorr = getContourList(observedContSLNoCorr,"observedSLNoCorr",contLevel);
  TList* expContourSLNoCorr = getContourList(expectedContSLNoCorr,"expectedSLNoCorr",contLevel);
  TList* obsContourSLSR = getContourList(observedContSLSR,"observedSLSR",contLevel);
  TList* expContourSLSR = getContourList(expectedContSLSR,"expectedSLSR",contLevel);
  TList* obsContourSLSRNoCorr = getContourList(observedContSLSRNoCorr,"observedSLSRNoCorr",contLevel);
  TList* expContourSLSRNoCorr = getContourList(expectedContSLSRNoCorr,"expectedSLSRNoCorr",contLevel);

  //////////////////
  TCanvas* c = new TCanvas("c","",750,650);
  c->cd();
  c->SetLogz();
  c->SetTickx();
  c->SetTicky();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.10);
  c->cd();
  frame->Draw();

  for(int i = 0; i < obsContourSL->GetSize(); i++){
    TGraph* gr = (TGraph*) obsContourSL->At(i);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < obsContourSLNoCorr->GetSize(); i++){
    TGraph* gr = (TGraph*) obsContourSLNoCorr->At(i);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < obsContourSLSR->GetSize(); i++){
    TGraph* gr = (TGraph*) obsContourSLSR->At(i);
    gr->SetLineColor(kGreen+1);
    gr->SetLineWidth(4);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < obsContourSLSRNoCorr->GetSize(); i++){
    TGraph* gr = (TGraph*) obsContourSLSRNoCorr->At(i);
    gr->SetLineColor(kOrange);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  c->RedrawAxis("sameaxis");

  TLatex * texLumi = new TLatex();
  texLumi->SetNDC();
  texLumi->SetTextFont(42);
  texLumi->SetLineWidth(2);
  texLumi->SetTextSize(0.036);
  texLumi->Draw();
  texLumi->DrawLatex(0.60,0.94,"12.9 fb^{-1} (13 TeV)");

  TLatex* texCMS = new TLatex(0.12,0.94,"#bf{CMS} #it{Simplfied Likelihood}");
  texCMS->SetNDC();
  texCMS->SetTextFont(42);
  texCMS->SetLineWidth(2);
  texCMS->SetTextSize(0.040); texCMS->Draw();

  TLatex *   texZaxis = new TLatex();
  texZaxis->SetNDC();
  texZaxis->SetTextFont(42);
  texZaxis->SetLineWidth(2);
  texZaxis->SetTextSize(0.042);
  texZaxis->SetTextAngle(270);
  texZaxis->DrawLatex(0.963,0.73,"Observed #sigma_{95% CL}/#sigma_{th}");

  TLegend leg (0.2,0.6,0.65,0.9);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  if(mediatorType == "vector")
    leg.AddEntry((TObject*)0,"Vector med., Dirac DM, g_{q}=0.25, g_{DM}=1","");
  else if(mediatorType == "axial")
    leg.AddEntry((TObject*)0,"Axial Vector med., Dirac DM, g_{q}=0.25, g_{DM}=1","");
  else if(mediatorType == "scalar")
    leg.AddEntry((TObject*)0,"Scalar med., Dirac DM, g_{q}=1, g_{DM}=1","");
  else if(mediatorType == "pseudoscalar")
    leg.AddEntry((TObject*)0,"PseudoScalar med., Dirac DM, g_{q}=1, g_{DM}=1","");

  leg.AddEntry(obsContourSL->At(0),"Observed 95% CL Simp. Likelihood","L"); 
  leg.AddEntry(obsContourSLNoCorr->At(0),"Observed 95% CL Simp. Likelihood w/o corr.","L"); 
  leg.AddEntry(obsContourSLSR->At(0),"Observed 95% CL Fit SR","L"); 
  leg.AddEntry(obsContourSLSRNoCorr->At(0),"Observed 95% CL Fit SR w/o corr.","L"); 
  leg.Draw("same");

  c->SaveAs((outputPlotDIR+"/scan_comparison_observed_"+category+"_"+mediatorType+".png").c_str());
  c->SaveAs((outputPlotDIR+"/scan_comparison_observed_"+category+"_"+mediatorType+".pdf").c_str());

  ////////////
  frame->Draw();

  for(int i = 0; i < expContourSL->GetSize(); i++){
    TGraph* gr = (TGraph*) expContourSL->At(i);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < expContourSLNoCorr->GetSize(); i++){
    TGraph* gr = (TGraph*) expContourSLNoCorr->At(i);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < expContourSLSR->GetSize(); i++){
    TGraph* gr = (TGraph*) expContourSLSR->At(i);
    gr->SetLineColor(kGreen+1);
    gr->SetLineWidth(4);
    gr->Draw("Lsame");  
  }

  for(int i = 0; i < expContourSLSRNoCorr->GetSize(); i++){
    TGraph* gr = (TGraph*) expContourSLSRNoCorr->At(i);
    gr->SetLineColor(kOrange);
    gr->SetLineWidth(2);
    gr->Draw("Lsame");  
  }

  c->RedrawAxis("sameaxis");
  texLumi->DrawLatex(0.60,0.94,"12.9 fb^{-1} (13 TeV)");
  texCMS->Draw();
  texZaxis->DrawLatex(0.963,0.73,"Expected #sigma_{95% CL}/#sigma_{th}");
 
  leg.Clear();
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  if(mediatorType == "vector")
    leg.AddEntry((TObject*)0,"Vector med., Dirac DM, g_{q}=0.25, g_{DM}=1","");
  else if(mediatorType == "axial")
    leg.AddEntry((TObject*)0,"Axial Vector med., Dirac DM, g_{q}=0.25, g_{DM}=1","");
  else if(mediatorType == "scalar")
    leg.AddEntry((TObject*)0,"Scalar med., Dirac DM, g_{q}=1, g_{DM}=1","");
  else if(mediatorType == "pseudoscalar")
    leg.AddEntry((TObject*)0,"PseudoScalar med., Dirac DM, g_{q}=1, g_{DM}=1","");

  leg.AddEntry(expContourSL->At(0),"Observed 95% CL Simp. Likelihood","L"); 
  leg.AddEntry(expContourSLNoCorr->At(0),"Observed 95% CL Simp. Likelihood w/o corr.","L"); 
  leg.AddEntry(expContourSLSR->At(0),"Observed 95% CL Fit SR","L"); 
  leg.AddEntry(expContourSLSRNoCorr->At(0),"Observed 95% CL Fit SR w/o corr.","L"); 
  leg.Draw("same");

  c->SaveAs((outputPlotDIR+"/scan_comparison_expected_"+category+"_"+mediatorType+".png").c_str());
  c->SaveAs((outputPlotDIR+"/scan_comparison_expected_"+category+"_"+mediatorType+".pdf").c_str());
}  
