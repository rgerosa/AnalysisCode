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
      if(postfix == "observed"){
	obs->SetLineColor(kRed+1);
	obs->SetLineWidth(3);
      }
      else if(postfix == "expected"){
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

void makeScanPlots (string inputDIR, string category, string mediatorType, string outputPlotDIR){

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


  TChain* limit = new TChain("limit","limit");
  limit->Add((inputDIR+"/*"+category+"*.root").c_str());

  TGraph2D* observedLimit    = new TGraph2D();
  TGraph2D* observedLimitAlt = new TGraph2D();
  observedLimit->SetNpx(500);
  observedLimit->SetNpy(500);
  observedLimitAlt->SetNpx(500);
  observedLimitAlt->SetNpy(500);
  TGraph2D* expectedLimit = new TGraph2D();
  expectedLimit->SetNpx(500);
  expectedLimit->SetNpy(500);
  TGraph2D* expectedLimit1sUp = new TGraph2D();
  expectedLimit1sUp->SetNpx(500);
  expectedLimit1sUp->SetNpy(500);
  TGraph2D* expectedLimit1sDw = new TGraph2D();
  expectedLimit1sDw->SetNpx(500);
  expectedLimit1sDw->SetNpy(500);
  TGraph2D* expectedLimit2sUp = new TGraph2D();
  expectedLimit2sUp->SetNpx(500);
  expectedLimit2sUp->SetNpy(500);
  TGraph2D* expectedLimit2sDw = new TGraph2D();
  expectedLimit2sDw->SetNpx(500);
  expectedLimit2sDw->SetNpy(500);
  
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

  double minObs = 1e06;
  double minExp = 1e06;
  double minExp1sUp = 1e06;
  double minExp1sDw = 1e06;
  double minExp2sUp = 1e06;
  double minExp2sDw = 1e06;
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

    if(**limitObs < minObs and (**limitObs) != bugLimitObs) minObs = **limitObs;    
    if(**limitObs < minExp and **limitExp != bugLimitExp) minExp = **limitExp;
    if(**limitObs < minExp1sUp and **limitExp1sUp != bugLimitExp) minExp1sUp = **limitExp1sUp;
    if(**limitObs < minExp1sDw and **limitExp1sDw != bugLimitExp) minExp1sDw = **limitExp1sDw;
    if(**limitObs < minExp2sUp and **limitExp2sUp != bugLimitExp) minExp2sUp = **limitExp2sUp;
    if(**limitObs < minExp2sDw and **limitExp2sDw != bugLimitExp) minExp2sDw = **limitExp2sDw;
    
    if(atof(mediatorMass.c_str()) > 2*atof(dmMass.c_str()) and (mediatorType == "scalar" or mediatorType  == "pseudoscalar"))
      observedLimitAlt->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
    else 
      observedLimitAlt->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);

    if(not TString(mediatorType.c_str()).Contains("scalar")){
      observedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
      expectedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
      expectedLimit1sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp1sUp);
      expectedLimit1sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp1sDw);
      expectedLimit2sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp2sUp);
      expectedLimit2sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp2sDw);
    }
    else{
      if(atof(mediatorMass.c_str()) > 2*atof(dmMass.c_str())){
	observedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);
	expectedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp);
	expectedLimit1sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp1sUp);
	expectedLimit1sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp1sDw);
	expectedLimit2sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp2sUp);
	expectedLimit2sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()), **limitExp2sDw);	
      }
      else{

	if(**limitObs <= 1) 
	  observedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  observedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitObs);

	if(**limitExp <= 1) 
	  expectedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  expectedLimit->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitExp);

	if(**limitExp1sUp <= 1) 
	  expectedLimit1sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  expectedLimit1sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitExp1sUp);

	if(**limitExp1sDw <= 1) 
	  expectedLimit1sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  expectedLimit1sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitExp1sDw);

	if(**limitExp2sUp <= 1) 
	  expectedLimit2sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  expectedLimit2sUp->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitExp2sUp);

	if(**limitExp2sDw <= 1) 
	  expectedLimit2sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),1.1);
	else
	  expectedLimit2sDw->SetPoint(nPoint,atof(mediatorMass.c_str()),atof(dmMass.c_str()),**limitExp2sDw);
      }
    }
    nPoint++;				  
  }

  // avoid problem with low mass points
  double* z = observedLimit->GetZ();
  double* x = observedLimit->GetX();
  double* y = observedLimit->GetY();

  for(int N = 0; N < observedLimit->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimit->SetPoint(N,x[N],y[N],z[N]);
  }

  z = observedLimitAlt->GetZ();
  x = observedLimitAlt->GetX();
  y = observedLimitAlt->GetY();

  for(int N = 0; N < observedLimitAlt->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    observedLimitAlt->SetPoint(N,x[N],y[N],z[N]);
  }

  
  z = expectedLimit->GetZ();
  x = expectedLimit->GetX();
  y = expectedLimit->GetY();

  for(int N = 0; N < expectedLimit->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimit->SetPoint(N,x[N],y[N],z[N]);
  }

  z = expectedLimit1sUp->GetZ();
  x = expectedLimit1sUp->GetX();
  y = expectedLimit1sUp->GetY();

  for(int N = 0; N < expectedLimit1sUp->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimit1sUp->SetPoint(N,x[N],y[N],z[N]);
  }

  z = expectedLimit1sDw->GetZ();
  x = expectedLimit1sDw->GetX();
  y = expectedLimit1sDw->GetY();

  for(int N = 0; N < expectedLimit1sDw->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimit1sDw->SetPoint(N,x[N],y[N],z[N]);
  }


  z = expectedLimit2sUp->GetZ();
  x = expectedLimit2sUp->GetX();
  y = expectedLimit2sUp->GetY();

  for(int N = 0; N < expectedLimit2sUp->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimit2sUp->SetPoint(N,x[N],y[N],z[N]);
  }

  z = expectedLimit2sDw->GetZ();
  x = expectedLimit2sDw->GetX();
  y = expectedLimit2sDw->GetY();

  for(int N = 0; N < expectedLimit2sDw->GetN(); N++){
    if(z[N] == bugLimitObs or z[N] == bugLimitExp) z[N] = minObs;
    expectedLimit2sDw->SetPoint(N,x[N],y[N],z[N]);
  }
    

  TCanvas* canvas = new TCanvas("canvas","",750,650);
  canvas->cd();
  canvas->SetLogz();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.18);
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
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);  

  
  //make interpolation of graphs
  TRandom3 rand;
  if(limit != NULL){
    for (Int_t N=0; N<nPointInterpolation; N++) {
      double x,y;
      x = rand.Uniform(medMin,medMax);
      y = rand.Uniform(dmMin,dmMax);
      observedLimit->Interpolate(x,y);
      observedLimitAlt->Interpolate(x,y);
      expectedLimit->Interpolate(x,y);
      expectedLimit1sUp->Interpolate(x,y);
      expectedLimit1sDw->Interpolate(x,y);
      expectedLimit2sUp->Interpolate(x,y);
      expectedLimit2sDw->Interpolate(x,y);
    }
  }

  TH2* observed = observedLimit->GetHistogram();
  TH2* observedAlt = observedLimitAlt->GetHistogram();
  TH2* expected = expectedLimit->GetHistogram();
  TH2* expected1sUp = expectedLimit1sUp->GetHistogram();
  TH2* expected1sDw = expectedLimit1sDw->GetHistogram();
  TH2* expected2sUp = expectedLimit2sUp->GetHistogram();
  TH2* expected2sDw = expectedLimit2sDw->GetHistogram();

  observed->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  observedAlt->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  observed->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/observedLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/observedLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);
  expected->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expected->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/expectedLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/expectedLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);
  expected1sUp->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expected1sUp->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/expected1sUpLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/expected1sUpLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);
  expected1sDw->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expected1sDw->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/expected1sDwLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/expected1sDwLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);
  expected2sUp->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expected2sUp->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/expected2sUpLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/expected2sUpLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
  frame->Draw();
  CMS_lumi(canvas,"12.9",true);
  expected2sDw->GetZaxis()->SetRangeUser(minYaxis,maxYaxis);
  expected2sDw->Draw("colz same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputPlotDIR+"/expected2sDwLimit_"+category+"_"+mediatorType+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/expected2sDwLimit_"+category+"_"+mediatorType+".pdf").c_str(),"pdf");
		   
  // contour lines  
  double contLevel[1];
  contLevel[0] = 1;
  
  TH2* observedCont = (TH2*) observed->Clone("observedCont");
  TH2* expectedCont = (TH2*) expected->Clone("expectedCont");
  TH2* expectedCont1sUp = (TH2*) expected1sUp->Clone("expectedCont1sUp");
  TH2* expectedCont1sDw = (TH2*) expected1sDw->Clone("expectedCont1sDw");
  TH2* expectedCont2sUp = (TH2*) expected2sUp->Clone("expectedCont2sUp");
  TH2* expectedCont2sDw = (TH2*) expected2sDw->Clone("expectedCont2sDw");

  TList* obsContour = getContourList(observedCont,"observed",contLevel);
  TList* expContour = getContourList(expectedCont,"expected",contLevel);
  TList* expContour1sUp = getContourList(expectedCont1sUp,"expected1sUp",contLevel);
  TList* expContour1sDw = getContourList(expectedCont1sDw,"expected1sDw",contLevel);
  TList* expContour2sUp = getContourList(expectedCont2sUp,"expected2sUp",contLevel);
  TList* expContour2sDw = getContourList(expectedCont2sDw,"expected2sDw",contLevel);

  TCanvas* c = new TCanvas("c","",750,650);
  c->cd();
  c->SetLogz();
  c->SetTickx();
  c->SetTicky();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.15);
 
  c->cd();
  frame->Draw();
  observedAlt->Draw("COLZ same");

  obsContour->Draw("Lsame");  
  expContour->Draw("Lsame");  
  if(mediatorType == "vector" or mediatorType == "axial"){
    expContour1sUp->Draw("Lsame");  
    expContour1sDw->Draw("Lsame");  
  }
  else{
    expContour1sDw->Draw("Lsame");
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

  leg.AddEntry(obsContour->At(0),"Observed 95% CL","L"); 
  leg.AddEntry(expContour->At(0),"Median Expected 95% CL","L"); 
  if(mediatorType == "vector" or mediatorType == "axial")
    leg.AddEntry(expContour1sUp->At(0),"Expected #pm 1#sigma","L"); 
  else
    leg.AddEntry(expContour1sUp->At(0),"Expected - 1#sigma","L"); 

  leg.Draw("same");

  c->SaveAs((outputPlotDIR+"/scan_"+category+"_"+mediatorType+".png").c_str());
  c->SaveAs((outputPlotDIR+"/scan_"+category+"_"+mediatorType+".pdf").c_str());

}  
