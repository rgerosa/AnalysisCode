#include "../CMS_lumi.h"
#include "../../../../HiggsAnalysis/CombinedLimit/interface/RooSplineND.h"

static int nbinsX = 75;
static int nbinsY = 75;

void makeKappaVKappaFLimit(string inputFileName, string outputDIR, string postfix, bool useObserved = true, bool useSpline = true){

  gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc530/libHiggsAnalysisCombinedLimit.so");

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  gStyle->SetPalette(kLightTemperature);

  // take the LHC best fit result
  TFile *bestfit_lhc = TFile::Open("externalFiles/scan2d_kappa_V_kappa_F_exp.root");
  TGraph *grlhc_68 = (TGraph*) bestfit_lhc->Get("graph68_comb_0");
  TGraph *grlhc_95 = (TGraph*) bestfit_lhc->Get("graph95_comb_0");
  TGraph *grlhc_bf = (TGraph*) bestfit_lhc->Get("comb_best");
  grlhc_68->SetLineColor(kBlack);
  grlhc_68->SetLineWidth(3);
  grlhc_95->SetLineColor(kBlack);
  grlhc_95->SetLineWidth(3);
  grlhc_95->SetLineStyle(7);
  grlhc_bf->SetMarkerStyle(34);
  grlhc_bf->SetMarkerSize(2.0);
  grlhc_bf->SetMarkerColor(kBlack);
  

  TGraph* sm_value = new TGraph();
  sm_value->SetPoint(0,1,1);
  sm_value->SetMarkerStyle(20);
  sm_value->SetMarkerSize(1.8);
  sm_value->SetMarkerColor(kWhite);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
  TTree* tree_limit = (TTree*) inputFile->Get("limit");

  TTreeReader reader (tree_limit);
  TTreeReaderValue<float> kappa_V (reader,"trackedParam_kappa_V");
  TTreeReaderValue<float> kappa_F (reader,"trackedParam_kappa_F");
  TTreeReaderValue<float> quantileExpected (reader,"quantileExpected");
  TTreeReaderValue<double> limit (reader,"limit");

  TGraph2D* graph_limit = new TGraph2D();
  int    npoints = 0;

  float  kappa_V_min = 100;
  float  kappa_V_max = 0;
  float  kappa_F_min = 100;
  float  kappa_F_max = 0;
  float  limit_min = 100;
  float  limit_max = 0;

  while(reader.Next()){    
    if(useObserved and *quantileExpected != -1) continue;
    else if(not useObserved and *quantileExpected != 0.5) continue;
    graph_limit->SetPoint(npoints,*kappa_V,*kappa_F,*limit);
    npoints++;   
    if(*kappa_V > kappa_V_max) kappa_V_max = *kappa_V;
    if(*kappa_V < kappa_V_min) kappa_V_min = *kappa_V;
    if(*kappa_F > kappa_F_max) kappa_F_max = *kappa_F;
    if(*kappa_F < kappa_F_min) kappa_F_min = *kappa_F;
    if(*limit < limit_min) limit_min = *limit;
    if(*limit > limit_max) limit_max = *limit;    
  }

  TH2F* histo_limit = NULL;
  if(useSpline)
    histo_limit = new TH2F("histo_limit","",nbinsX,kappa_V_min,kappa_V_max,nbinsY,kappa_F_min,kappa_F_max);
  else
    histo_limit = new TH2F("histo_limit","",nbinsX*3,kappa_V_min,kappa_V_max,nbinsY*3,kappa_F_min,kappa_F_max);

  for(int iBinX = 0; iBinX < histo_limit->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < histo_limit->GetNbinsY(); iBinY++){
      histo_limit->SetBinContent(iBinX+1,iBinY+1,graph_limit->Interpolate(histo_limit->GetXaxis()->GetBinCenter(iBinX+1),histo_limit->GetYaxis()->GetBinCenter(iBinY+1)));
    }
  }
  histo_limit->Smooth();

  RooSplineND *spline = NULL;
  RooRealVar* xvar = new RooRealVar("x","x",(kappa_V_max-kappa_V_min)/2,kappa_V_min,kappa_V_max);
  RooRealVar* yvar = new RooRealVar("y","y",(kappa_F_max-kappa_F_min)/2,kappa_F_min,kappa_F_max);

  TFile* file_temp = new TFile("file_temp.root","RECREATE");
  TTree *tree = new TTree("tree","tree"); 
  float  x,y,z;
  tree->Branch("x",&x,"x/F");
  tree->Branch("y",&y,"y/F");
  tree->Branch("z",&z,"z/F");
      
  if(useSpline){
    for(int i = 1; i <= histo_limit->GetNbinsX(); i++){
      for(int j = 1; j <= histo_limit->GetNbinsY(); j++){
	x = histo_limit->GetXaxis()->GetBinCenter(i);
	y = histo_limit->GetYaxis()->GetBinCenter(j);
	z = histo_limit->GetBinContent(i,j);
	tree->Fill();
      }
    }
    
    RooArgList list (*xvar,*yvar);
    spline = new RooSplineND(Form("spline_%s",postfix.c_str()),Form("spline_%s",postfix.c_str()),list,tree,"z",1,true);  
  }
  if(tree) delete tree;

  TH2F* histo_limit_ext = NULL;
  if(useSpline)
    histo_limit_ext = new TH2F("histo_limit_ext","",nbinsX*3,kappa_V_min,kappa_V_max,nbinsY*3,kappa_F_min,kappa_F_max);
  else
    histo_limit_ext  = histo_limit;

  if(useSpline){
    cout<<"Evaluating the splines "<<endl;
    for(int i = 1; i <= histo_limit_ext->GetNbinsX(); i++){
      for(int j = 1; j <= histo_limit_ext->GetNbinsY(); j++){
	xvar->setVal(histo_limit_ext->GetXaxis()->GetBinCenter(i));
	yvar->setVal(histo_limit_ext->GetYaxis()->GetBinCenter(j));
	histo_limit_ext->SetBinContent(i,j,spline->getVal());
      }
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas",700,650);
  canvas->SetRightMargin(0.17);
  canvas->SetLeftMargin(0.11);
  TH1* frame = canvas->DrawFrame(kappa_V_min,kappa_F_min,kappa_V_max,kappa_F_max, "");
  frame->GetXaxis()->SetTitle("#kappa_{V}");
  frame->GetYaxis()->SetTitle("#kappa_{F}");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetXaxis()->SetLabelSize(0.90*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetLabelSize(0.90*frame->GetYaxis()->GetLabelSize());
  frame->GetZaxis()->SetLabelSize(0.90*frame->GetZaxis()->GetLabelSize());
  frame->Draw();

  
  histo_limit_ext->GetZaxis()->SetRangeUser(0.05,0.7);
  histo_limit_ext->GetZaxis()->SetLabelSize(0.90*frame->GetZaxis()->GetLabelSize());
  histo_limit_ext->Draw("COLZ SAME");
  TLine* line_vertical = new TLine(1,kappa_F_min,1,kappa_F_max);
  line_vertical->SetLineWidth(2);
  line_vertical->SetLineStyle(7);
  line_vertical->SetLineColor(kWhite);
  line_vertical->Draw("L same");

  TLine* line_horizontal = new TLine(kappa_V_min,1,kappa_V_max,1);
  line_horizontal->SetLineWidth(2);
  line_horizontal->SetLineStyle(7);
  line_horizontal->SetLineColor(kWhite);
  line_horizontal->Draw("L same");

  sm_value->Draw("P same");
  grlhc_68->Draw("L same");
  grlhc_95->Draw("L same");
  grlhc_bf->Draw("Psame");

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.12);

  TLegend *leg =  new TLegend(0.15,0.61,0.39,0.80);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(grlhc_bf,"LHC best fit","P");
  leg->AddEntry(grlhc_68,"68% CL","L");
  leg->AddEntry(grlhc_95,"95% CL","L");
  leg->AddEntry(sm_value,"SM production","P");
  leg->Draw("same");

  TLatex *   tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.040);
  tex->SetTextAngle(90);
  tex->DrawLatex(0.978,0.20,"95% CL upper limit on #sigma #times BR(H#rightarrow inv.)/#sigma_{SM}");
  
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/kappa_limit_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/kappa_limit_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/kappa_limit_"+postfix+".root").c_str(),"root");

  system(("rm "+string(file_temp->GetName())).c_str());
  if(file_temp) delete file_temp;

}
