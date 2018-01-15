#include "../CMS_lumi.h"
#include "../../../../HiggsAnalysis/CombinedLimit/interface/RooSplineND.h"

static int nbinsX = 50;
static int nbinsY = 50;
static int reductionForContour = 10;

TGraph* produceContour (const int & reduction){

  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
        if(i2%reduction != 0) continue; // reduce number of points                                                                                                                                     
        lXE.push_back(pCurv->GetX()[i2]);
        lYE.push_back(pCurv->GetY()[i2]);
      }
      pCurv->SetLineColor(kRed);
      pCurv = (TGraph*)pContLevel->After(pCurv);
    }
  }
  if(lXE.size() == 0) {
    lXE.push_back(0);
    lYE.push_back(0);
  }

  TGraph *lTotalE = new TGraph(lXE.size(),&lXE[0],&lYE[0]);
  return lTotalE;
}

void makeKappaVKappaFBestFit(string inputFileName, string outputDIR, string postfix){

  gSystem->Load("$CMSSW_BASE/lib/slc6_amd64_gcc530/libHiggsAnalysisCombinedLimit.so");

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  //gStyle->SetPalette(kBeach);
  //gStyle->SetPalette(kBrownCyan);
  //gStyle->SetPalette(kLake);
  gStyle->SetPalette(kLightTemperature);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();
  TTree* tree_limit = (TTree*) inputFile->Get("limit");

  TTreeReader reader (tree_limit);
  TTreeReaderValue<float> kappa_V (reader,"kappa_V");
  TTreeReaderValue<float> kappa_F (reader,"kappa_F");
  TTreeReaderValue<float> deltaNLL (reader,"deltaNLL");

  TGraph2D* graph_limit = new TGraph2D();
  TGraph* bestfit = new TGraph();
  int    npoints = 0;

  float  kappa_V_min = 100;
  float  kappa_V_max = 0;
  float  kappa_F_min = 100;
  float  kappa_F_max = 0;

  while(reader.Next()){    
    if(npoints == 0) bestfit->SetPoint(npoints,*kappa_V,*kappa_F);
    graph_limit->SetPoint(npoints,*kappa_V,*kappa_F,2*(*deltaNLL));
    npoints++;   
    if(*kappa_V > kappa_V_max) kappa_V_max = *kappa_V;
    if(*kappa_V < kappa_V_min) kappa_V_min = *kappa_V;
    if(*kappa_F > kappa_F_max) kappa_F_max = *kappa_F;
    if(*kappa_F < kappa_F_min) kappa_F_min = *kappa_F;
  }

  TH2F* histo_limit = new TH2F("histo_limit","",nbinsX,kappa_V_min,kappa_V_max,nbinsY,kappa_F_min,kappa_F_max);
  for(int iBinX = 0; iBinX < histo_limit->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < histo_limit->GetNbinsY(); iBinY++){
      histo_limit->SetBinContent(iBinX+1,iBinY+1,graph_limit->Interpolate(histo_limit->GetXaxis()->GetBinCenter(iBinX+1),histo_limit->GetYaxis()->GetBinCenter(iBinY+1)));
    }
  }
  histo_limit->Smooth();

  RooRealVar* xvar = new RooRealVar("x","x",(kappa_V_max-kappa_V_min)/2,kappa_V_min,kappa_V_max);
  RooRealVar* yvar = new RooRealVar("y","y",(kappa_F_max-kappa_F_min)/2,kappa_F_min,kappa_F_max);

  TTree *tree = new TTree("tree","tree"); 
  float  x,y,z;
  tree->Branch("x",&x,"x/F");
  tree->Branch("y",&y,"y/F");
  tree->Branch("z",&z,"z/F");

  
  for(int i = 1; i <= histo_limit->GetNbinsX(); i++){
    for(int j = 1; j <= histo_limit->GetNbinsY(); j++){
      x = histo_limit->GetXaxis()->GetBinCenter(i);
      y = histo_limit->GetYaxis()->GetBinCenter(j);
      z = histo_limit->GetBinContent(i,j);
      tree->Fill();
    }
  }

  RooArgList list (*xvar,*yvar);
  RooSplineND *spline = new RooSplineND(Form("spline_%s",postfix.c_str()),Form("spline_%s",postfix.c_str()),list,tree,"z",1,true);  
  if(tree) delete tree;
  
  TH2F* histo_limit_ext = new TH2F("histo_limit_ext","",nbinsX*10,kappa_V_min,kappa_V_max,nbinsY*10,kappa_F_min,kappa_F_max);

  cout<<"Evaluating the splines "<<endl;
  for(int i = 1; i <= histo_limit_ext->GetNbinsX(); i++){
    for(int j = 1; j <= histo_limit_ext->GetNbinsY(); j++){
      xvar->setVal(histo_limit_ext->GetXaxis()->GetBinCenter(i));
      yvar->setVal(histo_limit_ext->GetYaxis()->GetBinCenter(j));
      histo_limit_ext->SetBinContent(i,j,spline->getVal());
    }
  }

  TH2* histo_cl68 = (TH2*) histo_limit_ext->Clone("histo_cl68");
  TH2* histo_cl95 = (TH2*) histo_limit_ext->Clone("histo_cl95");

  double contours[1]; contours[0]=1;
  histo_cl68->SetContour(1,contours);
  contours[0]=3.9; // 95% CL
  histo_cl95->SetContour(1,contours);

  TCanvas* canvas = new TCanvas("canvas", "canvas",700,650);
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.12);
  TH1* frame = canvas->DrawFrame(kappa_V_min,kappa_F_min,kappa_V_max,kappa_F_max, "");
  frame->GetXaxis()->SetTitle("#kappa_{V}");
  frame->GetYaxis()->SetTitle("#kappa_{F}");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetLabelSize(0.90*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetLabelSize(0.90*frame->GetYaxis()->GetLabelSize());
  frame->GetZaxis()->SetLabelSize(0.90*frame->GetZaxis()->GetLabelSize());
  frame->Draw();
  
  histo_cl68->GetZaxis()->SetLabelSize(0);
  histo_cl68->Draw("contz list same");
  canvas->Update();
  TGraph* contour_cl68 = produceContour(1);

  histo_cl95->GetZaxis()->SetLabelSize(0);
  histo_cl95->Draw("contz list same");
  canvas->Update();
  TGraph* contour_cl95 = produceContour(reductionForContour);

  frame->Draw();
  histo_limit_ext->GetZaxis()->SetRangeUser(0,50);
  histo_limit_ext->GetZaxis()->SetLabelSize(0.90*frame->GetZaxis()->GetLabelSize());
  histo_limit_ext->Draw("COLZ SAME");
  contour_cl95->SetLineColor(kBlack);
  contour_cl95->SetLineWidth(3);
  contour_cl95->SetLineStyle(7);
  contour_cl68->SetLineColor(kBlack);
  contour_cl68->SetLineWidth(3);
  contour_cl95->Draw("L same");
  contour_cl68->Draw("L same");
  bestfit->SetMarkerColor(kBlack);
  bestfit->SetMarkerSize(1.6);
  bestfit->SetMarkerStyle(34);
  bestfit->Draw("Psame");

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.1);

  TLegend *leg =  new TLegend(0.15,0.65,0.355,0.80);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(bestfit,"Best fit H_{inv}","P");
  leg->AddEntry(contour_cl68,"68% CL","L");
  leg->AddEntry(contour_cl95,"95% CL","L");
  leg->Draw("same");

  TLatex *   tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.040);
  tex->SetTextAngle(90);
  tex->DrawLatex(0.978,0.65,"-2 x #Delta Log(L)");
  
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/kappa_bestfit_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/kappa_bestfit_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/kappa_bestfit_"+postfix+".root").c_str(),"root");


}
