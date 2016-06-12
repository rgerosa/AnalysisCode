#include "../CMS_lumi.h"

void drawUpperPlot(TPad* pad1, TH1* histo_1, TH1* histo_2){
 
  TH1* frame = pad1->DrawFrame(histo_1->GetXaxis()->GetBinLowEdge(1),0.0001,histo_1->GetXaxis()->GetBinLowEdge(histo_1->GetNbinsX()+1),1.,"");

  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->GetYaxis()->SetTitle("arbitrary unit");
  frame ->Draw();
  CMS_lumi(pad1,"2.30",true);
  
  histo_1->SetLineColor(kBlack);
  histo_1->SetLineWidth(2);

  histo_2->SetLineColor(kRed);
  histo_2->SetMarkerColor(kRed);
  histo_2->SetMarkerStyle(20);
  histo_2->SetMarkerSize(0.8);
  
  histo_1->Draw("HIST SAME");
  histo_2->Draw("P SAME");

  TLegend* leg = new TLegend(0.5,0.65,0.93,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(histo_1,"Private  Production m_{DM} = 100 GeV","L");
  leg->AddEntry(histo_2,"Official Production m_{DM} = 100 GeV","L");
  leg->Draw("same");

  return ;
}

void drawDownPlot(TPad* pad2, TH1* histo_1, TH1* histo_2,string xAxisTitle){

  
  TH1* frame = pad2->DrawFrame(histo_1->GetXaxis()->GetBinLowEdge(1),0,histo_1->GetXaxis()->GetBinLowEdge(histo_1->GetNbinsX()+1),2.0,"");


  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.10);
  frame->GetXaxis()->SetLabelOffset(0.03);
  frame->GetXaxis()->SetTitleSize(0.13);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetLabelSize(0.08);
  frame->GetYaxis()->SetTitleSize(0.10);
  frame->GetYaxis()->SetTitle("Ratio");
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetYaxis()->SetNdivisions(504, false);
  frame->GetYaxis()->SetTitleOffset(0.5);
  frame->Draw();

  TH1* ratio = (TH1*) histo_1->Clone("ratio");
  ratio->Divide(histo_2);

  ratio->SetLineColor(kRed);
  ratio->SetMarkerColor(kRed);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.8);

  ratio->Draw("Psame");

  return ;

}

void makeValidationPrivateGeneration(string outputDirectory){


  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  system(("mkdir "+outputDirectory).c_str());

  TFile* inputFile_1 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoW_Scalar/tree_DM_ScalarWH_Mphi-100_Mchi-50_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
  TFile* inputFile_2 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-19-2-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-100_Mchi-50_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");

  TTree* tree_1 = (TTree*) inputFile_1->Get("tree/tree");
  TTree* tree_2 = (TTree*) inputFile_2->Get("tree/tree");

  TH1F* jetPt_1 = new TH1F("jetPt_1","",30,100,1000);
  TH1F* jetPt_2 = new TH1F("jetPt_2","",30,100,1000);
  jetPt_1->Sumw2();
  jetPt_2->Sumw2();

  TH1F* met_1 = new TH1F("met_1","",30,200,1000);
  TH1F* met_2 = new TH1F("met_2","",30,200,1000);
  met_1->Sumw2();
  met_2->Sumw2();

  TH1F* nvtx_1 = new TH1F("nvtx_1","",25,0,25);
  TH1F* nvtx_2 = new TH1F("nvtx_2","",25,0,25);
  nvtx_1->Sumw2();
  nvtx_2->Sumw2();

  TH1F* njets_1 = new TH1F("njets_1","",7,0,7);
  TH1F* njets_2 = new TH1F("njets_2","",7,0,7);
  njets_1->Sumw2();
  njets_2->Sumw2();

  TH1F* mpruned_1 = new TH1F("mpruned_1","",25,0,120);
  TH1F* mpruned_2 = new TH1F("mpruned_2","",25,0,120);
  mpruned_1->Sumw2();
  mpruned_2->Sumw2();

  TH1F* tau2tau1_1 = new TH1F("tau2tau1_1","",25,0,1);
  TH1F* tau2tau1_2 = new TH1F("tau2tau1_2","",25,0,1);
  tau2tau1_1->Sumw2();
  tau2tau1_2->Sumw2();

  string cutString = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && (hltmet90 > 0 || hltmet120 > 0 || hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0) && njets >= 1 && nbjetslowpt < 1 && t1pfmet > 200";

  string cutString_2 = "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && (hltmet90 > 0 || hltmet120 > 0 || hltmetwithmu120 > 0 || hltmetwithmu170 > 0 || hltmetwithmu300 > 0 || hltmetwithmu90 > 0) && njets >= 1 && t1pfmet > 200 && boostedJetpt[0] > 200";

  tree_1->Draw("combinejetpt[0] >> jetPt_1",cutString.c_str(),"goff");
  tree_2->Draw("combinejetpt[0] >> jetPt_2",cutString.c_str(),"goff");

  tree_1->Draw("t1pfmet >> met_1",cutString.c_str(),"goff");
  tree_2->Draw("t1pfmet >> met_2",cutString.c_str(),"goff");

  tree_1->Draw("nvtx >> nvtx_1",cutString.c_str(),"goff");
  tree_2->Draw("nvtx >> nvtx_2",cutString.c_str(),"goff");

  tree_1->Draw("njets >> njets_1",cutString.c_str(),"goff");
  tree_2->Draw("njets >> njets_2",cutString.c_str(),"goff");

  tree_1->Draw("prunedJetm[0] >> mpruned_1",cutString_2.c_str(),"goff");
  tree_2->Draw("prunedJetm[0] >> mpruned_2",cutString_2.c_str(),"goff");

  tree_1->Draw("boostedJettau2[0]/boostedJettau1[0] >> tau2tau1_1",cutString_2.c_str(),"goff");
  tree_2->Draw("boostedJettau2[0]/boostedJettau1[0] >> tau2tau1_2",cutString_2.c_str(),"goff");


  // normalize histo to compare shapes
  jetPt_1->Scale(1./jetPt_1->Integral());
  jetPt_2->Scale(1./jetPt_2->Integral());

  met_1->Scale(1./met_1->Integral());
  met_2->Scale(1./met_2->Integral());
  
  nvtx_1->Scale(1./nvtx_1->Integral());
  nvtx_2->Scale(1./nvtx_2->Integral());

  njets_1->Scale(1./njets_1->Integral());
  njets_2->Scale(1./njets_2->Integral());

  mpruned_1->Scale(1./mpruned_1->Integral());
  mpruned_2->Scale(1./mpruned_2->Integral());

  tau2tau1_1->Scale(1./tau2tau1_1->Integral());  
  tau2tau1_2->Scale(1./tau2tau1_2->Integral());  

  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.13);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.28);
  pad2->SetTickx();
  pad2->SetTicky();

  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  drawUpperPlot(pad1,jetPt_1,jetPt_2);
  pad1->SetLogy();
  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  drawDownPlot(pad2,jetPt_1,jetPt_2,"jet p_{T} [GeV]");

  canvas->SaveAs((outputDirectory+"/jetPt.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/jetPt.png").c_str(),"png");

  canvas->cd();
  pad1->cd();
  drawUpperPlot(pad1,met_1,met_2);
  canvas->cd();
  pad2->cd();
  drawDownPlot(pad2,met_1,met_2,"E_{T}^{miss} (GeV)");

  canvas->SaveAs((outputDirectory+"/met.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/met.png").c_str(),"png");

  canvas->cd();
  pad1->cd();
  drawUpperPlot(pad1,nvtx_1,nvtx_2);
  pad1->SetLogy(0);
  canvas->cd();
  pad2->cd();
  drawDownPlot(pad2,nvtx_1,nvtx_2,"N_{PV}");

  canvas->SaveAs((outputDirectory+"/nvtx.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/nvtx.png").c_str(),"png");

  canvas->cd();
  pad1->cd();
  drawUpperPlot(pad1,njets_1,njets_2);
  pad1->SetLogy(0);
  canvas->cd();
  pad2->cd();
  drawDownPlot(pad2,njets_1,njets_2,"N_{jets}");

  canvas->SaveAs((outputDirectory+"/njet.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/njet.png").c_str(),"png");

  canvas->cd();
  pad1->cd();
  pad1->SetLogy();
  drawUpperPlot(pad1,mpruned_1,mpruned_2);
  canvas->cd();
  pad2->cd();
  drawDownPlot(pad2,mpruned_1,mpruned_2,"m_{pruned} [GeV]");

  canvas->SaveAs((outputDirectory+"/mpruned.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/mpruned.png").c_str(),"png");

  canvas->cd();
  pad1->cd();
  drawUpperPlot(pad1,tau2tau1_1,tau2tau1_2);
  canvas->cd();
  pad2->cd();
  drawDownPlot(pad2,mpruned_1,mpruned_2,"#tau_{2}/#tau_{1}");

  canvas->SaveAs((outputDirectory+"/tau2tau1.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/tau2tau1.png").c_str(),"png");

}
