#include "CMS_lumi.h"

using namespace std;

void prepostGJ(string fitFilename, string templateFileName, string observable, int category) {

  gROOT->SetBatch(kTRUE);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.295);
  pad2->SetTickx();
  pad2->SetTicky();
  
  TFile* pfile = new TFile(fitFilename.c_str());
  TFile* dfile = new TFile(templateFileName.c_str());
  
  TH1* dthist = (TH1*)dfile->Get(("datahistgam_"+observable).c_str());
  TH1* wlhist = (TH1*)pfile->Get("shapes_fit_b/ch4/QCD_GJ");
  TH1* pohist = (TH1*)pfile->Get("shapes_fit_b/ch4/total_background");    
  TH1* prhist = (TH1*)pfile->Get("shapes_prefit/ch4/total_background");    

  dthist->Scale(1.0, "width");

  ofstream outputfile;
  outputfile.open("prepostGJ.txt");

  stringstream QCDRate;
  QCDRate << "Process: QCD";
  stringstream PreRate;
  PreRate << "Process: Pre-fit (total)";
  stringstream PostRate;
  PostRate << "Process: Post-fit (total)";
  stringstream DataRate;
  DataRate << "Process: Data";

  for(int iBin = 0; iBin < wlhist->GetNbinsX(); iBin++){
    QCDRate << "   ";
    QCDRate << wlhist->GetBinContent(iBin);
  }

  for(int iBin = 0; iBin < prhist->GetNbinsX(); iBin++){
    PreRate << "   ";
    PreRate << prhist->GetBinContent(iBin);
  }

  for(int iBin = 0; iBin < pohist->GetNbinsX(); iBin++){
    PostRate << "   ";
    PostRate << pohist->GetBinContent(iBin);
  }

  for(int iBin = 0; iBin < dthist->GetNbinsX(); iBin++){
    DataRate << "   ";
    DataRate << dthist->GetBinContent(iBin);
  }

  outputfile<<"######################"<<endl;
  outputfile<<QCDRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PreRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PostRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<DataRate.str()<<endl;
  outputfile<<"######################"<<endl;

  outputfile.close();

  
  prhist->SetLineColor(kRed);
  prhist->SetLineWidth(2);
  pohist->SetLineColor(kBlue);
  pohist->SetLineWidth(2);
  prhist->SetMarkerColor(kRed);
  pohist->SetMarkerColor(kBlue);
  
  wlhist->SetFillColor(kOrange+1);
  wlhist->SetLineColor(kBlack);

  pad1->SetRightMargin(0.075);
  pad1->SetLeftMargin(0.10);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
    
  TH1* frame = NULL;
  if(category <=1)
   frame = pad1->DrawFrame(200., 0.0005, 1250., 10000, "");
  else
   frame = pad1->DrawFrame(250., 0.0005, 1000., 10000, "");

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.05);
  
  frame ->Draw();
  
  CMS_lumi(pad1, 4, 0, true);
  prhist->Draw("HIST SAME");
  pohist->Draw("HIST SAME");
  wlhist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1);
  dthist->SetLineWidth(2);
  dthist->SetMarkerStyle(20);
  
  dthist->Draw("PE SAME");
  
  TLegend* leg = new TLegend(0.5, 0.65, 0.75, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","P");
  leg->AddEntry(prhist, "Pre-fit #gamma+jet","L");
  leg->AddEntry(pohist, "Post-fit #gamma+jet","L");
  leg->AddEntry(wlhist, "Background", "F");
  leg->Draw("SAME");
  
  pad1->RedrawAxis("sameaxis");
  pad1->SetLogy();
  
  canvas->cd();
  pad2->SetTopMargin(0.08);
  pad2->SetRightMargin(0.075);
  pad2->SetLeftMargin(0.10);    
  pad2->SetGridy();
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();

  TH1* frame2 = NULL;
  if(category <=1)
   frame2 = pad2->DrawFrame(200., 0., 1250., 2., "");
  else
   frame2 = pad2->DrawFrame(250., 0., 1000., 2., "");

  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetYaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->Draw();
  
  TH1* d1hist = (TH1*)dthist->Clone("d1hist");
  TH1* d2hist = (TH1*)dthist->Clone("d2hist");
  TH1* m1hist = (TH1*)prhist->Clone("m1hist");
  TH1* m2hist = (TH1*)pohist->Clone("m2hist");
  TH1* erhist = (TH1*)pohist->Clone("erhist");
  
  d1hist->SetLineColor(kRed);
  d1hist->SetLineWidth(2);
  d2hist->SetLineColor(kBlue);
  d2hist->SetLineWidth(2);
  d1hist->SetMarkerColor(kRed);
  d1hist->SetMarkerSize(1);
  d1hist->SetMarkerStyle(20);
  d2hist->SetMarkerColor(kBlue);
  d2hist->SetMarkerSize(1);
  d2hist->SetMarkerStyle(20);
  
  for (int i = 1; i <= m1hist->GetNbinsX(); i++) m1hist->SetBinError(i, 0);
  for (int i = 1; i <= m2hist->GetNbinsX(); i++) m2hist->SetBinError(i, 0);
  
  d1hist->Divide(m1hist);
  d2hist->Divide(m2hist);
  erhist->Divide(m2hist);
  erhist->SetLineColor(0);
  erhist->SetMarkerColor(0);
  erhist->SetMarkerSize(0);
  erhist->SetFillColor(kGray);
  
  d1hist->SetMarkerSize(0.7);
  d2hist->SetMarkerSize(0.7);
  d1hist->SetStats(kFALSE);
  
  d1hist->GetXaxis()->SetLabelOffset(999999);
  d1hist->GetXaxis()->SetLabelSize(0);
  d1hist->GetXaxis()->SetTitleOffset(999999);
  d1hist->GetXaxis()->SetTitleSize(0);
  

  d1hist->GetYaxis()->SetTitleOffset(0.3);
  d1hist->GetYaxis()->SetLabelSize(0.12);
  d1hist->GetYaxis()->SetTitleSize(0.15);
  d1hist->GetYaxis()->SetTitle("Data/Pred.");
  
  d1hist->Draw("PE SAME");    
  d2hist->Draw("PE SAME");
  erhist->Draw("E2 SAME");
  d1hist->Draw("PE SAME");
  d2hist->Draw("PE SAME");

  TH1* unhist = (TH1*)pohist->Clone("unhist");

  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetFillColor(0);
  unhist->Draw("hist same");  
  pad2->RedrawAxis("sameaxis");
  
  canvas->SaveAs("prepostfit_gam.pdf");
  canvas->SaveAs("prepostfit_gam.png");
  
  TH1* postcorr = (TH1*)m2hist->Clone("postcorr");
  postcorr->Divide(m1hist);
  
  TFile* pcfile = new TFile("postcorrgam.root", "RECREATE");
  postcorr->Write();
  pcfile->Close();
}

