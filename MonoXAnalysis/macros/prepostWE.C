#include "CMS_lumi.h"

void prepostWE(string fitFilename, string templateFileName, string observable, int category,bool plotSBFit = false) {

  gROOT->SetBatch(kTRUE); 
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 675);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
   
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();
  
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.295);
  pad2->SetTickx();
  pad2->SetTicky();
  
  setTDRStyle();

  TFile* pfile = new TFile(fitFilename.c_str());
  TFile* dfile = new TFile(templateFileName.c_str());

  TH1* dthist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* pohist = NULL;
  TH1* prhist = NULL;

  if(!plotSBFit){
    
    dthist = (TH1*)dfile->Get(("datahistwen_"+observable).c_str());
    wlhist = (TH1*)pfile->Get("shapes_fit_b/ch6/ZJets_WE");
    tthist = (TH1*)pfile->Get("shapes_fit_b/ch6/Top");
    dihist = (TH1*)pfile->Get("shapes_fit_b/ch6/Dibosons");
    pohist = (TH1*)pfile->Get("shapes_fit_b/ch6/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch6/total_background");

  }
  else{

    dthist = (TH1*)dfile->Get(("datahistwen_"+observable).c_str());
    wlhist = (TH1*)pfile->Get("shapes_fit_s/ch6/ZJets_WE");
    tthist = (TH1*)pfile->Get("shapes_fit_s/ch6/Top");
    dihist = (TH1*)pfile->Get("shapes_fit_s/ch6/Dibosons");
    pohist = (TH1*)pfile->Get("shapes_fit_s/ch6/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch6/total_background");

  }
  dthist->Scale(1.0, "width");

  ofstream  outputfile;
  outputfile.open("prepostWE.txt");
  stringstream TopRate;
  TopRate << "Process: Top";
  stringstream VVRate;
  VVRate << "Process: DiBoson";
  stringstream WJetRate;
  WJetRate << "Process: WJet";
  stringstream PreRate;
  PreRate << "Process: Pre-fit (total)";
  stringstream PostRate;
  PostRate << "Process: Post-fit (total)";
  stringstream DataRate;
  DataRate << "Process: Data";

  for(int iBin = 0; iBin < tthist->GetNbinsX(); iBin++){
    TopRate << "   ";
    TopRate << tthist->GetBinContent(iBin);
  }

  for(int iBin = 0; iBin < dihist->GetNbinsX(); iBin++){
    VVRate << "   ";
    VVRate << dihist->GetBinContent(iBin);
  }

  for(int iBin = 0; iBin < wlhist->GetNbinsX(); iBin++){
    WJetRate << "   ";
    WJetRate << wlhist->GetBinContent(iBin);
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
  outputfile<<TopRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<VVRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<WJetRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PreRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PostRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<DataRate.str()<<endl;
  
  outputfile.close();
  
  
  prhist->SetLineColor(kRed);
  prhist->SetLineWidth(2);
  pohist->SetLineColor(kBlue);
  pohist->SetLineWidth(2);
  prhist->SetMarkerColor(kRed);
  pohist->SetMarkerColor(kBlue);
  
  wlhist->SetFillColor(kOrange+1);
  wlhist->SetLineColor(kBlack);
  wlhist->Add(tthist);
  wlhist->Add(dihist);

  
  pad1->SetRightMargin(0.06);
  pad1->SetLeftMargin(0.12);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();

  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category <=1)
    frame->GetYaxis()->SetRangeUser(0.0005,6000);
  else
    frame->GetYaxis()->SetRangeUser(0.0005,100);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.055);  
  frame ->Draw();
  
  CMS_lumi(pad1,"2.30");
  prhist->Draw("HIST SAME");
  pohist->Draw("HIST SAME");
  wlhist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  
  dthist->Draw("EP SAME");
  
  TLegend* leg = new TLegend(0.55, 0.55, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","PEL");
  leg->AddEntry(pohist, "Post-fit (W #rightarrow e#nu)","L");
  leg->AddEntry(prhist, "Pre-fit (W #rightarrow e#nu)","L");
  leg->AddEntry(wlhist, "Other Backgrounds", "F");
  leg->Draw("SAME");
  
  pad1->RedrawAxis("sameaxis");
  pad1->SetLogy();
  
  canvas->cd();
  pad2->SetTopMargin(0.08);
  pad2->SetRightMargin(0.06);
  pad2->SetLeftMargin(0.12);    
  pad2->SetBottomMargin(0.35);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  (TH1*) dthist->Clone("frame");
  frame2->Reset("ICES");

  if(category <=1)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else
    frame2->GetYaxis()->SetRangeUser(0.25,1.75);

  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetXaxis()->SetLabelSize(0.11);
  frame2->GetYaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetTitleSize(0.135);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->GetXaxis()->SetNdivisions(510);
  frame2->Draw();

  TH1* d1hist = (TH1*)dthist->Clone("d1hist");
  TH1* d2hist = (TH1*)dthist->Clone("d2hist");
  TH1* m1hist = (TH1*)prhist->Clone("m1hist");
  TH1* m2hist = (TH1*)pohist->Clone("m2hist");
  TH1* erhist = (TH1*)pohist->Clone("erhist");
  
  d1hist->SetLineColor(kRed);
  d2hist->SetLineColor(kBlue);
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
  
  d1hist->SetMarkerSize(1.2);
  d2hist->SetMarkerSize(1.2);
  d1hist->SetStats(kFALSE);
  
  d1hist->GetXaxis()->SetLabelOffset(999999);
  d1hist->GetXaxis()->SetLabelSize(0);
  d1hist->GetXaxis()->SetTitleOffset(999999);
  d1hist->GetXaxis()->SetTitleSize(0);

  for(int iBin = 1; iBin <= erhist->GetNbinsX(); iBin++){
    if(erhist->GetBinError(iBin) > erhist->GetBinError(iBin+1) && iBin != erhist->GetNbinsX())
      erhist->SetBinError(iBin,erhist->GetBinError(iBin+1)*0.9);
  }

  d1hist->GetYaxis()->SetTitleOffset(0.3);
  d1hist->GetYaxis()->SetLabelSize(0.12);
  d1hist->GetYaxis()->SetTitleSize(0.15);
  d1hist->GetYaxis()->SetTitle("Data/Pred.");
  
  d1hist->Draw("PE1 SAME");    
  d2hist->Draw("PE1 SAME");
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
  unhist->SetLineWidth(2);
  unhist->Draw("hist same");  
  pad2->RedrawAxis("G sameaxis");
  
  canvas->SaveAs("prepostfit_wen.pdf");
  canvas->SaveAs("prepostfit_wen.png");

}

