#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void prepostGJ(string fitFilename, string templateFileName, string observable, Category category,bool plotSBFit = false,  bool dumpHisto = false) {

  gROOT->SetBatch(kTRUE); 
  setTDRStyle();
  

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  

  TFile* pfile = new TFile(fitFilename.c_str());
  TFile* dfile = new TFile(templateFileName.c_str());

  TH1* dthist = NULL;
  TH1* qchist = NULL;
  TH1* vghist = NULL;
  TH1* vlhist = NULL;
  TH1* pohist = NULL;
  TH1* prhist = NULL;

  if(!plotSBFit){
    
    dthist = (TH1*)dfile->FindObjectAny(("datahistgam_"+observable).c_str());
    qchist = (TH1*)pfile->Get("shapes_fit_b/ch4/QCD_GJ");
    vghist = (TH1*)pfile->Get("shapes_fit_b/ch4/VGamma_GJ");
    vlhist = (TH1*)pfile->Get("shapes_fit_b/ch4/WJets_GJ");
    pohist = (TH1*)pfile->Get("shapes_fit_b/ch4/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch4/total_background");

  }
  else{

    dthist = (TH1*)dfile->FindObjectAny(("datahistgam_"+observable).c_str());
    qchist = (TH1*)pfile->Get("shapes_fit_s/ch4/QCD_GJ");
    vghist = (TH1*)pfile->Get("shapes_fit_s/ch4/VGamma_GJ");
    vlhist = (TH1*)pfile->Get("shapes_fit_s/ch4/WJets_GJ");
    pohist = (TH1*)pfile->Get("shapes_fit_s/ch4/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch4/total_background");

  }
  dthist->Scale(1.0, "width");

  ofstream  outputfile;
  outputfile.open("prepostGJ.txt");
  stringstream QCDRate;
  QCDRate << "Process: QCD";
  stringstream PreRate;
  PreRate << "Process: Pre-fit (total)";
  stringstream PostRate;
  PostRate << "Process: Post-fit (total)";
  stringstream DataRate;
  DataRate << "Process: Data";

  for(int iBin = 0; iBin < qchist->GetNbinsX(); iBin++){
    QCDRate << "   ";
    QCDRate << qchist->GetBinContent(iBin);
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
  
  outputfile.close();
  
  
  prhist->SetLineColor(kRed);
  prhist->SetLineWidth(2);
  pohist->SetLineColor(kBlue);
  pohist->SetLineWidth(2);
  prhist->SetMarkerColor(kRed);
  pohist->SetMarkerColor(kBlue);
  
  qchist->Add(vghist);
  qchist->Add(vlhist);
  qchist->SetFillColor(kOrange+1);
  vghist->SetFillColor(kOrange+1);
  vlhist->SetFillColor(kOrange+1);
  qchist->SetLineColor(kBlack);
  
  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,prhist->GetMaximum()*100);
  else
    frame->GetYaxis()->SetRangeUser(0.0007,prhist->GetMaximum()*100);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);


  frame ->Draw();

  CMS_lumi(canvas,"12.9");
  
  prhist->Draw("HIST SAME");
  pohist->Draw("HIST SAME");
  qchist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetLineColor(kBlack);
  dthist->Draw("EP SAME");
  
  TLegend* leg = new TLegend(0.6, 0.60, 0.92, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","PEL");
  leg->AddEntry(pohist, "Post-fit (#gamma + jets)","L");
  leg->AddEntry(prhist, "Pre-fit (#gamma + jets)","L");
  leg->AddEntry(qchist, "Other Backgrounds", "F");
  leg->Draw("SAME");


  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  (TH1*) dthist->Clone("frame");
  frame2->Reset("ICES");

  if(category ==  Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(210);
  frame2->GetYaxis()->SetNdivisions(5);


  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
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
  unhist->SetLineWidth(2);
  unhist->SetFillColor(0);
  unhist->Draw("hist same");  
  pad2->RedrawAxis("G sameaxis");
  
  canvas->SaveAs("prepostfit_gam.pdf");
  canvas->SaveAs("prepostfit_gam.png");

  if(dumpHisto){

    TFile* outFile = new TFile("postfit_weights_GJ.root","RECREATE");
    outFile->cd();

    dthist->Write("data");
    qchist->Write("qcd_post_fit");
    TH1* qchist_prefit = (TH1*) pfile->Get("shapes_prefit/ch4/QCD_GJ");
    qchist_prefit->Write("wjets_pre_fit");
    pohist->Write("post_fit_post_fit");
    prhist->Write("pre_fit_post_fit");


    outFile->Close();
  }
}

