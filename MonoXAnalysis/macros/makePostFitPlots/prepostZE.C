#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void prepostZE(string fitFilename, string templateFileName, string observable, Category category,bool plotSBFit = false, bool dumpHisto = false) {

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

  TH1* dthist  = NULL;
  TH1* zllhist = NULL;
  TH1* wlhist  = NULL;
  TH1* tthist  = NULL;
  TH1* dihist  = NULL;
  TH1* ewkwhist  = NULL;
  TH1* ewkzhist  = NULL;
  TH1* pohist  = NULL;
  TH1* prhist  = NULL;

  if(!plotSBFit){
    
    dthist = (TH1*)dfile->FindObjectAny(("datahistzee_"+observable).c_str());
    zllhist = (TH1*)pfile->Get("shapes_fit_b/ch5/Znunu");
    wlhist = (TH1*)pfile->Get("shapes_fit_b/ch5/WJets_ZE");
    tthist = (TH1*)pfile->Get("shapes_fit_b/ch5/Top");
    dihist = (TH1*)pfile->Get("shapes_fit_b/ch5/Dibosons");
    ewkwhist = (TH1*)pfile->Get("shapes_fit_b/ch5/EWKW");
    ewkzhist = (TH1*)pfile->Get("shapes_fit_b/ch5/EWKZ");
    pohist = (TH1*)pfile->Get("shapes_fit_b/ch5/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch5/total_background");

  }
  else{

    dthist = (TH1*)dfile->FindObjectAny(("datahistzee_"+observable).c_str());
    zllhist = (TH1*)pfile->Get("shapes_fit_s/ch5/Znunu");
    wlhist = (TH1*)pfile->Get("shapes_fit_s/ch5/WJets_ZE");
    tthist = (TH1*)pfile->Get("shapes_fit_s/ch5/Top");
    ewkwhist = (TH1*)pfile->Get("shapes_fit_s/ch5/EWKW");
    ewkzhist = (TH1*)pfile->Get("shapes_fit_s/ch5/EWKZ");
    dihist = (TH1*)pfile->Get("shapes_fit_s/ch5/Dibosons");
    pohist = (TH1*)pfile->Get("shapes_fit_s/ch5/total_background");
    prhist = (TH1*)pfile->Get("shapes_prefit/ch5/total_background");

  }
  dthist->Scale(1.0, "width");
 
  ofstream  outputfile;
  outputfile.open("prepostZE.txt");
  stringstream TopRate;
  TopRate << "Process: Top";
  stringstream VVRate;
  VVRate << "Process: DiBoson";
  stringstream EWKWRate;
  EWKWRate << "Process: EWKW";
  stringstream EWKZRate;
  EWKZRate << "Process: EWKZ";
  stringstream WJetRate;
  WJetRate << "Process: WJet";
  stringstream ZJetRate;
  ZJetRate << "Process: ZJet";
  stringstream PreRate;
  PreRate << "Process: Pre-fit (total)";
  stringstream PostRate;
  PostRate << "Process: Post-fit (total)";
  stringstream DataRate;
  DataRate << "Process: Data";

  for(int iBin = 1; iBin <= tthist->GetNbinsX(); iBin++){
    TopRate << "   ";
    TopRate << tthist->GetBinContent(iBin) <<" pm "<<tthist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= dihist->GetNbinsX(); iBin++){
    VVRate <<"   ";
    VVRate <<dihist->GetBinContent(iBin) << " pm "<<dihist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= ewkwhist->GetNbinsX(); iBin++){
    EWKWRate <<"   ";
    EWKWRate <<ewkwhist->GetBinContent(iBin) << " pm "<<ewkwhist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= ewkzhist->GetNbinsX(); iBin++){
    EWKZRate <<"   ";
    EWKZRate <<ewkzhist->GetBinContent(iBin) << " pm "<<ewkzhist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= wlhist->GetNbinsX(); iBin++){
    WJetRate << "   ";
    WJetRate << wlhist->GetBinContent(iBin)<< " pm "<<wlhist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= zllhist->GetNbinsX(); iBin++){
    ZJetRate << "   ";
    ZJetRate << zllhist->GetBinContent(iBin)<< " pm "<<zllhist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= prhist->GetNbinsX(); iBin++){
    PreRate << "   ";
    PreRate << prhist->GetBinContent(iBin) << " pm "<<prhist->GetBinError(iBin);
  }

  for(int iBin = 1; iBin <= pohist->GetNbinsX(); iBin++){
    PostRate << "   ";
    PostRate << pohist->GetBinContent(iBin) <<" pm "<<pohist->GetBinError(iBin);
  }  

  for(int iBin = 1; iBin <= dthist->GetNbinsX(); iBin++){
    DataRate << "   ";
    DataRate << dthist->GetBinContent(iBin) << " pm "<<dthist->GetBinError(iBin);
  }

  outputfile<<"######################"<<endl;
  outputfile<<TopRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<VVRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<EWKWRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<EWKZRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<WJetRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<ZJetRate.str()<<endl;
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
  wlhist->Add(ewkwhist);
  wlhist->Add(ewkzhist);

  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.001,wlhist->GetMaximum()*100);
  else
    frame->GetYaxis()->SetRangeUser(0.0007,wlhist->GetMaximum()*100);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);  
  if(category  == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);


  frame ->Draw();
  
  CMS_lumi(canvas,"2.6");
  prhist->Draw("HIST SAME");
  pohist->Draw("HIST SAME");
  wlhist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  
  dthist->Draw("EP SAME");
  
  TLegend* leg = new TLegend(0.5, 0.65, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","PEL");
  leg->AddEntry(pohist, "Post-fit (Z #rightarrow ee)","L");
  leg->AddEntry(prhist, "Pre-fit (Z #rightarrow ee)","L");
  leg->AddEntry(wlhist, "Other Backgrounds", "F");
  leg->Draw("SAME");
  
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();  
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  (TH1*) dthist->Clone("frame");
  frame2->Reset("ICES");

  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.25,1.75);
  else
    frame2->GetYaxis()->SetRangeUser(0.25,1.75);

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

  d1hist->GetYaxis()->SetTitleOffset(1.5);
  d1hist->GetYaxis()->SetLabelSize(0.03);
  d1hist->GetYaxis()->SetTitleSize(0.04);
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

  canvas->Update();
  
  canvas->SaveAs("prepostfit_zee.pdf");
  canvas->SaveAs("prepostfit_zee.png");

  if(dumpHisto){

    TFile* outFile = new TFile("postfit_weights_ZE.root","RECREATE");
    outFile->cd();
    
    dthist->Write("data");
    zllhist->Write("zjets_post_fit");
    wlhist->Write("wjets_post_fit");
    tthist->Write("top_post_fit");
    dihist->Write("diboson_post_fit");

    TH1* zllhist_prefit = (TH1*) pfile->Get("shapes_prefit/ch5/Znunu");
    TH1* wlhist_prefit = (TH1*) pfile->Get("shapes_prefit/ch5/WJets_ZE");
    TH1* tthist_prefit = (TH1*) pfile->Get("shapes_prefit/ch5/Top");
    TH1* dihist_prefit = (TH1*) pfile->Get("shapes_prefit/ch5/Dibosons");
    
    zllhist_prefit->Write("zjets_pre_fit");
    wlhist_prefit->Write("wjets_pre_fit");
    tthist_prefit->Write("top_pre_fit");
    dihist_prefit->Write("diboson_pre_fit");

    pohist->Write("post_fit_post_fit");
    prhist->Write("pre_fit_post_fit");      


    outFile->Close();
  }
}

