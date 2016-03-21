#include "CMS_lumi.h"

void prepostSig(string fitFilename, string templateFileName, string observable, int category, 
		bool isHiggsInvisible, int scaleSig = 1, bool blind = true, bool plotSBFit = false, 
		string interaction = "Vector", string mediatorMass = "2000", string DMMass = "10") {
  

  gROOT->SetBatch(kTRUE);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
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

  // in case of b-only fit just dispaly three possible signal on the stack
  TH1* mjhist = NULL;
  TH1* mwhist = NULL;
  TH1* mzhist = NULL;

  TH1* ggHhist = NULL;
  TH1* vbfhist = NULL;
  TH1* wHhist = NULL;
  TH1* zHhist = NULL;
  TH1* zhhist = NULL;
  

  if(! isHiggsInvisible){
    mjhist = (TH1*) dfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    mwhist = (TH1*) dfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    mzhist = (TH1*) dfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());  
    mjhist->Scale(1.0, "width");
    mwhist->Scale(1.0, "width");
    mzhist->Scale(1.0, "width");
  }
  else{
    ggHhist = (TH1*) dfile->FindObjectAny(("ggHhist_"+observable).c_str());
    vbfhist = (TH1*) dfile->FindObjectAny(("vbfHhist_"+observable).c_str());
    wHhist  = (TH1*) dfile->FindObjectAny(("wHhist_"+observable).c_str());
    zHhist  = (TH1*) dfile->FindObjectAny(("zHhist_"+observable).c_str());
    ggHhist->Scale(1.0, "width");
    vbfhist->Scale(1.0, "width");
    wHhist->Scale(1.0, "width");
    zHhist->Scale(1.0, "width");
  }

  TH1* znhist = NULL;
  TH1* zlhist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* qchist = NULL;
  TH1* gmhist = NULL;
  TH1* tohist = NULL;
  TH1* tphist = NULL;
  TH1* sighist = NULL;

  if(!plotSBFit){

    znhist = (TH1*)pfile->Get("shapes_fit_b/ch1/Znunu");    
    zlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/ZJets");    
    wlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/WJets");    
    tthist = (TH1*)pfile->Get("shapes_fit_b/ch1/Top");    
    dihist = (TH1*)pfile->Get("shapes_fit_b/ch1/Dibosons");    
    qchist = (TH1*)pfile->Get("shapes_fit_b/ch1/QCD");    
    gmhist = (TH1*)pfile->Get("shapes_fit_b/ch1/GJets");    
    tohist = (TH1*)pfile->Get("shapes_fit_b/ch1/total_background");    
    tphist = (TH1*)pfile->Get("shapes_prefit/ch1/total_background");    
    
  }
  else{
    znhist = (TH1*)pfile->Get("shapes_fit_s/ch1/Znunu");    
    zlhist = (TH1*)pfile->Get("shapes_fit_s/ch1/ZJets");    
    wlhist = (TH1*)pfile->Get("shapes_fit_s/ch1/WJets");    
    tthist = (TH1*)pfile->Get("shapes_fit_s/ch1/Top");    
    dihist = (TH1*)pfile->Get("shapes_fit_s/ch1/Dibosons");    
    qchist = (TH1*)pfile->Get("shapes_fit_s/ch1/QCD");    
    gmhist = (TH1*)pfile->Get("shapes_fit_s/ch1/GJets");    
    tohist = (TH1*)pfile->Get("shapes_fit_s/ch1/total_background");    
    tphist = (TH1*)pfile->Get("shapes_prefit/ch1/total_background");      
    sighist = (TH1*)pfile->Get("shapes_fit_s/ch1/total_signal");
  }

  TH1* dthist = NULL;
  if(!blind){
    dthist = (TH1*)dfile->FindObjectAny(("datahist_"+observable).c_str());
    dthist->Scale(1.0,"width");
  }
  else{
    dthist = (TH1*)dfile->FindObjectAny(("datahist_"+observable).c_str());
    for (int i = 0; i <= dthist->GetNbinsX(); i++) {
      double yield = 0.0;
      yield += zlhist->GetBinContent(i);
      yield += gmhist->GetBinContent(i);
      yield += wlhist->GetBinContent(i);
      yield += tthist->GetBinContent(i);
      yield += dihist->GetBinContent(i);
      yield += qchist->GetBinContent(i);
      yield += znhist->GetBinContent(i);
      dthist->SetBinContent(i, yield);
      dthist->SetBinError(i, 0.);
    }
  }

  // print in text file yields
  ofstream outputfile;
  outputfile.open("prepostSR.txt");

  stringstream QCDRate;
  QCDRate << "Process: QCD";
  stringstream GJetsRate;
  GJetsRate << "Process: GJets";
  stringstream DiBosonRate;
  DiBosonRate << "Process: DiBoson";
  stringstream TopRate;
  TopRate << "Process: TopRate";
  stringstream ZJetsRate;
  ZJetsRate << "Process: ZJetsRate";
  stringstream WJetsRate;
  WJetsRate << "Process: WJetsRate";
  stringstream ZnunuRate;
  ZnunuRate << "Process: ZnunuRate";
  stringstream PreRate;
  PreRate << "Process: Pre-fit (total)";
  stringstream PostRate;
  PostRate << "Process: Post-fit (total)";
  stringstream PostRateUnc;
  PostRateUnc << "Process: Post-fit uncertainty (total)";
  stringstream DataRate;
  DataRate << "Process: Data";

  for(int iBin = 0; iBin < qchist->GetNbinsX(); iBin++){
    QCDRate << "   ";
    QCDRate << qchist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < gmhist->GetNbinsX(); iBin++){
    GJetsRate << "   ";
    GJetsRate << gmhist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < dihist->GetNbinsX(); iBin++){
    DiBosonRate << "   ";
    DiBosonRate << dihist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < tthist->GetNbinsX(); iBin++){
    TopRate << "   ";
    TopRate << tthist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < zlhist->GetNbinsX(); iBin++){
    ZJetsRate << "   ";
    ZJetsRate << zlhist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < wlhist->GetNbinsX(); iBin++){
    WJetsRate << "   ";
    WJetsRate << wlhist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < znhist->GetNbinsX(); iBin++){
    ZnunuRate << "   ";
    ZnunuRate << znhist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < tphist->GetNbinsX(); iBin++){
    PreRate << "   ";
    PreRate << tphist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < tohist->GetNbinsX(); iBin++){
    PostRate << "   ";
    PostRate << tohist->GetBinContent(iBin+1);
    PostRateUnc << "   ";
    PostRateUnc << tohist->GetBinError(iBin+1);
  }

  for(int iBin = 0; iBin < dthist->GetNbinsX(); iBin++){
    DataRate << "   ";
    DataRate << dthist->GetBinContent(iBin+1);
  }

  outputfile<<"######################"<<endl;
  outputfile<<QCDRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<GJetsRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<DiBosonRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<TopRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<ZJetsRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<WJetsRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<ZnunuRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PreRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PostRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<PostRateUnc.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<DataRate.str()<<endl;
  outputfile<<"######################"<<endl;

  outputfile.close();

  // pad 1
  pad1->SetRightMargin(0.06);
  pad1->SetLeftMargin(0.12);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();
  
  //signal style  
  if(mjhist){
    mjhist->SetFillColor(0);
    mjhist->SetFillStyle(0);
    mjhist->SetLineColor(kBlack);
    mjhist->SetLineWidth(3);
    mjhist->Scale(scaleSig);
  }

  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kBlack);
    ggHhist->SetLineWidth(3);
    ggHhist->Scale(scaleSig);
  }

  if(mwhist){
    mwhist->SetFillColor(0);
    mwhist->SetFillStyle(0);
    mwhist->SetLineColor(kViolet);
    mwhist->SetLineWidth(3);
    mwhist->SetLineStyle(7);
    mwhist->Scale(scaleSig);
  }

  if(vbfhist){
    vbfhist->SetFillColor(0);
    vbfhist->SetFillStyle(0);
    vbfhist->SetLineColor(kViolet);
    vbfhist->SetLineWidth(3);
    vbfhist->SetLineStyle(7);
    vbfhist->Scale(scaleSig);
  }

  if(mzhist){
    mzhist->SetFillColor(0);
    mzhist->SetFillStyle(0);
    mzhist->SetLineColor(kOrange+1);
    mzhist->SetLineWidth(3);
    mzhist->SetLineStyle(7);
    mzhist->Scale(scaleSig);
  }

  if(wHhist){
    wHhist->SetFillColor(0);
    wHhist->SetFillStyle(0);
    wHhist->SetLineColor(kOrange+1);
    wHhist->SetLineWidth(3);
    wHhist->SetLineStyle(7);
    wHhist->Scale(scaleSig);
  }

  if(zHhist){
    zHhist->SetFillColor(0);
    zHhist->SetFillStyle(0);
    zHhist->SetLineColor(kGreen);
    zHhist->SetLineWidth(3);
    zHhist->SetLineStyle(7);
    zHhist->Scale(scaleSig);
  }
  
  znhist->SetFillColor(kGreen+1);
  znhist->SetLineColor(kBlack);

  gmhist->SetFillColor(13);
  gmhist->SetLineColor(13);

  zlhist->SetFillColor(13);
  zlhist->SetLineColor(kBlack);

  wlhist->SetFillColor(kRed);
  wlhist->SetLineColor(kBlack);

  tthist->SetFillColor(kBlue);
  tthist->SetLineColor(kBlack);

  dihist->SetFillColor(13);
  dihist->SetLineColor(13);

  qchist->SetFillColor(kGray);
  qchist->SetLineColor(kBlack);

  if(sighist){
    sighist->SetFillColor(kBlack);
    sighist->SetLineColor(kBlack);
    sighist->SetFillStyle(3001);
  }

  // make the stack
  THStack* stack = new THStack("stack", "stack");
  stack->Add(gmhist);
  stack->Add(dihist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(qchist);
  stack->Add(wlhist);
  stack->Add(znhist);
  if(plotSBFit && sighist)
    stack->Add(sighist);


  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category <=1)
    frame->GetYaxis()->SetRangeUser(0.0005,1000000);
  else
    frame->GetYaxis()->SetRangeUser(0.002,9000);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame ->Draw();

  CMS_lumi(pad1,"2.30",false,true);

  stack ->Draw("HIST SAME");
  if(mjhist && !plotSBFit)
    mjhist->Draw("HIST SAME");
  if(mwhist && !plotSBFit)
    mwhist->Draw("HIST SAME");
  if(mzhist && !plotSBFit)
    mzhist->Draw("HIST SAME");

  if(ggHhist && !plotSBFit)
    ggHhist->Draw("HIST SAME");
  if(vbfhist && !plotSBFit)
    vbfhist->Draw("HIST SAME");
  if(wHhist && !plotSBFit)
    wHhist->Draw("HIST SAME");
  if(zHhist && !plotSBFit)
    zHhist->Draw("HIST SAME");

  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetLineWidth(1);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("PE SAME");
  
  TLegend* leg = new TLegend(0.38, 0.38, 0.93, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);  
  leg->SetNColumns(2);

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

  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
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
  
  // for post-fit pre-fit data/mc
  TH1* dphist = (TH1*)dthist->Clone("dahist");
  TH1* dahist = (TH1*)dthist->Clone("dahist");

  dphist->SetLineColor(kRed);
  dahist->SetLineColor(kBlue);
  dphist->SetMarkerColor(kRed);
  dphist->SetMarkerSize(1);
  dphist->SetMarkerStyle(20);
  dahist->SetMarkerColor(kBlue);
  dahist->SetMarkerSize(1);
  dahist->SetMarkerStyle(20);

  TH1* mphist = (TH1*)tphist->Clone("mphist");
  TH1* mchist = (TH1*)zlhist->Clone("mchist");
  TH1* unhist = (TH1*)zlhist->Clone("unhist");
  mchist->Add(wlhist);
  mchist->Add(gmhist);
  mchist->Add(tthist);
  mchist->Add(dihist);
  mchist->Add(qchist);
  mchist->Add(znhist);
  if(sighist && plotSBFit)
    mchist->Add(sighist);

  for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0);
  for (int i = 1; i <= mphist->GetNbinsX(); i++) mphist->SetBinError(i, 0);

  // ratio data/post-fit
  dahist->Divide(mchist);
  // ratio data/pre-fit
  dphist->Divide(mphist);

  // ratio post fit at 1 with uncertaitny
  tohist->Divide(mchist);
  tohist->SetLineColor(0);
  tohist->SetMarkerColor(0);
  tohist->SetMarkerSize(0);
  tohist->SetFillColor(kGray);

  dahist->SetMarkerSize(1);
  dphist->SetMarkerSize(1);
  dahist->SetMarkerStyle(20);
  dphist->SetMarkerStyle(20);

  dahist->SetStats(kFALSE);
  dphist->SetStats(kFALSE);

  // line at 1
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetLineWidth(2);
  unhist->SetFillColor(0);

  dahist->GetXaxis()->SetLabelOffset(999999);
  dahist->GetXaxis()->SetLabelSize(0);
  dahist->GetXaxis()->SetTitleOffset(999999);
  dahist->GetXaxis()->SetTitleSize(0);

  tohist->Draw("E2 SAME");
  unhist->Draw("SAME");
  dahist->Draw("PE1 SAME");
  if(!blind)
    dphist->Draw("PE1 SAME");
  
  pad2->RedrawAxis("G sameaxis");

  canvas->cd();
  pad1->cd();

  leg->AddEntry(dthist, "Data", "PEL");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted Total Mono-X Signal", "F");

  leg->AddEntry(znhist, "Z #rightarrow #nu#nu", "F");
  leg->AddEntry(wlhist, "W #rightarrow l#nu", "F");
  leg->AddEntry(qchist, "QCD", "F");
  leg->AddEntry(tthist, "Top Quark", "F");
  leg->AddEntry(zlhist, "Others", "F");

  if(mjhist && !plotSBFit)
    leg->AddEntry(mjhist, Form("Mono-J (V,2 TeV x%d)",scaleSig));

  if(mwhist && !plotSBFit)
    leg->AddEntry(mwhist, Form("Mono-W (V,2 TeV x%d)",scaleSig));

  if(mzhist && !plotSBFit)
    leg->AddEntry(mzhist, Form("Mono-Z (V,2 TeV x%d)",scaleSig));

  if(ggHhist && !plotSBFit)
    leg->AddEntry(ggHhist, "gg #rightarrow H (m_{H} = 125 GeV)");

  if(vbfhist && !plotSBFit)
    leg->AddEntry(vbfhist, "qq #rightarrow H (m_{H} = 125 GeV)");

  if(wHhist && !plotSBFit)
    leg->AddEntry(wHhist, "qq #rightarrow WH (m_{H} = 125 GeV)");

  if(zHhist && !plotSBFit)
    leg->AddEntry(zHhist, "qq #rightarrow ZH (m_{H} = 125 GeV)");

  leg->AddEntry(dahist,"Expected (post-fit)","PEL");
  leg->AddEntry(dphist,"Expected (pre-fit)","PEL");


  leg->Draw("SAME");
  
  pad1->RedrawAxis("sameaxis");
  pad1->SetLogy();


  if(blind){
    canvas->SaveAs("postfit_sig_blind.pdf");
    canvas->SaveAs("postfit_sig_blind.png");
  }
  else{
    canvas->SaveAs("postfit_sig.pdf");
    canvas->SaveAs("postfit_sig.png");
  }
}

