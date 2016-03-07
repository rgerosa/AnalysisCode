#include "CMS_lumi.h"

void prepostSig_fromScan(string fitFilename, 
			 string dataFile,
			 string templateFileMonoJ, 
			 string templateFileMonoW, 
			 string templateFileMonoZ, 
			 string observable, 
			 int category, 
			 string interaction,
			 string mediatorMass = "1000", 
			 string DMMass = "0050", 
			 bool blind = true, bool plotSBFit = false) {
  

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

  TFile* dfileMonoJ = TFile::Open(("root://eoscms.cern.ch//"+templateFileMonoJ).c_str());
  TFile* dfileMonoW = TFile::Open(("root://eoscms.cern.ch//"+templateFileMonoW).c_str());
  TFile* dfileMonoZ = TFile::Open(("root://eoscms.cern.ch//"+templateFileMonoZ).c_str());
  TFile* dfile = TFile::Open((dataFile).c_str());
  // in case of b-only fit just dispaly three possible signal on the stack
  string nameInteraction;
  if(interaction == "Vector")
    nameInteraction = "800";
  else if(interaction == "Axial")
    nameInteraction = "801";
  else if(interaction == "Scalar")
    nameInteraction = "805";
  else if(interaction == "Pseudoscalar")
    nameInteraction = "806";
  
  TH1* mjhist = (TH1*) dfileMonoJ->Get(("category_monojet/signal_signal_"+nameInteraction+mediatorMass+DMMass).c_str());
  TH1* mwhist = (TH1*) dfileMonoW->Get(("category_monojet/signal_signal_"+nameInteraction+mediatorMass+DMMass).c_str());
  TH1* mzhist = (TH1*) dfileMonoZ->Get(("category_monojet/signal_signal_"+nameInteraction+mediatorMass+DMMass).c_str());
  mjhist->Scale(1.0, "width");
  mwhist->Scale(1.0, "width");
  mzhist->Scale(1.0, "width");
  
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

  if(! plotSBFit){

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
    dthist = (TH1*)dfile->Get(("datahist_"+observable).c_str());
    dthist->Scale(1.0,"width");
  }
  else{
    dthist = (TH1*)dfile->Get(("datahist_"+observable).c_str());
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
  outputfile<<DataRate.str()<<endl;
  outputfile<<"######################"<<endl;

  outputfile.close();

  // pad 1
  pad1->SetRightMargin(0.075);
  pad1->SetLeftMargin(0.10);
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
  }

  if(mwhist){
    mwhist->SetFillColor(0);
    mwhist->SetFillStyle(0);
    mwhist->SetLineColor(kViolet);
    mwhist->SetLineWidth(3);
    mwhist->SetLineStyle(7);
  }

  if(mzhist){
    mzhist->SetFillColor(0);
    mzhist->SetFillStyle(0);
    mzhist->SetLineColor(kOrange+1);
    mzhist->SetLineWidth(3);
    mzhist->SetLineStyle(7);
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


  TH1* frame = NULL;
  if(category <=1)
    frame = pad1->DrawFrame(200., 0.0005, 1250., 100000, "");
  else
    frame = pad1->DrawFrame(250., 0.0002, 1000., 100000, "");

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame ->Draw();

  CMS_lumi(pad1, 4, 0,true);

  stack ->Draw("HIST SAME");
  if(mjhist && !plotSBFit)
    mjhist->Draw("HIST SAME");
  if(mwhist && !plotSBFit)
    mwhist->Draw("HIST SAME");
  if(mzhist && !plotSBFit)
    mzhist->Draw("HIST SAME");

  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetMarkerSize(1);
  dthist->SetLineWidth(2);
  dthist->SetMarkerStyle(20);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("P SAME");
  
  TLegend* leg = new TLegend(0.53, 0.32, 0.9, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  

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
    frame2 = pad2->DrawFrame(200., 0.5, 1250., 1.5, "");
  else
    frame2 = pad2->DrawFrame(250., 0.5, 1000., 1.5, "");

  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->Draw();
  
  // for post-fit pre-fit data/mc
  TH1* dphist = (TH1*)dthist->Clone("dahist");
  TH1* dahist = (TH1*)dthist->Clone("dahist");
  dphist->SetLineColor(kRed);
  dphist->SetMarkerColor(kRed);
  dphist->SetMarkerSize(0.7);
  dahist->SetLineColor(kBlue);
  dahist->SetMarkerColor(kBlue);
  dahist->SetMarkerSize(0.7);

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
  unhist->SetFillColor(0);

  dahist->GetXaxis()->SetLabelOffset(999999);
  dahist->GetXaxis()->SetLabelSize(0);
  dahist->GetXaxis()->SetTitleOffset(999999);
  dahist->GetXaxis()->SetTitleSize(0);

  tohist->Draw("E2 SAME");
  unhist->Draw("SAME");
  dahist->Draw("PE SAME");
  if(!blind)
    dphist->Draw("PE SAME");
  
  pad2->RedrawAxis("sameaxis");

  canvas->cd();
  pad1->cd();

  leg->AddEntry(dthist, "Data");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted Total Mono-X Signal", "F");
  leg->AddEntry(znhist, "Z(#nu#nu)", "F");
  leg->AddEntry(wlhist, "W(l#nu)", "F");
  leg->AddEntry(qchist, "QCD", "F");
  leg->AddEntry(tthist, "Top", "F");
  leg->AddEntry(zlhist, "Others", "F");
  if(mjhist && !plotSBFit)
    leg->AddEntry(mjhist, "Mono-J (S,150 GeV)");
  if(mwhist && !plotSBFit)
    leg->AddEntry(mwhist, "Mono-W (S,150 GeV)");
  if(mzhist && !plotSBFit)
    leg->AddEntry(mzhist, "Mono-Z (S,150 GeV)");
  leg->AddEntry(dphist,"Expected (post-fit)");
  leg->AddEntry(dahist,"Expected (pre-fit)");

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

