#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void prepostSig(string fitFilename, 
		string templateFileName, 
		string observable, 
		Category category, 
		bool   isHiggsInvisible, 
		int    scaleSig = 1, 
		bool   blind = true, 
		bool   plotSBFit = false, 
		string interaction = "Vector", string mediatorMass = "2000", string DMMass = "10") {
  

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
  pad2->SetFillStyle(0);
  pad2->SetLineColor(0);
  pad2->SetGridy();

  TColor *color; // for color definition with alpha                                                                                                                             
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

  if(not blind){
    if(!isHiggsInvisible){
      mjhist = (TH1*) dfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      mwhist = (TH1*) dfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      mzhist = (TH1*) dfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());  
      if(mjhist)
	mjhist->Scale(1.0, "width");
      if(mwhist)
	mwhist->Scale(1.0, "width");
      if(mzhist)
	mzhist->Scale(1.0, "width");
    }
    else{
      ggHhist = (TH1*) dfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
      vbfhist = (TH1*) dfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
      wHhist  = (TH1*) dfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str());
      zHhist  = (TH1*) dfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str());
      if(ggHhist)
	ggHhist->Scale(1.0, "width");
      if(vbfhist)
	vbfhist->Scale(1.0, "width");
      if(wHhist)
	wHhist->Scale(1.0, "width");
      if(zHhist)
	zHhist->Scale(1.0, "width");
    }
  }

  TH1* znhist = NULL;
  TH1* zlhist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* qchist = NULL;
  TH1* gmhist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* tohist = NULL;
  TH1* tphist = NULL;
  TH1* sighist = NULL;

  if(!plotSBFit){

    znhist = (TH1*)pfile->Get("shapes_fit_b/ch1/Znunu");    
    zlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/ZJets");    
    wlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/WJets");    
    tthist = (TH1*)pfile->Get("shapes_fit_b/ch1/Top");    
    dihist = (TH1*)pfile->Get("shapes_fit_b/ch1/Dibosons");    
    ewkwhist = (TH1*)pfile->Get("shapes_fit_b/ch1/EWKW");    
    ewkzhist = (TH1*)pfile->Get("shapes_fit_b/ch1/EWKZ");    
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
    ewkwhist = (TH1*)pfile->Get("shapes_fit_s/ch1/EWKW");    
    ewkzhist = (TH1*)pfile->Get("shapes_fit_s/ch1/EWKZ");    
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
    dthist = (TH1*) tohist->Clone(("datahist_"+observable).c_str());
    for (int i = 0; i <= dthist->GetNbinsX(); i++) {
      dthist->SetBinError(i, 0.);      
    }
  }

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
  stringstream EWKWRate;
  EWKWRate << "Process: EWKWRate";
  stringstream EWKZRate;
  EWKZRate << "Process: EWKZRate";
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

  for(int iBin = 0; iBin < ewkwhist->GetNbinsX(); iBin++){
    EWKWRate << "   ";
    EWKWRate << ewkwhist->GetBinContent(iBin+1);
  }

  for(int iBin = 0; iBin < ewkzhist->GetNbinsX(); iBin++){
    EWKZRate << "   ";
    EWKZRate << ewkzhist->GetBinContent(iBin+1);
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
  outputfile<<EWKWRate.str()<<endl;
  outputfile<<"######################"<<endl;
  outputfile<<EWKZRate.str()<<endl;
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

  //signal style  
  if(mjhist){
    mjhist->SetFillColor(0);
    mjhist->SetFillStyle(0);
    mjhist->SetLineColor(kBlack);
    mjhist->SetLineStyle(7);
    mjhist->SetLineWidth(3);
    mjhist->Scale(scaleSig);
    mjhist->SetMarkerSize(0);
  }
  
  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kBlack);
    ggHhist->SetLineStyle(7);
    ggHhist->SetLineWidth(3);
    ggHhist->Scale(scaleSig);
    ggHhist->SetMarkerSize(0);
  }
    
  if(mwhist){
    mwhist->SetFillColor(0);
    mwhist->SetFillStyle(0);
    mwhist->SetLineColor(kBlue);
    mwhist->SetLineWidth(3);
    mwhist->Scale(scaleSig);
    mwhist->SetMarkerSize(0);
  }
  
  if(vbfhist){
    vbfhist->SetFillColor(0);
    vbfhist->SetFillStyle(0);
    vbfhist->SetLineColor(kBlue);
    vbfhist->SetLineWidth(3);
    //    vbfhist->SetLineStyle(7);
    vbfhist->Scale(scaleSig);
    vbfhist->SetMarkerSize(0);
  }

  if(mzhist){
    mzhist->SetFillColor(0);
    mzhist->SetFillStyle(0);
    mzhist->SetLineColor(TColor::GetColor("#A2C523"));
    mzhist->SetLineWidth(3);
    mzhist->Scale(scaleSig);
    mzhist->SetMarkerSize(0);
  }

  if(wHhist){
    wHhist->SetFillColor(0);
    wHhist->SetFillStyle(0);
    wHhist->SetLineColor(TColor::GetColor("#A2C523"));
    wHhist->SetLineWidth(3);
    wHhist->Scale(scaleSig);
    wHhist->SetMarkerSize(0);
  }

  if(zHhist){
    zHhist->SetFillColor(0);
    zHhist->SetFillStyle(0);
    zHhist->SetLineColor(TColor::GetColor("#A2C523"));
    zHhist->SetLineWidth(3);
    zHhist->Scale(scaleSig);
    zHhist->SetMarkerSize(0);
  }

  if(wHhist and zHhist)
    wHhist->Add(zHhist);
  
  qchist->SetFillColor(TColor::GetColor("#F1F1F2"));
  qchist->SetLineColor(kBlack);

  gmhist->SetFillColor(TColor::GetColor("#9A9EAB"));
  gmhist->SetLineColor(TColor::GetColor("#9A9EAB"));
  zlhist->SetFillColor(TColor::GetColor("#9A9EAB"));  
  zlhist->SetLineColor(kBlack);
  zlhist->Add(gmhist);

  znhist->SetFillColor(TColor::GetColor("#258039"));
  znhist->SetLineColor(kBlack);

  wlhist->SetFillColor(TColor::GetColor("#FAAF08"));
  wlhist->SetLineColor(kBlack);

  dihist->SetFillColor(TColor::GetColor("#4897D8"));
  dihist->SetLineColor(kBlack);

  tthist->SetFillColor(TColor::GetColor("#CF3721"));
  tthist->SetLineColor(kBlack);

  ewkwhist->SetFillColor(kBlack);
  ewkwhist->SetLineColor(kBlack);
  ewkzhist->SetFillColor(kBlack);
  ewkzhist->SetLineColor(kBlack);

  if(sighist){
    sighist->SetFillColor(kBlack);
    sighist->SetLineColor(kBlack);
    sighist->SetFillStyle(3001);
  }

  // make the stack
  THStack* stack = new THStack("stack", "stack");
  stack->Add(qchist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(dihist);
  //  ewkzhist->Add(ewkwhist);
  //  stack->Add(ewkzhist);
  stack->Add(wlhist);
  stack->Add(znhist);
  if(plotSBFit && sighist)
    stack->Add(sighist);


  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*100);
  else
    frame->GetYaxis()->SetRangeUser(0.0005,wlhist->GetMaximum()*200);

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

  CMS_lumi(canvas,"2.6");

  stack ->Draw("HIST SAME");
  if(mwhist && !plotSBFit)
    mwhist->Draw("HIST SAME");
  if(mzhist && !plotSBFit)
    mzhist->Draw("HIST SAME");
  if(mjhist && !plotSBFit)
    mjhist->Draw("HIST SAME");

  if(vbfhist && !plotSBFit)
    vbfhist->Draw("HIST SAME");
  if(wHhist && !plotSBFit)
    wHhist->Draw("HIST SAME");
  if(ggHhist && !plotSBFit)
    ggHhist->Draw("HIST SAME");

  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetLineWidth(1);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("PE SAME");
  
  TLegend* leg = new TLegend(0.6, 0.55, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);  

  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  TH1* frame2 = (TH1*) dthist->Clone("frame2");
  frame2->Reset();
  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else
    frame2->GetYaxis()->SetRangeUser(-0.5,2);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(210);
  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->Draw("AXIS");
  frame2->Draw("AXIG same");

  // for post-fit pre-fit data/mc
  TH1* dphist = (TH1*)dthist->Clone("dahist");
  TH1* dahist = (TH1*)dthist->Clone("dahist");

  dphist->SetLineColor(kRed);
  dphist->SetMarkerColor(kRed);

  dahist->SetLineColor(kBlue);
  dahist->SetMarkerColor(kBlue);

  dphist->SetMarkerSize(1);
  dphist->SetMarkerStyle(20);
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
  TH1* htemp = (TH1*) tohist->Clone("postfit_over_prefit");
  if(plotSBFit && sighist)
    tohist->Add(sighist);

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
  if(!blind)
    dphist->Draw("PE1 SAME");
  dahist->Draw("PE1 SAME");

  TLegend* leg2 = new TLegend(0.14,0.14,0.65,0.21,NULL,"brNDC");
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(dahist,"Backgraound (Post-Fit)","PLE");
  leg2->AddEntry(dphist,"Backgraound (Pre-Fit)","PLE");
  leg2->Draw("same");

  
  pad2->RedrawAxis("G sameaxis");

  canvas->cd();
  leg->AddEntry(dthist, "Data", "PEL");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted Total Mono-X Signal", "F");

  leg->AddEntry(znhist, "Z #rightarrow #nu#nu", "F");
  leg->AddEntry(wlhist, "W #rightarrow l#nu", "F");
  leg->AddEntry(dihist, "WW/WZ/ZZ", "F");
  //  leg->AddEntry(ewkzhist, "EWK W/Z+2jets", "F");
  leg->AddEntry(tthist, "Top Quark", "F");
  leg->AddEntry(tthist, "Top Quark", "F");
  leg->AddEntry(zlhist, "Z/#gamma #rightarrow ll, #gamma+jets", "F");
  leg->AddEntry(qchist, "QCD", "F");

  if(mjhist && !plotSBFit)
    leg->AddEntry(mjhist, Form("Mono-J (V,2 TeV x%d)",scaleSig),"L");

  if(mwhist && !plotSBFit)
    leg->AddEntry(mwhist, Form("Mono-W (V,2 TeV x%d)",scaleSig),"L");

  if(mzhist && !plotSBFit)
    leg->AddEntry(mzhist, Form("Mono-Z (V,2 TeV x%d)",scaleSig),"L");

  if(ggHhist && !plotSBFit)
    leg->AddEntry(ggHhist, "ggH m_{H} = 125 GeV","L");

  if(vbfhist && !plotSBFit)
    leg->AddEntry(vbfhist, "qqH m_{H} = 125 GeV","L");

  if(wHhist && !plotSBFit)
    leg->AddEntry(wHhist, "VH m_{H} = 125 GeV","L");


  leg->Draw("SAME");  
  pad2->RedrawAxis("sameaxis");
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();

  if(blind){
    canvas->SaveAs("postfit_sig_blind.pdf");
    canvas->SaveAs("postfit_sig_blind.png");
  }
  else{
    canvas->SaveAs("postfit_sig.pdf");
    canvas->SaveAs("postfit_sig.png");
  }
  
  TH1* totalSignal = NULL;

  if(isHiggsInvisible and not plotSBFit){
    totalSignal = (TH1*) ggHhist->Clone("totalSignal");
    totalSignal->Add(vbfhist);
    totalSignal->Add(wHhist);
    totalSignal->Add(zHhist);
  }
  else if(not isHiggsInvisible and not plotSBFit){
    totalSignal = (TH1*) mjhist->Clone("totalSignal");
    totalSignal->Add(mwhist);
    totalSignal->Add(mzhist);
  }
  else if(plotSBFit)
    totalSignal = (TH1*) sighist->Clone("totalSignal");

  canvas->cd();
  pad2->Draw();
  pad2->cd();
  if(not plotSBFit)
    frame2->GetYaxis()->SetTitle("(S+B)/B");
  else
    frame2->GetYaxis()->SetTitle("(S_{fit}+B)/B");

  TH1* SoverB_prefit = (TH1*) totalSignal->Clone("SoverB_prefit");
  TH1* SoverB_postfit = (TH1*) totalSignal->Clone("SoverB_postfit");
  SoverB_prefit->SetLineColor(kRed);
  SoverB_prefit->SetMarkerColor(kRed);
  SoverB_prefit->SetMarkerSize(1);
  SoverB_prefit->SetMarkerStyle(20);
  SoverB_postfit->SetLineColor(kBlue);
  SoverB_postfit->SetMarkerColor(kBlue);
  SoverB_postfit->SetMarkerSize(1);
  SoverB_postfit->SetMarkerStyle(20);

  SoverB_prefit->Add(tphist);
  SoverB_prefit->Divide(tphist);
  SoverB_postfit->Add(htemp);
  SoverB_postfit->Divide(htemp);

  frame2->GetYaxis()->SetRangeUser(0.5,SoverB_postfit->GetMaximum()*1.2);
  frame2->Draw();

  SoverB_postfit->Draw("hist same");
  TH1* SoverB_postfit_d = (TH1*) SoverB_postfit->Clone("SoverB_postfit_d");
  for(int iBin = 0; iBin < SoverB_postfit_d->GetNbinsX(); iBin++)
    SoverB_postfit_d->SetBinContent(iBin+1,1);
  SoverB_postfit_d->SetLineColor(0);
  
  SoverB_postfit->Draw("hist same");
  SoverB_postfit_d->SetMarkerColor(0);
  SoverB_postfit_d->SetMarkerSize(0);
  SoverB_postfit_d->SetFillColor(kGray);
  SoverB_postfit_d->SetFillStyle(1001);
  SoverB_postfit_d->Draw("E2 SAME");
  unhist->Draw("SAME");
  SoverB_prefit->Draw("hist same");
  SoverB_postfit->Draw("hist same");
  pad2->RedrawAxis("sameaxis");
  
  canvas->SaveAs("postfit_sig_SoB.pdf");
  canvas->SaveAs("postfit_sig_SoB.png");
 
  TFile* outFile = new TFile("postfit_weights_Sig.root","RECREATE");
  outFile->cd();
  htemp->Divide(tphist);
  htemp->Write("postfit_over_prefit");
  outFile->Close();
}

