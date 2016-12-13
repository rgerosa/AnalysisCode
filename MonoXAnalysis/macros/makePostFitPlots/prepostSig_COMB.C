#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static bool saveTextFile = false;
static bool plot_significance = false;

void prepostSig_COMB(string fitFilename, 
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

  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF)
    postfix = "_VBF";

  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  string dir = "ch1_ch1";
  if(category == Category::monoV)
    dir = "ch2_ch1";
  else if(category == Category::VBF)
    dir = "ch3_ch1";


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
  
  // take signals used for the fit
  if(!isHiggsInvisible){
    mjhist = (TH1*) pfile->Get("shapes_prefit/ch1/MonoJ");
    mwhist = (TH1*) pfile->Get("shapes_prefit/ch1/MonoW");
    mzhist = (TH1*) pfile->Get("shapes_prefit/ch1/MonoZ");
  }
  else{
    ggHhist  = (TH1*) pfile->Get("shapes_prefit/ch1/ggH");
    vbfhist = (TH1*) pfile->Get("shapes_prefit/ch1/qqH");
    wHhist   = (TH1*) pfile->Get("shapes_prefit/ch1/WH");
    zHhist   = (TH1*) pfile->Get("shapes_prefit/ch1/ZH");
    ggZHhist = (TH1*) pfile->Get("shapes_prefit/ch1/ggZH");
  }

  // take backgrouds
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
    
    znhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Znunu").c_str());    
    zlhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/ZJets").c_str());    
    wlhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/WJets").c_str());    
    tthist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Top").c_str());    
    dihist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Dibosons").c_str());    
    ewkwhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/WJets_EWK").c_str());    
    ewkzhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/Znunu_EWK").c_str());    
    qchist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/QCD").c_str());    
    gmhist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/GJets").c_str());    
    tohist = (TH1*)pfile->Get(("shapes_fit_b/"+dir+"/total_background").c_str());    
    tphist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());    
  }
  else{
    znhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/Znunu").c_str());    
    zlhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/ZJets").c_str());    
    wlhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/WJets").c_str());    
    tthist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/Top").c_str());    
    dihist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/Dibosons").c_str());    
    ewkwhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/WJets_EWK").c_str());    
    ewkzhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/Znunu_EWK").c_str());    
    qchist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/QCD").c_str());    
    gmhist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/GJets").c_str());    
    tohist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/total_background").c_str());    
    tphist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());      
    sighist = (TH1*)pfile->Get(("shapes_fit_s/"+dir+"/total_signal").c_str());
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

  if(saveTextFile){
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
      QCDRate << qchist->GetBinContent(iBin+1)*qchist->GetBinWidth(iBin+1) << " \\pm "<<qchist->GetBinError(iBin+1)*qchist->GetBinWidth(iBin+1);
    }

    for(int iBin = 0; iBin < gmhist->GetNbinsX(); iBin++){
      GJetsRate << "   ";
      GJetsRate << gmhist->GetBinContent(iBin+1)*gmhist->GetBinWidth(iBin+1) << " \\pm "<<gmhist->GetBinError(iBin+1)*gmhist->GetBinWidth(iBin+1);
    }

    for(int iBin = 0; iBin < dihist->GetNbinsX(); iBin++){
      DiBosonRate << "   ";
      DiBosonRate << dihist->GetBinContent(iBin+1)*dihist->GetBinWidth(iBin+1) << " \\pm "<<dihist->GetBinError(iBin+1)*dihist->GetBinWidth(iBin+1);
    }
    
    for(int iBin = 0; iBin < tthist->GetNbinsX(); iBin++){
      TopRate << "   ";
      TopRate << tthist->GetBinContent(iBin+1)*tthist->GetBinWidth(iBin+1) << " \\pm "<<tthist->GetBinError(iBin+1)*tthist->GetBinWidth(iBin+1);
    }
    
    if(category == Category::VBF){
      for(int iBin = 0; iBin < ewkwhist->GetNbinsX(); iBin++){
	EWKWRate << "   ";
	EWKWRate << ewkwhist->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < ewkzhist->GetNbinsX(); iBin++){
	EWKZRate << "   ";
	EWKZRate << ewkzhist->GetBinContent(iBin+1);
      }
    }
    
    for(int iBin = 0; iBin < zlhist->GetNbinsX(); iBin++){
      ZJetsRate << "   ";
      ZJetsRate << zlhist->GetBinContent(iBin+1)*zlhist->GetBinWidth(iBin+1) << " \\pm "<<zlhist->GetBinError(iBin+1)*zlhist->GetBinWidth(iBin+1);
    }
    
    for(int iBin = 0; iBin < wlhist->GetNbinsX(); iBin++){
      WJetsRate << "   ";
      WJetsRate << wlhist->GetBinContent(iBin+1)*wlhist->GetBinWidth(iBin+1) << " \\pm "<<wlhist->GetBinError(iBin+1)*wlhist->GetBinWidth(iBin+1);
    }
    
    for(int iBin = 0; iBin < znhist->GetNbinsX(); iBin++){
      ZnunuRate << "   ";
      ZnunuRate << znhist->GetBinContent(iBin+1)*znhist->GetBinWidth(iBin+1) << " \\pm "<<znhist->GetBinError(iBin+1)*znhist->GetBinWidth(iBin+1);
    }
    
    for(int iBin = 0; iBin < tphist->GetNbinsX(); iBin++){
      PreRate << "   ";
      PreRate << tphist->GetBinContent(iBin+1)*tphist->GetBinWidth(iBin+1) << " \\pm "<<tphist->GetBinError(iBin+1)*tphist->GetBinWidth(iBin+1);
    }

    for(int iBin = 0; iBin < tohist->GetNbinsX(); iBin++){
      PostRate << "   ";
      PostRate << tohist->GetBinContent(iBin+1)*tohist->GetBinWidth(iBin+1) << " \\pm "<<tohist->GetBinError(iBin+1)*tohist->GetBinWidth(iBin+1);
    }
    
    for(int iBin = 0; iBin < dthist->GetNbinsX(); iBin++){
      DataRate << "   ";
      DataRate << dthist->GetBinContent(iBin+1)*dthist->GetBinWidth(iBin+1) << " \\pm "<<dthist->GetBinError(iBin+1)*dthist->GetBinWidth(iBin+1);
    }
    
    outputfile<<"######################"<<endl;
    outputfile<<QCDRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<GJetsRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<DiBosonRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<TopRate.str()<<endl;
    if(category == Category::VBF){
      outputfile<<"######################"<<endl;
      outputfile<<EWKWRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<EWKZRate.str()<<endl;
    }
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
  }


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

  if(qchist){
    qchist->SetFillColor(TColor::GetColor("#F1F1F2"));
    qchist->SetLineColor(kBlack);
  }

  gmhist->SetFillColor(TColor::GetColor("#9A9EAB"));
  gmhist->SetLineColor(TColor::GetColor("#9A9EAB"));
  zlhist->SetFillColor(TColor::GetColor("#9A9EAB"));  
  zlhist->SetLineColor(kBlack);
  zlhist->Add(gmhist);

  znhist->SetFillColor(TColor::GetColor("#3A8C4C"));
  znhist->SetLineColor(kBlack);

  wlhist->SetFillColor(TColor::GetColor("#FAAF08"));
  wlhist->SetLineColor(kBlack);

  if(category != Category::VBF)
    dihist->SetFillColor(TColor::GetColor("#4897D8"));
  else
    dihist->SetFillColor(kRed+3);
  dihist->SetLineColor(kBlack);

  tthist->SetFillColor(TColor::GetColor("#CF3721"));
  tthist->SetLineColor(kBlack);

  if(category == Category::VBF){
    ewkzhist->SetFillColor(kCyan+1);
    ewkzhist->SetLineColor(kBlack);
    ewkwhist->SetFillColor(kAzure+1);
    ewkwhist->SetLineColor(kBlack);
  }

  if(sighist){
    sighist->SetFillColor(kBlack);
    sighist->SetLineColor(kBlack);
    sighist->SetFillStyle(3001);
  }
  
  // make the stack
  THStack* stack = new THStack("stack", "stack");
  if(qchist)
    stack->Add(qchist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(dihist);
  if(category == Category::VBF)
    stack->Add(ewkwhist);
  if(category == Category::VBF)
    stack->Add(ewkzhist);
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


  CMS_lumi(canvas,"36.2");

  TLatex* categoryLabel = new TLatex();
  categoryLabel->SetNDC();
  categoryLabel->SetTextSize(0.5*canvas->GetTopMargin());
  categoryLabel->SetTextFont(42);
  categoryLabel->SetTextAlign(11);
  if(category == Category::monojet)
    categoryLabel ->DrawLatex(0.175,0.80,"monojet");
  else if(category == Category::monoV)
    categoryLabel ->DrawLatex(0.175,0.80,"mono-V");
  else if(category == Category::VBF)
    categoryLabel ->DrawLatex(0.175,0.80,"VBF");
  categoryLabel->Draw("same");

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
  
  TLegend* leg = new TLegend(0.6, 0.55, 0.92, 0.92);
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
    frame2->GetYaxis()->SetRangeUser(0.4,1.6);
  else
    frame2->GetYaxis()->SetRangeUser(0.4,1.6);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(510);

  if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi"))
    frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  else if(TString(observable).Contains("mjj"))
    frame2->GetXaxis()->SetTitle("m_{jj} [GeV]");
  else if(TString(observable).Contains("detajj"))
    frame2->GetXaxis()->SetTitle("#Delta#eta_{jj}");

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
  TH1* mchist = (TH1*)tphist->Clone("mchist");
  TH1* unhist = (TH1*)tphist->Clone("unhist");
  mchist->Reset();
  unhist->Reset();
  if(qchist)
    mchist->Add(qchist);
  mchist->Add(wlhist);
  mchist->Add(zlhist);
  mchist->Add(tthist);
  mchist->Add(dihist);
  mchist->Add(znhist);
  if(category == Category::VBF){
    mchist->Add(ewkwhist);
    mchist->Add(ewkzhist);
  }

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
    dphist->Draw("P0E1 SAME");
  dahist->Draw("P0E1 SAME");

  TLegend* leg2 = new TLegend(0.14,0.24,0.40,0.28,NULL,"brNDC");
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(dphist,"post-fit","PLE");
  leg2->AddEntry(dahist,"pre-fit","PLE");
  leg2->Draw("same");

  
  pad2->RedrawAxis("G sameaxis");

  canvas->cd();
  leg->AddEntry(dthist, "Data", "PEL");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted Total Mono-X Signal", "F");

  leg->AddEntry(znhist, "Z #rightarrow #nu#nu", "F");
  leg->AddEntry(wlhist, "W #rightarrow l#nu", "F");
  if(category == Category::VBF){
    leg->AddEntry(ewkzhist,"Z-EWK #rightarrow #nu#nu", "F");
    leg->AddEntry(ewkwhist,"W-EWK #rightarrow l#nu", "F");
  } 
  leg->AddEntry(dihist, "WW/WZ/ZZ", "F");
  leg->AddEntry(tthist, "Top Quark", "F");
  leg->AddEntry(zlhist, "Z/#gamma #rightarrow ll, #gamma+jets", "F");
  if(qchist)
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
    canvas->SaveAs(("postfit_sig_blind"+postfix+".pdf").c_str());
    canvas->SaveAs(("postfit_sig_blind"+postfix+".png").c_str());
  }
  else{
    canvas->SaveAs(("postfit_sig"+postfix+".pdf").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+".png").c_str());
  }

  if(plot_significance){
    
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
    
    canvas->SaveAs(("postfit_sig_SoB"+postfix+".pdf").c_str());
    canvas->SaveAs(("postfit_sig_SoB"+postfix+".png").c_str());
    
    TFile* outFile = new TFile("postfit_weights_Sig.root","RECREATE");
    outFile->cd();
    htemp->Divide(tphist);
    htemp->Write("postfit_over_prefit");
    outFile->Close();
  }
}

