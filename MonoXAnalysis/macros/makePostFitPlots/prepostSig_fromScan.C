#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static bool saveTextFile = true;
static bool dumpInfo     = false;
static bool plotSignificance = true;
static float lumiScale_Higgs = 15;
static float lumiScale_DM = 2.78;
static bool addStatUncPull = false;

void prepostSig_fromScan(string   fitFilename, 
			 string   observable, 
			 Category category, 
			 bool     isCombinedFit = false,
			 bool     plotSBFit   = false,
			 bool     addPullPlot = false,  
			 bool     addPreFitOnPull = false,
			 int      scaleSig = 1, 
			 bool     blind    = false){
  
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();


  TCanvas* canvas = NULL;
  TPad *pad2 = NULL;
  TPad *pad3 = NULL;

  if(not addPullPlot){

    canvas = new TCanvas("canvas", "canvas", 600, 700);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->cd();
    canvas->SetBottomMargin(0.3);
    canvas->SetRightMargin(0.06);

    pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
    pad2->SetTopMargin(0.7);
    pad2->SetRightMargin(0.06);
    pad2->SetFillColor(0);
    //    pad2->SetGridy(1);
    pad2->SetFillStyle(0);
  }
  else{

    canvas = new TCanvas("canvas", "canvas", 600, 800);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->cd();
    canvas->SetBottomMargin(0.38);
    canvas->SetRightMargin(0.06);

    pad2 = new TPad("pad2","pad2",0,0.,1,1.);
    pad2->SetTopMargin(0.63);
    pad2->SetBottomMargin(0.25);
    pad2->SetRightMargin(0.06);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    pad2->SetLineColor(0);
    //    pad2->SetGridy();

    pad3 = new TPad("pad3","pad3",0,0.,1,1.);
    pad3->SetTopMargin(0.76);
    pad3->SetRightMargin(0.06);
    pad3->SetFillColor(0);
    pad3->SetFillStyle(0);
    pad3->SetLineColor(0);
    //    pad3->SetGridy();
  }

  TColor *color; // for color definition with alpha                                                                                                                             
  TFile* pfile = new TFile(fitFilename.c_str());

  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  string dir;
  if(isCombinedFit){
    if(category == Category::monojet)
      dir = "ch1_ch1";
    else if(category == Category::monoV)
      dir = "ch2_ch1";
    else if(category == Category::VBF)
      dir = "ch3_ch1";
  }
  else if( category != Category::VBF)
    dir = "ch1";
  else
    dir = "ch1";

  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF)
    postfix = "_VBF";

  TFile*monoj_av = NULL, *monow_av = NULL, *monoz_av = NULL, *higgs = NULL;
  
  if(category == Category::monoV){
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonov_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonov_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonov_13TeV_v1.root","READ");
    higgs    = new TFile("~/work/MONOJET_ANALYSIS/CMSSW_7_4_16/src/AnalysisCode/MonoXAnalysis/macros/monoV_hinv_forCombination/templates_met_v2.root","READ");
  }
  else{
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");
    higgs    = new TFile("~/work/MONOJET_ANALYSIS/CMSSW_7_4_16/src/AnalysisCode/MonoXAnalysis/macros/monoj_hinv_forCombination/templates_met_v2.root","READ");
  }

  // in case of b-only fit just dispaly three possible signal on the stack
  TH1* mjhist_av = NULL;
  TH1* mwhist_av = NULL;
  TH1* mzhist_av = NULL;

  // in case of higgs invisible limits
  TH1* ggHhist = NULL;
  TH1* vbfhist = NULL;
  TH1* wHhist = NULL;
  TH1* zHhist = NULL;
  TH1* ggZHhist = NULL;  

  // signals for axial vector model 
  mjhist_av = (TH1*) monoj_av->FindObjectAny("signal_signal_80120000001");
  mwhist_av = (TH1*) monow_av->FindObjectAny("signal_signal_80120000001");
  mzhist_av = (TH1*) monoz_av->FindObjectAny("signal_signal_80120000001");
  mjhist_av->Scale(1.,"width");
  mwhist_av->Scale(1.,"width");
  mzhist_av->Scale(1.,"width");
  // summing all signals together
  mjhist_av->Add(mwhist_av);
  mjhist_av->Add(mzhist_av);
  mjhist_av->Scale(lumiScale_DM);

  // signals for higgs invisible
  ggHhist  = (TH1*) higgs->FindObjectAny("ggHhist_125_met");
  vbfhist  = (TH1*) higgs->FindObjectAny("vbfHhist_125_met");
  wHhist   = (TH1*) higgs->FindObjectAny("wHhist_125_met");
  zHhist   = (TH1*) higgs->FindObjectAny("zHhist_125_met");
  ggZHhist = (TH1*) higgs->FindObjectAny("ggzHhist_125_met");
  ///
  ggHhist->Scale(1.,"width");
  vbfhist->Scale(1.,"width");
  wHhist->Scale(1.,"width");
  zHhist->Scale(1.,"width");
  ggZHhist->Scale(1.,"width");
  ///
  ggHhist->Add(vbfhist);
  ggHhist->Add(wHhist);
  ggHhist->Add(zHhist);
  ggHhist->Add(ggZHhist);
  ggHhist->Scale(lumiScale_Higgs);

  // background
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

  znhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Znunu").c_str());    
  zlhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/ZJets").c_str());    
  wlhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/WJets").c_str());    
  tthist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Top").c_str());    
  dihist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Dibosons").c_str());    

  if(category == Category::VBF){
    ewkwhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/WJets_EWK").c_str());    
    ewkzhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Znunu_EWK").c_str());    
  }
  else{
    ewkwhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/WJets_EWK").c_str());    
    ewkzhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/ZJets_EWK").c_str());    
  }
  qchist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/QCD").c_str());    
  gmhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/GJets").c_str());    
  tohist = (TH1*) ((TH1*)pfile->Get((fit_dir+"/"+dir+"/total_background").c_str()))->Clone("tohist");    
  tphist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());    

  if(plotSBFit)
    sighist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/total_signal").c_str());

  ////////////////
  TGraphAsymmErrors* dthist = NULL;
  if(!blind)
    dthist = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir+"/data").c_str());
  else{

    dthist = new TGraphAsymmErrors();
    for(int iBin = 1; iBin < tohist->GetNbinsX()+1; iBin++){
      dthist->SetPoint(iBin,tohist->GetBinCenter(iBin),tohist->GetBinContent(iBin));
      dthist->SetPointError(iBin,tohist->GetBinWidth(iBin)/2,tohist->GetBinWidth(iBin)/2,tohist->GetBinError(iBin)/2,tohist->GetBinError(iBin)/2);
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
    
    for(int iBin = 0; iBin < dthist->GetN(); iBin++){
      double x,y;
      dthist->GetPoint(iBin,x,y);
      DataRate << "   ";
      DataRate << y*tohist->GetBinWidth(iBin+1);
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
  if(mjhist_av){
    mjhist_av->SetFillColor(0);
    mjhist_av->SetFillStyle(0);
    mjhist_av->SetLineColor(kBlue);
    mjhist_av->SetLineWidth(3);
    mjhist_av->Scale(scaleSig);
    mjhist_av->SetMarkerSize(0);
  }
  
  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kBlack);
    ggHhist->SetLineWidth(3);
    ggHhist->Scale(scaleSig);
    ggHhist->SetMarkerSize(0);
  }
    
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
  else{
    if(ewkwhist){
      ewkwhist->SetFillColor(kViolet+1);
      ewkwhist->SetLineColor(kBlack);
    }
    if(ewkzhist){
      ewkzhist->SetFillColor(kCyan+1);
      ewkzhist->SetLineColor(kBlack);
    }
  }

  
  if(sighist){
    sighist->SetFillColor(kBlack);
    sighist->SetLineColor(kBlack);
    sighist->SetLineWidth(3);
    sighist->SetFillColor(0);
    sighist->SetFillStyle(0);
  }


  // make the stack for backgrounds
  THStack* stack = new THStack("stack", "stack");
  if(qchist)
    stack->Add(qchist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(dihist);  
  if(ewkwhist)
    stack->Add(ewkwhist);
  if(ewkzhist)
    stack->Add(ewkzhist);
  stack->Add(wlhist);
  stack->Add(znhist);

  TH1* frame = (TH1*) tohist->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);

  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*500);
  else if(category == Category::monoV)
    frame->GetYaxis()->SetRangeUser(0.01,wlhist->GetMaximum()*500);
  else if(category == Category::VBF)
    frame->GetYaxis()->SetRangeUser(0.005,tphist->GetMaximum()*500);

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

  frame->Draw();

  CMS_lumi(canvas,"35.9");

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
  if(plotSBFit)
    sighist->Draw("HIST same");
  if(not plotSBFit and category != Category::VBF){
    mjhist_av->Draw("hist same");
    ggHhist->Draw("hist same"); 
  }

  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetLineWidth(1);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("PE SAME");

  TLegend* leg = new TLegend(0.50, 0.55, 0.92, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(dthist, "Data", "PEL");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted signal", "L");
  leg->AddEntry(znhist,  "Z(#nu#nu)+jets", "F");
  leg->AddEntry(wlhist,  "W(l#nu)+jets", "F");
  if(ewkzhist)  leg->AddEntry(ewkzhist,"Z(#nu#nu)+jets EWK", "F");
  if(ewkwhist)  leg->AddEntry(ewkwhist,"W(l#nu)+jets EWK", "F");
  leg->AddEntry(dihist,  "WW/WZ/ZZ", "F");
  leg->AddEntry(tthist,  "Top quark", "F");
  leg->AddEntry(zlhist,  "Z/#gamma(ll), #gamma+jets", "F");
  if(qchist)
    leg->AddEntry(qchist,  "QCD", "F");
  if(not plotSBFit and category != Category::VBF)
    leg->AddEntry(ggHhist,   "Higgs invisible, m_{H} = 125 GeV","L");
  if(not plotSBFit and category != Category::VBF)
    leg->AddEntry(mjhist_av, "Axial-vector, m_{med} = 2.0 TeV","L");

  leg->Draw("SAME");    
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();


  ////
  TH1* frame2 =  (TH1*) tohist->Clone("frame");
  frame2->Reset();
  frame2->SetLineColor(kBlack);
  frame2->SetLineWidth(1);

  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.90,1.10);
  else
    frame2->GetYaxis()->SetRangeUser(0.75,1.25);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(210);
  frame2->GetYaxis()->SetNdivisions(5);

  if(not addPullPlot){
    frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    if(category == Category::VBF and TString(observable).Contains("mjj"))
      frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");
    frame2->GetYaxis()->SetTitle("Data / Pred.");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetTitleOffset(1.5);
    frame2->GetYaxis()->SetLabelSize(0.04);
    frame2->GetYaxis()->SetTitleSize(0.04);
    frame2->GetXaxis()->SetLabelSize(0.04);
    frame2->GetXaxis()->SetTitleSize(0.05); 
    frame2->GetXaxis()->SetTitleOffset(1.1);
  }
  else{
    frame2->GetYaxis()->SetTitleOffset(1.9);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetXaxis()->SetLabelSize(0);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetYaxis()->SetTitle("Data / Pred.");
    frame2->GetYaxis()->CenterTitle();
  }
  frame2->GetXaxis()->SetTickLength(0.025);
  frame2->Draw();


  // for post-fit pre-fit data/mc
  TGraphAsymmErrors* dphist = (TGraphAsymmErrors*)dthist->Clone("dphist");
  TGraphAsymmErrors* dahist = (TGraphAsymmErrors*)dthist->Clone("dahist");
  
  dphist->SetLineColor(kRed);
  dphist->SetMarkerColor(kRed);

  dahist->SetLineColor(TColor::GetColor("#0066ff"));
  dahist->SetMarkerColor(TColor::GetColor("#0066ff"));

  dphist->SetMarkerSize(1);
  dphist->SetMarkerStyle(24);
  dahist->SetMarkerSize(1);
  dahist->SetMarkerStyle(20);

  TH1* mphist = (TH1*)tphist->Clone("mphist");
  TH1* mchist = (TH1*)tphist->Clone("mchist");
  TH1* unhist = (TH1*)tphist->Clone("unhist");
  mchist->Reset();
  unhist->Reset();
  mchist->Add(qchist);
  mchist->Add(wlhist);
  mchist->Add(zlhist);
  mchist->Add(tthist);
  mchist->Add(dihist);
  mchist->Add(znhist);
  if(ewkwhist) mchist->Add(ewkwhist);
  if(ewkzhist) mchist->Add(ewkzhist);

  //  if(plotSBFit and sighist){
  //    mchist->Add(sighist);
  //    tohist->Add(sighist);
  //  }

  for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0);
  for (int i = 1; i <= mphist->GetNbinsX(); i++) mphist->SetBinError(i, 0);

  for(int iPoint = 0; iPoint < dphist->GetN(); iPoint++){
    double x,y;
    dphist->GetPoint(iPoint,x,y);
    dphist->SetPoint(iPoint,x,y/mphist->GetBinContent(iPoint+1));
    dphist->SetPointError(iPoint,dphist->GetErrorXlow(iPoint),dphist->GetErrorXhigh(iPoint),
                          dphist->GetErrorYlow(iPoint)/mphist->GetBinContent(iPoint+1),dphist->GetErrorYhigh(iPoint)/mphist->GetBinContent(iPoint+1));
    dahist->GetPoint(iPoint,x,y);
    dahist->SetPoint(iPoint,x,y/mchist->GetBinContent(iPoint+1));
    dahist->SetPointError(iPoint,dahist->GetErrorXlow(iPoint),dahist->GetErrorXhigh(iPoint),
                          dahist->GetErrorYlow(iPoint)/mchist->GetBinContent(iPoint+1),dahist->GetErrorYhigh(iPoint)/mchist->GetBinContent(iPoint+1));
  }

  

  TH1F* band = (TH1F*) tohist->Clone("band");
  tohist->Divide(mchist);
  tohist->SetLineColor(0);
  tohist->SetMarkerColor(0);
  tohist->SetMarkerSize(0);
  tohist->SetFillColor(kGray);

  dahist->SetMarkerSize(1);
  dphist->SetMarkerSize(1);
  dahist->SetMarkerStyle(20);
  dphist->SetMarkerStyle(24);

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
  if(!blind){
    if(addPreFitOnPull and addPullPlot)
      dphist->Draw("P0E1 SAME");
    else if(not addPullPlot)
      dphist->Draw("P0E1 SAME");
  }

  dahist->Draw("P0E1 SAME");  
  pad2->RedrawAxis("G sameaxis");

  TLegend* leg2 = NULL;
  if(addPreFitOnPull  and addPullPlot)
    leg2 = new TLegend(0.14,0.32,0.40,0.36,NULL,"brNDC");
  else
    leg2 = new TLegend(0.14,0.24,0.40,0.28,NULL,"brNDC");

  leg2->SetFillColor(0);
  leg2->SetFillStyle(1);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(dahist,"Post-fit","PLE");
  leg2->AddEntry(dphist,"Pre-fit","PLE");
  if(addPullPlot and addPreFitOnPull)
    leg2->Draw("same");
  else if(not addPullPlot)
    leg2->Draw("same");

  canvas->cd();
  pad2->RedrawAxis("sameaxis");
  canvas->RedrawAxis("sameaxis");

  if(addPullPlot){

    pad3->Draw();
    pad3->cd();

    TH1* frame3 = (TH1*) tohist->Clone("frame2");
    frame3->Reset();
    frame3->SetLineColor(kBlack);
    frame3->SetLineWidth(1);
    frame3->GetYaxis()->SetRangeUser(-3.5,3.5);
    if(category == Category::monojet)
      frame3->GetXaxis()->SetNdivisions(510);
    else
      frame3->GetXaxis()->SetNdivisions(210);

    frame3->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    if(addStatUncPull)
      frame3->GetYaxis()->SetTitle("#frac{(Data-Pred.)}{#sigma}");
    else
      frame3->GetYaxis()->SetTitle("#frac{(Data-Pred.)}{#sigma_{pred}}");

    frame3->GetYaxis()->CenterTitle();
    frame3->GetYaxis()->SetTitleOffset(1.5);
    frame3->GetXaxis()->SetTitleOffset(0.9);
    frame3->GetYaxis()->SetLabelSize(0.03);
    frame3->GetYaxis()->SetTitleSize(0.03);
    frame3->GetXaxis()->SetLabelSize(0.04);
    frame3->GetXaxis()->SetTitleSize(0.05);
    frame3->GetYaxis()->SetNdivisions(504);
    frame3->GetXaxis()->SetTickLength(0.025);
    frame3->Draw("AXIS");
    frame3->Draw("AXIG same");

    TH1F* data_pull_post = (TH1F*) tohist->Clone("data_pull_post");
    data_pull_post->Reset();
    for(int iPoint = 0; iPoint < dthist->GetN(); iPoint++){
      double x,y;
      dthist->GetPoint(iPoint,x,y);
      data_pull_post->SetBinContent(iPoint+1,y);
      data_pull_post->SetBinError(iPoint+1,(dthist->GetErrorYlow(iPoint)+dthist->GetErrorYhigh(iPoint))/2);
    }
    data_pull_post->Add(mchist,-1);
    data_pull_post->SetMarkerColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetLineColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetFillColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetLineWidth(1);
    for(int iBin = 0; iBin < data_pull_post->GetNbinsX()+1; iBin++){
      if(addStatUncPull)
        data_pull_post->SetBinContent(iBin+1,data_pull_post->GetBinContent(iBin+1)/sqrt(pow(band->GetBinError(iBin+1),2)+pow((dthist->GetErrorYlow(iBin)+dthist->GetErrorYhigh(iBin))/2,2)));
      else
        data_pull_post->SetBinContent(iBin+1,data_pull_post->GetBinContent(iBin+1)/band->GetBinError(iBin+1));
      data_pull_post->SetBinError(iBin+1,+1); // divide by sigma data                                                                                                                                 
    }

    // line at 1                                                                                                                                                                                      
    TH1* unhist2 = (TH1*) tohist->Clone("unhist");
    unhist2->Reset();
    for (int i = 1; i <= unhist2->GetNbinsX(); i++) unhist2->SetBinContent(i, 0);
    for (int i = 1; i <= unhist2->GetNbinsX(); i++) unhist2->SetBinError(i, 0);
    unhist2->SetMarkerSize(0);
    unhist2->SetLineColor(kBlack);
    unhist2->SetLineStyle(2);
    unhist2->SetLineWidth(2);
    unhist2->SetFillColor(0);
    unhist2->Draw("SAME");
    data_pull_post->Draw("hist same");
    pad3->RedrawAxis("G sameaxis");
    pad3->Modified();

  }


  if(blind and not addPullPlot){
    canvas->SaveAs(("postfit_sig_blind"+postfix+".pdf").c_str());
    canvas->SaveAs(("postfit_sig_blind"+postfix+".png").c_str());
  }
  else if(not blind and not addPullPlot){
    canvas->SaveAs(("postfit_sig"+postfix+".pdf").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+".png").c_str());
  }
  else if(addPullPlot and blind and not addPreFitOnPull){
    canvas->SaveAs(("postfit_sig_blind"+postfix+"_pull.pdf").c_str());
    canvas->SaveAs(("postfit_sig_blind"+postfix+"_pull.png").c_str());
  }
  else if(addPullPlot and not blind and not addPreFitOnPull){
    canvas->SaveAs(("postfit_sig"+postfix+"_pull.pdf").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+"_pull.png").c_str());
  }
  else if(addPullPlot and blind and addPreFitOnPull){
    canvas->SaveAs(("postfit_sig_blind"+postfix+"_wprefit_pull.pdf").c_str());
    canvas->SaveAs(("postfit_sig_blind"+postfix+"_wprefit_pull.png").c_str());
  }
  else if(addPullPlot and not blind and addPreFitOnPull){
    canvas->SaveAs(("postfit_sig"+postfix+"_wprefit_pull.pdf").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+"_wprefit_pull.png").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+"_wprefit_pull.C").c_str());
  }
  
  
  if(plotSignificance and not addPullPlot) {

    // ratio post fit at 1 with uncertaitny
    TH1* htemp = (TH1*) ((TH1*)pfile->Get((fit_dir+"/"+dir+"/total_background").c_str()))->Clone("postfit_over_prefit");
    TH1* totalSignal   = NULL;
    TH1* totalSignal_s = NULL;
    TH1* totalSignal_av = NULL;

    if(not plotSBFit){
      totalSignal_s  = (TH1*) ggHhist->Clone("totalSignal_s");
      totalSignal_av = (TH1*) mjhist_av->Clone("totalSignal_av");
    }
    else if(plotSBFit)
      totalSignal = (TH1*) sighist->Clone("totalSignal");

    canvas->cd();
    pad2->cd();

    if(not plotSBFit)
      frame2->GetYaxis()->SetTitle("(S+B)/B");
    else
      frame2->GetYaxis()->SetTitle("(S_{fit}+B)/B");
    TH1* SoverB_s  = NULL; 
    TH1* SoverB_av = NULL;
    TH1* SoverB    = NULL;

    if(not plotSBFit){
      SoverB_s  = (TH1*) totalSignal_s->Clone("SoverB_s");
      SoverB_av = (TH1*) totalSignal_av->Clone("SoverB_av");
    }
    else
      SoverB = (TH1*) totalSignal->Clone("SoverB");

    if(not plotSBFit and SoverB_s and SoverB_av){
      SoverB_s->SetLineColor(kBlack);
      SoverB_s->SetMarkerColor(kBlack);
      SoverB_s->SetMarkerSize(1);
      SoverB_s->SetMarkerStyle(20);
      SoverB_av->SetLineColor(TColor::GetColor("#0066ff"));
      SoverB_av->SetMarkerColor(TColor::GetColor("#0066ff"));
      SoverB_av->SetMarkerSize(1);
      SoverB_av->SetMarkerStyle(20);
      SoverB_s->Add(htemp);
      SoverB_s->Divide(htemp);
      SoverB_av->Add(htemp);
      SoverB_av->Divide(htemp);
    }
    else{
      SoverB->SetLineColor(kBlack);
      SoverB->SetMarkerColor(kBlack);
      SoverB->SetMarkerSize(1);
      SoverB->SetMarkerStyle(20);      
      SoverB->Add(htemp);
      SoverB->Divide(htemp);
    }

    if(not plotSBFit)
      //      frame2->GetYaxis()->SetRangeUser(SoverB_av->GetMinimum()*0.9,SoverB_av->GetMaximum()*1.1);
      frame2->GetYaxis()->SetRangeUser(1,2);
    else
      frame2->GetYaxis()->SetRangeUser(0.5,SoverB->GetMaximum()*1.2);
    frame2->Draw();

    TH1* SoverB_postfit_d = NULL;
    if(not plotSBFit)
      SoverB_postfit_d = (TH1*) SoverB_s->Clone("SoverB_postfit_d");
    else
      SoverB_postfit_d = (TH1*) SoverB->Clone("SoverB_postfit_d");

    for(int iBin = 0; iBin < SoverB_postfit_d->GetNbinsX(); iBin++)
      SoverB_postfit_d->SetBinContent(iBin+1,1);
    SoverB_postfit_d->SetLineColor(0);
    
    SoverB_postfit_d->SetMarkerColor(0);
    SoverB_postfit_d->SetMarkerSize(0);
    SoverB_postfit_d->SetFillColor(kGray);
    SoverB_postfit_d->SetFillStyle(1001);
    SoverB_postfit_d->Draw("E2 SAME");
    unhist->Draw("SAME");
    if(not plotSBFit){
      //SoverB_s->Draw("hist same");
      SoverB_av->Draw("hist same");
    }
    else
      SoverB->Draw("hist same");

    pad2->RedrawAxis("sameaxis");

    canvas->SaveAs(("postfit_sig"+postfix+"_SoB.pdf").c_str());
    canvas->SaveAs(("postfit_sig"+postfix+"_SoB.png").c_str());
    
    if(plotSignificance){
      TFile* outFile = new TFile("postfit_weights_Sig.root","RECREATE");
      outFile->cd();
      htemp->Divide(tphist);
      htemp->Write("postfit_over_prefit");
      outFile->Close();
    }
  }
}

