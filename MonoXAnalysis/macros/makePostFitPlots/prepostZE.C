#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static bool saveTextFile = false;
static bool dumpInfo     = false;


void prepostZE(string fitFilename, string observable, Category category, bool isCombinedFit = false, bool plotSBFit = false, bool addPullPlot = false,  bool dumpHisto = false) {

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
  

  TFile* pfile = new TFile(fitFilename.c_str());

  TGraphAsymmErrors* dthist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* pohist = NULL;
  TH1* prhist = NULL;


  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  string dir;
  if(isCombinedFit){
    if(category == Category::monojet)
      dir = "ch1_ch5";
    else if(category == Category::monoV)
      dir = "ch2_ch5";
    else if(category == Category::VBF)
      dir = "ch3_ch4";
  }
  else if( category != Category::VBF)
    dir = "ch5";
  else
    dir = "ch4";


  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF)
    postfix = "_VBF";


  dthist = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir+"/data").c_str());
  wlhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/WJets_ZE").c_str());
  tthist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Top").c_str());
  dihist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/Dibosons").c_str());
  ewkwhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/WJets_EWK_ZE").c_str());
  ewkzhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/ZJets_EWK").c_str());
  pohist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/total_background").c_str());
  prhist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());
  
  if(saveTextFile){

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
    stringstream ZJetRate;
    ZJetRate << "Process: ZJet";
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
    
    if(category == Category::VBF){
      for(int iBin = 0; iBin < ewkwhist->GetNbinsX(); iBin++){
	EWKWRate << "   ";
	EWKWRate << ewkwhist->GetBinContent(iBin);
      }
      
      for(int iBin = 0; iBin < ewkzhist->GetNbinsX(); iBin++){
	EWKZRate << "   ";
	EWKZRate << ewkzhist->GetBinContent(iBin);
      }
    }
    
    for(int iBin = 0; iBin < wlhist->GetNbinsX(); iBin++){
      ZJetRate << "   ";
      ZJetRate << wlhist->GetBinContent(iBin);
    }
    
    for(int iBin = 0; iBin < prhist->GetNbinsX(); iBin++){
      PreRate << "   ";
      PreRate << prhist->GetBinContent(iBin);
    }
    
    for(int iBin = 0; iBin < pohist->GetNbinsX(); iBin++){
      PostRate << "   ";
      PostRate << pohist->GetBinContent(iBin);
    }  

    for(int iBin = 0; iBin < dthist->GetN(); iBin++){
      double x,y;
      dthist->GetPoint(iBin,x,y);
      DataRate << "   ";
      DataRate << y;
    }
          
    outputfile<<"######################"<<endl;
    outputfile<<TopRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<VVRate.str()<<endl;
    outputfile<<"######################"<<endl;
    if(category == Category::VBF){
      outputfile<<EWKWRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<EWKZRate.str()<<endl;
      outputfile<<"######################"<<endl;
    }
    outputfile<<ZJetRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<PreRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<PostRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<DataRate.str()<<endl;
    
    outputfile.close();
  }
  
  prhist->SetLineColor(kRed);
  prhist->SetLineStyle(7);
  prhist->SetLineWidth(2);
  pohist->SetLineColor(TColor::GetColor("#0066ff"));
  pohist->SetLineWidth(2);
  prhist->SetMarkerColor(kRed);
  pohist->SetMarkerColor(TColor::GetColor("#0066ff"));
  
  //  wlhist->SetFillColor(kOrange+1);
  wlhist->SetFillColor(33);
  wlhist->SetLineColor(kBlack);
  wlhist->Add(tthist);
  wlhist->Add(dihist);
  if(category == Category::VBF){
    wlhist->Add(ewkwhist);
    ewkzhist->SetFillColor(kCyan+1);
    ewkzhist->SetLineColor(kBlack);
  }
  else{
    if(ewkwhist) wlhist->Add(ewkwhist);
    if(ewkzhist) wlhist->Add(ewkzhist);
  }
  
  TH1* frame = (TH1*) pohist->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);

  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,prhist->GetMaximum()*500);
  else if(category == Category::monoV)
    frame->GetYaxis()->SetRangeUser(0.0007,prhist->GetMaximum()*500);
  else if(category == Category::VBF)
    frame->GetYaxis()->SetRangeUser(0.0007,prhist->GetMaximum()*500);

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

  prhist->Draw("HIST SAME");
  pohist->Draw("HIST SAME");
  wlhist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetLineColor(kBlack);
  dthist->Draw("EP SAME");
  
 
  TLegend* leg = new TLegend(0.50, 0.62, 0.95, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","PEL");
  leg->AddEntry(pohist, "Post-fit Z(ee)+jets","L");
  leg->AddEntry(prhist, "Pre-fit Z(ee)+jets","L");
  if(category == Category::VBF)
    leg->AddEntry(ewkwhist, "Z-EWK","F");
  leg->AddEntry(wlhist, "Other backgrounds", "F");
  leg->Draw("SAME");
  
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  
  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  (TH1*) pohist->Clone("frame");
  frame2->Reset();
  frame2->SetLineColor(kBlack);
  frame2->SetLineWidth(1);

  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.4,1.6);
  else
    frame2->GetYaxis()->SetRangeUser(0.4,1.6);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(210);
  //  frame2->GetYaxis()->SetNdivisions(5);
  frame2->GetYaxis()->SetNdivisions(3);

  if(not addPullPlot){
    frame2->GetXaxis()->SetTitle("Hadronic recoil p_{T} [GeV]");
    if(category == Category::VBF and TString(observable).Contains("mjj"))
      frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");
    frame2->GetYaxis()->SetTitle("Data / Pred.");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetTitleOffset(1.5);
    frame2->GetXaxis()->SetTitleOffset(1.1);
    frame2->GetYaxis()->SetLabelSize(0.04);
    frame2->GetYaxis()->SetTitleSize(0.04);
    frame2->GetXaxis()->SetLabelSize(0.04);
    frame2->GetXaxis()->SetTitleSize(0.05);
  }
  else{
    frame2->GetYaxis()->SetTitleOffset(1.9);
    frame2->GetXaxis()->SetTitleOffset(1.1);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetXaxis()->SetLabelSize(0);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetYaxis()->SetTitle("Data / Pred.");
    frame2->GetYaxis()->CenterTitle();
  }
  frame2->GetXaxis()->SetTickLength(0.025);
  frame2->Draw();

  TGraphAsymmErrors* d1hist = (TGraphAsymmErrors*)dthist->Clone("d1hist");
  TGraphAsymmErrors* d2hist = (TGraphAsymmErrors*)dthist->Clone("d2hist");
  TH1* m1hist = (TH1*)prhist->Clone("m1hist");
  TH1* m2hist = (TH1*)pohist->Clone("m2hist");
  TH1* erhist = (TH1*)pohist->Clone("erhist");
  
  d1hist->SetLineColor(kRed);
  d2hist->SetLineColor(TColor::GetColor("#0066ff"));
  d1hist->SetLineColor(kRed);
  d2hist->SetLineColor(TColor::GetColor("#0066ff"));
  d1hist->SetMarkerColor(kRed);
  d1hist->SetMarkerSize(1);
  d1hist->SetMarkerStyle(24);
  d2hist->SetMarkerColor(TColor::GetColor("#0066ff"));
  d2hist->SetMarkerSize(1);
  d2hist->SetMarkerStyle(20);
  
  for (int i = 1; i <= m1hist->GetNbinsX(); i++) m1hist->SetBinError(i, 0);
  for (int i = 1; i <= m2hist->GetNbinsX(); i++) m2hist->SetBinError(i, 0);

  for(int iPoint = 0; iPoint < d1hist->GetN(); iPoint++){
    double x,y;
    d1hist->GetPoint(iPoint,x,y);
    d1hist->SetPoint(iPoint,x,y/m1hist->GetBinContent(iPoint+1));
    d1hist->SetPointError(iPoint,d1hist->GetErrorXlow(iPoint),d1hist->GetErrorXhigh(iPoint),
                          d1hist->GetErrorYlow(iPoint)/m1hist->GetBinContent(iPoint+1),d1hist->GetErrorYhigh(iPoint)/m1hist->GetBinContent(iPoint+1));
    d2hist->GetPoint(iPoint,x,y);
    d2hist->SetPoint(iPoint,x,y/m2hist->GetBinContent(iPoint+1));
    d2hist->SetPointError(iPoint,d2hist->GetErrorXlow(iPoint),d2hist->GetErrorXhigh(iPoint),
                          d2hist->GetErrorYlow(iPoint)/m2hist->GetBinContent(iPoint+1),d2hist->GetErrorYhigh(iPoint)/m2hist->GetBinContent(iPoint+1));
  }
  
  erhist->Divide(m2hist);
  erhist->SetLineColor(0);
  erhist->SetMarkerColor(0);
  erhist->SetMarkerSize(0);
  erhist->SetFillColor(kGray);
  
  d1hist->SetMarkerSize(1.2);
  d2hist->SetMarkerSize(1.2);
  
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
  
  //  if(not addPullPlot)
  d1hist->Draw("PE1 SAME");    
  d2hist->Draw("PE1 SAME");
  erhist->Draw("E2 SAME");
  //if(not addPullPlot)
  d1hist->Draw("P0E1 SAME");
  d2hist->Draw("P0E1 SAME");

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

  TLegend* leg2 = new TLegend(0.14,0.24,0.40,0.28,NULL,"brNDC");
  //leg2->SetFillColor(0);
  //leg2->SetFillStyle(1);
  //leg2->SetLineColor(0);
  leg2->SetBorderSize(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(d2hist,"Post-fit","PLE");
  leg2->AddEntry(d1hist,"Pre-fit","PLE");
  //  if(not addPullPlot)
    leg2->Draw("same");


  if(addPullPlot){

    pad3->Draw();
    pad3->cd();

    TH1* frame3 = (TH1*) pohist->Clone("frame2");
    frame3->Reset();
    frame3->SetLineColor(kBlack);
    frame3->SetLineWidth(1);
    frame3->GetYaxis()->SetRangeUser(-3.5,3.5);
    if(category == Category::monojet)
      frame3->GetXaxis()->SetNdivisions(510);
    else
      frame3->GetXaxis()->SetNdivisions(210);
    frame3->GetXaxis()->SetTitle("Hadronic recoil p_{T} [GeV]");
    frame3->GetYaxis()->SetTitle("#frac{(Data-Pred.)}{#sigma_{pred}}");

    frame3->GetYaxis()->CenterTitle();
    frame3->GetYaxis()->SetTitleOffset(1.5);
    frame3->GetYaxis()->SetLabelSize(0.03);
    frame3->GetYaxis()->SetTitleSize(0.03);
    frame3->GetXaxis()->SetLabelSize(0.04);
    frame3->GetXaxis()->SetTitleSize(0.05);
    frame3->GetYaxis()->SetNdivisions(504);
    frame3->Draw("AXIS");
    frame3->Draw("AXIG same");
    
    TH1F* data_pull_post = (TH1F*) pohist->Clone("data_pull_post");
    data_pull_post->Reset();
    for(int iPoint = 0; iPoint < dthist->GetN(); iPoint++){
      double x,y;
      dthist->GetPoint(iPoint,x,y);
      data_pull_post->SetBinContent(iPoint+1,y);
      data_pull_post->SetBinError(iPoint+1,(dthist->GetErrorYlow(iPoint+1)+dthist->GetErrorYhigh(iPoint+1))/2);
    }
    data_pull_post->Add(pohist,-1);
    data_pull_post->SetMarkerColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetLineColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetFillColor(TColor::GetColor("#0066ff"));
    data_pull_post->SetLineWidth(1);
    for(int iBin = 0; iBin < data_pull_post->GetNbinsX()+1; iBin++){
      data_pull_post->SetBinContent(iBin+1,data_pull_post->GetBinContent(iBin+1)/pohist->GetBinError(iBin+1)); // divide by sigma data                                                               
      data_pull_post->SetBinError(iBin+1,+1); // divide by sigma data                                                                                                                                 
    }

    // line at 1                                                                                                                                                                                      
    TH1* unhist2 = (TH1*) pohist->Clone("unhist");
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

  if(not addPullPlot){
    canvas->SaveAs(("prepostfit_zee"+postfix+".pdf").c_str());
    canvas->SaveAs(("prepostfit_zee"+postfix+".png").c_str());
  }
  else{
    canvas->SaveAs(("prepostfit_zee"+postfix+"_pull.pdf").c_str());
    canvas->SaveAs(("prepostfit_zee"+postfix+"_pull.png").c_str());
  }

  if(dumpHisto){

    TFile* outFile = new TFile(("postfit_weights_ZE"+postfix+".root").c_str(),"RECREATE");
    outFile->cd();

    dthist->Write("data");
    wlhist->Write("zjets_post_fit");
    tthist->Write("top_post_fit");
    dihist->Write("diboson_post_fit");

    TH1* wlhist_prefit = (TH1*) pfile->Get(("shapes_prefit/"+dir+"/WJets_ZE").c_str());
    TH1* tthist_prefit = (TH1*) pfile->Get(("shapes_prefit/"+dir+"/Top").c_str());
    TH1* dihist_prefit = (TH1*) pfile->Get(("shapes_prefit/"+dir+"/Dibosons").c_str());
    
    wlhist_prefit->Write("zjets_pre_fit");
    tthist_prefit->Write("top_pre_fit");
    dihist_prefit->Write("diboson_pre_fit");

    pohist->Write("post_fit_post_fit");
    prhist->Write("pre_fit_post_fit");

    outFile->Close();
  }
}

