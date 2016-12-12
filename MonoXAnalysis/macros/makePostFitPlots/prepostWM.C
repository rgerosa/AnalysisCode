#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static bool saveTextFile = false;

void prepostWM(string fitFilename, string templateFileName, string observable, Category category, bool isLog = true, bool plotSBFit = false,  bool dumpHisto = false) {

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
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* qcdhist = NULL;
  TH1* pohist = NULL;
  TH1* prhist = NULL;

  string dir = "ch3";
  if(category == Category::VBF)
    dir = "ch3";

  string dir_fit = "shapes_fit_b";
  if(plotSBFit)
    dir_fit = "shapes_fit_s";

  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF)
    postfix = "_VBF";
    
  dthist = (TH1*)dfile->FindObjectAny(("datahistwmn_"+observable).c_str());
  dthist->Scale(1.0, "width");

  wlhist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/ZJets_WM").c_str());
  tthist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/Top").c_str());
  dihist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/Dibosons").c_str());
  qcdhist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/QCD_WM").c_str());
  if(category == Category::VBF){
    ewkwhist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/WJets_EWK").c_str());
    ewkzhist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/ZJets_EWK_WM").c_str());
  }
  pohist = (TH1*)pfile->Get((dir_fit+"/"+dir+"/total_background").c_str());
  prhist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());

  // write down rates
  if(saveTextFile){

    ofstream  outputfile;
    outputfile.open("prepostWM.txt");
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
    
    for(int iBin = 0; iBin < dthist->GetNbinsX(); iBin++){
      DataRate << "   ";
      DataRate << dthist->GetBinContent(iBin);
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

  // use colors for plotting  
  prhist->SetLineColor(kRed);
  prhist->SetLineWidth(2);
  prhist->SetMarkerColor(kRed);
  pohist->SetLineColor(kBlue);
  pohist->SetLineWidth(2);
  pohist->SetMarkerColor(kBlue);
  
  wlhist->SetFillColor(kOrange+1);
  wlhist->SetLineColor(kBlack);
  wlhist->Add(qcdhist);
  wlhist->Add(tthist);
  wlhist->Add(dihist);
  if(category == Category::VBF){
    wlhist->Add(ewkzhist);
    ewkwhist->SetFillColor(kCyan+1);
    ewkwhist->SetLineColor(kBlack);
  }

  
  TH1* frame = (TH1*) dthist->Clone("frame");
  frame->Reset();
  if(isLog){
    if(category == Category::monojet)
      frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*200);
    else if(category == Category::monoV)
      frame->GetYaxis()->SetRangeUser(0.002,wlhist->GetMaximum()*200);
    else if(category == Category::VBF){
      if(TString(observable).Contains("detajj"))
	frame->GetYaxis()->SetRangeUser(0.1,wlhist->GetMaximum()*300);
      else
	frame->GetYaxis()->SetRangeUser(0.01,wlhist->GetMaximum()*300);
    }
  }
  else 
    frame->GetYaxis()->SetRangeUser(0,dthist->GetMaximum()*2);
    
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else if(category == Category::monoV)
    frame->GetXaxis()->SetNdivisions(504);
  else if(category == Category::VBF)
    frame->GetXaxis()->SetNdivisions(504);

  frame->Draw();

  ///
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
  if(category ==  Category::VBF)
    ewkwhist->Draw("HIST SAME");
  wlhist->Draw("HIST SAME");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetLineColor(kBlack);
  dthist->Draw("EP SAME");
  
  ///
  TLegend* leg = new TLegend(0.50, 0.62, 0.95, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(dthist, "Data","PEL");
  leg->AddEntry(pohist, "Post-fit single-muon CR","L");
  leg->AddEntry(prhist, "Pre-fit single-muon CR","L");
  if(category == Category::VBF)
    leg->AddEntry(ewkwhist, "W-EWK","F");
  leg->AddEntry(wlhist, "Other Backgrounds", "F");
  leg->Draw("SAME");  
  canvas->RedrawAxis("sameaxis");
  if(isLog)
    canvas->SetLogy();

  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  (TH1*) dthist->Clone("frame");
  frame2->Reset("ICES");

  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if(category == Category::monoV)
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if(category == Category::VBF)
    frame2->GetYaxis()->SetRangeUser(0.0,2.0);
  
  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else if(category == Category::monoV)
    frame2->GetXaxis()->SetNdivisions(210);
  else
    frame2->GetXaxis()->SetNdivisions(505);  
  frame2->GetYaxis()->SetNdivisions(5);

  if(TString(observable).Contains("met") and not TString(observable).Contains("jetmetdphi"))
    frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  else if(TString(observable).Contains("mjj"))
    frame2->GetXaxis()->SetTitle("m_{jj} [GeV]");
  else if(TString(observable).Contains("detajj"))
    frame2->GetXaxis()->SetTitle("#Delta#eta_{jj}");

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
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1);
  leg2->SetBorderSize(0);
  leg2->SetLineColor(0);
  leg2->SetNColumns(2);
  leg2->AddEntry(d2hist,"post-fit","PLE");
  leg2->AddEntry(d1hist,"pre-fit","PLE");
  leg2->Draw("same");

  canvas->SaveAs(("prepostfit_wmn"+postfix+".pdf").c_str());
  canvas->SaveAs(("prepostfit_wmn"+postfix+".png").c_str());

  if(dumpHisto){

    TFile* outFile = new TFile("postfit_weights_WM.root","RECREATE");
    outFile->cd();

    dthist->Write("data");
    wlhist->Write("zjets_post_fit");
    tthist->Write("top_post_fit");
    dihist->Write("diboson_post_fit");

    TH1* wlhist_prefit = (TH1*) pfile->Get("shapes_prefit/ch3/ZJets_WM");
    TH1* tthist_prefit = (TH1*) pfile->Get("shapes_prefit/ch3/Top");
    TH1* dihist_prefit = (TH1*) pfile->Get("shapes_prefit/ch3/Dibosons");

    wlhist_prefit->Write("zjets_pre_fit");
    tthist_prefit->Write("top_pre_fit");
    dihist_prefit->Write("diboson_pre_fit");

    pohist->Write("post_fit_post_fit");
    prhist->Write("pre_fit_post_fit");

    outFile->Close();
  }
  
}

