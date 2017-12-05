1;95;0c#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void makePostFit_VBF(string   fitFilename, 
		     bool     plotSBFit  = false,
		     bool     blind = false){
  
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();


  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 800);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.38);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,1.);
  pad2->SetTopMargin(0.63);
  pad2->SetBottomMargin(0.25);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->SetLineColor(0);

  TPad*  pad3 = new TPad("pad3","pad3",0,0.,1,1.);
  pad3->SetTopMargin(0.76);
  pad3->SetRightMargin(0.06);
  pad3->SetFillColor(0);
  pad3->SetFillStyle(0);
  pad3->SetLineColor(0);

  // open the input file
  TFile* pfile = new TFile(fitFilename.c_str());  
  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  /// only signal region VBF only fit are produced
  string dir = "ch1";
  string postfix = "VBF";
  
  // in case of higgs invisible limits
  TH1* ggHhist = (TH1*) pfile->Get(("shapes_prefit/"+dir+"/ggH_hinv").c_str());;
  TH1* vbfhist = (TH1*) pfile->Get(("shapes_prefit/"+dir+"/qqH_hinv").c_str());;
  
  TH1* znhist = NULL;
  TH1* zlhist = NULL;
  TH1* wlhist = NULL;
  TH1* tthist = NULL;
  TH1* dihist = NULL;
  TH1* qchist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* tohist = NULL;
  TH1* tphist = NULL;
  TH1* sighist = NULL;

  znhist = (TH1*) pfile->Get((fit_dir+"/"+dir+"/qcd_znunu").c_str());    
  zlhist = (TH1*) pfile->Get((fit_dir+"/"+dir+"/ZJets").c_str());    
  wlhist = (TH1*) pfile->Get((fit_dir+"/"+dir+"/qcd_wjets").c_str());    
  tthist = (TH1*) pfile->Get((fit_dir+"/"+dir+"/Top").c_str());    
  dihist = (TH1*) pfile->Get((fit_dir+"/"+dir+"/Dibosons").c_str());    
  ewkwhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/ewk_wjets").c_str());    
  ewkzhist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/ewk_znunu").c_str());    
  qchist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/QCD").c_str());    

  ///
  tohist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/total_background").c_str());    
  tphist = (TH1*)pfile->Get(("shapes_prefit/"+dir+"/total_background").c_str());    

  if(plotSBFit)
    sighist = (TH1*)pfile->Get((fit_dir+"/"+dir+"/total_signal").c_str());

  ////////////////
  TGraphAsymmErrors* dthist = NULL;
  if(!blind)
    dthist = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir+"/data").c_str());
  else{
    dthist = new TGraphAsymmErrors();
    for(int iBin = 0; iBin < tohist->GetNbinsX()+1; iBin++){
      dthist->SetPoint(iBin,tohist->GetBinCenter(iBin+1),tohist->GetBinContent(iBin+1));
      dthist->SetPointError(iBin,tohist->GetBinWidth(iBin+1)/2,tohist->GetBinWidth(iBin+1)/2,tohist->GetBinError(iBin+1),tohist->GetBinError(iBin+1));
    }
  }
  
  // start to change style
  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kRed);
    ggHhist->SetLineWidth(3);
    ggHhist->SetMarkerSize(0);
  }
    
  if(vbfhist){
    vbfhist->SetFillColor(0);
    vbfhist->SetFillStyle(0);
    vbfhist->SetLineColor(kBlack);
    vbfhist->SetLineWidth(3);
    vbfhist->SetMarkerSize(0);
  }

  if(qchist){
    qchist->SetFillColor(TColor::GetColor("#F1F1F2"));
    qchist->SetLineColor(kBlack);
  }

  zlhist->SetFillColor(TColor::GetColor("#9A9EAB"));  
  zlhist->SetLineColor(kBlack);
  
  znhist->SetFillColor(TColor::GetColor("#3A8C4C"));
  znhist->SetLineColor(kBlack);

  wlhist->SetFillColor(TColor::GetColor("#FAAF08"));
  wlhist->SetLineColor(kBlack);

  dihist->SetFillColor(kRed+3);
  dihist->SetLineColor(kBlack);

  tthist->SetFillColor(TColor::GetColor("#CF3721"));
  tthist->SetLineColor(kBlack);

  if(ewkzhist){
    ewkzhist->SetFillColor(kCyan+1);
    ewkzhist->SetLineColor(kBlack);
  }
  if(ewkwhist){
    ewkwhist->SetFillColor(kAzure+1);
    ewkwhist->SetLineColor(kBlack);
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
  stack->Add(qchist);
  stack->Add(zlhist); 
  stack->Add(tthist);
  stack->Add(dihist);  
  stack->Add(ewkwhist);
  stack->Add(ewkzhist);
  stack->Add(wlhist);
  stack->Add(znhist);

  TH1* frame = (TH1*) tohist->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);
  frame->GetYaxis()->SetRangeUser(0.0007,tphist->GetMaximum()*1000);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->GetXaxis()->SetNdivisions(504);
  frame->Draw();
  
  CMS_lumi(canvas,"35.9");

  TLatex* categoryLabel = new TLatex();
  categoryLabel->SetNDC();
  categoryLabel->SetTextSize(0.5*canvas->GetTopMargin());
  categoryLabel->SetTextFont(42);
  categoryLabel->SetTextAlign(11);
  categoryLabel ->DrawLatex(0.175,0.81,"VBF-tagged");
  categoryLabel->Draw("same");
  stack ->Draw("HIST SAME");

  if(vbfhist and !plotSBFit){
    vbfhist->SetLineStyle(2);
    vbfhist->Draw("HIST SAME");
  }

  if(ggHhist and !plotSBFit){
    ggHhist->SetLineStyle(2);
    ggHhist->Draw("HIST SAME");
  }

  if(plotSBFit)
    sighist->Draw("HIST same");
  
  dthist->SetMarkerSize(1.2);
  dthist->SetMarkerStyle(20);
  dthist->SetFillStyle(0);
  dthist->SetFillColor(0);
  dthist->SetLineColor(kBlack);
  dthist->SetLineWidth(1);
  dthist->SetMarkerColor(kBlack);
  dthist->Draw("PE SAME");

  TLegend* leg = new TLegend(0.40, 0.60, 0.95, 0.90);
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  leg->AddEntry(dthist,"Data", "PEL");
  leg->AddEntry(znhist,"Z+jets (QCD)", "F");
  leg->AddEntry(wlhist,"W+jets (QCD)", "F");
  leg->AddEntry(ewkzhist,"Z+jets (EW)", "F");
  leg->AddEntry(ewkzhist,"W+jets (EW)", "F");
  leg->AddEntry(dihist,  "VV", "F");
  leg->AddEntry(tthist,  "Top quark", "F");
  leg->AddEntry(zlhist,  "Z/#gamma(ll), #gamma+jets", "F");
  leg->AddEntry(qchist,  "QCD", "F");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted signal", "L");
  else if(not plotSBFit){
    leg->AddEntry(vbfhist, "qqH(125) #rightarrow inv.", "L");
    leg->AddEntry(ggHhist, "ggH(125) #rightarrow inv.", "L");
  }

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
  frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  frame2->GetXaxis()->SetNdivisions(505);
  frame2->GetYaxis()->SetNdivisions(5);

  frame2->GetYaxis()->SetTitleOffset(1.9);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetXaxis()->SetLabelSize(0);
  frame2->GetXaxis()->SetTitleSize(0);
  frame2->GetYaxis()->SetTitleSize(0.03);
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
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
  mchist->Add(ewkwhist);
  mchist->Add(ewkzhist);

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
  if(!blind)
    dphist->Draw("P0E1 SAME");
  dahist->Draw("P0E1 SAME");  
  pad2->RedrawAxis("G sameaxis");

  canvas->cd();
  pad2->RedrawAxis("sameaxis");
  canvas->RedrawAxis("sameaxis");


  //////////////
  pad3->Draw();
  pad3->cd();

  TH1* frame3 = (TH1*) tohist->Clone("frame2");
  frame3->Reset();
  frame3->SetLineColor(kBlack);
  frame3->SetLineWidth(1);
  frame3->GetYaxis()->SetRangeUser(-3,3);
  frame3->GetXaxis()->SetNdivisions(505);
  frame3->GetXaxis()->SetTitle("M_{jj} [GeV]");
  frame3->GetYaxis()->SetTitle("#frac{(Data-Pred.)}{Unc.}");

  frame3->GetYaxis()->CenterTitle();
  frame3->GetYaxis()->SetTitleOffset(1.5);
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
    data_pull_post->SetBinContent(iBin+1,data_pull_post->GetBinContent(iBin+1)/sqrt(pow(band->GetBinError(iBin+1),2)+pow((dthist->GetErrorYlow(iBin)+dthist->GetErrorYhigh(iBin))/2,2)));
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

  canvas->SaveAs(("postfit_sig_"+postfix+".pdf").c_str());
  canvas->SaveAs(("postfit_sig_"+postfix+".png").c_str());
}

