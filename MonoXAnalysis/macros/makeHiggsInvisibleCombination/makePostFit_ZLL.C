#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void makePostFit_ZLL(string   fitFilename, 
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
  string postfix = "ZLL";
  
  // in case of higgs invisible limits
  TH1F* ZHhist = new TH1F("ZHhist","",12,0,12);
  TH1F* ggZHhist = new TH1F("ggZHhist","",12,0,12);
  TH1F* EMhist = new TH1F("EMhist","",12,0,12);
  TH1F* VVVhist = new TH1F("VVVhist","",12,0,12);
  TH1F* WZhist = new TH1F("WZhist","",12,0,12);
  TH1F* ZZhist = new TH1F("ZZhist","",12,0,12);
  TH1F* Zjethist = new TH1F("Zjethist","",12,0,12);
  TH1F* total_bkg_post = new TH1F("total_bkg_post","",12,0,12);
  TH1F* total_bkg_pre = new TH1F("total_bkg_pre","",12,0,12);

  TH1F* sighist = NULL;
  if(plotSBFit)
    sighist = new TH1F("sighist","",12,0,12);

  TGraphAsymmErrors* dthist = new TGraphAsymmErrors();

  for(int i = 1; i <= 12; i++){ // loop on bins
    
    if(pfile->Get(Form("shapes_prefit/ch1_ch%d/ZH_hinv",i))){
      ZHhist->SetBinContent(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/ZH_hinv",i)))->GetBinContent(1));
      ZHhist->SetBinError(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/ZH_hinv",i)))->GetBinError(1));
    }

    if(pfile->Get(Form("shapes_prefit/ch1_ch%d/ggZH_hinv",i))){
      ggZHhist->SetBinContent(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/ggZH_hinv",i)))->GetBinContent(1));
      ggZHhist->SetBinError(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/ggZH_hinv",i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/EM",fit_dir.c_str(),i))){
      EMhist->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/EM",fit_dir.c_str(),i)))->GetBinContent(1));
      EMhist->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/EM",fit_dir.c_str(),i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/VVV",fit_dir.c_str(),i))){
      VVVhist->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/VVV",fit_dir.c_str(),i)))->GetBinContent(1));
      VVVhist->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/VVV",fit_dir.c_str(),i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/Zjets",fit_dir.c_str(),i))){
      Zjethist->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/Zjets",fit_dir.c_str(),i)))->GetBinContent(1));
      Zjethist->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/Zjets",fit_dir.c_str(),i)))->GetBinError(1));
    }
    
    if(pfile->Get(Form("%s/ch1_ch%d/WZ",fit_dir.c_str(),i))){
      WZhist->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/WZ",fit_dir.c_str(),i)))->GetBinContent(1));
      WZhist->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/WZ",fit_dir.c_str(),i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/ZZ",fit_dir.c_str(),i))){
      ZZhist->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/ZZ",fit_dir.c_str(),i)))->GetBinContent(1));
      ZZhist->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/ZZ",fit_dir.c_str(),i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/total_background",fit_dir.c_str(),i))){
      total_bkg_post->SetBinContent(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/total_background",fit_dir.c_str(),i)))->GetBinContent(1));
      total_bkg_post->SetBinError(i,((TH1*) pfile->Get(Form("%s/ch1_ch%d/total_background",fit_dir.c_str(),i)))->GetBinError(1));
    }

    if(pfile->Get(Form("shapes_prefit/ch1_ch%d/total_background",i))){
      total_bkg_pre->SetBinContent(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/total_background",i)))->GetBinContent(1));
      total_bkg_pre->SetBinError(i,((TH1*) pfile->Get(Form("shapes_prefit/ch1_ch%d/total_background",i)))->GetBinError(1));
    }

    if(pfile->Get(Form("%s/ch1_ch%d/data",fit_dir.c_str(),i))){
      if(!blind){
	double x,y;
	double x_err_dw,x_err_up,y_err_dw,y_err_up;
	TGraphAsymmErrors* dt_temp = (TGraphAsymmErrors*) pfile->Get(Form("%s/ch1_ch%d/data",fit_dir.c_str(),i));
	dt_temp->GetPoint(0,x,y);
	x_err_dw = dt_temp->GetErrorXlow(0);
	x_err_up = dt_temp->GetErrorXhigh(0);
	y_err_dw = dt_temp->GetErrorYlow(0);
	y_err_up = dt_temp->GetErrorYhigh(0);  
	dthist->SetPoint(i-1,total_bkg_post->GetBinCenter(i),y);
	dthist->SetPointError(i-1,x_err_dw,x_err_up,y_err_dw,y_err_up);
      }    
      else{
	dthist->SetPoint(i-1,total_bkg_post->GetBinCenter(i),total_bkg_post->GetBinContent(i));
	dthist->SetPointError(i-1,total_bkg_post->GetBinWidth(i)/2,total_bkg_post->GetBinWidth(i)/2,total_bkg_post->GetBinError(i),total_bkg_post->GetBinError(i));
      }
    }    
  }

  // start to change style
  if(ZHhist){
    ZHhist->SetFillColor(0);
    ZHhist->SetFillStyle(0);
    ZHhist->SetLineColor(kBlack);
    ZHhist->SetLineWidth(3);
    ZHhist->SetMarkerSize(0);
  }
    
  if(ggZHhist){
    ggZHhist->SetFillColor(0);
    ggZHhist->SetFillStyle(0);
    ggZHhist->SetLineColor(kRed);
    ggZHhist->SetLineWidth(3);
    ggZHhist->SetMarkerSize(0);
  }
  
  EMhist->SetFillColor(TColor::GetColor("#F1F1F2"));
  EMhist->SetLineColor(kBlack);

  VVVhist->SetFillColor(TColor::GetColor("#9A9EAB"));  
  VVVhist->SetLineColor(kBlack);

  Zjethist->SetFillColor(TColor::GetColor("#4897D8"));
  Zjethist->SetLineColor(kBlack);
  
  ZZhist->SetFillColor(TColor::GetColor("#3A8C4C"));
  ZZhist->SetLineColor(kBlack);

  WZhist->SetFillColor(TColor::GetColor("#FAAF08"));
  WZhist->SetLineColor(kBlack);
  
  if(sighist){
    sighist->SetFillColor(kBlack);
    sighist->SetLineColor(kBlack);
    sighist->SetLineWidth(3);
    sighist->SetFillColor(0);
    sighist->SetFillStyle(0);
  }

  // make the stack for backgrounds
  THStack* stack = new THStack("stack", "stack");
  stack->Add(EMhist);
  stack->Add(VVVhist); 
  stack->Add(Zjethist);
  stack->Add(WZhist);  
  stack->Add(ZZhist);

  TH1* frame = (TH1*) total_bkg_pre->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);
  frame->GetYaxis()->SetRangeUser(0.1,total_bkg_pre->GetMaximum()*1000);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->Draw();
  
  CMS_lumi(canvas,"35.9");
  canvas->SetLogy();

  TLatex* categoryLabel = new TLatex();
  categoryLabel->SetNDC();
  categoryLabel->SetTextSize(0.5*canvas->GetTopMargin());
  categoryLabel->SetTextFont(42);
  categoryLabel->SetTextAlign(11);
  categoryLabel ->DrawLatex(0.175,0.81,"Z(ll)-tagged");
  categoryLabel->Draw("same");
  stack ->Draw("HIST SAME");

  if(ZHhist and !plotSBFit){
    ZHhist->SetLineStyle(2);
    ZHhist->Draw("HIST SAME");
  }

  if(ggZHhist and !plotSBFit){
    ggZHhist->SetLineStyle(2);
    ggZHhist->Draw("HIST SAME");
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

  TLegend* leg = new TLegend(0.45, 0.65, 0.95, 0.92);
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  leg->AddEntry(dthist,"Data", "PEL");
  leg->AddEntry(ZZhist,"ZZ", "F");
  leg->AddEntry(WZhist,"WZ", "F");
  leg->AddEntry(Zjethist,"Drell-Yan", "F");
  leg->AddEntry(VVVhist, "VVV", "F");
  leg->AddEntry(EMhist,  "Other bkg.", "F");
  if(sighist && plotSBFit)
    leg->AddEntry(sighist, "Fitted signal", "L");
  else if(not plotSBFit){
    leg->AddEntry(ZHhist, "ZH(125) #rightarrow inv.", "L");
    leg->AddEntry(ggZHhist, "ggZH(125) #rightarrow inv.", "L");
  }

  leg->Draw("SAME");    
  canvas->RedrawAxis("sameaxis");
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  ////
  TH1* frame2 =  (TH1*) total_bkg_pre->Clone("frame");
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

  TH1* mphist = (TH1*) total_bkg_pre->Clone("mphist");
  TH1* mchist = (TH1*) total_bkg_pre->Clone("mchist");
  TH1* unhist = (TH1*) total_bkg_pre->Clone("unhist");
  mchist->Reset();
  unhist->Reset();
  mchist->Add(EMhist);
  mchist->Add(VVVhist);
  mchist->Add(Zjethist);
  mchist->Add(WZhist);
  mchist->Add(ZZhist);

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

  TH1F* band = (TH1F*) total_bkg_post->Clone("band");

  total_bkg_post->Divide(mchist);
  total_bkg_post->SetLineColor(0);
  total_bkg_post->SetMarkerColor(0);
  total_bkg_post->SetMarkerSize(0);
  total_bkg_post->SetFillColor(kGray);

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

 
  total_bkg_post->Draw("E2 SAME");
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

  TH1* frame3 = (TH1*) total_bkg_pre->Clone("frame2");
  frame3->Reset();
  frame3->SetLineColor(kBlack);
  frame3->SetLineWidth(1);
  frame3->GetYaxis()->SetRangeUser(-3,3);
  frame3->GetXaxis()->SetNdivisions(505);
  frame3->GetXaxis()->SetTitle("BDT classifier");
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

  TH1F* data_pull_post = (TH1F*) total_bkg_post->Clone("data_pull_post");
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
  TH1* unhist2 = (TH1*) total_bkg_post->Clone("unhist");
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

