#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void plotComparison(TCanvas* canvas, TH1F* tf_prefit, TH1F* tf_postfit,const string & postfix, const string & observable,const string & outputDIR){

  canvas->cd();
  tf_prefit->SetLineColor(kBlack);
  tf_prefit->SetLineWidth(3);
  tf_prefit->SetMarkerColor(kBlack);
  tf_prefit->SetMarkerStyle(20);
  tf_prefit->SetMarkerSize(1);

  tf_postfit->SetLineColor(kRed);
  tf_postfit->SetLineWidth(3);
  tf_postfit->SetLineStyle(7);
  tf_postfit->SetMarkerColor(kRed);
  tf_postfit->SetMarkerStyle(24);
  tf_postfit->SetMarkerSize(1);

  tf_prefit->GetXaxis()->SetTitle(observable.c_str());
  tf_prefit->GetXaxis()->SetTitleOffset(1.10);

  if(postfix == "rzmm")
    tf_prefit->GetYaxis()->SetTitle("R(Z#mu#mu)");
  else if(postfix == "rzee")
    tf_prefit->GetYaxis()->SetTitle("R(Zee)");
  else if(postfix == "rwmn")
    tf_prefit->GetYaxis()->SetTitle("R(W#mu#nu)");
  else if(postfix == "rwen")
    tf_prefit->GetYaxis()->SetTitle("R(We#nu)");
  else if(postfix == "rzw")
    tf_prefit->GetYaxis()->SetTitle("R(ZW)");
  else if(postfix == "rzmm_ewk")
    tf_prefit->GetYaxis()->SetTitle("R(Z#mu#mu EW)");
  else if(postfix == "rzee_ewk")
    tf_prefit->GetYaxis()->SetTitle("R(Zee EW)");
  else if(postfix == "rwmn_ewk")
    tf_prefit->GetYaxis()->SetTitle("R(W#mu#nu EW)");
  else if(postfix == "rwen_ewk")
    tf_prefit->GetYaxis()->SetTitle("R(We#nu EW)");
  else if(postfix == "rzw_ewk")
    tf_prefit->GetYaxis()->SetTitle("R(ZW EW)");
  
  tf_prefit->GetYaxis()->SetRangeUser(tf_prefit->GetMinimum()*0.75,tf_prefit->GetMaximum()*1.25);
  
  tf_prefit->Draw("hist");
  tf_postfit->Draw("hist same");

  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(tf_prefit,"Pre-fit TF","L");
  leg.AddEntry(tf_postfit,"Post-fit TF","L");
  leg.Draw("same");

  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}

void makeTransferFactorPreFitVsPostFit(string inputFileName, string outputDIR, Category category){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  
  string sr_folder;
  string zmm_folder;
  string zee_folder;
  string wen_folder;
  string wmn_folder;
  string gam_folder;

  if(category == Category::monojet or category == Category::monoV){
    sr_folder  = "ch1";
    zmm_folder = "ch2";
    zee_folder = "ch5";
    wen_folder = "ch6";
    wmn_folder = "ch3";
    gam_folder = "ch4";
  }
  else if(category == Category::VBF or category == Category::VBFrelaxed){
    //    sr_folder  = "ch1";
    //    zmm_folder = "ch2";
    //    zee_folder = "ch4";
    //    wen_folder = "ch5";
    //    wmn_folder = "ch3";
    sr_folder  = "vbf_signal";
    zmm_folder = "vbf_dimuon";
    zee_folder = "vbf_dielec";
    wmn_folder = "vbf_singlemu";
    wen_folder = "vbf_singleel";
  }

  // pre-fit histos
  TH1F* znn_qcd_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+sr_folder+"/qcd_znunu").c_str());
  TH1F* wjet_qcd_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+sr_folder+"/qcd_wjets").c_str());
  TH1F* zmm_qcd_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+zmm_folder+"/qcd_zll").c_str());
  TH1F* zee_qcd_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+zee_folder+"/qcd_zll").c_str());
  TH1F* wmn_qcd_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+wmn_folder+"/qcd_wjets").c_str());
  TH1F* wen_qcd_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+wen_folder+"/qcd_wjets").c_str());

  // post-fit histos
  TH1F* znn_qcd_postfit  = (TH1F*) inputFile->Get(("shapes_fit_b/"+sr_folder+"/qcd_znunu").c_str());
  TH1F* wjet_qcd_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+sr_folder+"/qcd_wjets").c_str());
  TH1F* zmm_qcd_postfit  = (TH1F*) inputFile->Get(("shapes_fit_b/"+zmm_folder+"/qcd_zll").c_str());
  TH1F* zee_qcd_postfit  = (TH1F*) inputFile->Get(("shapes_fit_b/"+zee_folder+"/qcd_zll").c_str());
  TH1F* wmn_qcd_postfit  = (TH1F*) inputFile->Get(("shapes_fit_b/"+wmn_folder+"/qcd_wjets").c_str());
  TH1F* wen_qcd_postfit  = (TH1F*) inputFile->Get(("shapes_fit_b/"+wen_folder+"/qcd_wjets").c_str());

  TH1F* gamma_prefit  = NULL;
  TH1F* gamma_postfit = NULL;

  if(category == Category::monojet or category == Category::monoV){
    gamma_prefit  = (TH1F*) inputFile->Get(("shapes_prefit/"+gam_folder+"/gjets").c_str());
    gamma_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+gam_folder+"/gjets").c_str());
  }



  TH1F* znn_ewk_prefit =  NULL;
  TH1F* wjet_ewk_prefit =  NULL;
  TH1F* zmm_ewk_prefit =  NULL;
  TH1F* zee_ewk_prefit =  NULL;
  TH1F* wmn_ewk_prefit =  NULL;
  TH1F* wen_ewk_prefit =  NULL;

  TH1F* znn_ewk_postfit =  NULL;
  TH1F* wjet_ewk_postfit =  NULL;
  TH1F* zmm_ewk_postfit =  NULL;
  TH1F* zee_ewk_postfit =  NULL;
  TH1F* wmn_ewk_postfit =  NULL;
  TH1F* wen_ewk_postfit =  NULL;

  if(category == Category::VBF or category == Category::VBFrelaxed){
    znn_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+sr_folder+"/ewk_znunu").c_str());
    wjet_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+sr_folder+"/ewk_wjets").c_str());
    zmm_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+zmm_folder+"/ewk_zll").c_str());
    zee_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+zee_folder+"/ewk_zll").c_str());
    wmn_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+wmn_folder+"/ewk_wjets").c_str());
    wen_ewk_prefit = (TH1F*) inputFile->Get(("shapes_prefit/"+wen_folder+"/ewk_wjets").c_str());

    znn_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+sr_folder+"/ewk_znunu").c_str());
    wjet_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+sr_folder+"/ewk_wjets").c_str());
    zmm_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+zmm_folder+"/ewk_zll").c_str());
    zee_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+zee_folder+"/ewk_zll").c_str());
    wmn_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+wmn_folder+"/ewk_wjets").c_str());
    wen_ewk_postfit = (TH1F*) inputFile->Get(("shapes_fit_b/"+wen_folder+"/ewk_wjets").c_str());
  }

  // transfer factors
  TH1F* tf_zmm_qcd_prefit = (TH1F*) znn_qcd_prefit->Clone("tf_zmm_qcd_prefit");
  tf_zmm_qcd_prefit->Divide(zmm_qcd_prefit);
  TH1F* tf_zmm_qcd_postfit = (TH1F*) znn_qcd_postfit->Clone("tf_zmm_qcd_postfit");
  tf_zmm_qcd_postfit->Divide(zmm_qcd_postfit);

  TH1F* tf_zee_qcd_prefit = (TH1F*) znn_qcd_prefit->Clone("tf_zee_qcd_prefit");
  tf_zee_qcd_prefit->Divide(zee_qcd_prefit);
  TH1F* tf_zee_qcd_postfit = (TH1F*) znn_qcd_postfit->Clone("tf_zee_qcd_postfit");
  tf_zee_qcd_postfit->Divide(zee_qcd_postfit);

  TH1F* tf_wmn_qcd_prefit = (TH1F*) wjet_qcd_prefit->Clone("tf_wmn_qcd_prefit");
  tf_wmn_qcd_prefit->Divide(wmn_qcd_prefit);
  TH1F* tf_wmn_qcd_postfit = (TH1F*) wjet_qcd_postfit->Clone("tf_wmn_qcd_postfit");
  tf_wmn_qcd_postfit->Divide(wmn_qcd_postfit);

  TH1F* tf_wen_qcd_prefit = (TH1F*) wjet_qcd_prefit->Clone("tf_wen_qcd_prefit");
  tf_wen_qcd_prefit->Divide(wen_qcd_prefit);
  TH1F* tf_wen_qcd_postfit = (TH1F*) wjet_qcd_postfit->Clone("tf_wen_qcd_postfit");
  tf_wen_qcd_postfit->Divide(wen_qcd_postfit);
  TH1F* tf_zw_qcd_prefit = (TH1F*) znn_qcd_prefit->Clone("tf_zw_qcd_prefit");
  tf_zw_qcd_prefit->Divide(wjet_qcd_prefit);
  TH1F* tf_zw_qcd_postfit = (TH1F*) znn_qcd_postfit->Clone("tf_zw_qcd_postfit");
  tf_zw_qcd_postfit->Divide(wjet_qcd_postfit);

  TH1F* tf_zg_prefit = NULL;
  TH1F* tf_zg_postfit = NULL;

  if(category == Category::monojet or category == Category::monoV){
    tf_zg_prefit = (TH1F*) znn_qcd_prefit->Clone("tf_zg_prefit");
    tf_zg_prefit->Divide(gamma_prefit);
    tf_zg_postfit = (TH1F*) znn_qcd_postfit->Clone("tf_zg_postfit");
    tf_zg_postfit->Divide(gamma_postfit);
  }

  TH1F* tf_zmm_ewk_prefit = NULL;
  TH1F* tf_zee_ewk_prefit = NULL;
  TH1F* tf_wmn_ewk_prefit = NULL;
  TH1F* tf_wen_ewk_prefit = NULL;
  TH1F* tf_zw_ewk_prefit = NULL;

  TH1F* tf_zmm_ewk_postfit = NULL;
  TH1F* tf_zee_ewk_postfit = NULL;
  TH1F* tf_wmn_ewk_postfit = NULL;
  TH1F* tf_wen_ewk_postfit = NULL;
  TH1F* tf_zw_ewk_postfit = NULL;

  if(category == Category::VBF or category == Category::VBFrelaxed){

    tf_zmm_ewk_prefit = (TH1F*) znn_ewk_prefit->Clone("tf_zmm_ewk_prefit");
    tf_zmm_ewk_prefit->Divide(zmm_ewk_prefit);
    tf_zmm_ewk_postfit = (TH1F*) znn_ewk_postfit->Clone("tf_zmm_ewk_postfit");
    tf_zmm_ewk_postfit->Divide(zmm_ewk_postfit);

    tf_zee_ewk_prefit = (TH1F*) znn_ewk_prefit->Clone("tf_zee_ewk_prefit");
    tf_zee_ewk_prefit->Divide(zee_ewk_prefit);
    tf_zee_ewk_postfit = (TH1F*) znn_ewk_postfit->Clone("tf_zee_ewk_postfit");
    tf_zee_ewk_postfit->Divide(zee_ewk_postfit);
    
    tf_wmn_ewk_prefit = (TH1F*) wjet_ewk_prefit->Clone("tf_wmn_ewk_prefit");
    tf_wmn_ewk_prefit->Divide(wmn_ewk_prefit);
    tf_wmn_ewk_postfit = (TH1F*) wjet_ewk_postfit->Clone("tf_wmn_ewk_postfit");
    tf_wmn_ewk_postfit->Divide(wmn_ewk_postfit);
    
    tf_wen_ewk_prefit = (TH1F*) wjet_ewk_prefit->Clone("tf_wen_ewk_prefit");
    tf_wen_ewk_prefit->Divide(wen_ewk_prefit);
    tf_wen_ewk_postfit = (TH1F*) wjet_ewk_postfit->Clone("tf_wen_ewk_postfit");
    tf_wen_ewk_postfit->Divide(wen_ewk_postfit);
    
    tf_zw_ewk_prefit = (TH1F*) znn_ewk_prefit->Clone("tf_zw_ewk_prefit");
    tf_zw_ewk_prefit->Divide(wjet_ewk_prefit);
    tf_zw_ewk_postfit = (TH1F*) znn_ewk_postfit->Clone("tf_zw_ewk_postfit");
    tf_zw_ewk_postfit->Divide(wjet_ewk_postfit);
  }

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  
  string observable = "Recoil [GeV]";
  if(category == Category::VBF or category == Category::VBFrelaxed)
    observable = "M_{jj} [GeV]";

  plotComparison(canvas,tf_zmm_qcd_prefit,tf_zmm_qcd_postfit,"rzmm",observable,outputDIR);
  plotComparison(canvas,tf_zee_qcd_prefit,tf_zee_qcd_postfit,"rzee",observable,outputDIR);
  plotComparison(canvas,tf_wmn_qcd_prefit,tf_wmn_qcd_postfit,"rwmn",observable,outputDIR);
  plotComparison(canvas,tf_wen_qcd_prefit,tf_wen_qcd_postfit,"rwen",observable,outputDIR);
  plotComparison(canvas,tf_zw_qcd_prefit,tf_zw_qcd_postfit,"rzw",observable,outputDIR);
  
  if(category == Category::monojet or category == Category::monoV)
    plotComparison(canvas,tf_zg_prefit,tf_zg_postfit,"rgam",observable,outputDIR);

  if(category == Category::VBF or category == Category::VBFrelaxed){
    plotComparison(canvas,tf_zmm_ewk_prefit,tf_zmm_ewk_postfit,"rzmm_ewk",observable,outputDIR);
    plotComparison(canvas,tf_zee_ewk_prefit,tf_zee_ewk_postfit,"rzee_ewk",observable,outputDIR);
    plotComparison(canvas,tf_wmn_ewk_prefit,tf_wmn_ewk_postfit,"rwmn_ewk",observable,outputDIR);
    plotComparison(canvas,tf_wen_ewk_prefit,tf_wen_ewk_postfit,"rwen_ewk",observable,outputDIR);
    plotComparison(canvas,tf_zw_ewk_prefit,tf_zw_ewk_postfit,"rzw_ewk",observable,outputDIR);
  }
    
}
