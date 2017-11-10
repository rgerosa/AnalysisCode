#include "../CMS_lumi.h"

void plotHistogram(TCanvas* canvas,
		   const vector<TH1F*> & histo_1,
		   const vector<TH1F*> & histo_2,
		   const string & outputDIR,
		   const string & observable,
		   const string & leg_title,
		   const string & postfix,
		   const pair<string,string> label){


  TPad* pad2 = new TPad("pad2","",0,0,1.,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  TH1F* histo_total_1 = (TH1F*) histo_1.at(0)->Clone(Form("%s_total",histo_1.at(0)->GetName()));
  TH1F* histo_total_2 = (TH1F*) histo_2.at(0)->Clone(Form("%s_total",histo_2.at(0)->GetName()));

  TH1F* histo_error_1_up = (TH1F*) histo_1.at(0)->Clone(Form("%s_error_up",histo_1.at(0)->GetName()));
  TH1F* histo_error_2_up = (TH1F*) histo_2.at(0)->Clone(Form("%s_error_up",histo_2.at(0)->GetName()));
  histo_error_1_up->Reset();
  histo_error_2_up->Reset();
  TH1F* histo_error_1_dw = (TH1F*) histo_1.at(0)->Clone(Form("%s_error_dw",histo_1.at(0)->GetName()));
  TH1F* histo_error_2_dw = (TH1F*) histo_2.at(0)->Clone(Form("%s_error_dw",histo_2.at(0)->GetName()));
  histo_error_1_dw->Reset();
  histo_error_2_dw->Reset();

  for(int ibin = 0; ibin < histo_error_1_up->GetNbinsX(); ibin++){
    histo_error_1_up->SetBinContent(ibin+1,(histo_1.at(1)->GetBinContent(ibin+1)-histo_1.at(0)->GetBinContent(ibin+1))/(histo_total_1->GetBinContent(ibin+1)));
    histo_error_2_up->SetBinContent(ibin+1,(histo_2.at(1)->GetBinContent(ibin+1)-histo_2.at(0)->GetBinContent(ibin+1))/(histo_total_2->GetBinContent(ibin+1)));
    histo_error_1_dw->SetBinContent(ibin+1,(histo_1.at(2)->GetBinContent(ibin+1)-histo_1.at(0)->GetBinContent(ibin+1))/(histo_total_1->GetBinContent(ibin+1)));
    histo_error_2_dw->SetBinContent(ibin+1,(histo_2.at(2)->GetBinContent(ibin+1)-histo_2.at(0)->GetBinContent(ibin+1))/(histo_total_2->GetBinContent(ibin+1)));
  }

  histo_total_1->SetLineColor(kBlue);
  histo_total_1->SetLineWidth(2);

  histo_total_2->SetLineColor(kRed);
  histo_total_2->SetLineWidth(3);
  histo_total_2->SetLineStyle(7);

  histo_total_1->GetXaxis()->SetTitleSize(0);
  histo_total_1->GetXaxis()->SetLabelSize(0);
  histo_total_1->GetYaxis()->SetTitle("Entries");
  histo_total_1->GetYaxis()->SetRangeUser(min(histo_total_1->GetMinimum(),histo_total_2->GetMinimum())*0.1,max(histo_total_1->GetMaximum(),histo_total_2->GetMaximum())*10);

  histo_total_1->Draw("hist");
  histo_total_2->Draw("hist same");

  TLegend leg (0.6,0.68,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry((TObject*)(0),leg_title.c_str(),"");
  leg.AddEntry(histo_total_1,label.first.c_str(),"FL");
  leg.AddEntry(histo_total_2,label.second.c_str(),"FL");
  leg.Draw("same");
  
  CMS_lumi(canvas,"35.9");

  canvas->cd();
  pad2->Draw();
  pad2->cd();

  histo_error_1_up->SetLineColor(kBlue);
  histo_error_1_up->SetLineWidth(2);
  histo_error_1_dw->SetLineColor(kBlue);
  histo_error_1_dw->SetLineWidth(2);

  histo_error_2_up->SetLineColor(kRed);
  histo_error_2_up->SetLineWidth(2);
  histo_error_2_up->SetLineStyle(7);
  histo_error_2_dw->SetLineColor(kRed);
  histo_error_2_dw->SetLineWidth(2);
  histo_error_2_dw->SetLineStyle(7);

  histo_error_1_up->GetYaxis()->SetTitle("Uncertainty");

  histo_error_1_up->GetYaxis()->SetTitleOffset(1.45);
  histo_error_1_up->GetXaxis()->SetTitleOffset(1.1);
  histo_error_1_up->GetYaxis()->SetLabelSize(0.03);
  histo_error_1_up->GetYaxis()->SetTitleSize(0.04);
  histo_error_1_up->GetXaxis()->SetLabelSize(0.04);
  histo_error_1_up->GetXaxis()->SetTitleSize(0.05);
  histo_error_1_up->GetXaxis()->SetNdivisions(505);

  if(observable == "mjj")
    histo_error_1_up->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "detajj")
    histo_error_1_up->GetXaxis()->SetTitle("#Delta#eta_{jj}");
  else if(observable == "met_onebin")
    histo_error_1_up->GetXaxis()->SetTitle("Recoil [GeV]");

  if(TString(postfix).Contains("QCD") and TString(postfix).Contains("RenScale2")) histo_error_1_up->GetYaxis()->SetRangeUser(-0.15,0.15);
  else if(TString(postfix).Contains("EWK") and TString(postfix).Contains("RenScale2")) histo_error_1_up->GetYaxis()->SetRangeUser(-0.12,0.12);
  else histo_error_1_up->GetYaxis()->SetRangeUser(-0.08,0.08);

  histo_error_1_up->GetYaxis()->SetNdivisions(505);
  histo_error_1_up->Draw("hist");
  histo_error_1_dw->Draw("hist same");
  histo_error_2_up->Draw("hist same");
  histo_error_2_dw->Draw("hist same");
  
  canvas->cd();
  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  if(pad2) delete pad2;

}

////
void makeHistogramList(vector<TH1F*> & histogram_QCD, // variation of QCD nuisance 
		       vector<TH1F*> & histogram_EWK,  // variation of EWK nuisance
		       RooWorkspace* w, 
		       const string & controlRegion, 
		       const string & observable, 
		       const string & systematic){

  // take bin contents                                                                                                                                                                                 
  vector<RooFormulaVar*> WJets_QCD;
  vector<RooFormulaVar*> WJets_EWK;

  RooRealVar* var_1 = (RooRealVar*) w->obj(Form("%s_VBF",observable.c_str()));
  const RooAbsBinning & binning = var_1->getBinning();

  for(int ibin = 0; ibin < binning.numBins(); ibin++){
    if((RooFormulaVar*) w->obj(Form("WJets_%s_VBF_bin%d",controlRegion.c_str(),ibin+1)))
      WJets_QCD.push_back((RooFormulaVar*) w->obj(Form("WJets_%s_VBF_bin%d",controlRegion.c_str(),ibin+1)));
    if((RooFormulaVar*) w->obj(Form("WJets_EWK_%s_VBF_bin%d",controlRegion.c_str(),ibin+1)))
      WJets_EWK.push_back((RooFormulaVar*) w->obj(Form("WJets_EWK_%s_VBF_bin%d",controlRegion.c_str(),ibin+1)));
  }
  
  // Nuisances
  RooRealVar* ZW_QCD_sys = (RooRealVar*) w->obj(Form("ZW_QCD_SR_%s",systematic.c_str()));
  RooRealVar* ZW_EWK_sys = (RooRealVar*) w->obj(Form("ZW_EWK_SR_%s",systematic.c_str()));

  // histograms
  vector<TH1F*> histo_WJets_QCD ;
  vector<TH1F*> histo_WJets_EWK ;
  
  if(WJets_QCD.size() != 0 and WJets_EWK.size() != 0){ // split case
  
    for(int ierr = 0; ierr < 3; ierr++){

      histo_WJets_QCD.push_back(new TH1F(Form("histo_WJets_%s_QCD_%d",controlRegion.c_str(),ierr),"",binning.numBins(),binning.array()));
      histo_WJets_EWK.push_back(new TH1F(Form("histo_WJets_%s_EWK_%d",controlRegion.c_str(),ierr),"",binning.numBins(),binning.array()));      
      histo_WJets_QCD.back()->Sumw2();
      histo_WJets_EWK.back()->Sumw2();

      //////////////////
      if(ierr == 0){
	ZW_QCD_sys->setVal(0);
	ZW_EWK_sys->setVal(0);
      }
      else if(ierr == 1){
	ZW_QCD_sys->setVal(1);
	ZW_EWK_sys->setVal(1);
      }
      else if(ierr == 2){
	ZW_QCD_sys->setVal(-1);
	ZW_EWK_sys->setVal(-1);
      }

      for(int ibin = 0; ibin < binning.numBins(); ibin++){// set bin content      
	histo_WJets_QCD.back()->SetBinContent(ibin+1,WJets_QCD.at(ibin)->getVal());
	histo_WJets_EWK.back()->SetBinContent(ibin+1,WJets_EWK.at(ibin)->getVal());
      }

      // QCD histograms --> central value for both QCD and EWK
      histogram_QCD.push_back((TH1F*) histo_WJets_QCD.back()->Clone(Form("WJets_QCD_%s_split_%d",controlRegion.c_str(),ierr)));
      histogram_QCD.back()->Add(histo_WJets_EWK.at(0));

      histogram_EWK.push_back((TH1F*) histo_WJets_QCD.at(0)->Clone(Form("WJets_EWK_%s_split_%d",controlRegion.c_str(),ierr)));
      histogram_EWK.back()->Add(histo_WJets_EWK.back());
      
    } 
  }
  else{ // Merged case

    for(int ierr = 0; ierr < 3; ierr++){      
      histogram_QCD.push_back(new TH1F(Form("WJets_QCD_%s_merged_%d",controlRegion.c_str(),ierr),"",binning.numBins(),binning.array()));
      histogram_EWK.push_back(new TH1F(Form("WJets_EWK_%s_merged_%d",controlRegion.c_str(),ierr),"",binning.numBins(),binning.array()));
      
      histogram_QCD.back()->Sumw2();
      histogram_EWK.back()->Sumw2();
              
      for(int ibin = 0; ibin < binning.numBins(); ibin++){// set bin content
	
	//////////////////
	ZW_EWK_sys->setVal(0);
	if(ierr == 0)
	  ZW_QCD_sys->setVal(0);
	else if(ierr == 1)
	  ZW_QCD_sys->setVal(1);
	else if(ierr == 2)
	  ZW_QCD_sys->setVal(-1);
	
	histogram_QCD.back()->SetBinContent(ibin+1,WJets_QCD.at(ibin)->getVal());
	
	//////////////////
	ZW_QCD_sys->setVal(0);
	if(ierr == 0)
	  ZW_EWK_sys->setVal(0);
	else if(ierr == 1)
	  ZW_EWK_sys->setVal(1);
	else if(ierr == 2)
	  ZW_EWK_sys->setVal(-1);
	
	histogram_EWK.back()->SetBinContent(ibin+1,WJets_QCD.at(ibin)->getVal());
	
      }
    }    
  }
}

void makeWorkspaceComparisonVBF(string fileName1, string fileName2, string outputDIR, string observable, string systematic, pair<string,string> labels){

  if(observable != "mjj" and observable != "met_onebin" and observable != "detajj"){
    std::cerr<<"Unrecognized variable --> exit "<<std::endl;
    return;
  }

  if(systematic != "RenScale2" and systematic != "FactScale2"){
    std::cerr<<"Unrecognized systematic --> exit "<<std::endl;
    return;
  }

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* file1 = TFile::Open(fileName1.c_str(),"READ");
  TFile* file2  = TFile::Open(fileName2.c_str(),"READ");

  // splitting case --> take workspaces
  RooWorkspace* SR_VBF_1 = (RooWorkspace*) file1->Get("SR_VBF");
  RooWorkspace* WM_VBF_1 = (RooWorkspace*) file1->Get("WM_VBF");
  RooWorkspace* WE_VBF_1 = (RooWorkspace*) file1->Get("WE_VBF");

  /// final histograms used in the comparison
  vector<TH1F*> histo_WJets_SR_QCD_1 ;
  vector<TH1F*> histo_WJets_WM_QCD_1 ;
  vector<TH1F*> histo_WJets_WE_QCD_1 ;

  vector<TH1F*> histo_WJets_SR_EWK_1 ;
  vector<TH1F*> histo_WJets_WM_EWK_1 ;
  vector<TH1F*> histo_WJets_WE_EWK_1 ;

  makeHistogramList(histo_WJets_SR_QCD_1,histo_WJets_SR_EWK_1,SR_VBF_1,"SR",observable,systematic);
  makeHistogramList(histo_WJets_WM_QCD_1,histo_WJets_WM_EWK_1,WM_VBF_1,"WM",observable,systematic);
  makeHistogramList(histo_WJets_WE_QCD_1,histo_WJets_WE_EWK_1,WE_VBF_1,"WE",observable,systematic);

  /////////////////////
  RooWorkspace* SR_VBF_2 = (RooWorkspace*) file2->Get("SR_VBF");
  RooWorkspace* WM_VBF_2 = (RooWorkspace*) file2->Get("WM_VBF");
  RooWorkspace* WE_VBF_2 = (RooWorkspace*) file2->Get("WE_VBF");

  vector<TH1F*> histo_WJets_SR_QCD_2 ;
  vector<TH1F*> histo_WJets_WM_QCD_2 ;
  vector<TH1F*> histo_WJets_WE_QCD_2 ;

  vector<TH1F*> histo_WJets_SR_EWK_2 ;
  vector<TH1F*> histo_WJets_WM_EWK_2 ;
  vector<TH1F*> histo_WJets_WE_EWK_2 ;

  makeHistogramList(histo_WJets_SR_QCD_2,histo_WJets_SR_EWK_2,SR_VBF_2,"SR",observable,systematic);
  makeHistogramList(histo_WJets_WM_QCD_2,histo_WJets_WM_EWK_2,WM_VBF_2,"WM",observable,systematic);
  makeHistogramList(histo_WJets_WE_QCD_2,histo_WJets_WE_EWK_2,WE_VBF_2,"WE",observable,systematic);

  /// plotting final histograms
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->SetBottomMargin(0.3);
  canvas->cd();

  plotHistogram(canvas,histo_WJets_SR_QCD_1,histo_WJets_SR_QCD_2,outputDIR,observable,"W #rightarrow l#nu (SR)","WJet_SR_QCD_"+systematic,labels);
  plotHistogram(canvas,histo_WJets_WM_QCD_1,histo_WJets_WM_QCD_2,outputDIR,observable,"W #rightarrow #mu#nu (WM)","WJet_WM_QCD_"+systematic,labels);
  plotHistogram(canvas,histo_WJets_WE_QCD_1,histo_WJets_WE_QCD_2,outputDIR,observable,"W #rightarrow e#nu (WE)","WJet_WE_QCD_"+systematic,labels);
  
  plotHistogram(canvas,histo_WJets_SR_EWK_1,histo_WJets_SR_EWK_2,outputDIR,observable,"W #rightarrow l#nu (SR)","WJet_SR_EWK_"+systematic,labels);
  plotHistogram(canvas,histo_WJets_WM_EWK_1,histo_WJets_WM_EWK_2,outputDIR,observable,"W #rightarrow #mu#nu (WM)","WJet_WM_EWK_"+systematic,labels);
  plotHistogram(canvas,histo_WJets_WE_EWK_1,histo_WJets_WE_EWK_2,outputDIR,observable,"W #rightarrow e#nu (WE)","WJet_WE_EWK_"+systematic,labels);

  TFile* outputFile = new TFile((outputDIR+"/template_file.root").c_str(),"RECREATE");
  outputFile->cd();

  for(auto hist : histo_WJets_SR_QCD_1) hist->Write();
  for(auto hist : histo_WJets_SR_QCD_2) hist->Write();
  for(auto hist : histo_WJets_WM_QCD_1) hist->Write();
  for(auto hist : histo_WJets_WM_QCD_2) hist->Write();
  for(auto hist : histo_WJets_WE_QCD_1) hist->Write();
  for(auto hist : histo_WJets_WE_QCD_2) hist->Write();

  for(auto hist : histo_WJets_SR_EWK_1) hist->Write();
  for(auto hist : histo_WJets_SR_EWK_2) hist->Write();
  for(auto hist : histo_WJets_WM_EWK_1) hist->Write();
  for(auto hist : histo_WJets_WM_EWK_2) hist->Write();
  for(auto hist : histo_WJets_WE_EWK_1) hist->Write();
  for(auto hist : histo_WJets_WE_EWK_2) hist->Write();
  
  outputFile->Close();
}
