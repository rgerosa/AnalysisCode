#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float luminosity = 35.9;
static int reductionFactor = 1;

void plotDistributions(TH1* histo1, TH1* histo2, TH1* histo3, const string & outputDIR, const string & plot){

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  //                                                                                                                                                                                                  
  TH1* ratio_1 =  (TH1*) histo2->Clone("ratio_1");
  TH1* ratio_2 =  (TH1*) histo3->Clone("ratio_2");
  ratio_1->Divide(histo1);
  ratio_2->Divide(histo1);
  histo1->GetXaxis()->SetLabelSize(0);
  histo1->GetXaxis()->SetTitleSize(0);
  
  ratio_1->GetXaxis()->SetTitle("M_{jj} [GeV]");
  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.1,
				   max(histo1->GetMaximum(),max(histo2->GetMaximum(),histo3->GetMaximum()))*100);
  
  histo1->GetYaxis()->SetTitle("Events");
  histo1->GetYaxis()->SetTitleOffset(1.1);
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2);
  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2);
  histo3->SetLineColor(kBlue);
  histo3->SetLineWidth(2);

  histo1->Draw("hist");
  histo2->Draw("hist same");
  histo3->Draw("hist same");

  ///                                                                                                                                                                                                  
  CMS_lumi(canvas,Form("%.1f",luminosity));
  
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(plot == "fraction"){
    leg->AddEntry(histo1,"W #rightarrow #tau#nu","L");
    leg->AddEntry(histo2,"W #rightarrow #mu#nu","L");
    leg->AddEntry(histo3,"W #rightarrow e#nu","L");
  }
  else if(plot == "muonAcc"){
    leg->AddEntry(histo1,"W #rightarrow #mu#nu","L");
    leg->AddEntry(histo2,"W #rightarrow #mu#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow #mu#nu OUT","L");    
  }
  else if(plot == "electronAcc"){
    leg->AddEntry(histo1,"W #rightarrow e#nu","L");
    leg->AddEntry(histo2,"W #rightarrow e#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow e#nu OUT","L");    
  }
  else if(plot == "tauAcc"){
    leg->AddEntry(histo1,"W #rightarrow #tau#nu","L");
    leg->AddEntry(histo2,"W #rightarrow #tau#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow #tau#nu OUT","L");    
  }
  leg->Draw("same");

  canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_1->GetYaxis()->SetTitle("Ratio");
  ratio_1->GetYaxis()->SetTitleOffset(1.20);
  ratio_1->GetYaxis()->SetTitleSize(0.04);
  ratio_1->GetYaxis()->SetLabelSize(0.03);
  ratio_1->GetYaxis()->SetNdivisions(505);
  ratio_1->GetXaxis()->SetTitleOffset(1.10);
  ratio_1->GetXaxis()->SetNdivisions(505);
  ratio_1->SetLineColor(kRed);
  ratio_1->SetLineWidth(2);
  ratio_2->SetLineColor(kBlue);
  ratio_2->SetLineWidth(2);
  ratio_1->Draw("hist");
  ratio_2->Draw("hist same");
  ratio_1->GetYaxis()->SetRangeUser(min(ratio_1->GetMinimum(),ratio_2->GetMinimum())*0.8,
				    max(ratio_1->GetMaximum(),ratio_2->GetMaximum())*1.2);
  canvas->SaveAs((outputDIR+"/"+plot+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+plot+".pdf").c_str(),"pdf");

  if(canvas) delete canvas;

}



void makeWJetsFractions(string inputDIR, string outputDIR, string kfactorFile, Category category){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  initializeBinning();

  // k-factors
  cout<<"Load k-factors"<<endl;
  TFile* kffile = TFile::Open(kfactorFile.c_str());
  TH1* hist_nloqcdewk = NULL;
  TH1* hist_nloqcd    = NULL;
  TH1* hist_loqcd     = NULL;

  hist_nloqcdewk = (TH1*) kffile->Get("EWKcorr/W");
  hist_nloqcd    = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
  hist_loqcd     = (TH1*) kffile->Get("WJets_LO/inv_pt");

  hist_nloqcdewk->Divide(hist_nloqcd);
  hist_nloqcd->Divide(hist_loqcd);

  TFile* kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets_v2.root");
  ///////////////////                                                                                                                                                                                
  TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
  if(category == Category::VBFrelaxed)
    wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf_relaxed");
  if(category == Category::VBF)
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
  else if(category == Category::VBFrelaxed)
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf_relaxed"));

  TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
  wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));

  wjet_nlo_vbf->Divide(wjet_nlo_mj);

  // VBF k-factor
  vector<TH1*> khists; khists.push_back(hist_nloqcd); khists.push_back(hist_nloqcdewk); khists.push_back(wjet_nlo_vbf);

  // histograms
  cout<<"Book histograms"<<endl;
  vector<double> mjj_bin;
  if(category == Category::VBFrelaxed)
    mjj_bin = selectBinning("mjj",category);
  else if(category == Category::VBF)
    mjj_bin = {1300,2000,3500};

  TH1F* histo_mjj_total = new TH1F("histo_mjj_total","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon  = new TH1F("histo_mjj_muon","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele   = new TH1F("histo_mjj_ele","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau   = new TH1F("histo_mjj_tau","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon_inaccept  = new TH1F("histo_mjj_muon_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon_outaccept = new TH1F("histo_mjj_muon_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_inaccept  = new TH1F("histo_mjj_ele_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_outaccept = new TH1F("histo_mjj_ele_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_inaccept  = new TH1F("histo_mjj_tau_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_outaccept = new TH1F("histo_mjj_tau_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);

  histo_mjj_total->Sumw2();
  histo_mjj_muon->Sumw2();
  histo_mjj_ele->Sumw2();
  histo_mjj_tau->Sumw2();
  histo_mjj_muon_inaccept->Sumw2();
  histo_mjj_ele_inaccept->Sumw2();
  histo_mjj_tau_inaccept->Sumw2();
  histo_mjj_muon_outaccept->Sumw2();
  histo_mjj_ele_outaccept->Sumw2();
  histo_mjj_tau_outaccept->Sumw2();
  
  // tree reader
  TChain* tree_wjet = new TChain("tree/tree");
  tree_wjet->Add((inputDIR+"/*root").c_str());

  TTreeReader reader (tree_wjet);
  TTreeReaderValue<int>    putrue  (reader,"putrue");
  TTreeReaderValue<float>  xsec    (reader,"xsec");
  TTreeReaderValue<float>  wgt     (reader,"wgt");
  TTreeReaderValue<double> wgtsum  (reader,"wgtsum");
  TTreeReaderValue<float>  wgtpu   (reader,"wgtpileup");
  TTreeReaderValue<float>  wgtbtag (reader,"wgtbtag");
  TTreeReaderValue<UChar_t> hltm90      (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90    (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm120   (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<unsigned int> nbjets     (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> ntaus      (reader,"ntausold");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<unsigned int> nincjets  (reader,"njetsinc");
  TTreeReaderValue<float> met         (reader,"t1pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> l1eta   (reader,"l1eta");
  TTreeReaderValue<float> l1pt    (reader,"l1pt");
  TTreeReaderValue<int>   l1id    (reader,"l1id");
  TTreeReaderValue<float> l2eta   (reader,"l2eta");
  TTreeReaderValue<float> l2pt    (reader,"l2pt");
  TTreeReaderValue<int>   l2id    (reader,"l2id");
  TTreeReaderValue<float> wzpt    (reader,"wzpt");


  // Met trigger efficiency                                                                                                                                                                         
  vector<TFile*> triggerfile_MET_binned_Wmn;
  vector<TF1*> triggermet_func_binned_Wmn;
  if(category == Category::VBF or category == Category::VBFrelaxed){ // monojet                                                                                    
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_0.0_800.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_800.0_1200.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1200.0_1700.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1700.0_3000.0.root"));
    if(triggerfile_MET_binned_Wmn.size() != 0){
      for(auto ifile : triggerfile_MET_binned_Wmn)
	triggermet_func_binned_Wmn.push_back((TF1*) ifile->Get("efficiency_func"));
    }
  }

  // loop on events                                                                                                                                                                                   
  cout<<"Event loop"<<endl;
  long int nEvents = tree_wjet->GetEntries();
  long int iEvent = 0;

  while(reader.Next()){

    cout.flush();
    if(iEvent %10000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(double(nEvents)/reductionFactor)<<" % ";
    if(iEvent > double(nEvents)/reductionFactor) break;
    iEvent++;

    // check trigger depending on the sample                                                                                                                                                          
    Double_t hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm90+*hltmwm120+*hltmwm170+*hltmwm300;
    if(hlt == 0) continue;
    ///----
    if(*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fcsc == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    ///----
    if(*nbjets > 0) continue;
    if(*ntaus > 0) continue;
    ///----
    if(fabs(*met-*metcalo)/(*met) > 0.5) continue;
    if(*met < 250) continue;
    ///----
    if(*nincjets < 2) continue;
    if(jetpt->size() < 2) continue;
    if(jetpt->at(0) < 80) continue;
    if(jetpt->at(1) < 40) continue;
    if(fabs(jeteta->at(0)) > 4.7) continue;
    if(fabs(jeteta->at(1)) > 4.7) continue;
    if(jeteta->at(0)*jeteta->at(1) > 0) continue;
    ///----
    if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
    if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
    if(*jmmdphi < 0.5) continue;

    TLorentzVector jet1,jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));

    if(category == Category::VBF){
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 4) continue;
      if((jet1+jet2).M() < 1300) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
    }
    else if(category == Category::VBFrelaxed){
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1) continue;
      if((jet1+jet2).M() < 200) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.3) continue;
    }

    // met trigger scale factor                                                                                                                                                                    
    Double_t trig_wgt = 1.;
    if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed){
      double pfmet = *met;
      if((jet1+jet2).M() < 800)
	trig_wgt *= triggermet_func_binned_Wmn.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(0)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	trig_wgt *= triggermet_func_binned_Wmn.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(1)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	trig_wgt *= triggermet_func_binned_Wmn.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(2)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1700)
	trig_wgt *= triggermet_func_binned_Wmn.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(3)->GetXaxis()->GetXmax()));
    }

    //Gen level info --> NLO re-weight                                                                                                                                                                
    Double_t kwgt = 1.0;
    double genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
        if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
        if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
        kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    double puwgt = 1;
    if(fabs(*wgtpu) > 100)  puwgt = 1;
    else if(fabs(*wgtpu) < 0.01) puwgt = 1;
    else puwgt = *wgtpu;
 
    // Fill histograms
    histo_mjj_total->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum));

    if((fabs(*l1id) == 11 and fabs(*l2id) == 12) or (fabs(*l1id) == 12 and fabs(*l2id) == 11)) // electron decays
      histo_mjj_ele->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum));
    else if((fabs(*l1id) == 13 and fabs(*l2id) == 14) or (fabs(*l1id) == 14 and fabs(*l2id) == 13)) // muon decays
      histo_mjj_muon->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum));
    else if((fabs(*l1id) == 15 and fabs(*l2id) == 16) or (fabs(*l1id) == 15 and fabs(*l2id) == 16)) // tau decays
      histo_mjj_tau->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum));

    if((fabs(*l1id) == 11 and fabs(*l2id) == 12) or (fabs(*l1id) == 12 and fabs(*l2id) == 11)){ // electron decays 
      if(fabs(*l1id) == 11){
	if(*l1pt < 10 or fabs(*l1eta) > 2.5)
	  histo_mjj_ele_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_ele_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
      else if(fabs(*l2id) == 11){
	if(*l2pt < 10 or fabs(*l2eta) > 2.5)
	  histo_mjj_ele_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_ele_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
    }

    if((fabs(*l1id) == 13 and fabs(*l2id) == 14) or (fabs(*l1id) == 14 and fabs(*l2id) == 13)){ // muon decays 
      if(fabs(*l1id) == 13){
	if(*l1pt < 10 or fabs(*l1eta) > 2.4)
	  histo_mjj_muon_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_muon_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
      else if(fabs(*l2id) == 13){
	if(*l2pt < 10 or fabs(*l2eta) > 2.4)
	  histo_mjj_muon_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_muon_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
    }

    if((fabs(*l1id) == 15 and fabs(*l2id) == 16) or (fabs(*l1id) == 16 and fabs(*l2id) == 15)){ // tau decays 
      if(fabs(*l1id) == 15){
	if(*l1pt < 18 or fabs(*l1eta) > 2.3)
	  histo_mjj_tau_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_tau_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
      else if(fabs(*l2id) == 15){
	if(*l2pt < 18 or fabs(*l2eta) > 2.3)
	  histo_mjj_tau_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons
	else
	  histo_mjj_tau_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)/(*wgtsum)); // out of acceptance electrons	
      }
    }
  }
  cout<<endl;

  //// ----- ////
  cout<<"W+jets events in SR "<<histo_mjj_total->Integral()<<endl;
  cout<<"Wmn events in SR "<<histo_mjj_muon->Integral()<<" fraction "<<histo_mjj_muon->Integral()/histo_mjj_total->Integral()<<endl;
  cout<<"Wen events in SR "<<histo_mjj_ele->Integral()<<" fraction "<<histo_mjj_ele->Integral()/histo_mjj_total->Integral()<<endl;
  cout<<"Wtau events in SR "<<histo_mjj_tau->Integral()<<" fraction "<<histo_mjj_tau->Integral()/histo_mjj_total->Integral()<<endl;

  cout<<"Wmn events in SR: in acceptance "<<histo_mjj_muon_inaccept->Integral()<<"  out acceptance "<<histo_mjj_muon_outaccept->Integral()<<" fraction "<<histo_mjj_muon_inaccept->Integral()/(histo_mjj_muon_inaccept->Integral()+histo_mjj_muon_outaccept->Integral())<<endl;
  cout<<"Wen events in SR: in acceptance "<<histo_mjj_ele_inaccept->Integral()<<"  out acceptance "<<histo_mjj_ele_outaccept->Integral()<<" fraction "<<histo_mjj_ele_inaccept->Integral()/(histo_mjj_ele_inaccept->Integral()+histo_mjj_ele_outaccept->Integral())<<endl;
  cout<<"Wtau events in SR: in acceptance "<<histo_mjj_tau_inaccept->Integral()<<"  out acceptance "<<histo_mjj_tau_outaccept->Integral()<<" fraction "<<histo_mjj_tau_inaccept->Integral()/(histo_mjj_tau_inaccept->Integral()+histo_mjj_tau_outaccept->Integral())<<endl;

  ////// ------------ //////
  plotDistributions(histo_mjj_tau,histo_mjj_muon,histo_mjj_ele,outputDIR,"fraction");
  plotDistributions(histo_mjj_muon,histo_mjj_muon_inaccept,histo_mjj_muon_outaccept,outputDIR,"muonAcc");
  plotDistributions(histo_mjj_ele,histo_mjj_ele_inaccept,histo_mjj_ele_outaccept,outputDIR,"eleAcc");
  plotDistributions(histo_mjj_tau,histo_mjj_tau_inaccept,histo_mjj_tau_outaccept,outputDIR,"tauAcc");  
  

}
