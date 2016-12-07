#include "../CMS_lumi.h"

void drawPlot(TCanvas* canvas, TH1* sig1, TH1* sig2, TH1* sig3, TH1* bkg, string outputDIR, string xAxisLabel){

  sig1->Scale(1./sig1->Integral());
  sig2->Scale(1./sig2->Integral());
  sig3->Scale(1./sig3->Integral());
  bkg->Scale(1./bkg->Integral());

  sig1->GetYaxis()->SetTitleSize(0.050);
  sig1->GetXaxis()->SetTitle(xAxisLabel.c_str());
  sig1->GetYaxis()->SetTitle("a.u.");

  sig1->SetLineColor(kBlack);
  sig1->SetLineWidth(2);
  sig1->Draw("hist");
  bkg->SetLineColor(kBlack);
  bkg->SetFillStyle(1001);
  bkg->SetFillColor(kGray);
  bkg->Draw("HIST same");
  sig1->Draw("HIST same");
  
  sig2->SetLineColor(kBlue);
  sig2->SetLineWidth(2);
  sig3->SetLineColor(kRed);
  sig3->SetLineWidth(2);
  sig2->Draw("HIST same");
  sig3->Draw("HIST same");

  sig1->GetYaxis()->SetRangeUser(0,max(sig1->GetMaximum(),max(sig3->GetMaximum(),bkg->GetMaximum()))*1.25);
  CMS_lumi(canvas,"");
  
  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.AddEntry(sig1,"Mono-W m_{med} = 150 GeV","L");
  leg.AddEntry(sig2,"Mono-W m_{med} = 400 GeV","L");
  leg.AddEntry(sig3,"Mono-W m_{med} = 600 GeV","L");
  leg.AddEntry(bkg,"Z #rightarrow #nu#nu","F");

  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/"+string(sig1->GetName())+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+string(sig1->GetName())+".png").c_str(),"png");

}

void makeBoostedVariablesComparison(string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  
  TChain* signal_mass1 = new TChain("tree/tree","tree/tree");
  TChain* signal_mass2 = new TChain("tree/tree","tree/tree");;
  TChain* signal_mass3 = new TChain("tree/tree","tree/tree");;
  TChain* background = new TChain("tree/tree","tree/tree");

  signal_mass1->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WminusH_HToInvisible_WToQQ_M150_13TeV_powheg_pythia8.root");
  signal_mass1->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WplusH_HToInvisible_WToQQ_M150_13TeV_powheg_pythia8.root");
  signal_mass2->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WminusH_HToInvisible_WToQQ_M400_13TeV_powheg_pythia8.root");
  signal_mass2->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WplusH_HToInvisible_WToQQ_M400_13TeV_powheg_pythia8.root");
  signal_mass3->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WminusH_HToInvisible_WToQQ_M600_13TeV_powheg_pythia8.root");
  signal_mass3->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/HiggsInvisible/sigfilter/sig_WplusH_HToInvisible_WToQQ_M600_13TeV_powheg_pythia8.root");
  background->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_28_11_2016/ZJets/sigfilter/*root");


  TH1F* prunedmass_sig_mass1 = new TH1F("prunedmass_sig_mass1","",30,0,150);
  TH1F* prunedmass_sig_mass2 = new TH1F("prunedmass_sig_mass2","",30,0,150);
  TH1F* prunedmass_sig_mass3 = new TH1F("prunedmass_sig_mass3","",30,0,150);
  TH1F* prunedmass_bkg       = new TH1F("prunedmass_bkg","",30,0,150);
  prunedmass_sig_mass1->Sumw2();
  prunedmass_sig_mass2->Sumw2();
  prunedmass_sig_mass3->Sumw2();
  prunedmass_bkg->Sumw2();

  TH1F* softdropmass_sig_mass1 = new TH1F("softdropmass_sig_mass1","",30,0,150);
  TH1F* softdropmass_sig_mass2 = new TH1F("softdropmass_sig_mass2","",30,0,150);
  TH1F* softdropmass_sig_mass3 = new TH1F("softdropmass_sig_mass3","",30,0,150);
  TH1F* softdropmass_bkg       = new TH1F("softdropmass_bkg","",30,0,150);
  softdropmass_sig_mass1->Sumw2();
  softdropmass_sig_mass2->Sumw2();
  softdropmass_sig_mass3->Sumw2();
  softdropmass_bkg->Sumw2();

  TH1F* tau2tau1_sig_mass1 = new TH1F("tau2tau1_sig_mass1","",30,0.1,1);
  TH1F* tau2tau1_sig_mass2 = new TH1F("tau2tau1_sig_mass2","",30,0.1,1);
  TH1F* tau2tau1_sig_mass3 = new TH1F("tau2tau1_sig_mass3","",30,0.1,1);
  TH1F* tau2tau1_bkg       = new TH1F("tau2tau1_bkg","",30,0.1,1);
  tau2tau1_sig_mass1->Sumw2();
  tau2tau1_sig_mass2->Sumw2();
  tau2tau1_sig_mass3->Sumw2();
  tau2tau1_bkg->Sumw2();

  TH1F* tau2tau1puppi_sig_mass1 = new TH1F("tau2tau1puppi_sig_mass1","",30,0.1,1);
  TH1F* tau2tau1puppi_sig_mass2 = new TH1F("tau2tau1puppi_sig_mass2","",30,0.1,1);
  TH1F* tau2tau1puppi_sig_mass3 = new TH1F("tau2tau1puppi_sig_mass3","",30,0.1,1);
  TH1F* tau2tau1puppi_bkg       = new TH1F("tau2tau1puppi_bkg","",30,0.1,1);
  tau2tau1puppi_sig_mass1->Sumw2();
  tau2tau1puppi_sig_mass2->Sumw2();
  tau2tau1puppi_sig_mass3->Sumw2();
  tau2tau1puppi_bkg->Sumw2();

  // pileup weight
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");

  TFile* triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet.root");
  TEfficiency*       triggermet       = (TEfficiency*) triggerfile_MET->Get("efficiency");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();


  TTreeReader reader(signal_mass1);
  TTreeReaderValue<vector<float> > boostedJetpt    (reader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (reader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetphi   (reader,"boostedJetphi");
  TTreeReaderValue<vector<float> > boostedJetm     (reader,"boostedJetm");
  TTreeReaderValue<vector<float> > boostedJettau2  (reader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (reader,"boostedJettau1");
  TTreeReaderValue<vector<float> > prunedJetm      (reader,"prunedJetm");

  TTreeReaderValue<vector<float> > boostedPuppipt    (reader,"boostedPuppiJetpt");
  TTreeReaderValue<vector<float> > boostedPuppieta   (reader,"boostedPuppiJeteta");
  TTreeReaderValue<vector<float> > boostedPuppiphi   (reader,"boostedPuppiJetphi");
  TTreeReaderValue<vector<float> > boostedPuppim     (reader,"boostedPuppiJetm");
  TTreeReaderValue<vector<float> > boostedPuppitau2  (reader,"boostedPuppiJettau2");
  TTreeReaderValue<vector<float> > boostedPuppitau1  (reader,"boostedPuppiJettau1");
  TTreeReaderValue<vector<float> > softDropPuppim (reader,"softDropPuppiJetm");
  
  TTreeReaderValue<vector<float> > jetpt    (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta   (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi   (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetCHfrac (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > jetNHfrac (reader,"combinejetNHfrac");

  TTreeReaderValue<unsigned int> nbjets (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> ntaus  (reader,"ntausrawold");
  TTreeReaderValue<unsigned int> nphotons (reader,"nphotons");
  TTreeReaderValue<unsigned int> nvtx   (reader,"nvtx");

  TTreeReaderValue<float> mmet   (reader,"t1mumet");
  TTreeReaderValue<float> jmmetdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> wzpt   (reader,"wzpt");
  TTreeReaderValue<float> wgt    (reader,"wgt");

  cout<<"Loop on first signal sample "<<endl;
  while(reader.Next()){
    if(*ntaus != 0) continue;
    if(*nbjets != 0) continue;
    if(*nphotons !=0) continue;
    if(*jmmetdphi < 0.5) continue;
    if(*mmet < 250) continue;
    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetCHfrac->at(0) < 0.1) continue;
    if(jetNHfrac->at(0) > 0.1) continue;
    if(boostedJetpt->size() == 0  or boostedPuppipt->size() == 0) continue;
    
    float puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    float trgwgt = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));
    
    if(boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4){
      prunedmass_sig_mass1->Fill(prunedJetm->at(0),*wgt*puwgt*trgwgt);
      tau2tau1_sig_mass1->Fill(boostedJettau2->at(0)/boostedJettau1->at(0),*wgt*puwgt*trgwgt);
    }
    if(boostedPuppipt->at(0) > 250 and fabs(boostedPuppieta->at(0)) < 2.4){
      softdropmass_sig_mass1->Fill(softDropPuppim->at(0),*wgt*puwgt*trgwgt);
      tau2tau1puppi_sig_mass1->Fill(boostedPuppitau2->at(0)/boostedPuppitau1->at(0),*wgt*puwgt*trgwgt);
    }
  }

  // signal 2
  cout<<"Loop on second signal sample "<<endl;
  reader.SetTree(signal_mass2);
  reader.SetEntry(0);
  while(reader.Next()){
    if(*ntaus != 0) continue;
    if(*nbjets != 0) continue;
    if(*nphotons !=0) continue;
    if(*jmmetdphi < 0.5) continue;
    if(*mmet < 250) continue;
    if(jetpt->size() == 0) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetCHfrac->at(0) < 0.1) continue;
    if(jetNHfrac->at(0) > 0.1) continue;
    if(boostedJetpt->size() == 0 or boostedPuppipt->size() == 0) continue;
    
    float puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    float trgwgt = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));

    if(boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4){
      prunedmass_sig_mass2->Fill(prunedJetm->at(0),*wgt*puwgt*trgwgt);
      tau2tau1_sig_mass2->Fill(boostedJettau2->at(0)/boostedJettau1->at(0),*wgt*puwgt*trgwgt);
    }
    if(boostedPuppipt->at(0) > 250 and fabs(boostedPuppieta->at(0)) < 2.4){
      softdropmass_sig_mass2->Fill(softDropPuppim->at(0),*wgt*puwgt*trgwgt);
      tau2tau1puppi_sig_mass2->Fill(boostedPuppitau2->at(0)/boostedPuppitau1->at(0),*wgt*puwgt*trgwgt);
    }
    
  }

  // signal 3
  cout<<"Loop on third signal sample "<<endl;
  reader.SetTree(signal_mass3);
  reader.SetEntry(0);
  while(reader.Next()){
    if(*ntaus != 0) continue;
    if(*nbjets != 0) continue;
    if(*nphotons !=0) continue;
    if(*jmmetdphi < 0.5) continue;
    if(*mmet < 250) continue;
    if(jetpt->size() == 0) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetCHfrac->at(0) < 0.1) continue;
    if(jetNHfrac->at(0) > 0.1) continue;
    if(boostedJetpt->size() == 0 or boostedPuppipt->size() == 0) continue;

    float puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    float trgwgt = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));

    if(boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4){
      prunedmass_sig_mass3->Fill(prunedJetm->at(0),*wgt*puwgt*trgwgt);
      tau2tau1_sig_mass3->Fill(boostedJettau2->at(0)/boostedJettau1->at(0),*wgt*puwgt*trgwgt);
    }
    if(boostedPuppipt->at(0) > 250 and fabs(boostedPuppieta->at(0)) < 2.4){
      softdropmass_sig_mass3->Fill(softDropPuppim->at(0),*wgt*puwgt*trgwgt);
      tau2tau1puppi_sig_mass3->Fill(boostedPuppitau2->at(0)/boostedPuppitau1->at(0),*wgt*puwgt*trgwgt);
    }    
  }

  ////// background
  cout<<"Loop on background sample "<<endl; 
  reader.SetTree(background);
  reader.SetEntry(0);

  vector<TH1*> khists;
  TFile* kffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root","READ");

  TH1* znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
  TH1* zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
  TH1* zewkhist = (TH1*) kffile->Get("EWKcorr/Z");
  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  
  khists.push_back(zewkhist);
  khists.push_back(znlohist);
  

  while(reader.Next()){
    if(*ntaus != 0) continue;
    if(*nbjets != 0) continue;
    if(*nphotons !=0) continue;
    if(*jmmetdphi < 0.5) continue;
    if(*mmet < 250) continue;
    if(jetpt->size() == 0) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetCHfrac->at(0) < 0.1) continue;
    if(jetNHfrac->at(0) > 0.1) continue;
    if(boostedJetpt->size() == 0 or boostedPuppipt->size() == 0) continue;

    float puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    float trgwgt = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));
    Double_t kwgt  = 1.0;
    double   genpt = *wzpt;
    
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }      
    }

    if(boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4){
      prunedmass_bkg->Fill(prunedJetm->at(0),*wgt*puwgt*trgwgt*kwgt);
      tau2tau1_bkg->Fill(boostedJettau2->at(0)/boostedJettau1->at(0),*wgt*puwgt*trgwgt*kwgt);
    }
    if(boostedPuppipt->at(0) > 250 and fabs(boostedPuppieta->at(0)) < 2.4){
      softdropmass_bkg->Fill(softDropPuppim->at(0),*wgt*puwgt*trgwgt*kwgt);
      tau2tau1puppi_bkg->Fill(boostedPuppitau2->at(0)/boostedPuppitau1->at(0),*wgt*puwgt*trgwgt*kwgt);    
    }
  }

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  drawPlot(canvas,prunedmass_sig_mass1,prunedmass_sig_mass2,prunedmass_sig_mass3,prunedmass_bkg,outputDIR,"pruned mass [GeV]");
  drawPlot(canvas,softdropmass_sig_mass1,softdropmass_sig_mass2,softdropmass_sig_mass3,softdropmass_bkg,outputDIR,"softdrop+puppi mass [GeV]");
  drawPlot(canvas,tau2tau1_sig_mass1,tau2tau1_sig_mass2,tau2tau1_sig_mass3,tau2tau1_bkg,outputDIR,"#tau_{2}/#tau_{1}");
  drawPlot(canvas,tau2tau1puppi_sig_mass1,tau2tau1puppi_sig_mass2,tau2tau1puppi_sig_mass3,tau2tau1puppi_bkg,outputDIR,"#tau_{2}/#tau_{1} puppi");

}
