#include "../CMS_lumi.h"

void fillHisto(TH1* histo, float val, float weight){ // embed the overflow                                                                                                                       
 
  if(val < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    histo->Fill(val,weight);
  else
    histo->Fill(histo->GetXaxis()->GetBinCenter(histo->GetNbinsX()),weight);

}

void fillHisto2D(TH2F* histo, float valx, float valy, float weight){ // Embed- the overflow                                                                                                       
  
  float x = 0;
  if(valx < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    x = valx;
  else
    x =histo->GetXaxis()->GetBinCenter(histo->GetNbinsX());
  float y = 0;
  if(valy < histo->GetYaxis()->GetBinLowEdge(histo->GetNbinsY()+1))
    y = valy;
  else
    y = histo->GetYaxis()->GetBinCenter(histo->GetNbinsY());

  histo->Fill(x,y,weight);
}


//// ---
void createHistograms1D(vector<TH1F*> & histo, const string & postfix){

  histo.push_back(new TH1F(Form("%s_met",postfix.c_str()),"",20,250,1250));
  histo.push_back(new TH1F(Form("%s_jetpt1",postfix.c_str()),"",20,80,500));
  histo.push_back(new TH1F(Form("%s_jetpt2",postfix.c_str()),"",20,40,400));
  histo.push_back(new TH1F(Form("%s_jeteta1",postfix.c_str()),"",20,-4.5,4.5));
  histo.push_back(new TH1F(Form("%s_jeteta2",postfix.c_str()),"",20,-4.5,4.5));
  histo.push_back(new TH1F(Form("%s_dphijj",postfix.c_str()),"",20,0,3.14));
  histo.push_back(new TH1F(Form("%s_detajj",postfix.c_str()),"",20,0,9));
  histo.push_back(new TH1F(Form("%s_mjj",postfix.c_str()),"",20,0,4000));

  for(auto hist: histo){
    hist->Sumw2();
    hist->SetDirectory(0);
  }
 
}

//// ---
void createHistograms2D(vector<TH2F*> & histo, const string postfix){
  
  histo.push_back(new TH2F(Form("%s_mjj_detajj",postfix.c_str()),"",15,400,2500,10,2,8));
  histo.push_back(new TH2F(Form("%s_mjj_dphijj",postfix.c_str()),"",15,400,2500,10,0.,3.14));
  histo.push_back(new TH2F(Form("%s_mjj_met",postfix.c_str()),"",15,400,2500,10,0.,3.14));
  histo.push_back(new TH2F(Form("%s_detajj_dphijj",postfix.c_str()),"",10,2,8,10,0.,3.14));
  histo.push_back(new TH2F(Form("%s_detajj_met",postfix.c_str()),"",10,2,8,10,0.,3.14));

  for(auto hist: histo){
    hist->Sumw2();
    hist->SetDirectory(0);
  }

}


///// ---
void plotDistribution(TCanvas* canvas, TH1F* qqH, TH1F* ggH, TH1F* zjet, TH1F* wjet, TH1F* ewkz, TH1F* ewkw, const string & outputDIR){
  
  // normalize to 1
  canvas->cd();

  qqH->Scale(1./qqH->Integral());
  ggH->Scale(1./ggH->Integral());
  zjet->Scale(1./zjet->Integral());
  wjet->Scale(1./wjet->Integral());
  ewkz->Scale(1./ewkz->Integral());
  ewkw->Scale(1./ewkw->Integral());

  TString name (qqH->GetName());
  string label = "";
  bool isLog = false;
  if(name.Contains("detajj"))
    label = "#Delta#eta_{jj}";
  else if(name.Contains("mjj")){
    label = "M_{jj} [GeV]";
    isLog = true;
  }
  else if(name.Contains("dphijj"))
    label = "#Delta#phi_{jj}";
  else if(name.Contains("jeteta1"))
    label = "#eta_{j1}";
  else if(name.Contains("jeteta2"))
    label = "#eta_{j2}";
  else if(name.Contains("jetpt1")){
    label = "p_{T}^{j1} [GeV]";
    isLog = true;
  }
  else if(name.Contains("jetpt2")){
    label = "p_{T}^{j2} [GeV]";
    isLog = true;
  }
  else if(name.Contains("met")){
    label = "E_{T}^{miss} [GeV]";
    isLog = true;
  }

  zjet->GetXaxis()->SetTitle(label.c_str());
  zjet->GetYaxis()->SetTitle("a.u");

  zjet->SetLineColor(kRed+1);
  zjet->SetLineWidth(3);
  zjet->Draw("hist");

  wjet->SetLineColor(kOrange+8);
  wjet->SetLineWidth(3);
  wjet->Draw("hist same");
  
  ewkz->SetLineColor(kBlue);
  ewkz->SetLineWidth(3);
  ewkz->Draw("hist same");

  ewkw->SetLineColor(kAzure);
  ewkw->SetLineWidth(3);
  ewkw->Draw("hist same");

  ggH->SetLineColor(kBlack);
  ggH->SetLineWidth(5);
  ggH->SetLineStyle(2);
  ggH->Draw("hist same");

  qqH->SetLineColor(kBlack);
  qqH->SetLineWidth(5);
  qqH->Draw("hist same");
  
  zjet->GetYaxis()->SetRangeUser(0,max(qqH->GetMaximum(),max(ggH->GetMaximum(),max(zjet->GetMaximum(),max(wjet->GetMaximum(),max(ewkz->GetMaximum(),ewkw->GetMaximum())))))*1.5);
  if(isLog)
    zjet->GetYaxis()->SetRangeUser(0.00001,max(qqH->GetMaximum(),max(ggH->GetMaximum(),max(zjet->GetMaximum(),max(wjet->GetMaximum(),max(ewkz->GetMaximum(),ewkw->GetMaximum())))))*10);

  CMS_lumi(canvas,"");

  TLegend leg(0.7,0.7,0.9,0.9);
  leg.SetFillStyle(1001);
  leg.SetFillColor(0);
  leg.AddEntry(qqH,"VBF H = 125 GeV","L");
  leg.AddEntry(ggH,"gg H = 125 GeV","L");
  leg.AddEntry(zjet,"Z#nu#nu (QCD)","L");
  leg.AddEntry(wjet,"W+jets (QCD)","L");
  leg.AddEntry(ewkz,"Z#nu#nu (EW)","L");
  leg.AddEntry(ewkw,"W+jets (EW)","L");
  leg.Draw("same");

  if(isLog)
    canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  
}



void plotDistribution2D(TCanvas* canvas, TH2F* histo, const string & outputDIR){
  
  system(("mkdir -p "+outputDIR+"/correlation/").c_str());
  canvas->cd();
  canvas->SetRightMargin(0.18);

  ///// --- normalize to 1
  histo->Scale(1./histo->Integral());

  ///// ---
  TGraph2D* graph = new TGraph2D();
  graph->SetNpx(50);
  graph->SetNpy(50);
  int nPoint = 0;
  for(int iBinX = 0; iBinX < histo->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < histo->GetNbinsY() ; iBinY++){
      graph->SetPoint(nPoint,histo->GetXaxis()->GetBinCenter(iBinX+1),histo->GetYaxis()->GetBinCenter(iBinY+1),histo->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  std::stringstream name(histo->GetName());
  std::string segment;
  std::vector<std::string> seglist;
  
  while(std::getline(name, segment, '_')){
    seglist.push_back(segment);
  }

  string labelX;
  string labelY;

  TString nameVarX (seglist.at(1));
  TString nameVarY (seglist.at(2));

  if(nameVarX.Contains("detajj"))
    labelX = "#Delta#eta_{jj}";
  else if(nameVarX.Contains("mjj"))
    labelX = "M_{jj} [GeV]";
  else if(nameVarX.Contains("dphijj"))
    labelX = "#Delta#phi_{jj}";
  else if(nameVarX.Contains("jetpt1"))
    labelX = "p_{T}^{j1} [GeV]";
  else if(nameVarX.Contains("jetpt2"))
    labelX = "p_{T}^{j2} [GeV]";
  else if(nameVarX.Contains("met"))
    labelX = "E_{T}^{miss} [GeV]";

  if(nameVarY.Contains("detajj"))
    labelY = "#Delta#eta_{jj}";
  else if(nameVarY.Contains("mjj"))
    labelY = "M_{jj} [GeV]";
  else if(nameVarY.Contains("dphijj"))
    labelY = "#Delta#phi_{jj}";
  else if(nameVarY.Contains("jetpt1"))
    labelY = "p_{T}^{j1} [GeV]";
  else if(nameVarY.Contains("jetpt2"))
    labelY = "p_{T}^{j2} [GeV]";
  else if(nameVarY.Contains("met"))
    labelY = "E_{T}^{miss} [GeV]";

  /////
  TH2D* hist = graph->GetHistogram();
  hist->GetXaxis()->SetTitle(labelX.c_str());
  hist->GetYaxis()->SetTitle(labelY.c_str());
  hist->GetZaxis()->SetTitle("a.u");   
  hist->Draw("colz");

  TProfile* profile = hist->ProfileX(Form("%s_pfx",histo->GetName()));
  profile->SetMarkerColor(kBlack);
  profile->SetMarkerStyle(20);
  profile->SetMarkerSize(1);
  profile->Draw("EPsame");

  CMS_lumi(canvas,"");
  TLegend leg(0.55,0.7,0.8,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",hist->GetCorrelationFactor()),"");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/correlation/"+string(histo->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(histo->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();

}

//// -----
void makeAnalysis(TTree* tree, vector<TH1F*> hist_1D, vector<TH2F*> hist_2D, vector<TH1F*> khists){
  
  TTreeReader reader (tree);
  //// ----
  TTreeReaderValue<float> xsec         (reader,"xsec");
  TTreeReaderValue<float> wgt          (reader,"wgt");
  TTreeReaderValue<double> wgtsum       (reader,"wgtsum");
  //// ----
  TTreeReaderValue<unsigned int> njetsinc   (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nphotons   (reader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus   (reader,"ntausold");
  TTreeReaderValue<unsigned int> nmuons  (reader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets  (reader,"nbjetslowpt");
  //// ----
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  //// ----
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  //// ----
  TTreeReaderValue<vector<float> > chfrac (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac (reader,"combinejetNHfrac");
  //// ----
  TTreeReaderValue<float> wzpt        (reader,"wzpt");
  TTreeReaderValue<float> mmet        (reader,"t1mumet");
  TTreeReaderValue<float> jmmdphi4    (reader,"incjetmumetdphimin4");
  //// ----
  TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm100  (reader,"hltmetwithmu100");
  TTreeReaderValue<UChar_t> hltmwm110  (reader,"hltmetwithmu110");
  TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");

  //// ---
  while(reader.Next()){
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    // number of jets
    if (*njetsinc  < 2) continue;
    // b-veto
    if (*nbjets > 0) continue;
    // vetos
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;

    if (*hltm90 == 0 and *hltm100 == 0 and *hltm110 and *hltm120 == 0 and 
	*hltmwm90 == 0 and *hltmwm100 == 0 and *hltmwm110 == 0 and *hltmwm120 == 0) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 80) continue; 
    if (jetpt->at(1) < 40) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
    if (fabs(jeteta->at(0)) > 3 and fabs(jeteta->at(1)) > 3) continue; 
    if (jeteta->at(0)*jeteta->at(1) > 0) continue;
    // relaxed met cut
    if (*mmet     < 250) continue;
    if (*jmmdphi4 < 0.5) continue;
    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    
    /////
    Float_t kwgt = 1.0;
    float genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
        if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
        kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    float weight = *xsec*(*wgt)*kwgt/(*wgtsum);

    /////
    for(auto hist : hist_1D){
      TString name (hist->GetName());
      float variable = 0; 
      if(name.Contains("mjj"))
	fillHisto(hist,(jet1+jet2).M(),weight);
      else if(name.Contains("met"))
	fillHisto(hist,*mmet,weight);
      else if(name.Contains("detajj"))
	fillHisto(hist,fabs(jet1.Eta()-jet2.Eta()),weight);
      else if(name.Contains("dphijj"))
	fillHisto(hist,fabs(jet1.DeltaPhi(jet2)),weight);
      else if(name.Contains("jetpt1"))
	fillHisto(hist,jet1.Pt(),weight);
      else if(name.Contains("jetpt2"))
	fillHisto(hist,jet2.Pt(),weight);
      else if(name.Contains("jeteta1"))
	fillHisto(hist,jet1.Eta(),weight);
      else if(name.Contains("jeteta2"))
	fillHisto(hist,jet2.Eta(),weight);     
    }

    /////
    for(auto hist : hist_2D){
      TString name (hist->GetName());
      float variable = 0;
      if(name.Contains("mjj_detajj"))
	fillHisto2D(hist,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()),weight);
      else if(name.Contains("mjj_dphijj"))
	fillHisto2D(hist,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)),weight);
      else if(name.Contains("mjj_ptj1"))
	fillHisto2D(hist,(jet1+jet2).M(),jet1.Pt(),weight);
      else if(name.Contains("mjj_ptj2"))
	fillHisto2D(hist,(jet1+jet2).M(),jet2.Pt(),weight);
      else if(name.Contains("mjj_met"))
	fillHisto2D(hist,(jet1+jet2).M(),*mmet,weight);
      else if(name.Contains("detajj_dphijj"))
	fillHisto2D(hist,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)),weight);
      else if(name.Contains("detajj_ptj1"))
	fillHisto2D(hist,fabs(jet1.Eta()-jet2.Eta()),jet1.Pt(),weight);
      else if(name.Contains("detajj_ptj2"))
	fillHisto2D(hist,fabs(jet1.Eta()-jet2.Eta()),jet2.Pt(),weight);
      else if(name.Contains("detajj_met"))
	fillHisto2D(hist,fabs(jet1.Eta()-jet2.Eta()),*mmet,weight);
      
    }
  }
}


//// -- make plots of the distributions
void makeVBFSignalVSBackground(string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TChain* vbftree = new TChain("tree/tree");
  TChain* ggHtree = new TChain("tree/tree");
  TChain* znntree = new TChain("tree/tree");
  TChain* wjettree = new TChain("tree/tree");
  TChain* znnewktree = new TChain("tree/tree");
  TChain* wjetewktree = new TChain("tree/tree");

  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*125*root");
  ggHtree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/HiggsInvisible/sigfilter/sig_*GluGlu*_HToInvisible_M*125*root");
  znntree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/ZJets/sigfilter/*root");
  wjettree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/WJets/sigfilter/*root");
  znnewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/ZJetsToNuNuEWK/sigfilter/*root");
  wjetewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/WJetsEWK/sigfilter/*root");

  // get k-factors NLO                                                                                                                                                                                
  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_24bins.root");
  TH1F*  znlohist = (TH1F*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1F*  zlohist  = (TH1F*) kffile.Get("ZJets_LO/inv_pt");
  TH1F* zewkhist  = (TH1F*) kffile.Get("EWKcorr/Z");

  TH1F*  wnlohist = (TH1F*) kffile.Get("WJets_012j_NLO/nominal");
  TH1F*  wlohist  = (TH1F*) kffile.Get("WJets_LO/inv_pt");
  TH1F* wewkhist  = (TH1F*) kffile.Get("EWKcorr/W");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(znlohist)
    wnlohist->Divide(wlohist);
  
  vector<TH1F*> zhists;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);
  vector<TH1F*> whists;
  whists.push_back(wnlohist);
  whists.push_back(wewkhist);
  vector<TH1F*> ehists;

  // make histograms
  vector<TH1F*> qqH_histo;
  vector<TH1F*> ggH_histo;
  vector<TH1F*> zjet_histo;
  vector<TH1F*> wjet_histo;
  vector<TH1F*> ewkzjet_histo;
  vector<TH1F*> ewkwjet_histo;
  
  createHistograms1D(qqH_histo,"qqH");
  createHistograms1D(ggH_histo,"ggH");
  createHistograms1D(zjet_histo,"zjet");
  createHistograms1D(wjet_histo,"wjet");
  createHistograms1D(ewkzjet_histo,"ewkzjet");
  createHistograms1D(ewkwjet_histo,"ewkwjet");

  vector<TH2F*> qqH_histo_2D;
  vector<TH2F*> ggH_histo_2D;
  vector<TH2F*> zjet_histo_2D;
  vector<TH2F*> wjet_histo_2D;
  vector<TH2F*> ewkzjet_histo_2D;
  vector<TH2F*> ewkwjet_histo_2D;

  createHistograms2D(qqH_histo_2D,"qqH");
  createHistograms2D(ggH_histo_2D,"ggH");
  createHistograms2D(zjet_histo_2D,"zjet");
  createHistograms2D(wjet_histo_2D,"wjet");
  createHistograms2D(ewkzjet_histo_2D,"ewkzjet");
  createHistograms2D(ewkwjet_histo_2D,"ewkwjet");

  ////////////
  cout<<"Loop on VBF trees"<<endl;
  makeAnalysis(vbftree,qqH_histo,qqH_histo_2D,ehists);
  cout<<"Loop on ggH trees"<<endl;
  makeAnalysis(ggHtree,ggH_histo,ggH_histo_2D,ehists);
  cout<<"Loop on Z+jets trees"<<endl;
  makeAnalysis(znntree,zjet_histo,zjet_histo_2D,zhists);
  cout<<"Loop on W+jets trees"<<endl;
  makeAnalysis(wjettree,wjet_histo,wjet_histo_2D,whists);
  cout<<"Loop on Z+jets EW trees"<<endl;
  makeAnalysis(znnewktree,ewkzjet_histo,ewkzjet_histo_2D,ehists);
  cout<<"Loop on W+jets EW trees"<<endl;
  makeAnalysis(wjetewktree,ewkwjet_histo,ewkwjet_histo_2D,ehists);

  //////////////////////////
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  /////////////
  cout<<"Plot 1D distributions "<<endl;
  for(int ihist = 0; ihist < qqH_histo.size(); ihist++)
    plotDistribution(canvas,qqH_histo.at(ihist),ggH_histo.at(ihist),zjet_histo.at(ihist),wjet_histo.at(ihist),ewkzjet_histo.at(ihist),ewkwjet_histo.at(ihist),outputDIR);

  cout<<"Plot 2D distributions "<<endl;
  ///////////// qqH
  for(int ihist = 0; ihist < qqH_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,qqH_histo_2D.at(ihist),outputDIR);

  ///////////// ggH
  for(int ihist = 0; ihist < ggH_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,ggH_histo_2D.at(ihist),outputDIR);

  ///////////// Z+jets
  for(int ihist = 0; ihist < zjet_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,zjet_histo_2D.at(ihist),outputDIR);

  ///////////// W+jets
  for(int ihist = 0; ihist < wjet_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,wjet_histo_2D.at(ihist),outputDIR);

  ///////////// EWK-Z+jets
  for(int ihist = 0; ihist < ewkzjet_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,ewkzjet_histo_2D.at(ihist),outputDIR);

  ///////////// EWK-W+jets
  for(int ihist = 0; ihist < ewkwjet_histo_2D.size(); ihist++)
    plotDistribution2D(canvas,ewkwjet_histo_2D.at(ihist),outputDIR);

}
