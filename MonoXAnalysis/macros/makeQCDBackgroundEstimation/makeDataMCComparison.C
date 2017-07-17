#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float lowMetBound  = 100;
static float highMetBound = 250;
static float luminosity = 35.9;

void loadChain(const string & inputPath, TChain* chain, const bool & isEOS){

  if(isEOS)
    system(("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+inputPath+" -name \"*.root\" | grep -v failed > file.temp").c_str());
  else
    system(("find "+inputPath+" -name \"*.root\" | grep -v failed > file.temp").c_str());
  
  ifstream infile;
  string line;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      if(isEOS)
	chain->Add(("root://eoscms.cern.ch//"+line).c_str());
      else
	chain->Add(line.c_str());
    }
  }
  infile.close();
  system("rm file.temp");
}

void fillHistograms(TChain* chain, vector<TH1F*> & histoA, vector<TH1F*> & histoB, vector<TH1F*> & histoC, const bool & isMC, const Category & category){

  TFile* kffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_24bins.root");
  TH1* znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
  TH1* zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile->Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  
  TH1* wnlohist = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
  TH1* wlohist  = (TH1*) kffile->Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile->Get("EWKcorr/W");
  
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  vector<TH1*> khists;
  if(TString(histoA.at(0)->GetName()).Contains("ZJets")){
    khists.push_back(zewkhist);
    khists.push_back(znlohist);
  }
  else if(TString(histoA.at(0)->GetName()).Contains("WJets")){
    khists.push_back(wewkhist);
    khists.push_back(wnlohist);
  }
  
  TTreeReader reader(chain);

  TTreeReaderValue<UChar_t> hltPFHT350 (reader,"hltPFHT350");
  TTreeReaderValue<UChar_t> hltPFHT400 (reader,"hltPFHT400");
  TTreeReaderValue<UChar_t> hltPFHT475 (reader,"hltPFHT475");
  TTreeReaderValue<UChar_t> hltPFHT600 (reader,"hltPFHT600");
  TTreeReaderValue<UChar_t> hltPFHT650 (reader,"hltPFHT650");
  TTreeReaderValue<UChar_t> hltPFHT800 (reader,"hltPFHT800");
  TTreeReaderValue<UChar_t> hltPFHT900 (reader,"hltPFHT900");  

  TTreeReaderValue<float> pswgt_ht350 (reader,"pswgt_ht350");
  TTreeReaderValue<float> pswgt_ht400 (reader,"pswgt_ht400");
  TTreeReaderValue<float> pswgt_ht475 (reader,"pswgt_ht475");
  TTreeReaderValue<float> pswgt_ht600 (reader,"pswgt_ht600");
  TTreeReaderValue<float> pswgt_ht650 (reader,"pswgt_ht650");
  TTreeReaderValue<float> pswgt_ht800 (reader,"pswgt_ht800");
  TTreeReaderValue<float> pswgt_ht900 (reader,"pswgt_ht900");

  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<float> jmmdphi     (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> metcalo     (reader,"calomet");

  TTreeReaderValue<float> xsec        (reader,"xsec");
  TTreeReaderValue<float> met         (reader,"t1pfmet");

  string pubranch = "wgtpileup";
  string btagbranch = "wgtbtag";
  TTreeReaderValue<double>* wgtsum = NULL;
  if(not isMC){
    pubranch = "wgt";
    btagbranch = "wgt";
  }

  TTreeReaderValue<float> wgtpileup   (reader,pubranch.c_str());
  TTreeReaderValue<float> wgtbtag     (reader,btagbranch.c_str());
  if(isMC)
    wgtsum = new TTreeReaderValue<double> (reader,"wgtsum");
  TTreeReaderValue<float> wgt         (reader,"wgt");
  
  long int totalEntries = chain->GetEntries();
  long int nPart = 100000;
  long int nEvents = 0;
  double sumwgt  = 0;

  while(reader.Next()){

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/totalEntries*100<<" % ";
    nEvents++;

    // require one of the trigger to be accpted
    int passHT = *hltPFHT350+*hltPFHT400+*hltPFHT475+*hltPFHT600+*hltPFHT650+*hltPFHT800+*hltPFHT900;
    if( passHT == 0) continue;

    if(fabs(*met-*metcalo)/(*met) > 0.5) continue;

    TLorentzVector jet1, jet2;
    if(category == Category::VBF or category == Category::VBFrelaxed){
      if(jetpt->size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7) continue;
      if(fabs(jeteta->at(1)) > 4.7) continue;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    }
        
    if(category == Category::VBFrelaxed){
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1.0) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.3) continue;
      if((jet1+jet2).M() < 200) continue;
    }
    else if(category == Category::VBF){
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 4.0) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      if((jet1+jet2).M() < 1300) continue;
    }

    float ht = 0;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      if(jetpt->at(ijet) > 30 and fabs(jeteta->at(ijet)) < 4.7)
	ht += jetpt->at(ijet);
    }

    //if(ht < 350) continue;

    //Gen level info --> NLO re-weight                                                                                                                                                                
    Double_t kwgt = 1.0;
    double genpt = *met;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
        if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
        if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
        kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }


    double evtweight = 1 ;
    if(isMC)
      evtweight = *xsec*luminosity*(*wgt)*(*wgtpileup)*(*wgtbtag)*kwgt/(**wgtsum);
    else{ // check trigger prescales
      if(*hltPFHT900) evtweight = *pswgt_ht900;
      else if(not *hltPFHT900 and *hltPFHT800)  evtweight = *pswgt_ht800;
      else if(not *hltPFHT900 and not *hltPFHT800 and *hltPFHT650) evtweight = *pswgt_ht650;
      else if(not *hltPFHT900 and not *hltPFHT800 and not *hltPFHT650 and *hltPFHT600) evtweight = *pswgt_ht600;
      else if(not *hltPFHT900 and not *hltPFHT800 and not *hltPFHT650 and not *hltPFHT600 and *hltPFHT475) evtweight = *pswgt_ht475;
      else if(not *hltPFHT900 and not *hltPFHT800 and not *hltPFHT650 and not *hltPFHT600 and not *hltPFHT475 and *hltPFHT400) evtweight = *pswgt_ht400;
      else if(not *hltPFHT900 and not *hltPFHT800 and not *hltPFHT650 and not *hltPFHT600 and not *hltPFHT475 and not *hltPFHT400 and *hltPFHT350) evtweight = *pswgt_ht350;
    }
      

    // fill histo A
    if(*jmmdphi < 0.5 and *met > lowMetBound and *met < highMetBound){
      for(auto hist : histoA){
	TString name (hist->GetName());
	if(name.Contains("mjj"))
	  hist->Fill((jet1+jet2).M(),evtweight);
	else if(name.Contains("met"))
	  hist->Fill(*met,evtweight);
	else if(name.Contains("ht"))
	  hist->Fill(ht,evtweight);
	else if(name.Contains("jetpt2"))
	  hist->Fill(jet2.Pt(),evtweight);
	else if(name.Contains("jeteta2"))
	  hist->Fill(jet2.Eta(),evtweight);
	else if(name.Contains("jetpt"))
	  hist->Fill(jet1.Pt(),evtweight);
	else if(name.Contains("jeteta"))
	  hist->Fill(jet1.Eta(),evtweight);
      }
    }

    // fill histo N
    if(*jmmdphi < 0.5 and *met > highMetBound){
      for(auto hist : histoB){
	TString name (hist->GetName());
	if(name.Contains("mjj"))
	  hist->Fill((jet1+jet2).M(),evtweight);
	else if(name.Contains("met"))
	  hist->Fill(*met,evtweight);
	else if(name.Contains("ht"))
	  hist->Fill(ht,evtweight);
	else if(name.Contains("jetpt2"))
	  hist->Fill(jet2.Pt(),evtweight);
	else if(name.Contains("jeteta2"))
	  hist->Fill(jet2.Eta(),evtweight);
	else if(name.Contains("jetpt"))
	  hist->Fill(jet1.Pt(),evtweight);
	else if(name.Contains("jeteta"))
	  hist->Fill(jet1.Eta(),evtweight);
      }
    }
    
    // fill histo N
    if(*jmmdphi > 0.5 and *met > lowMetBound and *met < highMetBound){
      for(auto hist : histoC){
	TString name (hist->GetName());
	if(name.Contains("mjj"))
	  hist->Fill((jet1+jet2).M(),evtweight);
	else if(name.Contains("met"))
	  hist->Fill(*met,evtweight);
	else if(name.Contains("ht"))
	  hist->Fill(ht,evtweight);
	else if(name.Contains("jetpt2"))
	  hist->Fill(jet2.Pt(),evtweight);
	else if(name.Contains("jeteta2"))
	  hist->Fill(jet2.Eta(),evtweight);
	else if(name.Contains("jetpt"))
	  hist->Fill(jet1.Pt(),evtweight);
	else if(name.Contains("jeteta"))
	  hist->Fill(jet1.Eta(),evtweight);
      }
    }
  }
  cout<<endl;
  if(kffile) kffile->Close();
}

void plotDataMC(TH1* histoData, THStack* histoBkg, TCanvas* canvas, const string & outputDIR){

  canvas->cd();
  histoData->SetMarkerColor(kBlack);
  histoData->SetMarkerStyle(20);
  histoData->SetMarkerSize(1.0);
  histoData->Draw("EP");

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  // total background from the Stack
  TH1* totalBkg = (TH1*) histoBkg->GetStack()->At(histoBkg->GetNhists()-1);
  TString name (histoData->GetName());
  
  if(not name.Contains("jeteta"))
    histoData->GetYaxis()->SetRangeUser(min(totalBkg->GetMinimum(),histoData->GetMinimum())*0.01,max(totalBkg->GetMaximum(),histoData->GetMaximum())*100);
  else
    histoData->GetYaxis()->SetRangeUser(0,max(totalBkg->GetMaximum(),histoData->GetMaximum())*1.30);

  histoData->GetYaxis()->SetTitle("Events");
  histoData->GetYaxis()->SetTitleOffset(1.10);
  histoData->GetYaxis()->SetLabelSize(0.035);

  TH1* nhist   = (TH1*) histoData->Clone("nhist");
  TH1* dhist   = (TH1*) totalBkg->Clone("dhist");
  TH1* dhist_2 = (TH1*) totalBkg->Clone("dhist_2");
  TH1* unhist  = (TH1*) histoData->Clone("unhist");

  histoData->GetXaxis()->SetTitleSize(0);
  histoData->GetXaxis()->SetLabelSize(0);

  histoBkg->Draw("hist same");
  histoData->Draw("same");

  TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histoData,"Data","EP");
  TObjArray* listHisto = histoBkg->GetStack();
  TIter iObject(listHisto);
  while (TH1F* histo = dynamic_cast<TH1F*>(iObject())) {
    TString name (histo->GetName());
    if(name.Contains("QCD"))
      leg->AddEntry(histo,"QCD Multijet","F");
    else if(name.Contains("WJets"))
      leg->AddEntry(histo,"W+jets","F");
    else if(name.Contains("ZJets"))
      leg->AddEntry(histo,"Z+jets","F");    
  }
  
  leg->Draw("same");
  
  canvas->RedrawAxis("sameaxis");
  CMS_lumi(canvas,"");

  pad2->Draw();
  pad2->cd();
  
  nhist->GetYaxis()->SetTitle("Data/Pred");
  nhist->GetYaxis()->CenterTitle();
  nhist->GetYaxis()->SetTitleOffset(1.5);
  nhist->GetYaxis()->SetLabelSize(0.035);
  nhist->GetYaxis()->SetTitleSize(0.04);
  nhist->GetXaxis()->SetLabelSize(0.04);
  nhist->GetXaxis()->SetTitleSize(0.05);
  
  if(name.Contains("mjj"))
    nhist->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(name.Contains("jetpt2"))
    nhist->GetXaxis()->SetTitle("p_{T}^{j2} [GeV]");
  else if(name.Contains("jetpt"))
    nhist->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
  else if(name.Contains("jeteta2"))
    nhist->GetXaxis()->SetTitle("#eta^{j2}");
  else if(name.Contains("jeteta"))
    nhist->GetXaxis()->SetTitle("#eta^{j1}");
  else if(name.Contains("met"))
    nhist->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  else if(name.Contains("ht"))
    nhist->GetXaxis()->SetTitle("H_{T} [GeV]");

  nhist->GetXaxis()->SetTitleOffset(1.10);

  nhist->GetYaxis()->SetRangeUser(0,2);
  
  for (int i = 1; i <= dhist->GetNbinsX(); i++) dhist->SetBinError(i, 0);
  nhist->Divide(dhist);
  dhist_2->Divide(dhist);
  
  dhist_2->SetLineColor(0);
  dhist_2->SetMarkerColor(0);
  dhist_2->SetMarkerSize(0);
  dhist_2->SetFillColor(kGray);

  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetFillColor(0);
  unhist->SetLineWidth(2);
  
  nhist->GetYaxis()->SetNdivisions(505);

  nhist->Draw("PE1 SAME");
  dhist_2->Draw("E2 SAME");
  unhist->Draw("L SAME");
  nhist->Draw("PE SAME");

  string postfix = "DataMC_";
  if(name.Contains("regionA")) postfix += "RegionA_";
  if(name.Contains("regionB")) postfix += "RegionB_";
  if(name.Contains("regionC")) postfix += "RegionC_";
  if(name.Contains("mjj")) postfix += "mjj";
  else if(name.Contains("jetpt2")) postfix += "jetpt2";
  else if(name.Contains("jeteta2")) postfix += "jeteta2";
  else if(name.Contains("jetpt")) postfix += "jetpt";
  else if(name.Contains("jeteta")) postfix += "jeteta";
  else if(name.Contains("met")) postfix += "met";
  else if(name.Contains("ht")) postfix += "ht";
  
  pad2->RedrawAxis("sameaxis");

  if(not name.Contains("jeteta"))
    canvas->SetLogy();
  else
    canvas->SetLogy(0);


  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  if(leg) delete leg;
  if(nhist) delete nhist;
  if(dhist) delete dhist;
  if(dhist_2) delete dhist_2;
  if(unhist) delete unhist;
  if(pad2) delete pad2;
  
}

void bookHistogram(vector<TH1F*> & histoA, vector<TH1F*> & histoB, vector<TH1F*> & histoC, const string & postfix, const Category & category){

  // variable binning in the different regions (Region A)
  vector<double> bins_mjj_A     = {200.,350.,500.,650.,800.,1000.,1200.,1400.,1600.,1800.,2000.,2350.,2800.,3500,5000};
  vector<double> bins_jeteta_A  = selectBinning("jeteta",category);
  vector<double> bins_jeteta2_A = selectBinning("jeteta2",category);
  vector<double> bins_ht_A      = {250.,350.,400.,450.,500.,550.,600.,650,700.,750,800,850,950,1050,1150,1250};
  vector<double> bins_met_A     = {100.,115.,130.,145.,160.,175.,190.,210.,230.,250.};

  // variable binning in the different regions (Region B)
  vector<double> bins_mjj_B     = selectBinning("mjj",category);
  vector<double> bins_jeteta_B  = selectBinning("jeteta",category);
  vector<double> bins_jeteta2_B = selectBinning("jeteta2",category);
  vector<double> bins_ht_B      = {250.,350.,450.,600.,750,900,1050,1200,1400,1600,1800,2000,2500};
  vector<double> bins_met_B     = {250.,300.,350.,450.,600.,800.,1000.};

  // variable binning in the different regions (Region C)
  vector<double> bins_mjj_C     = {200.,450.,700.,1000.,1300.,1600.,2000.,2500.,3250.,5000.};
  vector<double> bins_jeteta_C  = selectBinning("jeteta",category);
  vector<double> bins_jeteta2_C = selectBinning("jeteta2",category);
  vector<double> bins_ht_C      = {250.,350.,400,450.,500,550.,600.,650.,700.,775.,850.,950.,1050.,1150,1250.};  
  vector<double> bins_met_C     = {100.,120.,140.,160.,180.,200.,225.,250.};
  
  histoA.push_back(new TH1F(Form("histo_%s_mjj_regionA",postfix.c_str()),"",bins_mjj_A.size()-1,&bins_mjj_A[0]));
  histoB.push_back(new TH1F(Form("histo_%s_mjj_regionB",postfix.c_str()),"",bins_mjj_B.size()-1,&bins_mjj_B[0]));
  histoC.push_back(new TH1F(Form("histo_%s_mjj_regionC",postfix.c_str()),"",bins_mjj_C.size()-1,&bins_mjj_C[0]));
  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();

  histoA.push_back(new TH1F(Form("histo_%s_ht_regionA",postfix.c_str()),"",bins_ht_A.size()-1,&bins_ht_A[0]));
  histoB.push_back(new TH1F(Form("histo_%s_ht_regionB",postfix.c_str()),"",bins_ht_B.size()-1,&bins_ht_B[0]));
  histoC.push_back(new TH1F(Form("histo_%s_ht_regionC",postfix.c_str()),"",bins_ht_C.size()-1,&bins_ht_C[0]));
  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();
  
  histoA.push_back(new TH1F(Form("histo_%s_met_regionA",postfix.c_str()),"",bins_met_A.size()-1,&bins_met_A[0]));
  histoB.push_back(new TH1F(Form("histo_%s_met_regionB",postfix.c_str()),"",bins_met_B.size()-1,&bins_met_B[0]));
  histoC.push_back(new TH1F(Form("histo_%s_met_regionC",postfix.c_str()),"",bins_met_C.size()-1,&bins_met_C[0]));
  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();

  histoA.push_back(new TH1F(Form("histo_%s_jeteta_regionA",postfix.c_str()),"",bins_jeteta_A.size()-1,&bins_jeteta_A[0]));
  histoB.push_back(new TH1F(Form("histo_%s_jeteta_regionB",postfix.c_str()),"",bins_jeteta_B.size()-1,&bins_jeteta_B[0]));
  histoC.push_back(new TH1F(Form("histo_%s_jeteta_regionC",postfix.c_str()),"",bins_jeteta_C.size()-1,&bins_jeteta_C[0]));
  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();

  histoA.push_back(new TH1F(Form("histo_%s_jeteta2_regionA",postfix.c_str()),"",bins_jeteta2_A.size()-1,&bins_jeteta2_A[0]));
  histoB.push_back(new TH1F(Form("histo_%s_jeteta2_regionB",postfix.c_str()),"",bins_jeteta2_B.size()-1,&bins_jeteta2_B[0]));
  histoC.push_back(new TH1F(Form("histo_%s_jeteta2_regionC",postfix.c_str()),"",bins_jeteta2_C.size()-1,&bins_jeteta2_C[0]));
  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();
}

vector<THStack*> createBackgroundStack(const vector<TH1F*> & histoWJets, const vector<TH1F*> & histoZJets, const vector<TH1F*> & histoQCD){

  // create new empty stacks
  vector<THStack*> bkgStack;
  for(auto hist : histoQCD){
    TString name = Form("%s",hist->GetName());
    name.ReplaceAll("QCD","Bkg");    
    bkgStack.push_back(new THStack(name,name));    
  }

  // Fill the stack and set the right colors
  for(size_t ihist = 0; ihist < histoQCD.size(); ihist++){

    if(histoWJets.size() != 0){
      histoWJets.at(ihist)->SetLineColor(kBlack);
      histoWJets.at(ihist)->SetFillColor(TColor::GetColor("#FAAF08"));
      bkgStack.at(ihist)->Add(histoWJets.at(ihist));
    }

    if(histoZJets.size() != 0){
      histoZJets.at(ihist)->SetLineColor(kBlack);
      histoZJets.at(ihist)->SetFillColor(TColor::GetColor("#4D975D"));
      bkgStack.at(ihist)->Add(histoZJets.at(ihist));
    }

    histoQCD.at(ihist)->SetLineColor(kBlack);
    histoQCD.at(ihist)->SetFillColor(TColor::GetColor("#F1F1F2"));
    bkgStack.at(ihist)->Add(histoQCD.at(ihist));
  }

  return bkgStack;
  
}

///////////
void makeDataMCComparison (Category category, string outputDIR, bool addEWKBackgrounds, bool isEOS){

  string inputPathData = "/home/rgerosa/MONOJET_ANALYSIS/QCDEstimation/JetHT_DATA/";
  string inputPathQCDMC = "/home/rgerosa/MONOJET_ANALYSIS/QCDEstimation/QCD_MC/";
  string inputPathWJetsMC = "/home/rgerosa/MONOJET_ANALYSIS/QCDEstimation/WJets_MC/";
  string inputPathZJetsMC = "/home/rgerosa/MONOJET_ANALYSIS/QCDEstimation/ZJets_MC/";

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  initializeBinning();
  gROOT->SetBatch(kTRUE);
  
  TChain* chainQCD  = new TChain("tree/tree");  
  TChain* chainData = new TChain("tree/tree");  
  TChain* chainWJets = NULL;
  TChain* chainZJets = NULL;

  loadChain(inputPathQCDMC,chainQCD,isEOS);
  loadChain(inputPathData,chainData,isEOS);
  
  if(addEWKBackgrounds){
    chainWJets = new TChain("tree/tree");
    chainZJets = new TChain("tree/tree");
    loadChain(inputPathWJetsMC,chainWJets,isEOS);
    loadChain(inputPathZJetsMC,chainZJets,isEOS);
  }
  
  // Book QCD MC and Fill
  vector<TH1F*> histoQCD_A, histoQCD_B, histoQCD_C;
  cout<<"Run on QCD-MC "<<endl;
  bookHistogram(histoQCD_A,histoQCD_B,histoQCD_C,"QCD",category);
  fillHistograms(chainQCD,histoQCD_A,histoQCD_B,histoQCD_C,true,category);

  // Book WJets MC and Fill
  vector<TH1F*> histoWJets_A, histoWJets_B, histoWJets_C;
  if(addEWKBackgrounds){
    cout<<"Run on WJets-MC "<<endl;
    bookHistogram(histoWJets_A,histoWJets_B,histoWJets_C,"WJets",category);
    fillHistograms(chainWJets,histoWJets_A,histoWJets_B,histoWJets_C,true,category);
  }

  // Book ZJets MC and Fill
  vector<TH1F*> histoZJets_A, histoZJets_B, histoZJets_C;
  if(addEWKBackgrounds){
    cout<<"Run on ZJets-MC "<<endl;
    bookHistogram(histoZJets_A,histoZJets_B,histoZJets_C,"ZJets",category);
    fillHistograms(chainZJets,histoZJets_A,histoZJets_B,histoZJets_C,true,category);
  }

  
  // Book Data and Fill
  vector<TH1F*> histoData_A, histoData_B, histoData_C;
  cout<<"Run on Data "<<endl;
  bookHistogram(histoData_A,histoData_B,histoData_C,"Data",category);
  fillHistograms(chainData,histoData_A,histoData_B,histoData_C,false,category);

  // Create THStack for the backgrounds
  vector<THStack*> histoBkg_A = createBackgroundStack(histoWJets_A,histoZJets_A,histoQCD_A);
  vector<THStack*> histoBkg_B = createBackgroundStack(histoWJets_B,histoZJets_B,histoQCD_B);
  vector<THStack*> histoBkg_C = createBackgroundStack(histoWJets_C,histoZJets_C,histoQCD_C);

  cout<<"W+jets rate in Region A: "<<histoWJets_A.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region A: "<<histoZJets_A.front()->Integral()<<endl;
  cout<<"QCD rate in Region A: "<<histoQCD_A.front()->Integral()<<endl;
  cout<<"Data rate in Region A: "<<histoData_A.front()->Integral()<<endl;

  cout<<"W+jets rate in Region B: "<<histoWJets_B.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region B: "<<histoZJets_B.front()->Integral()<<endl;
  cout<<"QCD rate in Region B: "<<histoQCD_B.front()->Integral()<<endl;
  cout<<"Data rate in Region B: "<<histoData_B.front()->Integral()<<endl;

  cout<<"W+jets rate in Region C: "<<histoWJets_C.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region C: "<<histoZJets_C.front()->Integral()<<endl;
  cout<<"QCD rate in Region C: "<<histoQCD_C.front()->Integral()<<endl;
  cout<<"Data rate in Region C: "<<histoData_C.front()->Integral()<<endl;
  
  //////  
  cout<<"Start plotting "<<endl;
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  for(size_t ihist = 0; ihist < histoData_A.size(); ihist++)
    plotDataMC(histoData_A.at(ihist),histoBkg_A.at(ihist),canvas,outputDIR);

  for(size_t ihist = 0; ihist < histoData_B.size(); ihist++)
    plotDataMC(histoData_B.at(ihist),histoBkg_B.at(ihist),canvas,outputDIR);

  for(size_t ihist = 0; ihist < histoData_C.size(); ihist++)
    plotDataMC(histoData_C.at(ihist),histoBkg_C.at(ihist),canvas,outputDIR);
  
}
