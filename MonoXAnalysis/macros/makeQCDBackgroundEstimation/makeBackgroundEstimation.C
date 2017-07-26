#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float lowMetBound  = 100;
static float highMetBound = 160;
static float luminosity   = 35.9;
static float scale_wjet   = 0.15;
static float scale_zjet   = 0.15;

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

void fillHistograms(TChain* chain, vector<TH1F*> & histoA, vector<TH1F*> & histoB, vector<TH1F*> & histoC, vector<TH1F*> & histoD, const bool & isMC, const Category & category){

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
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
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
      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
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

    // fill histo B
    if(*jmmdphi < 0.5 and *met > 250){
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
    
    // fill histo C
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

    // fill histo D
    if(*jmmdphi > 0.5 and *met > 250){
      for(auto hist : histoD){
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

void bookHistogram(vector<TH1F*> & histoA, vector<TH1F*> & histoB, vector<TH1F*> & histoC, vector<TH1F*> & histoD, const string & postfix, const Category & category, const int binning = 1){

  // variable binning in the different regions (Region A)
  vector<double> bins_mjj;
  if(binning == 1)
    bins_mjj = {200,400,600,900,1200,1500,2000,2750,3500,5000};
  else 
    bins_mjj  = {200.,500.,800.,1300.,2000.,2800,5000.};

  histoA.push_back(new TH1F(Form("histo_%s_mjj_regionA",postfix.c_str()),"",bins_mjj.size()-1,&bins_mjj[0]));
  histoB.push_back(new TH1F(Form("histo_%s_mjj_regionB",postfix.c_str()),"",bins_mjj.size()-1,&bins_mjj[0]));
  histoC.push_back(new TH1F(Form("histo_%s_mjj_regionC",postfix.c_str()),"",bins_mjj.size()-1,&bins_mjj[0]));
  histoD.push_back(new TH1F(Form("histo_%s_mjj_regionD",postfix.c_str()),"",bins_mjj.size()-1,&bins_mjj[0]));

  histoA.back()->Sumw2();
  histoB.back()->Sumw2();
  histoC.back()->Sumw2();
  histoD.back()->Sumw2();
}


///////////
static bool doSmoothing = true;
void makeBackgroundEstimation (Category category, string outputDIR, bool addEWKBackgrounds, bool isEOS){

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
  vector<TH1F*> histoQCD_A, histoQCD_B, histoQCD_C, histoQCD_D;
  cout<<"Run on QCD-MC "<<endl;
  bookHistogram(histoQCD_A,histoQCD_B,histoQCD_C,histoQCD_D,"QCD",category);
  fillHistograms(chainQCD,histoQCD_A,histoQCD_B,histoQCD_C,histoQCD_D,true,category);

  vector<TH1F*> histoQCD_A_v2, histoQCD_B_v2, histoQCD_C_v2, histoQCD_D_v2;
  cout<<"Run on QCD-MC: wider binning"<<endl;
  bookHistogram(histoQCD_A_v2,histoQCD_B_v2,histoQCD_C_v2,histoQCD_D_v2,"QCD_v2",category,2);
  fillHistograms(chainQCD,histoQCD_A_v2,histoQCD_B_v2,histoQCD_C_v2,histoQCD_D_v2,true,category);
  
  // Book WJets MC and Fill
  vector<TH1F*> histoWJets_A, histoWJets_B, histoWJets_C,histoWJets_D;
  vector<TH1F*> histoWJets_A_v2, histoWJets_B_v2, histoWJets_C_v2,histoWJets_D_v2;
  if(addEWKBackgrounds){
    cout<<"Run on WJets-MC "<<endl;
    bookHistogram(histoWJets_A,histoWJets_B,histoWJets_C,histoWJets_D,"WJets",category);
    fillHistograms(chainWJets,histoWJets_A,histoWJets_B,histoWJets_C,histoWJets_D,true,category);
    cout<<"Run on WJets-MC: wider binning "<<endl;
    bookHistogram(histoWJets_A_v2,histoWJets_B_v2,histoWJets_C_v2,histoWJets_D_v2,"WJets_v2",category,2);
    fillHistograms(chainWJets,histoWJets_A_v2,histoWJets_B_v2,histoWJets_C_v2,histoWJets_D_v2,true,category);
  }

  // Book ZJets MC and Fill
  vector<TH1F*> histoZJets_A, histoZJets_B, histoZJets_C, histoZJets_D;
  vector<TH1F*> histoZJets_A_v2, histoZJets_B_v2, histoZJets_C_v2, histoZJets_D_v2;
  if(addEWKBackgrounds){
    cout<<"Run on ZJets-MC "<<endl;
    bookHistogram(histoZJets_A,histoZJets_B,histoZJets_C,histoZJets_D,"ZJets",category);
    fillHistograms(chainZJets,histoZJets_A,histoZJets_B,histoZJets_C,histoZJets_D,true,category);
    cout<<"Run on ZJets-MC: wider binning "<<endl;
    bookHistogram(histoZJets_A_v2,histoZJets_B_v2,histoZJets_C_v2,histoZJets_D_v2,"ZJets_v2",category,2);
    fillHistograms(chainZJets,histoZJets_A_v2,histoZJets_B_v2,histoZJets_C_v2,histoZJets_D_v2,true,category);
  }
  
  
  // Book Data and Fill
  vector<TH1F*> histoData_A, histoData_B, histoData_C, histoData_D;
  cout<<"Run on Data "<<endl;
  bookHistogram(histoData_A,histoData_B,histoData_C,histoData_D,"Data",category);
  fillHistograms(chainData,histoData_A,histoData_B,histoData_C,histoData_D,false,category);

  cout<<"########################"<<endl;
  cout<<"W+jets rate in Region A: "<<histoWJets_A.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region A: "<<histoZJets_A.front()->Integral()<<endl;
  cout<<"QCD rate in Region A: "<<histoQCD_A.front()->Integral()<<endl;
  cout<<"Data rate in Region A: "<<histoData_A.front()->Integral()<<endl;

  cout<<"########################"<<endl;
  cout<<"W+jets rate in Region B: "<<histoWJets_B.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region B: "<<histoZJets_B.front()->Integral()<<endl;
  cout<<"QCD rate in Region B: "<<histoQCD_B.front()->Integral()<<endl;
  cout<<"Data rate in Region B: "<<histoData_B.front()->Integral()<<endl;

  cout<<"########################"<<endl;
  cout<<"W+jets rate in Region C: "<<histoWJets_C.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region C: "<<histoZJets_C.front()->Integral()<<endl;
  cout<<"QCD rate in Region C: "<<histoQCD_C.front()->Integral()<<endl;
  cout<<"Data rate in Region C: "<<histoData_C.front()->Integral()<<endl;

  cout<<"########################"<<endl;
  cout<<"W+jets rate in Region D: "<<histoWJets_D.front()->Integral()<<endl;
  cout<<"Z+jets rate in Region D: "<<histoZJets_D.front()->Integral()<<endl;
  cout<<"QCD rate in Region D: "<<histoQCD_D.front()->Integral()<<endl;
  cout<<"Data rate in Region D: "<<histoData_D.front()->Integral()<<endl;
  cout<<"########################"<<endl;

  if(doSmoothing){

    for(auto hist : histoQCD_A_v2) hist->Smooth();
    for(auto hist : histoQCD_B_v2) hist->Smooth();
    for(auto hist : histoQCD_C_v2) hist->Smooth();
    for(auto hist : histoQCD_D_v2) hist->Smooth();
    
  }

  // make transfer factors based on QCD MC: D/B and C/A
  cout<<"Make transfer factors "<<endl;

  vector<TH1F*> transferQCD_DB, transferQCD_CA;
  vector<TH1F*> transferQCD_DB_wjet_up; 
  vector<TH1F*> transferQCD_DB_wjet_dw; 
  vector<TH1F*> transferQCD_DB_zjet_up; 
  vector<TH1F*> transferQCD_DB_zjet_dw; 
  vector<TH1F*> transferQCD_CA_wjet_up; 
  vector<TH1F*> transferQCD_CA_wjet_dw; 
  vector<TH1F*> transferQCD_CA_zjet_up; 
  vector<TH1F*> transferQCD_CA_zjet_dw; 

  vector<TGraph*> transferQCD_DB_graph, transferQCD_CA_graph;

  vector<TGraph*> transferQCD_DB_graph_wjet_up; 
  vector<TGraph*> transferQCD_DB_graph_wjet_dw; 
  vector<TGraph*> transferQCD_DB_graph_zjet_up; 
  vector<TGraph*> transferQCD_DB_graph_zjet_dw; 
  vector<TGraph*> transferQCD_CA_graph_wjet_up; 
  vector<TGraph*> transferQCD_CA_graph_wjet_dw; 
  vector<TGraph*> transferQCD_CA_graph_zjet_up; 
  vector<TGraph*> transferQCD_CA_graph_zjet_dw; 

  vector<vector<TGraph*> > transferQCD_DB_graph_BinUp;
  vector<vector<TGraph*> > transferQCD_DB_graph_BinDw;
  vector<vector<TGraph*> > transferQCD_CA_graph_BinUp;
  vector<vector<TGraph*> > transferQCD_CA_graph_BinDw;


  for(size_t ihist = 0; ihist < histoQCD_B.size(); ihist++){
    transferQCD_DB_graph_BinUp.push_back(vector<TGraph*> ());
    transferQCD_DB_graph_BinDw.push_back(vector<TGraph*> ());
    transferQCD_CA_graph_BinUp.push_back(vector<TGraph*> ());
    transferQCD_CA_graph_BinDw.push_back(vector<TGraph*> ());
  }

  for(size_t ihist = 0; ihist < histoQCD_B.size(); ihist++){

    transferQCD_DB.push_back((TH1F*) histoQCD_D_v2.at(ihist)->Clone(Form("transferQCD_DB_%d",int(ihist))));
    transferQCD_CA.push_back((TH1F*) histoQCD_C.at(ihist)->Clone(Form("transferQCD_CA_%d",int(ihist))));

    TH1F* denominator_DB = (TH1F*) histoQCD_B_v2.at(ihist)->Clone(Form("denominator_DB_%d",int(ihist)));
    denominator_DB->Add(histoWJets_B_v2.at(ihist));
    denominator_DB->Add(histoZJets_B_v2.at(ihist));
    transferQCD_DB.back()->Divide(denominator_DB);

    TH1F* denominator_CA = (TH1F*) histoQCD_A.at(ihist)->Clone(Form("denominator_CA_%d",int(ihist)));
    denominator_CA->Add(histoWJets_A.at(ihist));
    denominator_CA->Add(histoZJets_A.at(ihist));
    transferQCD_CA.back()->Divide(denominator_CA);
    
    transferQCD_DB_graph.push_back(new TGraph(transferQCD_DB.back()));
    transferQCD_CA_graph.push_back(new TGraph(transferQCD_CA.back()));

    ///// ----------
    transferQCD_DB_wjet_up.push_back((TH1F*) histoQCD_D_v2.at(ihist)->Clone(Form("transferQCD_DB_wjet_up_%d",int(ihist))));
    transferQCD_DB_wjet_dw.push_back((TH1F*) histoQCD_D_v2.at(ihist)->Clone(Form("transferQCD_DB_wjet_dw_%d",int(ihist))));

    TH1F* denominator_DB_wjet_up = (TH1F*) histoQCD_B_v2.at(ihist)->Clone(Form("denominator_DB_wjet_up_%d",int(ihist)));    
    histoWJets_B_v2.at(ihist)->Scale(1+scale_wjet);
    denominator_DB_wjet_up->Add(histoWJets_B_v2.at(ihist));
    denominator_DB_wjet_up->Add(histoZJets_B_v2.at(ihist));
    transferQCD_DB_wjet_up.back()->Divide(denominator_DB_wjet_up);
    histoWJets_B_v2.at(ihist)->Scale(1./(1+scale_wjet));

    TH1F* denominator_DB_wjet_dw = (TH1F*) histoQCD_B_v2.at(ihist)->Clone(Form("denominator_DB_wjet_dw_%d",int(ihist)));
    histoWJets_B_v2.at(ihist)->Scale(1-scale_wjet);
    denominator_DB_wjet_dw->Add(histoWJets_B_v2.at(ihist));
    denominator_DB_wjet_dw->Add(histoZJets_B_v2.at(ihist));
    transferQCD_DB_wjet_dw.back()->Divide(denominator_DB_wjet_dw);
    histoWJets_B_v2.at(ihist)->Scale(1./(1-scale_wjet));

    ///// ----------
    transferQCD_DB_zjet_up.push_back((TH1F*) histoQCD_D_v2.at(ihist)->Clone(Form("transferQCD_DB_zjet_up_%d",int(ihist))));
    transferQCD_DB_zjet_dw.push_back((TH1F*) histoQCD_D_v2.at(ihist)->Clone(Form("transferQCD_DB_zjet_dw_%d",int(ihist))));

    TH1F* denominator_DB_zjet_up = (TH1F*) histoQCD_B_v2.at(ihist)->Clone(Form("denominator_DB_zjet_up_%d",int(ihist)));    
    histoZJets_B_v2.at(ihist)->Scale(1+scale_zjet);
    denominator_DB_zjet_up->Add(histoWJets_B_v2.at(ihist));
    denominator_DB_zjet_up->Add(histoZJets_B_v2.at(ihist));
    transferQCD_DB_zjet_up.back()->Divide(denominator_DB_wjet_up);
    histoZJets_B_v2.at(ihist)->Scale(1./(1+scale_zjet));

    TH1F* denominator_DB_zjet_dw = (TH1F*) histoQCD_B_v2.at(ihist)->Clone(Form("denominator_DB_zjet_dw_%d",int(ihist)));
    histoZJets_B_v2.at(ihist)->Scale(1-scale_zjet);
    denominator_DB_zjet_dw->Add(histoWJets_B_v2.at(ihist));
    denominator_DB_zjet_dw->Add(histoZJets_B_v2.at(ihist));
    transferQCD_DB_zjet_dw.back()->Divide(denominator_DB_wjet_dw);
    histoZJets_B_v2.at(ihist)->Scale(1./(1-scale_zjet));

    transferQCD_DB_graph_wjet_up.push_back(new TGraph(transferQCD_DB_wjet_up.back()));
    transferQCD_DB_graph_wjet_dw.push_back(new TGraph(transferQCD_DB_wjet_dw.back()));
    transferQCD_DB_graph_zjet_up.push_back(new TGraph(transferQCD_DB_zjet_up.back()));
    transferQCD_DB_graph_zjet_dw.push_back(new TGraph(transferQCD_DB_zjet_dw.back()));

    ///// ----------
    transferQCD_CA_wjet_up.push_back((TH1F*) histoQCD_C.at(ihist)->Clone(Form("transferQCD_CA_wjet_up_%d",int(ihist))));
    transferQCD_CA_wjet_dw.push_back((TH1F*) histoQCD_C.at(ihist)->Clone(Form("transferQCD_CA_wjet_dw_%d",int(ihist))));

    TH1F* denominator_CA_wjet_up = (TH1F*) histoQCD_A.at(ihist)->Clone(Form("denominator_CA_wjet_up_%d",int(ihist)));    
    histoWJets_A.at(ihist)->Scale(1+scale_wjet);
    denominator_CA_wjet_up->Add(histoWJets_A.at(ihist));
    denominator_CA_wjet_up->Add(histoZJets_A.at(ihist));
    transferQCD_CA_wjet_up.back()->Divide(denominator_CA_wjet_up);
    histoWJets_A.at(ihist)->Scale(1./(1+scale_wjet));

    TH1F* denominator_CA_wjet_dw = (TH1F*) histoQCD_A.at(ihist)->Clone(Form("denominator_CA_wjet_dw_%d",int(ihist)));
    histoWJets_A.at(ihist)->Scale(1-scale_wjet);
    denominator_CA_wjet_dw->Add(histoWJets_A.at(ihist));
    denominator_CA_wjet_dw->Add(histoZJets_A.at(ihist));
    transferQCD_CA_wjet_dw.back()->Divide(denominator_CA_wjet_dw);
    histoWJets_A.at(ihist)->Scale(1./(1-scale_wjet));

    ////
    transferQCD_CA_zjet_up.push_back((TH1F*) histoQCD_C.at(ihist)->Clone(Form("transferQCD_CA_zjet_up_%d",int(ihist))));
    transferQCD_CA_zjet_dw.push_back((TH1F*) histoQCD_C.at(ihist)->Clone(Form("transferQCD_CA_zjet_dw_%d",int(ihist))));

    TH1F* denominator_CA_zjet_up = (TH1F*) histoQCD_A.at(ihist)->Clone(Form("denominator_CA_zjet_up_%d",int(ihist)));    
    histoZJets_A.at(ihist)->Scale(1+scale_zjet);
    denominator_CA_zjet_up->Add(histoWJets_A.at(ihist));
    denominator_CA_zjet_up->Add(histoZJets_A.at(ihist));
    transferQCD_CA_zjet_up.back()->Divide(denominator_CA_wjet_up);
    histoZJets_A.at(ihist)->Scale(1./(1+scale_zjet));

    TH1F* denominator_CA_zjet_dw = (TH1F*) histoQCD_A.at(ihist)->Clone(Form("denominator_CA_zjet_dw_%d",int(ihist)));
    histoZJets_A.at(ihist)->Scale(1-scale_zjet);
    denominator_CA_zjet_dw->Add(histoWJets_A.at(ihist));
    denominator_CA_zjet_dw->Add(histoZJets_A.at(ihist));
    transferQCD_CA_zjet_dw.back()->Divide(denominator_CA_wjet_dw);
    histoZJets_A.at(ihist)->Scale(1./(1-scale_zjet));
    
    transferQCD_CA_graph_wjet_up.push_back(new TGraph(transferQCD_CA_wjet_up.back()));
    transferQCD_CA_graph_wjet_dw.push_back(new TGraph(transferQCD_CA_wjet_dw.back()));
    transferQCD_CA_graph_zjet_up.push_back(new TGraph(transferQCD_CA_zjet_up.back()));
    transferQCD_CA_graph_zjet_dw.push_back(new TGraph(transferQCD_CA_zjet_dw.back()));

    ////
    vector<TH1F*> transferDB_binUp;
    vector<TH1F*> transferDB_binDw;

    for(int iBin = 0; iBin < transferQCD_DB.back()->GetNbinsX(); iBin++){
      transferDB_binUp.push_back((TH1F*) transferQCD_DB.back()->Clone(Form("%s_bin%d_Up",transferQCD_DB.back()->GetName(),iBin)));
      transferDB_binUp.back()->SetBinContent(iBin+1,transferQCD_DB.back()->GetBinContent(iBin+1)+transferQCD_DB.back()->GetBinError(iBin+1));
      transferDB_binDw.push_back((TH1F*) transferQCD_DB.back()->Clone(Form("%s_bin%d_Dw",transferQCD_DB.back()->GetName(),iBin)));
      transferDB_binDw.back()->SetBinContent(iBin+1,transferQCD_DB.back()->GetBinContent(iBin+1)-transferQCD_DB.back()->GetBinError(iBin+1));      
    }

    for(auto hist : transferDB_binUp)
      transferQCD_DB_graph_BinUp.at(ihist).push_back(new TGraph(hist));
    for(auto hist : transferDB_binDw)
      transferQCD_DB_graph_BinDw.at(ihist).push_back(new TGraph(hist));

    vector<TH1F*> transferCA_binUp;
    vector<TH1F*> transferCA_binDw;
    for(int iBin = 0; iBin < transferQCD_CA.back()->GetNbinsX(); iBin++){
      transferCA_binUp.push_back((TH1F*) transferQCD_CA.back()->Clone(Form("%s_bin%d_Up",transferQCD_CA.back()->GetName(),iBin)));
      transferCA_binUp.back()->SetBinContent(iBin+1,transferQCD_CA.back()->GetBinContent(iBin+1)+transferQCD_CA.back()->GetBinError(iBin+1));
      transferCA_binDw.push_back((TH1F*) transferQCD_CA.back()->Clone(Form("%s_bin%d_Dw",transferQCD_CA.back()->GetName(),iBin)));
      transferCA_binDw.back()->SetBinContent(iBin+1,transferQCD_CA.back()->GetBinContent(iBin+1)+transferQCD_CA.back()->GetBinError(iBin+1));      
    }

    for(auto hist : transferCA_binUp)
      transferQCD_CA_graph_BinUp.at(ihist).push_back(new TGraph(hist));
    for(auto hist : transferCA_binDw)
      transferQCD_CA_graph_BinDw.at(ihist).push_back(new TGraph(hist));    
  }  

  // data TGraph
  vector<TGraph*> histoData_A_graph, histoData_B_graph, histoData_C_graph, histoData_D_graph;
  for(auto hist: histoData_A)
    histoData_A_graph.push_back(new TGraph(hist));
  for(auto hist: histoData_B)
    histoData_B_graph.push_back(new TGraph(hist));
  for(auto hist: histoData_C)
    histoData_C_graph.push_back(new TGraph(hist));
  for(auto hist: histoData_D)
    histoData_D_graph.push_back(new TGraph(hist));

  /// do the central value estimation
  cout<<"Perfom Background estimation "<<endl;
  vector<TH1F*> histoQCD_estimatedInC;
  vector<TH1F*> histoQCD_estimatedInD;  

  for(size_t ihist = 0; ihist < histoData_C.size(); ihist++){
    histoQCD_estimatedInC.push_back((TH1F*) histoData_C.at(ihist)->Clone(Form("histoQCD_estimatedInC_%d",int(ihist))));
    histoQCD_estimatedInD.push_back((TH1F*) histoData_D.at(ihist)->Clone(Form("histoQCD_estimatedInD_%d",int(ihist))));
    histoQCD_estimatedInC.back()->Reset();
    histoQCD_estimatedInD.back()->Reset();

    for(int iBin = 0; iBin < histoQCD_estimatedInC.back()->GetNbinsX(); iBin++){
      float binCenter = histoData_C.at(ihist)->GetBinCenter(iBin+1);
      histoQCD_estimatedInC.back()->SetBinContent(iBin+1,histoData_A_graph.at(ihist)->Eval(binCenter)*transferQCD_CA_graph.at(ihist)->Eval(binCenter));
    }
    for(int iBin = 0; iBin < histoQCD_estimatedInD.back()->GetNbinsX(); iBin++){
      float binCenter = histoData_D.at(ihist)->GetBinCenter(iBin+1);
      histoQCD_estimatedInD.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph.at(ihist)->Eval(binCenter));
    }
  }
  
  // Do bin-by-bin uncertainty on the total background from TFs
  vector< vector<TH1F*> > histoQCD_estimatedInD_binUp;
  vector< vector<TH1F*> > histoQCD_estimatedInD_binDw;

  for(size_t ihist = 0; ihist < histoQCD_estimatedInD.size(); ihist++){
    vector<TH1F*> vectorHistBinByBin_Up;
    vector<TH1F*> vectorHistBinByBin_Dw;
    for(int iBin = 0; iBin < histoQCD_estimatedInD.at(ihist)->GetNbinsX(); iBin++){
      vectorHistBinByBin_Up.push_back((TH1F*) histoQCD_estimatedInD.at(ihist)->Clone(Form("histoQCD_estimatedInD_%d_Bin_%d_Up",int(ihist),iBin)));
      vectorHistBinByBin_Dw.push_back((TH1F*) histoQCD_estimatedInD.at(ihist)->Clone(Form("histoQCD_estimatedInD_%d_Bin_%d_Dw",int(ihist),iBin)));

      float binCenter = histoData_D.at(ihist)->GetBinCenter(iBin+1);                                
      int binTransfer = transferQCD_DB.at(ihist)->FindBin(binCenter);
      vectorHistBinByBin_Up.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_BinUp.at(ihist).at(binTransfer-1)->Eval(binCenter));
      vectorHistBinByBin_Dw.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_BinDw.at(ihist).at(binTransfer-1)->Eval(binCenter));
    }
    histoQCD_estimatedInD_binUp.push_back(vectorHistBinByBin_Up);
    histoQCD_estimatedInD_binDw.push_back(vectorHistBinByBin_Dw);
  }

  vector<TH1F*> histoQCD_estimatedInD_zjet_up;
  vector<TH1F*> histoQCD_estimatedInD_zjet_dw;
  vector<TH1F*> histoQCD_estimatedInD_wjet_up;
  vector<TH1F*> histoQCD_estimatedInD_wjet_dw;

  //////////
  for(size_t ihist = 0; ihist < histoQCD_estimatedInD.size(); ihist++){

    histoQCD_estimatedInD_zjet_up.push_back((TH1F*) histoData_B.at(ihist)->Clone(Form("histoQCD_estimatedInD_zjet_up_%d",int(ihist))));
    histoQCD_estimatedInD_zjet_dw.push_back((TH1F*) histoData_B.at(ihist)->Clone(Form("histoQCD_estimatedInD_zjet_dw_%d",int(ihist))));
    histoQCD_estimatedInD_wjet_up.push_back((TH1F*) histoData_B.at(ihist)->Clone(Form("histoQCD_estimatedInD_wjet_up_%d",int(ihist))));
    histoQCD_estimatedInD_wjet_dw.push_back((TH1F*) histoData_B.at(ihist)->Clone(Form("histoQCD_estimatedInD_wjet_dw_%d",int(ihist))));

    histoQCD_estimatedInD_zjet_up.back()->Reset();
    histoQCD_estimatedInD_zjet_dw.back()->Reset();
    histoQCD_estimatedInD_wjet_up.back()->Reset();
    histoQCD_estimatedInD_wjet_dw.back()->Reset();

    for(int iBin = 0; iBin < histoQCD_estimatedInD.at(ihist)->GetNbinsX(); iBin++){
      float binCenter = histoData_B.at(ihist)->GetBinCenter(iBin+1);      
      histoQCD_estimatedInD_zjet_up.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_zjet_up.at(ihist)->Eval(binCenter));
      histoQCD_estimatedInD_zjet_dw.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_zjet_dw.at(ihist)->Eval(binCenter));
      histoQCD_estimatedInD_wjet_up.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_wjet_up.at(ihist)->Eval(binCenter));
      histoQCD_estimatedInD_wjet_dw.back()->SetBinContent(iBin+1,histoData_B_graph.at(ihist)->Eval(binCenter)*transferQCD_DB_graph_wjet_dw.at(ihist)->Eval(binCenter));
    }
  }

  
  cout<<"Save and close the file"<<endl;

  TFile* output = new TFile((outputDIR+"/distributions_qcd.root").c_str(),"RECREATE");
  output->mkdir("Data");
  output->cd("Data");
  for(auto hist : histoData_A) hist->Write();
  for(auto hist : histoData_B) hist->Write();
  for(auto hist : histoData_C) hist->Write();
  output->cd();
  output->mkdir("QCD_MC");
  output->cd("QCD_MC");
  for(auto hist : histoQCD_A) hist->Write();
  for(auto hist : histoQCD_B) hist->Write();
  for(auto hist : histoQCD_C) hist->Write();
  for(auto hist : histoQCD_D) hist->Write();
  output->cd();
  output->mkdir("QCD_MC_v2");
  output->cd("QCD_MC_v2");
  for(auto hist : histoQCD_A_v2) hist->Write();
  for(auto hist : histoQCD_B_v2) hist->Write();
  for(auto hist : histoQCD_C_v2) hist->Write();
  for(auto hist : histoQCD_D_v2) hist->Write();
  output->cd();
  output->mkdir("WJets_MC");
  output->cd("WJets_MC");
  for(auto hist : histoWJets_A) hist->Write();
  for(auto hist : histoWJets_B) hist->Write();
  for(auto hist : histoWJets_C) hist->Write();
  for(auto hist : histoWJets_D) hist->Write();
  output->cd();
  output->mkdir("WJets_MC_v2");
  output->cd("WJets_MC_v2");
  for(auto hist : histoWJets_A_v2) hist->Write();
  for(auto hist : histoWJets_B_v2) hist->Write();
  for(auto hist : histoWJets_C_v2) hist->Write();
  for(auto hist : histoWJets_D_v2) hist->Write();
  output->cd();
  output->mkdir("ZJets_MC");
  output->cd("ZJets_MC");
  for(auto hist : histoZJets_A) hist->Write();
  for(auto hist : histoZJets_B) hist->Write();
  for(auto hist : histoZJets_C) hist->Write();
  for(auto hist : histoZJets_D) hist->Write();
  output->cd();
  output->mkdir("ZJets_MC_v2");
  output->cd("ZJets_MC_v2");
  for(auto hist : histoZJets_A_v2) hist->Write();
  for(auto hist : histoZJets_B_v2) hist->Write();
  for(auto hist : histoZJets_C_v2) hist->Write();
  for(auto hist : histoZJets_D_v2) hist->Write();
  output->cd();
  output->mkdir("Transfer_QCD");
  output->cd("Transfer_QCD");
  for(auto hist : transferQCD_DB) hist->Write();
  for(auto hist : transferQCD_CA) hist->Write();
  for(auto hist : transferQCD_DB_graph) hist->Write(Form("%s_graph",hist->GetName()));
  for(auto hist : transferQCD_CA_graph) hist->Write(Form("%s_graph",hist->GetName()));
  output->mkdir("Transfer_DB_BinByBin");
  output->cd("Transfer_DB_BinByBin");
  for(auto vec : transferQCD_DB_graph_BinUp){ for(auto graph : vec) graph->Write(); }
  for(auto vec : transferQCD_DB_graph_BinDw){ for(auto graph : vec) graph->Write(); }
  for(auto vec : transferQCD_CA_graph_BinUp){ for(auto graph : vec) graph->Write(); }
  for(auto vec : transferQCD_CA_graph_BinDw){ for(auto graph : vec) graph->Write(); }
  output->cd();
  output->mkdir("Transfer_ZJets_up");
  output->cd("Transfer_ZJets_up");
  for(auto graph : transferQCD_CA_graph_zjet_up) graph->Write();
  for(auto graph : transferQCD_DB_graph_zjet_up) graph->Write();
  output->cd();
  output->mkdir("Transfer_ZJets_dw");
  output->cd("Transfer_ZJets_dw");
  for(auto graph : transferQCD_CA_graph_zjet_dw) graph->Write();
  for(auto graph : transferQCD_DB_graph_zjet_dw) graph->Write();
  output->cd();
  output->mkdir("Transfer_WJets_up");
  output->cd("Transfer_WJets_up");
  for(auto graph : transferQCD_CA_graph_wjet_up) graph->Write();
  for(auto graph : transferQCD_DB_graph_wjet_up) graph->Write();
  output->cd();
  output->mkdir("Transfer_WJets_dw");
  output->cd("Transfer_WJets_dw");
  for(auto graph : transferQCD_CA_graph_wjet_dw) graph->Write();
  for(auto graph : transferQCD_DB_graph_wjet_dw) graph->Write();
  output->cd();
  output->mkdir("Estimation");
  output->cd("Estimation");
  for(auto hist : histoQCD_estimatedInC) hist->Write(Form("%s_graph",hist->GetName()));
  for(auto hist : histoQCD_estimatedInD) hist->Write(Form("%s_graph",hist->GetName()));
  output->mkdir("BinByBin");
  output->cd("BinByBin");
  for(auto vec : histoQCD_estimatedInD_binUp){ for (auto hist : vec) hist->Write();} 
  for(auto vec : histoQCD_estimatedInD_binDw){ for (auto hist : vec) hist->Write();} 
  output->cd();  
  output->mkdir("ZJets");
  output->cd("ZJets");  
  for(auto hist : histoQCD_estimatedInD_zjet_up) hist->Write();
  for(auto hist : histoQCD_estimatedInD_zjet_dw) hist->Write();
  output->cd();
  output->mkdir("WJets");  
  output->cd("WJets");  
  for(auto hist : histoQCD_estimatedInD_wjet_up) hist->Write();
  for(auto hist : histoQCD_estimatedInD_wjet_dw) hist->Write();
  output->cd();  
  output->cd();
  output->Close();
  
}
