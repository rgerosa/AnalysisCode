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

void fillHistograms(TChain* chain, vector<TH1F*> histoA, vector<TH1F*> histoB, vector<TH1F*> histoC, const bool & isMC, const Category & category){


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

    if(ht < 350) continue;

    double evtweight = 1 ;
    if(isMC)
      evtweight = *xsec*luminosity*(*wgt)*(*wgtpileup)*(*wgtbtag)/(**wgtsum);
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
	  hist->Fill(jet2.Eta(),evtweight);
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
	  hist->Fill(jet2.Eta(),evtweight);
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
	  hist->Fill(jet2.Eta(),evtweight);
      }
    }
  }
  cout<<endl;
  
}

void plotDataMC(TH1* histoData, TH1* histoQCD, TCanvas* canvas, const string & outputDIR){

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

  histoData->GetYaxis()->SetRangeUser(min(histoQCD->GetMinimum(),histoData->GetMinimum())*0.1,max(histoQCD->GetMaximum(),histoData->GetMaximum())*100);
  histoData->GetYaxis()->SetTitle("Events");
  histoData->GetYaxis()->SetTitleOffset(1.10);
  histoData->GetYaxis()->SetLabelSize(0.035);

  TH1* nhist = (TH1*) histoData->Clone("nhist");
  TH1* dhist = (TH1*) histoQCD->Clone("dhist");
  TH1* dhist_2 = (TH1*) histoQCD->Clone("dhist_2");
  TH1* unhist = (TH1*) histoData->Clone("unhist");

  histoData->GetXaxis()->SetTitleSize(0);
  histoData->GetXaxis()->SetLabelSize(0);

  histoQCD->SetFillColor(kAzure+1);
  histoQCD->SetFillStyle(1001);
  histoQCD->SetLineColor(kBlack);
  histoQCD->Draw("hist same");
  histoData->Draw("same");

  TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histoData,"Data","EP");
  leg->AddEntry(histoQCD,"QCD MC","F");
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
 
  TString name (histoData->GetName());
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

  canvas->SetLogy();
  
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

  if(leg) delete leg;
  if(nhist) delete nhist;
  if(dhist) delete dhist;
  if(dhist_2) delete dhist_2;
  if(unhist) delete unhist;
  if(pad2) delete pad2;
  
}

///////////
void makeDataMCComparison (string inputPathMC, string inputPathData, string outputDIR, Category category, bool isEOS){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  initializeBinning();
  gROOT->SetBatch(kTRUE);
  
  TChain* chainMC = new TChain("tree/tree");  
  TChain* chainData = new TChain("tree/tree");  

  loadChain(inputPathMC,chainMC,isEOS);
  loadChain(inputPathData,chainData,isEOS);

  // variables --> we can just check agreement in region A, B and C
  vector<double> bins_mjj  = selectBinning("mjj",category);
  vector<double> bins_jeteta  = selectBinning("jeteta",category);
  vector<double> bins_jeteta2 = selectBinning("jeteta2",category);
  vector<double> bins_ht   = selectBinning("HT",category);
  vector<double> bins_lowMet  = {100.,120.,140.,160.,180,200.,225.,250.};
  vector<double> bins_highMet = {250.,300.,350.,500.,1000};


  vector<TH1F*> histoQCD_A, histoQCD_B, histoQCD_C;

  histoQCD_A.push_back(new TH1F("histoQCD_mjj_regionA","",bins_mjj.size()-1,&bins_mjj[0]));
  histoQCD_B.push_back(new TH1F("histoQCD_mjj_regionB","",bins_mjj.size()-1,&bins_mjj[0]));
  histoQCD_C.push_back(new TH1F("histoQCD_mjj_regionC","",bins_mjj.size()-1,&bins_mjj[0]));
  histoQCD_A.back()->Sumw2();
  histoQCD_B.back()->Sumw2();
  histoQCD_C.back()->Sumw2();

  histoQCD_A.push_back(new TH1F("histoQCD_ht_regionA","",bins_ht.size()-1,&bins_ht[0]));
  histoQCD_B.push_back(new TH1F("histoQCD_ht_regionB","",bins_ht.size()-1,&bins_ht[0]));
  histoQCD_C.push_back(new TH1F("histoQCD_ht_regionC","",bins_ht.size()-1,&bins_ht[0]));
  histoQCD_A.back()->Sumw2();
  histoQCD_B.back()->Sumw2();
  histoQCD_C.back()->Sumw2();
  
  histoQCD_A.push_back(new TH1F("histoQCD_met_regionA","",bins_lowMet.size()-1,&bins_lowMet[0]));
  histoQCD_B.push_back(new TH1F("histoQCD_met_regionB","",bins_highMet.size()-1,&bins_highMet[0]));
  histoQCD_C.push_back(new TH1F("histoQCD_met_regionC","",bins_lowMet.size()-1,&bins_lowMet[0]));
  histoQCD_A.back()->Sumw2();
  histoQCD_B.back()->Sumw2();
  histoQCD_C.back()->Sumw2();

  histoQCD_A.push_back(new TH1F("histoQCD_jeteta_regionA","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoQCD_B.push_back(new TH1F("histoQCD_jeteta_regionB","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoQCD_C.push_back(new TH1F("histoQCD_jeteta_regionC","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoQCD_A.back()->Sumw2();
  histoQCD_B.back()->Sumw2();
  histoQCD_C.back()->Sumw2();

  histoQCD_A.push_back(new TH1F("histoQCD_jeteta2_regionA","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoQCD_B.push_back(new TH1F("histoQCD_jeteta2_regionB","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoQCD_C.push_back(new TH1F("histoQCD_jeteta2_regionC","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoQCD_A.back()->Sumw2();
  histoQCD_B.back()->Sumw2();
  histoQCD_C.back()->Sumw2();

  // fill histograms
  fillHistograms(chainMC,histoQCD_A,histoQCD_B,histoQCD_C,true,category);

  vector<TH1F*> histoData_A, histoData_B, histoData_C;

  histoData_A.push_back(new TH1F("histoData_mjj_regionA","",bins_mjj.size()-1,&bins_mjj[0]));
  histoData_B.push_back(new TH1F("histoData_mjj_regionB","",bins_mjj.size()-1,&bins_mjj[0]));
  histoData_C.push_back(new TH1F("histoData_mjj_regionC","",bins_mjj.size()-1,&bins_mjj[0]));
  histoData_A.back()->Sumw2();
  histoData_B.back()->Sumw2();
  histoData_C.back()->Sumw2();

  histoData_A.push_back(new TH1F("histoData_ht_regionA","",bins_ht.size()-1,&bins_ht[0]));
  histoData_B.push_back(new TH1F("histoData_ht_regionB","",bins_ht.size()-1,&bins_ht[0]));
  histoData_C.push_back(new TH1F("histoData_ht_regionC","",bins_ht.size()-1,&bins_ht[0]));
  histoData_A.back()->Sumw2();
  histoData_B.back()->Sumw2();
  histoData_C.back()->Sumw2();
  
  histoData_A.push_back(new TH1F("histoData_met_regionA","",bins_lowMet.size()-1,&bins_lowMet[0]));
  histoData_B.push_back(new TH1F("histoData_met_regionB","",bins_highMet.size()-1,&bins_highMet[0]));
  histoData_C.push_back(new TH1F("histoData_met_regionC","",bins_lowMet.size()-1,&bins_lowMet[0]));
  histoData_A.back()->Sumw2();
  histoData_B.back()->Sumw2();
  histoData_C.back()->Sumw2();

  histoData_A.push_back(new TH1F("histoData_jeteta_regionA","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoData_B.push_back(new TH1F("histoData_jeteta_regionB","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoData_C.push_back(new TH1F("histoData_jeteta_regionC","",bins_jeteta.size()-1,&bins_jeteta[0]));
  histoData_A.back()->Sumw2();
  histoData_B.back()->Sumw2();
  histoData_C.back()->Sumw2();

  histoData_A.push_back(new TH1F("histoData_jeteta2_regionA","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoData_B.push_back(new TH1F("histoData_jeteta2_regionB","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoData_C.push_back(new TH1F("histoData_jeteta2_regionC","",bins_jeteta2.size()-1,&bins_jeteta2[0]));
  histoData_A.back()->Sumw2();
  histoData_B.back()->Sumw2();
  histoData_C.back()->Sumw2();

  // fill histograms
  fillHistograms(chainData,histoData_A,histoData_B,histoData_C,false,category);
  
  //////
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  for(size_t ihist = 0; ihist < histoData_A.size(); ihist++){
    plotDataMC(histoData_A.at(ihist),histoQCD_A.at(ihist),canvas,outputDIR);
  }

  for(size_t ihist = 0; ihist < histoData_B.size(); ihist++){
    plotDataMC(histoData_B.at(ihist),histoQCD_B.at(ihist),canvas,outputDIR);
  }

  for(size_t ihist = 0; ihist < histoData_C.size(); ihist++){
    plotDataMC(histoData_C.at(ihist),histoQCD_C.at(ihist),canvas,outputDIR);
  }
  
}
