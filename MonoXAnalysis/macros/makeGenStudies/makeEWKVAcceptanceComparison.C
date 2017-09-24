#include "../CMS_lumi.h"

static float bosonPt        = 200;
static float bosonPtMax     = 80000;
static float mjj            = 500;
static float detajj         = 2.5;
static float dphijj         = 1.5;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;
static float minDeltaPhiVBF = 0.5;
static bool  removetaus     = true;
enum class Sample {znn, zll, wjet};
float lumi_ = 12.9;

/////////////////////////
void plotResults(TH1* histo_ewk, TH1* histo_qcd,const string & postfix,const string & xAxisLabel, const string & outputDIR, const Sample & sample){

  
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,650);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  TPad *pad1 = NULL;
  TPad *pad2 = NULL;

  histo_ewk->SetTitle("");
  histo_qcd->GetYaxis()->SetTitle("");
  histo_ewk->GetYaxis()->SetTitle("a.u.");

  histo_ewk->Scale(1./histo_ewk->Integral());
  histo_qcd->Scale(1./histo_qcd->Integral());

  histo_ewk->SetLineColor(kRed);
  histo_ewk->SetLineWidth(2);
  histo_ewk->SetMarkerStyle(20);
  histo_ewk->SetMarkerSize(1);
  histo_ewk->SetMarkerColor(kRed);

  histo_qcd->SetLineColor(kBlack);
  histo_qcd->SetLineWidth(2);
  histo_qcd->SetMarkerColor(kBlack);
  histo_qcd->SetMarkerSize(0.8);
  histo_qcd->SetMarkerStyle(20);

  histo_ewk->GetYaxis()->SetRangeUser(min(histo_ewk->GetMinimum(),histo_qcd->GetMinimum())*0.75,max(histo_ewk->GetMaximum(),histo_qcd->GetMaximum())*1.25);

  histo_ewk->GetXaxis()->SetTitle(xAxisLabel.c_str());
  histo_ewk->Draw("hist");
  CMS_lumi(canvas,"",true);    
  histo_qcd->Draw("hist same");

  TLegend* leg = new TLegend(0.7,0.70,0.90,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(sample == Sample::znn or sample == Sample::zll){
    leg->AddEntry(histo_qcd, "Z+jets QCD","L");
    leg->AddEntry(histo_ewk, "Z+jets EWK","L");
  }
  else if(sample == Sample::wjet){
    leg->AddEntry(histo_qcd, "W+jets QCD","L");
    leg->AddEntry(histo_ewk, "W+jets EWK","L");
  }
  
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".pdf").c_str(),"pdf");

  canvas->SetLogy();

  histo_ewk->GetYaxis()->SetRangeUser(min(histo_ewk->GetMinimum(),histo_qcd->GetMinimum())*0.1,max(histo_ewk->GetMaximum(),histo_qcd->GetMaximum())*500);
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+"_log.pdf").c_str(),"pdf");
  
}


///////////////////
void makeEWKVAcceptanceComparison(string inputDIR_QCD, string inputDIR_EWK, Sample sample, string outputDIR, bool applyKfactor = false, bool applyVBFCuts = false){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  
  //////
  vector<TH1*> khists;  
  TFile* kffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root","READ");

  if(applyKfactor){
    TH1* znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
    TH1* zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
    TH1* zewkhist = (TH1*) kffile->Get("EWKcorr/Z");
    if(zewkhist)
      zewkhist->Divide(znlohist);
    if(znlohist)
      znlohist->Divide(zlohist);

    TH1* wnlohist = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
    TH1* wlohist  = (TH1*) kffile->Get("WJets_LO/inv_pt");
    TH1* wewkhist = (TH1*) kffile->Get("EWKcorr/W");
    if(wewkhist)
      wewkhist->Divide(wnlohist);
    if(wnlohist)
      wnlohist->Divide(wlohist);

    if(sample == Sample::zll){
      khists.push_back(znlohist);
      khists.push_back(zewkhist);
    }
    else if(sample == Sample::wjet){
      khists.push_back(wnlohist);
      khists.push_back(wewkhist);
    }
  }

  ///////
  system(("mkdir -p "+outputDIR).c_str());
  vector<TTree*> tree_QCD;
  vector<TTree*> tree_EWK;
  vector<TFile*> file_QCD;
  vector<TFile*> file_EWK;

  system(("ls "+inputDIR_QCD+" | grep root > list.txt").c_str());
  ifstream infileQCD ("list.txt");
  string line;
  if(infileQCD.is_open()){
    while(!infileQCD.eof()){
      getline(infileQCD,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_QCD.push_back(TFile::Open((inputDIR_QCD+"/"+line).c_str(),"READ"));
      tree_QCD.push_back((TTree*) file_QCD.back()->Get("gentree/tree"));
    }
  }
  infileQCD.close();

  system(("ls "+inputDIR_EWK+" | grep root > list.txt").c_str());
  ifstream infileEWK ("list.txt");
  if(infileEWK.is_open()){
    while(!infileEWK.eof()){
      getline(infileEWK,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_EWK.push_back(TFile::Open((inputDIR_EWK+"/"+line).c_str(),"READ"));
      tree_EWK.push_back((TTree*) file_EWK.back()->Get("gentree/tree"));
    }
  }
  infileEWK.close();
  
  system("rm list.txt");

  // calculate sumwgt                                                                                                                                                                                 
  vector<double> sumwgt_qcd;
  int ifile = 0;
  double eventsNoBoson = 0;
  for(auto tree : tree_QCD){
    eventsNoBoson = 0;
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    cout<<"Calculate sumwgt for QCD file "<<file_QCD.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
    while(reader.Next()){
      // filter away bad events with no matching                                                                                                                                                      
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23){
	eventsNoBoson++;
      }
      else if(sample == Sample::wjet and fabs(*wzid) != 24){
	eventsNoBoson++;
      }
      
      sumwgt += *wgt;
    }
    cout<<"Tree QCD with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
    sumwgt_qcd.push_back(sumwgt);
    ifile++;
  }

  vector<double> sumwgt_ewk;
  ifile = 0;
  for(auto tree : tree_EWK){
    eventsNoBoson = 0;
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    cout<<"Calculate sumwgt for EWK file "<<file_EWK.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
    while(reader.Next()){

      // filter away bad events with no matching                                                                                                                                                      
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23){
	eventsNoBoson++;
      }
      else if(sample == Sample::wjet and fabs(*wzid) != 24){
	eventsNoBoson++;
      }

      sumwgt += *wgt;
    }
    cout<<"Tree EWK with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
    sumwgt_ewk.push_back(sumwgt);
    ifile++;
  }

  // distributions                                                                                                                                                                                    
  TH1F* bosonPt_QCD = new TH1F("bosonPt_QCD","",30,bosonPt,1300);
  TH1F* bosonPt_EWK = new TH1F("bosonPt_EWK","",30,bosonPt,1300);
  bosonPt_QCD->Sumw2();
  bosonPt_EWK->Sumw2();

  TH1F* leptonEta_QCD = new TH1F("leptonEta_QCD","",30,-4,4);
  TH1F* leptonEta_EWK = new TH1F("leptonEta_EWK","",30,-4,4);
  leptonEta_QCD->Sumw2();
  leptonEta_EWK->Sumw2();

  TH1F* leptonEta_2_QCD = new TH1F("leptonEta_2_QCD","",30,-4,4);
  TH1F* leptonEta_2_EWK = new TH1F("leptonEta_2_EWK","",30,-4,4);
  leptonEta_2_QCD->Sumw2();
  leptonEta_2_EWK->Sumw2();

  TH1F* leptonPt_QCD = new TH1F("leptonPt_QCD","",30,0,400);
  TH1F* leptonPt_EWK = new TH1F("leptonPt_EWK","",30,0,400);
  leptonPt_QCD->Sumw2();
  leptonPt_EWK->Sumw2();

  TH1F* leptonPt_2_QCD = new TH1F("leptonPt_2_QCD","",30,0,400);
  TH1F* leptonPt_2_EWK = new TH1F("leptonPt_2_EWK","",30,0,400);
  leptonPt_2_QCD->Sumw2();
  leptonPt_2_EWK->Sumw2();

  TH1F* mjj_QCD = NULL;
  TH1F* mjj_EWK = NULL;
  if(applyVBFCuts){
    mjj_EWK = new TH1F("mjj_EWK","",30,mjj,3000);
    mjj_QCD = new TH1F("mjj_QCD","",30,mjj,3000);
  }
  else{
    mjj_EWK = new TH1F("mjj_EWK","",30,0,3000);
    mjj_QCD = new TH1F("mjj_QCD","",30,0,3000);
  }
  mjj_QCD->Sumw2();
  mjj_EWK->Sumw2();

  TH1F* detajj_QCD = NULL;
  TH1F* detajj_EWK = NULL;
  if(applyVBFCuts){
    detajj_EWK = new TH1F("detajj_EWK","",15,detajj,10);
    detajj_QCD = new TH1F("detajj_QCD","",15,detajj,10);
  }
  else{
    detajj_EWK = new TH1F("detajj_EWK","",15,0,10);
    detajj_QCD = new TH1F("detajj_QCD","",15,0,10);
  }
  detajj_QCD->Sumw2();
  detajj_EWK->Sumw2();

  TH1F* jetmetdphi_QCD = new TH1F("jetmetdphi_QCD","",20,1,3.14);
  TH1F* jetmetdphi_EWK = new TH1F("jetmetdphi_EWK","",20,1,3.14);
  jetmetdphi_QCD->Sumw2();
  jetmetdphi_EWK->Sumw2();

  TH1F* jetpt1_QCD = new TH1F("jetpt1_QCD","",30,leadingJetVBF,800);
  TH1F* jetpt1_EWK = new TH1F("jetpt1_EWK","",30,leadingJetVBF,800);
  jetpt1_QCD->Sumw2();
  jetpt1_EWK->Sumw2();

  TH1F* jetpt2_QCD = new TH1F("jetpt2_QCD","",25,trailingJetVBF,500);
  TH1F* jetpt2_EWK = new TH1F("jetpt2_EWK","",25,trailingJetVBF,500);
  jetpt2_QCD->Sumw2();
  jetpt2_EWK->Sumw2();

  /// start making selections
  ifile = 0;
  for(auto tree: tree_QCD){

    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    TTreeReaderValue<int>   l1id  (reader,"l1id");
    TTreeReaderValue<int>   l2id  (reader,"l2id");
    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");

    TTreeReaderValue<float> l1pt   (reader,"l1pt");
    TTreeReaderValue<float> l1eta  (reader,"l1eta");
    TTreeReaderValue<float> l1phi  (reader,"l1phi");
    TTreeReaderValue<float> l2pt   (reader,"l2pt");
    TTreeReaderValue<float> l2eta  (reader,"l2eta");
    TTreeReaderValue<float> l2phi  (reader,"l2phi");

    TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
    TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");
    TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");
    TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");


    cout<<"Loop on QCD file "<<file_QCD.at(ifile)->GetName()<<endl;
    while(reader.Next()){

      float mindphi = 100;
      for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
        if(ijet > 3) break; // limiting min dphi to first 4 leading jets                                                                                                                     
        float dphi = fabs(*wzphi-jetphi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          mindphi = dphi;
      }

      if(*wzpt < bosonPt) continue;
      if(*wzpt > bosonPtMax) continue;

      Double_t kwgt  = 1.0;
      double   genpt = *wzpt;
      
      for (size_t i = 0; i < khists.size(); i++) {
	if (khists[i]) {
	  if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	  if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	  kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
	}
      }

      if(removetaus and fabs(*l1id) == 15) continue;

      // apply dijet selections
      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > minDeltaPhiVBF){

        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));

	if(fabs(jeteta->at(0)) > 3 and fabs(jeteta->at(1)) > 3) continue;
	if(fabs(jet1.DeltaPhi(jet2) < dphijj) continue;

	if(applyVBFCuts){
	  if((jet1+jet2).M() < mjj) continue;
	  if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
	}

	//////////////////
	bosonPt_QCD->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetpt1_QCD->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetpt2_QCD->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	mjj_QCD->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	detajj_QCD->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetmetdphi_QCD->Fill(mindphi,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	if(sample == Sample::wjet){
	  leptonPt_QCD->Fill(*l1pt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	  leptonEta_QCD->Fill(*l1eta,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	}
	else{
	  if((*l1pt) > (*l2pt)){
	    leptonPt_QCD->Fill(*l1pt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonEta_QCD->Fill(*l1eta,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonPt_2_QCD->Fill(*l2pt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonEta_2_QCD->Fill(*l2eta,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	  }
	  else{
	    leptonPt_QCD->Fill(*l2pt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonEta_QCD->Fill(*l2eta,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonPt_2_QCD->Fill(*l1pt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	    leptonEta_2_QCD->Fill(*l1eta,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	  }
	}
      }
    }
    ifile++;
  }

  /// start making selections
  ifile = 0;
  for(auto tree: tree_EWK){

    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    TTreeReaderValue<int>   l1id  (reader,"l1id");
    TTreeReaderValue<int>   l2id  (reader,"l2id");
    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");

    TTreeReaderValue<float> l1pt   (reader,"l1pt");
    TTreeReaderValue<float> l1eta  (reader,"l1eta");
    TTreeReaderValue<float> l1phi  (reader,"l1phi");
    TTreeReaderValue<float> l2pt   (reader,"l2pt");
    TTreeReaderValue<float> l2eta  (reader,"l2eta");
    TTreeReaderValue<float> l2phi  (reader,"l2phi");

    TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
    TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");
    TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");
    TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");

    cout<<"Loop on EWK file "<<file_EWK.at(ifile)->GetName()<<endl;
    while(reader.Next()){
      
      float mindphi = 100;
      for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
        if(ijet > 3) break; // limiting min dphi to first 4 leading jets                                                                                                                     
        float dphi = fabs(*wzphi-jetphi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          mindphi = dphi;
      }

      if(*wzpt < bosonPt) continue;
      if(*wzpt > bosonPtMax) continue;

      if(removetaus and fabs(*l1id) == 15) continue;

      // apply dijet selections
      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > minDeltaPhiVBF){
	
        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));

	
	if(applyVBFCuts){
	  if((jet1+jet2).M() < mjj) continue;
	  if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
	}

	//////////////////
	bosonPt_EWK->Fill(*wzpt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetpt1_EWK->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetpt2_EWK->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	mjj_EWK->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	detajj_EWK->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetmetdphi_EWK->Fill(mindphi,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	if(sample == Sample::wjet){
	  leptonPt_EWK->Fill(*l1pt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	  leptonEta_EWK->Fill(*l1eta,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	}
	else{
	  if((*l1pt) > (*l2pt)){
	    leptonPt_EWK->Fill(*l1pt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonEta_EWK->Fill(*l1eta,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonPt_2_EWK->Fill(*l2pt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonEta_2_EWK->Fill(*l2eta,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	  }
	  else{
	    leptonPt_EWK->Fill(*l2pt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonEta_EWK->Fill(*l2eta,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonPt_2_EWK->Fill(*l1pt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	    leptonEta_2_EWK->Fill(*l1eta,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	  }
	}
	///////////////////
      }
    }
    ifile++;
  }
  ///////// plotting results
  plotResults(bosonPt_EWK,bosonPt_QCD,"bosonPt","Boson p_{T} [GeV]",outputDIR,sample);
  ///////// plotting results
  plotResults(jetpt1_EWK,jetpt1_QCD,"jetpt1","p_{T}^{j1} [GeV]",outputDIR,sample);
  ///////// plotting results
  plotResults(jetpt2_EWK,jetpt2_QCD,"jetpt2","p_{T}^{j2} [GeV]",outputDIR,sample);
  ///////// plotting results
  plotResults(mjj_EWK,mjj_QCD,"mjj","m_{jj} [GeV]",outputDIR,sample);
  ///////// plotting results
  plotResults(detajj_EWK,detajj_QCD,"detajj","#Delta#eta_{jj}",outputDIR,sample);
  ///////// plotting results
  plotResults(jetmetdphi_EWK,jetmetdphi_QCD,"jetmetdphi","#Delta#phi(jet,met)",outputDIR,sample);
  /////////
  plotResults(leptonPt_EWK,leptonPt_QCD,"leptonpt","p_{T}^{lep} [GeV]",outputDIR,sample);
  /////////
  plotResults(leptonEta_EWK,leptonEta_QCD,"leptoneta","#eta_{lep}",outputDIR,sample);
  ////
  if(sample != Sample::wjet){
    plotResults(leptonPt_2_EWK,leptonPt_2_QCD,"leptonpt_2","p_{T}^{lep} [GeV]",outputDIR,sample);
    plotResults(leptonEta_2_EWK,leptonEta_2_QCD,"leptoneta_2","#eta_{lep}",outputDIR,sample);    
  }
  
  
  // acceptance:
  double integralInAcceptance  = 0;
  double integralOutAcceptance = 0;

  for(int iBin = 0; iBin < leptonEta_QCD->GetNbinsX()+1; iBin++){
    if(fabs(leptonEta_QCD->GetBinLowEdge(iBin+1)) > 3.2)
      integralOutAcceptance += leptonEta_QCD->GetBinContent(iBin+1);
    else
      integralInAcceptance += leptonEta_QCD->GetBinContent(iBin+1);
  }
  
  cout<<"Integral in acceptance for Lepton 1  V-QCD "<<integralInAcceptance/leptonEta_QCD->Integral()<<endl;
  cout<<"Integral out acceptance for Lepton 1 V-QCD "<<integralOutAcceptance/leptonEta_QCD->Integral()<<endl;

  integralInAcceptance  = 0;
  integralOutAcceptance = 0;

  for(int iBin = 0; iBin < leptonEta_EWK->GetNbinsX()+1; iBin++){
    if(fabs(leptonEta_EWK->GetBinLowEdge(iBin+1)) > 3.2)
      integralOutAcceptance += leptonEta_EWK->GetBinContent(iBin+1);
    else
      integralInAcceptance += leptonEta_EWK->GetBinContent(iBin+1);
  }
  
  cout<<"Integral in acceptance for Lepton 1  V-EWK "<<integralInAcceptance/leptonEta_EWK->Integral()<<endl;
  cout<<"Integral out acceptance for Lepton 1 V-EWK "<<integralOutAcceptance/leptonEta_EWK->Integral()<<endl;


  integralInAcceptance  = 0;
  integralOutAcceptance = 0;

  for(int iBin = 0; iBin < leptonEta_2_QCD->GetNbinsX()+1; iBin++){
    if(fabs(leptonEta_2_QCD->GetBinLowEdge(iBin+1)) > 3.2)
      integralOutAcceptance += leptonEta_2_QCD->GetBinContent(iBin+1);
    else
      integralInAcceptance += leptonEta_2_QCD->GetBinContent(iBin+1);
  }
  
  cout<<"Integral in acceptance for Lepton 2  V-QCD "<<integralInAcceptance/leptonEta_2_QCD->Integral()<<endl;
  cout<<"Integral out acceptance for Lepton 2 V-QCD "<<integralOutAcceptance/leptonEta_2_QCD->Integral()<<endl;

  integralInAcceptance  = 0;
  integralOutAcceptance = 0;

  for(int iBin = 0; iBin < leptonEta_2_EWK->GetNbinsX()+1; iBin++){
    if(fabs(leptonEta_2_EWK->GetBinLowEdge(iBin+1)) > 3.2)
      integralOutAcceptance += leptonEta_2_EWK->GetBinContent(iBin+1);
    else
      integralInAcceptance += leptonEta_2_EWK->GetBinContent(iBin+1);
  }
  
  cout<<"Integral in acceptance for Lepton 2  V-EWK "<<integralInAcceptance/leptonEta_2_EWK->Integral()<<endl;
  cout<<"Integral out acceptance for Lepton 2 V-EWK "<<integralOutAcceptance/leptonEta_2_EWK->Integral()<<endl;


}  
