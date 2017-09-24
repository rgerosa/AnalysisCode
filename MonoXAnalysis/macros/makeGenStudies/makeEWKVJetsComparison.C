#include "../CMS_lumi.h"

static float bosonPt        = 200;
static float mjj            = 450;
static float detajj         = 2.5;
static float dphijj         = 1.5;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;

enum class Sample {znn, zll, wjet, gam};
float lumi_ = 12.9;

vector<float> bin_bosonPt      {200.,250.,300.,350.,400.,450.,500.,550.,600.,700.,850.,1000.,1300.};
vector<float> bin_bosonPt_vbf  {250.,300.,350.,450.,550.,750.,1000.,1300.};
vector<float> bin_jetpt1       {70.,85.,100.,125,150,175,200,250.,300.,350,400,450.,500.,600.,700.,850.,1000.};
vector<float> bin_jetpt1_vbf   {70.,85.,100.,125,150,175,200,250.,300.,350,400,450.,500.,600.,700.,850.,1000.};
vector<float> bin_jetpt2       {50,60,70.,80.,90,100.,110,125,150,175,200,250.,300.,};
vector<float> bin_jetpt2_vbf   {50,60,70.,80.,90,100.,110,125,150,175,200,250.,300.,};
vector<float> bin_mjj          {200.,300.,400.,500.,600.,700.,800.,900.,1000.,1250.,1500.,1750.,2000.,2250.,2500.,3000.};
vector<float> bin_mjj_vbf      {450.,550.,650.,750.,850.,950.,1050.,1250.,1500.,1750.,2000.,2250.,2500.,3000.};
vector<float> bin_detajj       {0.,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,9};
vector<float> bin_detajj_vbf   {2.5,3,3.5,4,4.5,5,5.5,6,9};
vector<float> bin_jetmetdphi   {0.5,0.8,1.2,1.5,1.8,2.1,2.4,2.7,3.14};
vector<float> bin_jetmetdphi_vbf {1.0,1.4,1.8,2.2,2.6,3,3.14};

/////////////////////////
void plotResults(TH1* histo_ewk, TH1* histo_qcd,const string & postfix,const string & xAxisLabel, const string & outputDIR){

  
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,650);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  TPad *pad1 = new TPad(("pad1_"+postfix).c_str(),"",0,0.20,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad(("pad2_"+postfix).c_str(),"",0,0.05,1,0.25);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->SetGridy();
  
  canvas->cd();
  pad1->Draw();
  canvas->cd();
  pad2->Draw();
  canvas->cd();
  pad1->cd();
      
  histo_ewk->SetTitle("");
  histo_qcd->GetYaxis()->SetTitle("");
  histo_ewk->GetYaxis()->SetTitle("Entries / Bin width");
  
  histo_ewk->Scale(1.,"width");
  histo_qcd->Scale(1.,"width");


  histo_ewk->SetLineColor(kRed);
  histo_ewk->SetLineWidth(2);
  histo_ewk->SetMarkerStyle(20);
  histo_ewk->SetMarkerSize(1);
  histo_ewk->SetMarkerColor(kRed);

  TH1* histo_qcd_band = (TH1*) histo_qcd->Clone("histo_qcd_band");
  histo_qcd_band->SetFillColor(kGray);
  histo_qcd_band->SetFillStyle(3001);
  
  histo_qcd->SetLineColor(kBlack);
  histo_qcd->SetLineWidth(2);
  histo_qcd->SetMarkerColor(kBlack);
  histo_qcd->SetMarkerSize(0.8);
  histo_qcd->SetMarkerStyle(20);

  histo_ewk->GetYaxis()->SetRangeUser(min(histo_ewk->GetMinimum(),histo_qcd->GetMinimum())*0.1,max(histo_ewk->GetMaximum(),histo_qcd->GetMaximum())*500);
  histo_ewk->Draw("PE");
  CMS_lumi(pad1,Form("%.1f",lumi_),true);

  histo_qcd_band->Draw("E2 same");
  histo_qcd->Draw("hist same");
  histo_ewk->Draw("PEsame");

  TLegend* leg = new TLegend(0.5,0.70,0.90,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_qcd_band, "V+jets QCD","FL");
  leg->AddEntry(histo_ewk, "V+jets EWK","PEL");
  leg->Draw("same");
  pad1->RedrawAxis("sameaxis");
  pad1->SetLogy();

  canvas->cd();
  pad2->cd();



  TH1F* ratio = (TH1F*) histo_ewk->Clone("ratio");
  ratio->Reset();
  ratio->Add(histo_ewk);
  ratio->Divide(histo_qcd);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerSize(0.8);
  ratio->SetMarkerStyle(20);
  ratio->GetYaxis()->SetRangeUser(ratio->GetMinimum()*0.75,ratio->GetMaximum()*1.25);
  ratio->GetXaxis()->SetLabelSize(0.10);
  ratio->GetXaxis()->SetLabelOffset(0.03);
  ratio->GetXaxis()->SetTitleSize(0.13);
  ratio->GetXaxis()->SetTitleOffset(1.05);
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitle(xAxisLabel.c_str());
  ratio->GetYaxis()->SetTitle("Ratio EWK/QCD");
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetNdivisions(504, false);
  ratio->Draw("PEsame");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparison_"+postfix+".pdf").c_str(),"pdf");

}


///////////////////
void makeEWKVComparison(string inputDIR_QCD, string inputDIR_EWK, Sample sample, string outputDIR, bool applyKfactor = false){

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

    TH1* anlohist = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    TH1* alohist  = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    TH1* aewkhist = (TH1*) kffile->Get("EWKcorr/photon");
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);

    if(sample == Sample::zll or sample == Sample::znn){
      khists.push_back(znlohist);
      khists.push_back(zewkhist);
    }
    else if(sample == Sample::wjet){
      khists.push_back(wnlohist);
      khists.push_back(wewkhist);
    }
    else if(sample == Sample::gam){
      khists.push_back(aewkhist);
      khists.push_back(anlohist);
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
  double eventsNoBoson;
  for(auto tree : tree_QCD){
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
	continue;
      }
      else if(sample == Sample::wjet and fabs(*wzid) != 24){
	eventsNoBoson++;
	continue;
      }
      else if(sample == Sample::gam and fabs(*wzid) != 22){
	eventsNoBoson++;
	continue;
      }

      sumwgt += *wgt;
    }
    cout<<"Tree QCD with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
    sumwgt_qcd.push_back(sumwgt);
    ifile++;
  }

  vector<double> sumwgt_ewk;
  ifile = 0;
  eventsNoBoson = 0;
  for(auto tree : tree_EWK){
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
	continue;
      }
      else if(sample == Sample::wjet and fabs(*wzid) != 24){
	eventsNoBoson++;
	continue;
      }
      else if(sample == Sample::gam and fabs(*wzid) != 22){
	eventsNoBoson++;
	continue;
      }

      sumwgt += *wgt;
    }
    cout<<"Tree EWK with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
    sumwgt_ewk.push_back(sumwgt);
    ifile++;
  }

  // distributions                                                                                                                                                                                    
  TH1F* bosonPt_QCD = new TH1F("bosonPt_QCD","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* bosonPt_EWK = new TH1F("bosonPt_EWK","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  bosonPt_QCD->Sumw2();
  bosonPt_EWK->Sumw2();

  TH1F* bosonPt_QCD_vbf = new TH1F("bosonPt_QCD_vbf","",bin_bosonPt_vbf.size()-1,&bin_bosonPt_vbf[0]);
  TH1F* bosonPt_EWK_vbf = new TH1F("bosonPt_EWK_vbf","",bin_bosonPt_vbf.size()-1,&bin_bosonPt_vbf[0]);
  bosonPt_QCD_vbf->Sumw2();
  bosonPt_EWK_vbf->Sumw2();

  TH1F* mjj_QCD = new TH1F("mjj_QCD","",bin_mjj.size()-1,&bin_mjj[0]);
  TH1F* mjj_EWK = new TH1F("mjj_EWK","",bin_mjj.size()-1,&bin_mjj[0]);
  mjj_QCD->Sumw2();
  mjj_EWK->Sumw2();

  TH1F* mjj_QCD_vbf = new TH1F("mjj_QCD_vbf","",bin_mjj_vbf.size()-1,&bin_mjj_vbf[0]);
  TH1F* mjj_EWK_vbf = new TH1F("mjj_EWK_vbf","",bin_mjj_vbf.size()-1,&bin_mjj_vbf[0]);
  mjj_QCD_vbf->Sumw2();
  mjj_EWK_vbf->Sumw2();

  TH1F* detajj_QCD = new TH1F("detajj_QCD","",bin_detajj.size()-1,&bin_detajj[0]);
  TH1F* detajj_EWK = new TH1F("detajj_EWK","",bin_detajj.size()-1,&bin_detajj[0]);
  detajj_QCD->Sumw2();
  detajj_EWK->Sumw2();

  TH1F* detajj_QCD_vbf = new TH1F("detajj_QCD_vbf","",bin_detajj_vbf.size()-1,&bin_detajj_vbf[0]);
  TH1F* detajj_EWK_vbf = new TH1F("detajj_EWK_vbf","",bin_detajj_vbf.size()-1,&bin_detajj_vbf[0]);
  detajj_QCD_vbf->Sumw2();
  detajj_EWK_vbf->Sumw2();

  TH1F* jetmetdphi_QCD = new TH1F("jetmetdphi_QCD","",bin_jetmetdphi.size()-1,&bin_jetmetdphi[0]);
  TH1F* jetmetdphi_EWK = new TH1F("jetmetdphi_EWK","",bin_jetmetdphi.size()-1,&bin_jetmetdphi[0]);
  jetmetdphi_QCD->Sumw2();
  jetmetdphi_EWK->Sumw2();

  TH1F* jetmetdphi_QCD_vbf = new TH1F("jetmetdphi_QCD_vbf","",bin_jetmetdphi_vbf.size()-1,&bin_jetmetdphi_vbf[0]);
  TH1F* jetmetdphi_EWK_vbf = new TH1F("jetmetdphi_EWK_vbf","",bin_jetmetdphi_vbf.size()-1,&bin_jetmetdphi_vbf[0]);
  jetmetdphi_QCD_vbf->Sumw2();
  jetmetdphi_EWK_vbf->Sumw2();

  TH1F* jetpt1_QCD = new TH1F("jetpt1_QCD","",bin_jetpt1.size()-1,&bin_jetpt1[0]);
  TH1F* jetpt1_EWK = new TH1F("jetpt1_EWK","",bin_jetpt1.size()-1,&bin_jetpt1[0]);
  jetpt1_QCD->Sumw2();
  jetpt1_EWK->Sumw2();

  TH1F* jetpt1_QCD_vbf = new TH1F("jetpt1_QCD_vbf","",bin_jetpt1_vbf.size()-1,&bin_jetpt1_vbf[0]);
  TH1F* jetpt1_EWK_vbf = new TH1F("jetpt1_EWK_vbf","",bin_jetpt1_vbf.size()-1,&bin_jetpt1_vbf[0]);
  jetpt1_QCD_vbf->Sumw2();
  jetpt1_EWK_vbf->Sumw2();

  TH1F* jetpt2_QCD = new TH1F("jetpt2_QCD","",bin_jetpt2.size()-1,&bin_jetpt2[0]);
  TH1F* jetpt2_EWK = new TH1F("jetpt2_EWK","",bin_jetpt2.size()-1,&bin_jetpt2[0]);
  jetpt2_QCD->Sumw2();
  jetpt2_EWK->Sumw2();

  TH1F* jetpt2_QCD_vbf = new TH1F("jetpt2_QCD_vbf","",bin_jetpt2_vbf.size()-1,&bin_jetpt2_vbf[0]);
  TH1F* jetpt2_EWK_vbf = new TH1F("jetpt2_EWK_vbf","",bin_jetpt2_vbf.size()-1,&bin_jetpt2_vbf[0]);
  jetpt2_QCD_vbf->Sumw2();
  jetpt2_EWK_vbf->Sumw2();

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

      // filter away bad events with no matching                                                                                                                                              
      // if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;                                                                                                 
      // else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;                                                                                                                      
      // else if(sample == Sample::gam and fabs(*wzid) != 22) continue;                                                                                                                            

      float mindphi = 100;
      for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
        if(ijet > 3) break; // limiting min dphi to first 4 leading jets                                                                                                                     
        float dphi = fabs(*wzphi-jetphi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          mindphi = dphi;
      }

      if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus                                                      


      Double_t kwgt = 1.0;
      double genpt = *wzpt;

      if(*wzpt < bosonPt) continue;

      for (size_t i = 0; i < khists.size(); i++) {
	if (khists[i]) {
	  if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	  if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	  kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
	}
      }

      // apply dijet selections
      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > 0.5){

        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));

	if(fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3) continue;
	if(fabs(jet1.DeltaPhi(jet2)) > dphijj) continue;

	//////////////////
	bosonPt_QCD->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetpt1_QCD->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetpt2_QCD->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	mjj_QCD->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	detajj_QCD->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	jetmetdphi_QCD->Fill(mindphi,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	///////////////////
	
        if((jet1+jet2).M() < mjj) continue;
        if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
        if(mindphi < 1.0) continue;

        bosonPt_QCD_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
        jetpt1_QCD_vbf->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
        jetpt2_QCD_vbf->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
	mjj_QCD_vbf->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
        detajj_QCD_vbf->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
        jetmetdphi_QCD_vbf->Fill(mindphi,lumi_*(*wgt)*(*xsec)*kwgt/sumwgt_qcd.at(ifile));
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
      
      // filter away bad events with no matching                                                                                                                                                      
      /*                                                                                                                                                                                             
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;                                                                                                            
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;                                                                                                                                 
      else if(sample == Sample::gam and fabs(*wzid) != 22) continue;                                                                                                                                  
      */

      float mindphi = 100;
      for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
        if(ijet > 3) break; // limiting min dphi to first 4 leading jets                                                                                                                     
        float dphi = fabs(*wzphi-jetphi->at(ijet));
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        if(dphi < mindphi)
          mindphi = dphi;
      }

      if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus                                                      

      if(*wzpt < bosonPt) continue;

      // apply dijet selections
      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > 0.5){

        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));

	//////////////////
	bosonPt_EWK->Fill(*wzpt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetpt1_EWK->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetpt2_EWK->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	mjj_EWK->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	detajj_EWK->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	jetmetdphi_EWK->Fill(mindphi,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	///////////////////
	
        if((jet1+jet2).M() < mjj) continue;
        if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
        if(mindphi < 1.0) continue;

        bosonPt_EWK_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
        jetpt1_EWK_vbf->Fill(jetpt->at(0),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
        jetpt2_EWK_vbf->Fill(jetpt->at(1),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
	mjj_EWK_vbf->Fill((jet1+jet2).M(),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
        detajj_EWK_vbf->Fill(fabs(jeteta->at(0)-jeteta->at(1)),lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
        jetmetdphi_EWK_vbf->Fill(mindphi,lumi_*(*wgt)*(*xsec)/sumwgt_ewk.at(ifile));
      }
    }
    ifile++;
  }

  ///////// plotting results
  plotResults(bosonPt_EWK,bosonPt_QCD,"bosonPt","Boson p_{T} [GeV]",outputDIR);
  plotResults(bosonPt_EWK_vbf,bosonPt_QCD_vbf,"bosonPt_vbf","Boson p_{T} [GeV]",outputDIR);
  ///////// plotting results
  plotResults(jetpt1_EWK,jetpt1_QCD,"jetpt1","p_{T}^{j1} [GeV]",outputDIR);
  plotResults(jetpt1_EWK_vbf,jetpt1_QCD_vbf,"jetpt1_vbf","p_{T}^{j1} [GeV]",outputDIR);
  ///////// plotting results
  plotResults(jetpt2_EWK,jetpt2_QCD,"jetpt2","p_{T}^{j2} [GeV]",outputDIR);
  plotResults(jetpt2_EWK_vbf,jetpt2_QCD_vbf,"jetpt2_vbf","p_{T}^{j2} [GeV]",outputDIR);
  ///////// plotting results
  plotResults(mjj_EWK,mjj_QCD,"mjj","m_{jj} [GeV]",outputDIR);
  plotResults(mjj_EWK_vbf,mjj_QCD_vbf,"mjj_vbf","m_{jj} [GeV]",outputDIR);
  ///////// plotting results
  plotResults(detajj_EWK,detajj_QCD,"detajj","#Delta#eta_{jj}",outputDIR);
  plotResults(detajj_EWK_vbf,detajj_QCD_vbf,"detajj_vbf","#Delta#eta_{jj}",outputDIR);
  ///////// plotting results
  plotResults(jetmetdphi_EWK,jetmetdphi_QCD,"jetmetdphi","#Delta#phi(jet,met)",outputDIR);
  plotResults(jetmetdphi_EWK_vbf,jetmetdphi_QCD_vbf,"jetmetdphi_vbf","#Delta#phi(jet,met)",outputDIR);
  /////////

}  
