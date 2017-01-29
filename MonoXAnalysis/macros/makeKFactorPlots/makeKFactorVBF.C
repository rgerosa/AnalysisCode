#include "../CMS_lumi.h"

vector<float> bosonPt          {150.,200.,250.,350,500.,700.,1000};
vector<float> bosonPt_vbf      {150.,200.,250.,350,500.,700.,1000};

static float mjj            = 800;
static float mjjrelaxed     = 350;
static float detajj         = 2.75;
static float detajjrelaxed  = 1.75;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;
static float dphijj         = 2.5;

enum class Sample {znn, zll, wjet, gam};
float lumi_ = 1;

void makeKFactorVBF(string inputDIR_LO, string inputDIR_NLO, string outputDIR, Sample sample){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  float scale_lo = 1;
  float scale_nlo = 1;
  if(sample == Sample::znn)
    scale_nlo = 3;

  system(("mkdir -p "+outputDIR).c_str());
  vector<TTree*> tree_LO;
  vector<TTree*> tree_NLO;
  vector<TFile*> file_LO;
  vector<TFile*> file_NLO;
  
  system(("ls "+inputDIR_LO+" | grep root > list.txt").c_str());
  ifstream infileLO ("list.txt");
  string line;
  if(infileLO.is_open()){
    while(!infileLO.eof()){
      getline(infileLO,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_LO.push_back(TFile::Open((inputDIR_LO+"/"+line).c_str(),"READ"));
      tree_LO.push_back((TTree*) file_LO.back()->FindObjectAny("tree"));
    }
  }
  infileLO.close();

  system(("ls "+inputDIR_NLO+" | grep root > list.txt").c_str());
  ifstream infileNLO ("list.txt");
  if(infileNLO.is_open()){
    while(!infileNLO.eof()){
      getline(infileNLO,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_NLO.push_back(TFile::Open((inputDIR_NLO+"/"+line).c_str(),"READ"));
      tree_NLO.push_back((TTree*) file_NLO.back()->FindObjectAny("tree"));
    }
  }
  infileNLO.close();

  system("rm list.txt");

  // calculate sumwgt
  vector<double> sumwgt_lo;
  int ifile = 0;
  for(auto tree : tree_LO){
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");  
    TTreeReaderValue<float> wgt   (reader,"wgt");  
    TTreeReaderValue<int>   wzid  (reader,"wzid");  

    cout<<"Calculate sumwgt for LO file "<<file_LO.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
    while(reader.Next()){

      // filter away bad events
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      else if(sample == Sample::gam and fabs(*wzid) != 22) continue;

      sumwgt += *wgt;
    }
    cout<<"Tree LO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
    sumwgt_lo.push_back(sumwgt);
    ifile++;
  }


  vector<double> sumwgt_nlo;
  ifile = 0;
  for(auto tree : tree_NLO){
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");  
    TTreeReaderValue<int>   wzid  (reader,"wzid");  
    TTreeReaderValue<float> wgt   (reader,"wgt");  

    cout<<"Calculate sumwgt for NLO file "<<file_NLO.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
    while(reader.Next()){
      sumwgt += *wgt;
    }
    cout<<"Tree NLO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
    sumwgt_nlo.push_back(sumwgt);
    ifile++;
  }

  // histograms
  TH1F* bosonPt_LO_monojet  = new TH1F("bosonPt_LO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_monojet = new TH1F("bosonPt_NLO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_LO_twojet   = new TH1F("bosonPt_LO_twojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_twojet  = new TH1F("bosonPt_NLO_twojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_LO_vbf_relaxed      = new TH1F("bosonPt_LO_vbf_relaxed","", bosonPt_vbf.size()-1,&bosonPt_vbf[0]);
  TH1F* bosonPt_NLO_vbf_relaxed     = new TH1F("bosonPt_NLO_vbf_relaxed","", bosonPt_vbf.size()-1,&bosonPt_vbf[0]);
  TH1F* bosonPt_LO_vbf      = new TH1F("bosonPt_LO_vbf","", bosonPt_vbf.size()-1,&bosonPt_vbf[0]);
  TH1F* bosonPt_NLO_vbf     = new TH1F("bosonPt_NLO_vbf","", bosonPt_vbf.size()-1,&bosonPt_vbf[0]);

  bosonPt_LO_monojet->Sumw2();
  bosonPt_NLO_monojet->Sumw2();
  bosonPt_LO_twojet->Sumw2();
  bosonPt_NLO_twojet->Sumw2();
  bosonPt_LO_vbf_relaxed->Sumw2();
  bosonPt_NLO_vbf_relaxed->Sumw2();
  bosonPt_LO_vbf->Sumw2();
  bosonPt_NLO_vbf->Sumw2();

  // Loop on LO trees
  ifile = 0;
  for(auto tree: tree_LO){

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

  
    cout<<"Loop on LO file "<<file_LO.at(ifile)->GetName()<<endl;
    while(reader.Next()){

      // filter away bad events
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      else if(sample == Sample::gam and fabs(*wzid) != 22) continue;
      
      // calculate min-dphi at gen level where met is boson 4V
      float mindphi   = 100;
      for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
	if(ijet > 3) break; // limiting min dphi to first 4 leading jets
	float dphi = fabs(*wzphi-jetphi->at(ijet));
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(dphi < mindphi)
	  mindphi = dphi;	
      }
      
      if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus

      // apply monojet-like selections
      if(*njets >= 1 and jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and *wzpt >= 150 and mindphi > 0.5)
	bosonPt_LO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));

      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > 0.5){      

	// fill two jet phase space
	bosonPt_LO_twojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));

	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));
	if((jet1+jet2).M() < mjjrelaxed) continue;
	if(jet1.Eta()*jet2.Eta()  > 0) continue; 
	if(fabs(jeteta->at(0)-jeteta->at(1)) < detajjrelaxed) continue;
	float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	if(deltaPhi > dphijj) continue;
	bosonPt_LO_vbf_relaxed->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      

	if((jet1+jet2).M() < mjj) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
	bosonPt_LO_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      

      }
    }
    ifile++;
  }

  ifile = 0;
  for(auto tree: tree_NLO){

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

  
    cout<<"Loop on NLO file "<<file_NLO.at(ifile)->GetName()<<endl;

    if(sample == Sample::zll)
      scale_nlo = 1.07;
    
    while(reader.Next()){
      
      // filter away bad events with no matching
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      else if(sample == Sample::gam and fabs(*wzid) != 22) continue;

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
      
      // apply monojet-like selections
      if(*njets >= 1 and jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and *wzpt >= 150 and mindphi > 0.5){
	  bosonPt_NLO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
      }

      if(*njetsinc >= 2 and jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)) < 4.7 and  fabs(jeteta->at(1)) < 4.7 and mindphi > 0.5){      
	bosonPt_NLO_twojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));
	if((jet1+jet2).M() < mjjrelaxed) continue;
	if(jet1.Eta()*jet2.Eta()  > 0) continue; 
	if(fabs(jeteta->at(0)-jeteta->at(1)) < detajjrelaxed) continue;
	float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	if(deltaPhi > dphijj) continue;
	bosonPt_NLO_vbf_relaxed->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));      

	if((jet1+jet2).M() < mjj) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
	bosonPt_NLO_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));	
      }
    }
    ifile++;
  }
  
  string postfix;
  if(sample == Sample::znn)       postfix = "znn";
  else if(sample == Sample::zll)  postfix = "zll";
  else if(sample == Sample::wjet) postfix = "wjet";
  else if(sample == Sample::gam)  postfix = "gam";


  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  canvas->cd();


  TH1F* kfactor_monojet = (TH1F*) bosonPt_NLO_monojet->Clone("kfactor_monojet");
  kfactor_monojet->Divide(bosonPt_LO_monojet);
  TH1F* kfactor_twojet = (TH1F*) bosonPt_NLO_twojet->Clone("kfactor_twojet");
  kfactor_twojet->Divide(bosonPt_LO_twojet);
  TH1F* kfactor_vbf_relaxed = (TH1F*) bosonPt_NLO_vbf_relaxed->Clone("kfactor_vbf_relaxed");
  kfactor_vbf_relaxed->Divide(bosonPt_LO_vbf_relaxed);
  TH1F* kfactor_vbf = (TH1F*) bosonPt_NLO_vbf->Clone("kfactor_vbf");
  kfactor_vbf->Divide(bosonPt_LO_vbf);
  
  kfactor_monojet->GetYaxis()->SetTitleOffset(1.1);
  kfactor_monojet->GetXaxis()->SetTitleOffset(1.1);
  kfactor_monojet->GetYaxis()->SetRangeUser(min(kfactor_monojet->GetMinimum(),min(kfactor_twojet->GetMinimum(),min(kfactor_vbf_relaxed->GetMinimum(),kfactor_vbf->GetMinimum())))*0.7,
					    max(kfactor_monojet->GetMaximum(),max(kfactor_twojet->GetMaximum(),max(kfactor_vbf_relaxed->GetMaximum(),kfactor_vbf->GetMaximum())))*1.3);

  kfactor_monojet->SetLineColor(kBlack);
  kfactor_monojet->SetMarkerColor(kBlack);
  kfactor_monojet->SetLineWidth(2);
  kfactor_monojet->SetMarkerSize(1);

  kfactor_twojet->SetLineColor(kBlue);
  kfactor_twojet->SetMarkerColor(kBlue);
  kfactor_twojet->SetLineWidth(2);
  kfactor_twojet->SetMarkerSize(1);

  kfactor_vbf_relaxed->SetLineColor(kOrange+1);
  kfactor_vbf_relaxed->SetMarkerColor(kOrange+1);
  kfactor_vbf_relaxed->SetLineWidth(2);
  kfactor_vbf_relaxed->SetMarkerSize(1);

  kfactor_vbf->SetLineColor(kRed);
  kfactor_vbf->SetMarkerColor(kRed);
  kfactor_vbf->SetLineWidth(2);
  kfactor_vbf->SetMarkerSize(1);
  
  TH1F* kfactor_monojet_band = (TH1F*) kfactor_monojet->Clone("kfactor_monojet_band");
  kfactor_monojet_band->SetFillColor(kBlack);
  kfactor_monojet_band->SetFillStyle(3001);

  TH1F* kfactor_twojet_band = (TH1F*) kfactor_twojet->Clone("kfactor_twojet_band");
  kfactor_twojet_band->SetFillColor(kBlue);
  kfactor_twojet_band->SetFillStyle(3001);

  TH1F* kfactor_vbf_relaxed_band = (TH1F*) kfactor_vbf_relaxed->Clone("kfactor_vbf_relaxed_band");
  kfactor_vbf_relaxed_band->SetFillColor(kOrange+1);
  kfactor_vbf_relaxed_band->SetFillStyle(3001);

  TH1F* kfactor_vbf_band = (TH1F*) kfactor_vbf->Clone("kfactor_vbf_band");
  kfactor_vbf_band->SetFillColor(kRed);
  kfactor_vbf_band->SetFillStyle(3001);
  
  kfactor_monojet->GetYaxis()->SetTitle("K-factor");
  kfactor_monojet->GetXaxis()->SetLabelSize(0);
  kfactor_monojet->GetXaxis()->SetTitleSize(0);
  kfactor_monojet->Draw("hist");
  kfactor_monojet_band->Draw("E2 same");
  kfactor_twojet_band->Draw("E2 same");
  kfactor_vbf_relaxed_band->Draw("E2 same");
  kfactor_vbf_band->Draw("E2 same");

  kfactor_monojet->Draw("hist same");
  kfactor_twojet->Draw("hist same");
  kfactor_vbf_relaxed->Draw("hist same");
  kfactor_vbf->Draw("hist same");

  TLegend leg (0.65,0.65,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(kfactor_monojet_band,"K-factor monojet","FL");
  leg.AddEntry(kfactor_twojet_band, "K-factor 2-jet","FL");
  leg.AddEntry(kfactor_vbf_relaxed_band,    "K-factor VBF relaxed","FL");
  leg.AddEntry(kfactor_vbf_band,    "K-factor VBF","FL");
  leg.Draw("same");

  TLatex* latex = new TLatex ();
  latex->SetNDC();
  latex->SetTextSize(0.6*gPad->GetTopMargin());
  latex->SetTextFont(62);
  latex->SetTextAlign(11);
  if(sample == Sample::znn)
    latex->DrawLatex(0.15,0.95,"CMS Simulation Z+jets");
  else if(sample == Sample::zll)
    latex->DrawLatex(0.15,0.95,"CMS Simulation DY+jets");
  else if(sample == Sample::wjet)
    latex->DrawLatex(0.15,0.95,"CMS Simulation W+jets");
  else if(sample == Sample::gam)
    latex->DrawLatex(0.15,0.95,"CMS Simulation #gamma+jets");

  pad2->Draw();
  pad2->cd();

  //ratio plot
  TH1* ratio_1 = (TH1*) kfactor_vbf->Clone("ratio_1");
  ratio_1->Divide(kfactor_monojet);
  TH1* ratio_2 = (TH1*) kfactor_vbf_relaxed->Clone("ratio_2");
  ratio_2->Divide(kfactor_monojet);
  TH1* ratio_3 = (TH1*) kfactor_twojet->Clone("ratio_3");
  ratio_3->Divide(kfactor_monojet);

  ratio_3->SetLineColor(kBlue);
  ratio_3->SetMarkerColor(kBlue);
  ratio_3->GetYaxis()->SetTitle("Ratio");
  if(sample == Sample::gam)
    ratio_3->GetYaxis()->SetRangeUser(0.85,1.5);
  else
    ratio_3->GetYaxis()->SetRangeUser(0.85,1.25);
  
  ratio_3->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  ratio_3->GetYaxis()->CenterTitle();
  ratio_3->GetYaxis()->SetTitleOffset(1.5);
  ratio_3->GetYaxis()->SetLabelSize(0.04);
  ratio_3->GetYaxis()->SetTitleSize(0.04);
  ratio_3->GetXaxis()->SetLabelSize(0.04);
  ratio_3->GetXaxis()->SetTitleSize(0.05);  
  ratio_3->GetYaxis()->SetNdivisions(5);
  ratio_3->Draw("hist");

  ratio_2->SetLineColor(kOrange+1);
  ratio_2->SetMarkerColor(kOrange+1);
  ratio_2->Draw("hist same");

  ratio_1->SetLineColor(kRed);
  ratio_1->SetMarkerColor(kRed);
  ratio_1->Draw("hist same");

  TH1F* ratio_3_band = (TH1F*) ratio_3->Clone("ratio_3_band");
  ratio_3_band->SetFillColor(kBlue);
  ratio_3_band->SetFillStyle(3001);
  ratio_3_band->Draw("E2 same");

  TH1F* ratio_2_band = (TH1F*) ratio_2->Clone("ratio_2_band");
  ratio_2_band->SetFillColor(kOrange+1);
  ratio_2_band->SetFillStyle(3001);
  ratio_2_band->Draw("E2 same");

  TH1F* ratio_1_band = (TH1F*) ratio_1->Clone("ratio_1_band");
  ratio_1_band->SetFillColor(kRed);
  ratio_1_band->SetFillStyle(3001);
  ratio_1_band->Draw("E2 same");

  ratio_3->Draw("hist same");
  ratio_2->Draw("hist same");
  ratio_1->Draw("hist same");

  if(sample == Sample::znn){
    canvas->SaveAs((outputDIR+"/kfactor_znn.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_znn.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::zll){
    canvas->SaveAs((outputDIR+"/kfactor_zll.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_zll.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDIR+"/kfactor_wjet.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_wjet.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDIR+"/kfactor_gam.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_gam.pdf").c_str(),"pdf");
  }

  TFile* output = new TFile((outputDIR+"/kfactor_"+postfix+".root").c_str(),"RECREATE");
  output->cd();
  bosonPt_NLO_monojet->Write();
  bosonPt_NLO_twojet->Write();
  bosonPt_NLO_vbf_relaxed->Write();
  bosonPt_NLO_vbf->Write();
  bosonPt_LO_monojet->Write();
  bosonPt_LO_twojet->Write();
  bosonPt_LO_vbf_relaxed->Write();
  bosonPt_LO_vbf->Write();
  output->Close();

}
