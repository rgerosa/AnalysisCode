#include "../CMS_lumi.h"

vector<float> bosonPt  {150.,170.,200.,230.,260.,290.0,320.0,350.0,390.0,430.0,470.0,510.0,550.0,590.0,640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0};

enum class Sample {znn, zll, wjet, gam};
float lumi_ = 1;

void makeKFactorMonojet(string inputDIR_LO, string inputDIR_NLO, string outputDIR, Sample sample){

  // to clone the EWK corrections:
  TFile* inputEWK = TFile::Open("../../data/kFactors/uncertainties_EWK_24bins.root","READ");
  TH1F* ewkCorr =  NULL;
  if(sample == Sample::znn or sample == Sample::zll){
    ewkCorr = (TH1F*) inputEWK->Get("EWKcorr/Z");
    ewkCorr->Divide((TH1F*) inputEWK->Get("ZJets_012j_NLO/nominal"));
  }
  if(sample == Sample::wjet){
    ewkCorr = (TH1F*) inputEWK->Get("EWKcorr/W");
    ewkCorr->Divide((TH1F*) inputEWK->Get("WJets_012j_NLO/nominal"));
  }
  if(sample == Sample::gam){
    ewkCorr = (TH1F*) inputEWK->Get("EWKcorr/photon");
    ewkCorr->Divide((TH1F*) inputEWK->Get("GJets_1j_NLO/nominal_G"));
  }

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
  string linea;
  if(infileLO.is_open()){
    while(!infileLO.eof()){
      getline(infileLO,linea);
      if(not TString(linea).Contains("root") or linea == "") continue;
      file_LO.push_back(TFile::Open((inputDIR_LO+"/"+linea).c_str(),"READ"));
      if((TTree*) file_LO.back()->Get("tree/tree"))
	tree_LO.push_back((TTree*) file_LO.back()->Get("tree/tree"));
      else
	tree_LO.push_back((TTree*) file_LO.back()->Get("gentree/tree"));

    }
  }
  infileLO.close();

  system(("ls "+inputDIR_NLO+" | grep root > list.txt").c_str());
  ifstream infileNLO ("list.txt");
  if(infileNLO.is_open()){
    while(!infileNLO.eof()){
      getline(infileNLO,linea);
      if(not TString(linea).Contains("root") or linea == "") continue;
      file_NLO.push_back(TFile::Open((inputDIR_NLO+"/"+linea).c_str(),"READ"));
      if((TTree*) file_NLO.back()->FindObjectAny("tree/tree"))
	tree_NLO.push_back((TTree*) file_NLO.back()->Get("tree"));
      else
	tree_NLO.push_back((TTree*) file_NLO.back()->Get("gentree/tree"));
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

      // filter away bad events in which we don't have a good weak boson candidate
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

      // filter away bad events in which we don't have a good weak boson candidate
      if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      else if(sample == Sample::gam and fabs(*wzid) != 22) continue;

      sumwgt += *wgt;
    }
    cout<<"Tree NLO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
    sumwgt_nlo.push_back(sumwgt);
    ifile++;
  }

  // histograms
  TH1F* bosonPt_LO_monojet  = new TH1F("bosonPt_LO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_monojet = new TH1F("bosonPt_NLO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  bosonPt_LO_monojet->Sumw2();
  bosonPt_NLO_monojet->Sumw2();

  // Loop on LO trees
  ifile = 0;
  for(auto tree: tree_LO){

    TTreeReader reader (tree);

    TTreeReaderValue<float> xsec       (reader,"xsec");
    TTreeReaderValue<float> wgt        (reader,"wgt");
    TTreeReaderValue<int>   l1id       (reader,"l1id");
    TTreeReaderValue<int>   l2id       (reader,"l2id");
    // post-showered i.e. status 1                                                                                                                                                                   
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<int>   wzid    (reader,"wzid");

    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");  
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");  
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

      if(*wzpt < bosonPt.front()) continue;

      // filter out taus
      if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus
      if(sample == Sample::zll  and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus

      // for photons only look at the central region
      if(sample == Sample::gam and fabs(*wzeta) > 1.449) continue;
      if(sample == Sample::zll and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::zll and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::zll and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::zll and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::wjet and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::wjet and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;

      if(sample == Sample::wjet and fabs(*l1id) == 11 and *l1pt < 40) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 13 and *l1pt < 20) continue;

      TLorentzVector lepton1, lepton2;
      lepton1.SetPtEtaPhiM(*l1pt,*l1eta,*l1phi,0.);
      lepton2.SetPtEtaPhiM(*l2pt,*l2eta,*l2phi,0.);

      if(sample == Sample::zll and max(*l1pt,*l2pt) < 20) continue;
      if(sample == Sample::zll and min(*l1pt,*l2pt) < 10) continue;
      if(sample == Sample::zll and (lepton1+lepton2).M() < 60) continue;
      if(sample == Sample::zll and (lepton1+lepton2).M() > 120) continue;

      // jets not in overlap with leptions
      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
	if(sample == Sample::wjet or sample == Sample::znn or sample == Sample::zll){
	  float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	  float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	  if(dphi1 > TMath::Pi())
	    dphi1 = 2*TMath::Pi()-dphi1;
	  if(dphi2 > TMath::Pi())
	    dphi2 = 2*TMath::Pi()-dphi2;
	  if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	     sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4
	     ){ // check cleaning
	    jets.push_back(jet4V);
	  }
	}
 	else if(sample == Sample::gam){
	  float dphi1 = fabs(jetphi->at(ijet)-*wzphi);
	  if(dphi1 > TMath::Pi())
	    dphi1 = 2*TMath::Pi()-dphi1;
	  if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*wzeta)*fabs(jeteta->at(ijet)-*wzeta)) > 0.4)
	    jets.push_back(jet4V);
	}
      }
      
      if(jets.size() < 1) continue;
      if(jets.at(0).Pt()  < 100) continue;
      if(fabs(jets.at(0).Eta()) > 2.5) continue;
      // jet met dphi cut
      // calculate min-dphi at gen level where met is boson 4V
      float mindphi   = 100;
      for(size_t ijet = 0; ijet < jets.size(); ijet++){
	if(ijet > 3) break; // limiting min dphi to first 4 leading jets
	float dphi = fabs(*wzphi-jets.at(ijet).Phi());
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(dphi < mindphi)
	  mindphi = dphi;	
      }

      if(mindphi < 0.5) continue;

      bosonPt_LO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));
    }
    ifile++;
  }

  ifile = 0;
  for(auto tree: tree_NLO){

    TTreeReader reader (tree);

    TTreeReaderValue<float> xsec       (reader,"xsec");
    TTreeReaderValue<float> wgt        (reader,"wgt");
    TTreeReaderValue<int>   l1id       (reader,"l1id");
    TTreeReaderValue<int>   l2id       (reader,"l2id");
    // post-showered i.e. status 1                                                                                                                                                                   
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<int>   wzid    (reader,"wzid");

    TTreeReaderValue<float> wzmass_lhe (reader,"mvmass");
    TTreeReaderValue<float> wzeta_lhe  (reader,"mveta");
    TTreeReaderValue<float> wzphi_lhe  (reader,"mvphi");
    TTreeReaderValue<float> wzpt_lhe   (reader,"mvpt");
    TTreeReaderValue<int>   wzid_lhe   (reader,"mvid");

    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");  
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");  
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

      if(*wzpt < bosonPt.front()) continue;

      if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus
      if(sample == Sample::zll  and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus

      // for photons only look at the central region
      if(sample == Sample::gam and fabs(*wzeta) > 1.449) continue;
      if(sample == Sample::zll and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::zll and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::zll and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::zll and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::wjet and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
      if(sample == Sample::wjet and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;

      if(sample == Sample::wjet and fabs(*l1id) == 11 and *l1pt < 40) continue;
      if(sample == Sample::wjet and fabs(*l1id) == 13 and *l1pt < 20) continue;

      TLorentzVector lepton1, lepton2;
      lepton1.SetPtEtaPhiM(*l1pt,*l1eta,*l1phi,0.);
      lepton2.SetPtEtaPhiM(*l2pt,*l2eta,*l2phi,0.);

      if(sample == Sample::zll and max(*l1pt,*l2pt) < 20) continue;
      if(sample == Sample::zll and min(*l1pt,*l2pt) < 10) continue;
      if(sample == Sample::zll and (lepton1+lepton2).M() < 60) continue;
      if(sample == Sample::zll and (lepton1+lepton2).M() > 120) continue;

      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
	if(sample == Sample::wjet or sample == Sample::znn or sample == Sample::zll){
	  float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	  float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	  if(dphi1 > TMath::Pi())
	    dphi1 = 2*TMath::Pi()-dphi1;
	  if(dphi2 > TMath::Pi())
	    dphi2 = 2*TMath::Pi()-dphi2;
	  if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	     sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4){ // check cleaning
	    jets.push_back(jet4V);
	  }
	}
	else if(sample == Sample::gam){
	  float dphi1 = fabs(jetphi->at(ijet)-*wzphi);
	  if(dphi1 > TMath::Pi())
	    dphi1 = 2*TMath::Pi()-dphi1;
	  if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*wzeta)*fabs(jeteta->at(ijet)-*wzeta)) > 0.4)
	    jets.push_back(jet4V);
	}
      }

      if(jets.size() < 1) continue;
      if(jets.at(0).Pt()  < 100) continue;
      if(fabs(jets.at(0).Eta()) > 2.5) continue;
      
      // calculate min-dphi at gen level where met is boson 4V
      float mindphi   = 100;
      for(size_t ijet = 0; ijet < jets.size(); ijet++){
	if(ijet > 3) break; // limiting min dphi to first 4 leading jets
	float dphi = fabs(*wzphi-jets.at(ijet).Phi());
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(dphi < mindphi)
	  mindphi = dphi;	
      }

      if(sample == Sample::wjet and *wzpt_lhe >= 400 and fabs(*wgt) > 100) continue;
      if(sample == Sample::wjet and *wzpt_lhe >= 600 and fabs(*wgt) > 10)  continue;
      if(sample == Sample::znn  and *wzpt_lhe >= 400 and fabs(*wgt) > 10)  continue;
      if(sample == Sample::zll  and *wzpt_lhe >= 400 and fabs(*wgt) > 10)  continue;
      if(sample == Sample::znn  and *wzpt_lhe >= 650 and fabs(*wgt) > 1)   continue;
      if(sample == Sample::zll  and *wzpt_lhe >= 650 and fabs(*wgt) > 1)   continue;
      if(mindphi < 0.5) continue;
      
      // apply monojet-like selections
      bosonPt_NLO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
    }
    ifile++;
  }
  
  TH1F* bosonPt_NLOEWK = (TH1F*) bosonPt_NLO_monojet->Clone("bosonPt_NLOEWK");
  bosonPt_NLOEWK->Multiply(ewkCorr);


  string postfix;
  if(sample == Sample::znn)       postfix = "znn";
  else if(sample == Sample::zll)  postfix = "zll";
  else if(sample == Sample::wjet) postfix = "wjet";
  else if(sample == Sample::gam)  postfix = "gam";

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetLogy();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.27);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->SetBottomMargin(0.35);

  TH1F* frame =  canvas->DrawFrame(bosonPt_NLO_monojet->GetXaxis()->GetXmin(),bosonPt_NLO_monojet->GetMinimum()*0.1,1090,bosonPt_NLO_monojet->GetMaximum()*100, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV]");
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame ->Draw();

  bosonPt_LO_monojet->SetLineColor(kBlack);
  bosonPt_LO_monojet->SetMarkerColor(kBlack);
  bosonPt_LO_monojet->SetLineWidth(2);  
  bosonPt_LO_monojet->SetMarkerStyle(20);
  bosonPt_LO_monojet->SetMarkerSize(1);
  bosonPt_NLO_monojet->SetLineColor(kRed);
  bosonPt_NLO_monojet->SetMarkerColor(kRed);
  bosonPt_NLO_monojet->SetLineWidth(2);
  bosonPt_NLO_monojet->SetMarkerStyle(20);
  bosonPt_NLO_monojet->SetMarkerSize(1);
  bosonPt_NLOEWK->SetLineColor(kBlue);
  bosonPt_NLOEWK->SetMarkerColor(kBlue);
  bosonPt_NLOEWK->SetMarkerStyle(20);
  bosonPt_NLOEWK->SetMarkerSize(1);
  bosonPt_NLOEWK->SetLineWidth(2);

  bosonPt_LO_monojet->Draw("hist same");
  bosonPt_NLO_monojet->Draw("hist same");
  bosonPt_NLOEWK->Draw("hist same");

  CMS_lumi(canvas,"");

  TLegend* leg = new TLegend(0.48, 0.66, 0.90, 0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  if(sample == Sample::znn)
    leg->AddEntry((TObject*)(0),"Z+jets","");
  else if(sample == Sample::wjet)
    leg->AddEntry((TObject*)(0),"W+jets","");
  else if(sample == Sample::gam)
    leg->AddEntry((TObject*)(0),"#gamma+jets","");

  leg->AddEntry(bosonPt_NLOEWK,"NLO QCD + EWK","FL");
  leg->AddEntry(bosonPt_NLO_monojet,"NLO-QCD","FL");
  leg->AddEntry(bosonPt_LO_monojet,"LO","FL");
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();


  canvas->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1F* frame2 =  pad2->DrawFrame(bosonPt_NLO_monojet->GetXaxis()->GetXmin(),0.5,bosonPt_NLO_monojet->GetXaxis()->GetXmax(),2.0, "");

  frame2->GetXaxis()->SetLabelSize(0.12);
  frame2->GetXaxis()->SetTitleSize(0.15);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.12);
  frame2->GetYaxis()->SetTitleSize(0.15);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetXaxis()->SetTitle("boson p_{T} [GeV]");
  frame2->GetYaxis()->SetTitle("#sigma_{NLO}/#sigma_{LO}");
  frame2->GetYaxis()->SetNdivisions(505);
  frame2->Draw();

  TH1* ratio = (TH1*) bosonPt_NLO_monojet->Clone("ratio");
  ratio->Divide(bosonPt_LO_monojet);

  TH1* ratio2 = (TH1*) bosonPt_NLOEWK->Clone("ratio2");
  ratio2->Divide(bosonPt_LO_monojet);

  TF1* line = new TF1("line","1",bosonPt_NLO_monojet->GetXaxis()->GetXmin(),bosonPt_NLO_monojet->GetXaxis()->GetXmax());
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  ratio->SetLineWidth(1);
  ratio->SetMarkerSize(1);
  ratio->SetMarkerStyle(20);

  ratio->Draw("Psame");
  ratio2->Draw("Psame");
  line->Draw("same");
  pad2->RedrawAxis("sameaxis");

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

  if(sample == Sample::znn){
    output->mkdir("ZJets_012j_NLO");
    output->cd("ZJets_012j_NLO");
    bosonPt_NLO_monojet->Write("nominal");
    output->cd();
    output->mkdir("ZJets_LO");
    output->cd("ZJets_LO");
    bosonPt_LO_monojet->Write("inv_pt");
    output->cd();
    output->mkdir("EWKcorr");
    output->cd("EWKcorr");
    bosonPt_NLOEWK->Write("Z");
    output->cd();
  }
  else if(sample == Sample::zll){
    output->mkdir("DYJets_012j_NLO");
    output->cd("DYJets_012j_NLO");
    bosonPt_NLO_monojet->Write("nominal");
    output->cd();
    output->mkdir("DYJets_LO");
    output->cd("DYJets_LO");
    bosonPt_LO_monojet->Write("inv_pt");
    output->cd();
    output->mkdir("EWKcorr");
    output->cd("EWKcorr");
    bosonPt_NLOEWK->Write("DY");
    output->cd();
  }
  else if(sample == Sample::wjet){
    output->mkdir("WJets_012j_NLO");
    output->cd("WJets_012j_NLO");
    bosonPt_NLO_monojet->Write("nominal");
    output->cd();
    output->mkdir("WJets_LO");
    output->cd("WJets_LO");
    bosonPt_LO_monojet->Write("inv_pt");
    output->cd();
    output->mkdir("EWKcorr");
    output->cd("EWKcorr");
    bosonPt_NLOEWK->Write("W");
    output->cd();
  }
  else if(sample == Sample::gam){
    output->mkdir("GJets_1j_NLO");
    output->cd("GJets_1j_NLO");
    bosonPt_NLO_monojet->Write("nominal_G");
    output->cd();
    output->mkdir("GJets_LO");
    output->cd("GJets_LO");
    bosonPt_LO_monojet->Write("inv_pt_G");
    output->cd();
    output->mkdir("EWKcorr");
    output->cd("EWKcorr");
    bosonPt_NLOEWK->Write("photon");
    output->cd();
  }
  
  output->Close();

}
