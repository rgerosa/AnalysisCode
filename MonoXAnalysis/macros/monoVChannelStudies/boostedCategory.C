#include "../CMS_lumi.h"

void boostedCategory(string inputDIR, string signalHypothesis, string outputDir, bool isWboson = true){

  vector<double> bins = {200,225,250,300,350,400,450,500,550,600,650,700,750,800.,850.,900.,950.,1000.,1050,1100,1150,1200};

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);


  TChain* chain = new TChain("tree/tree");  
  if(signalHypothesis == "Scalar" and isWboson){
    chain->Add((inputDIR+"/sigfilter/*WminusH*root").c_str());
    chain->Add((inputDIR+"/sigfilter/*WplusH*root").c_str());
  }
  else if(signalHypothesis == "Scalar" and not isWboson)
    chain->Add((inputDIR+"/sigfilter/*ZH*root").c_str());
  else if(signalHypothesis == "Axial" and isWboson)
    chain->Add((inputDIR+"/sigfilter/*AxialMonoW*root").c_str());
  else if(signalHypothesis == "Axial" and not isWboson)
    chain->Add((inputDIR+"/sigfilter/*AxialMonoZ*root").c_str());
  else if(signalHypothesis == "Vector" and isWboson)
    chain->Add((inputDIR+"/sigfilter/*VectorMonoW*root").c_str());
  else if(signalHypothesis == "Vector" and not isWboson)
    chain->Add((inputDIR+"/sigfilter/*VectorMonoZ*root").c_str());
  
  cout<<"Make weights for flat V-bsoon pt "<<endl;
  TH1F* bosonPt_D = new TH1F("bosonPt_D","",120,200,1200);
  bosonPt_D->Sumw2();

  chain->Draw("wzpt_h >> bosonPt_D","wzpt_h > 150 && abs(wzeta_h) < 2.4 && njets >= 1 && ntausrawold == 0 && nbjetslowpt == 0 && combinejetCHfrac[0] > 0.1 && combinejetNHfrac[0] < 0.8&& incjetmumetdphimin4 > 0.5 && combinejetpt[0] > 100 && t1pfmet > 200 && abs(combinejeteta[0]) < 2.5","goff");
  bosonPt_D->Scale(1./bosonPt_D->Integral());
  
  TH1F* bosonPt_W = new TH1F("bosonPt_W","",120,200,1200);
  bosonPt_W->Sumw2();
  
  for(int iBin = 0 ; iBin < bosonPt_W->GetNbinsX(); iBin++)
    bosonPt_W->SetBinContent(iBin+1,1);
  bosonPt_W->Scale(1./bosonPt_W->Integral());
  bosonPt_W->Divide(bosonPt_D);

  TTreeReader myReader(chain);
  TTreeReaderValue<UChar_t> hlt       (myReader,"hltmet90");
  TTreeReaderValue<unsigned int> nbjets    (myReader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> njets     (myReader,"njets");
  TTreeReaderValue<float>       jmmdphi   (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<float>       met       (myReader,"t1pfmet");

  TTreeReaderValue<vector<float> > jetpt  (myReader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta (myReader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi (myReader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm (myReader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac (myReader,"combinejetNHfrac");
  TTreeReaderValue<vector<float> > genjetpt  (myReader,"combinejetGenpt");
  TTreeReaderValue<vector<float> > genjeteta (myReader,"combinejetGeneta");
  TTreeReaderValue<vector<float> > genjetphi (myReader,"combinejetGenphi");
  TTreeReaderValue<vector<float> > genjetm (myReader,"combinejetGenm");
  

  TTreeReaderValue<vector<float> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<float> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<float> > boostedJetGenpt    (myReader,"boostedJetGenpt");
  TTreeReaderValue<vector<float> > boostedJetGeneta   (myReader,"boostedJetGeneta");
  TTreeReaderValue<vector<float> > boostedJetGenphi   (myReader,"boostedJetGenphi");
  TTreeReaderValue<vector<float> > boostedJetGenm     (myReader,"boostedJetGenm");
  TTreeReaderValue<vector<float> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (myReader,"boostedJettau1");
  
  TTreeReaderValue<float > wzpt_h (myReader,"wzpt_h");
  TTreeReaderValue<float > wzeta_h (myReader,"wzeta_h");
  TTreeReaderValue<float > wzphi_h (myReader,"wzphi_h");
  TTreeReaderValue<float > wzmass_h (myReader,"wzmass_h");
  TTreeReaderValue<float > q1pt_h (myReader,"q1pt");
  TTreeReaderValue<float > q1eta_h (myReader,"q1eta");
  TTreeReaderValue<float > q1phi_h (myReader,"q1phi");
  TTreeReaderValue<float > q2pt_h (myReader,"q2pt");
  TTreeReaderValue<float > q2eta_h (myReader,"q2eta");
  TTreeReaderValue<float > q2phi_h (myReader,"q2phi");

  TH1F* denominator        = new TH1F("den_boosted","",bins.size()-1,&bins[0]);
  TH1F* numerator_boosted  = new TH1F("numerator_boosted","",bins.size()-1,&bins[0]);
  TH1F* numerator_resolved = new TH1F("numerator_resolved","",bins.size()-1,&bins[0]);

  denominator->Sumw2();
  numerator_boosted->Sumw2();
  numerator_resolved->Sumw2();
  cout<<"Loop on the event "<<endl;
  while(myReader.Next()){
    if(*njets  < 1) continue;
    if(*nbjets > 0) continue;
    if(*jmmdphi < 0.5) continue;
    if(*met < 200) continue;
    if(jetpt->size() <= 0) continue;
    if(nhfrac->at(0) > 0.8) continue;
    if(chfrac->at(0) < 0.1) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;

    // ensuring one hadronic W in the eta acceptance, with a GenPt > 200 GeV
    if(*wzpt_h <= 150. or fabs(*wzeta_h) > 2.4 or fabs(*q1eta_h) > 2.4 or fabs(*q2eta_h) > 2.4 or fabs(*q2pt_h) < 30 or fabs(*q1pt_h) < 30) continue;

    denominator->Fill(*wzpt_h,bosonPt_W->GetBinContent(bosonPt_W->FindBin(*wzpt_h)));

    TLorentzVector jetAK8, genJetAK8, boson;
    if(boostedJetpt->size() > 0)
      jetAK8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));
    if(boostedJetGenpt->size() > 0)
      genJetAK8.SetPtEtaPhiM(boostedJetGenpt->at(0),boostedJetGeneta->at(0),boostedJetGenphi->at(0),boostedJetGenm->at(0));
    boson.SetPtEtaPhiM(*wzpt_h,*wzeta_h,*wzphi_h,*wzmass_h);

    // look for: reco jet pt > 200, mathing with ak8 Gen Jet    
    if(boostedJetpt->size() > 0 and boostedJetpt->at(0) > 200 and jetAK8.DeltaR(genJetAK8) < 0.1 and jetAK8.DeltaR(boson) < 0.1)
      numerator_boosted->Fill(*wzpt_h,bosonPt_W->GetBinContent(bosonPt_W->FindBin(*wzpt_h)));


    //////
    TLorentzVector jetAK4_1, jetAK4_2, genJetAK4;
    bool foundW = false;
    for(size_t iJet = 0; iJet < jetpt->size(); iJet++){

      if(foundW) continue;
      jetAK4_1.SetPtEtaPhiM(jetpt->at(iJet),jeteta->at(iJet),jetphi->at(iJet),jetm->at(iJet));

      // find a gen jet matching
      bool findGoodGenJet = false;
      int genJetPos = 0.;
      for(size_t iGenJet = 0; iGenJet < genjetpt->size(); iGenJet++){
	genJetAK4.SetPtEtaPhiM(genjetpt->at(iGenJet),genjeteta->at(iGenJet),genjetphi->at(iGenJet),genjetm->at(iGenJet));
	if(jetAK4_1.DeltaR(genJetAK4) < 0.1){
	  findGoodGenJet = true;
	  genJetPos = iGenJet;
	}
      }
      
      if(not findGoodGenJet) continue;
      
      for(size_t jJet = iJet; jJet < jetpt->size(); jJet++){
	if(foundW) continue;
	jetAK4_2.SetPtEtaPhiM(jetpt->at(jJet),jeteta->at(jJet),jetphi->at(jJet),jetm->at(jJet));

	// find a gen jet matching
	bool findGoodGenJet = false;
	for(size_t iGenJet = 0; iGenJet < genjetpt->size(); iGenJet++){
	  genJetAK4.SetPtEtaPhiM(genjetpt->at(iGenJet),genjeteta->at(iGenJet),genjetphi->at(iGenJet),genjetm->at(iGenJet));
	  if(jetAK4_2.DeltaR(genJetAK4) < 0.1 and iGenJet != genJetPos)
	    findGoodGenJet = true;
	}
	
	if(not findGoodGenJet) continue;
	
	float deltaPhi_1 = fabs(jetAK4_1.Phi()-(*q1phi_h)); 
	if(deltaPhi_1 > TMath::Pi()) deltaPhi_1 = 2*TMath::Pi()-fabs(jetAK4_1.Phi()-(*q1phi_h)); 

	float deltaPhi_2 = fabs(jetAK4_1.Phi()-(*q2phi_h)); 
	if(deltaPhi_2 > TMath::Pi()) deltaPhi_2 = 2*TMath::Pi()-fabs(jetAK4_1.Phi()-(*q2phi_h)); 
	
	float deltaR_1 = sqrt(deltaPhi_1*deltaPhi_1+fabs(jetAK4_1.Eta()-*q1eta_h)*fabs(jetAK4_1.Eta()-*q1eta_h));
	float deltaR_2 = sqrt(deltaPhi_2*deltaPhi_2+fabs(jetAK4_1.Eta()-*q2eta_h)*fabs(jetAK4_1.Eta()-*q2eta_h));

	float deltaRmin = std::min(deltaR_1,deltaR_2);
	
	if(deltaRmin > 0.1) continue;

	float deltaRmin_2 = 0.;

	if(deltaRmin == deltaR_1){
	  
	  deltaPhi_2 = fabs(jetAK4_2.Phi()-*q2phi_h);	   
	  if(deltaPhi_2 > TMath::Pi()) deltaPhi_2 = 2*TMath::Pi()-fabs(jetAK4_2.Phi()-(*q2phi_h));

	  deltaRmin_2 = sqrt(deltaPhi_2*deltaPhi_2+fabs(jetAK4_2.Eta()-*q2eta_h)*fabs(jetAK4_2.Eta()-*q2eta_h));
	}
	else if(deltaRmin == deltaR_2){
	  deltaPhi_2 = fabs(jetAK4_2.Phi()-*q1phi_h);
          if(deltaPhi_2 > TMath::Pi()) deltaPhi_2 = 2*TMath::Pi()-fabs(jetAK4_2.Phi()-(*q1phi_h));

          deltaRmin_2 = sqrt(deltaPhi_2*deltaPhi_2+fabs(jetAK4_2.Eta()-*q1eta_h)*fabs(jetAK4_2.Eta()-*q1eta_h));
	}
	  
       if(deltaRmin_2 > 0.1) continue;

	numerator_resolved->Fill(*wzpt_h,bosonPt_W->GetBinContent(bosonPt_W->FindBin(*wzpt_h)));
	foundW = true;
	
      }
    }
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.12);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.28);
  pad2->SetTickx();
  pad2->SetTicky();


  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();

  TH1* frame = pad1->DrawFrame(numerator_boosted->GetXaxis()->GetBinLowEdge(1),0.,numerator_boosted->GetXaxis()->GetBinLowEdge(numerator_boosted->GetNbinsX()+1),1.2,"");
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.1);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleOffset(0.85);
  frame->Draw();

  CMS_lumi(canvas,"");

  TGraphAsymmErrors* boostedEff = new TGraphAsymmErrors();
  boostedEff->BayesDivide(numerator_boosted,denominator);
  boostedEff->SetMarkerColor(kBlack);
  boostedEff->SetLineColor(kBlack);
  boostedEff->SetMarkerStyle(20);
  boostedEff->SetMarkerSize(1.);
  boostedEff->Draw("PEsame");

  TGraphAsymmErrors* resolvedEff = new TGraphAsymmErrors();
  resolvedEff->BayesDivide(numerator_resolved,denominator);
  resolvedEff->SetMarkerColor(kRed);
  resolvedEff->SetLineColor(kRed);
  resolvedEff->SetMarkerStyle(20);
  resolvedEff->SetMarkerSize(1.);
  resolvedEff->Draw("PEsame");

  TLegend* leg = new TLegend(0.58, 0.36, 0.85, 0.62);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(boostedEff,"#varepsilon_{AK8} #DeltaR(W,jet) < 0.1","PLE");
  leg->AddEntry(resolvedEff,"#varepsilon_{AK4} #DeltaR(q_{i},jet) < 0.1","PLE");
  leg->Draw("same");

  pad1->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 =  pad2->DrawFrame(numerator_boosted->GetXaxis()->GetBinLowEdge(1),0.001,numerator_boosted->GetXaxis()->GetBinLowEdge(numerator_boosted->GetNbinsX()+1),1.5,"");
  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleSize(0.13);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleSize(0.10);
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->GetYaxis()->SetTitle("Resolved/Boosted");
  frame2->GetYaxis()->SetTitleOffset(0.5);
  frame2->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  frame2->Draw();

  numerator_boosted->Divide(denominator);
  numerator_resolved->Divide(denominator);
  numerator_resolved->Divide(numerator_boosted);

  numerator_resolved->SetStats(kFALSE);
  numerator_resolved->SetLineColor(kBlack);
  numerator_resolved->SetMarkerColor(kBlack);
  numerator_resolved->SetMarkerStyle(20);
  numerator_resolved->SetMarkerSize(1.0);
  numerator_resolved->SetLineWidth(2);

  numerator_resolved->Draw("hist same");
  pad2->SetLogy();

  pad2->RedrawAxis("sameaxis");
  

  system(("mkdir -p "+outputDir).c_str());
  
  if(isWboson){
    canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+"_Wboson.png").c_str(),"png");
    canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+"_Wboson.pdf").c_str(),"pdf");
  }
  else{
    canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+"_Zboson.png").c_str(),"png");
    canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+"_Zboson.pdf").c_str(),"pdf");
  }
}
