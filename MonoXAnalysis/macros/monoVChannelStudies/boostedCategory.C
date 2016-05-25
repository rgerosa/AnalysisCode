#include "../CMS_lumi.h"

void boostedCategory(string signalHypothesis, string outputDir){

  vector<double> bins = {200,225,250,300,350,400,450,500,550,600,650,700,750,800.,850.,900.,950.,1000.,1050,1100,1150,1200};

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  TFile* inputSignal_1 = NULL;
  TFile* inputSignal_2 = NULL;
  TFile* inputSignal_3 = NULL;
  TFile* inputSignal_4 = NULL;
  TFile* inputSignal_5 = NULL;
  
  if(signalHypothesis == "Scalar"){
    inputSignal_1 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-100_Mchi-50_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_2 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-1000_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_3 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-2000_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_4 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-10000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
  }
  else if(signalHypothesis == "Pseudoscalar"){
    inputSignal_1 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-100_Mchi-1_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_2 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_3 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
    inputSignal_4 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-2000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-JHUGen.root");
  }
  else if(signalHypothesis == "Axial"){
    inputSignal_1 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-300_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_2 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_3 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_4 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-2000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_5 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-10000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
  }
  else if(signalHypothesis == "Vector"){
    inputSignal_1 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-300_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_2 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_3 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_4 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-2000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
    inputSignal_5 = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-10000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph.root");
  }

  TH1F* bosonPt_D = new TH1F("bosonPt_D","",120,200,1200);
  bosonPt_D->Sumw2();

  TChain* chain = new TChain("tree/tree");
  chain->AddFile(inputSignal_1->GetName());
  chain->AddFile(inputSignal_2->GetName());
  chain->AddFile(inputSignal_3->GetName());
  chain->AddFile(inputSignal_4->GetName());
  if(inputSignal_5)
    chain->AddFile(inputSignal_5->GetName());

  chain->Draw("wzpt_h >> bosonPt_D","wzpt_h > 200 && abs(wzeta_h) < 2.4 && hltmet90 && njets >= 1 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8&& incjetmumetdphimin4 > 0.5 && centraljetpt[0] > 100 && t1pfmet > 200 ","goff");
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
  TTreeReaderValue<double>       jmmdphi   (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double>       met       (myReader,"t1pfmet");

  TTreeReaderValue<vector<double> > jetpt  (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jeteta (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetm (myReader,"centraljetm");
  TTreeReaderValue<vector<double> > chfrac (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac (myReader,"centraljetNHfrac");
  TTreeReaderValue<vector<double> > genjetpt  (myReader,"centraljetGenpt");
  TTreeReaderValue<vector<double> > genjeteta (myReader,"centraljetGeneta");
  TTreeReaderValue<vector<double> > genjetphi (myReader,"centraljetGenphi");
  TTreeReaderValue<vector<double> > genjetm (myReader,"centraljetGenm");
  

  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > boostedJetGenpt    (myReader,"boostedJetGenpt");
  TTreeReaderValue<vector<double> > boostedJetGeneta   (myReader,"boostedJetGeneta");
  TTreeReaderValue<vector<double> > boostedJetGenphi   (myReader,"boostedJetGenphi");
  TTreeReaderValue<vector<double> > boostedJetGenm     (myReader,"boostedJetGenm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<double> > boostedJetBosoneta  (myReader,"boostedJetBosoneta");
  TTreeReaderValue<vector<double> > boostedJetBosonphi  (myReader,"boostedJetBosonphi");
  TTreeReaderValue<vector<double> > boostedJetBosonpt   (myReader,"boostedJetBosonpt");
  TTreeReaderValue<vector<double> > boostedJetBosonm    (myReader,"boostedJetBosonm");
  

  TTreeReaderValue<double > wzpt_h (myReader,"wzpt_h");
  TTreeReaderValue<double > wzeta_h (myReader,"wzeta_h");
  TTreeReaderValue<double > wzphi_h (myReader,"wzphi_h");
  TTreeReaderValue<double > wzmass_h (myReader,"wzmass_h");
  TTreeReaderValue<double > q1pt_h (myReader,"q1pt");
  TTreeReaderValue<double > q1eta_h (myReader,"q1eta");
  TTreeReaderValue<double > q1phi_h (myReader,"q1phi");
  TTreeReaderValue<double > q2pt_h (myReader,"q2pt");
  TTreeReaderValue<double > q2eta_h (myReader,"q2eta");
  TTreeReaderValue<double > q2phi_h (myReader,"q2phi");

  TH1F* denominator        = new TH1F("den_boosted","",bins.size()-1,&bins[0]);
  TH1F* numerator_boosted  = new TH1F("numerator_boosted","",bins.size()-1,&bins[0]);
  TH1F* numerator_resolved = new TH1F("numerator_resolved","",bins.size()-1,&bins[0]);

  denominator->Sumw2();
  numerator_boosted->Sumw2();
  numerator_resolved->Sumw2();

  while(myReader.Next()){
    if(*hlt == 0) continue;
    if(*njets < 1) continue;
    if(*nbjets > 0) continue;
    if(*jmmdphi < 0.5) continue;
    if(*met < 200) continue;
    if(jetpt->size() <= 0) continue;
    if(nhfrac->at(0) > 0.8) continue;
    if(chfrac->at(0) < 0.1) continue;
    if(jetpt->at(0) < 100) continue;

    // ensuring one hadronic W in the eta acceptance, with a GenPt > 200 GeV
    if(*wzpt_h <= 200. or fabs(*wzeta_h) > 2.4 or fabs(*q1eta_h) > 2.4 or fabs(*q2eta_h) > 2.4 or fabs(*q2pt_h) < 30 or fabs(*q1pt_h) < 30) continue;

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

  TCanvas* canvas = new TCanvas("canvas","canvas",600,700);
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

  TH1* frame = pad1->DrawFrame(numerator_boosted->GetXaxis()->GetBinLowEdge(1),0.,numerator_boosted->GetXaxis()->GetBinLowEdge(numerator_boosted->GetNbinsX()+1),1.,"");
  frame->GetYaxis()->SetTitle("efficiency");
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->Draw();

  CMS_lumi(pad1,"2.30",true);

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

  TH1* frame2 =  pad2->DrawFrame(numerator_boosted->GetXaxis()->GetBinLowEdge(1),0.,numerator_boosted->GetXaxis()->GetBinLowEdge(numerator_boosted->GetNbinsX()+1),1.,"");
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

  numerator_resolved->Draw("PEsame");

  pad2->RedrawAxis("sameaxis");
  

  system(("mkdir -p "+outputDir).c_str());
  
  canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+".png").c_str(),"png");
  canvas->SaveAs((outputDir+"/efficiency_"+signalHypothesis+".pdf").c_str(),"pdf");
}
