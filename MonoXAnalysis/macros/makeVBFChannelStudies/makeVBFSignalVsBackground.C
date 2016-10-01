#include "../CMS_lumi.h"

void fillHisto(TH1* histo, double val, double weight = 1){ // embed the overflow
  
  if(val < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    histo->Fill(val,weight);
  else
    histo->Fill(histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1)-1,weight);

}

void fillHisto(TH2* histo, double valx, double valy, double weight = 1){ // Embed- the overflow
  
  double x = 0;
  if(valx < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    x = valx;
  else
    x = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1)-1;
  double y = 0;
  if(valy < histo->GetYaxis()->GetBinLowEdge(histo->GetNbinsY()+1))
    y = valy;
  else
    y = histo->GetYaxis()->GetBinLowEdge(histo->GetNbinsY()+1)-1;

  histo->Fill(x,y,weight);
}


void plotDistribution(TCanvas* canvas, TH1* qqH, TH1* znn, TH1* znn_ewk, const string & label, const string & outputDIR){
  
  // normalize to 1
  canvas->cd();
  qqH->Scale(1./qqH->Integral());
  znn->Scale(1./znn->Integral());
  znn_ewk->Scale(1./znn_ewk->Integral());
  qqH->GetXaxis()->SetTitle(label.c_str());
  qqH->GetYaxis()->SetTitle("a.u");
  qqH->SetLineColor(kBlack);
  qqH->SetLineWidth(2);
  qqH->Draw("hist");
  
  znn->SetLineColor(kRed);
  znn->SetLineWidth(2);
  znn->Draw("hist same");
  znn_ewk->SetLineColor(kBlue);
  znn_ewk->SetLineWidth(2);
  znn_ewk->Draw("hist same");

  qqH->GetYaxis()->SetRangeUser(0.,max(qqH->GetMaximum(),max(znn->GetMaximum(),znn_ewk->GetMaximum()))*1.2);
  CMS_lumi(canvas,"");
  TLegend leg(0.7,0.7,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(qqH,"VBF H = 125 GeV","L");
  leg.AddEntry(znn,"Z#nu#nu","L");
  leg.AddEntry(znn_ewk,"Z#nu#nu (EWK)","L");
  leg.Draw("same");
  
  //  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".pdf").c_str(),"pdf");

  qqH->GetYaxis()->SetRangeUser(max(0.0001,min(qqH->GetMinimum(),znn->GetMinimum())*0.8),max(qqH->GetMaximum(),znn->GetMaximum())*1000);
  canvas->SetLogy();
  //  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  
}



void plotDistribution(TCanvas* canvas, TH2* qqH, TH2* znn, TH2* znn_ewk, const string & labelX, const string & labelY, const string & outputDIR){
  
  system(("mkdir -p "+outputDIR+"/correlation/").c_str());
  // normalize to 1
  canvas->cd();
  canvas->SetRightMargin(0.18);
  qqH->Scale(1./qqH->Integral());
  znn->Scale(1./znn->Integral());
  znn_ewk->Scale(1./znn_ewk->Integral());

  TGraph2D* qqHGraph = new TGraph2D();
  qqHGraph->SetNpx(300);
  qqHGraph->SetNpy(300);
  int nPoint = 0;
  for(int iBinX = 0; iBinX < qqH->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < qqH->GetNbinsY() ; iBinY++){
      qqHGraph->SetPoint(nPoint,qqH->GetXaxis()->GetBinCenter(iBinX+1),qqH->GetYaxis()->GetBinCenter(iBinY+1),qqH->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  TH2* qqHPlot = qqHGraph->GetHistogram();
  qqHPlot->GetXaxis()->SetTitle(labelX.c_str());
  qqHPlot->GetYaxis()->SetTitle(labelY.c_str());
  qqHPlot->GetZaxis()->SetTitle("a.u");   
  qqHPlot->Draw("colz");

  TProfile* qqHProfile = qqH->ProfileX(Form("%s_pfx",qqH->GetName()));
  qqHProfile->SetMarkerColor(kBlack);
  qqHProfile->SetMarkerStyle(20);
  qqHProfile->SetMarkerSize(1);
  qqHProfile->Draw("EPsame");

  CMS_lumi(canvas,"");
  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,"VBF H = 125 GeV","");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",qqHPlot->GetCorrelationFactor()),"");
  leg.Draw("same");

  //  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  //  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);

  TGraph2D* znnGraph = new TGraph2D();
  znnGraph->SetNpx(300);
  znnGraph->SetNpy(300);
  nPoint = 0;
  for(int iBinX = 0; iBinX < znn->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < znn->GetNbinsY() ; iBinY++){
      znnGraph->SetPoint(nPoint,znn->GetXaxis()->GetBinCenter(iBinX+1),znn->GetYaxis()->GetBinCenter(iBinY+1),znn->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  TH2* znnPlot = znnGraph->GetHistogram();
  znnPlot->GetXaxis()->SetTitle(labelX.c_str());
  znnPlot->GetYaxis()->SetTitle(labelY.c_str());
  znnPlot->GetZaxis()->SetTitle("a.u");   
  znnPlot->Draw("colz");
  CMS_lumi(canvas,"");

  TProfile* znnProfile = znn->ProfileX(Form("%s_pfx",znn->GetName()));
  znnProfile->SetMarkerColor(kBlack);
  znnProfile->SetMarkerStyle(20);
  znnProfile->SetMarkerSize(1);
  znnProfile->Draw("EPsame");
 
  leg.Clear();
  leg.AddEntry((TObject*)0,"Z#nu#nu","");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",znnPlot->GetCorrelationFactor()),"");
  leg.Draw();

  //  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(znn->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  //  canvas->SaveAs((outputDIR+"/correlation/"+string(znn->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(znn->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);

  TGraph2D* znn_ewkGraph = new TGraph2D();
  znn_ewkGraph->SetNpx(300);
  znn_ewkGraph->SetNpy(300);
  nPoint = 0;
  for(int iBinX = 0; iBinX < znn_ewk->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < znn_ewk->GetNbinsY() ; iBinY++){
      znn_ewkGraph->SetPoint(nPoint,znn_ewk->GetXaxis()->GetBinCenter(iBinX+1),znn_ewk->GetYaxis()->GetBinCenter(iBinY+1),znn_ewk->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  TH2* znn_ewkPlot = znn_ewkGraph->GetHistogram();
  znn_ewkPlot->GetXaxis()->SetTitle(labelX.c_str());
  znn_ewkPlot->GetYaxis()->SetTitle(labelY.c_str());
  znn_ewkPlot->GetZaxis()->SetTitle("a.u");   
  znn_ewkPlot->Draw("colz");
  CMS_lumi(canvas,"");

  TProfile* znn_ewkProfile = znn_ewk->ProfileX(Form("%s_pfx",znn_ewk->GetName()));
  znn_ewkProfile->SetMarkerColor(kBlack);
  znn_ewkProfile->SetMarkerStyle(20);
  znn_ewkProfile->SetMarkerSize(1);
  znn_ewkProfile->Draw("EPsame");
 
  leg.Clear();
  leg.AddEntry((TObject*)0,"Z#nu#nu EWK","");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",znn_ewkPlot->GetCorrelationFactor()),"");
  leg.Draw();

  //  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(znn_ewk->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  //  canvas->SaveAs((outputDIR+"/correlation/"+string(znn_ewk->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(znn_ewk->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);
  
}

void makeVBFSignalVSBackground(string outputPlot){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputPlot).c_str());

  TChain* vbftree = new TChain("tree/tree");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*110*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*125*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*150*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*200*root");

  TChain* znntree = new TChain("tree/tree");
  znntree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-*root");

  TChain* znnewktree = new TChain("tree/tree");
  znnewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJetsToNuNuEWK/sigfilter/sig_EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8.root");

  // get k-factors NLO                                                                                                                                                                                
  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  vector<TH1*> khists;
  khists.push_back(znlohist);
  khists.push_back(zewkhist);
  
  TTreeReader myReader (vbftree);
  TTreeReaderValue<double> xsec         (myReader,"xsec");
  TTreeReaderValue<double> wgt          (myReader,"wgt");
  TTreeReaderValue<double> wgtsum       (myReader,"wgtsum");
  TTreeReaderValue<unsigned int> njetsinc  (myReader,"njetsinc");
  TTreeReaderValue<unsigned int> nphotons  (myReader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons   (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus   (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets  (myReader,"nbjetslowpt");
  TTreeReaderValue<UChar_t> fhbhe   (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso  (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc    (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb    (myReader,"flageebadsc");
  TTreeReaderValue<vector<double> > jetpt    (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetbtag  (myReader,"combinejetbtag");
  TTreeReaderValue<vector<double> > jetm     (myReader,"combinejetm");
  TTreeReaderValue<vector<double> > chfrac   (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac   (myReader,"combinejetNHfrac");
  TTreeReaderValue<double > mediatorMass     (myReader,"dmmass");
  TTreeReaderValue<double > mediatorPt       (myReader,"dmpt");
  TTreeReaderValue<double > mediatorEta      (myReader,"dmeta");
  TTreeReaderValue<double > mediatorPhi      (myReader,"dmphi");
  TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> wzpt        (myReader,"wzpt");
  TTreeReaderValue<double> wzeta       (myReader,"wzeta");
  TTreeReaderValue<double> wzphi       (myReader,"wzphi");
  TTreeReaderValue<double> wzmass      (myReader,"wzmass");
  TTreeReaderValue<double> genmet      (myReader,"genmet");
  TTreeReaderValue<double> genmetphi   (myReader,"genmetphi");
  TTreeReaderValue<double> jmmdphi4    (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jmmdphi     (myReader,"incjetmumetdphimin");  
  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  
  TH1F* qqH_bosonPt = new TH1F("qqH_bosonPt","qqH_bosonPt",25,150,1250);
  TH1F* qqH_bosonEta = new TH1F("qqH_bosonEta","qqH_bosonEta",25,-3.5,3.5);
  TH1F* qqH_bosonPhi = new TH1F("qqH_bosonPhi","qqH_bosonPhi",25,0.,3.14);
  TH1F* qqH_bosonPz    = new TH1F("qqH_bosonPz","qqH_bosonPz",25,150,1250);
  TH1F* qqH_bosonPzOPt = new TH1F("qqH_bosonPzOPt","qqH_bosonPzOPt",20,0,20);
  TH1F* qqH_jetpt1   = new TH1F("qqH_jetpt1","qqH_jetpt1",25,50,500);
  TH1F* qqH_jetpt2   = new TH1F("qqH_jetpt2","qqH_jetpt2",25,50,500);
  TH1F* qqH_jetpz1   = new TH1F("qqH_jetpz1","qqH_jetpz1",25,50,500);  
  TH1F* qqH_jetpz2   = new TH1F("qqH_jetpz2","qqH_jetpz2",25,50,500);
  TH1F* qqH_jetpz1Opt   = new TH1F("qqH_jetpz1Opt","qqH_jetpz1",20,0,20);  
  TH1F* qqH_jetpz2Opt   = new TH1F("qqH_jetpz2OPt","qqH_jetpz2",20,0,20);
  TH1F* qqH_jeteta1  = new TH1F("qqH_jeteta1","qqH_jeteta1",25,-3.5,3.5);
  TH1F* qqH_jeteta2  = new TH1F("qqH_jeteta2","qqH_jeteta2",25,-3.5,3.5);
  TH1F* qqH_jetphi1  = new TH1F("qqH_jetphi1","qqH_jetphi1",25,0,3.14);
  TH1F* qqH_jetphi2  = new TH1F("qqH_jetphi2","qqH_jetphi2",25,0,3.14);
  TH1F* qqH_jetDeltaPhi  = new TH1F("qqH_jetDeltaPhi","qqH_jetDeltaPhi",25,0,3.14);
  TH1F* qqH_jetDeltaEta  = new TH1F("qqH_jetDeltaEta","qqH_jetDeltaEta",25,0,9);
  TH1F* qqH_Mjj     = new TH1F("qqH_Mjj","qqH_Mjj",40,400,2500);
  TH1F* qqH_jetmetDphi4 = new TH1F("qqH_jetmetdphi4","qqH_jetmetdphi4",25,0,3.14);
  TH1F* qqH_jetmetDphi  = new TH1F("qqH_jetmetdphi","qqH_jetmetdphi",25,0,3.14);
  TH1F* qqH_jetmediatorDphi = new TH1F("qqH_jetmediatorDphi","qqH_jetmediatorDphi",25,0,3.14);
  TH1F* qqH_jetmediatorDeta = new TH1F("qqH_jetmediatorDeta","qqH_jetmediatorDeta",25,0,9);
  TH1F* qqH_MT = new TH1F("qqH_MT","qqH_MT",30,0,600);
  TH1F* qqH_MetMT = new TH1F("qqH_MetMT","qqH_MetMT",30,100,1000);

  qqH_bosonPt->Sumw2();
  qqH_bosonEta->Sumw2();
  qqH_bosonPhi->Sumw2();
  qqH_bosonPz->Sumw2();
  qqH_bosonPzOPt->Sumw2();
  qqH_jetpt1->Sumw2();
  qqH_jetpt2->Sumw2();
  qqH_jetpz1->Sumw2();
  qqH_jetpz2->Sumw2();
  qqH_jetpz1Opt->Sumw2();
  qqH_jetpz2Opt->Sumw2();
  qqH_jeteta1->Sumw2();
  qqH_jeteta2->Sumw2();
  qqH_jetphi1->Sumw2();
  qqH_jetphi2->Sumw2();
  qqH_jetDeltaPhi->Sumw2();
  qqH_jetDeltaEta->Sumw2();
  qqH_Mjj->Sumw2();
  qqH_jetmetDphi->Sumw2();
  qqH_jetmetDphi4->Sumw2();
  qqH_jetmediatorDphi->Sumw2();
  qqH_jetmediatorDeta->Sumw2();
  qqH_MT->Sumw2();
  qqH_MetMT->Sumw2();

  TH2F* qqH_mjj_detajj = new TH2F("qqH_mjj_detajj","qqH_mjj_detajj",15,400,2500,10,2,8);
  TH2F* qqH_mjj_jetmetDphi = new TH2F("qqH_mjj_jetmetDphi","qqH_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_mjj_jetDphi = new TH2F("qqH_mjj_jetDphi","qqH_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_mjj_MT = new TH2F("qqH_mjj_MT","qqH_mjj_MT",15,400,2500,10,0,500);
  TH2F* qqH_detajj_jetmetDphi = new TH2F("qqH_detajj_jetmetDphi","qqH_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_detajj_jetDphi = new TH2F("qqH_detajj_jetDphi","qqH_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_detajj_MT = new TH2F("qqH_detajj_MT","qqH_detajj_MT",10,2,8,10,0,500);

  qqH_mjj_detajj->Sumw2();
  qqH_mjj_jetmetDphi->Sumw2();
  qqH_mjj_jetDphi->Sumw2();
  qqH_mjj_MT->Sumw2();
  qqH_detajj_jetmetDphi->Sumw2();
  qqH_detajj_jetDphi->Sumw2();
  qqH_detajj_MT->Sumw2();

  ////////////

  while(myReader.Next()){
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets
    if (*njetsinc  < 2) continue;
    // b-veto
    if (*nbjets > 0) continue;
    // vetos
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 50) continue; 
    if (jetpt->at(1) < 50) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut
    if (*mmet     < 200) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2) continue;
    if ((jet1+jet2).M() < 450) continue;
    
    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*mediatorPt,*mediatorEta,*mediatorPhi,*mediatorMass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);
    double weight = *xsec*(*wgt)/(*wgtsum);
    // fill some hisotgrams
    fillHisto(qqH_bosonPt,*mediatorPt,weight);
    fillHisto(qqH_bosonEta,*mediatorEta,weight);
    fillHisto(qqH_bosonPhi,fabs(*mediatorPhi),weight);
    fillHisto(qqH_bosonPz,mediator.Pz(),weight);
    fillHisto(qqH_bosonPzOPt,mediator.Pz()/mediator.Pt(),weight);
    fillHisto(qqH_jetpt1,jetpt->at(0),weight);
    fillHisto(qqH_jetpt2,jetpt->at(1),weight);
    fillHisto(qqH_jetpz1,jet1.Pz(),weight);
    fillHisto(qqH_jetpz2,jet2.Pz(),weight);
    fillHisto(qqH_jetpz1Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(qqH_jetpz2Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(qqH_jeteta1,jeteta->at(0),weight);
    fillHisto(qqH_jeteta2,jeteta->at(1),weight);
    fillHisto(qqH_jetphi1,fabs(jetphi->at(0)),weight);
    fillHisto(qqH_jetphi2,fabs(jetphi->at(1)),weight);
    fillHisto(qqH_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(qqH_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(qqH_Mjj,(jet1+jet2).M(),weight);
    fillHisto(qqH_jetmetDphi,*jmmdphi,weight);
    fillHisto(qqH_jetmetDphi4,*jmmdphi4,weight);
    fillHisto(qqH_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)),weight);
    fillHisto(qqH_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()),weight);
    fillHisto(qqH_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(qqH_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))),weight);

    // 2D part
    fillHisto(qqH_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(qqH_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi,weight);
    fillHisto(qqH_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(qqH_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(qqH_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi,weight);
    fillHisto(qqH_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(qqH_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);

  }

  ////  
  TH1F* znnewk_bosonPt = new TH1F("znnewk_bosonPt","znnewk_bosonPt",25,150,1250);
  TH1F* znnewk_bosonEta = new TH1F("znnewk_bosonEta","znnewk_bosonEta",25,-3.5,3.5);
  TH1F* znnewk_bosonPhi = new TH1F("znnewk_bosonPhi","znnewk_bosonPhi",25,0.,3.14);
  TH1F* znnewk_bosonPz    = new TH1F("znnewk_bosonPz","znnewk_bosonPz",25,150,1250);
  TH1F* znnewk_bosonPzOPt = new TH1F("znnewk_bosonPzOPt","znnewk_bosonPzOPt",20,0,20);
  TH1F* znnewk_jetpt1   = new TH1F("znnewk_jetpt1","znnewk_jetpt1",25,50,500);
  TH1F* znnewk_jetpt2   = new TH1F("znnewk_jetpt2","znnewk_jetpt2",25,50,500);
  TH1F* znnewk_jetpz1   = new TH1F("znnewk_jetpz1","znnewk_jetpz1",25,50,500);  
  TH1F* znnewk_jetpz2   = new TH1F("znnewk_jetpz2","znnewk_jetpz2",25,50,500);
  TH1F* znnewk_jetpz1Opt   = new TH1F("znnewk_jetpz1Opt","znnewk_jetpz1",20,0,20);  
  TH1F* znnewk_jetpz2Opt   = new TH1F("znnewk_jetpz2OPt","znnewk_jetpz2",20,0,20);
  TH1F* znnewk_jeteta1  = new TH1F("znnewk_jeteta1","znnewk_jeteta1",25,-3.5,3.5);
  TH1F* znnewk_jeteta2  = new TH1F("znnewk_jeteta2","znnewk_jeteta2",25,-3.5,3.5);
  TH1F* znnewk_jetphi1  = new TH1F("znnewk_jetphi1","znnewk_jetphi1",25,0,3.14);
  TH1F* znnewk_jetphi2  = new TH1F("znnewk_jetphi2","znnewk_jetphi2",25,0,3.14);
  TH1F* znnewk_jetDeltaPhi  = new TH1F("znnewk_jetDeltaPhi","znnewk_jetDeltaPhi",25,0,3.14);
  TH1F* znnewk_jetDeltaEta  = new TH1F("znnewk_jetDeltaEta","znnewk_jetDeltaEta",25,0,9);
  TH1F* znnewk_Mjj     = new TH1F("znnewk_Mjj","znnewk_Mjj",40,400,2500);
  TH1F* znnewk_jetmetDphi4 = new TH1F("znnewk_jetmetdphi4","znnewk_jetmetdphi4",25,0,3.14);
  TH1F* znnewk_jetmetDphi  = new TH1F("znnewk_jetmetdphi","znnewk_jetmetdphi",25,0,3.14);
  TH1F* znnewk_jetmediatorDphi = new TH1F("znnewk_jetmediatorDphi","znnewk_jetmediatorDphi",25,0,3.14);
  TH1F* znnewk_jetmediatorDeta = new TH1F("znnewk_jetmediatorDeta","znnewk_jetmediatorDeta",25,0,9);
  TH1F* znnewk_MT = new TH1F("znnewk_MT","znnewk_MT",30,0,600);
  TH1F* znnewk_MetMT = new TH1F("znnewk_MetMT","znnewk_MetMT",30,100,1000);

  znnewk_bosonPt->Sumw2();
  znnewk_bosonEta->Sumw2();
  znnewk_bosonPhi->Sumw2();
  znnewk_bosonPz->Sumw2();
  znnewk_bosonPzOPt->Sumw2();
  znnewk_jetpt1->Sumw2();
  znnewk_jetpt2->Sumw2();
  znnewk_jetpz1->Sumw2();
  znnewk_jetpz2->Sumw2();
  znnewk_jetpz1Opt->Sumw2();
  znnewk_jetpz2Opt->Sumw2();
  znnewk_jeteta1->Sumw2();
  znnewk_jeteta2->Sumw2();
  znnewk_jetphi1->Sumw2();
  znnewk_jetphi2->Sumw2();
  znnewk_jetDeltaPhi->Sumw2();
  znnewk_jetDeltaEta->Sumw2();
  znnewk_Mjj->Sumw2();
  znnewk_jetmetDphi->Sumw2();
  znnewk_jetmetDphi4->Sumw2();
  znnewk_jetmediatorDphi->Sumw2();
  znnewk_jetmediatorDeta->Sumw2();
  znnewk_MT->Sumw2();
  znnewk_MetMT->Sumw2();

  TH2F* znnewk_mjj_detajj = new TH2F("znnewk_mjj_detajj","znnewk_mjj_detajj",15,400,2500,10,2,8);
  TH2F* znnewk_mjj_jetmetDphi = new TH2F("znnewk_mjj_jetmetDphi","znnewk_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* znnewk_mjj_jetDphi = new TH2F("znnewk_mjj_jetDphi","znnewk_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* znnewk_mjj_MT = new TH2F("znnewk_mjj_MT","znnewk_mjj_MT",15,400,2500,10,0,500);
  TH2F* znnewk_detajj_jetmetDphi = new TH2F("znnewk_detajj_jetmetDphi","znnewk_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* znnewk_detajj_jetDphi = new TH2F("znnewk_detajj_jetDphi","znnewk_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* znnewk_detajj_MT = new TH2F("znnewk_detajj_MT","znnewk_detajj_MT",10,2,8,10,0,500);

  znnewk_mjj_detajj->Sumw2();
  znnewk_mjj_jetmetDphi->Sumw2();
  znnewk_mjj_jetDphi->Sumw2();
  znnewk_mjj_MT->Sumw2();
  znnewk_detajj_jetmetDphi->Sumw2();
  znnewk_detajj_jetDphi->Sumw2();
  znnewk_detajj_MT->Sumw2();


  ////////////
  myReader.SetTree(znnewktree);
  myReader.SetEntry(0);
  while(myReader.Next()){
    
    // met filters                                                                                                                                                                                    
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets                                                                                                                                                                                 
    if (*njetsinc  < 2) continue;
    // b-veto                                                                                                                                                                                         
    if (*nbjets > 0) continue;
    // vetos                                                                                                                                                                                          
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts                                                                                                                                                                        
    if (jetpt->at(0) < 50) continue;
    if (jetpt->at(1) < 50) continue;
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut                                                                                                                                                                                
    if (*mmet     < 130) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts                                                                                                                                                                               
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2) continue;
    if ((jet1+jet2).M() < 450) continue;

    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*wzpt,*wzeta,*wzphi,*wzmass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);

    double weight = *xsec*(*wgt)/(*wgtsum);
    // fill some hisotgrams
    fillHisto(znnewk_bosonPt,*wzpt,weight);
    fillHisto(znnewk_bosonEta,*wzeta,weight);
    fillHisto(znnewk_bosonPhi,fabs(*wzphi),weight);
    fillHisto(znnewk_bosonPz,mediator.Pz(),weight);
    fillHisto(znnewk_bosonPzOPt,mediator.Pz()/mediator.Pt(),weight);
    fillHisto(znnewk_jetpt1,jetpt->at(0),weight);
    fillHisto(znnewk_jetpt2,jetpt->at(1),weight);
    fillHisto(znnewk_jetpz1,jet1.Pz(),weight);
    fillHisto(znnewk_jetpz2,jet2.Pz(),weight);
    fillHisto(znnewk_jetpz1Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(znnewk_jetpz2Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(znnewk_jeteta1,jeteta->at(0),weight);
    fillHisto(znnewk_jeteta2,jeteta->at(1),weight);
    fillHisto(znnewk_jetphi1,fabs(jetphi->at(0)),weight);
    fillHisto(znnewk_jetphi2,fabs(jetphi->at(1)),weight);
    fillHisto(znnewk_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znnewk_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(znnewk_Mjj,(jet1+jet2).M(),weight);
    fillHisto(znnewk_jetmetDphi,*jmmdphi,weight);
    fillHisto(znnewk_jetmetDphi4,*jmmdphi4,weight);
    fillHisto(znnewk_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)),weight);
    fillHisto(znnewk_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()),weight);
    fillHisto(znnewk_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(znnewk_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))),weight);

    fillHisto(znnewk_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(znnewk_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi,weight);
    fillHisto(znnewk_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znnewk_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(znnewk_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi,weight);
    fillHisto(znnewk_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znnewk_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);  
  }

  /// Znunu QCD
  myReader.SetTree(znntree);
  myReader.SetEntry(0);
  
  TH1F* znn_bosonPt = new TH1F("znn_bosonPt","znn_bosonPt",25,150,1250);
  TH1F* znn_bosonEta = new TH1F("znn_bosonEta","znn_bosonEta",25,-3.5,3.5);
  TH1F* znn_bosonPhi = new TH1F("znn_bosonPhi","znn_bosonPhi",25,0.,3.14);
  TH1F* znn_bosonPz    = new TH1F("znn_bosonPz","znn_bosonPz",25,150,1250);
  TH1F* znn_bosonPzOPt = new TH1F("znn_bosonPzOPt","znn_bosonPzOPt",20,0,20);
  TH1F* znn_jetpt1   = new TH1F("znn_jetpt1","znn_jetpt1",25,50,500);
  TH1F* znn_jetpt2   = new TH1F("znn_jetpt2","znn_jetpt2",25,50,500);
  TH1F* znn_jetpz1   = new TH1F("znn_jetpz1","znn_jetpz1",25,50,500);  
  TH1F* znn_jetpz2   = new TH1F("znn_jetpz2","znn_jetpz2",25,50,500);
  TH1F* znn_jetpz1Opt   = new TH1F("znn_jetpz1Opt","znn_jetpz1",20,0,20);  
  TH1F* znn_jetpz2Opt   = new TH1F("znn_jetpz2OPt","znn_jetpz2",20,0,20);
  TH1F* znn_jeteta1  = new TH1F("znn_jeteta1","znn_jeteta1",25,-3.5,3.5);
  TH1F* znn_jeteta2  = new TH1F("znn_jeteta2","znn_jeteta2",25,-3.5,3.5);
  TH1F* znn_jetphi1  = new TH1F("znn_jetphi1","znn_jetphi1",25,0,3.14);
  TH1F* znn_jetphi2  = new TH1F("znn_jetphi2","znn_jetphi2",25,0,3.14);
  TH1F* znn_jetDeltaPhi  = new TH1F("znn_jetDeltaPhi","znn_jetDeltaPhi",25,0,3.14);
  TH1F* znn_jetDeltaEta  = new TH1F("znn_jetDeltaEta","znn_jetDeltaEta",25,0,9);
  TH1F* znn_Mjj     = new TH1F("znn_Mjj","znn_Mjj",40,400,2500);
  TH1F* znn_jetmetDphi4 = new TH1F("znn_jetmetdphi4","znn_jetmetdphi4",25,0,3.14);
  TH1F* znn_jetmetDphi  = new TH1F("znn_jetmetdphi","znn_jetmetdphi",25,0,3.14);
  TH1F* znn_jetmediatorDphi = new TH1F("znn_jetmediatorDphi","znn_jetmediatorDphi",25,0,3.14);
  TH1F* znn_jetmediatorDeta = new TH1F("znn_jetmediatorDeta","znn_jetmediatorDeta",25,0,9);
  TH1F* znn_MT = new TH1F("znn_MT","znn_MT",30,0,600);
  TH1F* znn_MetMT = new TH1F("znn_MetMT","znn_MetMT",30,100,1000);

  znn_bosonPt->Sumw2();
  znn_bosonEta->Sumw2();
  znn_bosonPhi->Sumw2();
  znn_bosonPz->Sumw2();
  znn_bosonPzOPt->Sumw2();
  znn_jetpt1->Sumw2();
  znn_jetpt2->Sumw2();
  znn_jetpz1->Sumw2();
  znn_jetpz2->Sumw2();
  znn_jetpz1Opt->Sumw2();
  znn_jetpz2Opt->Sumw2();
  znn_jeteta1->Sumw2();
  znn_jeteta2->Sumw2();
  znn_jetphi1->Sumw2();
  znn_jetphi2->Sumw2();
  znn_jetDeltaPhi->Sumw2();
  znn_jetDeltaEta->Sumw2();
  znn_Mjj->Sumw2();
  znn_jetmetDphi->Sumw2();
  znn_jetmetDphi4->Sumw2();
  znn_jetmediatorDphi->Sumw2();
  znn_jetmediatorDeta->Sumw2();
  znn_MT->Sumw2();
  znn_MetMT->Sumw2();

  TH2F* znn_mjj_detajj = new TH2F("znn_mjj_detajj","znn_mjj_detajj",15,400,2500,10,2,8);
  TH2F* znn_mjj_jetmetDphi = new TH2F("znn_mjj_jetmetDphi","znn_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* znn_mjj_jetDphi = new TH2F("znn_mjj_jetDphi","znn_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* znn_mjj_MT = new TH2F("znn_mjj_MT","znn_mjj_MT",15,400,2500,10,0,500);
  TH2F* znn_detajj_jetmetDphi = new TH2F("znn_detajj_jetmetDphi","znn_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* znn_detajj_jetDphi = new TH2F("znn_detajj_jetDphi","znn_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* znn_detajj_MT = new TH2F("znn_detajj_MT","znn_detajj_MT",10,2,8,10,0,500);

  znn_mjj_detajj->Sumw2();
  znn_mjj_jetmetDphi->Sumw2();
  znn_mjj_jetDphi->Sumw2();
  znn_mjj_MT->Sumw2();
  znn_detajj_jetmetDphi->Sumw2();
  znn_detajj_jetDphi->Sumw2();
  znn_detajj_MT->Sumw2();

  ////////////
  while(myReader.Next()){

    // met filters                                                                                                                                                                                    
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets                                                                                                                                                                                 
    if (*njetsinc  < 2) continue;
    // b-veto                                                                                                                                                                                         
    if (*nbjets > 0) continue;
    // vetos                                                                                                                                                                                          
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts                                                                                                                                                                        
    if (jetpt->at(0) < 50) continue;
    if (jetpt->at(1) < 50) continue;
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut                                                                                                                                                                                
    if (*mmet     < 130) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts                                                                                                                                                                               
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2.0) continue;
    if ((jet1+jet2).M() < 450) continue;

    double kwgt = 1;
    double genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
        if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
        if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
        kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }    

    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*wzpt,*wzeta,*wzphi,*wzmass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);
    double weight = *xsec*(*wgt)*(kwgt)/(*wgtsum);

    // fill some hisotgrams
    fillHisto(znn_bosonPt,*wzpt,weight);
    fillHisto(znn_bosonEta,*wzeta,weight);
    fillHisto(znn_bosonPhi,fabs(*wzphi),weight);
    fillHisto(znn_bosonPz,mediator.Pz(),weight);
    fillHisto(znn_bosonPzOPt,mediator.Pz()/mediator.Pt(),weight);
    fillHisto(znn_jetpt1,jetpt->at(0),weight);
    fillHisto(znn_jetpt2,jetpt->at(1),weight);
    fillHisto(znn_jetpz1,jet1.Pz(),weight);
    fillHisto(znn_jetpz2,jet2.Pz(),weight);
    fillHisto(znn_jetpz1Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(znn_jetpz2Opt,jet1.Pz()/jet1.Pt(),weight);
    fillHisto(znn_jeteta1,jeteta->at(0),weight);
    fillHisto(znn_jeteta2,jeteta->at(1),weight);
    fillHisto(znn_jetphi1,fabs(jetphi->at(0)),weight);
    fillHisto(znn_jetphi2,fabs(jetphi->at(1)),weight);
    fillHisto(znn_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znn_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(znn_Mjj,(jet1+jet2).M(),weight);
    fillHisto(znn_jetmetDphi,*jmmdphi,weight);
    fillHisto(znn_jetmetDphi4,*jmmdphi4,weight);
    fillHisto(znn_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)),weight);
    fillHisto(znn_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()),weight);
    fillHisto(znn_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(znn_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))),weight);

    fillHisto(znn_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()),weight);
    fillHisto(znn_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi,weight);
    fillHisto(znn_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znn_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);    
    fillHisto(znn_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi,weight);
    fillHisto(znn_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)),weight);
    fillHisto(znn_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))),weight);  
  }

  cout<<"Expected signal    rate normalized to 1 fb-1: "<<qqH_bosonPt->Integral()<<endl;
  cout<<"Expected EWK Z     rate normalized to 1 fb-1: "<<znnewk_bosonPt->Integral()<<endl;
  cout<<"Expected Znunu QCD rate normalized to 1 fb-1: "<<znn_bosonPt->Integral()<<endl;
  
  //////////////////////////

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  /////////////
  plotDistribution(canvas,qqH_bosonPt,znn_bosonPt,znnewk_bosonPt,"Mediator p_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_bosonEta,znn_bosonEta,znnewk_bosonEta,"Mediator #eta",outputPlot);
  plotDistribution(canvas,qqH_bosonPhi,znn_bosonPhi,znnewk_bosonPhi,"Mediator #phi",outputPlot);
  plotDistribution(canvas,qqH_bosonPz,znn_bosonPz,znnewk_bosonPz,"Mediator p_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_bosonPzOPt,znn_bosonPzOPt,znnewk_bosonPzOPt,"Mediator p_{z}/p_{T}",outputPlot);
  plotDistribution(canvas,qqH_jetpt1,znn_jetpt1,znnewk_jetpt1,"p^{j1}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpt2,znn_jetpt2,znnewk_jetpt2,"p^{j2}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz1,znn_jetpz1,znnewk_jetpz1,"p^{j1}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz2,znn_jetpz2,znnewk_jetpz2,"p^{j2}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz1Opt,znn_jetpz1Opt,znnewk_jetpz1Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_jetpz2Opt,znn_jetpz2Opt,znnewk_jetpz2Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_jeteta1,znn_jeteta1,znnewk_jeteta1,"#eta_{j1}",outputPlot);
  plotDistribution(canvas,qqH_jeteta2,znn_jeteta2,znnewk_jeteta2,"#eta_{j2}",outputPlot);
  plotDistribution(canvas,qqH_jetphi1,znn_jetphi1,znnewk_jetphi1,"#phi_{j1}",outputPlot);
  plotDistribution(canvas,qqH_jetphi2,znn_jetphi2,znnewk_jetphi2,"#phi_{j2}",outputPlot);
  plotDistribution(canvas,qqH_jetDeltaPhi,znn_jetDeltaPhi,znnewk_jetDeltaPhi,"#Delta#phi_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_jetDeltaEta,znn_jetDeltaEta,znnewk_jetDeltaEta,"#Delta#eta_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_Mjj,znn_Mjj,znnewk_Mjj,"M_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_jetmetDphi,znn_jetmetDphi,znnewk_jetmetDphi,"#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_jetmetDphi4,znn_jetmetDphi4,znnewk_jetmetDphi4,"#Delta#phi_{4-jet,met}",outputPlot);
  plotDistribution(canvas,qqH_jetmediatorDeta,znn_jetmediatorDeta,znnewk_jetmediatorDeta,"#Delta#eta_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_jetmediatorDphi,znn_jetmediatorDphi,znnewk_jetmediatorDphi,"#Delta#phi_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_MetMT,znn_MetMT,znnewk_MetMT,"m_{T}(jj,met) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_MT,znn_MT,znnewk_MT,"m_{T}(jj) [GeV]",outputPlot);
  /////////////
  plotDistribution(canvas,qqH_mjj_detajj,znn_mjj_detajj,znnewk_mjj_detajj,"m_{jj} [GeV]","#Delta#eta_{jj}",outputPlot);
  plotDistribution(canvas,qqH_mjj_jetmetDphi,znn_mjj_jetmetDphi,znnewk_mjj_jetmetDphi,"m_{jj} [GeV]","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_mjj_jetDphi,znn_mjj_jetDphi,znnewk_mjj_jetDphi,"m_{jj} [GeV]","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_mjj_MT,znn_mjj_MT,znnewk_mjj_MT,"m_{jj} [GeV]","m_{T}(jj) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_detajj_jetmetDphi,znn_detajj_jetmetDphi,znnewk_detajj_jetmetDphi,"#Delta#eta_{jj}","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_detajj_jetDphi,znn_detajj_jetDphi,znnewk_detajj_jetDphi,"#Delta#eta_{jj}","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_detajj_MT,znn_detajj_MT,znnewk_detajj_MT,"#Delta#eta_{jj}","m_{T}(jj) [GeV]",outputPlot);
}
