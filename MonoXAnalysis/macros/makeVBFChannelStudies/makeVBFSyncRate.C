#include "../CMS_lumi.h"

void fillHisto(TH1* histo, double val, double weight){ // embed the overflow                                                                                                                       
 
  if(val < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    histo->Fill(val,weight);
  else
    histo->Fill(histo->GetXaxis()->GetBinCenter(histo->GetNbinsX()),weight);

}

void makeSRSelections(TH1F* histo, TChain* chain, vector<TH1*> khists, bool isKfactor, double lumi){

  // pileup re-weight                                                                                                                                                                                 
  TFile* pileupFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_12p9fb.root");
  TH1*   puhist = (TH1*) pileupFile->Get("puhist");
  // trigger MET                                                                                                                                                                                      
  TFile* triggerfile_MET   =  TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency_12p9.root");
  TEfficiency*  triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();

  TTreeReader myReader (chain);
  TTreeReaderValue<double> xsec         (myReader,"xsec");
  TTreeReaderValue<double> wgt          (myReader,"wgt");
  TTreeReaderValue<double> wgtsum       (myReader,"wgtsum");
  TTreeReaderValue<unsigned int> njetsinc   (myReader,"njetsinc");
  TTreeReaderValue<unsigned int> nphotons   (myReader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus   (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons  (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets  (myReader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nvtx    (myReader,"nvtx");
  TTreeReaderValue<UChar_t> fhbhe   (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso  (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc    (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb    (myReader,"flageebadsc");
  TTreeReaderValue<vector<double> > jetpt   (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<double> > jetm    (myReader,"combinejetm");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"combinejetNHfrac");
  TTreeReaderValue<double > mediatorMass    (myReader,"dmmass");
  TTreeReaderValue<double > mediatorPt      (myReader,"dmpt");
  TTreeReaderValue<double > mediatorEta     (myReader,"dmeta");
  TTreeReaderValue<double > mediatorPhi     (myReader,"dmphi");
  TTreeReaderValue<double>  mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double>  mmet        (myReader,"t1mumet");
  TTreeReaderValue<double>  wzpt        (myReader,"wzpt");
  TTreeReaderValue<double>  wzeta       (myReader,"wzeta");
  TTreeReaderValue<double>  wzphi       (myReader,"wzphi");
  TTreeReaderValue<double>  wzmass      (myReader,"wzmass");
  TTreeReaderValue<double>  genmet      (myReader,"genmet");
  TTreeReaderValue<double>  genmetphi   (myReader,"genmetphi");
  TTreeReaderValue<double>  jmmdphi4    (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double>  jmmdphi     (myReader,"incjetmumetdphimin");  
  TTreeReaderValue<UChar_t> hltm90      (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120     (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120   (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90    (myReader,"hltmetwithmu90");

  while(myReader.Next()){

    if (*hltm90      == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets
    if (*njetsinc   < 2) continue;
    // b-veto
    if (*nbjets     > 0) continue;
    // vetos
    if (*nmuons     > 0) continue;
    if (*nelectrons > 0) continue;
    if (*ntaus      > 0) continue;
    if (*nphotons   > 0) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 80) continue; 
    if (jetpt->at(1) < 70) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut
    if (*mmet     < 150) continue;
    if (*jmmdphi4 < 2.3) continue;
    if (jeteta->at(0)*jeteta->at(1) > 0) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 3.6) continue;
    if ((jet1+jet2).M() < 1100) continue;
    
    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*mediatorPt,*mediatorEta,*mediatorPhi,*mediatorMass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);

    // pileup weight                                                                                                                                                                                  
    double puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    // met trigger                                                                                                                                                                                    
    double trgwgt = triggermet_graph->Eval(min(*mmet,triggermet_graph->GetXaxis()->GetXmax()));
    // kfactor                                                                                                                                                                                        
    double kwgt = 1.0;
    double genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i] and isKfactor) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    double weight = lumi*(*xsec)*(*wgt)*puwgt*trgwgt*kwgt/(*wgtsum); 
    fillHisto(histo,*mediatorPt,weight);
  }
}

void makeVBFSyncRate(double lumi){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TChain* vbftree = new TChain("tree/tree");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*125*root");

  TChain* gghtree = new TChain("tree/tree");
  gghtree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/sig_*GluGlu*_M*125*root");

  TChain* znntree = new TChain("tree/tree");
  znntree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-*root");

  TChain* wjettree = new TChain("tree/tree");
  wjettree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJets/sigfilter/*root");

  TChain* znnewktree = new TChain("tree/tree");
  znnewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJetsToNuNuEWK/sigfilter/sig_EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8.root");

  TChain* wjetewktree = new TChain("tree/tree");
  wjetewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJetsEWK/sigfilter/*root");

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

  TH1F* qqH_bosonPt = new TH1F("qqH_bosonPt","qqH_bosonPt",30,0,1250);
  TH1F* ggH_bosonPt = new TH1F("ggH_bosonPt","ggH_bosonPt",30,0,1250);
  TH1F* znn_bosonPt = new TH1F("znn_bosonPt","znn_bosonPt",30,0,1250);
  TH1F* wjet_bosonPt = new TH1F("wjet_bosonPt","wjet_bosonPt",30,0,1250);
  TH1F* znnewk_bosonPt = new TH1F("znnewk_bosonPt","znnewk_bosonPt",30,0,1250);
  TH1F* wjetewk_bosonPt = new TH1F("wjetewk_bosonPt","wjetewk_bosonPt",30,0,1250);

  makeSRSelections(qqH_bosonPt,vbftree,khists,false,lumi);
  double error = 0;
  double integral = qqH_bosonPt->IntegralAndError(0,qqH_bosonPt->GetNbinsX()+1,error);
  cout<<" VBF rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  makeSRSelections(ggH_bosonPt,gghtree,khists,false,lumi);
  integral = ggH_bosonPt->IntegralAndError(0,ggH_bosonPt->GetNbinsX()+1,error);
  cout<<" ggH rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  makeSRSelections(znn_bosonPt,znntree,khists,true,lumi);
  integral = znn_bosonPt->IntegralAndError(0,znn_bosonPt->GetNbinsX()+1,error);
  cout<<" Znn rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  makeSRSelections(wjet_bosonPt,wjettree,khists,true,lumi);
  integral = wjet_bosonPt->IntegralAndError(0,wjet_bosonPt->GetNbinsX()+1,error);
  cout<<" WJet rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  makeSRSelections(znnewk_bosonPt,znnewktree,khists,false,lumi);
  integral = znnewk_bosonPt->IntegralAndError(0,znnewk_bosonPt->GetNbinsX()+1,error);
  cout<<" Znn EWK rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  makeSRSelections(wjetewk_bosonPt,wjetewktree,khists,false,lumi);
  integral = wjetewk_bosonPt->IntegralAndError(0,wjetewk_bosonPt->GetNbinsX()+1,error);
  cout<<" WJet EWK rate in 1fb-1 : "<<integral<<" error "<<error<<endl;
  
}
