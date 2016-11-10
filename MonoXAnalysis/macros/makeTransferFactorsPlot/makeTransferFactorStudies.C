#include "../CMS_lumi.h"

static float bosonPt        = 200;
static float mjj            = 450;
static float detajj         = 2.5;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;

vector<float> bin_bosonPt {200.,250,350,800};
float lumi_ = 12.9;

enum class Sample {sig, zmm, zee, wmn, wen};

void makePlot(TH1F* tf_qcd, TH1F* tf_ewk, string xAxisName, string outputPlotDIR, string postfix){

  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,625);

  tf_qcd->GetXaxis()->SetTitle(xAxisName.c_str());
  tf_qcd->GetYaxis()->SetTitle("Transfer factor");
  tf_qcd->GetYaxis()->SetTitleOffset(1.1);
  tf_qcd->GetXaxis()->SetTitleOffset(1.1);
  tf_qcd->GetYaxis()->SetTitleSize(tf_qcd->GetYaxis()->GetTitleSize()*1.2);
  tf_qcd->GetXaxis()->SetTitleSize(tf_qcd->GetXaxis()->GetTitleSize()*1.2);
  tf_qcd->GetYaxis()->CenterTitle();

  tf_qcd->SetLineColor(kRed);
  tf_qcd->SetLineWidth(2);
  if(tf_ewk != NULL){
    tf_ewk->SetLineColor(kBlue);
    tf_ewk->SetLineWidth(2);
  }
  TH1F* temp = (TH1F*) tf_qcd->Clone(Form("%s_temp",tf_qcd->GetName()));
  temp->SetFillColor(kRed);
  temp->SetFillStyle(3001);
  temp->SetMarkerSize(0);

  TH1F* temp2 = NULL;
  if(tf_ewk != NULL){
    temp2 = (TH1F*) tf_ewk->Clone(Form("%s_temp",tf_ewk->GetName()));
    temp2->SetFillColor(kBlue);
    temp2->SetFillStyle(3004);
    temp2->SetMarkerSize(0);
  }
  if(tf_ewk != NULL){
    tf_qcd->GetYaxis()->SetRangeUser(min(tf_qcd->GetMinimum(),tf_ewk->GetMinimum())*0.7,max(tf_qcd->GetMaximum(),tf_ewk->GetMaximum())*1.3);
  }
  else
    tf_qcd->GetYaxis()->SetRangeUser(tf_qcd->GetMinimum()*0.5, tf_qcd->GetMaximum()*2);

  if(tf_qcd->Integral() != 0){
    tf_qcd->Draw("hist");
    temp->Draw("E2 same");
    tf_qcd->Draw("hist same");  
    if(tf_ewk != NULL and tf_ewk->Integral() != 0){
      temp2->Draw("E2 same");
      tf_qcd->Draw("hist same");
      tf_ewk->Draw("hist same");
    }
  }
  CMS_lumi(canvas,Form("%.1f",lumi_),true);

  canvas->SaveAs((outputPlotDIR+"/transfer"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputPlotDIR+"/transfer"+postfix+".pdf").c_str(),"pdf");  

}

void applySelection(TTree* tree, TH1F* histo_mu, TH1F* histo_el, TH1F* histo_tau, vector<TH1F*> kfact, const Sample & sample, const bool & applyVBFselection){

  TTreeReader reader(tree);
  TTreeReaderValue<double>  mu1pt     (reader,"mu1pt");
  TTreeReaderValue<double>  mu1eta    (reader,"mu1eta");
  TTreeReaderValue<double>  mu1phi    (reader,"mu1phi");
  TTreeReaderValue<int>     mu1id     (reader,"mu1id");
  TTreeReaderValue<int>     mu1pid    (reader,"mu1pid");
  TTreeReaderValue<double>  mu2pt     (reader,"mu2pt");
  TTreeReaderValue<double>  mu2eta    (reader,"mu2eta");
  TTreeReaderValue<double>  mu2phi    (reader,"mu2phi");
  TTreeReaderValue<int>     mu2id     (reader,"mu2id");
  TTreeReaderValue<int>     mu2pid    (reader,"mu2pid");
  TTreeReaderValue<double>  el1pt     (reader,"el1pt");
  TTreeReaderValue<double>  el1eta    (reader,"el1eta");
  TTreeReaderValue<double>  el1phi    (reader,"el1phi");
  TTreeReaderValue<int>     el1id     (reader,"el1id");
  TTreeReaderValue<int>     el1pid    (reader,"el1pid");
  TTreeReaderValue<double>  el2pt     (reader,"el2pt");
  TTreeReaderValue<double>  el2eta    (reader,"el2eta");
  TTreeReaderValue<double>  el2phi    (reader,"el2phi");
  TTreeReaderValue<int>     el2id     (reader,"el2id");
  TTreeReaderValue<int>     el2pid    (reader,"el2pid");
  TTreeReaderValue<unsigned int> ntaus       (reader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
  TTreeReaderValue<vector<double> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<double> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<double> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<double> > jetchfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > jetnhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<double> met         (reader,"t1pfmet");
  TTreeReaderValue<double> metphi      (reader,"t1pfmetphi");
  TTreeReaderValue<double> mmet        (reader,"t1mumet");
  TTreeReaderValue<double> mmetphi     (reader,"t1mumetphi");
  TTreeReaderValue<double> emet        (reader,"t1elmet");
  TTreeReaderValue<double> emetphi     (reader,"t1elmetphi");
  TTreeReaderValue<double> metpf       (reader,"pfmet");
  TTreeReaderValue<double> metcalo     (reader,"calomet");
  TTreeReaderValue<double> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi (reader,"incjetelmetdphimin4");

  TTreeReaderValue<double> zmass     (reader,"zmass");
  TTreeReaderValue<double> zeemass   (reader,"zeemass");

  TTreeReaderValue<double> wmt    (reader,"wmt");
  TTreeReaderValue<double> wemt   (reader,"wemt");
  TTreeReaderValue<double> wgt    (reader,"wgt");
  TTreeReaderValue<double> xsec   (reader,"xsec");
  TTreeReaderValue<double> wgtsum (reader,"wgtsum");

  TTreeReaderValue<double> wzpt    (reader,"wzpt");
  TTreeReaderValue<int> l1id    (reader,"l1id");
  TTreeReaderValue<int> l2id    (reader,"l2id");

 
  
  //////////////////////////
  while(reader.Next()){

    if(*nphotons != 0) continue;
    if(*nbjets   != 0) continue;
    if(*nincjets  < 2) continue;
    if(*ntaus    != 0) continue;

    if(jetpt->at(0) < leadingJetVBF)  continue;
    if(jetpt->at(1) < trailingJetVBF) continue;
    if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;

    if(fabs(jeteta->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;
    if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;

    ///////////////
    if(sample == Sample::sig){
      if(*nmuons != 0)       continue; 
      if(*nelectrons != 0)   continue;
      if(*met < bosonPt)     continue;
      if(fabs(*jmmdphi) < 1) continue;
    }
    else if(sample == Sample::zmm){
      if (*nelectrons != 0)    continue;
      if (*nmuons != 2)        continue;
      if (fabs(*mu1eta) > 2.4) continue;
      if (fabs(*mu2eta) > 2.4) continue;
      if (*mu1pid == *mu2pid)  continue;
      if (not ((*mu1pt > 20 and *mu1id == 1) or (*mu2pt > 20 and *mu2id == 1))) continue;
      if (*zmass < 60 or *zmass > 120) continue;
      if (fabs(*jmmdphi) < 1)  continue;
      if (*mmet < bosonPt)     continue;
    }
    else if(sample == Sample::wmn){
      if(*nelectrons != 0) continue;
      if(*nmuons != 1)     continue;
      if(*mu1pt < 20)      continue;
      if(fabs(*mu1eta) > 2.4) continue;
      if(*mu1id != 1) continue;
      if(*wmt > 160) continue;
      if(fabs(*jmmdphi) < 1) continue;
      if(*mmet < bosonPt) continue;
    }

    else if(sample == Sample::wen){
      if(*nmuons != 0) continue;
      if(*nelectrons != 1 ) continue;      
      if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
      if(*el1id != 1) continue;
      if(*met < 50) continue;
      if(*wemt > 160) continue;
      if(fabs(*jemdphi) < 1) continue;
      if(*emet < bosonPt) continue;
    }
    else if(sample == Sample::zee){
      if(*nmuons != 0) continue;
      if (*nelectrons != 2 ) continue;
      if (fabs(*el1eta) > 2.5) continue;
      if (fabs(*el2eta) > 2.5) continue;
      if (*el1pid == *el2pid) continue;
      if (not ((*el1pt > 40 and *el1id == 1) or (*el2pt > 40 and *el2id == 1))) continue;
      if (*zeemass < 60 or *zeemass > 120) continue;
      if(fabs(*jemdphi) < 1) continue;
      if(*emet < bosonPt) continue;  
    }

    ////
    if(applyVBFselection){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      if(leadingJet.Eta()*subleadingJet.Eta() > 0) continue;
      if(fabs(leadingJet.Eta()-subleadingJet.Eta()) < detajj) continue;
      if((leadingJet+subleadingJet).M() < mjj) continue;
    }
 

    // k-factord
    Double_t kwgt = 1.0;
    double genpt = *wzpt;
    for (size_t i = 0; i < kfact.size(); i++) {
      if (kfact[i]) {
	if(genpt <= kfact[i]->GetXaxis()->GetBinLowEdge(1)) genpt = kfact[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= kfact[i]->GetXaxis()->GetBinLowEdge(kfact[i]->GetNbinsX()+1)) genpt = kfact[i]->GetXaxis()->GetBinLowEdge(kfact[i]->GetNbinsX()+1)-1;
	kwgt *= kfact[i]->GetBinContent(kfact[i]->FindBin(genpt));
      }
    }
    
    if(sample == Sample::zmm){
      if(fabs(*l1id) == 13 and fabs(*l2id) == 13)
	histo_mu->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 11 and fabs(*l2id) == 11)
	histo_el->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 15 and fabs(*l2id) == 15)
	histo_tau->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else
	cerr<<"Problem z->mm events : gen lep 1 pdg "<<*l1id<<" gen lep 2 pdg "<<*l2id<<endl;
    }
    else if(sample == Sample::zee){
      if(fabs(*l1id) == 11 and fabs(*l2id) == 11)
	histo_el->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 13 and fabs(*l2id) == 13)
	histo_mu->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 15 and fabs(*l2id) == 15)
	histo_tau->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else
	cerr<<"Problem z->ee events : gen lep 1 pdg "<<*l1id<<" gen lep 2 pdg "<<*l2id<<endl;
    }
    else if(sample == Sample::wen){
      if(fabs(*l1id) == 11 or fabs(*l2id) == 11)
	histo_el->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 13 or fabs(*l2id) == 13)
	histo_mu->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 15 or fabs(*l2id) == 15)
	histo_tau->Fill(*emet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else
	cerr<<"Problem w->en events : gen lep 1 pdg "<<*l1id<<" gen lep 2 pdg "<<*l2id<<endl;
    }
    else if(sample == Sample::wmn){    
      if(fabs(*l1id) == 11 or fabs(*l2id) == 11)
	histo_el->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 13 or fabs(*l2id) == 13)
	histo_mu->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 15 or fabs(*l2id) == 15)
	histo_tau->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else
	cerr<<"Problem w->mn events : gen lep 1 pdg "<<*l1id<<" gen lep 2 pdg "<<*l2id<<endl;
    }
    else if(sample == Sample::sig){
      if(fabs(*l1id) == 11 or fabs(*l2id) == 11 or (fabs(*l1id) == 12 and fabs(*l2id) == 12))
	histo_el->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 13 or fabs(*l2id) == 13 or (fabs(*l1id) == 14 and fabs(*l2id) == 14))
	histo_mu->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else if(fabs(*l1id) == 15 or fabs(*l2id) == 15 or (fabs(*l1id) == 16 and fabs(*l2id) == 16))
	histo_tau->Fill(*mmet,lumi_*(*wgt)*(*xsec)*kwgt/(*wgtsum));
      else
	cerr<<"Problem sig events : gen lep 1 pdg "<<*l1id<<" gen lep 2 pdg "<<*l2id<<endl;
    }
  }
}


void makeTransferFactorStudies(Sample sample, string outputPlotDIR, bool applyVBFselection = true, bool applyKfactor = true){

  string postfix;
  if(sample == Sample::wmn)
    postfix = "_wm";
  else if(sample == Sample::wen)
    postfix = "_wen";
  else if(sample == Sample::zmm)
    postfix = "_zmm";
  else if(sample == Sample::zee)
    postfix = "_zee";


  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  
  system(("mkdir -p "+outputPlotDIR).c_str());

  TChain* tree_vjet_qcd_den = new TChain("tree/tree"); 
  TChain* tree_vjet_ewk_den = new TChain("tree/tree");
  TChain* tree_vjet_qcd_num = new TChain("tree/tree");
  TChain* tree_vjet_ewk_num = new TChain("tree/tree");

  if(sample == Sample::zmm){
    tree_vjet_qcd_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/DYJets/zmmfilter/*root");
    tree_vjet_qcd_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJets/sigfilter/*root");
    tree_vjet_ewk_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJetsToLLEWK/zmmfilter/*root");
    tree_vjet_ewk_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJetsToNuNuEWK/sigfilter/*root");
  }
  else if(sample == Sample::zee){
    tree_vjet_qcd_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/DYJets/zeefilter/*root");
    tree_vjet_qcd_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJets/sigfilter/*root");
    tree_vjet_ewk_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJetsToLLEWK/zeefilter/*root");
    tree_vjet_ewk_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/ZJetsToNuNuEWK/sigfilter/*root");
  }
  else if(sample == Sample::wmn){
    tree_vjet_qcd_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJets/wmnfilter/*root");
    tree_vjet_qcd_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJets/sigfilter/*root");
    tree_vjet_ewk_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJetsEWK/wmnfilter/*root");
    tree_vjet_ewk_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJetsEWK/sigfilter/*root");
  }
  else if(sample == Sample::wen){
    tree_vjet_qcd_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJets/wenfilter/*root");
    tree_vjet_qcd_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJets/sigfilter/*root");
    tree_vjet_ewk_den->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJetsEWK/wenfilter/*root");
    tree_vjet_ewk_num->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/WJetsEWK/sigfilter/*root");
 }

  ////// k-factors                                                                                                                                                                      
  vector<TH1F*> khists;
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
  }

  // declare histograms
  cout<<"Loop vjet qcd denominator "<<endl;
  TH1F* vjet_qcd_den_muon = new TH1F("vjet_qcd_den_muon","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_qcd_den_ele  = new TH1F("vjet_qcd_den_ele","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_qcd_den_tau  = new TH1F("vjet_qcd_den_tau","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  vjet_qcd_den_muon->Sumw2();
  vjet_qcd_den_ele->Sumw2();
  vjet_qcd_den_tau->Sumw2();
  applySelection(tree_vjet_qcd_den,vjet_qcd_den_muon,vjet_qcd_den_ele,vjet_qcd_den_tau,khists,sample,applyVBFselection);
  cout<<"Rate: vjet_qcd_den_muon = "<<vjet_qcd_den_muon->Integral()<<endl;
  cout<<"Rate: vjet_qcd_den_ele  = "<<vjet_qcd_den_ele->Integral()<<endl;
  cout<<"Rate: vjet_qcd_den_tau  = "<<vjet_qcd_den_tau->Integral()<<endl;

  cout<<"Loop vjet ewk denominator "<<endl;
  TH1F* vjet_ewk_den_muon = new TH1F("vjet_ewk_den_muon","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_ewk_den_ele  = new TH1F("vjet_ewk_den_ele","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_ewk_den_tau  = new TH1F("vjet_ewk_den_tau","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  vjet_ewk_den_muon->Sumw2();
  vjet_ewk_den_ele->Sumw2();
  vjet_ewk_den_tau->Sumw2();
  applySelection(tree_vjet_ewk_den,vjet_ewk_den_muon,vjet_ewk_den_ele,vjet_ewk_den_tau,khists,sample,applyVBFselection);
  cout<<"Rate: vjet_ewk_den_muon = "<<vjet_ewk_den_muon->Integral()<<endl;
  cout<<"Rate: vjet_ewk_den_ele  = "<<vjet_ewk_den_ele->Integral()<<endl;
  cout<<"Rate: vjet_ewk_den_tau  = "<<vjet_ewk_den_tau->Integral()<<endl;

  cout<<"Loop vjet qcd numerator "<<endl;
  TH1F* vjet_qcd_num_muon = new TH1F("vjet_qcd_num_muon","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_qcd_num_ele  = new TH1F("vjet_qcd_num_ele","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_qcd_num_tau  = new TH1F("vjet_qcd_num_tau","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  vjet_qcd_num_muon->Sumw2();
  vjet_qcd_num_ele->Sumw2();
  vjet_qcd_num_tau->Sumw2();
  applySelection(tree_vjet_qcd_num,vjet_qcd_num_muon,vjet_qcd_num_ele,vjet_qcd_num_tau,khists,Sample::sig,applyVBFselection);
  cout<<"Rate: vjet_qcd_num_muon = "<<vjet_qcd_num_muon->Integral()<<endl;
  cout<<"Rate: vjet_qcd_num_ele  = "<<vjet_qcd_num_ele->Integral()<<endl;
  cout<<"Rate: vjet_qcd_num_tau  = "<<vjet_qcd_num_tau->Integral()<<endl;

  TH1F* vjet_qcd_num_total = (TH1F*) vjet_qcd_num_muon->Clone("vjet_qcd_num_total");
  if(sample == Sample::zmm or sample == Sample::zee){
    vjet_qcd_num_total->Add(vjet_qcd_num_ele);
    vjet_qcd_num_total->Add(vjet_qcd_num_tau);

  }

  cout<<"Loop vjet ewk numerator "<<endl;
  TH1F* vjet_ewk_num_muon = new TH1F("vjet_ewk_num_muon","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_ewk_num_ele  = new TH1F("vjet_ewk_num_ele","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  TH1F* vjet_ewk_num_tau  = new TH1F("vjet_ewk_num_tau","",bin_bosonPt.size()-1,&bin_bosonPt[0]);
  vjet_ewk_num_muon->Sumw2();
  vjet_ewk_num_ele->Sumw2();
  vjet_ewk_num_tau->Sumw2();
  applySelection(tree_vjet_ewk_num,vjet_ewk_num_muon,vjet_ewk_num_ele,vjet_ewk_num_tau,khists,Sample::sig,applyVBFselection);
  cout<<"Rate: vjet_ewk_num_muon = "<<vjet_ewk_num_muon->Integral()<<endl;
  cout<<"Rate: vjet_ewk_num_ele  = "<<vjet_ewk_num_ele->Integral()<<endl;
  cout<<"Rate: vjet_ewk_num_tau  = "<<vjet_ewk_num_tau->Integral()<<endl;

  TH1F* vjet_ewk_num_total = (TH1F*) vjet_ewk_num_muon->Clone("vjet_ewk_num_total");
  if(sample == Sample::zmm or sample == Sample::zee){
    vjet_ewk_num_total->Add(vjet_ewk_num_ele);
    vjet_ewk_num_total->Add(vjet_ewk_num_tau);

  }
  
  if(sample != Sample::zmm and sample != Sample::zee){
    TH1F* tf_qcd_muon = (TH1F*) vjet_qcd_num_muon->Clone("tf_qcd_muon");
    tf_qcd_muon->Divide(vjet_qcd_den_muon);
    TH1F* tf_ewk_muon = (TH1F*) vjet_ewk_num_muon->Clone("tf_ewk_muon");
    tf_ewk_muon->Divide(vjet_ewk_den_muon);
    
    makePlot(tf_qcd_muon,tf_ewk_muon,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_muon");
    
    TH1F* tf_qcd_ele = (TH1F*) vjet_qcd_num_ele->Clone("tf_qcd_ele");
    TH1F* tf_ewk_ele = (TH1F*) vjet_ewk_num_ele->Clone("tf_ewk_ele");
    tf_qcd_ele->Divide(vjet_qcd_den_ele);
    tf_ewk_ele->Divide(vjet_ewk_den_ele);
    
    makePlot(tf_qcd_ele,tf_ewk_ele,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_ele");
    
    TH1F* tf_qcd_tau = (TH1F*) vjet_qcd_num_tau->Clone("tf_qcd_tau");
    TH1F* tf_ewk_tau = (TH1F*) vjet_ewk_num_tau->Clone("tf_ewk_tau");
    tf_qcd_tau->Divide(vjet_qcd_den_tau);
    tf_ewk_tau->Divide(vjet_ewk_den_tau);
    
    makePlot(tf_qcd_tau,tf_ewk_tau,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_tau");
  }
  else{
    TH1F* tf_qcd_muon = (TH1F*) vjet_qcd_num_total->Clone("tf_qcd_muon");
    tf_qcd_muon->Divide(vjet_qcd_den_muon);
    TH1F* tf_ewk_muon = (TH1F*) vjet_ewk_num_total->Clone("tf_ewk_muon");
    tf_ewk_muon->Divide(vjet_ewk_den_muon);
    
    makePlot(tf_qcd_muon,tf_ewk_muon,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_muon");
    
    TH1F* tf_qcd_ele = (TH1F*) vjet_qcd_num_total->Clone("tf_qcd_ele");
    TH1F* tf_ewk_ele = (TH1F*) vjet_ewk_num_total->Clone("tf_ewk_ele");
    tf_qcd_ele->Divide(vjet_qcd_den_ele);
    tf_ewk_ele->Divide(vjet_ewk_den_ele);
    
    makePlot(tf_qcd_ele,tf_ewk_ele,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_ele");
    
    TH1F* tf_qcd_tau = (TH1F*) vjet_qcd_num_total->Clone("tf_qcd_tau");
    TH1F* tf_ewk_tau = (TH1F*) vjet_ewk_num_total->Clone("tf_ewk_tau");
    tf_qcd_tau->Divide(vjet_qcd_den_tau);
    tf_ewk_tau->Divide(vjet_ewk_den_tau);
    
    makePlot(tf_qcd_tau,tf_ewk_tau,"boson p_{T} [GeV]",outputPlotDIR,postfix+"_tau");
    
  }

}
