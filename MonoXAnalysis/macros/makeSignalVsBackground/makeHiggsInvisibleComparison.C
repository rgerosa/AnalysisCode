#include "../CMS_lumi.h"

// all the observables asaf of different selections applied: baseline is object veto, then leading jet requirement, then jet met dphi, then recoil cut, then VBF selections
vector<TH1F*> ggH_leadjetpt;
vector<TH1F*> ggH_traijetpt;
vector<TH1F*> ggH_leadjeteta;
vector<TH1F*> ggH_traijeteta;
vector<TH1F*> ggH_met;
vector<TH1F*> ggH_leadQGL;
vector<TH1F*> ggH_traiQGL;
vector<TH1F*> ggH_njet;
vector<TH1F*> ggH_ncjet;
vector<TH1F*> ggH_mT;
vector<TH1F*> ggH_HT;
vector<TH1F*> ggH_detajj;
vector<TH1F*> ggH_dphijj;
vector<TH1F*> ggH_mjj;
vector<TH1F*> ggH_etaj1etaj2;
vector<TH1F*> ggH_cjv;
vector<TH1F*> ggH_dphijmet;

vector<TH2F*> ggH_jetpt1jetpt2;
vector<TH2F*> ggH_jeteta1jeteta2;
vector<TH2F*> ggH_metnjet;
vector<TH2F*> ggH_metQGL;
vector<TH2F*> ggH_metmjj;
vector<TH2F*> ggH_metdetajj;
vector<TH2F*> ggH_detajjmjj;


void makeHiggsInvisibleComparison(string ggHDIR, string qqHDIR, string vjetsQCD, string vjetsEWK, string masspoint, string outputDIR){
  
  TChain* ggH     = new TChain("tree/tree");
  TChain* qqH     = new TChain("tree/tree");
  TChain* vjetQCD = new TChain("tree/tree");
  TChain* vjetEWK = new TChain("tree/tree");
  
  ggH->Add((ggHDIR+"/*150*root").c_str());
  qqH->Add((qqHDIR+"/*150*root").c_str());
  
  vjetQCD->Add((vjetsQCD+"/*root").c_str());
  vjetEWK->Add((vjetsEWK+"/*root").c_str());

  

}
