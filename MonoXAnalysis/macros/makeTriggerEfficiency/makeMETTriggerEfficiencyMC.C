#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// recoil binning for monojet                                                                                                                                                                    
vector <float> bins_monojet_recoil = {50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 265., 280., 300, 320., 340., 360., 380., 400., 430., 460., 490., 520, 550., 580., 610., 650., 700., 740., 800., 900., 1000.,1250};

enum class Sample {sig,wmn,wen,zmm,zee};

static bool makeSelectionGen = false;
static bool applyJetSelections = true;
static bool askForTriggerDenominator = true;
static bool computeAsFunctionOfMET = false;

void calculateSumWeight(vector<TTree*> gentree, vector<double> & wgtsum){

  for(auto tree: gentree){
    wgtsum.push_back(0);
  }

  int itree = 0;
  for(auto tree: gentree){
    TTreeReader reader(tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    //////////////////                                                                                                                                                                           
    long int nTotal = tree->GetEntries();
    long int nEvents = 0;
    long int nPart = 100000;
    cout<<"Looping on itree "<<itree<<" of "<<gentree.size()<<" Total number of events: "<<nTotal<<endl;
    while(reader.Next()){
      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;
      wgtsum.at(itree) += *wgt;
    }
    cout<<endl;
    cout<<"Sum of weigths for this tree "<<wgtsum.at(itree)<<endl;
    itree++;
  }
}

////////////////
TH1F* metCalo_SR_pas = new TH1F("metCalo_SR_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metCalo_WM_pas = new TH1F("metCalo_WM_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metCalo_ZM_pas = new TH1F("metCalo_ZM_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_SR_pas = new TH1F("metPF_SR_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_WM_pas = new TH1F("metPF_WM_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_ZM_pas = new TH1F("metPF_ZM_pas","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);

TH1F* metCalo_SR_fail = new TH1F("metCalo_SR_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metCalo_WM_fail = new TH1F("metCalo_WM_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metCalo_ZM_fail = new TH1F("metCalo_ZM_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_SR_fail = new TH1F("metPF_SR_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_WM_fail = new TH1F("metPF_WM_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);
TH1F* metPF_ZM_fail = new TH1F("metPF_ZM_fail","",bins_monojet_recoil.size()-1,&bins_monojet_recoil[0]);

//////////
void makeTriggerAnalysis(vector<TTree*> trees, TH1F* hnum, TH1F* hden, const Sample & sample, vector<double> wgtsum, float luminosity){
  
  cout<<"Loop on trees "<<endl;

  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");

  int itree = 0;
  for(auto tree : trees){    

    TTreeReader reader(tree);
    TTreeReaderValue<float> xsec         (reader,"xsec");
    TTreeReaderValue<float> wgt          (reader,"wgt");
    TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
    TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
    TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
    TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
    TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
    TTreeReaderValue<UChar_t> hltmwm100  (reader,"hltmetwithmu100");
    TTreeReaderValue<UChar_t> hltmwm110  (reader,"hltmetwithmu110");
    TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
    TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");
    TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
    TTreeReaderValue<UChar_t> hltjm      (reader,"hltjetmet");
    TTreeReaderValue<UChar_t> hltsingleel (reader,"hltsingleel");
    TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
    TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
    TTreeReaderValue<float>   mu1pt      (reader,"mu1pt");
    TTreeReaderValue<float>   mu1eta     (reader,"mu1eta");
    TTreeReaderValue<float>   mu1phi     (reader,"mu1phi");
    TTreeReaderValue<int>     mu1id      (reader,"mu1id");
    TTreeReaderValue<int>     mu1pid     (reader,"mu1pid");
    TTreeReaderValue<float>   el1pt      (reader,"el1pt");
    TTreeReaderValue<float>   el1eta     (reader,"el1eta");
    TTreeReaderValue<float>   el1phi     (reader,"el1phi");
    TTreeReaderValue<int>     el1id      (reader,"el1id");
    TTreeReaderValue<float>   mu2pt      (reader,"mu2pt");
    TTreeReaderValue<float>   mu2eta     (reader,"mu2eta");
    TTreeReaderValue<float>   mu2phi     (reader,"mu2phi");
    TTreeReaderValue<float>   zmass      (reader,"zmass");
    TTreeReaderValue<int>     mu2id      (reader,"mu2id");
    TTreeReaderValue<int>     mu2pid     (reader,"mu2pid");
    TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
    TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
    TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
    TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
    TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
    TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
    TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
    TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
    TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
    TTreeReaderValue<unsigned int> nvtx        (reader,"nvtx");
    TTreeReaderValue<unsigned int> ntaus       (reader,"ntaus");
    TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
    TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
    TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
    TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
    TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
    TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
    TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
    TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
    TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
    TTreeReaderValue<vector<float> > jetchfrac  (reader,"combinejetCHfrac");
    TTreeReaderValue<vector<float> > jetnhfrac  (reader,"combinejetNHfrac");
    TTreeReaderValue<float> met         (reader,"t1pfmet");
    TTreeReaderValue<float> metphi      (reader,"t1pfmetphi");
    TTreeReaderValue<float> mmet        (reader,"t1mumet");
    TTreeReaderValue<float> mmetphi     (reader,"t1mumetphi");
    TTreeReaderValue<float> metpf       (reader,"pfmet");
    TTreeReaderValue<float> metcalo     (reader,"calomet");
    TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
    TTreeReaderValue<int>   l1id (reader,"l1id");
    TTreeReaderValue<int>   l2id (reader,"l2id");

    //////////////////                                                                                                                                                                              
    long int nTotal = tree->GetEntries();
    cout<<"Looping on itree "<<itree<<" of "<<trees.size()<<" Total number of events: "<<nTotal<<endl;
    long int nEvents = 0;    
    long int nPart = 100000;

    while(reader.Next()){
      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;

      if(*nbjets   != 0)    continue;
      if(*ntaus    != 0)    continue;
      if(*nphotons  != 0)   continue;
      if(sample != Sample::wen){
	if(*nelectrons != 0 ) continue;
      }
      else
	if(*nmuons != 0 ) continue;

      if(not *fcsc)  continue;
      if(not *fcsct) continue;
      if(not *feeb)  continue;
      if(not *fetp)  continue;
      if(not *fvtx)  continue;
      if(not *fbadmu) continue;
      if(not *fbadch) continue;
      if(not *fhbhe)  continue;
      if(not *fhbiso) continue;

      // apply calo met cleaning
      if(not makeSelectionGen and fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;

      if(askForTriggerDenominator and sample == Sample::wmn and not *hltsinglemu) continue;
      else if(askForTriggerDenominator and sample == Sample::zmm and (not *hltdoublemu or not *hltsinglemu)) continue;
      else if(askForTriggerDenominator and sample == Sample::wen and not *hltsingleel) continue;
      // apply jet pt selections
      if(applyJetSelections){
	if(*nincjets < 1) continue;
	if(jetpt->size() == 0) continue;
	if(jetpt->at(0) < 100) continue;
	if(fabs(jeteta->at(0)) > 2.5) continue;
	if(jetchfrac->at(0) < 0.1) continue;
	if(jetnhfrac->at(0) > 0.8) continue;

	// apply jet-met dphi--> not for gen level analysis
	if(not makeSelectionGen and *jmmdphi < 0.5) continue;
      }

      // apply standard reco-level cuts
      if(not makeSelectionGen){	
	if(sample == Sample::wmn){
	  if(*mu1pt < 20) continue;
	  if(fabs(*mu1eta) > 2.4) continue;
	  if(*mu1id  !=1) continue;
	  if(*nmuons !=1) continue;
	  float dphi = fabs(*mu1phi-*metphi);
	  if(dphi > TMath::Pi())
	    dphi = 2*TMath::Pi()-dphi;
	  float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
	  if(mtw > 160) continue;
	}
	if(sample == Sample::wen){
	  if(*el1pt < 40) continue;
	  if(fabs(*el1eta) > 2.5) continue;
	  if(*el1id  !=1) continue;
	  if(*nelectrons !=1) continue;
	  float dphi = fabs(*el1phi-*metphi);
	  if(dphi > TMath::Pi())
	    dphi = 2*TMath::Pi()-dphi;
	  float mtw = sqrt(2*(*el1pt)*(*met)*(1-cos(dphi)));
	  if(mtw > 160) continue;
	}
	else if(sample == Sample::zmm){
	  if(*nmuons !=2) continue;
	  if(*mu1pt < 20) continue;
	  if(fabs(*mu1eta) > 2.4) continue;
	  if(fabs(*mu2eta) > 2.4) continue;
	  if(*mu1pid == *mu2pid) continue;
	  if((*mu1pt > 20 and *mu1id != 1) or (*mu2pt > 20 and *mu2id != 1)) continue;
	  if(*zmass < 60 or *zmass > 120) continue;
	}
	else if(sample == Sample::sig){
	  if(*nmuons != 0) continue;
	}
      }
      else{// gen level analysis
	if(sample == Sample::wmn){
	  bool goodEvent = false;
	  if((fabs(*l1id) == 13 and fabs(*l2id) == 14) or (fabs(*l1id) == 14 and fabs(*l2id) == 13)) goodEvent = true;
	  if(not goodEvent) continue;
	}
	else if(sample == Sample::zmm){
	  if(fabs(*l1id) != 13) continue;
	  if(fabs(*l2id) != 13) continue;
	  if(*l1id == *l2id) continue;
	}
	else if(sample == Sample::sig){
	  if(fabs(*l1id) != 12 and fabs(*l1id) != 14 and fabs(*l1id) != 16) continue;
	  if(fabs(*l2id) != 12 and fabs(*l2id) != 14 and fabs(*l2id) != 16) continue;
	}
      }

      // pileup re-weight
      double puwgt = 1;
      if(*nvtx < 60)
	puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
      if(sample != Sample::wen)
	hden->Fill(*mmet,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
      else
	hden->Fill(*met,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));

      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm){
	if(sample != Sample::wen)
	  hnum->Fill(*mmet,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));     
	  else
	  hnum->Fill(*met,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));     

	if(sample == Sample::sig){
	  metCalo_SR_pas->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_SR_pas->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}
	else if(sample == Sample::wmn){
	  metCalo_WM_pas->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_WM_pas->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}
	else if(sample == Sample::zmm){
	  metCalo_ZM_pas->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_ZM_pas->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}	
      }
      else{
	if(sample == Sample::sig){
	  metCalo_SR_fail->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_SR_fail->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}
	else if(sample == Sample::wmn){
	  metCalo_WM_fail->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_WM_fail->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}
	else if(sample == Sample::zmm){
	  metCalo_ZM_fail->Fill(*metcalo,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	  metPF_ZM_fail->Fill(*metpf,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
	}
      }
    }
    cout<<endl;
    itree++;
  }
  
  if(pufile) pufile->Close();

}

void makeMETTriggerEfficiencyMC(string inputDIR, string outputDIR, float luminosity = 35.9, bool applyGenSelection = false){

  makeSelectionGen = applyGenSelection;

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
 
  // input tree                                                                                                                                                                                      
  vector<TFile*> filelist_zvv;
  vector<TFile*> filelist_wjet;
  vector<TFile*> filelist_zll;
  vector<TTree*> tree_zvv;
  vector<TTree*> tree_wjet;
  vector<TTree*> tree_zll;
  vector<TTree*> gentree_zvv;
  vector<TTree*> gentree_wjet;
  vector<TTree*> gentree_zll;

  cout<<"Read list of files for Zvv process "<<endl;
  system(("ls "+inputDIR+" | grep ZJets > list_dir.txt").c_str());
  ifstream file_dir_zvv("list_dir.txt");
  if(file_dir_zvv.is_open()){
    string line;
    while(!file_dir_zvv.eof()){
      getline(file_dir_zvv,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_zvv("list.txt");
      if(file_zvv.is_open()){
	string line2;
	while(!file_zvv.eof()){
	  getline(file_zvv,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open Zvv file with name: "<<line2<<endl;
	  filelist_zvv.push_back(TFile::Open(line2.c_str()));
	  tree_zvv.push_back((TTree*) filelist_zvv.back()->Get("tree/tree"));
	  gentree_zvv.push_back((TTree*) filelist_zvv.back()->Get("gentree/gentree"));
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  cout<<"Read list of files for Wjet process "<<endl;
  system(("ls "+inputDIR+" | grep WJets > list_dir.txt").c_str());
  ifstream file_dir_wjet("list_dir.txt");
  if(file_dir_wjet.is_open()){
    string line;
    while(!file_dir_wjet.eof()){
      getline(file_dir_wjet,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_wjet("list.txt");
      if(file_wjet.is_open()){
	string line2;
	while(!file_wjet.eof()){
	  getline(file_wjet,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open Wjet file with name: "<<line2<<endl;
	  filelist_wjet.push_back(TFile::Open(line2.c_str()));
	  tree_wjet.push_back((TTree*) filelist_wjet.back()->Get("tree/tree"));
	  gentree_wjet.push_back((TTree*) filelist_wjet.back()->Get("gentree/gentree"));
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  cout<<"Read list of files for Zll process "<<endl;
  system(("ls "+inputDIR+" | grep DYJets > list_dir.txt").c_str());
  ifstream file_dir_zll("list_dir.txt");
  if(file_dir_zll.is_open()){
    string line;
    while(!file_dir_zll.eof()){
      getline(file_dir_zll,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_zll("list.txt");
      if(file_zll.is_open()){
	string line2;
	while(!file_zll.eof()){
	  getline(file_zll,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open Zll file with name: "<<line2<<endl;
	  filelist_zll.push_back(TFile::Open(line2.c_str()));
	  tree_zll.push_back((TTree*) filelist_zll.back()->Get("tree/tree"));
	  gentree_zll.push_back((TTree*) filelist_zll.back()->Get("gentree/gentree"));
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  
  // sum of weights
  vector<double> wgtsum_zvv;
  vector<double> wgtsum_wjet;
  vector<double> wgtsum_zll;
  cout<<"######### Calculate Weights for : Zvv "<<endl;
  calculateSumWeight(gentree_zvv,wgtsum_zvv);
  cout<<"######### Calculate Weights for : W+jets "<<endl;
  calculateSumWeight(gentree_wjet,wgtsum_wjet);
  cout<<"######### Calculate Weights for : Zll "<<endl;
  if(not computeAsFunctionOfMET)
    calculateSumWeight(gentree_zll,wgtsum_zll);

  /////// start analysis
  TF1 *fitfunc_monojet_recoil_zvv = new TF1("fitfunc_monojet_recoil_zvv",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_zvv->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_zvv = new TH1F("hnum_monojet_recoil_zvv", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_zvv = new TH1F("hden_monojet_recoil_zvv", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_zvv->Sumw2();
  hden_monojet_recoil_zvv->Sumw2();

  /////// start analysis
  TF1 *fitfunc_monojet_recoil_wjet = new TF1("fitfunc_monojet_recoil_wjet",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_wjet->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_wjet = new TH1F("hnum_monojet_recoil_wjet", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_wjet = new TH1F("hden_monojet_recoil_wjet", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_wjet->Sumw2();
  hden_monojet_recoil_wjet->Sumw2();

  /////// start analysis
  TF1 *fitfunc_monojet_recoil_wmn = new TF1("fitfunc_monojet_recoil_wmn",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_wmn->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_wmn = new TH1F("hnum_monojet_recoil_wmn", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_wmn = new TH1F("hden_monojet_recoil_wmn", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_wmn->Sumw2();
  hden_monojet_recoil_wmn->Sumw2();

  TF1 *fitfunc_monojet_recoil_wen = new TF1("fitfunc_monojet_recoil_wen",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_wen->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_wen = new TH1F("hnum_monojet_recoil_wen", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_wen = new TH1F("hden_monojet_recoil_wen", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_wen->Sumw2();
  hden_monojet_recoil_wen->Sumw2();

  /////// start analysis
  TF1 *fitfunc_monojet_recoil_zmm = new TF1("fitfunc_monojet_recoil_zmm",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_zmm->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_zmm = new TH1F("hnum_monojet_recoil_zmm", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_zmm = new TH1F("hden_monojet_recoil_zmm", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_zmm->Sumw2();
  hden_monojet_recoil_zmm->Sumw2();

  ////
  cout<<"######### Loop on Zvv trees for SR selection "<<endl;
  makeTriggerAnalysis(tree_zvv,hnum_monojet_recoil_zvv,hden_monojet_recoil_zvv,Sample::sig,wgtsum_zvv,luminosity);
  if(not applyGenSelection){
    cout<<"######### Loop on W+jets trees for SR selection "<<endl;
    makeTriggerAnalysis(tree_wjet,hnum_monojet_recoil_wjet,hden_monojet_recoil_wjet,Sample::sig,wgtsum_wjet,luminosity);
  }
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_wjet,hnum_monojet_recoil_wmn,hden_monojet_recoil_wmn,Sample::wmn,wgtsum_wjet,luminosity);
  if(computeAsFunctionOfMET){
    cout<<"######### Loop on W+jets trees for Wen selection "<<endl;
    makeTriggerAnalysis(tree_wjet,hnum_monojet_recoil_wen,hden_monojet_recoil_wen,Sample::wen,wgtsum_wjet,luminosity);
  }
  else{
    cout<<"######### Loop on Zll+jets trees for Zmm selection "<<endl;
    makeTriggerAnalysis(tree_zll,hnum_monojet_recoil_zmm,hden_monojet_recoil_zmm,Sample::zmm,wgtsum_zll,luminosity);
  }

  // Make efficiencies
  TEfficiency* eff_monojet_recoil_zvv = new TEfficiency(*hnum_monojet_recoil_zvv,*hden_monojet_recoil_zvv);
  eff_monojet_recoil_zvv->SetMarkerColor(kBlack);
  eff_monojet_recoil_zvv->SetLineColor(kBlack);
  eff_monojet_recoil_zvv->SetMarkerStyle(20);
  eff_monojet_recoil_zvv->SetMarkerSize(1);
  fitfunc_monojet_recoil_zvv->SetLineColor(kBlack);
  fitfunc_monojet_recoil_zvv->SetLineWidth(2);
  fitfunc_monojet_recoil_zvv->SetLineStyle(7);
  eff_monojet_recoil_zvv->Fit(fitfunc_monojet_recoil_zvv,"RE");

  TEfficiency* eff_monojet_recoil_wjet = new TEfficiency(*hnum_monojet_recoil_wjet,*hden_monojet_recoil_wjet);
  eff_monojet_recoil_wjet->SetMarkerColor(kRed);
  eff_monojet_recoil_wjet->SetLineColor(kRed);
  eff_monojet_recoil_wjet->SetMarkerStyle(20);
  eff_monojet_recoil_wjet->SetMarkerSize(1);
  fitfunc_monojet_recoil_wjet->SetLineColor(kRed);
  fitfunc_monojet_recoil_wjet->SetLineWidth(2);
  fitfunc_monojet_recoil_wjet->SetLineStyle(7);
  eff_monojet_recoil_wjet->Fit(fitfunc_monojet_recoil_wjet,"RE");

  TEfficiency* eff_monojet_recoil_wmn = new TEfficiency(*hnum_monojet_recoil_wmn,*hden_monojet_recoil_wmn);
  eff_monojet_recoil_wmn->SetMarkerColor(kBlue);
  eff_monojet_recoil_wmn->SetLineColor(kBlue);
  eff_monojet_recoil_wmn->SetMarkerStyle(20);
  eff_monojet_recoil_wmn->SetMarkerSize(1);
  fitfunc_monojet_recoil_wmn->SetLineColor(kBlue);
  fitfunc_monojet_recoil_wmn->SetLineWidth(2);
  fitfunc_monojet_recoil_wmn->SetLineStyle(7);
  eff_monojet_recoil_wmn->Fit(fitfunc_monojet_recoil_wmn,"RE");

  TEfficiency* eff_monojet_recoil_wen = new TEfficiency(*hnum_monojet_recoil_wen,*hden_monojet_recoil_wen);
  eff_monojet_recoil_wen->SetMarkerColor(kCyan);
  eff_monojet_recoil_wen->SetLineColor(kCyan);
  eff_monojet_recoil_wen->SetMarkerStyle(20);
  eff_monojet_recoil_wen->SetMarkerSize(1);
  fitfunc_monojet_recoil_wen->SetLineColor(kCyan);
  fitfunc_monojet_recoil_wen->SetLineWidth(2);
  fitfunc_monojet_recoil_wen->SetLineStyle(7);
  eff_monojet_recoil_wen->Fit(fitfunc_monojet_recoil_wen,"RE");

  TEfficiency* eff_monojet_recoil_zmm = new TEfficiency(*hnum_monojet_recoil_zmm,*hden_monojet_recoil_zmm);
  eff_monojet_recoil_zmm->SetMarkerColor(kCyan);
  eff_monojet_recoil_zmm->SetLineColor(kCyan);
  eff_monojet_recoil_zmm->SetMarkerStyle(20);
  eff_monojet_recoil_zmm->SetMarkerSize(1);
  fitfunc_monojet_recoil_zmm->SetLineColor(kCyan);
  fitfunc_monojet_recoil_zmm->SetLineStyle(7);
  fitfunc_monojet_recoil_zmm->SetLineWidth(2);
  eff_monojet_recoil_zmm->Fit(fitfunc_monojet_recoil_zmm,"RE");

  // Plotting final result for MC turn ons
  TH1* frame = canvas->DrawFrame(bins_monojet_recoil.front(),0.,bins_monojet_recoil.back(), 1.1, "");
  frame->GetXaxis()->SetTitle("Recoil [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);

  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);
  
  fitfunc_monojet_recoil_zvv->Draw("Lsame");
  if(not applyGenSelection)
    fitfunc_monojet_recoil_wjet->Draw("Lsame");
  fitfunc_monojet_recoil_wmn->Draw("Lsame");
  if(not computeAsFunctionOfMET)
    fitfunc_monojet_recoil_zmm->Draw("Lsame");
  else
    fitfunc_monojet_recoil_wen->Draw("Lsame");

  eff_monojet_recoil_zvv->Draw("EPsame");
  if(not applyGenSelection)
    eff_monojet_recoil_wjet->Draw("EPsame"); 
  eff_monojet_recoil_wmn->Draw("EPsame");
  if(not computeAsFunctionOfMET)
    eff_monojet_recoil_zmm->Draw("EPsame");
  else
    eff_monojet_recoil_wen->Draw("EPsame");
  
  TLegend leg (0.6,0.3,0.9,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(eff_monojet_recoil_zvv,"Z #rightarrow #nu#nu","EPL");
  if(not applyGenSelection)
    leg.AddEntry(eff_monojet_recoil_wjet,"W-jet SR","EPL");
  leg.AddEntry(eff_monojet_recoil_wmn,"W #rightarrow #mu#nu","EPL");
  if(not computeAsFunctionOfMET)
    leg.AddEntry(eff_monojet_recoil_zmm,"Z #rightarrow #mu#mu","EPL");
  else
    leg.AddEntry(eff_monojet_recoil_wen,"W #rightarrow e#nu","EPL");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil.pdf").c_str(),"pdf");

  // Plotting final result for MC turn ons
  TH1* frame2 = canvas->DrawFrame(150,0.85,600,1.05, "");
  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Trigger Efficiency");
  frame2->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame2->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame2->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame2->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame2->GetXaxis()->SetTitleOffset(1.0);
  
  canvas->cd();
  frame2->Draw();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);
  
  fitfunc_monojet_recoil_zvv->Draw("Lsame");
  if(not applyGenSelection)
    fitfunc_monojet_recoil_wjet->Draw("Lsame");
  fitfunc_monojet_recoil_wmn->Draw("Lsame");
  if(not computeAsFunctionOfMET)
    fitfunc_monojet_recoil_zmm->Draw("Lsame");
  else
    fitfunc_monojet_recoil_wen->Draw("Lsame");

  eff_monojet_recoil_zvv->Draw("EPsame");
  if(not applyGenSelection)
    eff_monojet_recoil_wjet->Draw("EPsame");
  eff_monojet_recoil_wmn->Draw("EPsame");
  if(not computeAsFunctionOfMET)
    eff_monojet_recoil_zmm->Draw("EPsame");
  else
    eff_monojet_recoil_wen->Draw("EPsame");
  
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil_zoom.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil_zoom.pdf").c_str(),"pdf");

  TFile* outputFile = new TFile((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil.root").c_str(),"RECREATE");
  outputFile->cd();

  fitfunc_monojet_recoil_zvv->Write("func_zvv");
  if(not applyGenSelection)
    fitfunc_monojet_recoil_wjet->Write("func_wjet");
  fitfunc_monojet_recoil_wmn->Write("func_wmn");
  if(not computeAsFunctionOfMET)
    fitfunc_monojet_recoil_zmm->Write("func_zmm");
  else
    fitfunc_monojet_recoil_wen->Write("func_zee");

  eff_monojet_recoil_zvv->Write("efficiency_zvv");
  if(not applyGenSelection)
    eff_monojet_recoil_wjet->Write("efficiency_wjet");
  eff_monojet_recoil_wmn->Write("efficiency_wmn");
  if(not computeAsFunctionOfMET)
    eff_monojet_recoil_zmm->Write("efficiency_zmm");
  else
    eff_monojet_recoil_wen->Write("efficiency_wen");

  outputFile->mkdir("Distributions");
  outputFile->cd("Distributions");

  metCalo_SR_pas->Write("metCalo_SR_pas");
  metCalo_WM_pas->Write("metCalo_WM_pas");
  metCalo_ZM_pas->Write("metCalo_ZM_pas");
  metCalo_SR_fail->Write("metCalo_SR_fail");
  metCalo_WM_fail->Write("metCalo_WM_fail");
  metCalo_ZM_fail->Write("metCalo_ZM_fail");
  metPF_SR_pas->Write("metPF_SR_pas");
  metPF_WM_pas->Write("metPF_WM_pas");
  metPF_ZM_pas->Write("metPF_ZM_pas");
  metPF_SR_fail->Write("metPF_SR_fail");
  metPF_WM_fail->Write("metPF_WM_fail");
  metPF_ZM_fail->Write("metPF_ZM_fail");
  hnum_monojet_recoil_zvv->Write("recoil_Zvv_pass");
  hnum_monojet_recoil_wjet->Write("recoil_WJet_pass");
  hnum_monojet_recoil_wmn->Write("recoil_WM_pass");
  hnum_monojet_recoil_zmm->Write("recoil_ZM_pass");
  hden_monojet_recoil_zvv->Write("recoil_Zvv_fail");
  hden_monojet_recoil_wjet->Write("recoil_WJet_fail");
  hden_monojet_recoil_wmn->Write("recoil_WM_fail");
  hden_monojet_recoil_zmm->Write("recoil_ZM_fail");
  
  outputFile->Close();
  
}
