#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// recoil binning for monojet                                                                                                                                                                    
vector <float> bins_monojet_recoil    = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBFrelaxed_recoil = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBF_recoil        = {150,175,200,225,230,250,275,300,350,400,450,500,600,700,850,1000};
vector <float> bins_VBFrelaxed_mjj    = {200,400,600,800,1000,1250,1500,1750,2000,2500,3500,5000};
vector <float> bins_VBF_mjj           = {1300,1500,1750,2000,3500,3500,5000};

// cut for trigger values
static float L1_ETM_CUT  = 50;
static float caloMET_CUT = 70;
static float caloMETClean_CUT = 60;
static float caloMHT_CUT = 70;
static float PFMHT_CUT   = 90;
static float PFMET_CUT   = 90;
static float PFMHTNoMu_CUT = 90;
static float PFMETNoMu_CUT = 90;

// possible samples to be selected
enum class Sample {sig,wmn};
enum class Category {monojet,VBFrelaxed,VBF};


// k-factor file
static string kfactorFile = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_24bins.root";

// calculate the sum of weights for each file from the gentree
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

//////////
void makeTriggerAnalysis(vector<TTree*> & trees, 
			 vector<TH1F*> &  histograms, // depends on the selection one wants to apply
			 const Sample &   sample, // sample to be selected
			 const Category & category,
			 vector<double> & wgtsum, // sum of weights
			 vector<TH1*> &   khists, // NLO k-factors			 
			 const float &    luminosity, // luminosity
			 const TString &  triggerPath, // trigger Path to be studied
			 const string &   observable
			 ){
  

  string triggerString;

  // fix the threshold according to the trigger path
  if(triggerPath.Contains("PFMET90_PFMHT90")){
    caloMET_CUT = 70;
    caloMETClean_CUT = 60;
    caloMHT_CUT = 70;
    PFMHT_CUT   = 90;
    PFMET_CUT   = 90;
    triggerString = "hltmetwithmu90";
  }
  else if(triggerPath.Contains("PFMETNoMu90_PFMHTNoMu90")){
    caloMET_CUT = 80;
    caloMETClean_CUT = 70;
    caloMHT_CUT = 80;
    PFMHTNoMu_CUT = 90;
    PFMETNoMu_CUT = 90;
    triggerString = "hltmet90";
  }
  else if(triggerPath.Contains("PFMET100_PFMHT100") or triggerPath.Contains("PFMETNoMu100_PFMHTNoMu100")){
    caloMET_CUT = 80;
    caloMETClean_CUT = 70;
    caloMHT_CUT = 80;
    PFMHT_CUT   = 100;
    PFMET_CUT   = 100;
    PFMHTNoMu_CUT = 100;
    PFMETNoMu_CUT = 100;
    if(triggerPath.Contains("PFMETNoMu100_PFMHTNoMu100"))
      triggerString = "hltmet100";
    else
      triggerString = "hltmetwithmu100";
  }
  else if(triggerPath.Contains("PFMET110_PFMHT110") or triggerPath.Contains("PFMETNoMu110_PFMHTNoMu110")){
    caloMET_CUT = 80;
    caloMETClean_CUT = 70;
    caloMHT_CUT = 80;
    PFMHT_CUT   = 110;
    PFMET_CUT   = 110;
    PFMHTNoMu_CUT   = 110;
    PFMETNoMu_CUT   = 110;
    if(triggerPath.Contains("PFMETNoMu110_PFMHTNoMu110"))
      triggerString = "hltmet110";
    else
      triggerString = "hltmetwithmu110";
      
  }
  else if(triggerPath.Contains("PFMET120_PFMHT120") or triggerPath.Contains("PFMETNoMu120_PFMHTNoMu120")){
    caloMET_CUT = 90;
    caloMETClean_CUT = 80;
    caloMHT_CUT = 90;
    PFMHT_CUT   = 120;
    PFMET_CUT   = 120;
    PFMHTNoMu_CUT   = 120;
    PFMETNoMu_CUT   = 120;
    if(triggerPath.Contains("PFMET120_PFMHT120"))
      triggerString = "hltmetwithmu120";
    else
      triggerString = "hltmet120";
  }

  cout<<"Loop on trees "<<endl;

  /// pileup weight
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");
  puhist->SetDirectory(0);

  /// muon scale factors
  TFile* muonSF_file  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/muon_scalefactors.root");
  TH2* msfloose_id  = (TH2*) muonSF_file->Get("scalefactors_MuonLooseId_Muon");
  TH2* msfloose_iso = (TH2*) muonSF_file->Get("scalefactors_Iso_MuonLooseId");
  TH2* msftight_id  = (TH2*) muonSF_file->Get("scalefactors_TightId_Muon");
  TH2* msftight_iso = (TH2*) muonSF_file->Get("scalefactors_Iso_MuonTightId");
  msfloose_id->SetDirectory(0);
  msfloose_iso->SetDirectory(0);
  msftight_id->SetDirectory(0);
  msftight_iso->SetDirectory(0);

  int itree = 0;

  for(auto tree : trees){    

    TTreeReader reader(tree);
    TTreeReaderValue<float> xsec     (reader,"xsec");
    TTreeReaderValue<float> wgt      (reader,"wgt");
    TTreeReaderValue<UChar_t> hltmettrigger (reader,triggerString.c_str());
    TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
    TTreeReaderValue<float>   mu1pt  (reader,"mu1pt");
    TTreeReaderValue<float>   mu1eta (reader,"mu1eta");
    TTreeReaderValue<float>   mu1phi (reader,"mu1phi");
    TTreeReaderValue<int>     mu1id  (reader,"mu1id");
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
    TTreeReaderValue<vector<float> > jetpt     (reader,"combinejetpt");
    TTreeReaderValue<vector<float> > jeteta    (reader,"combinejeteta");
    TTreeReaderValue<vector<float> > jetphi    (reader,"combinejetphi");
    TTreeReaderValue<vector<float> > jetm      (reader,"combinejetm");
    TTreeReaderValue<vector<float> > jetchfrac (reader,"combinejetCHfrac");
    TTreeReaderValue<vector<float> > jetnhfrac (reader,"combinejetNHfrac");
    TTreeReaderValue<float> mmet        (reader,"t1mumet");
    TTreeReaderValue<float> mmetphi     (reader,"t1mumetphi");
    TTreeReaderValue<float> metpf       (reader,"pfmet");
    TTreeReaderValue<float> metcalo     (reader,"calomet");
    TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> trig_L1ETM_pt (reader,"trig_L1ETM_pt");
    TTreeReaderValue<vector<float> > trig_obj_pt (reader,"trig_obj_pt");
    TTreeReaderValue<vector<string> > trig_obj_col (reader,"trig_obj_col");

    //////////////////                                                                                                                                                                              
    long int nTotal = tree->GetEntries();
    cout<<"Looping on itree "<<itree<<" of "<<trees.size()<<" Total number of events: "<<nTotal<<endl;
    long int nEvents = 0;    
    long int nPart = 100000;

    while(reader.Next()){/////

      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;
      
      // to speed-up
      if(category == Category::monojet and *mmet < bins_monojet_recoil.front()) continue;
      else if(category == Category::VBFrelaxed and *mmet < bins_VBFrelaxed_recoil.front()) continue;
      else if(category == Category::VBF and *mmet < bins_VBF_recoil.front()) continue;

      if(observable == "mjj" and *mmet < 250) continue; // cut harder in met
 
      // object veto
      if(*nbjets   != 0)    continue;
      if(*ntaus    != 0)    continue;
      if(*nphotons != 0)    continue;
      if(*nelectrons != 0 ) continue;
      
      // muons
      if(sample == Sample::sig and *nmuons != 0) continue;
      else if(sample == Sample::wmn and *nmuons !=1) continue;
      
      // met filters
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
      if(fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;

      // ask single muon trigger in case of selecting wmn events
      if(sample == Sample::wmn and not *hltsinglemu) continue;

      // apply standard reco-level cuts for the control region
      if(sample == Sample::wmn){
	if(*mu1pt < 20) continue;
	if(fabs(*mu1eta) > 2.4) continue;
	if(*mu1id  !=1) continue;
	float dphi = fabs(*mu1phi-*mmetphi);
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	float mtw = sqrt(2*(*mu1pt)*(*mmet)*(1-cos(dphi)));
	if(mtw > 160) continue;
      }

      // apply jet-met dphi--> not for gen level analysis
      if(*jmmdphi < 0.5) continue;

      TLorentzVector jet1,jet2;
	  
      //// monojet category
      if(category == Category::monojet){
	// apply jet pt selections
	if(*nincjets < 1) continue;
	if(jetpt->size() == 0) continue;
	if(jetpt->at(0) < 100) continue;
	if(fabs(jeteta->at(0)) > 2.5) continue;
	if(jetchfrac->at(0) < 0.1) continue;
	if(jetnhfrac->at(0) > 0.8) continue;
	
      }
      else if(category == Category::VBFrelaxed){
	if(*nincjets < 1) continue;
	if(jetpt->size() <= 1) continue;
	if(jetpt->at(0) < 80) continue;
	if(jetpt->at(1) < 40) continue;
	if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
	if(fabs(jeteta->at(0)) < 2.4 and jetchfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.4 and jetnhfrac->at(0) > 0.8) continue;
	if(jeteta->at(0)*jeteta->at(1) > 0) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < 1) continue;

 	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
 	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	
	if((jet1+jet2).M() < 200) continue;
	if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      }
      else if(category == Category::VBF){
	if(*nincjets < 1) continue;
	if(jetpt->size() <= 1) continue;
	if(jetpt->at(0) < 80) continue;
	if(jetpt->at(1) < 40) continue;
	if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
	if(fabs(jeteta->at(0)) < 2.4 and jetchfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.4 and jetnhfrac->at(0) > 0.8) continue;
	if(jeteta->at(0)*jeteta->at(1) > 0) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < 3) continue;

 	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
 	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	
	if((jet1+jet2).M() < 1300) continue;
	if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      }

      // pileup re-weight
      double puwgt = 1;
      if(*nvtx < 60)
	puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
      
      // apply muon scale factors
      double sfwgt = 1.;
      if(*mu1pt > 0. and sample == Sample::wmn){
	float ptValue = *mu1pt;
	if(ptValue < msftight_id->GetYaxis()->GetBinLowEdge(1)) 
	  ptValue =  msftight_id->GetYaxis()->GetBinLowEdge(1)+1;
	else if(ptValue > msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)) 
	  ptValue = msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)-1;
	
	if(*mu1id == 1)
	  sfwgt *= msftight_id->GetBinContent(msftight_id->FindBin(fabs(*mu1eta),ptValue))*msftight_iso->GetBinContent(msftight_iso->FindBin(fabs(*mu1eta),ptValue));
	else
	  sfwgt *= msfloose_id->GetBinContent(msfloose_id->FindBin(fabs(*mu1eta),ptValue))*msfloose_iso->GetBinContent(msfloose_iso->FindBin(fabs(*mu1eta),ptValue));
      }
      
      // calculate NLO k-factors
      Double_t kwgt = 1.0;
      double   genpt = *wzpt;
      for (size_t i = 0; i < khists.size(); i++) {
	if (khists[i]) {
	  if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	  if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	  kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));	  
	}
      }
	
      // found numerical values
      float caloMET = 0;
      float caloMETClean = 0;
      float caloMHT = 0;
      float PFMHT = 0;
      float PFMET = 0;
      float PFMHTNoMu = 0;
      float PFMETNoMu = 0;
      
      for(size_t itrig = 0; itrig < trig_obj_col->size(); itrig++){
	TString name (trig_obj_col->at(itrig));
	if(name.Contains("hltMetClean")) caloMETClean = trig_obj_pt->at(itrig);
	else if(name.Contains("hltMet")) caloMET = trig_obj_pt->at(itrig);
	else if(name.Contains("hltMht")) caloMHT = trig_obj_pt->at(itrig);
	else if(name.Contains("hltPFMHTNoMuTightID")) PFMHTNoMu = trig_obj_pt->at(itrig);
	else if(name.Contains("hltPFMETNoMuProducer")) PFMETNoMu = trig_obj_pt->at(itrig);
	else if(name.Contains("hltPFMHTTightID")) PFMHT = trig_obj_pt->at(itrig);
	else if(name.Contains("hltPFMETProducer")) PFMET = trig_obj_pt->at(itrig);
      }

      // to include the overflow
      float value = 0;
      if(observable == "met"){
	value = *mmet;
	if(category == Category::monojet and value > bins_monojet_recoil.back())
	  value = bins_monojet_recoil.back()-1;	
	else if(category == Category::VBFrelaxed and value > bins_VBFrelaxed_recoil.back())
	  value = bins_VBFrelaxed_recoil.back()-1;
	else if(category == Category::VBF and value > bins_VBF_recoil.back())
	  value = bins_VBF_recoil.back()-1;
      }
      else if(observable == "mjj"){
	value = (jet1+jet2).M();
	if(category == Category::VBFrelaxed and value > bins_VBFrelaxed_mjj.back())
	  value = bins_VBFrelaxed_mjj.back()-1;
	else if(category == Category::VBF and value > bins_VBF_mjj.back())
	  value = bins_VBF_mjj.back()-1;
      }

      /// Fill histograms
      for(auto hist : histograms){
	
	TString name (hist->GetName());
	if(name.Contains("Inclusive")) // fill with all the events i.e. numerator
	  hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	
	// L1 efficiency
	else if(name.Contains("L1ETM")){ // fill with all the events passing L1 selection
	  if(*trig_L1ETM_pt > L1_ETM_CUT) 
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	
	// L1 + caloMET + caloMETClean
	else if(name.Contains("caloMETClean")){// fill with all the events passing calo-MET + calo-MET clean cut
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT) 
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	
	// L1 + caloMET
	else if(name.Contains("caloMET")){ // fill with all the events passing calo-MET 
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT)
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	
	// L1 + caloMET + caloMETClean + caloMHT
	else if(name.Contains("caloMHT")){// fill with all the events passing calo MHT
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT)
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu
	else if(name.Contains("PFMHTNoMu")){
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and 
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT and
	     PFMHTNoMu > PFMHTNoMu_CUT)
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	  
	// L1 + caloMET + caloMETClean + caloMHT + PFMHT
        else if(name.Contains("PFMHT")){  // fill with all the events passing PFMHT
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT and
	     PFMHT > PFMHT_CUT)
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu + PFMETNoMu
	else if(name.Contains("PFMETNoMu")){
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT   > caloMHT_CUT and
	     PFMHTNoMu > PFMHTNoMu_CUT and
	     PFMETNoMu > PFMETNoMu_CUT)
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));	    
	}	  
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHT + PFMET
	else if(name.Contains("PFMET")){ // fill with all the events passing PFMET
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT and
	     PFMHT > PFMHT_CUT and
	     PFMET > PFMET_CUT) 
	    hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));	  
        } 
	
	// total path
        else if(name.Contains("LastFilter")){
 	 if(*trig_L1ETM_pt > L1_ETM_CUT and
	    *hltmettrigger != 0) 
	   hist->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(itree));	  
        } 
      }            
    }
    cout<<endl;
    itree++;
  }
  
  if(pufile) pufile->Close();
  if(muonSF_file) muonSF_file->Close();
}


//// -----
void fillTreeList(vector<TFile*> & file, vector<TTree*> & tree, vector<TTree*> & gentree, const string & inputDIR, const string & postfix){

  cout<<"Read list of files for "<<postfix<<" process "<<endl;
  system(("ls "+inputDIR+" | grep "+postfix+" > list_dir.txt").c_str());
  ifstream file_dir("list_dir.txt");
  if(file_dir.is_open()){
    string line;
    while(!file_dir.eof()){
      getline(file_dir,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_("list.txt");
      if(file_.is_open()){
	string line2;
	while(!file_.eof()){
	  getline(file_,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open "<<postfix<<" file with name: "<<line2<<endl;
	  file.push_back(TFile::Open(line2.c_str()));
	  tree.push_back((TTree*) file.back()->Get("tree/tree"));
	  gentree.push_back((TTree*) file.back()->Get("gentree/gentree"));
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");
}

/// plotting efficiency
void plotEfficiency(TCanvas* canvas, TEfficiency* histo_zjet, TEfficiency* histo_wjet, TEfficiency* histo_wmn, const string & outputDIR, const float &  luminosity){
  
  // Make efficiencies as TGraph
  TGraphAsymmErrors* graph_zjet = histo_zjet->CreateGraph();
  TGraphAsymmErrors* graph_wjet = histo_wjet->CreateGraph();
  TGraphAsymmErrors* graph_wmn = histo_wmn->CreateGraph();

  graph_zjet->SetMarkerColor(kBlack);
  graph_zjet->SetLineColor(kBlack);
  graph_zjet->SetMarkerStyle(20);
  graph_zjet->SetMarkerSize(0.75);
  graph_zjet->SetLineWidth(1);

  graph_wjet->SetMarkerColor(kRed);
  graph_wjet->SetLineColor(kRed);
  graph_wjet->SetMarkerStyle(24);
  graph_wjet->SetMarkerSize(0.75);
  graph_wjet->SetLineWidth(1);

  graph_wmn->SetMarkerColor(kBlue);
  graph_wmn->SetLineColor(kBlue);
  graph_wmn->SetMarkerStyle(20);
  graph_wmn->SetMarkerStyle(26);
  graph_wmn->SetMarkerSize(0.75);
  graph_wmn->SetLineWidth(1);
    
  // Plotting final result for MC turn ons
  TH1* frame = canvas->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.6,
				 histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.1, "");
  frame->GetXaxis()->SetTitle("Recoil [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.30);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());

  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.3);
  canvas->Draw();
  canvas->cd();
  frame->Draw();

  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  graph_zjet->Draw("EPLsame");
  graph_wjet->Draw("EPLsame");
  graph_wmn->Draw("EPLsame");
  
  TLegend leg (0.6,0.4,0.9,0.6);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(graph_zjet,"Z(#nu#nu)+jets SR","EPL");
  leg.AddEntry(graph_wjet,"W+jets SR","EPL");
  leg.AddEntry(graph_wmn,"W+jets W(#mu#nu)-CR","EPL");
  leg.Draw("same");

  TString plotName(histo_zjet->GetName());
  plotName.ReplaceAll("_zvv","");

  //// ----
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();

  // Plotting final result for MC turn ons
  TH1* frame2 = pad2->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.9,
				histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.05, "");
  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.32);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetNdivisions(504);
  frame2->GetYaxis()->SetLabelSize(0.85*frame2->GetYaxis()->GetLabelSize());
  frame2->Draw();

  /// ------
  TH1F* ratio_wjet_over_zjet = (TH1F*) histo_wjet->GetPassedHistogram();
  ratio_wjet_over_zjet->Divide(histo_wjet->GetTotalHistogram());
  TH1F* ratio_temp = (TH1F*) histo_zjet->GetPassedHistogram();
  ratio_temp->Divide(histo_zjet->GetTotalHistogram());
  ratio_wjet_over_zjet->Divide(ratio_temp);

  /// ------
  TH1F* ratio_wmn_over_zjet = (TH1F*) histo_wmn->GetPassedHistogram();
  ratio_wmn_over_zjet->Divide(histo_wmn->GetTotalHistogram());
  ratio_wmn_over_zjet->Divide(ratio_temp);
  
  ratio_wjet_over_zjet->SetMarkerColor(kRed);
  ratio_wjet_over_zjet->SetLineColor(kRed);
  ratio_wjet_over_zjet->SetMarkerStyle(24);
  ratio_wjet_over_zjet->SetMarkerSize(0.75);
  ratio_wjet_over_zjet->SetLineWidth(1);
  ratio_wjet_over_zjet->Draw("EPsame");

  ratio_wmn_over_zjet->SetMarkerColor(kBlue);
  ratio_wmn_over_zjet->SetLineColor(kBlue);
  ratio_wmn_over_zjet->SetMarkerStyle(20);
  ratio_wmn_over_zjet->SetMarkerStyle(26);
  ratio_wmn_over_zjet->SetMarkerSize(0.75);
  ratio_wmn_over_zjet->SetLineWidth(1);
  ratio_wmn_over_zjet->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".root").c_str(),"root");

  if(frame2) delete frame2;
  if(pad2) delete pad2;
  
}


/// Main function 
void makeMETTriggerEfficiencyMC_ZvsW(string inputDIR, 
				     string outputDIR, 
				     string triggerPath, /// PFMET110_PFMHT110, PFMET100_PFMHT100, PFMET120_PFMHT120, 
				     Category category,
				     string observable, // study the trigger efficiency vs mjj or recoil
				     bool   isRelativeEff = false,
				     float  luminosity = 35.9){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  // check trigger path
  if(triggerPath != "PFMET90_PFMHT90" and
     triggerPath != "PFMET100_PFMHT100" and 
     triggerPath != "PFMET110_PFMHT110" and 
     triggerPath != "PFMET120_PFMHT120" and 
     triggerPath != "PFMETNoMu90_PFMHTNoMu90" and
     triggerPath != "PFMETNoMu100_PFMHTNoMu100" and 
     triggerPath != "PFMETNoMu110_PFMHTNoMu110" and 
     triggerPath != "PFMETNoMu120_PFMHTNoMu120"){
    cerr<<"######### Trigger path not known --> exit "<<endl;
    return;
  }

  // check the observable
  if(observable != "met" and
     observable != "mjj"){
    cerr<<"Not recognized observable --> exit "<<endl;
    return;
  }
   
  // further checks
  if(category == Category::monojet and observable == "mjj"){
    cerr<<"Cannot always plot mjj in monojet events --> exit "<<endl;
    return;
  }
  
  // input tree                                                                                                                                                                                       
  vector<TFile*> filelist_zjet;
  vector<TTree*> tree_zjet;
  vector<TTree*> gentree_zjet;

  vector<TFile*> filelist_wjet;
  vector<TTree*> tree_wjet;
  vector<TTree*> gentree_wjet;

  fillTreeList(filelist_zjet,tree_zjet,gentree_zjet,inputDIR,"ZJets");
  fillTreeList(filelist_wjet,tree_wjet,gentree_wjet,inputDIR,"WJets");
 
  // sum of weights
  vector<double> wgtsum_zjet;
  vector<double> wgtsum_wjet;
  cout<<"######### Calculate Weights for : Z+jets "<<endl;
  calculateSumWeight(gentree_zjet,wgtsum_zjet);
  cout<<"######### Calculate Weights for : W+jets "<<endl;
  calculateSumWeight(gentree_zjet,wgtsum_wjet);

  /////// start analysis  
  // Distributions for Z+jets

  vector<float> binning;
  if(category == Category::monojet and observable == "met") binning = bins_monojet_recoil;
  else if(category == Category::VBFrelaxed and observable == "met") binning = bins_VBFrelaxed_recoil;
  else if(category == Category::VBF and observable == "met") binning = bins_VBF_recoil;
  
  if(category == Category::VBFrelaxed and observable == "mjj") binning = bins_VBFrelaxed_mjj;
  else if(category == Category::VBF and observable == "mjj") binning = bins_VBF_mjj;
  
  TH1F* h_inclusive_zjet = new TH1F("h_Inclusive_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_L1ETM_zjet     = new TH1F("h_L1ETM_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMET_zjet   = new TH1F("h_caloMET_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMETClean_zjet     = new TH1F("h_caloMETClean_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMHT_zjet   = new TH1F("h_caloMHT_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHT_zjet     = new TH1F("h_PFMHT_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMET_zjet     = new TH1F("h_PFMET_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHTNoMu_zjet = new TH1F("h_PFMHTNoMu_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMETNoMu_zjet = new TH1F("h_PFMETNoMu_zjet", "", binning.size()-1, &binning[0]);
  TH1F* h_LastFilter_zjet = new TH1F("h_LastFilter_zjet", "", binning.size()-1, &binning[0]);

  h_inclusive_zjet->Sumw2();
  h_L1ETM_zjet->Sumw2();
  h_caloMET_zjet->Sumw2();
  h_caloMETClean_zjet->Sumw2();
  h_caloMHT_zjet->Sumw2();
  h_PFMHT_zjet->Sumw2();
  h_PFMET_zjet->Sumw2();
  h_PFMHTNoMu_zjet->Sumw2();
  h_PFMETNoMu_zjet->Sumw2();
  h_LastFilter_zjet->Sumw2();

  h_inclusive_zjet->SetDirectory(0);
  h_L1ETM_zjet->SetDirectory(0);
  h_caloMET_zjet->SetDirectory(0);
  h_caloMETClean_zjet->SetDirectory(0);
  h_caloMHT_zjet->SetDirectory(0);
  h_PFMHT_zjet->SetDirectory(0);
  h_PFMET_zjet->SetDirectory(0);
  h_PFMHTNoMu_zjet->SetDirectory(0);
  h_PFMETNoMu_zjet->SetDirectory(0);
  h_LastFilter_zjet->SetDirectory(0);

  // Distributions for W+jets SR
  TH1F* h_inclusive_wjet = new TH1F("h_Inclusive_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_L1ETM_wjet     = new TH1F("h_L1ETM_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMET_wjet   = new TH1F("h_caloMET_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMETClean_wjet     = new TH1F("h_caloMETClean_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMHT_wjet   = new TH1F("h_caloMHT_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHT_wjet     = new TH1F("h_PFMHT_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMET_wjet     = new TH1F("h_PFMET_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHTNoMu_wjet = new TH1F("h_PFMHTNoMu_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMETNoMu_wjet = new TH1F("h_PFMETNoMu_wjet", "", binning.size()-1, &binning[0]);
  TH1F* h_LastFilter_wjet = new TH1F("h_LastFilter_wjet", "", binning.size()-1, &binning[0]);

  h_inclusive_wjet->Sumw2();
  h_L1ETM_wjet->Sumw2();
  h_caloMET_wjet->Sumw2();
  h_caloMETClean_wjet->Sumw2();
  h_caloMHT_wjet->Sumw2();
  h_PFMHT_wjet->Sumw2();
  h_PFMET_wjet->Sumw2();
  h_PFMHTNoMu_wjet->Sumw2();
  h_PFMETNoMu_wjet->Sumw2();
  h_LastFilter_wjet->Sumw2();

  h_inclusive_wjet->SetDirectory(0);
  h_L1ETM_wjet->SetDirectory(0);
  h_caloMET_wjet->SetDirectory(0);
  h_caloMETClean_wjet->SetDirectory(0);
  h_caloMHT_wjet->SetDirectory(0);
  h_PFMHT_wjet->SetDirectory(0);
  h_PFMET_wjet->SetDirectory(0);
  h_PFMHTNoMu_wjet->SetDirectory(0);
  h_PFMETNoMu_wjet->SetDirectory(0);
  h_LastFilter_wjet->SetDirectory(0);

  // Distributions for W+jets Wmn
  TH1F* h_inclusive_wmn = new TH1F("h_Inclusive_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_L1ETM_wmn     = new TH1F("h_L1ETM_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMET_wmn   = new TH1F("h_caloMET_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMETClean_wmn = new TH1F("h_caloMETClean_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_caloMHT_wmn   = new TH1F("h_caloMHT_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHT_wmn     = new TH1F("h_PFMHT_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMET_wmn     = new TH1F("h_PFMET_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMHTNoMu_wmn = new TH1F("h_PFMHTNoMu_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_PFMETNoMu_wmn = new TH1F("h_PFMETNoMu_wmn", "", binning.size()-1, &binning[0]);
  TH1F* h_LastFilter_wmn = new TH1F("h_LastFilter_wmn", "", binning.size()-1, &binning[0]);

  h_inclusive_wmn->Sumw2();
  h_L1ETM_wmn->Sumw2();
  h_caloMET_wmn->Sumw2();
  h_caloMETClean_wmn->Sumw2();
  h_caloMHT_wmn->Sumw2();
  h_PFMHT_wmn->Sumw2();
  h_PFMET_wmn->Sumw2();
  h_PFMHTNoMu_wmn->Sumw2();
  h_PFMETNoMu_wmn->Sumw2();
  h_LastFilter_wmn->Sumw2();

  h_inclusive_wmn->SetDirectory(0);
  h_L1ETM_wmn->SetDirectory(0);
  h_caloMET_wmn->SetDirectory(0);
  h_caloMETClean_wmn->SetDirectory(0);
  h_caloMHT_wmn->SetDirectory(0);
  h_PFMHT_wmn->SetDirectory(0);
  h_PFMET_wmn->SetDirectory(0);
  h_PFMHTNoMu_wmn->SetDirectory(0);
  h_PFMETNoMu_wmn->SetDirectory(0);
  h_LastFilter_wmn->SetDirectory(0);

  //// ----
  vector<TH1F*> histo_zjet;
  vector<TH1F*> histo_wjet;
  vector<TH1F*> histo_wmn;

  histo_zjet.push_back(h_inclusive_zjet);
  histo_zjet.push_back(h_L1ETM_zjet);
  histo_zjet.push_back(h_caloMET_zjet);
  histo_zjet.push_back(h_caloMETClean_zjet);
  histo_zjet.push_back(h_caloMHT_zjet);

  histo_wjet.push_back(h_inclusive_wjet);
  histo_wjet.push_back(h_L1ETM_wjet);
  histo_wjet.push_back(h_caloMET_wjet);
  histo_wjet.push_back(h_caloMETClean_wjet);
  histo_wjet.push_back(h_caloMHT_wjet);

  histo_wmn.push_back(h_inclusive_wmn);
  histo_wmn.push_back(h_L1ETM_wmn);
  histo_wmn.push_back(h_caloMET_wmn);
  histo_wmn.push_back(h_caloMETClean_wmn);
  histo_wmn.push_back(h_caloMHT_wmn);

  if(TString(triggerPath).Contains("PFMETNoMu")){	
    histo_zjet.push_back(h_PFMHTNoMu_zjet);
    histo_zjet.push_back(h_PFMETNoMu_zjet);

    histo_wjet.push_back(h_PFMHTNoMu_wjet);
    histo_wjet.push_back(h_PFMETNoMu_wjet);

    histo_wmn.push_back(h_PFMHTNoMu_wmn);
    histo_wmn.push_back(h_PFMETNoMu_wmn);
  }
  else{
    histo_zjet.push_back(h_PFMHT_zjet);
    histo_zjet.push_back(h_PFMET_zjet);

    histo_wjet.push_back(h_PFMHT_wjet);
    histo_wjet.push_back(h_PFMET_wjet);

    histo_wmn.push_back(h_PFMHT_wmn);
    histo_wmn.push_back(h_PFMET_wmn);

  }

  histo_zjet.push_back(h_LastFilter_zjet);
  histo_wjet.push_back(h_LastFilter_wjet);
  histo_wmn.push_back(h_LastFilter_wmn);

  // k-factors
  TFile* kffile = TFile::Open(kfactorFile.c_str());
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
 
  vector<TH1*> khists_zjet;
  vector<TH1*> khists_wjet;

  khists_zjet.push_back(zewkhist);
  khists_zjet.push_back(znlohist);
  khists_wjet.push_back(wewkhist);
  khists_wjet.push_back(wnlohist);

  //// ----- 
  cout<<"######### Loop on Z+jets trees for SR selection "<<endl;
  makeTriggerAnalysis(tree_zjet,histo_zjet,Sample::sig,category,wgtsum_zjet,khists_zjet,luminosity,triggerPath,observable);
  cout<<"######### Loop on W+jets trees for SR selection "<<endl;
  makeTriggerAnalysis(tree_wjet,histo_wjet,Sample::sig,category,wgtsum_wjet,khists_wjet,luminosity,triggerPath,observable);
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_wjet,histo_wmn,Sample::wmn,category,wgtsum_wjet,khists_wjet,luminosity,triggerPath,observable);

  
  ////---
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 650);
  canvas->cd();

  // calculate the efficiency
  vector<TEfficiency*> eff_zjet;
  vector<TEfficiency*> eff_wjet;
  vector<TEfficiency*> eff_wmn;

  if(isRelativeEff){ // relative efficiency i.e. filter wrt the previous one

    for(size_t ihist = 0; ihist < histo_zjet.size()-1; ihist++){
      
      TString name (histo_zjet.at(ihist+1)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_zjet.push_back(new TEfficiency(*histo_zjet.at(ihist+1),*histo_zjet.at(ihist)));
      eff_zjet.back()->SetName(name.Data());

      name  = TString(histo_wjet.at(ihist+1)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_wjet.push_back(new TEfficiency(*histo_wjet.at(ihist+1),*histo_wjet.at(ihist)));
      eff_wjet.back()->SetName(name.Data());

      name = TString(histo_wmn.at(ihist+1)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_wmn.push_back(new TEfficiency(*histo_wmn.at(ihist+1),*histo_wmn.at(ihist)));
      eff_wmn.back()->SetName(name.Data());
    }
  }
  else{ // absolute efficiency of each trigger filter

    for(size_t ihist = 1; ihist < histo_zjet.size(); ihist++){

      TString name (histo_zjet.at(ihist)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_zjet.push_back(new TEfficiency(*histo_zjet.at(ihist),*histo_zjet.at(0)));
      eff_zjet.back()->SetName(name.Data());

      name = TString(histo_wjet.at(ihist)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_wjet.push_back(new TEfficiency(*histo_wjet.at(ihist),*histo_wjet.at(0)));
      eff_wjet.back()->SetName(name.Data());

      name = TString(histo_wmn.at(ihist)->GetName());
      name.ReplaceAll("h_","eff_");
      eff_wmn.push_back(new TEfficiency(*histo_wmn.at(ihist),*histo_wmn.at(0)));
      eff_wmn.back()->SetName(name.Data());
      
    }
  }

  // plot Efficiency
  for(size_t ihist = 0; ihist < eff_zjet.size(); ihist++)
    plotEfficiency(canvas,eff_zjet.at(ihist),eff_wjet.at(ihist),eff_wmn.at(ihist),outputDIR,luminosity);

}
