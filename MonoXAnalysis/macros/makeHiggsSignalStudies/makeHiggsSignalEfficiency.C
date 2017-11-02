#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "../makeTemplates/histoUtils.h"

static float luminosity = 35.9;

//////////////////////////
void makeHiggsSignalEfficiency(string inputFileName, Category category, bool isGen = false){

  // Define cut values depending on the category 
  float jetmetdphi_cut = 0.5;
  float jetpt1_cut     = 0;
  float jetpt2_cut     = 40;
  float jeteta1_cut    = 4.7;
  float jeteta2_cut    = 4.7;
  float met_cut        = 0;
  float dphijj_cut     = 0;
  float detajj_cut     = 0;
  float mjj_cut        = 0;

  if(category == Category::monojet or category == Category::monoV){
    jetpt1_cut  = 100;
    jeteta1_cut = 2.5;
    met_cut     = 250;
  }
  else if(category == Category::VBFrelaxed or category == Category::VBF){
    met_cut    = 200;
    jetpt1_cut = 80;
    if(category == Category::VBF){
      detajj_cut = 4.0;
      mjj_cut    = 1300;
      dphijj_cut  = 1.5;
    }
    else if(category == Category::VBFrelaxed){
      detajj_cut = 1.0;
      mjj_cut    = 200;
      dphijj_cut = 1.5;
    }
  }

  ////////////////////////  
  
  TFile* inputFile = TFile::Open(inputFileName.c_str());  
  TTree* tree = (TTree*) inputFile->Get("tree/tree");
  
  TTreeReader reader(tree);
  TTreeReaderValue<unsigned int> run    (reader,"run");
  TTreeReaderValue<unsigned int> event  (reader,"event");
  TTreeReaderValue<unsigned int> nvtx   (reader,"nvtx");
  TTreeReaderValue<float> xsec          (reader,"xsec");
  TTreeReaderValue<float> wgt           (reader,"wgt");

  /////////////// triggers 
  TTreeReaderValue<UChar_t> hltm90      (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90    (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm120   (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (reader,"hltmetwithmu300");

  /////////////// met filters
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");

  // njets
  TTreeReaderValue<unsigned int> njets       (reader,"njets");
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjetslowpt (reader,"nbjetslowpt");

  // AK4 jets
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<vector<float> > emfrac  (reader,"combinejetEMfrac");

  // Boosted jets
  TTreeReaderValue<vector<float> > boostedJetpt    (reader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (reader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetm     (reader,"boostedJetm");
  TTreeReaderValue<vector<float> > prunedJetm      (reader,"prunedJetm");
  TTreeReaderValue<vector<float> > boostedJettau2  (reader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (reader,"boostedJettau1");

  // MET
  TTreeReaderValue<float> met  (reader,"t1pfmet");
  TTreeReaderValue<float> metcalo (reader,"calomet");
  TTreeReaderValue<float> metphi  (reader,"t1pfmetphi");

  // Dphi
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<unsigned int> nmuons     (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (reader,"ntausold");
  TTreeReaderValue<unsigned int> nphotons   (reader,"nphotons");

  // gen jets
  TTreeReaderValue<vector<float> > genjetpt   (reader,"combinejetGenpt");
  TTreeReaderValue<vector<float> > genjeteta  (reader,"combinejetGeneta");
  TTreeReaderValue<vector<float> > genjetphi  (reader,"combinejetGenphi");
  TTreeReaderValue<vector<float> > genjetm    (reader,"combinejetGenm");

  // -----------
  TTreeReaderValue<float> genmet    (reader,"genmet");
  TTreeReaderValue<float> genmetphi (reader,"genmetphi");

  // -----------
  TTreeReaderValue<vector<float> > boostedJetGenpt    (reader,"boostedJetGenpt");
  TTreeReaderValue<vector<float> > boostedJetGeneta   (reader,"boostedJetGeneta");
  TTreeReaderValue<vector<float> > boostedJetGenm     (reader,"boostedJetGenm");
  TTreeReaderValue<vector<float> > prunedJetGenm      (reader,"prunedJetGenm");
  TTreeReaderValue<vector<float> > boostedJetGentau2  (reader,"boostedJetGentau2");
  TTreeReaderValue<vector<float> > boostedJetGentau1  (reader,"boostedJetGentau1");

  // calculate sum of weight
  double sumwgt = 0;
  while(reader.Next()){
    sumwgt += *wgt;
  }
  cout<<"Sum of weights "<<sumwgt<<endl;


  // -----------
  double n_total = 0;
  double n_genlevel = 0;
  double n_recolevel = 0;
  reader.SetEntry(0);

  // trigger weight file
  TFile* triggerfile_MET     = NULL;
  vector<TFile*> triggerfile_MET_binned_Wmn;
  if(category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed){ // monojet                                                                                     
      // Use Wmn or Wen for W+jets and SR                                                                                                                                                              
    triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet.root");
  }
  else{
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_0.0_800.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_800.0_1200.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1200.0_1700.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1700.0_3000.0.root"));
  }

  // single turn on                                                                                                                                                                                    
  TEfficiency*       triggermet       = NULL;
  TGraphAsymmErrors* triggermet_graph = NULL;
  
  if(triggerfile_MET != NULL){
    triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency_vbf_loose");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency");
    triggermet_graph = triggermet->CreateGraph();
  }

  vector<TF1*> triggermet_func_binned_Wmn;
  if(triggerfile_MET_binned_Wmn.size() != 0){
    for(auto ifile : triggerfile_MET_binned_Wmn)
      triggermet_func_binned_Wmn.push_back((TF1*) ifile->Get("efficiency_func"));
  }


  // pileup weight file
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1*  puhist = (TH1*) pufile->Get("puhist");

  ////////////////////
  while(reader.Next()){    

    double puwgt = 1;
    if(*nvtx < 60) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    
    double trigeff = 1;
    double pfmet = *met;

    // met trigger scale factor                                                                                                                                                                      
    if (not isGen) {
      // single trigger turn on to be applied                                                                                                                                                          
      if(triggermet_graph)
        trigeff *= triggermet_graph->Eval(min(pfmet,triggermet_graph->GetXaxis()->GetXmax()));
      // for VBF                                                                                                                                                                                       
      else if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed){
        if(jetpt->size() >= 2){
          TLorentzVector jet1 ;
          TLorentzVector jet2 ;
          jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
          jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	  if((jet1+jet2).M() < 800)
	    trigeff *= triggermet_func_binned_Wmn.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(0)->GetXaxis()->GetXmax()));
	  else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	    trigeff *= triggermet_func_binned_Wmn.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(1)->GetXaxis()->GetXmax()));
	  else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	    trigeff *= triggermet_func_binned_Wmn.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(2)->GetXaxis()->GetXmax()));
	  else if((jet1+jet2).M() >= 1700)
              trigeff *= triggermet_func_binned_Wmn.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(3)->GetXaxis()->GetXmax()));
	}
      }
    }
    
    // no experimental weights
    n_total += *wgt*(*xsec)*luminosity/(sumwgt);
        
    //////////////////////
    if(isGen){

      ///////////
      if(category == Category::monojet){       
	if(genjetpt->size() <= 0) continue;
	if(genjetpt->at(0)  < jetpt1_cut) continue;
	if(fabs(genjeteta->at(0)) > jeteta1_cut) continue;
	if(*genmet < met_cut) continue;

	n_genlevel += *wgt*(*xsec)*luminosity/(sumwgt);	
      }
      ///////////
      else if(category == Category::monoV){
	// still apply the monojet id on the leading jet
	if(genjetpt->size() <= 0) continue;
	if(genjetpt->at(0)  < jetpt1_cut) continue;
	if(fabs(genjeteta->at(0)) > jeteta1_cut) continue;
	if(*genmet < met_cut) continue;
	if(boostedJetGenpt->size() <= 0) continue;
	if(boostedJetGenpt->at(0)  < 250) continue;
	if(fabs(boostedJetGeneta->at(0)) > 2.4 ) continue;
	if(boostedJetGentau2->at(0)/boostedJetGentau1->at(0) > 0.6) continue;
	if(prunedJetGenm->at(0) < 65 or prunedJetGenm->at(0) > 105 ) continue;
	
	n_genlevel += *wgt*(*xsec)*luminosity/(sumwgt);	
      }
      ///////////
      else if(category == Category::VBF or category == Category::VBFrelaxed){
      
      if(genjetpt->size() < 2) continue;
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(genjetpt->at(0),genjeteta->at(0),genjetphi->at(0),genjetm->at(0));
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(genjetpt->at(1),genjeteta->at(1),genjetphi->at(1),genjetm->at(1));      
      if(leadingJet.Pt() < jetpt1_cut) continue;
      if(fabs(leadingJet.Eta()) > jeteta1_cut) continue;
      if(subleadingJet.Pt() < jetpt2_cut) continue;
      if(fabs(subleadingJet.Eta()) > jeteta2_cut) continue;
      if(fabs(subleadingJet.Eta()) > 3 and fabs(leadingJet.Eta()) > 3 ) continue;
      if(*genmet < met_cut) continue;
      if(leadingJet.Eta()*subleadingJet.Eta() > 0) continue;

      if(fabs(leadingJet.Eta()-subleadingJet.Eta()) < detajj_cut) continue;
      if((leadingJet+subleadingJet).M() < mjj_cut) continue; 
      if(fabs(leadingJet.DeltaPhi(subleadingJet)) > dphijj_cut) continue;

      n_genlevel += *wgt*(*xsec)*luminosity/(sumwgt);
      
      }    
    }
    else{

      // met filters
      if(*fhbhe == 0 or *fhbiso == 0 or *fcsc == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0) continue; 
      if(*nmuons != 0) continue;
      if(*nelectrons != 0) continue;
      if(*nphotons != 0) continue;
      if(*ntaus != 0) continue;
      ///////////
      if(*nbjetslowpt != 0) continue;
      
      ///////////
      if(category == Category::monojet){
	
	if(jetpt->size() <= 0) continue;
	if(jetpt->at(0)  < jetpt1_cut) continue;
	if(fabs(jeteta->at(0)) > jeteta1_cut) continue;
	if(chfrac->at(0) < 0.1) continue;
	if(nhfrac->at(0) > 0.8) continue;      
	if(fabs(*jmmdphi) < jetmetdphi_cut) continue;
	if(*met < met_cut) continue;
	if(fabs(*met-*metcalo)/(*met) > 0.5) continue;
	if(boostedJetpt->size() > 0 and fabs(boostedJeteta->at(0)) < 2.4 and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and boostedJetpt->at(0) > 250 and prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105) continue;

	n_recolevel += *wgt*(*xsec)*luminosity*puwgt*trigeff/(sumwgt);
	
      }
      ///////////
      else if(category == Category::monoV){
	
	// still apply the monojet id on the leading jet
	if(*njets < 1) continue;
	if(jetpt->size() <= 0) continue;
	if(jetpt->at(0)  < jetpt1_cut) continue;
	if(fabs(jeteta->at(0)) > jeteta1_cut) continue;
	if(chfrac->at(0) < 0.1) continue;
	if(nhfrac->at(0) > 0.8) continue;
	if(fabs(*jmmdphi) < jetmetdphi_cut) continue;
	if(*met < met_cut) continue;
	if(fabs(*met-*metcalo)/(*met) > 0.5) continue;
	if(boostedJetpt->size() <= 0) continue;
	if(boostedJetpt->at(0)  < 250) continue;
	if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
	if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
	if(prunedJetm->at(0) < 65 or prunedJetm->at(0) > 105 ) continue;
	
	n_recolevel += *wgt*(*xsec)*luminosity*puwgt*trigeff/(sumwgt);
	
      }
      ///////////
      else if(category == Category::VBF or category == Category::VBFrelaxed){
	
	if(jetpt->size() < 2) continue;
	if(*nincjets < 2) continue;
	TLorentzVector leadingJet;
	leadingJet.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	TLorentzVector subleadingJet;
	subleadingJet.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));      
	if(leadingJet.Pt() < jetpt1_cut) continue;
	if(fabs(leadingJet.Eta()) > jeteta1_cut) continue;
	if(subleadingJet.Pt() < jetpt2_cut) continue;
	if(fabs(subleadingJet.Eta()) > jeteta2_cut) continue;
	
	if(fabs(leadingJet.Eta()) < 2.4 and chfrac->at(0) < 0.1) continue;
	if(fabs(leadingJet.Eta()) < 2.4 and nhfrac->at(0) > 0.8) continue;
	if(*met < met_cut) continue;
	if(fabs(*jmmdphi) < jetmetdphi_cut) continue;
	if(fabs(*met-*metcalo)/(*met) > 0.5) continue;
	if(leadingJet.Eta()*subleadingJet.Eta() > 0) continue;

	if(fabs(leadingJet.Eta()-subleadingJet.Eta()) < detajj_cut) continue;
	if((leadingJet+subleadingJet).M() < mjj_cut) continue; 
	if(fabs(leadingJet.DeltaPhi(subleadingJet)) > dphijj_cut) continue;
	
	n_recolevel += *wgt*(*xsec)*luminosity*puwgt*trigeff/(sumwgt);
      }
    }    
  }  

  if(not isGen)
    cout<<"ntotal "<<n_total<<" selected events "<<n_recolevel<<" efficiency "<<double(n_recolevel)/(double(n_total))<<endl;
  else
    cout<<"ntotal "<<n_total<<" selected events "<<n_genlevel<<" efficiency "<<double(n_genlevel)/(double(n_total))<<endl;
}
