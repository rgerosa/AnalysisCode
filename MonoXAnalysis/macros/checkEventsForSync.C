#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReader.h"

void checkEventsForSync(string inputFile, string outputDir){

  system(("mkdir -p "+outputDir).c_str());
  system(("rm "+outputDir+"/*").c_str());

  TFile* input = TFile::Open(inputFile.c_str());
  TTree* tree  = (TTree*) input->Get("tree/tree");
  
  // declare branches
  TTreeReader myReader(tree);
  TTreeReaderValue<unsigned int> run    (myReader,"run");
  TTreeReaderValue<unsigned int> lumi   (myReader,"lumi");
  TTreeReaderValue<unsigned int> event  (myReader,"event");

  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte       (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltp165    (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (myReader,"hltphoton175");

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");
  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nbjets (myReader,"nbjets");
  TTreeReaderValue<unsigned int> nbjetslowpt (myReader,"nbjetslowpt");
  
  TTreeReaderValue<vector<double> > jetpt   (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"centraljetbtag");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"centraljetNHfrac");
  TTreeReaderValue<vector<double> > emfrac  (myReader,"centraljetEMfrac");

  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > prunedJetm_v2   (myReader,"prunedJetm_v2");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");

  TTreeReaderValue<double> met (myReader,"t1pfmet");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> emet        (myReader,"t1elmet");
  TTreeReaderValue<double> pmet        (myReader,"t1phmet");
  TTreeReaderValue<double> jmmdphi (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi (myReader,"incjetelmetdphimin4");
  TTreeReaderValue<double> jpmdphi (myReader,"incjetphmetdphimin4");

  TTreeReaderValue<unsigned int> nmuons     (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (myReader,"ntaus");
  TTreeReaderValue<unsigned int> nphotons   (myReader,"nphotons");

  // output text files for the comparison
  ofstream leptonVeto((outputDir+"/leptonVeto_SR.txt").c_str());
  ofstream leptonPhtonVeto((outputDir+"/leptonPhtonVeto_SR.txt").c_str());
  ofstream leptonPhtonTauVeto((outputDir+"/leptonPhtonTauVeto_SR.txt").c_str());
  ofstream ak4JetSelections((outputDir+"/ak4JetSelections_SR.txt").c_str());
  ofstream metSelections((outputDir+"/metSelections_SR.txt").c_str());
  ofstream btagVetoSelections((outputDir+"/btagVeto_SR.txt").c_str());
  ofstream VtaggingSelections((outputDir+"/Vtagging_SR.txt").c_str());
  
  /// event loop
  while(myReader.Next()){

    if(*nmuons != 0) continue;
    if(*nelectrons != 0) continue;
    leptonVeto << *run << " "<<*lumi<<" "<<*event<<"\n";
    if(*nphotons != 0) continue;    
    leptonPhtonVeto << *run << " "<<*lumi<<" "<<*event<<"\n";
    if(*ntaus != 0) continue;
    leptonPhtonTauVeto << *run << " "<<*lumi<<" "<<*event<<"\n";
        
    if(*nbjetslowpt > 0) continue;
    btagVetoSelections << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;


    ak4JetSelections << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(*jmmdphi < 0.5) continue;
    if(*met < 200) continue;

    metSelections << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;

    VtaggingSelections << *run << " "<<*lumi<<" "<<*event<<"\n";
   

  }

  leptonVeto.close();
  leptonPhtonVeto.close();
  leptonPhtonTauVeto.close();
  ak4JetSelections.close();
  metSelections.close();
  btagVetoSelections.close();
  VtaggingSelections.close();

}
