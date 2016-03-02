#include <string>
#include <vector>
#include <sstream>

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLorentzVector.h"

void makeSmallGenTree(string interaction, string signalType, string outputDirectory){

  system(("mkdir -p "+outputDirectory).c_str());

  TChain* chain = new TChain("tree/tree");
  
  if(interaction == "Vector" and signalType == "MonoJ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/DMV_Vector/*root");
  else if(interaction == "Axial" and signalType == "MonoJ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/DMV_Axial/*root");
  else if(interaction == "Scalar" and signalType == "MonoJ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/DMS_Scalar/*root");
  else if(interaction == "Pseudoscalar" and signalType == "MonoJ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/DMS_Pseudoscalar/*root");

  if(interaction == "Vector" and signalType == "MonoW")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoW_Vector/*root");
  else if(interaction == "Axial" and signalType == "MonoW")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoW_Axial/*root");
  else if(interaction == "Scalar" and signalType == "MonoW")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoW_Scalar/*root");
  else if(interaction == "Pseudoscalar" and signalType == "MonoW")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoW_Pseudoscalar/*root");

  if(interaction == "Vector" and signalType == "MonoZ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoZ_Vector/*root");
  else if(interaction == "Axial" and signalType == "MonoZ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoZ_Axial/*root");
  else if(interaction == "Scalar" and signalType == "MonoZ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoZ_Scalar/*root");
  else if(interaction == "Pseudoscalar" and signalType == "MonoZ")
    chain->Add("/home/rgerosa/MONOJET_ANALYSIS/InterpolationFiles/MonoZ_Pseudoscalar/*root");

  // load all the re-weight files
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");

  // trigger efficiency for met trigger                                                                                                                                       
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
  TH1*   trmhist = (TH1*) trmfile->Get("mettrigSF");


  TFile* outputFile = new TFile((outputDirectory+"/tree_"+interaction+"_"+signalType+".root").c_str(),"RECREATE");
  outputFile->cd();

  TTree* outputTree = new TTree("tree","tree");

  int id;
  double  genVBosonPt, genVBosonEta, genVBosonPhi, genVBosonMass;
  double  genAK8JetPt, genAK8JetEta, genAK8JetPhi, genAK8JetMass, genAK8JetPrunedMass, genAK8JetTau2Tau1;
  double  genAK4JetPt, genAK4JetEta, genAK4JetPhi, genAK4JetMass;
  double  genMediatorPt, genMediatorEta, genMediatorPhi, genMediatorMass;
  double  genX1Pt, genX1Eta, genX1Phi, genX1Mass;
  double  genX2Pt, genX2Eta, genX2Phi, genX2Mass;
  double  genMetPt, genMetPhi;
  double  weight;
  double  pfMetPt, pfMetPhi;
  double  recoAK8JetPt,recoAK8JetPrunedMass, recoAK8JetTau2Tau1;

  outputTree->Branch("id",&id,"id/I");  
  outputTree->Branch("genVBosonPt", &genVBosonPt, "genVBosonPt/D");  
  outputTree->Branch("genVBosonEta", &genVBosonEta, "genVBosonEta/D");  
  outputTree->Branch("genVBosonPhi", &genVBosonPhi, "genVBosonPhi/D");  
  outputTree->Branch("genVBosonMass", &genVBosonMass, "genVBosonMass/D");  

  outputTree->Branch("recoAK8JetPt", &recoAK8JetPt, "recoAK8JetPt/D");  
  outputTree->Branch("recoAK8JetPrunedMass", &recoAK8JetPrunedMass, "recoAK8JetPrunedMass/D");  
  outputTree->Branch("recoAK8JetTau2Tau1", &recoAK8JetTau2Tau1, "recoAK8JetTau2Tau1/D");  

  outputTree->Branch("genAK8JetPt", &genAK8JetPt, "genAK8JetPt/D");  
  outputTree->Branch("genAK8JetEta", &genAK8JetEta, "genAK8JetEta/D");  
  outputTree->Branch("genAK8JetPhi", &genAK8JetPhi, "genAK8JetPhi/D");  
  outputTree->Branch("genAK8JetMass", &genAK8JetMass, "genAK8JetMass/D");  
  outputTree->Branch("genAK8JetPrunedMass", &genAK8JetPrunedMass, "genAK8JetPrunedMass/D");  
  outputTree->Branch("genAK8JetTau2Tau1", &genAK8JetTau2Tau1, "genAK8JetTau2Tau1/D");  
 
  outputTree->Branch("genAK4JetPt", &genAK4JetPt, "genAK4JetPt/D");  
  outputTree->Branch("genAK4JetEta", &genAK4JetEta, "genAK4JetEta/D");  
  outputTree->Branch("genAK4JetPhi", &genAK4JetPhi, "genAK4JetPhi/D");  
  outputTree->Branch("genAK4JetMass", &genAK4JetMass, "genAK4JetMass/D");  

  outputTree->Branch("genMediatorPt", &genMediatorPt, "genMediatorPt/D");  
  outputTree->Branch("genMediatorEta", &genMediatorEta, "genMediatorEta/D");  
  outputTree->Branch("genMediatorPhi", &genMediatorPhi, "genMediatorPhi/D");  
  outputTree->Branch("genMediatorMass", &genMediatorMass, "genMediatorMass/D");  

  outputTree->Branch("genX1Pt", &genX1Pt, "genX1Pt/D");  
  outputTree->Branch("genX1Eta", &genX1Eta, "genX1Eta/D");  
  outputTree->Branch("genX1Phi", &genX1Phi, "genX1Phi/D");  
  outputTree->Branch("genX1Mass", &genX1Mass, "genX1Mass/D");  

  outputTree->Branch("genX2Pt", &genX2Pt, "genX2Pt/D");  
  outputTree->Branch("genX2Eta", &genX2Eta, "genX2Eta/D");  
  outputTree->Branch("genX2Phi", &genX2Phi, "genX2Phi/D");  
  outputTree->Branch("genX2Mass", &genX2Mass, "genX2Mass/D");  

  outputTree->Branch("genMetPt", &genMetPt, "genMetPt/D");  
  outputTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/D");  
  outputTree->Branch("pfMetPt", &pfMetPt, "pfMetPt/D");  
  outputTree->Branch("pfMetPhi", &pfMetPhi, "pfMetPhi/D");  

  outputTree->Branch("weight", &weight, "weight/D");  
  

  // setup read input file                                                                                                                                                  
  TTreeReader myReader(chain);

  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");

  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbheloose");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nbjets (myReader,"nbjetslowpt");
  TTreeReaderValue<double> j1pt         (myReader,"leadingjetpt");
  TTreeReaderValue<vector<double> > jetpt   (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"centraljetbtag");
  TTreeReaderValue<vector<double> > jetm    (myReader,"centraljetm");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"centraljetNHfrac");  

  TTreeReaderValue<vector<double> > jetGenpt   (myReader,"centraljetGenpt");
  TTreeReaderValue<vector<double> > jetGeneta  (myReader,"centraljetGeneta");
  TTreeReaderValue<vector<double> > jetGenphi  (myReader,"centraljetGenphi");
  TTreeReaderValue<vector<double> > jetGenm    (myReader,"centraljetGenm");
  
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > boostedJetGenpt (myReader,"boostedJetGenpt");
  TTreeReaderValue<vector<double> > boostedJetGeneta   (myReader,"boostedJetGeneta");
  TTreeReaderValue<vector<double> > boostedJetGenphi   (myReader,"boostedJetGenphi");
  TTreeReaderValue<vector<double> > boostedJetGenm     (myReader,"boostedJetGenm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > prunedJetGenm      (myReader,"prunedJetGenm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<double> > boostedJetGentau2  (myReader,"boostedJetGentau2");
  TTreeReaderValue<vector<double> > boostedJetGentau1  (myReader,"boostedJetGentau1");
  TTreeReaderValue<double > bosonMass    (myReader,"wzmass_h");
  TTreeReaderValue<double > bosonPt      (myReader,"wzpt_h");
  TTreeReaderValue<double > bosonEta     (myReader,"wzeta_h");
  TTreeReaderValue<double > bosonPhi     (myReader,"wzphi_h");
  TTreeReaderValue<double > mediatorMass    (myReader,"dmmass");
  TTreeReaderValue<double > mediatorPt      (myReader,"dmpt");
  TTreeReaderValue<double > mediatorEta     (myReader,"dmeta");
  TTreeReaderValue<double > mediatorPhi     (myReader,"dmphi");
  TTreeReaderValue<double > x1Mass    (myReader,"dmX1mass");
  TTreeReaderValue<double > x1Pt      (myReader,"dmX1pt");
  TTreeReaderValue<double > x1Eta     (myReader,"dmX1eta");
  TTreeReaderValue<double > x1Phi     (myReader,"dmX1phi");
  TTreeReaderValue<double > x2Mass    (myReader,"dmX2mass");
  TTreeReaderValue<double > x2Pt      (myReader,"dmX2pt");
  TTreeReaderValue<double > x2Eta     (myReader,"dmX2eta");
  TTreeReaderValue<double > x2Phi     (myReader,"dmX2phi");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double> genmet        (myReader,"genmet");
  TTreeReaderValue<double> genmetphi     (myReader,"genmetphi");
  TTreeReaderValue<double> jmmdphi (myReader,"incjetmumetdphimin4");

  // loop on the event and apply selections
  cout<<"Events in the chain "<<chain->GetEntries()<<endl;
  TString fileName;
  vector<string> seglist;
  while(myReader.Next()){


    TString name_tmp(myReader.GetTree()->GetCurrentFile()->GetName());

    if(fileName != name_tmp){
      seglist.clear();
      fileName = name_tmp;
      fileName.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-madgraph.root","");
      fileName.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-powheg.root","");
      fileName.ReplaceAll("_gSM-1p0_gDM-1p0_13TeV-JHUGen.root","");
      
      stringstream name(fileName.Data());
      string segment;
    
      while(getline(name, segment, '-')){
	seglist.push_back(segment);
      }
    }

    // Set Branches
    id = -1;
    
    genVBosonPt   = 0.;
    genVBosonEta  = 0.;
    genVBosonPhi  = 0.;
    genVBosonMass = 0.;
    recoAK8JetPt   = 0.;
    genAK8JetPt   = 0.;
    genAK8JetEta  = 0.; 
    genAK8JetPhi  = 0.; 
    genAK8JetMass = 0.; 
    genAK8JetPrunedMass = 0.;
    recoAK8JetPrunedMass = 0.;
    genAK8JetTau2Tau1   = 0.;
    recoAK8JetTau2Tau1   = 0.;
    genAK4JetPt     = 0.; 
    genAK4JetEta    = 0.;
    genAK4JetPhi    = 0.; 
    genAK4JetMass   = 0.;
    genMediatorPt   = 0.;
    genMediatorEta  = 0.; 
    genMediatorPhi  = 0.;
    genMediatorMass = 0.;
    genX1Pt  = 0.;
    genX1Eta = 0.;
    genX1Phi = 0.;
    genX1Mass = 0.;
    genX2Pt  = 0.; 
    genX2Eta = 0.; 
    genX2Phi = 0.; 
    genX2Mass = 0.;
    genMetPt = 0.; 
    genMetPhi = 0.;
    pfMetPt = 0.;
    pfMetPhi = 0.;
    weight = 1.;

    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) id = 0;
    if (*fhbhe  == 0 or *fhbiso == 0) id = 0;
    //    if (*fcsc   == 0 or not *feeb) id = 0;
    if (*njets  < 1) id = 0;
    if (*nbjets > 0) id = 0;
    if (chfrac->size() == 0 or nhfrac->size() == 0 or jetpt->size() == 0) id = 0;                                                                
    if (chfrac->size() > 0 and chfrac->at(0) < 0.1)   id = 0;
    if (nhfrac->size() > 0 and nhfrac->at(0) > 0.8)   id = 0;
    if (jetpt->size()  > 0 and jetpt->at(0)  < 100.)  id = 0;
    if (jetpt->size()  > 0 and jetpt->at(0)  < *j1pt) id = 0;
    if (*jmmdphi < 0.5) id = 0;
    if (*mmet < 200)    id = 0;    

    if(id == -1 and boostedJetpt->size()  == 0) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and fabs(boostedJeteta->at(0)) > 2.4) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and boostedJetpt->at(0) < 250.) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) id = 1;
    if(id == -1 and boostedJetpt->size()  > 0 and (prunedJetm->at(0) < 65 or prunedJetm->at(0) > 105)) id = 1;

    if(id == -1 and prunedJetm->size() > 0 and prunedJetm->at(0) > 0 and boostedJetpt->at(0) > 250. and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105 and *mmet > 250 and fabs(boostedJeteta->at(0)) < 2.4) id = 2;
    
    if (id == -1 and *mmet > 200 and *mmet < 250) id = 0;

    if(id == -1) cout<<"Huston we have a problem ... "<<*hltm90<<" "<<*fhbhe<<" "<<*fhbiso<<" "<<*feeb<<" "<<*fcsc<<" "<<*njets<<" "<<*nbjets<<" "<<chfrac->at(0)<<" "<<nhfrac->at(0) << " "<<jetpt->at(0)<<" "<<*jmmdphi<<" "<<*mmet<<" "<<boostedJetpt->at(0)<<" "<<boostedJettau2->at(0)/boostedJettau1->at(0)<<" "<<prunedJetm->at(0)<<endl;
      
    // fill branches
    genVBosonPt   = *bosonPt;  
    genVBosonPhi  = *bosonPhi;  
    genVBosonEta  = *bosonEta;  
    genVBosonMass = *bosonMass;  
    
    if(boostedJetGenpt->size() > 0){
      genAK8JetPt   = boostedJetGenpt->at(0);
      genAK8JetEta  = boostedJetGeneta->at(0);
      genAK8JetPhi  = boostedJetGenphi->at(0);
      genAK8JetMass = boostedJetGenm->at(0);
      genAK8JetPrunedMass = prunedJetGenm->at(0);
      genAK8JetTau2Tau1   = boostedJetGentau2->at(0)/boostedJetGentau1->at(0);
    }
    
    if(jetGenpt->size() > 0){
      genAK4JetPt   = jetGenpt->at(0);
      genAK4JetEta  = jetGeneta->at(0);
      genAK4JetPhi  = jetGenphi->at(0);
      genAK4JetMass = jetGenm->at(0);
    }
    
    genX1Pt   = *x1Pt;
    genX1Eta  = *x1Eta;
    genX1Phi  = *x1Phi;
    genX1Mass = stod(seglist.back());
    genX2Pt   = *x2Pt;
    genX2Eta  = *x2Eta;
    genX2Phi  = *x2Phi;
    genX2Mass = stod(seglist.back());

    TLorentzVector dmX1, dmX2, dmMED;
    dmX1.SetPtEtaPhiM(genX1Pt,genX1Eta,genX1Phi,genX1Mass);
    dmX2.SetPtEtaPhiM(genX2Pt,genX2Eta,genX2Phi,genX2Mass);
    dmMED = dmX1+dmX2;

    genMediatorPt   = *mediatorPt;
    genMediatorEta  = *mediatorEta;
    genMediatorPhi  = *mediatorPhi;
    genMediatorMass = stod(seglist.at(seglist.size()-2));
    
    if(*nvtx<= 35)
      weight *= trmhist->GetBinContent(trmhist->FindBin(*mmet))*puhist->GetBinContent(*nvtx);
    else
      weight *= trmhist->GetBinContent(trmhist->FindBin(*mmet));


    pfMetPt = *mmet;
    pfMetPhi = *mmetphi;

    genMetPt  = *genmet;
    genMetPhi = *genmetphi;

    if(boostedJetpt->size() > 0){
      recoAK8JetPt = boostedJetpt->at(0);
      recoAK8JetPrunedMass = prunedJetm->at(0);
      recoAK8JetTau2Tau1 = boostedJettau2->at(0)/boostedJettau1->at(0);
    }

    outputTree->Fill();

  }

  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}
