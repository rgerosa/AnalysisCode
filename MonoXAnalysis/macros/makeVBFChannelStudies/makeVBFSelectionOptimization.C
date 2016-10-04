#include "../CMS_lumi.h"

static float leadingJetPt  = 50;
static float trailingJetPt = 50;
static float missingEnergy = 150;
static float jetMetDPhi    = 0.5;
static float deltaEtaJJ    = 2;
static float Mjj           = 400;
static int   reduceStatistics = 2;
////////////////////// 
float mjj, detajj, jetpt1, jeteta1, jeteta2, jetpt2, dphijj, jetmetdphi, mTjj, njet, weight;
//////////////////////
void makeSimpleTreeForTraining(TTree* outputTree, vector<TTree*> chain, vector<TH1*> khists, const float & luminosity, const bool & isKfactor = false, const bool & reduceStat = false){

  // pileup re-weight
  TFile* pileupFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/pileupWeight/puweight_12p9fb.root");
  TH1*   puhist = (TH1*) pileupFile->Get("puhist");
  // trigger MET
  TFile* triggerfile_MET   =  TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency_12p9.root");
  TEfficiency*  triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();

  // optput tree
  outputTree->Branch("mjj",&mjj,"mjj/F");
  outputTree->Branch("detajj",&detajj,"detajj/F");
  outputTree->Branch("jetpt1",&jetpt1,"jetpt1/F");
  outputTree->Branch("jeteta1",&jeteta1,"jeteta1/F");
  outputTree->Branch("jetpt2",&jetpt2,"jetpt2/F");
  outputTree->Branch("jeteta2",&jeteta2,"jeteta2/F");
  outputTree->Branch("dphijj",&dphijj,"dphijj/F");
  outputTree->Branch("jetmetdphi",&jetmetdphi,"jetmetdphi/F");
  outputTree->Branch("mTjj",&mTjj,"mTjj/F");
  outputTree->Branch("njet",&njet,"njet/F");
  outputTree->Branch("weight",&weight,"weight/F");  
  mjj = -1; detajj = -1; jetpt1 = -1; jeteta1 = -1; jeteta2 = -1; jetpt2 = -1; dphijj = -1; jetmetdphi = -1; mTjj= -1; njet = -1; weight = -1;

  for(auto tree : chain){ 

    TTreeReader myReader (tree);
    
    TTreeReaderValue<double> xsec         (myReader,"xsec");  
    // weights
    TTreeReaderValue<double> wgt          (myReader,"wgt");
    TTreeReaderValue<double> wgtsum       (myReader,"wgtsum");
    TTreeReaderValue<double> wgtbtag      (myReader,"wgtbtag");
    TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
  // met filters 
    TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
    TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
    TTreeReaderValue<UChar_t> fcsc   (myReader,"flagcsctight");
    TTreeReaderValue<UChar_t> fcsct  (myReader,"flagcsctight");
    TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
    TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
    TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");
    // Object count
    TTreeReaderValue<unsigned int> njetsinc  (myReader,"njetsinc");
    TTreeReaderValue<unsigned int> nphotons  (myReader,"nphotons");
    TTreeReaderValue<unsigned int> nelectrons   (myReader,"nelectrons");
    TTreeReaderValue<unsigned int> ntaus   (myReader,"ntausraw");
    TTreeReaderValue<unsigned int> nmuons (myReader,"nmuons");
    TTreeReaderValue<unsigned int> nbjets  (myReader,"nbjetslowpt");
    // jets
    TTreeReaderValue<vector<double> > jetpt    (myReader,"combinejetpt");
    TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
    TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
    TTreeReaderValue<vector<double> > jetbtag  (myReader,"combinejetbtag");
    TTreeReaderValue<vector<double> > jetm     (myReader,"combinejetm");
    TTreeReaderValue<vector<double> > chfrac   (myReader,"combinejetCHfrac");
    TTreeReaderValue<vector<double> > nhfrac   (myReader,"combinejetNHfrac");
    TTreeReaderValue<double> jmmdphi4    (myReader,"incjetmumetdphimin4");
    TTreeReaderValue<double> jmmdphi     (myReader,"incjetmumetdphimin");  
    // met information
    TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
    TTreeReaderValue<double> mmet        (myReader,"t1mumet");
    // gen level boson for background
    TTreeReaderValue<double> wzpt        (myReader,"wzpt");
    TTreeReaderValue<double> wzeta       (myReader,"wzeta");
    TTreeReaderValue<double> wzphi       (myReader,"wzphi");
    TTreeReaderValue<double> wzmass      (myReader,"wzmass");
    // triggers
    TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
    TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
    TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
    TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
    TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
    TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
    

    long int nEvent = 0;
    
    // Loop on events 
    while(myReader.Next()){
      
      if(nEvent > tree->GetEntries()/reduceStatistics and reduceStat) continue;
      nEvent++;

      mjj = -1; detajj = -1; jetpt1 = -1; jeteta1 = -1; jeteta2 = -1; jetpt2 = -1; dphijj = -1; jetmetdphi = -1; mTjj= -1; njet = -1; weight = -1;
      
      // met trigger requirement
      if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;    
      // met filters                                                                                                                                                                                  
      if (*fhbhe == 0 || *fhbiso == 0 || *feeb == 0 || *fetp == 0 || *fvtx == 0 || *fcsc == 0 || *fcsct == 0) continue;
      // number of jets                                                                                                                                                                               
      if (*njetsinc < 2) continue;
      // b-veto                                                                                                                                                                                       
      if (*nbjets   > 0) continue;
      // vetos                                                                                                                                                                                       
      if (*nmuons > 0)     continue;
      if (*nelectrons > 0) continue;
      if (*ntaus > 0)      continue;
      if (*nphotons > 0)   continue;
      // relaxed vbf jet pt cuts                                                                                                                                                                     
      if (jetpt->at(0) < leadingJetPt) continue;
      if (jetpt->at(1) < trailingJetPt) continue;
      if (fabs(jeteta->at(0)) > 4.7) continue;
      if (fabs(jeteta->at(1)) > 4.7) continue;
      if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
      if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) > 0.8) continue;
      
      // relaxed met cut                                                                                                                                                                               
      if (*mmet     < missingEnergy) continue;
      if (*jmmdphi4 < jetMetDPhi) continue;
      
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      // relaxed VBF cuts                                                                                                                                                                             
      if (fabs(jeteta->at(0)-jeteta->at(1)) < deltaEtaJJ) continue;
      if ((jet1+jet2).M() < Mjj) continue;
      
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
      
      //filling variables
      mjj = (jet1+jet2).M(); 
      detajj = fabs(jeteta->at(0)-jeteta->at(1)); 
      jetpt1 = jetpt->at(0); 
      jeteta1 = jeteta->at(0); 
      jeteta2 = jeteta->at(1); 
      jetpt2 = jetpt->at(1); 
      dphijj = fabs(jet1.DeltaPhi(jet2)); 
      jetmetdphi = *jmmdphi4; 
      mTjj   = sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))); 
      njet   = *njetsinc; 
      if(not reduceStat)
	weight = *xsec*luminosity*(*wgt)*kwgt*puwgt*trgwgt/(*wgtsum);
      else
	weight = *xsec*luminosity*(*wgt)*kwgt*puwgt*trgwgt/(*wgtsum/reduceStatistics);

      
      outputTree->Fill();
      
    }
  }
  
  pileupFile->Close();
  triggerfile_MET->Close();
}


////////////////////
// on single variable do rectangular cut optimization
// on pairs variable do rectangular cut optimization and Likelihood
// on combo do Likelihood, MLP and BDT
void makeVBFSelectionOptimization(float luminosity, string outputDirectory, string outputName, bool doSingleVariables, bool doVariablesPair, bool doCombo){

  //  vector<TString> inputVariablesSingle = {"mjj","detajj","dphijj","mTjj","jetmetdphi","jetpt1","jetpt2"};
  vector<TString> inputVariablesSingle = {"jetmetdphi","jetpt1","jetpt2"};
  vector<TString> inputVariablesPair   = {"mjj-detajj","mjj-dphijj","mjj-mTjj","mjj-njet","mjj-jetmetdphi","detajj-dphijj","detajj-njet","detajj-mTjj","detajj-jetmetdphi"};
  vector<TString> inputVariablesComb   = {"mjj","detajj","dphijj","mTjj","jetmetdphi","jetpt1","jetpt2","njet"};

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());

  // take trees for VBF signal
  vector<TFile*> vbftreeFile;
  vector<TTree*> vbftree;
  vbftreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M110_13TeV_powheg_pythia8.root"));
  vbftree.push_back((TTree*) vbftreeFile.back()->Get("tree/tree"));
  vbftreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M125_13TeV_powheg_pythia8.root"));
  vbftree.push_back((TTree*) vbftreeFile.back()->Get("tree/tree"));
  vbftreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M150_13TeV_powheg_pythia8.root"));
  vbftree.push_back((TTree*) vbftreeFile.back()->Get("tree/tree"));
  vbftreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M200_13TeV_powheg_pythia8.root"));
  vbftree.push_back((TTree*) vbftreeFile.back()->Get("tree/tree"));

  // Znn backgrounds : QCD and EWK
  vector<TFile*>  znntreeFile ;
  vector<TTree*>  znntree ;
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-100To200_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-200To400_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-400To600_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-600To800_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-800To1200_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-1200To2500_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));
  znntreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph.root"));
  znntree.push_back((TTree*) znntreeFile.back()->Get("tree/tree"));

  vector<TFile*> znnewktreeFile;
  znnewktreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJetsToNuNuEWK/sigfilter/sig_EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8.root"));
  vector<TTree*> znnewktree;
  znnewktree.push_back((TTree*) znnewktreeFile.back()->Get("tree/tree"));
  
  // WJet backgrounds : QCD and EWK
  vector<TFile*> wjettreeFile;
  vector<TTree*> wjettree;

  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));
  wjettreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJets/sigfilter/sig_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"));
  wjettree.push_back((TTree*) wjettreeFile.back()->Get("tree/tree"));

  vector<TFile*> wjetewktreeFile;
  vector<TTree*> wjetewktree;
  wjetewktreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJetsEWK/sigfilter/sig_EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8.root"));
  wjetewktree.push_back((TTree*) wjetewktreeFile.back()->Get("tree/tree"));
  wjetewktreeFile.push_back(TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJetsEWK/sigfilter/sig_EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8.root"));
  wjetewktree.push_back((TTree*) wjetewktreeFile.back()->Get("tree/tree"));
  
  // get k-factors NLO                                                                                                                                                                                
  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1*  zewkhist = (TH1*) kffile.Get("EWKcorr/Z");
  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1*  wewkhist = (TH1*) kffile.Get("EWKcorr/W");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  vector<TH1*> khists;
  khists.push_back(znlohist);
  khists.push_back(zewkhist);

  // make small trees for traning
  TFile* outputTemp = new TFile("/tmp/rgerosa/outputTemp.root","RECREATE");
  outputTemp->cd();
  
  TTree* vbftree_training = new TTree("vbftree","vbftree");
  TTree* znntree_training = new TTree("znntree","znntree");
  TTree* znnewktree_training = new TTree("znnewktree","znnewktree");
  TTree* wjettree_training = new TTree("wjettree","wjettree");
  TTree* wjetewktree_training = new TTree("wjetewktree","wjetewktree");

  cout<<"Fill VBF tree for training/testing "<<endl;
  makeSimpleTreeForTraining(vbftree_training,vbftree,khists,luminosity,false,false);
  cout<<"Fill Znn QCD tree for training/testing "<<endl;
  if(not doCombo)
    makeSimpleTreeForTraining(znntree_training,znntree,khists,luminosity,true,true);  
  else
    makeSimpleTreeForTraining(znntree_training,znntree,khists,luminosity,true,false);  
  cout<<"Fill Znn EWK tree for training/testing "<<endl;
  makeSimpleTreeForTraining(znnewktree_training,znnewktree,khists,luminosity,false,false);
  
  khists.clear();
  khists.push_back(wnlohist);
  khists.push_back(wewkhist);

  cout<<"Fill W+jets QCD tree for training/testing "<<endl;
  if(not doCombo)
    makeSimpleTreeForTraining(wjettree_training,wjettree,khists,luminosity,true,true);
  else
    makeSimpleTreeForTraining(wjettree_training,wjettree,khists,luminosity,true,false);
  cout<<"Fill W+jets EWK tree for training/testing "<<endl;
  makeSimpleTreeForTraining(wjetewktree_training,wjetewktree,khists,luminosity,false,false);
  

  // setup the MVA factory
  TMVA::Tools::Instance();
  
  //output file
  //
  if(doSingleVariables){
    for(auto var :inputVariablesSingle){
      // create a factory
      TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_"+string(var)+"_cuts.root").c_str(),"RECREATE"); 
      outputFile->cd();
      TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification:Transformations=I");
      system(("mkdir -p "+outputDirectory+"/WeightCut_"+string(var)).c_str());
      (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory+"/WeightCut_"+var;
      factory->AddVariable(var,'F');
      factory->AddSignalTree(vbftree_training,1.0);
      factory->AddBackgroundTree(znntree_training,1.0);
      factory->AddBackgroundTree(znnewktree_training,1.0);
      factory->AddBackgroundTree(wjettree_training,1.0);
      factory->AddBackgroundTree(wjetewktree_training,1.0);      
      factory->SetSignalWeightExpression("weight");
      factory->SetBackgroundWeightExpression("weight");
      factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
      factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:VarProp=FSmart" );
      
      factory->TrainAllMethods();
      factory->TestAllMethods();
      factory->EvaluateAllMethods();
      
      outputFile->Close();
      delete factory;
    }
  }
  
  //////////////////
  if(doVariablesPair){
    for(auto var :inputVariablesPair){
      TObjArray* array = var.Tokenize("-");
      if(array->GetEntries() != 2){
	cout<<"Problem in dividing the string "<<var<<endl;
	continue;
      }

      TString var1 (((TObjString *) array->At(0))->String());
      TString var2 (((TObjString *) array->At(1))->String());

      TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_"+string(var1)+"_"+string(var2)+".root").c_str(),"RECREATE");
      outputFile->cd();
      TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification:Transformations=I");
      system(("mkdir -p "+outputDirectory+"/WeightCut_"+string(var1)+"_"+string(var2)).c_str());
      (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory+"/WeightCut_"+var1+"_"+var2;
      factory->AddVariable(var1,'F');
      factory->AddVariable(var2,'F');
      factory->AddSignalTree(vbftree_training,1.0);
      factory->AddBackgroundTree(znntree_training,1.0);
      factory->AddBackgroundTree(znnewktree_training,1.0);
      factory->AddBackgroundTree(wjettree_training,1.0);
      factory->AddBackgroundTree(wjetewktree_training,1.0);
      factory->SetSignalWeightExpression("weight");
      factory->SetBackgroundWeightExpression("weight");      
      factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
      factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:VarProp=FSmart" );
      factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V:!TransformOutput:PDFInterpol=Spline3:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50:VarTransform=D");
      factory->TrainAllMethods();
      factory->TestAllMethods();
      factory->EvaluateAllMethods();
      
      outputFile->Close();
      delete factory;
    }
  }

  if(doCombo){

    TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_combo.root").c_str(),"RECREATE");
      outputFile->cd();
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification:Transformations=I");
    (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory;
    
    for(auto var : inputVariablesComb)
      factory->AddVariable(var,'F');
    
    factory->AddSignalTree(vbftree_training,1.0);
    factory->AddBackgroundTree(znntree_training,1.0);
    factory->AddBackgroundTree(znnewktree_training,1.0);
    factory->AddBackgroundTree(wjettree_training,1.0);
    factory->AddBackgroundTree(wjetewktree_training,1.0);
    factory->SetSignalWeightExpression("weight");
    factory->SetBackgroundWeightExpression("weight");
    factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
    factory->BookMethod(TMVA::Types::kBDT, "BDTG","!H:!V:CreateMVAPdfs:NTrees=1000:BoostType=Grad:!UseBaggedGrad:GradBaggingFraction=0.5:PruneMethod=NoPruning:PruneStrength=5:MaxDepth=5:SeparationType=GiniIndex:Shrinkage=0.1:nCuts=200:VarTransform=D,P,G");
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();
    delete factory;
  }
    
  outputTemp->Close();
  system("rm -r /tmp/rgerosa/outputTemp.root");

}
