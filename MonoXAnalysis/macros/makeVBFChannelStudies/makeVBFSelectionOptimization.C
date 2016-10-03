#include "../CMS_lumi.h"

static float leadingJetPt  = 50;
static float trailingJetPt = 50;
static float missingEnergy = 150;
static float jetMetDPhi    = 0.5;
static float deltaEtaJJ    = 2;
static float Mjj           = 400;
static 
////////////////////// 
TTree* makeSimpleTreeForTraining(TChain* chain, vector<TH1*> khists, const float & luminosity, const bool & isKfactor = false){

  // pileup re-weight
  TFile* pileupFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/pileupWeight/puweight_12p9fb.root");
  TH1*   puhist = (TH1*) pileupFile->Get("puhist");
  // trigger MET
  TFile* triggerfile_MET   =  TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency_12p9.root");
  TEfficiency*  triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();
  
  TTreeReader myReader (chain);

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

  // optput tree
  TTree* outputTree = new TTree("tree","");
  float mjj, detajj, jetpt1, jeteta1, jeteta2, jetpt2, dphijj, jetmetdphi, mTjj, njet, weight;
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

  // Loop on events 
  while(myReader.Next()){

    mjj = -1; detajj = -1; jetpt1 = -1; jeteta1 = -1; jeteta2 = -1; jetpt2 = -1; dphijj = -1; jetmetdphi = -1; mTjj= -1; njet = -1; weight = -1;

    // met trigger requirement
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;    
    // met filters                                                                                                                                                                                    
    if (*fhbhe == 0 || *fhbiso == 0 || *feeb == 0 || *fetp == 0 || *fvtx == 0 || *fcsc == 0 || *fcsct == 0) continue;
    // number of jets                                                                                                                                                                                 
    if (*njetsinc  < 2) continue;
    // b-veto                                                                                                                                                                                         
    if (*nbjets    > 0) continue;
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
    mTjj   = sqrt(2*jetpt->at(0)*jetpt->at(1)*(1-cos(fabs(jet1.DeltaPhi(jet2))))); 
    njet   = *njetsinc; 
    weight = *xsec*luminosity*(*wgt)*kwgt*puwgt*trgwgt/(*wgtsum);
    
    outputTree->Fill();

  }

  pileupFile->Close();
  triggerfile_MET->Close();

  return outputTree;
}


////////////////////
// on single variable do rectangular cut optimization
// on pairs variable do rectangular cut optimization and Likelihood
// on combo do Likelihood, MLP and BDT
void makeVBFSelectionOptimization(float luminosity, string outputDirectory, string outputName, bool doSingleVariables, bool doVariablesPair, bool doCombo){

  vector<TString> inputVariablesSingle = {"mjj","detajj","dphijj","mTjj","njet","jetmetdphi","ptj1","ptj2"};
  vector<TString> inputVariablesPair   = {"mjj-detajj","mjj-dphijj","mjj-mTjj","mjj-njet","mjj-jetmetdphi","detajj-dphijj","detajj-njet","detajj-mTjj","detajj-jmetmetdphi"};

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());

  // take trees for VBF signal
  TChain* vbftree = new TChain("tree/tree");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*110*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*125*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M*150*root");

  // Znn backgrounds : QCD and EWK
  TChain* znntree = new TChain("tree/tree");
  znntree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJets/sigfilter/sig_ZJetsToNuNu_HT-*root");
  TChain* znnewktree = new TChain("tree/tree");
  znnewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/ZJetsToNuNuEWK/sigfilter/sig_EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8.root");

  // WJet backgrounds : QCD and EWK
  TChain* wjettree = new TChain("tree/tree");
  wjettree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJetsNLO/sigfilter/*root");
  TChain* wjetewktree = new TChain("tree/tree");
  wjetewktree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_16_07_2016/WJetsEWK/sigfilter/*root");
  
  // get k-factors NLO                                                                                                                                                                                
  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1*  zewkhist = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  vector<TH1*> khists;
  khists.push_back(znlohist);
  khists.push_back(zewkhist);

  // make small trees for traning
  TTree* vbftree_training = makeSimpleTreeForTraining(vbftree,khists,luminosity,false);
  TTree* znntree_training = makeSimpleTreeForTraining(znntree,khists,luminosity,true);
  TTree* znnewktree_training = makeSimpleTreeForTraining(znnewktree,khists,luminosity,true);
  TTree* wjettree_training = makeSimpleTreeForTraining(wjettree,khists,luminosity,true);
  TTree* wjetewktree_training = makeSimpleTreeForTraining(wjetewktree,khists,luminosity,true);

  // setup the MVA factory
  TMVA::Tools::Instance();

  //output file
  //
  if(doSingleVariables){
    for(auto var :inputVariablesSingle){
      // create a factory
      TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_"+string(var)+"_cuts.root").c_str(),"RECREATE"); 
      TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification:Transformations=I");
      (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory;
      factory->AddVariable(var,'F');
      factory->AddSignalTree(vbftree_training,1.0);
      factory->AddBackgroundTree(znntree_training,1.0);
      factory->AddBackgroundTree(znnewktree_training,1.0);
      factory->AddBackgroundTree(wjettree_training,1.0);
      factory->AddBackgroundTree(wjetewktree_training,1.0);
      factory->SetSignalWeightExpression("weight");
      factory->SetBackgroundWeightExpression("weight");
      factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
      factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
      //factory->BookMethod( TMVA::Types::kCuts, "Cuts", "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
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
      if(array->GetSize() != 2){
	cout<<"Problem in dividing the string "<<var<<endl;
	continue;
      }

      TString var1 (((TObjString *) array->At(0))->GetString());
      TString var2 (((TObjString *) array->At(1))->GetString());

      TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_"+string(var1)+"_"+string(var2)+".root").c_str(),"RECREATE");
      TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification,Transformations=I");
      (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory;
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
      factory->BookMethod( TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
      //factory->BookMethod( TMVA::Types::kCuts, "Cuts", "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
      factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V:!TransformOutput:PDFInterpol=Spline3:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:CreateMVAPdfs:VarTransform=D;P;G;N");

      factory->TrainAllMethods();
      factory->TestAllMethods();
      factory->EvaluateAllMethods();
      outputFile->Close();
      delete factory;
    }
  }
  
  if(doCombo){

    TFile* outputFile = new TFile((outputDirectory+"/"+outputName+"_combo.root").c_str(),"RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification,Transformations=I");
    (TMVA::gConfig().GetIONames()).fWeightFileDir = outputDirectory;
    
    for(auto var : inputVariablesSingle)
      factory->AddVariable(var,'F');
    
    factory->AddSignalTree(vbftree_training,1.0);
    factory->AddBackgroundTree(znntree_training,1.0);
    factory->AddBackgroundTree(znnewktree_training,1.0);
    factory->AddBackgroundTree(wjettree_training,1.0);
    factory->AddBackgroundTree(wjetewktree_training,1.0);
    factory->SetSignalWeightExpression("weight");
    factory->SetBackgroundWeightExpression("weight");
    factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
        
    factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood","!H:!V:!TransformOutput:PDFInterpol=Spline3:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:CreateMVAPdfs:VarTransform=D;P;G;N");
    
    factory->BookMethod(TMVA::Types::kBDT, "BDTG","!H:!V:CreateMVAPdfs:NTrees=500:BoostType=Grad:!UseBaggedGrad:GradBaggingFraction=0.5:PruneMethod=NoPruning:PruneStrength=5:MaxDepth=5:SeparationType=GiniIndex:Shrinkage=0.1:NNodesMax=100000:UseYesNoLeaf=F:nCuts=2000:IgnoreNegWeightsInTraining:VarTransform=D;P;G;N");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outputFile->Close();
    delete factory;
  }
  
}
