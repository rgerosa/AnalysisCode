#include "makehist.h"
#include "TChain.h"
#include "TF1.h"

using namespace std;

// make histograms for Z->mumu to signal region correction                                                                                                                   
void makezmmcorhist( const string &   signalRegionFile, 
		     const string &   zmumuFile,  
		     const Category   & category, 
		     const SamplesNLO & nloSamples,
		     vector<string>   observables, 
		     vector<string>   observables_2D, 
		     const double &   lumi, 
		     const string &   outDir  = "", 
		     const string &   sysName = "", 
		     const bool &     isHiggsInvisible = false,
		     const bool &     isEWK = false,
		     const string &   ext = "") {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((zmumuFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zmm_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zmm_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zmm_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zmm_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }
  
  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zmm_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zmm_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zmm_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zmm_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }


  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");    

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  
  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> dyhists;
  
  if(nloSamples.useZJetsNLO){
    zhists.push_back(zewkhist);
  }
  else{
    zhists.push_back(zewkhist);
    zhists.push_back(znlohist);
  }

  if(nloSamples.useDYJetsNLO){
    dyhists.push_back(zewkhist);
  }
  else{
    dyhists.push_back(zewkhist);
    dyhists.push_back(znlohist);
  }


  TFile* kfactzjet_vbf = NULL;
  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                               
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets.root");
    TH1* zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* zjet_nlo_mj  = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");
    zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    zjet_nlo_mj->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);
    if(not nloSamples.useZJetsNLO)
      zhists.push_back(zjet_nlo_vbf);
    if(not nloSamples.useDYJetsNLO)
      dyhists.push_back(zjet_nlo_vbf);    
  }


  if(not isEWK){
    // NLO Znunu or LO
    if(nloSamples.useZJetsNLO)
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible); 
    else
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);     
    makehist4(dtree, dhist, dhist_2D,  true, Sample::zmm, category, false, 1.00, lumi, dyhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);  
    makehist4(dtree, dhist, dhist_2D,  true, Sample::zmm, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("zmmcor")+ext;
  if(isEWK)
    name = string("zewkmmcor")+ext;

  // Make 1D transfer factor
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2); // smooth numerator in case
      smoothEmptyBins(dhist.at(ihist),2); // smooth denominator in case
    }
  }
  
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  // Make 2D transfer factor
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }    
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){    
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2); // smooth ratio
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }  
 
  outfile.cd();
  outfile.Close();
  kffile.Close();
  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();

  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();

  cout << "Z(mumu)->Z(inv) transfer factor computed ..." << endl;
}

// make histograms for Z->ee to signal region correction                                                                                                                   
void makezeecorhist( const string &   signalRegionFile, 
		     const string &   zeeFile,  
		     const Category   & category, 
		     const SamplesNLO & nloSamples,
		     vector<string>   observables, 
		     vector<string>   observables_2D, 
		     const double &   lumi, 
		     const string &   outDir  = "", 
		     const string &   sysName = "", 
		     const bool &     isHiggsInvisible = false,
		     const bool &     isEWK = false,
		     const string &   ext = "") {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((zeeFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zee_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zee_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zee_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zee_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }
  
  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zee_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zee_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zee_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zee_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }


  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");    

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  
  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> dyhists;
  
  if(nloSamples.useZJetsNLO){
    zhists.push_back(zewkhist);
  }
  else{
    zhists.push_back(zewkhist);
    zhists.push_back(znlohist);
  }

  if(nloSamples.useDYJetsNLO){
    dyhists.push_back(zewkhist);
  }
  else{
    dyhists.push_back(zewkhist);
    dyhists.push_back(znlohist);
  }

  TFile* kfactzjet_vbf = NULL;
  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                               
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets.root");
    TH1* zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* zjet_nlo_mj  = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");
    zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    zjet_nlo_mj->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);
    if(not nloSamples.useZJetsNLO)
      zhists.push_back(zjet_nlo_vbf);
    if(not nloSamples.useDYJetsNLO)
      dyhists.push_back(zjet_nlo_vbf);
  }


  if(not isEWK){
    // NLO Znunu or LO
    if(nloSamples.useZJetsNLO)
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible); 
    else
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);     
    makehist4(dtree, dhist, dhist_2D,  true, Sample::zee, category, false, 1.00, lumi, dyhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);  
    makehist4(dtree, dhist, dhist_2D,  true, Sample::zee, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("zeecor")+ext;
  if(isEWK)
    name = string("zewkeecor")+ext;

  // Make 1D transfer factor
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2); // smooth numerator in case
      smoothEmptyBins(dhist.at(ihist),2); // smooth denominator in case
    }
  }
  
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  // Make 2D transfer factor
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }    
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){    
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2); // smooth ratio
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }  
 
  outfile.cd();
  outfile.Close();
  kffile.Close();
  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();
  
  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();

  cout << "Z(ee)->Z(inv) transfer factor computed ..." << endl;
}





// make histograms for W->mnu to signal region correction                                                                                                                   
void makewmncorhist( const string &  signalRegionFile,  
		     const string &  wmnFile,   
		     const Category & category, 
                     const SamplesNLO & nloSamples,
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double &  lumi, 
		     const string &  outDir = "", 
		     const string &  sysName = "", 
		     const bool &    isHiggsInvisible = false,
		     const bool &    isEWK = false,
		     const string &  ext = "") {

  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((wmnFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wmn_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wmn_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wmn_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wmn_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wmn_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wmn_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wmn_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wmn_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> whists;

  if(nloSamples.useWJetsNLO)
    whists.push_back(wewkhist);
  else{
    whists.push_back(wnlohist);
    whists.push_back(wewkhist);
  }

  TFile* kfactwjet_vbf = NULL;
  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                               
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets.root");
    TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);
  }

  if(not isEWK){
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D, true, Sample::sig, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D, true, Sample::wmn, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("wmncor")+ext;
  if(isEWK)
    name = string("wewkmncor")+ext;

  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }
  
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }


  outfile.Close();
  kffile.Close();
  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();

  nhist.clear();
  dhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  unrolled.clear();

  cout << "W(mnu)->W+Jets transfer factor computed ..." << endl;
}

// make histograms for W->mnu to signal region correction                                                                                                                   
void makewencorhist( const string &  signalRegionFile,  
		     const string &  wenFile,   
		     const Category & category, 
                     const SamplesNLO & nloSamples,
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double &  lumi, 
		     const string &  outDir = "", 
		     const string &  sysName = "", 
		     const bool &    isHiggsInvisible = false,
		     const bool &    isEWK = false,
		     const string &  ext = "") {

  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((wenFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wen_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wen_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wen_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wen_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wen_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wen_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wen_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wen_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> whists;

  if(nloSamples.useWJetsNLO)
    whists.push_back(wewkhist);
  else{
    whists.push_back(wnlohist);
    whists.push_back(wewkhist);
  }

  TFile* kfactwjet_vbf = NULL;
  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                               
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets.root");
    TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);
  }
  
  if(not isEWK){
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::wen, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D, true, Sample::sig, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D, true, Sample::wen, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("wencor")+ext;
  if(isEWK)
    name = string("wewkencor")+ext;

  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }
  
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }


  outfile.Close();
  kffile.Close();
  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();

  nhist.clear();
  dhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  unrolled.clear();

  cout << "W(mnu)->W+Jets transfer factor computed ..." << endl;
}



// make Z/W ratio
void  makezwjcorhist(const string & znunuFile,  
		     const string & wlnuFile,   
		     const Category & category, 
                     const SamplesNLO & nloSamples,
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &   isHiggsInvisible = false,
		     const bool &   isEWK = false,
		     const string & ext = "",
		     int    kfact = 0) {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((znunuFile+"/*root").c_str());
  dtree->Add((wlnuFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zwj_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zwj_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"zwj_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"zwj_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zwj_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zwj_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"zwj_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"zwj_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");
  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  // in order to make uncertainties use the old file
  TFile kffileUnc (kfactorFileUnc.c_str());
  TH1* zpdfhist = (TH1*) kffileUnc.Get("znlo012/znlo012_pdfUp");
  TH1* wpdfhist = (TH1*) kffileUnc.Get("wnlo012/wnlo012_pdfUp");
  TH1* nomhist  = (TH1*) kffileUnc.Get("znlo1_over_wnlo1/znlo1_over_wnlo1");
  TH1* re1hist  = (TH1*) kffileUnc.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renCorrUp");
  TH1* re2hist  = (TH1*) kffileUnc.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renAcorrUp");
  TH1* fa1hist  = (TH1*) kffileUnc.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facCorrUp");
  TH1* fa2hist  = (TH1*) kffileUnc.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facAcorrUp");

  TH1* znloOrig = (TH1*) kffileUnc.Get("znlo012/znlo012_nominal");
  TH1* wnloOrig = (TH1*) kffileUnc.Get("wnlo012/wnlo012_nominal");

  zpdfhist->Divide(znloOrig);
  wpdfhist->Divide(wnloOrig);

  // Z/W NLO QCD re up / Z/W NLO QCD                                                                                                                                          
  re1hist->Divide(nomhist);
  // Z/W NLO QCD re EWK up / Z/W NLO QCD                                                                                                                                       
  re2hist->Divide(nomhist);
  // Z/W NLO QCD fac  up / Z/W NLO QCD                                                                                                                                         
  fa1hist->Divide(nomhist);
  // Z/W NLO QCD fac EWK up / Z/W NLO QCD                                                                                                                                       
  fa2hist->Divide(nomhist);

  vector<TH1*> zhists;
  vector<TH1*> whists;
  vector<TH1*> ehists;

  //kfact == 1 --> Znunu corrected for by NLO QCD, Wlnu by NLO QCD                                                                                                              
  if (kfact == 1 and not nloSamples.useZJetsNLO) zhists.push_back(znlohist);
  if (kfact == 1 and not nloSamples.useWJetsNLO) whists.push_back(wnlohist);
    
  //kfact == 2 --> Znunu corrected for by NLO QCD+EWK, Wlnu by NLO QCD+EWK                                                                                                      
  if (kfact == 2 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  else if(kfact == 2 and nloSamples.useZJetsNLO) zhists.push_back(zewkhist);

  if (kfact == 2 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(wewkhist);}
  else if (kfact == 2 and nloSamples.useWJetsNLO) {whists.push_back(wewkhist);}

  //kfact == 3 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up QCD                                                                                                        
  if (kfact == 3 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  else if(kfact == 3 and nloSamples.useZJetsNLO) zhists.push_back(re1hist);

  if (kfact == 3 and not nloSamples.useWJetsNLO)  whists.push_back(wnlohist);

  //kfact == 4 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up QCD                                                                                                        
  if (kfact == 4 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  else if(kfact == 4 and nloSamples.useZJetsNLO){ zhists.push_back(fa1hist) ;}
  
  if (kfact == 4 and not nloSamples.useWJetsNLO)  whists.push_back(wnlohist);

  //kfact == 5 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up EWK                                                                                                        
  if (kfact == 5 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  else if(kfact == 5 and nloSamples.useZJetsNLO) {zhists.push_back(re2hist) ;}
  if (kfact == 5 and not nloSamples.useWJetsNLO)  whists.push_back(wnlohist);
  
  //kfact == 6 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up EWK                                                                                                        
  if (kfact == 6 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  else if(kfact == 6 and nloSamples.useZJetsNLO) zhists.push_back(fa2hist);

  if (kfact == 6 and not nloSamples.useWJetsNLO)  whists.push_back(wnlohist);

  //kfact == 7 --> Znunu corrected for by NLO NLO PDF, Wlnu by NLO                                                                                                              
  if (kfact == 7 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  else if(kfact == 7 and nloSamples.useZJetsNLO) zhists.push_back(zpdfhist);

  if (kfact == 7 and not nloSamples.useWJetsNLO)  {whists.push_back(wnlohist); whists.push_back(wpdfhist);}
  else if (kfact == 7 and nloSamples.useWJetsNLO) {whists.push_back(wpdfhist);}

  TFile* kfactzjet_vbf = NULL;
  TFile* kfactwjet_vbf = NULL;

  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                                
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets.root");
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets.root");

    TH1* zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* zjet_nlo_mj  = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");
    zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    zjet_nlo_mj->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);
    if(not nloSamples.useZJetsNLO)
      zhists.push_back(zjet_nlo_vbf);


    TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);
  }



  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon,   
  if(not isEWK){
    if(nloSamples.useZJetsNLO)
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);  
    else
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);    
    makehist4(dtree, dhist, dhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("zwjcor")+ext;
  if(isEWK)
    name = string("zwjewkcor")+ext;

  // divide the two                                                                                                                                                          
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){    
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);  

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }


  outfile.Close();
  kffile.Close();
  kffileUnc.Close();
  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();
  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();
  
  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();
  
  cout << "W+Jets->Z+inv transfer factor computed ..." << endl;
}


// make Z/gamma ratio
void makegamcorhist( const string & znunuFile,  
		     const string & photonFile,  
		     const string & fPfile,  
		     const Category & category, 
		     const SamplesNLO & nloSamples,
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &   isHiggsInvisible = false,
		     const bool &   isEWK = false,
		     const string & ext = "",
		     int    kfact = 0) {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((znunuFile+"/*root").c_str());
  dtree->Add((photonFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"gam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"gam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"gam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"gam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"gam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"gam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"gam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"gam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);
  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  TFile kffileUnc (kfactorFileUnc.c_str());
  TH1* nomhist  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1");
  TH1* zpdfhist = (TH1*) kffileUnc.Get("znlo012/znlo012_pdfUp");
  TH1* apdfhist = (TH1*) kffileUnc.Get("anlo1/anlo1_pdfUp");
  TH1* re1hist  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_renCorrUp");
  TH1* re2hist  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_renAcorrUp");
  TH1* fa1hist  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_facCorrUp");
  TH1* fa2hist  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_facAcorrUp");

  TH1* znloOrig = (TH1*) kffileUnc.Get("znlo012/znlo012_nominal");
  TH1* anloOrig = (TH1*) kffileUnc.Get("anlo1/anlo1_nominal");


  zpdfhist->Divide(znloOrig);
  apdfhist->Divide(anloOrig);

  // Z/gam NLO re QCD Up / Z/gamma NLO                                                                                                                                       
  re1hist->Divide(nomhist);
  // Z/gam NLO re EWK Up / Z/gamma NLO                                                                                                                                        
  re2hist->Divide(nomhist);
  // Z/gam NLO fact QCD Up / Z/gamma NLO                                                                                                                                      
  fa1hist->Divide(nomhist);
  // Z/gam NLO fact EWK Up / Z/gamma NLO                                                                                                                                      
  fa2hist->Divide(nomhist);

  TFile fpfile(fPfile.c_str());
  TH1* afpchist = (TH1*)fpfile.Get("FP_Down");

  vector<TH1*> zhists;
  vector<TH1*> ahists;
  vector<TH1*> ehists;
  // ZNLO QCD and Gamma NLO QCD                                                                                                                                                
  if (kfact == 1 and not nloSamples.useZJetsNLO)      zhists.push_back(znlohist);
  if (kfact == 1 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);

  //ZNLO QCD+EWK and Gamma NLO QCD+EWK                                                                                                                                         
  if (kfact == 2 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  else if (kfact == 2 and nloSamples.useZJetsNLO) {zhists.push_back(zewkhist);}
  if (kfact == 2 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); ahists.push_back(aewkhist);}
  else if(kfact == 2 and nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist);}

  // ZNLO QCD+Re up and Gamma NLO QCD                                                                                                                                          
  if (kfact == 3 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  else if (kfact == 3 and nloSamples.useZJetsNLO) {zhists.push_back(re1hist) ;}

  if(kfact == 3 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist);}

  // ZNLO QCD + fact Up and Gamma NLO QCD                                                                                                                                      
  if (kfact == 4 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  else if (kfact == 4 and nloSamples.useZJetsNLO) {zhists.push_back(fa1hist) ;}

  if(kfact == 4 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist);}

  // ZNLO QCD + re EWK up and Gamma NLO QCD                                                                                                                                   
  if (kfact == 5 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  else if (kfact == 5 and nloSamples.useZJetsNLO) {zhists.push_back(re2hist) ;}

  if(kfact == 5 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist);}

  // ZNLO QCD + fact EWK up and Gamma NLO QCD                                                                                                                                  
  if (kfact == 6 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  else if (kfact == 6 and nloSamples.useZJetsNLO) {zhists.push_back(fa2hist) ;}

  if(kfact == 6 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist);}

  // ZNLO QCD + PDF up and Gamma NLO QCD + PDF Up                                                                                                                               
  if (kfact == 7 and not nloSamples.useZJetsNLO) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  else if (kfact == 7 and nloSamples.useZJetsNLO) {zhists.push_back(zpdfhist);}
  
  if(kfact == 7 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); ahists.push_back(apdfhist);}
  else if(kfact == 7  and nloSamples.usePhotonJetsNLO) ahists.push_back(apdfhist);

  // ZNLO QCD and Gamma NLO QCD + FP                                                                                                                                            
  if (kfact == 8 and not nloSamples.useZJetsNLO) zhists.push_back(znlohist);

  if (kfact == 8 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); zhists.push_back(afpchist);}
  else if(kfact == 8 and nloSamples.usePhotonJetsNLO) {ahists.push_back(afpchist);}


  TFile* kfactzjet_vbf = NULL;
  TFile* kfactgjet_vbf = NULL;

  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                           
  
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets.root");
    kfactgjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_gjets.root");

    TH1* zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* zjet_nlo_mj  = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");
    zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    zjet_nlo_mj->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);
    if(not nloSamples.useZJetsNLO)
      zhists.push_back(zjet_nlo_vbf);
    
    TH1* gjet_nlo_vbf = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* gjet_nlo_mj  = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_monojet");
    gjet_nlo_vbf->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_vbf"));
    gjet_nlo_mj->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_monojet"));
    gjet_nlo_vbf->Divide(gjet_nlo_mj);
    if(not nloSamples.usePhotonJetsNLO)
      ahists.push_back(gjet_nlo_vbf);
  }


  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  if(not isEWK){
    if(not nloSamples.useZJetsNLO)
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, "", false, reweightNVTX, 0, isHiggsInvisible);
    else
      makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhists, "", false, reweightNVTX, 0, isHiggsInvisible);      
    makehist4(dtree, dhist, dhist_2D,  true, Sample::gam, category, false, 1.00, lumi, ahists, "", false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::gam, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  string name = string("gamcor")+ext;
  if(isEWK)
    name = string("gamewkcor")+ext;
  
  // divide the two                                                                                                                                                          
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
  }
  
  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }

  outfile.Close();
  kffile.Close();
  kffileUnc.Close();
  fpfile.Close();
  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();
  if(kfactgjet_vbf)
    kfactgjet_vbf->Close();

  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();

  cout << "Gamma+Jets->Z+inv transfer factor computed ..." << endl;
}


void makewgamcorhist( const string & wlnuFile,  
		      const string & photonFile,  
		      const string & fPfile,  
		      const Category & category, 
		      const SamplesNLO & nloSamples,
		      vector<string> observables, 
		      vector<string> observables_2D, 
		      const double & lumi, 
		      const string & outDir = "", 
		      const string & sysName = "", 
		      const bool &  isHiggsInvisible = false,
		      const bool &  isEWK = false,
		      const string & ext = "",
		      int    kfact = 0) {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((wlnuFile+"/*root").c_str());
  dtree->Add((photonFile+"/*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  string postfix = "_";
  if(isEWK)
    postfix = "_ewk_";

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist"+postfix+"wgam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist"+postfix+"wgam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wgam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wgam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist"+postfix+"wgam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist"+postfix+"wgam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);
  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  // take open loop hitogram                                                                                                                                                    
  TFile kffileUnc (kfactorFileUnc.c_str());
  TH1* nomhist    = (TH1*) kffileUnc.Get("wnlo1_over_znlo1/wnlo1_over_znlo1");
  TH1* nomhist_2  = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1");
  nomhist->Multiply(nomhist_2);

  TH1* wpdfhist = (TH1*) kffileUnc.Get("wnlo012/wnlo012_pdfUp");
  TH1* apdfhist = (TH1*) kffileUnc.Get("anlo1/anlo1_pdfUp");
  
  TH1* re1hist   = (TH1*) kffileUnc.Get("wnlo1_over_znlo1/wnlo1_over_znlo1_renCorrUp");
  TH1* re1hist_2 = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_renCorrUp");
  re1hist->Multiply(re1hist_2);

  TH1* re2hist   = (TH1*) kffileUnc.Get("wnlo1_over_znlo1/wnlo1_over_znlo1_renAcorrUp");
  TH1* re2hist_2 = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_renAcorrUp");
  re2hist->Multiply(re2hist_2);

  TH1* fa1hist   = (TH1*) kffileUnc.Get("wnlo1_over_znlo1/wnlo1_over_znlo1_facCorrUp");
  TH1* fa1hist_2 = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_facCorrUp");
  fa1hist->Multiply(fa1hist_2);


  TH1* fa2hist   = (TH1*) kffileUnc.Get("wnlo1_over_znlo1/wnlo1_over_znlo1_facAcorrUp");
  TH1* fa2hist_2 = (TH1*) kffileUnc.Get("znlo1_over_anlo1/znlo1_over_anlo1_facAcorrUp");
  fa2hist->Multiply(fa2hist_2);

  TH1* wnloOrig = (TH1*) kffileUnc.Get("wnlo012/wnlo012_nominal");
  TH1* anloOrig = (TH1*) kffileUnc.Get("anlo1/anlo1_nominal");
  
  wpdfhist->Divide(wnloOrig);
  apdfhist->Divide(anloOrig);

  // Z/gam NLO re QCD Up / Z/gamma NLO                                                                                                                                       
  re1hist->Divide(nomhist);  
  // Z/gam NLO re EWK Up / Z/gamma NLO                                                                                                                                        
  re2hist->Divide(nomhist);
  // Z/gam NLO fact QCD Up / Z/gamma NLO                                                                                                                                      
  fa1hist->Divide(nomhist);
  // Z/gam NLO fact EWK Up / Z/gamma NLO                                                                                                                                      
  fa2hist->Divide(nomhist);

  TFile fpfile(fPfile.c_str());
  TH1* afpchist = (TH1*)fpfile.Get("FP_Down");

  vector<TH1*> whists;
  vector<TH1*> ahists;
  vector<TH1*> ehists;

  // ZNLO QCD and Gamma NLO QCD                                                                                                                                                
  if (kfact == 1 and not nloSamples.useWJetsNLO) whists.push_back(wnlohist);
  if (kfact == 1 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);

  //ZNLO QCD+EWK and Gamma NLO QCD+EWK                                                                                                                                         
  if (kfact == 2 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(wewkhist);}
  else if (kfact == 2 and nloSamples.useWJetsNLO) {whists.push_back(wewkhist);}

  if (kfact == 2 and nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); ahists.push_back(aewkhist);}
  else if (kfact == 2 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(aewkhist);}
  
  // ZNLO QCD+Re up and Gamma NLO QCD                                                                                                                                          
  if (kfact == 3 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(re1hist) ;}
  else if (kfact == 3 and nloSamples.useWJetsNLO) {whists.push_back(re1hist) ;}

  if (kfact == 3 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);
  
  // ZNLO QCD + fact Up and Gamma NLO QCD                                                                                                                                      
  if (kfact == 4 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(fa1hist) ;}
  else if (kfact == 4 and nloSamples.useWJetsNLO) {whists.push_back(fa1hist) ;}
  if (kfact == 4 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);

  // ZNLO QCD + re EWK up and Gamma NLO QCD                                                                                                                                   
  if (kfact == 5 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(re2hist) ;}
  else if (kfact == 5 and nloSamples.useWJetsNLO) {whists.push_back(re2hist) ;}
  if (kfact == 5 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);

  // ZNLO QCD + fact EWK up and Gamma NLO QCD                                                                                                                                  
  if (kfact == 6 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(fa2hist) ;}
  else if (kfact == 6 and nloSamples.useWJetsNLO) {whists.push_back(fa2hist) ;}

  if (kfact == 6 and not nloSamples.usePhotonJetsNLO) ahists.push_back(anlohist);

  // ZNLO QCD + PDF up and Gamma NLO QCD + PDF Up                                                                                                                               
  if (kfact == 7 and not nloSamples.useWJetsNLO) {whists.push_back(wnlohist); whists.push_back(wpdfhist);}
  else if (kfact == 7 and nloSamples.useWJetsNLO) {whists.push_back(wpdfhist);}

  if (kfact == 7 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); ahists.push_back(apdfhist);}
  else if(kfact == 7 and nloSamples.usePhotonJetsNLO) {ahists.push_back(apdfhist);}

  // ZNLO QCD and Gamma NLO QCD + FP                                                                                                                                            
  if (kfact == 8 and not nloSamples.useWJetsNLO) whists.push_back(wnlohist);
  if (kfact == 8 and not nloSamples.usePhotonJetsNLO) {ahists.push_back(anlohist); ahists.push_back(afpchist);}
  else if(kfact == 8 and nloSamples.usePhotonJetsNLO) {ahists.push_back(afpchist);}
  

  TFile* kfactwjet_vbf = NULL;
  TFile* kfactgjet_vbf = NULL;

  if(category == Category::VBF){ // apply further k-factors going to the VBF selections                                                                                                                

    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets.root");
    kfactgjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_gjets.root");

    TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);

    TH1* gjet_nlo_vbf = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_vbf");
    TH1* gjet_nlo_mj  = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_monojet");
    gjet_nlo_vbf->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_vbf"));
    gjet_nlo_mj->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_monojet"));
    gjet_nlo_vbf->Divide(gjet_nlo_mj);
    if(not nloSamples.usePhotonJetsNLO)
      ahists.push_back(gjet_nlo_vbf);
  }


  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  if(isEWK){
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, "", false, reweightNVTX, 0, isHiggsInvisible);  
    makehist4(dtree, dhist, dhist_2D,  true, Sample::gam, category, false, 1.00, lumi, ahists, "", false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::gam, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }  

  string name = string("wgamcor")+ext;
  if(isEWK)
    name = string("wgamewkcor")+ext;

  // divide the two                                                                                                                                                          
  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist.size(); ihist++){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  if(doSmoothing){
    for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
  }

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }

  outfile.Close();
  kffile.Close();
  kffileUnc.Close();
  fpfile.Close();
  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();

  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();
  if(kfactgjet_vbf)
    kfactgjet_vbf->Close();

  cout << "Gamma+Jets->W+lnu transfer factor computed ..." << endl;
}


// correction for top
void maketopmucorhist( const string & signalRegionFile,  
		       const string & topFile,
		       const Category & category, 
		       vector<string> observables, 
		       vector<string> observables_2D, 
		       const double & lumi,
		       const string & signalRegionFile_alt = "", 
		       const string & topFile_alt = "", 
		       const string & outDir  = "", 
		       const string & sysName = "", 
		       const bool &   isHiggsInvisible = false,
		       const string & ext = ""){

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((topFile+"/*root").c_str());

  TChain* ntree_alt = NULL;
  TChain* dtree_alt = NULL;

  if(signalRegionFile_alt != ""){
    ntree_alt = new TChain("tree/tree");
    dtree_alt->Add((signalRegionFile_alt+"/*root").c_str());
  }
  if(topFile_alt != ""){
    dtree_alt = new TChain("tree/tree");
    dtree_alt->Add((topFile_alt+"/*root").c_str());
  }
  
  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH1*> nhist_alt;
  vector<TH1*> dhist_alt;
  vector<TH1*> tfhist_alt;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH2*> nhist_2D_alt;
  vector<TH2*> dhist_2D_alt;
  vector<TH2*> tfhist_2D_alt;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_topmu_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_topmu_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* nhist_alt_temp = new TH1F(("nhist_alt_topmu_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_alt_temp = new TH1F(("dhist_alt_topmu_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
      nhist_alt.push_back(dynamic_cast<TH1*>(nhist_alt_temp));
      dhist_alt.push_back(dynamic_cast<TH1*>(dhist_alt_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_topmu_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_topmu_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* nhist_alt_temp = new TH1F(("nhist_alt_topmu_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_alt_temp = new TH1F(("dhist_alt_topmu_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
      nhist_alt.push_back(dynamic_cast<TH1*>(nhist_alt_temp));
      dhist_alt.push_back(dynamic_cast<TH1*>(dhist_alt_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_topmu_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_topmu_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* nhist_alt_temp = new TH2F(("nhist_alt_topmu_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_alt_temp = new TH2F(("dhist_alt_topmu_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
      nhist_2D_alt.push_back(dynamic_cast<TH2*>(nhist_alt_temp));
      dhist_2D_alt.push_back(dynamic_cast<TH2*>(dhist_alt_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_topmu_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_topmu_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* nhist_alt_temp = new TH2F(("nhist_alt_topmu_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_alt_temp = new TH2F(("dhist_alt_topmu_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
      nhist_2D_alt.push_back(dynamic_cast<TH2*>(nhist_alt_temp));
      dhist_2D_alt.push_back(dynamic_cast<TH2*>(dhist_alt_temp));      
    }
  }
  vector<TH1*> ehists;
  vector<TH1*> zhists;

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 7 == b-tagged region, 
  makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(dtree, dhist, dhist_2D,  true, Sample::topmu, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(ntree_alt, nhist_alt, nhist_2D_alt,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(dtree_alt, dhist_alt, dhist_2D_alt,  true, Sample::topmu, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);

  string name = string("topmucor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    if(doSmoothing){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
    if(nhist_alt.size() >= ihist){
      if(doSmoothing){
	smoothEmptyBins(nhist_alt.at(ihist),2);
	smoothEmptyBins(dhist_alt.at(ihist),2);
      }
      tfhist_alt.push_back((TH1*) nhist_alt.at(ihist)->Clone(Form("%s_temp",nhist_alt.at(ihist)->GetName())));
      tfhist_alt.back()->Divide(dhist_alt.at(ihist));
    }
  }


  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    if(doSmoothing){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
    if(nhist_2D_alt.size() >= ihist){
      if(doSmoothing){
	smoothEmptyBins(nhist_2D_alt.at(ihist),2);
	smoothEmptyBins(dhist_2D_alt.at(ihist),2);
      }
      tfhist_2D_alt.push_back((TH2*) nhist_2D_alt.at(ihist)->Clone(Form("%s_temp",nhist_2D_alt.at(ihist)->GetName())));
      tfhist_2D_alt.back()->Divide(dhist_2D_alt.at(ihist));
    }
  }


  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    for(size_t ihist = 0; ihist < tfhist_alt.size(); ihist++)
      smoothEmptyBins(tfhist_alt.at(ihist),2);
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
    for(size_t ihist = 0; ihist < tfhist_2D_alt.size(); ihist++)
      smoothEmptyBins(tfhist_2D_alt.at(ihist),1);
  }

  // make average
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    makeAverage(tfhist.at(ihist),tfhist_alt.at(ihist));

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    makeAverage(tfhist_2D.at(ihist),tfhist_2D_alt.at(ihist));

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }

  outfile.Close();
  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  nhist_alt.clear();
  dhist_alt.clear();
  tfhist_alt.clear();
  nhist_2D_alt.clear();
  dhist_2D_alt.clear();
  tfhist_2D_alt.clear();
  unrolled.clear();

  cout << "Top(b-tag,mu)->Top(b-veto) transfer factor computed ..." << endl;
}


// correction for top
void maketopelcorhist( const string & signalRegionFile,  
		       const string & topFile,  
		       const Category &    category, 
		       vector<string> observables, 
		       vector<string> observables_2D, 
		       const double & lumi, 
		       const string & signalRegionFile_alt = "", 
		       const string & topFile_alt = "", 
		       const string & outDir      = "", 
		       const string & sysName = "", 
		       const bool   & isHiggsInvisible = false,
		       const string & ext = "") {
  
  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((topFile+"/*root").c_str());

  TChain* ntree_alt = NULL;
  TChain* dtree_alt = NULL;

  if(signalRegionFile_alt != ""){
    ntree_alt = new TChain("tree/tree");
    dtree_alt->Add((signalRegionFile_alt+"/*root").c_str());
  }
  if(topFile_alt != ""){
    dtree_alt = new TChain("tree/tree");
    dtree_alt->Add((topFile_alt+"/*root").c_str());
  }

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH1*> nhist_alt;
  vector<TH1*> dhist_alt;
  vector<TH1*> tfhist_alt;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH2*> nhist_2D_alt;
  vector<TH2*> dhist_2D_alt;
  vector<TH2*> tfhist_2D_alt;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_topel_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_topel_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* nhist_alt_temp = new TH1F(("nhist_alt_topel_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_alt_temp = new TH1F(("dhist_alt_topel_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
      nhist_alt.push_back(dynamic_cast<TH1*>(nhist_alt_temp));
      dhist_alt.push_back(dynamic_cast<TH1*>(dhist_alt_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_topel_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_topel_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* nhist_alt_temp = new TH1F(("nhist_alt_topel_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_alt_temp = new TH1F(("dhist_alt_topel_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
      nhist_alt.push_back(dynamic_cast<TH1*>(nhist_alt_temp));
      dhist_alt.push_back(dynamic_cast<TH1*>(dhist_alt_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_topel_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_topel_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* nhist_alt_temp = new TH2F(("nhist_alt_topel_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_alt_temp = new TH2F(("dhist_alt_topel_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
      nhist_2D_alt.push_back(dynamic_cast<TH2*>(nhist_alt_temp));
      dhist_2D_alt.push_back(dynamic_cast<TH2*>(dhist_alt_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_topel_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_topel_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* nhist_alt_temp = new TH2F(("nhist_alt_topel_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_alt_temp = new TH2F(("dhist_alt_topel_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
      nhist_2D_alt.push_back(dynamic_cast<TH2*>(nhist_alt_temp));
      dhist_2D_alt.push_back(dynamic_cast<TH2*>(dhist_alt_temp));
    }
  }

  vector<TH1*> ehists;
  vector<TH1*> zhists;

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 7 == b-tagged region,
  makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(dtree, dhist, dhist_2D,  true, Sample::topel, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(ntree_alt, nhist_alt, nhist_2D_alt,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(dtree_alt, dhist_alt, dhist_2D_alt,  true, Sample::topel, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);

  string name = string("topelcor")+ext;


  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    if(doSmoothing){
      smoothEmptyBins(nhist.at(ihist),2);
      smoothEmptyBins(dhist.at(ihist),2);
    }
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
    if(nhist_alt.size() >= ihist){
      if(doSmoothing){
	smoothEmptyBins(nhist_alt.at(ihist),2);
	smoothEmptyBins(dhist_alt.at(ihist),2);
      }
      tfhist_alt.push_back((TH1*) nhist_alt.at(ihist)->Clone(Form("%s_temp",nhist_alt.at(ihist)->GetName())));
      tfhist_alt.back()->Divide(dhist_alt.at(ihist));
    }
  }


  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    if(doSmoothing){
      smoothEmptyBins(nhist_2D.at(ihist),2);
      smoothEmptyBins(dhist_2D.at(ihist),2);
    }
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
    if(nhist_2D_alt.size() >= ihist){
      if(doSmoothing){
	smoothEmptyBins(nhist_2D_alt.at(ihist),2);
	smoothEmptyBins(dhist_2D_alt.at(ihist),2);
      }
      tfhist_2D_alt.push_back((TH2*) nhist_2D_alt.at(ihist)->Clone(Form("%s_temp",nhist_2D_alt.at(ihist)->GetName())));
      tfhist_2D_alt.back()->Divide(dhist_2D_alt.at(ihist));
    }
  }


  //check for empty bins and apply smoothing
  if(doSmoothing){
    for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
      smoothEmptyBins(tfhist.at(ihist),2);
    for(size_t ihist = 0; ihist < tfhist_alt.size(); ihist++)
      smoothEmptyBins(tfhist_alt.at(ihist),2);
    for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
      smoothEmptyBins(tfhist_2D.at(ihist),1);
    for(size_t ihist = 0; ihist < tfhist_2D_alt.size(); ihist++)
      smoothEmptyBins(tfhist_2D_alt.at(ihist),1);
  }

  // make average
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    makeAverage(tfhist.at(ihist),tfhist_alt.at(ihist));

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    makeAverage(tfhist_2D.at(ihist),tfhist_2D_alt.at(ihist));

  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");

  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist.size(); ihist++)
    dhist.at(ihist)->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist.size(); ihist++){
    tfhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    tfhist.at(ihist)->Write("",TObject::kOverwrite);
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    unroll2DHistograms(nhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < dhist_2D.size(); ihist++)
    unroll2DHistograms(dhist_2D.at(ihist))->Write("",TObject::kOverwrite);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++){
    tfhist_2D.at(ihist)->SetName((name+"hist_"+observables_2D.at(ihist)+"_2D").c_str());
    unrolled.push_back(unroll2DHistograms(tfhist_2D.at(ihist)));
    unrolled.back()->Write("",TObject::kOverwrite);
  }

  outfile.Close();
  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  nhist_alt.clear();
  dhist_alt.clear();
  tfhist_alt.clear();
  nhist_2D_alt.clear();
  dhist_2D_alt.clear();
  tfhist_2D_alt.clear();
  unrolled.clear();

  cout << "Top(b-tag,el)->Top(b-veto) transfer factor computed ..." << endl;
}

//  LocalWords:  signalRegionFile
