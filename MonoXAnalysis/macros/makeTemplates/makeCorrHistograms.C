#include "makehist.h"
#include "TChain.h"
#include "TF1.h"

using namespace std;

// make histograms for Z->mumu to signal region correction                                                                                                                   
void makezmmcorhist( const string &   signalRegionFile, 
		     const string &   signalRegionEWKFile,
		     const string &   zmumuFile,  
		     const string &   zmumuEWKFile,  
		     const Category & category, 
		     vector<string>   observables, 
		     vector<string>   observables_2D, 
		     const double &   lumi, 
		     const string &   outDir = "", 
		     const string &   sysName = "", 
		     const bool &     isHiggsInvisible = false,
		     const bool &     useNLOSamples = false,
		     const string &   ext = "") {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((zmumuFile+"/*root").c_str());
  if(signalRegionEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((signalRegionEWKFile+"/*root").c_str());
  }
  if(zmumuEWKFile != ""){
    dtreeEWK = new TChain("tree/tree");
    dtreeEWK->Add((zmumuEWKFile+"/*root").c_str());
  }

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_zmm_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zmm_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_zmm_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zmm_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }
  
  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_zmm_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zmm_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_zmm_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zmm_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  vector<TH1*> zhistsnlo;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);
  zhistsnlo.push_back(zewkhist);

  if(useNLOSamples){
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible,false,-1, NULL, NULL, true);
  }
  makehist4(dtree, dhist, dhist_2D,  true, Sample::zmm, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);

  if(ntreeEWK != NULL and ntreeEWK != 0){
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  if(dtreeEWK != NULL and dtreeEWK != 0){
    makehist4(dtreeEWK, dhist, dhist_2D,  true, Sample::zmm, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  string name = string("zmmcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2); // smooth numerator in case
    smoothEmptyBins(dhist.at(ihist),2); // smooth denominator in case
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2); // smooth ratio
  
  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);

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
void makezeecorhist( const string &  signalRegionFile,  
		     const string &  signalRegionEWKFile,  
		     const string &  zeeFile,   
		     const string &  zeeEWKFile,   
		     const Category & category, 
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &   isHiggsInvisible = false,
		     const bool &   useNLOSamples = false,
		     string  ext = "") {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((zeeFile+"/*root").c_str());
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  if(signalRegionEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((signalRegionEWKFile+"/*root").c_str());
  }
  if(zeeEWKFile != ""){
    dtreeEWK = new TChain("tree/tree");
    dtreeEWK->Add((zeeEWKFile+"/*root").c_str());
  }


  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_zee_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zee_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_zee_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zee_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_zee_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zee_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_zee_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zee_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  vector<TH1*> zhistsnlo;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);
  zhistsnlo.push_back(zewkhist);
  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  if(useNLOSamples){
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhistsnlo, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName,false, reweightNVTX, 0, isHiggsInvisible,false,-1, NULL, NULL, true);
  }
  makehist4(dtree, dhist, dhist_2D,  true, Sample::zee, category, false, 1.00, lumi, zhists, sysName,false, reweightNVTX, 0, isHiggsInvisible);

  if(ntreeEWK != NULL and ntreeEWK != 0)
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  if(dtreeEWK != NULL and dtreeEWK != 0)  
    makehist4(dtreeEWK, dhist, dhist_2D,  true, Sample::zee, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);

  string name = string("zeecor")+ext;

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);

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
		     const string &  signalRegionEWKFile,  
		     const string &  wmnFile,   
		     const string &  wmnEWKFile,   
		     const Category & category, 
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double &  lumi, 
		     const string &  outDir = "", 
		     const string &  sysName = "", 
		     const bool &    isHiggsInvisible = false,
		     const bool &    useNLOSamples = false,
		     const string &  ext = "") {

  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((wmnFile+"/*root").c_str());
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  if(signalRegionEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((signalRegionEWKFile+"/*root").c_str());
  }
  if(wmnEWKFile != ""){
    dtreeEWK = new TChain("tree/tree");
    dtreeEWK->Add((wmnEWKFile+"/*root").c_str());
  }

  // open files --> temp fix for W+jets xs                                                                                                                                      
  TChain* ntree_nlo1 = new TChain("tree/tree");
  TChain* ntree_nlo2 = new TChain("tree/tree");
  TChain* ntree_nlo3 = new TChain("tree/tree");
  TChain* ntree_nlo4 = new TChain("tree/tree");
  TChain* dtree_nlo1 = new TChain("tree/tree");
  TChain* dtree_nlo2 = new TChain("tree/tree");
  TChain* dtree_nlo3 = new TChain("tree/tree");
  TChain* dtree_nlo4 = new TChain("tree/tree");
                          
  ntree_nlo1->Add((signalRegionFile+"/*Pt-100To250*root").c_str());
  ntree_nlo2->Add((signalRegionFile+"/*Pt-250To400*root").c_str());
  ntree_nlo3->Add((signalRegionFile+"/*Pt-400To600*root").c_str());
  ntree_nlo4->Add((signalRegionFile+"/*Pt-600ToInf*root").c_str());

  dtree_nlo1->Add((wmnFile+"/*Pt-100To250*root").c_str());
  dtree_nlo2->Add((wmnFile+"/*Pt-250To400*root").c_str());
  dtree_nlo3->Add((wmnFile+"/*Pt-400To600*root").c_str());
  dtree_nlo4->Add((wmnFile+"/*Pt-600ToInf*root").c_str());


  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_wmn_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_wmn_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_wmn_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_wmn_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_wmn_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_wmn_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_wmn_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_wmn_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  vector<TH1*> whistsnlo;

  whists.push_back(wnlohist);
  whists.push_back(wewkhist);
  whistsnlo.push_back(wewkhist);
  
  if(useNLOSamples){
    makehist4(ntree_nlo1, nhist, nhist_2D,  true, Sample::sig, category, false, 1.10, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo2, nhist, nhist_2D,  true, Sample::sig, category, false, 1.11, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo3, nhist, nhist_2D,  true, Sample::sig, category, false, 1.15, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo4, nhist, nhist_2D,  true, Sample::sig, category, false, 1.13, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo1, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.10, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo2, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.11, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo3, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.15, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo4, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.13, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  if(ntreeEWK != NULL and ntreeEWK != 0)
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  if(dtreeEWK != NULL and dtreeEWK != 0)
    makehist4(dtreeEWK, dhist, dhist_2D,  true, Sample::wmn, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);

  string name = string("wmncor")+ext;

  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);

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

  nhist.clear();
  dhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  unrolled.clear();

  cout << "W(mnu)->W+Jets transfer factor computed ..." << endl;
}


// make histograms for W->enu to signal region correction                                                                                                                   
void makewencorhist( const string & signalRegionFile,  
		     const string & signalRegionEWKFile,  
		     const string & wenFile,   
		     const string & wenEWKFile,   
		     const Category & category, 
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &   isHiggsInvisible = false,
		     const bool &   useNLOSamples = false,
		     string ext = "") {

  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((signalRegionFile+"/*root").c_str());
  dtree->Add((wenFile+"/*root").c_str());
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  if(signalRegionEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((signalRegionEWKFile+"/*root").c_str());
  }
  if(wenEWKFile != ""){
    dtreeEWK = new TChain("tree/tree");
    dtreeEWK->Add((wenEWKFile+"/*root").c_str());
  }

  // open files --> temp fix for NLO xs                                                                                                                                 
  TChain* ntree_nlo1 = new TChain("tree/tree");
  TChain* ntree_nlo2 = new TChain("tree/tree");
  TChain* ntree_nlo3 = new TChain("tree/tree");
  TChain* ntree_nlo4 = new TChain("tree/tree");
  TChain* dtree_nlo1 = new TChain("tree/tree");
  TChain* dtree_nlo2 = new TChain("tree/tree");
  TChain* dtree_nlo3 = new TChain("tree/tree");
  TChain* dtree_nlo4 = new TChain("tree/tree");

  ntree_nlo1->Add((signalRegionFile+"/*Pt-100To250*root").c_str());
  ntree_nlo2->Add((signalRegionFile+"/*Pt-250To400*root").c_str());
  ntree_nlo3->Add((signalRegionFile+"/*Pt-400To600*root").c_str());
  ntree_nlo4->Add((signalRegionFile+"/*Pt-600ToInf*root").c_str());

  dtree_nlo1->Add((wenFile+"/*Pt-100To250*root").c_str());
  dtree_nlo2->Add((wenFile+"/*Pt-250To400*root").c_str());
  dtree_nlo3->Add((wenFile+"/*Pt-400To600*root").c_str());
  dtree_nlo4->Add((wenFile+"/*Pt-600ToInf*root").c_str());
  

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_wen_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_wen_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_wen_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_wen_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
    
    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_wen_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_wen_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_wen_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_wen_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  vector<TH1*> whistsnlo;

  whists.push_back(wnlohist);
  whists.push_back(wewkhist);
  whistsnlo.push_back(wewkhist);

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  if(useNLOSamples){
    makehist4(ntree_nlo1, nhist, nhist_2D,  true, Sample::sig, category, false, 1.10, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo2, nhist, nhist_2D,  true, Sample::sig, category, false, 1.11, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo3, nhist, nhist_2D,  true, Sample::sig, category, false, 1.15, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(ntree_nlo4, nhist, nhist_2D,  true, Sample::sig, category, false, 1.13, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo1, dhist, dhist_2D,  true, Sample::wen, category, false, 1.10, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo2, dhist, dhist_2D,  true, Sample::wen, category, false, 1.11, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo3, dhist, dhist_2D,  true, Sample::wen, category, false, 1.15, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo4, dhist, dhist_2D,  true, Sample::wen, category, false, 1.13, lumi, whistsnlo, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }
  else{
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree, dhist, dhist_2D,  true, Sample::wen, category, false, 1, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  if(ntreeEWK != NULL and ntreeEWK != 0)
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);
  if(dtreeEWK != NULL and dtreeEWK != 0)
    makehist4(dtreeEWK, dhist, dhist_2D,  true, Sample::wen, category, false, 1.00, lumi, ehists, sysName,false, reweightNVTX, 0, isHiggsInvisible);


  string name = string("wencor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);

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

  nhist.clear();
  dhist.clear();
  tfhist.clear();
  nhist_2D.clear();
  dhist_2D.clear();
  tfhist_2D.clear();
  unrolled.clear();

  cout << "W(enu)->W+Jets transfer factor computed ..." << endl;
}


// make Z/W ratio
void  makezwjcorhist(const string & znunuFile,  
		     const string & znunuEWKFile,  
		     const string & wlnuFile,   
		     const string & wlnuEWKFile,   
		     const Category & category, 
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &   isHiggsInvisible = false,
		     const bool &   useNLOSamples = false,
		     const string & ext = "",
		     int    kfact = 0) {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((znunuFile+"/*root").c_str());
  dtree->Add((wlnuFile+"/*root").c_str());
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  if(znunuEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((znunuEWKFile+"/*root").c_str());
  }
  if(wlnuEWKFile != ""){
    dtreeEWK = new TChain("tree/tree");
    dtreeEWK->Add((wlnuEWKFile+"/*root").c_str());
  }

  // open files --> temp fix for NLO xs                                                                                                                                         
  TChain* dtree_nlo1 = new TChain("tree/tree");
  TChain* dtree_nlo2 = new TChain("tree/tree");
  TChain* dtree_nlo3 = new TChain("tree/tree");
  TChain* dtree_nlo4 = new TChain("tree/tree");

  dtree_nlo1->Add((wlnuFile+"/*Pt-100To250*root").c_str());
  dtree_nlo2->Add((wlnuFile+"/*Pt-250To400*root").c_str());
  dtree_nlo3->Add((wlnuFile+"/*Pt-400To600*root").c_str());
  dtree_nlo4->Add((wlnuFile+"/*Pt-600ToInf*root").c_str());

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_zwj_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zwj_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_zwj_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_zwj_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }

  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_zwj_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zwj_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_zwj_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_zwj_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  if (kfact == 1  and not useNLOSamples) zhists.push_back(znlohist);
  if (kfact == 1  and not useNLOSamples) whists.push_back(wnlohist);
    
  //kfact == 2 --> Znunu corrected for by NLO QCD+EWK, Wlnu by NLO QCD+EWK                                                                                                      
  if (kfact == 2 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  else if (kfact == 2 and useNLOSamples) {zhists.push_back(zewkhist);}
  if (kfact == 2 and not useNLOSamples) {whists.push_back(wnlohist); whists.push_back(wewkhist);}
  else if (kfact == 2 and useNLOSamples) {whists.push_back(wewkhist);}
  //kfact == 3 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up QCD                                                                                                        
  if (kfact == 3 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  else if (kfact == 3 and useNLOSamples) {zhists.push_back(re1hist) ;}
  if (kfact == 3 and not useNLOSamples)  whists.push_back(wnlohist);
  //kfact == 4 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up QCD                                                                                                        
  if (kfact == 4 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  else if (kfact == 4 and useNLOSamples) {zhists.push_back(fa1hist) ;}
  if (kfact == 4 and not useNLOSamples)  whists.push_back(wnlohist);
  //kfact == 5 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up EWK                                                                                                        
  if (kfact == 5 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  else if (kfact == 5 and useNLOSamples) {zhists.push_back(re2hist) ;}
  if (kfact == 5 and not useNLOSamples)  whists.push_back(wnlohist);
  //kfact == 6 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up EWK                                                                                                        
  if (kfact == 6 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  else if (kfact == 6 and useNLOSamples) {zhists.push_back(fa2hist) ;}
  if (kfact == 6 and not useNLOSamples)  whists.push_back(wnlohist);
  //kfact == 7 --> Znunu corrected for by NLO NLO PDF, Wlnu by NLO                                                                                                              
  if (kfact == 7 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  else if (kfact == 7 and useNLOSamples) {zhists.push_back(zpdfhist);}
  if (kfact == 7 and not useNLOSamples) {whists.push_back(wnlohist); whists.push_back(wpdfhist);}
  else if (kfact == 7 and useNLOSamples) {whists.push_back(wpdfhist);}
  
  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  if(not useNLOSamples)
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible, false,-1, NULL, NULL, true);
  else
    makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 3.00, lumi, zhists, sysName, false, reweightNVTX, 0, isHiggsInvisible);

  if(not useNLOSamples)
    makehist4(dtree, dhist, dhist_2D,  true, Sample::sig, category, false, 1.00, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  else{
    makehist4(dtree_nlo1, dhist, dhist_2D,  true, Sample::sig, category, false, 1.10, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo2, dhist, dhist_2D,  true, Sample::sig, category, false, 1.11, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo3, dhist, dhist_2D,  true, Sample::sig, category, false, 1.15, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
    makehist4(dtree_nlo4, dhist, dhist_2D,  true, Sample::sig, category, false, 1.13, lumi, whists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  }

  if(ntreeEWK != NULL and ntreeEWK != 0)
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  if(dtreeEWK != NULL and dtreeEWK != 0)
    makehist4(dtreeEWK, dhist, dhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  
  string name = string("zwjcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);
  
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
		     const string & znunuEWKFile,  
		     const string & photonFile,  
		     const string & fPfile,  
		     const Category & category, 
		     vector<string> observables, 
		     vector<string> observables_2D, 
		     const double & lumi, 
		     const string & outDir = "", 
		     const string & sysName = "", 
		     const bool &  isHiggsInvisible = false,
		     const bool &  useNLOSamples = false,
		     const string & ext = "",
		     int    kfact = 0) {

  // open files                                                                                                                                                                
  TChain* ntree = new TChain("tree/tree");
  TChain* dtree = new TChain("tree/tree");
  ntree->Add((znunuFile+"/*root").c_str());
  dtree->Add((photonFile+"/*root").c_str());
  TChain* ntreeEWK = NULL;
  TChain* dtreeEWK = NULL;
  if(znunuEWKFile != ""){
    ntreeEWK = new TChain("tree/tree"); 
    ntreeEWK->Add((znunuEWKFile+"/*root").c_str());
  }

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH1*> tfhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;
  vector<TH2*> tfhist_2D;
  vector<TH1*> unrolled;

  vector<double> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH1F* nhist_temp = new TH1F(("nhist_gam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_gam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
    else{
      TH1F* nhist_temp = new TH1F(("nhist_gam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_gam_"+ext+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
    }
  }

  for(auto obs : observables_2D){
    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(ext == ""){
      TH2F* nhist_temp = new TH2F(("nhist_gam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_gam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      nhist_2D.push_back(dynamic_cast<TH2*>(nhist_temp));
      dhist_2D.push_back(dynamic_cast<TH2*>(dhist_temp));
    }
    else{
      TH2F* nhist_temp = new TH2F(("nhist_gam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dhist_temp = new TH2F(("dhist_gam_"+ext+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
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
  if (kfact == 1 and not useNLOSamples) zhists.push_back(znlohist);
  if (kfact == 1) ahists.push_back(anlohist);

  //ZNLO QCD+EWK and Gamma NLO QCD+EWK                                                                                                                                         
  if (kfact == 2 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  else if (kfact == 2 and useNLOSamples) {zhists.push_back(zewkhist);}
  if (kfact == 2) {ahists.push_back(anlohist); ahists.push_back(aewkhist);}

  // ZNLO QCD+Re up and Gamma NLO QCD                                                                                                                                          
  if (kfact == 3 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  else if (kfact == 3 and useNLOSamples) {zhists.push_back(re1hist) ;}
  if (kfact == 3) ahists.push_back(anlohist);

  // ZNLO QCD + fact Up and Gamma NLO QCD                                                                                                                                      
  if (kfact == 4 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  else if (kfact == 4 and useNLOSamples) {zhists.push_back(fa1hist) ;}
  if (kfact == 4) ahists.push_back(anlohist);

  // ZNLO QCD + re EWK up and Gamma NLO QCD                                                                                                                                   
  if (kfact == 5 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  else if (kfact == 5 and useNLOSamples) {zhists.push_back(re2hist) ;}
  if (kfact == 5) ahists.push_back(anlohist);

  // ZNLO QCD + fact EWK up and Gamma NLO QCD                                                                                                                                  
  if (kfact == 6 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  else if (kfact == 6 and useNLOSamples) {zhists.push_back(fa2hist) ;}
  if (kfact == 6) ahists.push_back(anlohist);

  // ZNLO QCD + PDF up and Gamma NLO QCD + PDF Up                                                                                                                               
  if (kfact == 7 and not useNLOSamples) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  else if (kfact == 7 and useNLOSamples) {zhists.push_back(zpdfhist);}
  if (kfact == 7) {ahists.push_back(anlohist); ahists.push_back(apdfhist);}

  // ZNLO QCD and Gamma NLO QCD + FP                                                                                                                                            
  if (kfact == 8 and not useNLOSamples) zhists.push_back(znlohist);
  if (kfact == 8) {ahists.push_back(anlohist); zhists.push_back(afpchist);}


  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, zhists, "", false, reweightNVTX, 0, isHiggsInvisible, false,-1, NULL, NULL, true);
  if(ntreeEWK != NULL and ntreeEWK != 0)
    makehist4(ntreeEWK, nhist, nhist_2D,  true, Sample::sig, category, false, 1.00, lumi, ehists, sysName, false, reweightNVTX, 0, isHiggsInvisible);
  makehist4(dtree, dhist, dhist_2D,  true, Sample::gam, category, false, 1.00, lumi, ahists, "", false, reweightNVTX, 0, isHiggsInvisible);

  string name = string("gamcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
  }

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);

  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);
  
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

  cout << "Gamma+Jets->Z+inv transfer factor computed ..." << endl;
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
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
    if(nhist_alt.size() >= ihist){
      smoothEmptyBins(nhist_alt.at(ihist),2);
      smoothEmptyBins(dhist_alt.at(ihist),2);
      tfhist_alt.push_back((TH1*) nhist_alt.at(ihist)->Clone(Form("%s_temp",nhist_alt.at(ihist)->GetName())));
      tfhist_alt.back()->Divide(dhist_alt.at(ihist));
    }
  }


  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
    if(nhist_2D_alt.size() >= ihist){
      smoothEmptyBins(nhist_2D_alt.at(ihist),2);
      smoothEmptyBins(dhist_2D_alt.at(ihist),2);
      tfhist_2D_alt.push_back((TH2*) nhist_2D_alt.at(ihist)->Clone(Form("%s_temp",nhist_2D_alt.at(ihist)->GetName())));
      tfhist_2D_alt.back()->Divide(dhist_2D_alt.at(ihist));
    }
  }


  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);
  for(size_t ihist = 0; ihist < tfhist_alt.size(); ihist++)
    smoothEmptyBins(tfhist_alt.at(ihist),2);
  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);
  for(size_t ihist = 0; ihist < tfhist_2D_alt.size(); ihist++)
    smoothEmptyBins(tfhist_2D_alt.at(ihist),1);

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
    smoothEmptyBins(nhist.at(ihist),2);
    smoothEmptyBins(dhist.at(ihist),2);
    tfhist.push_back((TH1*) nhist.at(ihist)->Clone(Form("%s_temp",nhist.at(ihist)->GetName())));
    tfhist.back()->Divide(dhist.at(ihist));
    if(nhist_alt.size() >= ihist){
      smoothEmptyBins(nhist_alt.at(ihist),2);
      smoothEmptyBins(dhist_alt.at(ihist),2);
      tfhist_alt.push_back((TH1*) nhist_alt.at(ihist)->Clone(Form("%s_temp",nhist_alt.at(ihist)->GetName())));
      tfhist_alt.back()->Divide(dhist_alt.at(ihist));
    }
  }


  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    smoothEmptyBins(nhist_2D.at(ihist),2);
    smoothEmptyBins(dhist_2D.at(ihist),2);
    tfhist_2D.push_back((TH2*) nhist_2D.at(ihist)->Clone(Form("%s_temp",nhist_2D.at(ihist)->GetName())));
    tfhist_2D.back()->Divide(dhist_2D.at(ihist));
    if(nhist_2D_alt.size() >= ihist){
      smoothEmptyBins(nhist_2D_alt.at(ihist),2);
      smoothEmptyBins(dhist_2D_alt.at(ihist),2);
      tfhist_2D_alt.push_back((TH2*) nhist_2D_alt.at(ihist)->Clone(Form("%s_temp",nhist_2D_alt.at(ihist)->GetName())));
      tfhist_2D_alt.back()->Divide(dhist_2D_alt.at(ihist));
    }
  }


  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < tfhist.size(); ihist++)
    smoothEmptyBins(tfhist.at(ihist),2);
  for(size_t ihist = 0; ihist < tfhist_alt.size(); ihist++)
    smoothEmptyBins(tfhist_alt.at(ihist),2);
  for(size_t ihist = 0; ihist < tfhist_2D.size(); ihist++)
    smoothEmptyBins(tfhist_2D.at(ihist),1);
  for(size_t ihist = 0; ihist < tfhist_2D_alt.size(); ihist++)
    smoothEmptyBins(tfhist_2D_alt.at(ihist),1);

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
