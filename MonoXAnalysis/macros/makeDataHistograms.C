#include "makehist.h"

using namespace std;


// Build templates for the signal region                                                                                                                                       
void sigdatamchist(TFile* outfile,
                   int category,
                   string kFactorFile,
                   vector<string> observables,
                   double lumi              = 2.24,
                   bool applyQGLReweight    = false,
		   bool doShapeSystematics  = false,
		   bool doAlternativeTop   = false,
                   string interaction       = "Vector",
                   vector<pair<string,string> > massPoint = {make_pair("100","1")},
                   bool blind = false) {

  // Files for Znunu, Wlnu, Zll, top, qcd , diboson, signal, data                                                                                                            
  TFile* znfile  = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root");
  TFile* wlfile  = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root");
  TFile* zlfile  = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/sigfilter/sig_tree_DYJetsToLL_M-50.root");
  TFile* ttfile  = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root");
  TFile* qcdfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/sigfilter/sig_tree_QCD.root");
  TFile* dbfile  = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/sigfilter/sig_tree_DiBoson.root");
  
  // additional top sample
  TFile* ttfile_alt  = NULL;
  if(doAlternativeTop)
    ttfile_alt = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root");

  vector<TFile*> monoJfile;
  vector<TFile*> monoWfile;
  vector<TFile*> monoZfile;

  if(interaction == "Vector"){
    for(size_t iPoint = 0; iPoint < massPoint.size(); iPoint++){
      monoJfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMV_Vector/sigfilter/sig_tree_DMV_NNPDF30_Vector_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Vector/sigfilter/sig_tree_VectorMonoZ_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
    }
  }

  else if(interaction == "Axial"){
    for(size_t iPoint = 0; iPoint < massPoint.size(); iPoint++){
      monoJfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMV_Axial/sigfilter/sig_tree_DMV_NNPDF30_Axial_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Axial/sigfilter/sig_tree_AxialMonoZ_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
    }
  }
  else if(interaction == "Scalar"){
    for(size_t iPoint = 0; iPoint < massPoint.size(); iPoint++){
      monoJfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMS_Scalar/sigfilter/sig_tree_DMS_NNPDF30_Scalar_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+massPoint.at(iPoint).first+
					"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));
    }
  }
  else if(interaction == "Pseudoscalar"){
    for(size_t iPoint = 0; iPoint < massPoint.size(); iPoint++){
      monoJfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMS_Pseudoscalar/sigfilter/sig_tree_DMS_NNPDF30_Pseudoscalar_Mphi-"+
					massPoint.at(iPoint).first+"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-"+
					massPoint.at(iPoint).first+"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));
      monoZfile .push_back(TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarZH_Mphi-"+
					massPoint.at(iPoint).first+"_Mchi-"+massPoint.at(iPoint).second+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));
    }
  }


  //data                                                                                                                                                                        
  TFile* dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/sigfilter/sig_tree_crab_MET-Run2015.root");

  // make met histograms                                                                                                                                                        
  vector<TH1*> znhist;
  vector<TH1*> wlhist;
  vector<TH1*> zlhist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_bUp;
  vector<TH1*> tthist_bDw;
  vector<TH1*> tthist_metJetUp;
  vector<TH1*> tthist_metJetDw;
  vector<TH1*> tthist_metResUp;
  vector<TH1*> tthist_metResDw;
  vector<TH1*> tthist_metUncUp;
  vector<TH1*> tthist_metUncDw;

  vector<TH1*> tthist_alt;
  vector<TH1*> tthist_alt_bUp;
  vector<TH1*> tthist_alt_bDw;
  vector<TH1*> tthist_alt_metJetUp;
  vector<TH1*> tthist_alt_metJetDw;
  vector<TH1*> tthist_alt_metResUp;
  vector<TH1*> tthist_alt_metResDw;
  vector<TH1*> tthist_alt_metUncUp;
  vector<TH1*> tthist_alt_metUncDw;

  vector<TH1*> dihist;
  vector<TH1*> dihist_bUp;
  vector<TH1*> dihist_bDw;
  vector<TH1*> dihist_metJetUp;
  vector<TH1*> dihist_metJetDw;
  vector<TH1*> dihist_metResUp;
  vector<TH1*> dihist_metResDw;
  vector<TH1*> dihist_metUncUp;
  vector<TH1*> dihist_metUncDw;

  vector<TH1*> qcdhist;

  vector< vector<TH1*> > monoJhist;
  vector< vector<TH1*> > monoWhist;
  vector< vector<TH1*> > monoZhist;
  vector< vector<TH1*> > monoJhist_bUp;
  vector< vector<TH1*> > monoWhist_bUp;
  vector< vector<TH1*> > monoZhist_bUp;
  vector< vector<TH1*> > monoJhist_bDw;
  vector< vector<TH1*> > monoWhist_bDw;
  vector< vector<TH1*> > monoZhist_bDw;
  vector< vector<TH1*> > monoJhist_metJetUp;
  vector< vector<TH1*> > monoWhist_metJetUp;
  vector< vector<TH1*> > monoZhist_metJetUp;
  vector< vector<TH1*> > monoJhist_metJetDw;
  vector< vector<TH1*> > monoWhist_metJetDw;
  vector< vector<TH1*> > monoZhist_metJetDw;
  vector< vector<TH1*> > monoJhist_metResUp;
  vector< vector<TH1*> > monoWhist_metResUp;
  vector< vector<TH1*> > monoZhist_metResUp;
  vector< vector<TH1*> > monoJhist_metResDw;
  vector< vector<TH1*> > monoWhist_metResDw;
  vector< vector<TH1*> > monoZhist_metResDw;
  vector< vector<TH1*> > monoJhist_metUncUp;
  vector< vector<TH1*> > monoWhist_metUncUp;
  vector< vector<TH1*> > monoZhist_metUncUp;
  vector< vector<TH1*> > monoJhist_metUncDw;
  vector< vector<TH1*> > monoWhist_metUncDw;
  vector< vector<TH1*> > monoZhist_metUncDw;

  vector<TH1*> dthist;

  vector<float> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    TH1F* znhist_temp = new TH1F(("zinvhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* wlhist_temp = new TH1F(("wjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* zlhist_temp = new TH1F(("zjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* tthist_temp = new TH1F(("tbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dihist_temp = new TH1F(("dbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qcdhist_temp = new TH1F(("qbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dthist_temp = new TH1F(("datahist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    znhist.push_back(dynamic_cast<TH1*>(znhist_temp));
    wlhist.push_back(dynamic_cast<TH1*>(wlhist_temp));
    zlhist.push_back(dynamic_cast<TH1*>(zlhist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
    dihist.push_back(dynamic_cast<TH1*>(dihist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));

    if(doAlternativeTop){
      TH1F* tthist_alt_temp = new TH1F(("tbkghist_alt_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }

    if(doShapeSystematics){

      // b-tagging
      TH1F* tthist_bUp_temp = new TH1F(("tbkghist_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_bDw_temp = new TH1F(("tbkghist_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_bUp_temp = new TH1F(("dbkghist_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_bDw_temp = new TH1F(("dbkghist_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      tthist_bUp.push_back(dynamic_cast<TH1*>(tthist_bUp_temp));
      tthist_bDw.push_back(dynamic_cast<TH1*>(tthist_bDw_temp));
      dihist_bUp.push_back(dynamic_cast<TH1*>(dihist_bUp_temp));
      dihist_bDw.push_back(dynamic_cast<TH1*>(dihist_bDw_temp));

      TH1F* tthist_metJetUp_temp = new TH1F(("tbkghist_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metJetDw_temp = new TH1F(("tbkghist_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metResUp_temp = new TH1F(("tbkghist_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metResDw_temp = new TH1F(("tbkghist_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metUncUp_temp = new TH1F(("tbkghist_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metUncDw_temp = new TH1F(("tbkghist_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      tthist_metJetUp.push_back(dynamic_cast<TH1*>(tthist_metJetUp_temp));
      tthist_metJetDw.push_back(dynamic_cast<TH1*>(tthist_metJetDw_temp));
      tthist_metResUp.push_back(dynamic_cast<TH1*>(tthist_metResUp_temp));
      tthist_metResDw.push_back(dynamic_cast<TH1*>(tthist_metResDw_temp));
      tthist_metUncUp.push_back(dynamic_cast<TH1*>(tthist_metUncUp_temp));
      tthist_metUncDw.push_back(dynamic_cast<TH1*>(tthist_metUncDw_temp));

      TH1F* dihist_metJetUp_temp = new TH1F(("dbkghist_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_metJetDw_temp = new TH1F(("dbkghist_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_metResUp_temp = new TH1F(("dbkghist_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_metResDw_temp = new TH1F(("dbkghist_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_metUncUp_temp = new TH1F(("dbkghist_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_metUncDw_temp = new TH1F(("dbkghist_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dihist_metJetUp.push_back(dynamic_cast<TH1*>(dihist_metJetUp_temp));
      dihist_metJetDw.push_back(dynamic_cast<TH1*>(dihist_metJetDw_temp));
      dihist_metResUp.push_back(dynamic_cast<TH1*>(dihist_metResUp_temp));
      dihist_metResDw.push_back(dynamic_cast<TH1*>(dihist_metResDw_temp));
      dihist_metUncUp.push_back(dynamic_cast<TH1*>(dihist_metUncUp_temp));
      dihist_metUncDw.push_back(dynamic_cast<TH1*>(dihist_metUncDw_temp));
      
      if(doAlternativeTop){

	TH1F* tthist_alt_bUp_temp = new TH1F(("tbkghist_alt_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_bDw_temp = new TH1F(("tbkghist_alt_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	tthist_alt_bUp.push_back(dynamic_cast<TH1*>(tthist_alt_bUp_temp));
	tthist_alt_bDw.push_back(dynamic_cast<TH1*>(tthist_alt_bDw_temp));	

	TH1F* tthist_alt_metJetUp_temp = new TH1F(("tbkghist_alt_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_metJetDw_temp = new TH1F(("tbkghist_alt_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_metResUp_temp = new TH1F(("tbkghist_alt_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_metResDw_temp = new TH1F(("tbkghist_alt_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_metUncUp_temp = new TH1F(("tbkghist_alt_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* tthist_alt_metUncDw_temp = new TH1F(("tbkghist_alt_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	
	tthist_alt_metJetUp.push_back(dynamic_cast<TH1*>(tthist_alt_metJetUp_temp));
	tthist_alt_metJetDw.push_back(dynamic_cast<TH1*>(tthist_alt_metJetDw_temp));
	tthist_alt_metResUp.push_back(dynamic_cast<TH1*>(tthist_alt_metResUp_temp));
	tthist_alt_metResDw.push_back(dynamic_cast<TH1*>(tthist_alt_metResDw_temp));
	tthist_alt_metUncUp.push_back(dynamic_cast<TH1*>(tthist_alt_metUncUp_temp));
	tthist_alt_metUncDw.push_back(dynamic_cast<TH1*>(tthist_alt_metUncDw_temp));

      }
    }
  }
  // signal                                                                                                                                                                 
  monoJhist.assign(massPoint.size(),vector<TH1*>());
  monoWhist.assign(massPoint.size(),vector<TH1*>());
  monoZhist.assign(massPoint.size(),vector<TH1*>());
  if(doShapeSystematics){
    monoJhist_bUp.assign(massPoint.size(),vector<TH1*>());
    monoWhist_bUp.assign(massPoint.size(),vector<TH1*>());
    monoZhist_bUp.assign(massPoint.size(),vector<TH1*>());
    monoJhist_bDw.assign(massPoint.size(),vector<TH1*>());
    monoWhist_bDw.assign(massPoint.size(),vector<TH1*>());
    monoZhist_bDw.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metJetUp.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metJetUp.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metJetUp.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metJetDw.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metJetDw.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metJetDw.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metResUp.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metResUp.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metResUp.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metResDw.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metResDw.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metResDw.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metUncUp.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metUncUp.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metUncUp.assign(massPoint.size(),vector<TH1*>());
    monoJhist_metUncDw.assign(massPoint.size(),vector<TH1*>());
    monoWhist_metUncDw.assign(massPoint.size(),vector<TH1*>());
    monoZhist_metUncDw.assign(massPoint.size(),vector<TH1*>());
  }

  int imass = 0;

  for( auto massp : massPoint){
    for(auto obs : observables){

      bins = selectBinning(obs,category);
      if(bins.empty())
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* monoJhist_temp = new TH1F(("monoJhist_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoWhist_temp = new TH1F(("monoWhist_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoZhist_temp = new TH1F(("monoZhist_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      monoJhist.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_temp));
      monoWhist.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_temp));
      monoZhist.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_temp));

      if(doShapeSystematics){

	TH1F* monoJhist_bUp_temp = new TH1F(("monoJhist_bUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_bUp_temp = new TH1F(("monoWhist_bUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_bUp_temp = new TH1F(("monoZhist_bUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bUp_temp));
	monoWhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bUp_temp));
	monoZhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bUp_temp));

	TH1F* monoJhist_bDw_temp = new TH1F(("monoJhist_bDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_bDw_temp = new TH1F(("monoWhist_bDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_bDw_temp = new TH1F(("monoZhist_bDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bDw_temp));
	monoWhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bDw_temp));
	monoZhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bDw_temp));


	TH1F* monoJhist_metJetUp_temp = new TH1F(("monoJhist_metJetUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metJetUp_temp = new TH1F(("monoWhist_metJetUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metJetUp_temp = new TH1F(("monoZhist_metJetUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetUp_temp));
	monoWhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetUp_temp));
	monoZhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetUp_temp));

	TH1F* monoJhist_metJetDw_temp = new TH1F(("monoJhist_metJetDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metJetDw_temp = new TH1F(("monoWhist_metJetDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metJetDw_temp = new TH1F(("monoZhist_metJetDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetDw_temp));
	monoWhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetDw_temp));
	monoZhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetDw_temp));

	TH1F* monoJhist_metResUp_temp = new TH1F(("monoJhist_metResUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metResUp_temp = new TH1F(("monoWhist_metResUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metResUp_temp = new TH1F(("monoZhist_metResUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResUp_temp));
	monoWhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResUp_temp));
	monoZhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResUp_temp));

	TH1F* monoJhist_metResDw_temp = new TH1F(("monoJhist_metResDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metResDw_temp = new TH1F(("monoWhist_metResDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metResDw_temp = new TH1F(("monoZhist_metResDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResDw_temp));
	monoWhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResDw_temp));
	monoZhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResDw_temp));

	TH1F* monoJhist_metUncUp_temp = new TH1F(("monoJhist_metUncUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metUncUp_temp = new TH1F(("monoWhist_metUncUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metUncUp_temp = new TH1F(("monoZhist_metUncUp_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metUncUp_temp));
	monoWhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metUncUp_temp));
	monoZhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metUncUp_temp));

	TH1F* monoJhist_metUncDw_temp = new TH1F(("monoJhist_metUncDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metUncDw_temp = new TH1F(("monoWhist_metUncDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metUncDw_temp = new TH1F(("monoZhist_metUncDw_"+massp.first+"_"+massp.second+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metUncDw_temp));
	monoWhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metUncDw_temp));
	monoZhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metUncDw_temp));

      }

    }
    imass++;
  }

  vector<TH2*> znhist_2D;
  vector<TH2*> wlhist_2D;
  vector<TH2*> zlhist_2D;
  vector<TH2*> tthist_2D;

  vector<TH2*> tthist_bUp_2D;
  vector<TH2*> tthist_bDw_2D;
  vector<TH2*> tthist_metJetUp_2D;
  vector<TH2*> tthist_metJetDw_2D;
  vector<TH2*> tthist_metResUp_2D;
  vector<TH2*> tthist_metResDw_2D;
  vector<TH2*> tthist_metUncUp_2D;
  vector<TH2*> tthist_metUncDw_2D;

  vector<TH2*> tthist_alt_2D;
  vector<TH2*> tthist_alt_bUp_2D;
  vector<TH2*> tthist_alt_bDw_2D;
  vector<TH2*> tthist_alt_metJetUp_2D;
  vector<TH2*> tthist_alt_metJetDw_2D;
  vector<TH2*> tthist_alt_metResUp_2D;
  vector<TH2*> tthist_alt_metResDw_2D;
  vector<TH2*> tthist_alt_metUncUp_2D;
  vector<TH2*> tthist_alt_metUncDw_2D;

  vector<TH2*> dihist_2D;
  vector<TH2*> dihist_bUp_2D;
  vector<TH2*> dihist_bDw_2D;
  vector<TH2*> dihist_metJetUp_2D;
  vector<TH2*> dihist_metJetDw_2D;
  vector<TH2*> dihist_metResUp_2D;
  vector<TH2*> dihist_metResDw_2D;
  vector<TH2*> dihist_metUncUp_2D;
  vector<TH2*> dihist_metUncDw_2D;

  vector<TH2*> qcdhist_2D;

  vector< vector<TH2*> > monoJhist_2D;
  vector< vector<TH2*> > monoWhist_2D;
  vector< vector<TH2*> > monoZhist_2D;
  vector< vector<TH2*> > monoJhist_bUp_2D;
  vector< vector<TH2*> > monoWhist_bUp_2D;
  vector< vector<TH2*> > monoZhist_bUp_2D;
  vector< vector<TH2*> > monoJhist_bDw_2D;
  vector< vector<TH2*> > monoWhist_bDw_2D;
  vector< vector<TH2*> > monoZhist_bDw_2D;

  vector< vector<TH2*> > monoJhist_metJetUp_2D;
  vector< vector<TH2*> > monoWhist_metJetUp_2D;
  vector< vector<TH2*> > monoZhist_metJetUp_2D;
  vector< vector<TH2*> > monoJhist_metJetDw_2D;
  vector< vector<TH2*> > monoWhist_metJetDw_2D;
  vector< vector<TH2*> > monoZhist_metJetDw_2D;
  vector< vector<TH2*> > monoJhist_metResUp_2D;
  vector< vector<TH2*> > monoWhist_metResUp_2D;
  vector< vector<TH2*> > monoZhist_metResUp_2D;
  vector< vector<TH2*> > monoJhist_metResDw_2D;
  vector< vector<TH2*> > monoWhist_metResDw_2D;
  vector< vector<TH2*> > monoZhist_metResDw_2D;
  vector< vector<TH2*> > monoJhist_metUncUp_2D;
  vector< vector<TH2*> > monoWhist_metUncUp_2D;
  vector< vector<TH2*> > monoZhist_metUncUp_2D;
  vector< vector<TH2*> > monoJhist_metUncDw_2D;
  vector< vector<TH2*> > monoWhist_metUncDw_2D;
  vector< vector<TH2*> > monoZhist_metUncDw_2D;

  monoJhist_2D.assign(massPoint.size(),vector<TH2*>());
  monoWhist_2D.assign(massPoint.size(),vector<TH2*>());
  monoZhist_2D.assign(massPoint.size(),vector<TH2*>());

  if(doShapeSystematics){
    monoJhist_bUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_bUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_bUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoJhist_bDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_bDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_bDw_2D.assign(massPoint.size(),vector<TH2*>());

    monoJhist_metJetUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metJetUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metJetUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoJhist_metJetDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metJetDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metJetDw_2D.assign(massPoint.size(),vector<TH2*>());

    monoJhist_metResUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metResUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metResUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoJhist_metResDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metResDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metResDw_2D.assign(massPoint.size(),vector<TH2*>());

    monoJhist_metUncUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metUncUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metUncUp_2D.assign(massPoint.size(),vector<TH2*>());
    monoJhist_metUncDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoWhist_metUncDw_2D.assign(massPoint.size(),vector<TH2*>());
    monoZhist_metUncDw_2D.assign(massPoint.size(),vector<TH2*>());
  }


  vector<TH2*> dthist_2D;

  TTree* zntree = (TTree*)znfile->Get("tree/tree");
  TTree* wltree = (TTree*)wlfile->Get("tree/tree");
  TTree* zltree = (TTree*)zlfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* tttree_alt = NULL;
  if(ttfile_alt and doAlternativeTop)
    tttree_alt = (TTree*)ttfile_alt->Get("tree/tree");

  TTree* ditree = (TTree*)dbfile->Get("tree/tree");
  TTree* qcdtree = (TTree*)qcdfile->Get("tree/tree");
  vector<TTree* > monoJtree;
  vector<TTree* > monoWtree;
  vector<TTree* > monoZtree;

  for(auto file : monoJfile){
    if(file)
      monoJtree.push_back((TTree*) file->Get("tree/tree"));
    else
      monoJtree.push_back(0);
  }

  for(auto file : monoWfile){
    if(file)
      monoWtree.push_back((TTree*) file->Get("tree/tree"));
    else
      monoWtree.push_back(0);
  }

  for(auto file : monoZfile){
    if(file)
      monoZtree.push_back((TTree*) file->Get("tree/tree"));
    else
      monoZtree.push_back(0);
  }


  TTree* dttree = (TTree*)dtfile->Get("tree");

  // get k-factors NLO                                                                                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  // NLO corrections for Z and W base mc                                                                                                                                       
  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> whists;
  // apply EWK and QCD corrections                                                                                                                                              
  zhists.push_back(znlohist); zhists.push_back(zewkhist);
  whists.push_back(wnlohist); whists.push_back(wewkhist);

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  // make histograms for the signal region                                                                                                                                      
  int QGLZ_index = 1;
  int QGLW_index = 2;
  int QGLG_index = 3;
  int QGLT_index = 4;

  if(not applyQGLReweight){
    QGLZ_index = 0;
    QGLW_index = 0;
    QGLG_index = 0;
    QGLT_index = 0;
  }

  cout<<"signal region: Z->nunu sample "<<endl;
  makehist4(zntree, znhist,  znhist_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "", true, NULL);
  cout<<"signal region: W+jets sample "<<endl;
  makehist4(wltree, wlhist,  wlhist_2D,  true, 0, category, false, 1.00, lumi,    QGLW_index, whists, "", true, NULL);
  cout<<"signal region: Z+jets sample "<<endl;
  makehist4(zltree, zlhist,  zlhist_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "", true, NULL);
  cout<<"signal region: TTbar sample "<<endl;
  makehist4(tttree, tthist,  tthist_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "", true, NULL);

    //alternative ttbar             
  if(doAlternativeTop){
    cout<<"signal region: TTbar alternative sample "<<endl;
    makehist4(tttree_alt, tthist_alt,  tthist_alt_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "", true, NULL);
  }

  cout<<"signal region: Diboson sample "<<endl;
  makehist4(ditree, dihist,  dihist_2D,  true, 0, category, isWJet, 1.00, lumi,   0, ehists, "", true, NULL);
  cout<<"signal region: QCD sample "<<endl;
  makehist4(qcdtree, qcdhist,  qcdhist_2D,  true, 0, category, false, 1.00, lumi, 0,  ehists, "", true, NULL);

  if(doShapeSystematics){
    cout<<"signal region analysis --> do top shape sys "<<endl;
    makehist4(tttree, tthist_bUp,  tthist_bUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagUp", true, NULL);
    makehist4(tttree, tthist_bDw,  tthist_bDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagDw", true, NULL);
    makehist4(tttree, tthist_metJetUp,  tthist_metJetUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesUp", true, NULL);
    makehist4(tttree, tthist_metJetDw,  tthist_metJetDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesDw", true, NULL);
    makehist4(tttree, tthist_metResUp,  tthist_metResUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerUp", true, NULL);
    makehist4(tttree, tthist_metResDw,  tthist_metResDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerDw", true, NULL);
    makehist4(tttree, tthist_metUncUp,  tthist_metUncUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncUp", true, NULL);
    makehist4(tttree, tthist_metUncDw,  tthist_metUncDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncDw", true, NULL);
    
    if(doAlternativeTop){
      cout<<"signal region analysis --> do top alternative shape sys "<<endl;
      makehist4(tttree_alt, tthist_alt_bUp,       tthist_alt_bUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagUp", true, NULL);
      makehist4(tttree_alt, tthist_alt_bDw,       tthist_alt_bDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagDw", true, NULL);
      makehist4(tttree_alt, tthist_alt_metJetUp,  tthist_alt_metJetUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesUp", true, NULL);
      makehist4(tttree_alt, tthist_alt_metJetDw,  tthist_alt_metJetDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesDw", true, NULL);
      makehist4(tttree_alt, tthist_alt_metResUp,  tthist_alt_metResUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerUp", true, NULL);
      makehist4(tttree_alt, tthist_alt_metResDw,  tthist_alt_metResDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerDw", true, NULL);
      makehist4(tttree_alt, tthist_alt_metUncUp,  tthist_alt_metUncUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncUp", true, NULL);
      makehist4(tttree_alt, tthist_alt_metUncDw,  tthist_alt_metUncDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncDw", true, NULL);
    }


    cout<<"signal region analysis --> do diboson shape sys "<<endl;
    makehist4(tttree, dihist_bUp,  dihist_bUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "btagUp", true, NULL);
    makehist4(tttree, dihist_bDw,  dihist_bDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "btagDw", true, NULL);
    makehist4(tttree, dihist_metJetUp,  dihist_metJetUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jesUp", true, NULL);
    makehist4(tttree, dihist_metJetDw,  dihist_metJetDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jesDw", true, NULL);
    makehist4(tttree, dihist_metResUp,  dihist_metResUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jerUp", true, NULL);
    makehist4(tttree, dihist_metResDw,  dihist_metResDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jerDw", true, NULL);
    makehist4(tttree, dihist_metUncUp,  dihist_metUncUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "uncUp", true, NULL);
    makehist4(tttree, dihist_metUncDw,  dihist_metUncDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "uncDw", true, NULL);
    
  }

  // take average of ttbar                                                                                                                                                     
  if(doAlternativeTop){

    for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
      makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));    
    for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
      makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));

    if(doShapeSystematics){
      for(size_t iHisto = 0; iHisto < tthist_bUp.size(); iHisto++)
	makeAverage(tthist_bUp.at(iHisto),tthist_alt_bUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_bDw.size(); iHisto++)
	makeAverage(tthist_bDw.at(iHisto),tthist_alt_bDw.at(iHisto)); 
      for(size_t iHisto = 0; iHisto < tthist_metJetUp.size(); iHisto++)
	makeAverage(tthist_metJetUp.at(iHisto),tthist_alt_metJetUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metJetDw.size(); iHisto++)
	makeAverage(tthist_metJetDw.at(iHisto),tthist_alt_metJetDw.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metResUp.size(); iHisto++)
	makeAverage(tthist_metResUp.at(iHisto),tthist_alt_metResUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metResDw.size(); iHisto++)
	makeAverage(tthist_metResDw.at(iHisto),tthist_alt_metResDw.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metUncUp.size(); iHisto++)
	makeAverage(tthist_metUncUp.at(iHisto),tthist_alt_metUncUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metUncDw.size(); iHisto++)
	makeAverage(tthist_metUncDw.at(iHisto),tthist_alt_metUncDw.at(iHisto));	
    }
  }

  // Signals
  int itree = 0;
  for(auto tree : monoJtree){
    cout<<"signal region analysis --> Signal monoJet "<<endl;
    // signals                                                                                                                                                                 
    if(tree){
      makehist4(tree, monoJhist.at(itree),  monoJhist_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "",  true, NULL);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoJ sys "<<endl;
	makehist4(tree, monoJhist_bUp.at(itree),  monoJhist_bUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bUp",  true, NULL);
	makehist4(tree, monoJhist_bDw.at(itree),  monoJhist_bDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bDw",  true, NULL);
	makehist4(tree, monoJhist_metJetUp.at(itree),  monoJhist_metJetUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecUp",  true, NULL);
	makehist4(tree, monoJhist_metJetDw.at(itree),  monoJhist_metJetDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecDw",  true, NULL);
	makehist4(tree, monoJhist_metResUp.at(itree),  monoJhist_metResUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerUp",  true, NULL);
	makehist4(tree, monoJhist_metResDw.at(itree),  monoJhist_metResDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerDw",  true, NULL);
	makehist4(tree, monoJhist_metUncUp.at(itree),  monoJhist_metUncUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncUp",  true, NULL);
	makehist4(tree, monoJhist_metUncDw.at(itree),  monoJhist_metUncDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncDw",  true, NULL);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoWtree){
    cout<<"signal region analysis --> Signal monoW "<<endl;
    if(tree){
      makehist4(tree, monoWhist.at(itree),  monoWhist_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "", true, NULL);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoW sys "<<endl;
	makehist4(tree, monoWhist_bUp.at(itree),  monoWhist_bUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bUp",  true, NULL);
	makehist4(tree, monoWhist_bDw.at(itree),  monoWhist_bDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bDw",  true, NULL);
	makehist4(tree, monoWhist_metJetUp.at(itree),  monoWhist_metJetUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecUp",  true, NULL);
	makehist4(tree, monoWhist_metJetDw.at(itree),  monoWhist_metJetDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecDw",  true, NULL);
	makehist4(tree, monoWhist_metResUp.at(itree),  monoWhist_metResUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerUp",  true, NULL);
	makehist4(tree, monoWhist_metResDw.at(itree),  monoWhist_metResDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerDw",  true, NULL);
	makehist4(tree, monoWhist_metUncUp.at(itree),  monoWhist_metUncUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncUp",  true, NULL);
	makehist4(tree, monoWhist_metUncDw.at(itree),  monoWhist_metUncDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncDw",  true, NULL);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoZtree){
    cout<<"signal region analysis --> Signal monoZ "<<endl;
    if(tree){
      makehist4(tree, monoZhist.at(itree),  monoZhist_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "", true, NULL);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoZ sys "<<endl;
	makehist4(tree, monoZhist_bUp.at(itree),  monoZhist_bUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bUp",  true, NULL);
	makehist4(tree, monoZhist_bDw.at(itree),  monoZhist_bDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "bDw",  true, NULL);
	makehist4(tree, monoZhist_metJetUp.at(itree),  monoZhist_metJetUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecUp",  true, NULL);
	makehist4(tree, monoZhist_metJetDw.at(itree),  monoZhist_metJetDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jecDw",  true, NULL);
	makehist4(tree, monoZhist_metResUp.at(itree),  monoZhist_metResUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerUp",  true, NULL);
	makehist4(tree, monoZhist_metResDw.at(itree),  monoZhist_metResDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerDw",  true, NULL);
	makehist4(tree, monoZhist_metUncUp.at(itree),  monoZhist_metUncUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncUp",  true, NULL);
	makehist4(tree, monoZhist_metUncDw.at(itree),  monoZhist_metUncDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncDw",  true, NULL);
      }
    }
    itree++;
  }

  // data                                                
  cout<<"signal region analysis --> loop on data "<<endl;
  makehist4(dttree, dthist,  dthist_2D, false, 0, category, false, 1.00, lumi, 0, ehists, "", true, NULL);
    
  if (blind) {
    for( size_t ihist = 0; ihist < dthist.size(); ihist++){
      for (int i = 1; i <= dthist.at(ihist)->GetNbinsX(); i++) {
        double binval = 0.0;
        binval += znhist.at(ihist)->GetBinContent(i);
	binval += wlhist.at(ihist)->GetBinContent(i);
        binval += zlhist.at(ihist)->GetBinContent(i);
        binval += tthist.at(ihist)->GetBinContent(i);
        binval += dihist.at(ihist)->GetBinContent(i);
        binval += qcdhist.at(ihist)->GetBinContent(i);
        dthist.at(ihist)->SetBinContent(i, int(binval));
      }
    }

    for( size_t ihist = 0; ihist < dthist_2D.size(); ihist++){
      for (int iX = 1; iX <= dthist_2D.at(ihist)->GetNbinsX(); iX++) {
        for (int iY = 1; iY <= dthist_2D.at(ihist)->GetNbinsY(); iY++) {
          double binval = 0.0;
          binval += znhist_2D.at(ihist)->GetBinContent(iX,iY);
          binval += wlhist_2D.at(ihist)->GetBinContent(iX,iY);
          binval += zlhist_2D.at(ihist)->GetBinContent(iX,iY);
          binval += tthist_2D.at(ihist)->GetBinContent(iX,iY);
          binval += dihist_2D.at(ihist)->GetBinContent(iX,iY);
          binval += qcdhist_2D.at(ihist)->GetBinContent(iX,iY);
          dthist_2D.at(ihist)->SetBinContent(iX,iY,int(binval));
        }
      }
    }
  }

  
  outfile->cd();
  // store histograms                                                                                                                                                        
  for(auto hist : znhist) hist->Write();
  for(auto hist : wlhist) hist->Write();
  for(auto hist : zlhist) hist->Write();
  for(auto hist : tthist) hist->Write();
  for(auto hist : dihist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto sample : monoJhist){ for(auto hist : sample) hist->Write();}
  for(auto sample : monoWhist){ for(auto hist : sample) hist->Write();}  
  for(auto sample : monoZhist){ for(auto hist : sample) hist->Write();}
  for(auto hist : dthist) hist->Write();

  //
  if(doShapeSystematics){
    for(auto hist : tthist_bUp) hist->Write();
    for(auto hist : tthist_bDw) hist->Write();
    for(auto hist : tthist_metJetUp) hist->Write();
    for(auto hist : tthist_metJetDw) hist->Write();
    for(auto hist : tthist_metResUp) hist->Write();
    for(auto hist : tthist_metResDw) hist->Write();
    for(auto hist : tthist_metUncUp) hist->Write();
    for(auto hist : tthist_metUncDw) hist->Write();

    for(auto hist : dihist_bUp) hist->Write();
    for(auto hist : dihist_bDw) hist->Write();
    for(auto hist : dihist_metJetUp) hist->Write();
    for(auto hist : dihist_metJetDw) hist->Write();
    for(auto hist : dihist_metResUp) hist->Write();
    for(auto hist : dihist_metResDw) hist->Write();
    for(auto hist : dihist_metUncUp) hist->Write();
    for(auto hist : dihist_metUncDw) hist->Write();
    
    for(auto sample : monoJhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncDw){ for(auto hist : sample) hist->Write();}

    for(auto sample : monoWhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncDw){ for(auto hist : sample) hist->Write();}

    for(auto sample : monoZhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncDw){ for(auto hist : sample) hist->Write();}
  }

  // store hist_2Dograms                                                                                                                                                     
  for(auto hist_2D : znhist_2D) hist_2D->Write();
  for(auto hist_2D : wlhist_2D) hist_2D->Write();
  for(auto hist_2D : zlhist_2D) hist_2D->Write();
  for(auto hist_2D : tthist_2D) hist_2D->Write();
  for(auto hist_2D : dihist_2D) hist_2D->Write();
  for(auto hist_2D : qcdhist_2D) hist_2D->Write();
  for(auto sample : monoJhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}
  for(auto sample : monoWhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}
  for(auto sample : monoZhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}
  for(auto hist_2D : dthist_2D) hist_2D->Write();

  znfile->Close();
  wlfile->Close();
  zlfile->Close();
  ttfile->Close();
  dbfile->Close();
  qcdfile->Close();
  for(auto file : monoJfile){ if(file) file->Close();}
  for(auto file : monoWfile){ if(file) file->Close();}
  for(auto file : monoZfile){ if(file) file->Close();}
  dtfile->Close();
  kffile.Close();
  if(ttfile_alt) ttfile_alt->Close();

  cout << "Templates for the signal region computed ..." << endl;
}


// build templates for photon+jets control region                                                                                                                           
void gamdatamchist(TFile* outfile,
                   string kFactorFile,
                   int category,
                   vector<string> observables,
                   double lumi           = 2.24,
                   bool applyQGLReweight = false
                   ) {


  TFile* dtfile =  TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SinglePhoton/gamfilter/gam_tree_crab_SinglePhoton-Run2015.root");
  TFile* gmfile =  TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root");

  vector<TH1*> dthist;
  vector<TH1*> qcdhist;
  vector<TH1*> gmhist;
  vector<TH2*> dthist_2D;
  vector<TH2*> qcdhist_2D;
  vector<TH2*> gmhist_2D;

  vector<float> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    TH1F* gmhist_temp = new TH1F(("gbkghistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qchist_temp = new TH1F(("qbkghistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dthist_temp = new TH1F(("datahistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    qcdhist.push_back(dynamic_cast<TH1*>(qchist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
  }

  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* gmtree = (TTree*)gmfile->Get("tree/tree");

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                           
  TFile kffile(kFactorFile.c_str());
  TH1* anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
  TH1*  alohist = (TH1*)kffile.Get("alo/alo_nominal");
  TH1* aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");

  vector<TH1*> ahists;
  vector<TH1*> ehists;

  anlohist->Divide(alohist);
  ahists.push_back(anlohist);
  ahists.push_back(aewkhist);

  if(applyQGLReweight){
    cout<<"gamma+jets control region --> data"<<endl;
    makehist4(dttree, dthist, dthist_2D, false, 5, category, false, 1.00, lumi, 0, ehists, "", true, NULL);
    cout<<"gamma+jets control region --> gamma+jets"<<endl;
    makehist4(gmtree, gmhist, gmhist_2D, true,  5, category, false, 1.00, lumi, 3, ahists, "", true, NULL);
    cout<<"gamma+jets control region: QCD"<<endl;
    makehist4(dttree, qcdhist, qcdhist_2D, false, 6, category, false, 1.00, lumi, 0, ehists, "", true, NULL);
  }
  else{
    cout<<"gamma+jets control region --> data"<<endl;
    makehist4(dttree, dthist, dthist_2D, false, 5, category, false, 1.00, lumi, 0, ehists, "", true, NULL);
    cout<<"gamma+jets control region --> gamma+jets"<<endl;
    makehist4(gmtree, gmhist, gmhist_2D, true,  5, category, false, 1.00, lumi, 0, ahists, "", true, NULL);
    cout<<"gamma+jets control region --> QCD"<<endl;
    makehist4(dttree, qcdhist, qcdhist_2D, false, 6, category, false, 1.00, lumi, 0, ehists, "", true, NULL);
  }
			    
  outfile->cd();

  for(auto hist : dthist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto hist : dthist_2D) hist->Write();
  for(auto hist : gmhist_2D) hist->Write();
  for(auto hist : qcdhist_2D) hist->Write();

  dtfile->Close();
  gmfile->Close();

  cout << "Templates for the gamma+jets control region computed ..." << endl;
}


//build templates for Zmumu, Zee, Wenu, Wmunu                                                                                                                                  
void lepdatamchist(TFile* outfile, string kFactorFile, int sample, int category, vector<string> observables, double lumi=2.24, bool applyQGLReweight=false) {

  if (sample != 1 && sample != 2 && sample != 3 && sample != 4) return;

  TFile* ttfile  = NULL;
  TFile* dbfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;

  string suffix;

  if(sample == 1){
    suffix = "zmm";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zmmfilter/zmm_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/zmmfilter/zmm_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/zmmfilter/zmm_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/zmmfilter/zmm_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/zmmfilter/zmm_tree_Top_amc.root");
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/zmmfilter/zmm_tree_crab_MET-Run2015.root");
  }
  else if(sample == 2){

    suffix = "wmn";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/wmnfilter/wmn_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wmnfilter/wmn_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/wmnfilter/wmn_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/wmnfilter/wmn_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/wmnfilter/wmn_tree_Top_amc.root");
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/wmnfilter/wmn_tree_crab_MET-Run2015.root");
  }
  else if(sample == 3){

    suffix = "zee";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zeefilter/zee_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/zeefilter/zee_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/zeefilter/zee_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/zeefilter/zee_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/zeefilter/zee_tree_Top_amc.root");
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/zeefilter/zee_tree_crab_SingleEle-Run2015.root");
  }
  else if(sample == 4){

    suffix = "wen";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/wenfilter/wen_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wenfilter/wen_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/wenfilter/wen_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/wenfilter/wen_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/wenfilter/wen_tree_Top_amc.root");
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/wenfilter/wen_tree_crab_SingleEle-Run2015.root");
  }


  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;

  vector<float> bins;
  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;


    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;

  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* vltree = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* dbtree = (TTree*)dbfile->Get("tree/tree");
  TTree* qctree = (TTree*)qcfile->Get("tree/tree");

  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> vllhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu                                                                                                                      
  vlhists.push_back(wnlohist);
  vlhists.push_back(wewkhist);
  vllhists.push_back(znlohist);
  vllhists.push_back(zewkhist);

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  if(applyQGLReweight){
    cout<<"lepton+jets control region --> top"<<endl;
    makehist4(tttree, tthist,  tthist_2D,  true,  sample, category, false,  1.00, lumi, 4, ehists,   "", true, NULL);
    cout<<"lepton+jets control region --> W+jets"<<endl;
    makehist4(vltree, vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, 2, vlhists,  "", true, NULL);
    cout<<"lepton+jets control region --> Z+jets"<<endl;
    makehist4(vlltree,vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, 1, vllhists, "", true, NULL);
  }
  else{
    cout<<"lepton+jets control region --> top"<<endl;
    makehist4(tttree, tthist,  tthist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists,   "", true, NULL);
    cout<<"lepton+jets control region --> W+jets"<<endl;
    makehist4(vltree, vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, 0, vlhists,  "", true, NULL);
    cout<<"lepton+jets control region --> Z+jets"<<endl;
    makehist4(vlltree,vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, 0, vllhists, "", true, NULL);
  }
  
  cout<<"lepton+jets control region --> Diboson"<<endl;
  makehist4(dbtree, dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, 0, ehists, "", true, NULL);
  cout<<"lepton+jets control region --> QCD"<<endl;
  makehist4(qctree, qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
  cout<<"lepton+jets control region --> Data"<<endl;
  makehist4(dttree, dthist, dthist_2D,   false, sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);

  outfile->cd();
  for(auto hist :  dthist) hist->Write();
  for(auto hist :  tthist) hist->Write();
  for(auto hist :  dbhist) hist->Write();
  for(auto hist :  qchist) hist->Write();
  for(auto hist :  vlhist) hist->Write();
  for(auto hist :  vllhist) hist->Write();

  dtfile->Close();
  vlfile->Close();
  vllfile->Close();
  ttfile->Close();
  dbfile->Close();
  qcfile->Close();

  cout << "Templates for the lepton control region computed ..." << endl;
}


//build template for top
void topdatamchist(TFile* outfile, string kFactorFile, int sample, int category, vector<string> observables, double lumi = 2.24, bool applyQGLReweight = false,
		   bool makeResonantSelection = false) {

  if (sample != 7 && sample != 8) return;

  TFile* ttfile      = NULL;
  TFile* ttfile_alt  = NULL;
  TFile* dbfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;

  string suffix;

  if(sample == 7)
    suffix = "topmu";
  else if(sample == 8)
    suffix = "topel";

  vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/topfilter/top_tree_DYJetsToLL_M-50.root");
  vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/topfilter/top_tree_WJetsToLNu.root");
  qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/topfilter/top_tree_QCD.root");
  dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/topfilter/top_tree_DiBoson.root");
  ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root");
  ttfile_alt = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root");

  if(sample == 7)
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/topfilter/top_tree_crab_MET-Run2015.root");
  else if(sample == 8)
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/topfilter/top_tree_crab_SingleEle-Run2015.root");

  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_alt;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;

  vector<TH1*> tthist_matched;
  vector<TH1*> tthist_matched_alt;
  vector<TH1*> tthist_unmatched;
  vector<TH1*> tthist_unmatched_alt;

  vector<float> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

    if(makeResonantSelection){
      TH1F* tthist_matched_temp = new TH1F((string("tbkghist_matched")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_matched_alt_temp = new TH1F((string("tbkghist_matched_alt")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_unmatched_temp = new TH1F((string("tbkghist_unmatched")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_unmatched_alt_temp = new TH1F((string("tbkghist_unmatched_alt")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      tthist_matched.push_back(dynamic_cast<TH1*>(tthist_matched_temp));
      tthist_matched_alt.push_back(dynamic_cast<TH1*>(tthist_matched_alt_temp));
      tthist_unmatched.push_back(dynamic_cast<TH1*>(tthist_unmatched_temp));
      tthist_unmatched_alt.push_back(dynamic_cast<TH1*>(tthist_unmatched_alt_temp));

    }
    else{
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_alt_temp = new TH1F((string("tbkghist_alt")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }
  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;

  vector<TH2*> tthist_matched_2D;
  vector<TH2*> tthist_matched_alt_2D;
  vector<TH2*> tthist_unmatched_2D;
  vector<TH2*> tthist_unmatched_alt_2D;


  TTree* dttree  = (TTree*)dtfile->Get("tree");
  TTree* vltree  = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree  = (TTree*)ttfile->Get("tree/tree");
  TTree* tttree_alt = NULL;
  if(ttfile_alt)
    tttree_alt = (TTree*) ttfile_alt->Get("tree/tree");

  TTree* dbtree = (TTree*)dbfile->Get("tree/tree");
  TTree* qctree = (TTree*)qcfile->Get("tree/tree");

  vector<TH1*> ehists;

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  if(applyQGLReweight){
    if(not makeResonantSelection){
      cout<<"top+jets control region --> Top"<<endl;
      makehist4(tttree,     tthist,     tthist_2D,      true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL);
      cout<<"top+jets control region --> Top alternative"<<endl;
      makehist4(tttree_alt, tthist_alt, tthist_alt_2D,  true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL);
    }
    else{
      cout<<"top+jets control region --> Top"<<endl;
      makehist4(tttree,     tthist_matched,     tthist_matched_2D,      true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL,1);
      makehist4(tttree,     tthist_unmatched,   tthist_unmatched_2D,    true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL,2);
      cout<<"top+jets control region --> Top alternative"<<endl;
      makehist4(tttree_alt, tthist_matched_alt,   tthist_matched_alt_2D,   true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL,1);
      makehist4(tttree_alt, tthist_unmatched_alt, tthist_unmatched_alt_2D, true,  sample, category, false,  1.00, lumi, 4, ehists, "", true, NULL,2);
    }
    cout<<"top+jets control region --> W+jets"<<endl;
    makehist4(vltree,  vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, 2, ehists, "", true, NULL);
    cout<<"top+jets control region --> Z+jets"<<endl;
    makehist4(vlltree, vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, 1, ehists, "", true, NULL);
  }
  else{
    if(not makeResonantSelection){
      cout<<"top+jets control region --> Top"<<endl;
      makehist4(tttree,     tthist,     tthist_2D,      true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
      cout<<"top+jets control region --> Top alternative"<<endl;
      makehist4(tttree_alt, tthist_alt, tthist_alt_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
    }
    else{
      cout<<"top+jets control region --> Top"<<endl;
      makehist4(tttree,     tthist_matched,     tthist_matched_2D,      true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL,1);
      makehist4(tttree,     tthist_unmatched,   tthist_unmatched_2D,    true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL,2);
      cout<<"top+jets control region --> Top alternative"<<endl;
      makehist4(tttree_alt, tthist_matched_alt,   tthist_matched_alt_2D,   true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL,1);
      makehist4(tttree_alt, tthist_unmatched_alt, tthist_unmatched_alt_2D, true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL,2);
    }
    cout<<"top+jets control region: W+jets"<<endl;
    makehist4(vltree,  vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
    cout<<"top+jets control region: Z+jets"<<endl;
    makehist4(vlltree, vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
  }
  
  cout<<"top+jets control region: Diboson"<<endl;
  makehist4(dbtree,  dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, 0, ehists, "", true, NULL);
  cout<<"top+jets control region: QCD"<<endl;
  makehist4(qctree,  qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
  cout<<"top+jets control region: Data"<<endl;
  makehist4(dttree,  dthist,  dthist_2D,  false, sample, category, false,  1.00, lumi, 0, ehists, "", true, NULL);
  
  // take average of ttbar                                                                                                                                                     
  for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
    makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));

  for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
    makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));

  for(size_t iHisto = 0; iHisto < tthist_matched.size(); iHisto++)
    makeAverage(tthist_matched.at(iHisto),tthist_matched_alt.at(iHisto));

  for(size_t iHisto = 0; iHisto < tthist_matched_2D.size(); iHisto++)
    makeAverage(tthist_matched_2D.at(iHisto),tthist_matched_alt_2D.at(iHisto));

  for(size_t iHisto = 0; iHisto < tthist_unmatched.size(); iHisto++)
    makeAverage(tthist_unmatched.at(iHisto),tthist_unmatched_alt.at(iHisto));

  for(size_t iHisto = 0; iHisto < tthist_unmatched_2D.size(); iHisto++)
    makeAverage(tthist_unmatched_2D.at(iHisto),tthist_unmatched_alt_2D.at(iHisto));

  outfile->cd();
  for(auto hist :  dthist) hist->Write();
  for(auto hist :  tthist) hist->Write();
  for(auto hist :  tthist_matched) hist->Write();
  for(auto hist :  tthist_unmatched) hist->Write();
  for(auto hist :  dbhist) hist->Write();
  for(auto hist :  qchist) hist->Write();
  for(auto hist :  vlhist) hist->Write();
  for(auto hist :  vllhist) hist->Write();


  dtfile->Close();
  vlfile->Close();
  vllfile->Close();
  ttfile->Close();
  if(ttfile_alt) ttfile_alt->Close();
  dbfile->Close();
  qcfile->Close();
  
  cout << "Templates for the top control region computed ..." << endl;
}

