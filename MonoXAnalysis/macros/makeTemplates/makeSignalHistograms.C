#include "makehist.h"
#include "TChain.h"

using namespace std;

void signalHiggshist(TFile* outfile,
		     const Category &category,
		     vector<string>  observables,
		     vector<string>  observables_2D,
		     const double & lumi        = 2.30,
		     const bool   & doShapeSystematics  = false,
		     const string & mH                  = "125",
		     vector<double> xs                  = {4.198E+04,3.925E+03,1.475E+03,9.095E+02},
		     const int & typeOfHiggsSignal      = 0){

  if(xs.size() != 4 or xs.size() != 5)
    cerr<<"signalHiggshist: xs size wrong, should be 4 numbers for each mass points"<<endl;

  cout<<"Start HiggsInvisible: signalHiggshist --> "<<mH<<endl;

  TChain* ggHTree  = new TChain("tree/tree");
  TChain* vbfHTree = new TChain("tree/tree"); 
  TChain* wHTree   = new TChain("tree/tree");
  TChain* zHTree   = new TChain("tree/tree");

  if(typeOfHiggsSignal == 0){
    ggHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_GluGlu_HToInvisible_M"+mH+"*root").c_str());
    vbfHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M"+mH+"*root").c_str());
    wHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_DM_ScalarWH_Mphi-"+mH+"*root").c_str());
    zHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_DM_ScalarZH_Mphi-"+mH+"*.root").c_str());
  }
  else if(typeOfHiggsSignal == 1){ // fermion only
    ggHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_GluGlu_HToInvisible_M"+mH+"*.root").c_str());
    vbfHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_VBF_HToInvisible_M"+mH+"*.root").c_str());
  }
  else if(typeOfHiggsSignal == 2){  // boson only
    wHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_DM_ScalarWH_Mphi-"+mH+"*.root").c_str());
    zHTree->Add((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_DM_ScalarZH_Mphi-"+mH+"*.root").c_str());
  }

  vector<TH1*>  ggHhist;
  vector<TH1*>  vbfHhist;
  vector<TH1*>  wHhist;
  vector<TH1*>  zHhist;  
  vector<TH1*>  ggZHhist;
  vector<TH1*>  ggHhist_renUp, ggHhist_renDw, ggHhist_facUp, ggHhist_facDw;
  vector<TH1*>  ggHhist_bUp, ggHhist_bDw, ggHhist_metJetUp, ggHhist_metJetDw, ggHhist_metResUp, ggHhist_metResDw, ggHhist_metUncUp, ggHhist_metUncDw;
  vector<TH1*>  vbfHhist_bUp, vbfHhist_bDw, vbfHhist_metJetUp, vbfHhist_metJetDw, vbfHhist_metResUp, vbfHhist_metResDw, vbfHhist_metUncUp, vbfHhist_metUncDw;
  vector<TH1*>  wHhist_bUp, wHhist_bDw, wHhist_metJetUp, wHhist_metJetDw, wHhist_metResUp, wHhist_metResDw, wHhist_metUncUp, wHhist_metUncDw;
  vector<TH1*>  zHhist_bUp, zHhist_bDw, zHhist_metJetUp, zHhist_metJetDw, zHhist_metResUp, zHhist_metResDw, zHhist_metUncUp, zHhist_metUncDw;
  vector<TH1*>  ggZHhist_bUp, ggZHhist_bDw, ggZHhist_metJetUp, ggZHhist_metJetDw, ggZHhist_metResUp, ggZHhist_metResDw, ggZHhist_metUncUp, ggZHhist_metUncDw;
  vector<double> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(typeOfHiggsSignal <= 1){
      TH1F* ggHhist_temp  = new TH1F(("ggHhist_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vbfHhist_temp = new TH1F(("vbfHhist_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      ggHhist.push_back(dynamic_cast<TH1*>(ggHhist_temp));
      vbfHhist.push_back(dynamic_cast<TH1*>(vbfHhist_temp));
    }
    if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
      TH1F* wHhist_temp   = new TH1F(("wHhist_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zHhist_temp   = new TH1F(("zHhist_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* ggZHhist_temp = new TH1F(("ggZHhist_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      wHhist.push_back(dynamic_cast<TH1*>(wHhist_temp));
      zHhist.push_back(dynamic_cast<TH1*>(zHhist_temp));
      ggZHhist.push_back(dynamic_cast<TH1*>(ggZHhist_temp));
    }

    if(doShapeSystematics){

      if(typeOfHiggsSignal <= 1){
	TH1F* ggHhist_renUp_temp   = new TH1F(("ggHhist_renUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggHhist_renDw_temp   = new TH1F(("ggHhist_renDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggHhist_facUp_temp   = new TH1F(("ggHhist_facUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggHhist_facDw_temp   = new TH1F(("ggHhist_facDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);      
	ggHhist_renUp.push_back(dynamic_cast<TH1*>(ggHhist_renUp_temp));
	ggHhist_renDw.push_back(dynamic_cast<TH1*>(ggHhist_renDw_temp));
	ggHhist_facUp.push_back(dynamic_cast<TH1*>(ggHhist_facUp_temp));
	ggHhist_facDw.push_back(dynamic_cast<TH1*>(ggHhist_facDw_temp));

	TH1F* ggHhist_bUp_temp  = new TH1F(("ggHhist_bUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_bUp_temp = new TH1F(("vbfHhist_bUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_bUp.push_back(dynamic_cast<TH1*>(ggHhist_bUp_temp));
	vbfHhist_bUp.push_back(dynamic_cast<TH1*>(vbfHhist_bUp_temp));
      
	TH1F* ggHhist_bDw_temp  = new TH1F(("ggHhist_bDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_bDw_temp = new TH1F(("vbfHhist_bDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_bDw.push_back(dynamic_cast<TH1*>(ggHhist_bDw_temp));
	vbfHhist_bDw.push_back(dynamic_cast<TH1*>(vbfHhist_bDw_temp));

	TH1F* ggHhist_metJetUp_temp  = new TH1F(("ggHhist_metJetUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metJetUp_temp = new TH1F(("vbfHhist_metJetUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metJetUp.push_back(dynamic_cast<TH1*>(ggHhist_metJetUp_temp));
	vbfHhist_metJetUp.push_back(dynamic_cast<TH1*>(vbfHhist_metJetUp_temp));

	TH1F* ggHhist_metJetDw_temp  = new TH1F(("ggHhist_metJetDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metJetDw_temp = new TH1F(("vbfHhist_metJetDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metJetDw.push_back(dynamic_cast<TH1*>(ggHhist_metJetDw_temp));
	vbfHhist_metJetDw.push_back(dynamic_cast<TH1*>(vbfHhist_metJetDw_temp));

	TH1F* ggHhist_metResUp_temp  = new TH1F(("ggHhist_metResUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metResUp_temp = new TH1F(("vbfHhist_metResUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metResUp.push_back(dynamic_cast<TH1*>(ggHhist_metResUp_temp));
	vbfHhist_metResUp.push_back(dynamic_cast<TH1*>(vbfHhist_metResUp_temp));

	TH1F* ggHhist_metResDw_temp  = new TH1F(("ggHhist_metResDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metResDw_temp = new TH1F(("vbfHhist_metResDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metResDw.push_back(dynamic_cast<TH1*>(ggHhist_metResDw_temp));
	vbfHhist_metResDw.push_back(dynamic_cast<TH1*>(vbfHhist_metResDw_temp));

	TH1F* ggHhist_metUncUp_temp  = new TH1F(("ggHhist_metUncUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metUncUp_temp = new TH1F(("vbfHhist_metUncUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metUncUp.push_back(dynamic_cast<TH1*>(ggHhist_metUncUp_temp));
	vbfHhist_metUncUp.push_back(dynamic_cast<TH1*>(vbfHhist_metUncUp_temp));

	TH1F* ggHhist_metUncDw_temp  = new TH1F(("ggHhist_metUncDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* vbfHhist_metUncDw_temp = new TH1F(("vbfHhist_metUncDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	ggHhist_metUncDw.push_back(dynamic_cast<TH1*>(ggHhist_metUncDw_temp));
	vbfHhist_metUncDw.push_back(dynamic_cast<TH1*>(vbfHhist_metUncDw_temp));
      }

      if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
	
	TH1F* wHhist_bUp_temp   = new TH1F(("wHhist_bUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_bUp_temp   = new TH1F(("zHhist_bUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_bUp_temp = new TH1F(("ggZHhist_bUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_bUp.push_back(dynamic_cast<TH1*>(wHhist_bUp_temp));
	zHhist_bUp.push_back(dynamic_cast<TH1*>(zHhist_bUp_temp));
	ggZHhist_bUp.push_back(dynamic_cast<TH1*>(ggZHhist_bUp_temp));

	TH1F* wHhist_bDw_temp   = new TH1F(("wHhist_bDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_bDw_temp   = new TH1F(("zHhist_bDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_bDw_temp = new TH1F(("ggZHhist_bDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_bDw.push_back(dynamic_cast<TH1*>(wHhist_bDw_temp));
	zHhist_bDw.push_back(dynamic_cast<TH1*>(zHhist_bDw_temp));
	ggZHhist_bDw.push_back(dynamic_cast<TH1*>(ggZHhist_bDw_temp));
	
	TH1F* wHhist_metJetUp_temp   = new TH1F(("wHhist_metJetUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metJetUp_temp   = new TH1F(("zHhist_metJetUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metJetUp_temp = new TH1F(("ggZHhist_metJetUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	zHhist_metJetUp.push_back(dynamic_cast<TH1*>(zHhist_metJetUp_temp));
	wHhist_metJetUp.push_back(dynamic_cast<TH1*>(wHhist_metJetUp_temp));
	ggZHhist_metJetUp.push_back(dynamic_cast<TH1*>(ggZHhist_metJetUp_temp));
	
	TH1F* wHhist_metJetDw_temp   = new TH1F(("wHhist_metJetDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metJetDw_temp   = new TH1F(("zHhist_metJetDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metJetDw_temp = new TH1F(("ggZHhist_metJetDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_metJetDw.push_back(dynamic_cast<TH1*>(wHhist_metJetDw_temp));
	zHhist_metJetDw.push_back(dynamic_cast<TH1*>(zHhist_metJetDw_temp));
	ggZHhist_metJetDw.push_back(dynamic_cast<TH1*>(ggZHhist_metJetDw_temp));
	
	TH1F* wHhist_metResUp_temp   = new TH1F(("wHhist_metResUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metResUp_temp   = new TH1F(("zHhist_metResUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metResUp_temp = new TH1F(("ggZHhist_metResUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_metResUp.push_back(dynamic_cast<TH1*>(wHhist_metResUp_temp));
	zHhist_metResUp.push_back(dynamic_cast<TH1*>(zHhist_metResUp_temp));
	ggZHhist_metResUp.push_back(dynamic_cast<TH1*>(ggZHhist_metResUp_temp));
	
	TH1F* wHhist_metResDw_temp   = new TH1F(("wHhist_metResDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metResDw_temp   = new TH1F(("zHhist_metResDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metResDw_temp = new TH1F(("ggZHhist_metResDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_metResDw.push_back(dynamic_cast<TH1*>(wHhist_metResDw_temp));
	zHhist_metResDw.push_back(dynamic_cast<TH1*>(zHhist_metResDw_temp));
	ggZHhist_metResDw.push_back(dynamic_cast<TH1*>(ggZHhist_metResDw_temp));
	
	TH1F* wHhist_metUncUp_temp   = new TH1F(("wHhist_metUncUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metUncUp_temp   = new TH1F(("zHhist_metUncUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metUncUp_temp = new TH1F(("ggZHhist_metUncUp_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_metUncUp.push_back(dynamic_cast<TH1*>(wHhist_metUncUp_temp));
	zHhist_metUncUp.push_back(dynamic_cast<TH1*>(zHhist_metUncUp_temp));
	ggZHhist_metUncUp.push_back(dynamic_cast<TH1*>(ggZHhist_metUncUp_temp));
	
	TH1F* wHhist_metUncDw_temp   = new TH1F(("wHhist_metUncDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* zHhist_metUncDw_temp   = new TH1F(("zHhist_metUncDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* ggZHhist_metUncDw_temp = new TH1F(("ggZHhist_metUncDw_"+mH+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	wHhist_metUncDw.push_back(dynamic_cast<TH1*>(wHhist_metUncDw_temp));
	zHhist_metUncDw.push_back(dynamic_cast<TH1*>(zHhist_metUncDw_temp));
	ggZHhist_metUncDw.push_back(dynamic_cast<TH1*>(ggZHhist_metUncDw_temp));
	
      }  
    }
  }

  vector<TH2*> ggHhist_2D; 
  vector<TH2*> vbfHhist_2D;
  vector<TH2*> wHhist_2D;
  vector<TH2*> zHhist_2D;
  vector<TH2*> ggZHhist_2D;

  vector<TH2*> ggHhist_renUp_2D, ggHhist_renDw_2D, ggHhist_facUp_2D, ggHhist_facDw_2D;
  vector<TH2*> ggHhist_bUp_2D, ggHhist_bDw_2D, ggHhist_metJetUp_2D, ggHhist_metJetDw_2D, ggHhist_metResUp_2D, ggHhist_metResDw_2D, ggHhist_metUncUp_2D, ggHhist_metUncDw_2D;
  vector<TH2*> vbfHhist_bUp_2D, vbfHhist_bDw_2D, vbfHhist_metJetUp_2D, vbfHhist_metJetDw_2D, vbfHhist_metResUp_2D, vbfHhist_metResDw_2D, vbfHhist_metUncUp_2D, vbfHhist_metUncDw_2D;
  vector<TH2*> wHhist_bUp_2D, wHhist_bDw_2D, wHhist_metJetUp_2D, wHhist_metJetDw_2D, wHhist_metResUp_2D, wHhist_metResDw_2D, wHhist_metUncUp_2D, wHhist_metUncDw_2D;
  vector<TH2*> zHhist_bUp_2D, zHhist_bDw_2D, zHhist_metJetUp_2D, zHhist_metJetDw_2D, zHhist_metResUp_2D, zHhist_metResDw_2D, zHhist_metUncUp_2D, zHhist_metUncDw_2D;
  vector<TH2*> ggZHhist_bUp_2D, ggZHhist_bDw_2D, ggZHhist_metJetUp_2D, ggZHhist_metJetDw_2D, ggZHhist_metResUp_2D, ggZHhist_metResDw_2D, ggZHhist_metUncUp_2D, ggZHhist_metUncDw_2D;

  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(typeOfHiggsSignal <= 1){
      TH2F* ggHhist_temp  = new TH2F(("ggHhist_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vbfHhist_temp = new TH2F(("vbfHhist_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      ggHhist_2D.push_back(dynamic_cast<TH2*>(ggHhist_temp));
      vbfHhist_2D.push_back(dynamic_cast<TH2*>(vbfHhist_temp));
    }
    if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
      TH2F* wHhist_temp   = new TH2F(("wHhist_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zHhist_temp   = new TH2F(("zHhist_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* ggZHhist_temp = new TH2F(("ggZHhist_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      wHhist_2D.push_back(dynamic_cast<TH2*>(wHhist_temp));
      zHhist_2D.push_back(dynamic_cast<TH2*>(zHhist_temp));
      ggZHhist_2D.push_back(dynamic_cast<TH2*>(ggZHhist_temp));
    }
    
    if(doShapeSystematics){

      if(typeOfHiggsSignal <= 1){

	TH2F* ggHhist_renUp_temp   = new TH2F(("ggHhist_renUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggHhist_renDw_temp   = new TH2F(("ggHhist_renDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggHhist_facUp_temp   = new TH2F(("ggHhist_facUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggHhist_facDw_temp   = new TH2F(("ggHhist_facDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_renUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_renUp_temp));
	ggHhist_renDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_renDw_temp));
	ggHhist_facUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_facUp_temp));
	ggHhist_facDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_facDw_temp));
      
	TH2F* ggHhist_bUp_temp  = new TH2F(("ggHhist_bUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_bUp_temp = new TH2F(("vbfHhist_bUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_bUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_bUp_temp));
	vbfHhist_bUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_bUp_temp));

	TH2F* ggHhist_bDw_temp  = new TH2F(("ggHhist_bDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_bDw_temp = new TH2F(("vbfHhist_bDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_bDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_bDw_temp));
	vbfHhist_bDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_bDw_temp));

	TH2F* ggHhist_metJetUp_temp  = new TH2F(("ggHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metJetUp_temp = new TH2F(("vbfHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metJetUp_temp));
	vbfHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metJetUp_temp));

	TH2F* ggHhist_metJetDw_temp  = new TH2F(("ggHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metJetDw_temp = new TH2F(("vbfHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metJetDw_temp));
	vbfHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metJetDw_temp));

	TH2F* ggHhist_metResUp_temp  = new TH2F(("ggHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metResUp_temp = new TH2F(("vbfHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metResUp_temp));
	vbfHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metResUp_temp));
	
	TH2F* ggHhist_metResDw_temp  = new TH2F(("ggHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metResDw_temp = new TH2F(("vbfHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metResDw_temp));
	vbfHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metResDw_temp));
	
	TH2F* ggHhist_metUncUp_temp  = new TH2F(("ggHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metUncUp_temp = new TH2F(("vbfHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metUncUp_temp));
	vbfHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metUncUp_temp));

	TH2F* ggHhist_metUncDw_temp  = new TH2F(("ggHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* vbfHhist_metUncDw_temp = new TH2F(("vbfHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	ggHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metUncDw_temp));
	vbfHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metUncDw_temp));
      }

      if(typeOfHiggsSignal == 0 or typeOfHiggsSignal>= 2){

	TH2F* wHhist_bUp_temp   = new TH2F(("wHhist_bUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_bUp_temp   = new TH2F(("zHhist_bUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_bUp_temp = new TH2F(("ggZHhist_bUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_bUp_2D.push_back(dynamic_cast<TH2*>(wHhist_bUp_temp));
	zHhist_bUp_2D.push_back(dynamic_cast<TH2*>(zHhist_bUp_temp));
	ggZHhist_bUp_2D.push_back(dynamic_cast<TH2*>(ggZHhist_bUp_temp));
	
	TH2F* wHhist_bDw_temp   = new TH2F(("wHhist_bDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_bDw_temp   = new TH2F(("zHhist_bDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_bDw_temp = new TH2F(("ggZHhist_bDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_bDw_2D.push_back(dynamic_cast<TH2*>(wHhist_bDw_temp));
	zHhist_bDw_2D.push_back(dynamic_cast<TH2*>(zHhist_bDw_temp));
	ggZHhist_bDw_2D.push_back(dynamic_cast<TH2*>(ggZHhist_bDw_temp));

	TH2F* wHhist_metJetUp_temp   = new TH2F(("wHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metJetUp_temp   = new TH2F(("zHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metJetUp_temp = new TH2F(("ggZHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	zHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metJetUp_temp));
	wHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metJetUp_temp));
	ggZHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metJetUp_temp));
	
	TH2F* wHhist_metJetDw_temp   = new TH2F(("wHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metJetDw_temp   = new TH2F(("zHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metJetDw_temp = new TH2F(("ggZHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metJetDw_temp));
	zHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metJetDw_temp));
	ggZHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metJetDw_temp));
	
	TH2F* wHhist_metResUp_temp   = new TH2F(("wHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metResUp_temp   = new TH2F(("zHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metResUp_temp = new TH2F(("ggZHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metResUp_temp));
	zHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metResUp_temp));
	ggZHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metResUp_temp));

	TH2F* wHhist_metResDw_temp   = new TH2F(("wHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metResDw_temp   = new TH2F(("zHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metResDw_temp = new TH2F(("ggZHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metResDw_temp));
	zHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metResDw_temp));
	ggZHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metResDw_temp));
	
	TH2F* wHhist_metUncUp_temp   = new TH2F(("wHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metUncUp_temp   = new TH2F(("zHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metUncUp_temp = new TH2F(("ggZHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metUncUp_temp));
	zHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metUncUp_temp));
	ggZHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metUncUp_temp));
	
	TH2F* wHhist_metUncDw_temp   = new TH2F(("wHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* zHhist_metUncDw_temp   = new TH2F(("zHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* ggZHhist_metUncDw_temp = new TH2F(("ggZHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	wHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metUncDw_temp));
	zHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metUncDw_temp));
	ggZHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(ggZHhist_metUncDw_temp));
      }
    }  
  }
  
  // start running on signal samples                                                                                                                                 
  vector<TH1*> ehists;
  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;

  //Take histograms for Higgs
  TFile* higgsPT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/HiggsPT/hpt_ggH_125_13TeV_POWHEG_weights.root");

  TH1* renUp = (TH1*) higgsPT->FindObjectAny("weight_pT_rup");
  TH1* renDw = (TH1*) higgsPT->FindObjectAny("weight_pT_rdn");
  TH1* facUp = (TH1*) higgsPT->FindObjectAny("weight_pT_fup");
  TH1* facDw = (TH1*) higgsPT->FindObjectAny("weight_pT_fdn");
  for(int iBinX = 1; iBinX <= facDw->GetNbinsX(); iBinX++)
    facDw->SetBinContent(iBinX,1./facUp->GetBinContent(iBinX));


  // ggZH re-weight
  TFile* ggZHPT  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/HiggsPT/ggZH_weight.root");
  TH2* ggZH_weight = (TH2*) ggZHPT->FindObjectAny("Graph2D_from_weight");

  

  // higgs re-weight
  if(typeOfHiggsSignal <=1){
    makehist4(ggHTree,ggHhist,ggHhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(0));
    makehist4(vbfHTree,vbfHhist,vbfHhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(1));
  }
  if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
    makehist4(wHTree,wHhist,wHhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(2));
    makehist4(zHTree,zHhist,zHhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(3));
    if(xs.size() > 4) // apply ZH re-weight to estimated ggZH                                                                                                               
      makehist4(zHTree,ggZHhist,ggZHhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);

  }

  if(doShapeSystematics){
 
    if(typeOfHiggsSignal <= 1){
      makehist4(ggHTree,ggHhist_renUp,ggHhist_renUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(0),renUp);
      makehist4(ggHTree,ggHhist_renDw,ggHhist_renDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(0),renDw);
      makehist4(ggHTree,ggHhist_facUp,ggHhist_facUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(0),facUp);
      makehist4(ggHTree,ggHhist_facDw,ggHhist_facDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,true,false,xs.at(0),facDw);
    
      makehist4(ggHTree,ggHhist_bUp,ggHhist_bUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_bUp,vbfHhist_bUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_bDw,ggHhist_bDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_bDw,vbfHhist_bDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metJetUp,ggHhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metJetUp,vbfHhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metJetDw,ggHhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metJetDw,vbfHhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metResUp,ggHhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metResUp,vbfHhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metResDw,ggHhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metResDw,vbfHhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metUncUp,ggHhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metUncUp,vbfHhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,true,false,xs.at(1));
      makehist4(ggHTree,ggHhist_metUncDw,ggHhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,true,false,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metUncDw,vbfHhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,true,false,xs.at(1));
    }
    if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >= 2){

      makehist4(wHTree,wHhist_bUp,wHhist_bUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_bUp,zHhist_bUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,reweightNVTX,0,true,false,xs.at(3));      
      makehist4(wHTree,wHhist_bDw,wHhist_bDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_bDw,zHhist_bDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,reweightNVTX,0,true,false,xs.at(3));   
      makehist4(wHTree,wHhist_metJetUp,wHhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metJetUp,zHhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,true,false,xs.at(3));    
      makehist4(wHTree,wHhist_metJetDw,wHhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metJetDw,zHhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,true,false,xs.at(3));    
      makehist4(wHTree,wHhist_metResUp,wHhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metResUp,zHhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,true,false,xs.at(3));
      makehist4(wHTree,wHhist_metResDw,wHhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metResDw,zHhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,true,false,xs.at(3));
      makehist4(wHTree,wHhist_metUncUp,wHhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metUncUp,zHhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,true,false,xs.at(3));
      makehist4(wHTree,wHhist_metUncDw,wHhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,true,false,xs.at(2));
      makehist4(zHTree,zHhist_metUncDw,zHhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,true,false,xs.at(3));

      if(xs.size() >4){ // re-weight qqZH to get ggZH contribution                                                                                                           
        makehist4(zHTree,ggZHhist_bUp,ggZHhist_bUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_bDw,ggZHhist_bDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metJetUp,ggZHhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metJetDw,ggZHhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metResUp,ggZHhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metResDw,ggZHhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metUncUp,ggZHhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);
        makehist4(zHTree,ggZHhist_metUncDw,ggZHhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,true,0,true,false,xs.at(4),NULL,ggZH_weight);

      }
    }
  }
  

  //smooth
  for(auto hist: ggHhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: ggHhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: vbfHhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: vbfHhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist: wHhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: wHhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist: zHhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: zHhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist: ggZHhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist: ggZHhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}


  if(doShapeSystematics){

    for(auto hist: ggHhist_renUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_renUp_2D){if(TString(hist->GetName()).Contains("_met"))smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_renDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_renDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_facUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_facUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_facDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_facDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_bDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggHhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggHhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    ////
    for(auto hist: vbfHhist_bUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_bUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_bDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_bDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metJetUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metJetUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metJetDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metJetDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metResUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metResUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metResDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metResDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metUncUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metUncUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: vbfHhist_metUncDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: vbfHhist_metUncDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    ////
    for(auto hist: wHhist_bUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_bUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_bDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_bDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metJetUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metJetUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metJetDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metJetDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metResUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metResUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metResDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metResDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metUncUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metUncUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: wHhist_metUncDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: wHhist_metUncDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    ////
    for(auto hist: zHhist_bUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_bUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_bDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_bDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metJetUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metJetUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metJetDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metJetDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metResUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metResUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metResDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metResDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metUncUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metUncUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: zHhist_metUncDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: zHhist_metUncDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    //                                                                                                                                                                           
    for(auto hist: ggZHhist_bUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_bUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_bDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_bDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metJetUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metJetUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metJetDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metJetDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metResUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metResUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metResDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metResDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metUncUp){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metUncUp_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist: ggZHhist_metUncDw){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist: ggZHhist_metUncDw_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  }

  // store tempaltes in the output file                                                                                                                                     
  outfile->cd();
  if(not outfile->GetDirectory("ggH"))
    outfile->mkdir("ggH");
  outfile->cd("ggH");
  for(auto hist : ggHhist)  hist->Write();
  for(auto hist : ggHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

  outfile->cd();
  if(not outfile->GetDirectory("vbfH"))
    outfile->mkdir("vbfH");
  outfile->cd("vbfH");  
  for(auto hist : vbfHhist) hist->Write();
  for(auto hist : vbfHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

  outfile->cd();
  if(not outfile->GetDirectory("wH"))
    outfile->mkdir("wH");
  outfile->cd("wH");  
  for(auto hist : wHhist)   hist->Write();
  for(auto hist : wHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

  outfile->cd();
  if(not outfile->GetDirectory("zH"))
    outfile->mkdir("zH");
  outfile->cd("zH");  
  for(auto hist : zHhist)   hist->Write();
  for(auto hist : zHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

  outfile->cd();
  if(not outfile->GetDirectory("ggZH"))
    outfile->mkdir("ggZH");
  outfile->cd("ggZH");
  for(auto hist : ggZHhist)   hist->Write();
  for(auto hist : ggZHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}



  if(doShapeSystematics){

    outfile->cd();
    if(not outfile->GetDirectory("ggH/sysShape"))
      outfile->mkdir("ggH/sysShape");
    outfile->cd("ggH/sysShape");

    for(auto hist : ggHhist_renUp) hist->Write();
    for(auto hist : ggHhist_renDw) hist->Write();
    for(auto hist : ggHhist_facUp) hist->Write();
    for(auto hist : ggHhist_facDw) hist->Write();

    for(auto hist : ggHhist_bUp) hist->Write();
    for(auto hist : ggHhist_bDw) hist->Write();
    for(auto hist : ggHhist_metJetUp) hist->Write();
    for(auto hist : ggHhist_metJetDw) hist->Write();
    for(auto hist : ggHhist_metResUp) hist->Write();
    for(auto hist : ggHhist_metResDw) hist->Write();
    for(auto hist : ggHhist_metUncUp) hist->Write();
    for(auto hist : ggHhist_metUncDw) hist->Write();

    for(auto hist : ggHhist_renUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_renDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_facUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_facDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist : ggHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("vbfH/sysShape"))
      outfile->mkdir("vbfH/sysShape");
    outfile->cd("vbfH/sysShape");

    for(auto hist : vbfHhist_bUp) hist->Write();
    for(auto hist : vbfHhist_bDw) hist->Write();
    for(auto hist : vbfHhist_metJetUp) hist->Write();
    for(auto hist : vbfHhist_metJetDw) hist->Write();
    for(auto hist : vbfHhist_metResUp) hist->Write();
    for(auto hist : vbfHhist_metResDw) hist->Write();
    for(auto hist : vbfHhist_metUncUp) hist->Write();
    for(auto hist : vbfHhist_metUncDw) hist->Write();

    for(auto hist :vbfHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist :vbfHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("wH/sysShape"))
      outfile->mkdir("wH/sysShape");
    outfile->cd("wH/sysShape");

    for(auto hist : wHhist_bUp) hist->Write();
    for(auto hist : wHhist_bDw) hist->Write();
    for(auto hist : wHhist_metJetUp) hist->Write();
    for(auto hist : wHhist_metJetDw) hist->Write();
    for(auto hist : wHhist_metResUp) hist->Write();
    for(auto hist : wHhist_metResDw) hist->Write();
    for(auto hist : wHhist_metUncUp) hist->Write();
    for(auto hist : wHhist_metUncDw) hist->Write();

    for(auto hist : wHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metUncUp_2D){  TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("zH/sysShape"))
      outfile->mkdir("zH/sysShape");
    outfile->cd("zH/sysShape");

    for(auto hist : zHhist_bUp) hist->Write();
    for(auto hist : zHhist_bDw) hist->Write();
    for(auto hist : zHhist_metJetUp) hist->Write();
    for(auto hist : zHhist_metJetDw) hist->Write();
    for(auto hist : zHhist_metResUp) hist->Write();
    for(auto hist : zHhist_metResDw) hist->Write();
    for(auto hist : zHhist_metUncUp) hist->Write();
    for(auto hist : zHhist_metUncDw) hist->Write();


    for(auto hist : zHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("ggZH/sysShape"))
      outfile->mkdir("ggZH/sysShape");
    outfile->cd("ggZH/sysShape");

    for(auto hist : ggZHhist_bUp) hist->Write();
    for(auto hist : ggZHhist_bDw) hist->Write();
    for(auto hist : ggZHhist_metJetUp) hist->Write();
    for(auto hist : ggZHhist_metJetDw) hist->Write();
    for(auto hist : ggZHhist_metResUp) hist->Write();
    for(auto hist : ggZHhist_metResDw) hist->Write();
    for(auto hist : ggZHhist_metUncUp) hist->Write();
    for(auto hist : ggZHhist_metUncDw) hist->Write();
  
    for(auto hist : ggZHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggZHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}
  

  }

  outfile->cd();

  ggHhist.clear();
  vbfHhist.clear();
  wHhist.clear();
  zHhist.clear();
  ggZHhist.clear();
  ggHhist_renUp.clear(); 
  ggHhist_renDw.clear(); 
  ggHhist_facUp.clear(); 
  ggHhist_facDw.clear();
  ggHhist_bUp.clear(); 
  ggHhist_bDw.clear(); 
  ggHhist_metJetUp.clear(); 
  ggHhist_metJetDw.clear(); 
  ggHhist_metResUp.clear(); 
  ggHhist_metResDw.clear(); 
  ggHhist_metUncUp.clear(); 
  ggHhist_metUncDw.clear();
  vbfHhist_bUp.clear(); 
  vbfHhist_bDw.clear(); 
  vbfHhist_metJetUp.clear(); 
  vbfHhist_metJetDw.clear(); 
  vbfHhist_metResUp.clear(); 
  vbfHhist_metResDw.clear(); 
  vbfHhist_metUncUp.clear(); 
  vbfHhist_metUncDw.clear();
  wHhist_bUp.clear(); 
  wHhist_bDw.clear(); 
  wHhist_metJetUp.clear(); 
  wHhist_metJetDw.clear(); 
  wHhist_metResUp.clear(); 
  wHhist_metResDw.clear(); 
  wHhist_metUncUp.clear(); 
  wHhist_metUncDw.clear();
  zHhist_bUp.clear(); 
  zHhist_bDw.clear(); 
  zHhist_metJetUp.clear(); 
  zHhist_metJetDw.clear(); 
  zHhist_metResUp.clear(); 
  zHhist_metResDw.clear(); 
  zHhist_metUncUp.clear(); 
  zHhist_metUncDw.clear();
  ggZHhist_bUp.clear(); 
  ggZHhist_bDw.clear(); 
  ggZHhist_metJetUp.clear(); 
  ggZHhist_metJetDw.clear(); 
  ggZHhist_metResUp.clear(); 
  ggZHhist_metResDw.clear(); 
  ggZHhist_metUncUp.clear(); 
  ggZHhist_metUncDw.clear();

  ggHhist_2D.clear();
  vbfHhist_2D.clear();
  wHhist_2D.clear();
  zHhist_2D.clear();
  ggHhist_renUp_2D.clear(); 
  ggHhist_renDw_2D.clear(); 
  ggHhist_facUp_2D.clear(); 
  ggHhist_facDw_2D.clear();
  ggHhist_bUp_2D.clear(); 
  ggHhist_bDw_2D.clear(); 
  ggHhist_metJetUp_2D.clear(); 
  ggHhist_metJetDw_2D.clear(); 
  ggHhist_metResUp_2D.clear(); 
  ggHhist_metResDw_2D.clear(); 
  ggHhist_metUncUp_2D.clear(); 
  ggHhist_metUncDw_2D.clear();
  vbfHhist_bUp_2D.clear(); 
  vbfHhist_bDw_2D.clear(); 
  vbfHhist_metJetUp_2D.clear(); 
  vbfHhist_metJetDw_2D.clear(); 
  vbfHhist_metResUp_2D.clear(); 
  vbfHhist_metResDw_2D.clear(); 
  vbfHhist_metUncUp_2D.clear(); 
  vbfHhist_metUncDw_2D.clear();
  wHhist_bUp_2D.clear(); 
  wHhist_bDw_2D.clear(); 
  wHhist_metJetUp_2D.clear(); 
  wHhist_metJetDw_2D.clear(); 
  wHhist_metResUp_2D.clear(); 
  wHhist_metResDw_2D.clear(); 
  wHhist_metUncUp_2D.clear(); 
  wHhist_metUncDw_2D.clear();
  zHhist_bUp_2D.clear(); 
  zHhist_bDw_2D.clear(); 
  zHhist_metJetUp_2D.clear(); 
  zHhist_metJetDw_2D.clear(); 
  zHhist_metResUp_2D.clear(); 
  zHhist_metResDw_2D.clear(); 
  zHhist_metUncUp_2D.clear(); 
  zHhist_metUncDw_2D.clear();
  ggZHhist_bUp_2D.clear(); 
  ggZHhist_bDw_2D.clear(); 
  ggZHhist_metJetUp_2D.clear(); 
  ggZHhist_metJetDw_2D.clear(); 
  ggZHhist_metResUp_2D.clear(); 
  ggZHhist_metResDw_2D.clear(); 
  ggZHhist_metUncUp_2D.clear(); 
  ggZHhist_metUncDw_2D.clear();

  if(higgsPT) higgsPT->Close();
  if(ggZHPT) ggZHPT->Close();

  cout << "Templates for the Higgs invisible ..." << endl;
}


// Build templates for the signal region                                                                                                                                       
void signalmchist(TFile* outfile,
		  const Category & category,
		  vector<string> observables,
		  vector<string> observables_2D,
		  vector<signalSample> massPoint,
		  const string & interaction       = "Vector",
		  const double & lumi              = 2.24,
		  const bool & doShapeSystematics  = false){

  
  vector<TFile*> monoJfile;
  vector<TFile*> monoWfile;
  vector<TFile*> monoZfile;

  for(auto iPoint : massPoint){
    if(iPoint.interaction == interaction and interaction == "Vector"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMV_Vector/sigfilter/sig_tree_DMV_NNPDF30_Vector_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Vector/sigfilter/sig_tree_VectorMonoZ_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root").c_str()));
    }
    else if(iPoint.interaction == interaction and interaction == "Axial"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMV_Axial/sigfilter/sig_tree_DMV_NNPDF30_Axial_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Axial/sigfilter/sig_tree_AxialMonoZ_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root").c_str()));
    }
    else if(iPoint.interaction == interaction and interaction == "Scalar"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMS_Scalar/sigfilter/sig_tree_DMS_NNPDF30_Scalar_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_v2_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root").c_str()));      
    }
    else if(iPoint.interaction == interaction and interaction == "Pseudoscalar"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMS_Pseudoscalar/sigfilter/sig_tree_DMS_NNPDF30_Pseudoscalar_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_v2_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Pseudoscalar/sigfilter/sig_tree_DMS_PseudoscalarWH_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Pseudoscalar/sigfilter/sig_tree_DMS_PseudoscalarZH_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root").c_str()));
    }
  }

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

  // create 1D histgrams for each process and shape uncertainty
  int imass = 0;
  vector<double> bins;
  for(auto iPoint : massPoint){
    if(iPoint.interaction != interaction) continue;
    
    for(auto obs : observables){
      
      bins = selectBinning(obs,category);
      if(bins.empty())
        cout<<"No binning for this observable --> please define it"<<endl;
 
      TH1F* monoJhist_temp = new TH1F(("monoJhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* monoWhist_temp = new TH1F(("monoWhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* monoZhist_temp = new TH1F(("monoZhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      monoJhist.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_temp));
      monoWhist.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_temp));
      monoZhist.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_temp));

      if(doShapeSystematics){

	TH1F* monoJhist_bUp_temp = new TH1F(("monoJhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_bUp_temp = new TH1F(("monoWhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_bUp_temp = new TH1F(("monoZhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bUp_temp));
	monoWhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bUp_temp));
	monoZhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bUp_temp));

	TH1F* monoJhist_bDw_temp = new TH1F(("monoJhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_bDw_temp = new TH1F(("monoWhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_bDw_temp = new TH1F(("monoZhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bDw_temp));
	monoWhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bDw_temp));
	monoZhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bDw_temp));


	TH1F* monoJhist_metJetUp_temp = new TH1F(("monoJhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metJetUp_temp = new TH1F(("monoWhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metJetUp_temp = new TH1F(("monoZhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetUp_temp));
	monoWhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetUp_temp));
	monoZhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetUp_temp));

	TH1F* monoJhist_metJetDw_temp = new TH1F(("monoJhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metJetDw_temp = new TH1F(("monoWhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metJetDw_temp = new TH1F(("monoZhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetDw_temp));
	monoWhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetDw_temp));
	monoZhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetDw_temp));

	TH1F* monoJhist_metResUp_temp = new TH1F(("monoJhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metResUp_temp = new TH1F(("monoWhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metResUp_temp = new TH1F(("monoZhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResUp_temp));
	monoWhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResUp_temp));
	monoZhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResUp_temp));

	TH1F* monoJhist_metResDw_temp = new TH1F(("monoJhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metResDw_temp = new TH1F(("monoWhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metResDw_temp = new TH1F(("monoZhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResDw_temp));
	monoWhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResDw_temp));
	monoZhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResDw_temp));

	TH1F* monoJhist_metUncUp_temp = new TH1F(("monoJhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metUncUp_temp = new TH1F(("monoWhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metUncUp_temp = new TH1F(("monoZhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metUncUp_temp));
	monoWhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metUncUp_temp));
	monoZhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metUncUp_temp));

	TH1F* monoJhist_metUncDw_temp = new TH1F(("monoJhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoWhist_metUncDw_temp = new TH1F(("monoWhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* monoZhist_metUncDw_temp = new TH1F(("monoZhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	monoJhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metUncDw_temp));
	monoWhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metUncDw_temp));
	monoZhist_metUncDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metUncDw_temp));
	
      }      
    }
    imass++;
  }

  // allocate 2D histo
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

  imass = 0;

  for(auto iPoint : massPoint){
    if(iPoint.interaction != interaction) continue;
    for(auto obs : observables_2D){

      bin2D bins = selectBinning2D(obs,category);      
     if(bins.binX.empty() or bins.binY.empty())
        cout<<"No binning for this observable --> please define it"<<endl;

      TH2F* monoJhist_temp = new TH2F(("monoJhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* monoWhist_temp = new TH2F(("monoWhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* monoZhist_temp = new TH2F(("monoZhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      monoJhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_temp));
      monoWhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_temp));
      monoZhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_temp));

      if(doShapeSystematics){

	TH2F* monoJhist_bUp_temp = new TH2F(("monoJhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_bUp_temp = new TH2F(("monoWhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_bUp_temp = new TH2F(("monoZhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_bUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_bUp_temp));
	monoWhist_bUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_bUp_temp));
	monoZhist_bUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_bUp_temp));

	TH2F* monoJhist_bDw_temp = new TH2F(("monoJhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_bDw_temp = new TH2F(("monoWhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_bDw_temp = new TH2F(("monoZhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_bDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_bDw_temp));
	monoWhist_bDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_bDw_temp));
	monoZhist_bDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_bDw_temp));


	TH2F* monoJhist_metJetUp_temp = new TH2F(("monoJhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metJetUp_temp = new TH2F(("monoWhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metJetUp_temp = new TH2F(("monoZhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metJetUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metJetUp_temp));
	monoWhist_metJetUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metJetUp_temp));
	monoZhist_metJetUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metJetUp_temp));

	TH2F* monoJhist_metJetDw_temp = new TH2F(("monoJhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metJetDw_temp = new TH2F(("monoWhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metJetDw_temp = new TH2F(("monoZhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metJetDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metJetDw_temp));
	monoWhist_metJetDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metJetDw_temp));
	monoZhist_metJetDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metJetDw_temp));

	TH2F* monoJhist_metResUp_temp = new TH2F(("monoJhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metResUp_temp = new TH2F(("monoWhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metResUp_temp = new TH2F(("monoZhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metResUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metResUp_temp));
	monoWhist_metResUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metResUp_temp));
	monoZhist_metResUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metResUp_temp));

	TH2F* monoJhist_metResDw_temp = new TH2F(("monoJhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metResDw_temp = new TH2F(("monoWhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metResDw_temp = new TH2F(("monoZhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metResDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metResDw_temp));
	monoWhist_metResDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metResDw_temp));
	monoZhist_metResDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metResDw_temp));

	TH2F* monoJhist_metUncUp_temp = new TH2F(("monoJhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metUncUp_temp = new TH2F(("monoWhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metUncUp_temp = new TH2F(("monoZhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metUncUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metUncUp_temp));
	monoWhist_metUncUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metUncUp_temp));
	monoZhist_metUncUp_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metUncUp_temp));

	TH2F* monoJhist_metUncDw_temp = new TH2F(("monoJhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoWhist_metUncDw_temp = new TH2F(("monoWhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* monoZhist_metUncDw_temp = new TH2F(("monoZhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	monoJhist_metUncDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metUncDw_temp));
	monoWhist_metUncDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metUncDw_temp));
	monoZhist_metUncDw_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metUncDw_temp));

      }
    }
    imass++;
  }


  // getting trees
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


  // start running on signal samples
  vector<TH1*> ehists;
  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;

  int itree = 0;
  for(auto tree : monoJtree){
    // signals                                                                                                                                                                 
    if(tree){
      cout<<"signal region analysis --> Signal monoJet "<<endl;
      makehist4(tree,monoJhist.at(itree),monoJhist_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoJ sys "<<endl;
	makehist4(tree,monoJhist_bUp.at(itree),monoJhist_bUp_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"btagUp",false,reweightNVTX);
	makehist4(tree,monoJhist_bDw.at(itree),monoJhist_bDw_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"btagDown",false,reweightNVTX);
	makehist4(tree,monoJhist_metJetUp.at(itree),monoJhist_metJetUp_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX);
	makehist4(tree,monoJhist_metJetDw.at(itree),monoJhist_metJetDw_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX);
	makehist4(tree,monoJhist_metResUp.at(itree),monoJhist_metResUp_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX);
	makehist4(tree,monoJhist_metResDw.at(itree),monoJhist_metResDw_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX);
	makehist4(tree,monoJhist_metUncUp.at(itree),monoJhist_metUncUp_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX);
	makehist4(tree,monoJhist_metUncDw.at(itree),monoJhist_metUncDw_2D.at(itree),true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoWtree){
    if(tree){
      cout<<"signal region analysis --> Signal monoW "<<endl;
      makehist4(tree,monoWhist.at(itree),monoWhist_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoW sys "<<endl;
	makehist4(tree,monoWhist_bUp.at(itree),monoWhist_bUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"btagUp",false,reweightNVTX);
	makehist4(tree,monoWhist_bDw.at(itree),monoWhist_bDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"btagDown",false,reweightNVTX);
	makehist4(tree,monoWhist_metJetUp.at(itree),monoWhist_metJetUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesUp",false,reweightNVTX);
	makehist4(tree,monoWhist_metJetDw.at(itree),monoWhist_metJetDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesDw",false,reweightNVTX);
	makehist4(tree,monoWhist_metResUp.at(itree),monoWhist_metResUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerUp",false,reweightNVTX);
	makehist4(tree,monoWhist_metResDw.at(itree),monoWhist_metResDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerDw",false,reweightNVTX);
	makehist4(tree,monoWhist_metUncUp.at(itree),monoWhist_metUncUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncUp",false,reweightNVTX);
	makehist4(tree,monoWhist_metUncDw.at(itree),monoWhist_metUncDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncDw",false,reweightNVTX);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoZtree){
    if(tree){
      cout<<"signal region analysis --> Signal monoZ "<<endl;
      makehist4(tree,monoZhist.at(itree),monoZhist_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoZ sys "<<endl;
	makehist4(tree,monoZhist_bUp.at(itree),monoZhist_bUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"btagUp",false,reweightNVTX);
	makehist4(tree,monoZhist_bDw.at(itree),monoZhist_bDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"btagDown",false,reweightNVTX);
	makehist4(tree,monoZhist_metJetUp.at(itree),monoZhist_metJetUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesUp",false,reweightNVTX);
	makehist4(tree,monoZhist_metJetDw.at(itree),monoZhist_metJetDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesDw",false,reweightNVTX);
	makehist4(tree,monoZhist_metResUp.at(itree),monoZhist_metResUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerUp",false,reweightNVTX);
	makehist4(tree,monoZhist_metResDw.at(itree),monoZhist_metResDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerDw",false,reweightNVTX);
	makehist4(tree,monoZhist_metUncUp.at(itree),monoZhist_metUncUp_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncUp",false,reweightNVTX);
	makehist4(tree,monoZhist_metUncDw.at(itree),monoZhist_metUncDw_2D.at(itree),true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncDw",false,reweightNVTX);
      }
    }
    itree++;
  }
  

  //smooth                                                                                                                                                                      
  for(auto sample : monoJhist){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  for(auto sample : monoWhist){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  for(auto sample : monoZhist){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  for(auto sample : monoJhist_2D){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  for(auto sample : monoWhist_2D){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  for(auto sample : monoZhist_2D){ for (auto histo : sample){ if(TString(histo->GetName()).Contains("_met")) smoothEmptyBins(histo,2);}}
  

  if(doShapeSystematics){

    for(auto sample : monoJhist_bUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_bUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_bUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_bUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_bUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_bUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
   
    for(auto sample : monoJhist_bDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_bDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_bDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_bDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_bDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_bDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    
    for(auto sample : monoJhist_metJetUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metJetUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metJetUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metJetUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metJetUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metJetUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}

    for(auto sample : monoJhist_metJetDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metJetDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metJetDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metJetDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metJetDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metJetDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}

    for(auto sample : monoJhist_metResUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metResUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metResUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metResUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metResUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metResUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}

    for(auto sample : monoJhist_metResDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metResDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metResDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metResDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metResDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metResDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}

    for(auto sample : monoJhist_metUncUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metUncUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metUncUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metUncUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metUncUp){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metUncUp_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    
    for(auto sample : monoJhist_metUncDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoJhist_metUncDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metUncDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoWhist_metUncDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metUncDw){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
    for(auto sample : monoZhist_metUncDw_2D){for (auto histo : sample){if(TString(histo->GetName()).Contains("_met"))smoothEmptyBins(histo,2);}}
   
  }
  
  // store tempaltes in the output file
  outfile->cd();
  if(not outfile->GetDirectory("monoJ"))
    outfile->mkdir("monoJ");
  outfile->cd("monoJ");
  for(auto sample : monoJhist){ for(auto hist : sample) hist->Write();}
  outfile->cd();
  if(not outfile->GetDirectory("monoW"))
    outfile->mkdir("monoW");
  outfile->cd("monoW");
  for(auto sample : monoWhist){ for(auto hist : sample) hist->Write();}  
  outfile->cd();
  if(not outfile->GetDirectory("monoZ"))
    outfile->mkdir("monoZ");
  outfile->cd("monoZ");
  for(auto sample : monoZhist){ for(auto hist : sample) hist->Write();}
  outfile->cd();

  if(doShapeSystematics){
    if(not outfile->GetDirectory("monoJ/sysShape"))
      outfile->mkdir("monoJ/sysShape");
    outfile->cd("monoJ/sysShape");
    for(auto sample : monoJhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncDw){ for(auto hist : sample) hist->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("monoW/sysShape"))
      outfile->mkdir("monoW/sysShape");
    outfile->cd("monoW/sysShape");
    for(auto sample : monoWhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncDw){ for(auto hist : sample) hist->Write();}

    outfile->cd();
    if(not outfile->GetDirectory("monoZ/sysShape"))
      outfile->mkdir("monoZ/sysShape");
    outfile->cd("monoZ/sysShape");
    for(auto sample : monoZhist_bUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_bDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResDw){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncUp){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncDw){ for(auto hist : sample) hist->Write();}
  }

  outfile->cd();
  outfile->cd("monoJ");
  for(auto sample : monoJhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
  outfile->cd();
  outfile->cd("monoW");
  for(auto sample : monoWhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
  outfile->cd();
  outfile->cd("monoZ");
  for(auto sample : monoZhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}

  if(doShapeSystematics){

    outfile->cd();
    outfile->cd("monoJ/sysShape/");
    for(auto sample : monoJhist_bUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_bDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metJetUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metJetDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metResUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metResDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metUncUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metUncDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}


    outfile->cd();
    outfile->cd("monoW/sysShape/");
    for(auto sample : monoWhist_bUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_bDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metJetUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metJetDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metResUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metResDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metUncUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metUncDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}

    outfile->cd();
    outfile->cd("monoZ/sysShape/");
    for(auto sample : monoZhist_bUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_bDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metJetUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metJetDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metResUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metResDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metUncUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoZhist_metUncDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}    
  }

  for(auto file : monoJfile){ if(file) file->Close();}
  for(auto file : monoWfile){ if(file) file->Close();}
  for(auto file : monoZfile){ if(file) file->Close();}
  
  monoJhist.clear();
  monoWhist.clear();
  monoZhist.clear();

  monoJhist_bUp.clear();
  monoJhist_bUp.clear();
  monoWhist_bUp.clear();
  monoZhist_bUp.clear();
  monoJhist_bDw.clear();
  monoWhist_bDw.clear();
  monoZhist_bDw.clear();
  monoJhist_metJetUp.clear();
  monoWhist_metJetUp.clear();
  monoZhist_metJetUp.clear();
  monoJhist_metJetDw.clear();
  monoWhist_metJetDw.clear();
  monoZhist_metJetDw.clear();
  monoJhist_metResUp.clear();
  monoWhist_metResUp.clear();
  monoZhist_metResUp.clear();
  monoJhist_metResDw.clear();
  monoWhist_metResDw.clear();
  monoZhist_metResDw.clear();
  monoJhist_metUncUp.clear();
  monoWhist_metUncUp.clear();
  monoZhist_metUncUp.clear();
  monoJhist_metUncDw.clear();
  monoWhist_metUncDw.clear();
  monoZhist_metUncDw.clear();


  monoJhist_2D.clear();
  monoWhist_2D.clear();
  monoZhist_2D.clear();

  monoJhist_bUp_2D.clear();
  monoJhist_bUp_2D.clear();
  monoWhist_bUp_2D.clear();
  monoZhist_bUp_2D.clear();
  monoJhist_bDw_2D.clear();
  monoWhist_bDw_2D.clear();
  monoZhist_bDw_2D.clear();
  monoJhist_metJetUp_2D.clear();
  monoWhist_metJetUp_2D.clear();
  monoZhist_metJetUp_2D.clear();
  monoJhist_metJetDw_2D.clear();
  monoWhist_metJetDw_2D.clear();
  monoZhist_metJetDw_2D.clear();
  monoJhist_metResUp_2D.clear();
  monoWhist_metResUp_2D.clear();
  monoZhist_metResUp_2D.clear();
  monoJhist_metResDw_2D.clear();
  monoWhist_metResDw_2D.clear();
  monoZhist_metResDw_2D.clear();
  monoJhist_metUncUp_2D.clear();
  monoWhist_metUncUp_2D.clear();
  monoZhist_metUncUp_2D.clear();
  monoJhist_metUncDw_2D.clear();
  monoWhist_metUncDw_2D.clear();
  monoZhist_metUncDw_2D.clear();


  cout << "Templates for monoJ/monoV signal computed ..."<<interaction<< endl;

}

/////////////
void findAllPossibleMassPoints(vector<signalSample> & signalMassPoint, string interaction, int typeOfDMSignal){

  string baseDirMonoJet;
  string baseDirMonoW;
  string baseDirMonoZ;

  cout<<"interaction "<<interaction<<endl;

  if(interaction == "Vector"){
    baseDirMonoJet = baseInputTreePath+"/DMV_Vector/sigfilter/";
    baseDirMonoW   = baseInputTreePath+"/MonoW_Vector/sigfilter/";
    baseDirMonoZ   = baseInputTreePath+"/MonoZ_Vector/sigfilter/";
  }
  else if(interaction == "Axial"){
    baseDirMonoJet = baseInputTreePath+"/DMV_Axial/sigfilter/";
    baseDirMonoW   = baseInputTreePath+"/MonoW_Axial/sigfilter/";
    baseDirMonoZ   = baseInputTreePath+"/MonoZ_Axial/sigfilter/";
  }
  else if(interaction == "Scalar"){
    baseDirMonoJet = baseInputTreePath+"/DMS_Scalar/sigfilter/";
    baseDirMonoW   = baseInputTreePath+"/MonoW_Scalar/sigfilter/";
    baseDirMonoZ   = baseInputTreePath+"/MonoZ_Scalar/sigfilter/";
  }
  else if(interaction == "Pseudoscalar"){
    baseDirMonoJet = baseInputTreePath+"/DMS_Pseudoscalar/sigfilter/";
    baseDirMonoW   = baseInputTreePath+"/MonoW_Pseudoscalar/sigfilter/";
    baseDirMonoZ   = baseInputTreePath+"/MonoZ_Pseudoscalar/sigfilter/";
  }
  else{
    cout<<"[findAllPossibleMassPoints]: interaction type not found --> exit "<<endl;
    exit (EXIT_FAILURE);
  }

  vector<string> monoJetSamples;
  vector<string> monoWSamples;
  vector<string> monoZSamples;

  ifstream infile;
  string line;
  vector<string> seglist;

  //make the lsit of mono-jet samples
  if(typeOfDMSignal <= 1){

    system(("ls "+baseDirMonoJet+" | grep root > file_list.tmp").c_str());
    infile.open("file_list.tmp",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
        getline(infile,line);
        if(line == "") continue;
        TString fileName (line.c_str());
        seglist.clear();
	if(interaction == "Vector" or interaction == "Axial"){
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
	else{
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
        stringstream name(fileName.Data());
        string segment;
        while(getline(name, segment, '_')){
          seglist.push_back(segment);
        }
        monoJetSamples.push_back(seglist.at(seglist.size()-2)+"_"+seglist.back());
      }
    }
    infile.close();
  }

  // make the list of mono-W and mono-Z samples
  if(typeOfDMSignal == 0 or typeOfDMSignal >= 2){

    system(("ls "+baseDirMonoW+" | grep root > file_list.tmp").c_str());
    infile.open("file_list.tmp",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
        getline(infile,line);
        if(line == "") continue;
        TString fileName (line.c_str());
        seglist.clear();
	if(interaction == "Vector" or interaction == "Axial"){
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
	else{
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
        stringstream name(fileName.Data());
        string segment;
        while(getline(name, segment, '_')){
          seglist.push_back(segment);
        }
	monoWSamples.push_back(seglist.at(seglist.size()-2)+"_"+seglist.back());
      }
    }

    infile.close();

    system(("ls "+baseDirMonoZ+" | grep root > file_list.tmp").c_str());
    infile.open("file_list.tmp",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
        getline(infile,line);
        if(line == "") continue;
        TString fileName (line.c_str());
        seglist.clear();
	if(interaction == "Vector" or interaction == "Axial"){
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-0p25_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
	else{
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-madgraph.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-powheg.root","");
	  fileName.ReplaceAll("_gSM-1p0_gDM-1p0_v2_13TeV-JHUGen.root","");
	}
        stringstream name(fileName.Data());
        string segment;
        while(getline(name, segment, '_')){
          seglist.push_back(segment);
        }
        monoZSamples.push_back(seglist.at(seglist.size()-2)+"_"+seglist.back());
      }
    }

    infile.close();
  }

  // common point  and fill output vector                                                                                                                                       
  vector<string> commonPoint;
  if(typeOfDMSignal == 1)
    commonPoint = monoJetSamples;
  else if (typeOfDMSignal == 0){
    for(auto point : monoJetSamples){
      if(find(monoWSamples.begin(),monoWSamples.end(),point) != monoWSamples.end() and find(monoZSamples.begin(),monoZSamples.end(),point) != monoZSamples.end())
        commonPoint.push_back(point);
      else if(find(monoWSamples.begin(),monoWSamples.end(),point) != monoWSamples.end())
        commonPoint.push_back(point);
    }
  }

  else if(typeOfDMSignal >= 2){
    for(auto point : monoWSamples){
      if(find(monoZSamples.begin(),monoZSamples.end(),point) != monoZSamples.end())
        commonPoint.push_back(point);
      else if(interaction == "Scalar")
        commonPoint.push_back(point);
    }
  }

  system("rm file_list.tmp");

  // string manipulation
  for(auto point : commonPoint){
    stringstream name(point.c_str());
    string segment;
    vector<string> seglist;
    while(getline(name, segment, '-')){
      seglist.push_back(segment);
    }
    TString tmp (seglist.at(1));
    tmp.ReplaceAll("_Mchi","");
    signalMassPoint.push_back(signalSample(interaction,string(tmp.Data()),seglist.back()));
  }

  return;
}
